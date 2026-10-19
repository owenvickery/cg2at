#!/usr/bin/env python3
"""
Trajectory-guided backbone dihedral restraints (TDC method) built into the
CG2AT pipeline.  Enabled with -refine <CG trajectory>.

    1. compute_tdc()       φ/ψ distributions are read from the CG trajectory
                           (new Martini 3 N/CA/C/O backbone), split into the
                           same protein chains as CG2AT, and scored by
                           Trajectory Dihedral Confidence (TDC).
    2. write_restraints()  after pdb2gmx, one PROTEIN_<chain>_dihre.itp is
                           written per chain with molecule-local atom indices.
                           PROTEIN_<chain>.itp includes it under
                           #ifdef DIHEDRALPOSRES (see gro.convert_topology).
    3. gro.run_nvt()       the de novo NVT is run with -DDIHEDRALPOSRES.
    4. report()            Ramachandran quality, plots and a TDC B-factor PDB
                           are written to FINAL.

The TDC scoring, restraint writer and Ramachandran analysis are reused from
traj_dih_guide.py and cg2at_refine.py (Hafez Razmazma).
"""

import os
import contextlib
import numpy as np
import g_var
from exceptions import InputError

mda, tdg, refine = None, None, None   ## heavy optional dependencies, imported in check_input()
tdc_chains = {}  ## tdc_chains[chain] = {'dist':..., 'tdc':..., 'stats':...}


def _abs(path):
    return path if os.path.isabs(path) else os.path.join(g_var.start_dir, path)


def check_input():
#### validates -refine flags and imports MDAnalysis before any work is done
    global mda, tdg, refine
    args = g_var.args
    if args.refine is None:
        return
    if args.o not in ['all', 'de_novo']:
        raise InputError('-refine applies dihedral restraints during the de novo NVT, '
                         'so requires -o all or -o de_novo (not -o '+args.o+')')
    try:
        import MDAnalysis as mda
        import traj_dih_guide as tdg
        import cg2at_refine as refine
    except ImportError as exc:
        raise InputError('-refine requires MDAnalysis and matplotlib ('+str(exc)+').\n'
                         'Install them with: pip install MDAnalysis matplotlib '
                         '(or install CG2AT with the [refine] extra)')
    args.refine = _abs(args.refine)
    args.refine_top = _abs(args.refine_top if args.refine_top is not None else args.c)
    for flag, file_name in [('-refine', args.refine), ('-refine_top', args.refine_top)]:
        if not os.path.exists(file_name):
            raise InputError('cannot find '+flag+' file: '+file_name)


def _cg_chains():
#### CG resnames per protein chain, in the order CG2AT built them
    chains = {}
    for chain in sorted(g_var.coord_atomistic):
        chains[chain] = []
        for resid in g_var.coord_atomistic[chain]:
            beads = g_var.cg_residues['PROTEIN'][resid]
            chains[chain].append(beads[next(iter(beads))]['residue_name'])
    return chains


def _backbone_residues(universe):
#### residues carrying the new Martini 3 N, CA and C backbone sites
    residues = [res for res in universe.residues
                if all(name in res.atoms.names for name in ['N', 'CA', 'C'])]
    if len(residues) == 0:
        if 'BB' in universe.atoms.names:
            raise InputError('-refine requires the new Martini 3 backbone (N, CA, C, O beads); '
                             'the trajectory topology only contains single BB beads')
        raise InputError('no protein backbone (N, CA, C) found in: '+g_var.args.refine_top)
    return residues


def _extract_phi_psi(universe, residues):
#### φ/ψ per frame for one chain, no dihedrals across chain breaks
    from MDAnalysis.lib.distances import calc_dihedrals
    def index(name, res_list):
        return np.array([res.atoms.select_atoms('name '+name)[0].index for res in res_list])
    n_res = len(residues)
    dist = {i+1: {'resname': res.resname, 'phi': [], 'psi': []} for i, res in enumerate(residues)}
    if n_res < 2:
        return dist, 0
    c_prev, n_cur, ca_cur, c_cur = index('C', residues[:-1]), index('N', residues[1:]), index('CA', residues[1:]), index('C', residues[1:])
    n_psi, ca_psi, c_psi, n_next = index('N', residues[:-1]), index('CA', residues[:-1]), index('C', residues[:-1]), index('N', residues[1:])
    args = g_var.args
    times = np.array([ts.time for ts in universe.trajectory])
    start = 0 if args.refine_b <= 0 else int(np.searchsorted(times, args.refine_b))
    stop = None if args.refine_e < 0 else int(np.searchsorted(times, args.refine_e, side='right'))
    phi, psi = [], []
    for ts in universe.trajectory[start:stop:args.refine_step]:
        pos = universe.atoms.positions
        phi.append(calc_dihedrals(pos[c_prev], pos[n_cur], pos[ca_cur], pos[c_cur], box=ts.dimensions))
        psi.append(calc_dihedrals(pos[n_psi], pos[ca_psi], pos[c_psi], pos[n_next], box=ts.dimensions))
    phi, psi = np.degrees(np.array(phi)), np.degrees(np.array(psi))
    for i in range(n_res):
        dist[i+1]['phi'] = phi[:, i-1] if i > 0 else np.array([])
        dist[i+1]['psi'] = psi[:, i] if i < n_res-1 else np.array([])
    return dist, len(phi)


def compute_tdc():
#### reads the CG trajectory and scores every protein chain
    args = g_var.args
    print('\nReading CG trajectory for dihedral restraints: '+args.refine)
    try:
        universe = mda.Universe(args.refine_top, args.refine)
    except Exception as exc:
        raise InputError('cannot load CG trajectory with topology '+args.refine_top+': '+str(exc)+
                         '\nIf -c is not the trajectory topology, supply it with -refine_top')
    residues = _backbone_residues(universe)
    cg_chains = _cg_chains()
    n_cg = sum(len(resnames) for resnames in cg_chains.values())
    if len(residues) != n_cg:
        raise InputError('CG trajectory has '+str(len(residues))+' backbone residues but the CG input has '+
                         str(n_cg)+' protein residues. Use -refine_top to supply the trajectory topology')
    log = g_var.working_dir+'PROTEIN/TDC_scores.dat'
    pos, n_conf = 0, {'HIGH': 0, 'MEDIUM': 0, 'LOW': 0}
    with open(log, 'w') as log_out:
        for chain, resnames in cg_chains.items():
            chain_res = residues[pos:pos+len(resnames)]
            pos += len(resnames)
            mismatch = sum(1 for res, resname in zip(chain_res, resnames) if res.resname != resname)
            if mismatch and args.v >= 1:
                print('Chain '+str(chain)+': '+str(mismatch)+' residue names differ between trajectory and CG input')
            dist, n_frames = _extract_phi_psi(universe, chain_res)
            global_tdc, global_weight, tdc, weight, stats = tdg.compute_TDC(dist)
            tdc_chains[chain] = {'dist': dist, 'tdc': tdc, 'stats': stats}
            for value in tdc.values():
                n_conf['HIGH' if value >= 70 else 'MEDIUM' if value >= args.refine_tdc else 'LOW'] += 1
            log_out.write('\nChain '+str(chain)+'\n')
            with contextlib.redirect_stdout(log_out):
                tdg.print_TDC_table(dist, tdc, weight, stats, global_tdc, global_weight)
    all_tdc = [value for chain in tdc_chains.values() for value in chain['tdc'].values()]
    print('{0} frames, {1} chains, global TDC {2:.1f}/100'.format(n_frames, len(tdc_chains), np.mean(all_tdc)))
    print('Residues HIGH (TDC>=70): {0}  MEDIUM: {1}  LOW (TDC<{2:.0f}, unrestrained): {3}'.format(
          n_conf['HIGH'], n_conf['MEDIUM'], args.refine_tdc, n_conf['LOW']))
    print('Per-residue TDC scores saved in: '+log)


def write_restraints():
#### one dihedral restraint itp per chain, indexed on the de novo pdb2gmx output
    args = g_var.args
    log = g_var.working_dir+'PROTEIN/dihedral_restraints.log'
    total = 0
    with open(log, 'w') as log_out, contextlib.redirect_stdout(log_out):
        for chain, data in tdc_chains.items():
            print('\nChain '+str(chain))
            total += tdg.write_dihedral_restraints(
                data['dist'], data['tdc'], data['stats'],
                aa_pdb=g_var.working_dir+'PROTEIN/PROTEIN_de_novo_'+str(chain)+'_gmx.pdb',
                output_itp=g_var.working_dir+'PROTEIN/PROTEIN_'+str(chain)+'_dihre.itp',
                k_max=args.refine_kmax, sigma=args.refine_sigma,
                tdc_low=args.refine_tdc, fixed_dphi=args.refine_dphi)
    print('Wrote '+str(total)+' backbone dihedral restraints across '+str(len(tdc_chains))+' chains (details: '+log+')')


def report():
#### Ramachandran quality of the minimised de novo vs the restrained NVT output
    final = g_var.final_dir
    before = g_var.merged_directory+'checked_ringed_lipid_de_novo.pdb'
    after = final+'final_cg2at_de_novo.pdb'
    aligned = final+'final_cg2at_aligned.pdb' if os.path.exists(final+'final_cg2at_aligned.pdb') else None
    if not os.path.exists(after):
        print('Cannot find '+after+', skipping dihedral restraint report')
        return
    print('\n{0:-<100}\n{1:^100}\n{0:-<100}'.format('', 'Trajectory-guided dihedral restraints (TDC)'))
    print('"De novo" = minimised fragment fit before NVT; "Restrained" = final_cg2at_de_novo.pdb\n')

    combined, pos = {}, 0
    for data in tdc_chains.values():
        for resid in sorted(data['tdc']):
            combined[pos] = data['tdc'][resid]
            pos += 1
    refine.write_tdc_bfactor_pdb(after, combined, final+'final_cg2at_de_novo_TDC.pdb')

    cg_ref = g_var.input_directory+'CG_INPUT.pdb'
    rmsd_data = refine.print_rmsd_report(before, aligned, after, cg_ref if os.path.exists(cg_ref) else None)

    Z, phi_bins, psi_bins = refine._load_mda_reference()
    rama500 = {}
    rama500_dir = refine._find_rama500_dir()
    if rama500_dir:
        for key, cfg in refine._RAMA500_SETTINGS.items():
            if not os.path.exists(os.path.join(rama500_dir, cfg['filename'])):
                continue
            Z_r, mid = refine._load_rama500_file(os.path.join(rama500_dir, cfg['filename']))
            rama500[key] = (Z_r, mid, cfg['levels'])

    def assess(pdb):
        if pdb is None:
            return None, None, None
        return (refine.assess_ramachandran(pdb, Z, phi_bins, psi_bins, rama500.get('general')),
                refine.assess_ramachandran_gly(pdb, rama500.get('gly')),
                refine.assess_ramachandran_pro(pdb, rama500.get('pro')))

    dn, rest, al = assess(before), assess(after), assess(aligned)
    for label, idx in [('General', 0), ('GLY', 1), ('PRO', 2)]:
        refine.print_rama_summary('De novo — '+label, dn[idx])
        refine.print_rama_summary('Restrained — '+label, rest[idx])
        if al[idx]:
            refine.print_rama_summary('Aligned — '+label, al[idx])
    refine.write_quality_dat(final+'structure_quality_restrained.dat', dn[0], rest[0], al[0],
                             dn[1], rest[1], al[1], stats_pro_dn=dn[2], stats_pro_rest=rest[2],
                             stats_pro_al=al[2], rmsd_data=rmsd_data)
    refine.plot_ramachandran(before, after, aligned, final+'ramachandran_restrained.png',
                             Z=Z, phi_bins=phi_bins, psi_bins=psi_bins,
                             rama500_general_data=rama500.get('general'))
    refine.plot_ramachandran_gly_pro(before, after, aligned, final+'ramachandran_gly_pro_restrained.png',
                                     rama500_gly_data=rama500.get('gly'), rama500_pro_data=rama500.get('pro'))
    print('\nDihedral restraint outputs saved in '+final+':')
    print('  structure_quality_restrained.dat, ramachandran_restrained.png,')
    print('  ramachandran_gly_pro_restrained.png, final_cg2at_de_novo_TDC.pdb (TDC in B-factor column)')
