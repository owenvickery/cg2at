#!/usr/bin/env python3

import os, sys
from time import gmtime, strftime
import argparse
from pathlib import Path

# Single source of truth for the version string — referenced by argparse
# below and by cg2at_header() in gen.py.
VERSION: str = '1.0.0'

if __name__ == "g_var" and 'start_dir' not in locals():
    parser = argparse.ArgumentParser(
        prog='CG2AT2',
        description=(
            'CG2AT2: fragment-based back-mapping of coarse-grained systems to atomistic resolution.\n'
            'Requires a GROMACS installation. Documentation: https://github.com/pstansfeld/cg2at\n'
            'Original CG2AT by Owen Vickery: https://github.com/owenvickery/cg2at'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'Examples:\n'
            '  Minimal:     cg2at -c system.gro\n'
            '  Recommended: cg2at -c system.gro -a protein.pdb\n'
            '  Automated:   cg2at -c system.gro -a protein.pdb -ff charmm36-cg2at -fg martini_3-0_charmm36 -w tip3p\n'
            '  Database:    cg2at -info   (lists every force field and fragment library)\n'
        ),
    )

    # --- Information ---
    parser.add_argument('-info',
        help='Print available force fields and fragment libraries, then exit.',
        action='store_true')
    parser.add_argument('-version',
        action='version',
        version=f'%(prog)s {VERSION}')

    # --- Input files ---
    input_grp = parser.add_argument_group('input files')
    input_grp.add_argument('-c',
        help='Coarse-grained input structure (required).',
        metavar='FILE', type=str)
    input_grp.add_argument('-a',
        help='Atomistic input structure(s). Supplying the original atomistic file improves conversion quality. '
             'Multiple files are accepted.',
        metavar='FILE', type=str, nargs='*')

    # --- Force field and fragment library ---
    ff_grp = parser.add_argument_group('force field and fragment library')
    ff_grp.add_argument('-ff',
        help='Force field to use, e.g. charmm36-cg2at, amber14sb or oplsaa (see -info). '
             'If omitted, CG2AT2 will prompt interactively.',
        metavar='NAME', type=str)
    ff_grp.add_argument('-fg',
        help='Fragment library (or libraries, in priority order), e.g. martini_2-2_charmm36. '
             'If omitted, CG2AT2 will prompt interactively.',
        metavar='NAME', type=str, nargs='*')
    ff_grp.add_argument('-w',
        help='Water model to use, e.g. tip3p, tip4p, spc, spce. '
             'If omitted, CG2AT2 will prompt interactively.',
        metavar='MODEL', type=str)

    # --- Output ---
    out_grp = parser.add_argument_group('output')
    out_grp.add_argument('-o',
        help=(
            'Controls which final structures are produced (default: all).\n'
            '  none    – minimised de novo structure only, no MD\n'
            '  de_novo – as "none" plus a short 5 ps NVT run\n'
            '  align   – minimised de novo morphed to the user structure via steered MD\n'
            '  all     – both de_novo and align outputs'
        ),
        default='all', type=str, choices=['all', 'align', 'de_novo', 'none'])
    out_grp.add_argument('-loc',
        help='Name of the output folder. Default: CG2AT_<timestamp>. '
             'Re-using an existing folder lets you resume a failed run.',
        metavar='FOLDER', type=str,
        default='CG2AT_' + strftime("%Y-%m-%d_%H-%M-%S", gmtime()))

    # --- Protein chain handling ---
    prot_grp = parser.add_argument_group('protein chains')
    prot_grp.add_argument('-d',
        help='Duplicate atomistic chains. Format: CHAIN:COPIES [CHAIN:COPIES ...]. '
             'Example: "0:3 1:3" creates 3 copies of chains 0 and 1.',
        type=str, nargs='*', default=[], metavar='CHAIN:COPIES')
    prot_grp.add_argument('-group',
        help='Treat atomistic chains as rigid bodies during fitting. '
             'Provide comma-separated chain indices per group, separated by spaces, '
             'e.g. "0,2 1,3". Use "all" to fit the entire structure as one body, '
             'or "chain" to group by CG chain.',
        type=str, nargs='*', metavar='0,1')
    prot_grp.add_argument('-ter',
        help='Interactively choose N- and C-terminal capping groups for each chain.',
        action='store_true')
    prot_grp.add_argument('-nt',
        help='Apply a neutral N-terminus to all chains.',
        action='store_true')
    prot_grp.add_argument('-ct',
        help='Apply a neutral C-terminus to all chains.',
        action='store_true')

    # --- Disulphide bonds ---
    cys_grp = parser.add_argument_group('disulphide bonds')
    cys_grp.add_argument('-cys',
        help='Maximum S–S distance (Å) for disulphide bond detection (default: 7). '
             'Increase if bonds are missed; Martini S–S distances can reach ~10 Å.',
        metavar='DIST', type=float, default=7)
    cys_grp.add_argument('-silent',
        help='Automatically accept all detected disulphide bonds without prompting.',
        action='store_true')

    # --- Conversion tuning ---
    tune_grp = parser.add_argument_group('conversion tuning')
    tune_grp.add_argument('-mod',
        help='Treat every fragment bead independently rather than fitting grouped beads together. '
             'Grouping is on by default and generally gives better geometry.',
        action='store_true')
    tune_grp.add_argument('-sf',
        help='Scale factor applied to fragments before fitting (default: 0.9, i.e. 90%%). '
             'Reducing this lowers the chance of atom clashes before minimisation.',
        metavar='FACTOR', type=float, default=0.9)
    tune_grp.add_argument('-ov',
        help='Minimum allowed distance between atoms (Å) before overlap is flagged (default: 0.3).',
        metavar='DIST', type=float, default=0.3)
    tune_grp.add_argument('-swap',
        help='Swap residues or beads during conversion. '
             'Format: FROM[,BEAD]:TO[,BEAD][:RESID_RANGE]. '
             'Example: "POPC,NC3:POPG,GL0" or "ASP:ASN:0-10,30-40". '
             'Use "skip" as the target to remove a residue or bead entirely.',
        metavar='RULE', type=str, nargs='*')
    tune_grp.add_argument('-box',
        help='Override the periodic box dimensions (Å, cubic boxes only). '
             'Supply three values; use 0 to preserve the original size on that axis. '
             'Example: "-box 100 100 0" resizes X and Y only.',
        metavar='LEN', type=float, nargs=3)
    tune_grp.add_argument('-vs',
        help='Add virtual sites to protein hydrogen atoms.',
        action='store_true')
    tune_grp.add_argument('-disre',
        help='Apply a backbone H-bond distance restraint matrix during the de novo NVT. '
             'Off by default; applied only when -disre is given and -a is supplied.',
        action='store_true')

    # --- Trajectory-guided dihedral restraints ---
    ref_grp = parser.add_argument_group('trajectory-guided dihedral restraints (requires the new Martini 3 N/CA/C/O backbone)')
    ref_grp.add_argument('-refine',
        help='CG trajectory (e.g. .xtc). Backbone φ/ψ distributions are scored by Trajectory Dihedral '
             'Confidence (TDC) and applied as dihedral restraints during the de novo NVT. '
             'Requires the optional refine dependencies (MDAnalysis, matplotlib).',
        metavar='TRAJ', type=str)
    ref_grp.add_argument('-refine_top',
        help='Topology for the -refine trajectory (.tpr/.gro/.pdb). Default: the -c file.',
        metavar='FILE', type=str)
    ref_grp.add_argument('-refine_step',
        help='Use every Nth trajectory frame (default: 10).',
        metavar='N', type=int, default=10)
    ref_grp.add_argument('-refine_b',
        help='First trajectory time to use, in ps (default: 0).',
        metavar='PS', type=float, default=0)
    ref_grp.add_argument('-refine_e',
        help='Last trajectory time to use, in ps (default: -1, end of trajectory).',
        metavar='PS', type=float, default=-1)
    ref_grp.add_argument('-refine_kmax',
        help='Maximum dihedral force constant, kJ/mol/rad² (default: 500).',
        metavar='K', type=float, default=500)
    ref_grp.add_argument('-refine_sigma',
        help='Force constant decay scale k = kmax·exp(-σ/sigma), in degrees (default: 30).',
        metavar='DEG', type=float, default=30)
    ref_grp.add_argument('-refine_tdc',
        help='Minimum TDC score for a dihedral to be restrained (default: 40).',
        metavar='TDC', type=float, default=40)
    ref_grp.add_argument('-refine_dphi',
        help='Fixed restraint tolerance in degrees (default: adaptive, 0.3·σ clipped to 5–20°).',
        metavar='DEG', type=float)
    ref_grp.add_argument('-refine_nvt',
        help='Length of the restrained de novo NVT in 1 fs steps (default: 50000; unrestrained NVT is 5000).',
        metavar='STEPS', type=int, default=50000)

    # --- Utilities ---
    util_grp = parser.add_argument_group('utilities')
    util_grp.add_argument('-gmx',
        help='GROMACS executable to use (default: gmx). '
             'Set this if you have multiple GROMACS versions installed, e.g. gmx_avx.',
        metavar='EXE', type=str)
    util_grp.add_argument('-posre',
        help='Generate a position restraint .itp file for the named non-protein residue, then exit.',
        metavar='RESNAME', type=str)
    util_grp.add_argument('-compare',
        help='Compare all residues in an .itp file against the fragment database, then exit.',
        metavar='FILE.itp', type=str)
    util_grp.add_argument('-messy',
        help='Keep all temporary files after the run (useful for debugging).',
        action='store_true')
    util_grp.add_argument('-v',
        help='Increase output verbosity. Stack up to three times (-vvv) for maximum detail.',
        action='count', default=0)
    util_grp.add_argument('-ncpus',
        help=argparse.SUPPRESS,  # Deprecated; hidden from help
        type=int)
    args = parser.parse_args()
    opt = vars(args)
    opt['input']=os.path.abspath(sys.argv[0])+' '+''.join([ i+' ' for i in sys.argv[1:]])+'\n'

    #### virtual site information
    if args.vs:
        vs = '-vsite h'
        sf = args.sf-0.1
    else:
        vs = ''
        sf=args.sf

    ### hardcoded variables for use elsewhere in the script

    topology = {'HYDRATION':'','ALT_RES':'', 'C_TERMINAL':'default', 'N_TERMINAL':'default', 
                'CHIRAL':{'atoms':[]}, 'GROUPS':{'group_max':1}, 'CONNECT':{'atoms':{}}}

    box_line="CRYST1 %8.3f %8.3f %8.3f %8.2f %8.2f %8.2f P 1           1\n"

    pdbline = "ATOM  %5d %4s %4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f"

    os.environ['GMX_SUPPRESS_DUMP'] = '1'  ## prevent gromacs filling the file system with step files

    aas = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 
                'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 
                'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

    ### CG2AT folder locations

    start_dir       = os.getcwd()+'/'  ### initial working directory
    working_dir     = os.path.join(os.getcwd(), args.loc, '')   ### working directory (-loc may be relative or absolute)
    final_dir       = working_dir+'FINAL/'  ### final directory for run files
    input_directory = working_dir+'INPUT/'  ### contains input run files
    merged_directory = working_dir+'MERGED/'  ### contains run files
    scripts_dir     = os.path.dirname(os.path.realpath(__file__))+'/' ### contains script files
    database_dir    = str(Path(*Path(scripts_dir).parts[:-1]))+'/' ### contains database files
    box_vec = ''
    user_at_input = False
    ter_res = {}  ## contains the chain termini residue info e.g.  0:[ARG, PRO]
    system = {}  ## number of system components e.g. PROTEIN:2 POPE:10, POPG:20
    backbone_coords = {} ## CG coordinates of the backbone beads
    coord_atomistic = {} ## de_novo atomisitic information e.g. coord_atomistic[chain_count][residue_number][atom][info....]
    user_cys_bond = {} ## contains resid of disulphide bonds e.g. user_cys_bond[chain][[cys_resid,cys_resid], [cys_resid,cys_resid]])
    cg_residues = {} ## dictionary of CG beads eg cg_residues[residue type(POPE)][resid(1)][bead name(BB)][residue_name(PO4)/coordinates(coord)]
    seq_cg = {} ## CG sequence e.g. seq_cg[chain][sequence]
    seq_at = {} ## user AT sequence e.g. seq_at[chain][sequence]
    seq_cg_other = {} ## CG seq for linked non protein residues e.g. DNA
    tc = {} ## contains script timings
    atomistic_protein_input_raw = {} ## Raw user atomistic info coord_atomistic[chain_count][residue_number][atom][info....]
    atomistic_protein_input_aligned = {}
    chain_count = 0 ## number of user atomistic chains
    other_atomistic ={} ## 
    forcefield_available, fragments_available = '',''
    forcefield_location, forcefield = '','' ## forcefield info
    np_residues, p_residues, mod_residues, o_residues, sol_residues, ion_residues = [],[],[],[],[],[] ## fragment names
    np_directories, p_directories, mod_directories, o_directories, sol_directories, ion_directories = [],[],[],[],[],[]  ## fragment locations
    database_locations = [] ## grouped directories
    water, water_info = [],[] ## water information
    swap_dict ={} ## CG residue swap
    res_top, sorted_connect, hydrogen, heavy_bond, ions, at_mass = {},{},{},{},[],{} ### topology information
    group_chains = None
    np_blocks = {}
    skip_disul = {}  ## contains info on whether user supplied complete chain information   
    alt_res_name = {} ## alternative residue names
    hydration = {}
    get_forcefield = True
    cg_chain_group = {}
