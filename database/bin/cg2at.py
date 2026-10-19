#!/usr/bin/env python3

import os, sys

# ---------------------------------------------------------------------------
# Path bootstrap: the CG2AT modules in database/bin import each other by
# their plain names (gen, g_var, ...).  Put this directory first on sys.path
# so that works however CG2AT is started:
#   ./cg2at  (symlink in a git clone)  -> database/bin/cg2at.py
#   python database/bin/cg2at.py
#   conda install (symlink into $CONDA_PREFIX/bin)
#   pip install (console script `cg2at` -> cg2at.bin.cg2at:main)
# ---------------------------------------------------------------------------
_bin_dir = os.path.dirname(os.path.realpath(__file__))
if _bin_dir not in sys.path:
    sys.path.insert(0, _bin_dir)

import numpy as np
import multiprocessing as mp
import time
import gen, gro, at_mod, at_mod_p, at_mod_np, read_in, g_var, check_library, dihres
from exceptions import CG2ATError


def _run() -> None:
    """Core pipeline.  Separated from main() so it can be wrapped cleanly."""

    g_var.script_update = '17-10-2026'
    g_var.other = {'DA': 'A', 'DG': 'G', 'DC': 'C', 'DT': 'T',
                   'DAX': 'A', 'DGX': 'G', 'DCX': 'C', 'DTX': 'T'}

    g_var.tc['i_t'] = time.time()

    # --- Initialise ---
    gen.cg2at_header()
    gen.fetch_forcefield_water_info()
    gen.check_input_flag()
    gen.correct_number_cpus()
    if g_var.get_forcefield:   ## -info, -posre and -compare do not need GROMACS
        gen.find_gromacs()
    gen.read_database_directories()
    gen.forcefield_selection()
    gen.fragment_selection()
    gen.check_water_molecules()
    dihres.check_input()

    if g_var.args.posre is not None and len(g_var.np_directories) > 0:
        check_library.add_posres_file()
        return   ## -posre is a utility mode: write the posre file, then stop
    if g_var.args.compare is not None and len(g_var.np_directories) > 0:
        check_library.compare_forcefield_to_database()
        return   ## -compare is a utility mode: report differences, then stop
    if g_var.args.info:
        gen.database_information()
    if g_var.args.v >= 1:
        print(gen.fragments_in_use())

    gen.get_termini_selections()
    gen.fetch_fragment_multi()
    gen.fetch_fragment_single()
    gen.fetch_chain_groups()
    gen.sort_swap_group()
    print(gen.print_swap_residues())

    # --- Collect inputs ---
    gro.collect_input()
    gen.flags_used()
    g_var.tc['i_t_e'] = time.time()

    # --- Read CG file ---
    box_vec_initial = read_in.read_initial_cg_pdb()
    if g_var.args.box is not None:
        print('box cutting only works for cubic boxes currently')
        g_var.box_vec, box_shift = gen.new_box_vec(box_vec_initial, g_var.args.box)
    else:
        g_var.box_vec = box_vec_initial
        box_shift = np.array([0, 0, 0])
    read_in.real_box_vectors(g_var.box_vec)
    read_in.fix_pbc(box_vec_initial, g_var.box_vec, box_shift)
    at_mod.sanity_check()

    # --- Protein ---
    g_var.tc['r_i_t'] = time.time()
    if 'PROTEIN' in g_var.cg_residues:
        g_var.coord_atomistic = at_mod_p.build_multi_residue_atomistic_system(
            g_var.cg_residues, 'PROTEIN'
        )
        if not g_var.user_at_input and g_var.args.v >= 1:
            print(gen.print_sequence_info('PROTEIN'))
        if g_var.args.refine is not None:
            dihres.compute_tdc()

        g_var.tc['p_d_n_t'] = time.time()
        if g_var.user_at_input:
            for file_num, file_name in enumerate(g_var.args.a):
                atomistic_protein_input_raw, g_var.chain_count = read_in.read_in_atomistic(
                    g_var.input_directory + 'AT_INPUT_' + str(file_num) + '.pdb'
                )
                g_var.atomistic_protein_input_raw.update(atomistic_protein_input_raw)
            read_in.duplicate_chain()
            at_mod_p.check_sequence()
            at_mod_p.align_chain_sequence('PROTEIN')
            at_mod_p.find_disulphide_bonds_user_sup()

        at_mod_p.find_disulphide_bonds_de_novo()
        g_var.coord_atomistic = at_mod_p.correct_disulphide_bonds(g_var.coord_atomistic)
        final_coordinates_atomistic_de_novo = at_mod_p.finalise_novo_atomistic(
            g_var.coord_atomistic, 'PROTEIN'
        )
        if g_var.user_at_input:
            at_mod_p.align_user_chains(final_coordinates_atomistic_de_novo)

        if not os.path.exists(g_var.working_dir + 'PROTEIN/PROTEIN_de_novo_merged.pdb'):
            gro.run_parallel_pdb2gmx_min('PROTEIN', g_var.ter_res['PROTEIN'])
            print('Merging de_novo protein chains')
            at_mod.merge_individual_chain_pdbs( 
                g_var.working_dir + 'PROTEIN/MIN/PROTEIN_de_novo', '.pdb', 'PROTEIN'
            )
        if g_var.args.refine is not None:
            dihres.write_restraints()

        if g_var.user_at_input and not os.path.exists(
                g_var.working_dir + 'PROTEIN/PROTEIN_aligned_merged.pdb'):
            print('Merging aligned protein chains')
            if g_var.args.o not in ['none', 'align']:
                at_mod.merge_individual_chain_pdbs(
                    g_var.working_dir + 'PROTEIN/MIN/PROTEIN_aligned', '.pdb', 'PROTEIN'
                )
            else:
                at_mod.merge_individual_chain_pdbs(
                    g_var.working_dir + 'PROTEIN/PROTEIN_aligned', '_gmx_checked.pdb', 'PROTEIN'
                )

    # --- Other linked residues ---
    g_var.tc['f_p_t'] = time.time()
    if 'OTHER' in g_var.cg_residues:
        g_var.other_atomistic = at_mod_p.build_multi_residue_atomistic_system(
            g_var.cg_residues, 'OTHER'
        )
        if g_var.args.v >= 1:
            print(gen.print_sequence_info('OTHER'))
        at_mod_p.finalise_novo_atomistic(g_var.other_atomistic, 'OTHER')
        gro.run_parallel_pdb2gmx_min('OTHER', g_var.ter_res['OTHER'])
        if not os.path.exists(g_var.working_dir + 'OTHER/OTHER_de_novo_merged.pdb'):
            at_mod.merge_individual_chain_pdbs(    
                g_var.working_dir + 'OTHER/MIN/OTHER_de_novo', '.pdb', 'OTHER'
            )
    g_var.tc['f_o_t'] = time.time()

    # --- Non-protein residues ---
    non_protein_types = [k for k in g_var.cg_residues if k not in ['PROTEIN', 'OTHER']]
    if non_protein_types:
        print('\nConverting the following residues: \n')
        for residue_type in non_protein_types:
            number = at_mod_np.build_atomistic_system(residue_type)
            g_var.system.update(number)

        print('\nThis may take some time...\n')
        for residue_type in non_protein_types:
            if not os.path.exists(
                    g_var.working_dir + residue_type + '/' + residue_type + '_merged.pdb'):
                print('Minimising: ' + residue_type)
                error = gro.minimise_merged(
                    residue_type,
                    g_var.working_dir + residue_type + '/' + residue_type + '_all.pdb',
                )
                if error and residue_type not in ['SOL']:
                    print('Failed to minimise as a group: ' + residue_type)
                    print('please check your input file as there is likely something wrong')

    # --- Merge system ---
    g_var.tc['n_p_t'] = time.time()
    print('Merging all residue types into single file.\n')
    gro.write_merged_topol()

    for file_name in os.listdir(g_var.merged_directory):
        if not any(f in file_name for f in [
                'steered_posre.itp', 'low_posre.itp', 'mid_posre.itp', 'high_posre.itp']):
            if file_name.endswith('.itp') or file_name.endswith('final.top'):
                gen.file_copy_and_check(
                    g_var.merged_directory + file_name, g_var.final_dir + file_name
                )

    if not os.path.exists(g_var.merged_directory + 'merged_cg2at_de_novo.pdb'):
        at_mod.merge_system_pdbs('_de_novo')

    if not os.path.exists(g_var.merged_directory + 'MIN/merged_cg2at_de_novo_minimised.pdb'):
        gro.make_min('merged_cg2at')
        gro.minimise_merged_pdbs('_de_novo')

    g_var.tc['m_t'] = time.time()
    if not os.path.exists(g_var.merged_directory + 'checked_ringed_lipid_de_novo.pdb'):
        at_mod.check_ringed_lipids(
            g_var.merged_directory + 'MIN/merged_cg2at_de_novo_minimised.pdb'
        )

    if g_var.args.o not in ['none', 'align']:
        gro.run_nvt(g_var.merged_directory + 'checked_ringed_lipid_de_novo')
    else:
        print('Completed initial minimisation, please find final de_novo system: \n'
              + g_var.final_dir + 'final_cg2at_de_novo.pdb')
        gen.file_copy_and_check(
            g_var.merged_directory + 'checked_ringed_lipid_de_novo.pdb',
            g_var.final_dir + 'final_cg2at_de_novo.pdb',
        )

    g_var.tc['eq_t'] = time.time()

    if (g_var.user_at_input and 'PROTEIN' in g_var.cg_residues
            and g_var.args.o in ['all', 'align']):
        g_var.tc['a_s'] = time.time()
        gro.create_aligned()
        g_var.tc['a_e'] = time.time()

    if not g_var.args.messy:
        gen.clean()

    print(gen.write_system_components())

    if 'PROTEIN' in g_var.cg_residues:
        at_mod_p.write_RMSD()
        if g_var.args.refine is not None:
            dihres.report()

    g_var.tc['f_t'] = time.time()
    gen.print_script_timings()


def main() -> None:
    mp.freeze_support()
    try:
        _run()
    except CG2ATError as exc:
        sys.exit(f'\nError: {exc}\n')


if __name__ == "__main__":
    main()
