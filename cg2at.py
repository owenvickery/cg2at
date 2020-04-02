#!/usr/bin/env python3

import os, sys
import numpy as np
from shutil import copyfile
import time
import re
import multiprocessing as mp
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/database/script_files')
import gen, gro, at_mod, at_mod_p, at_mod_np, read_in, g_var, f_loc, at2cg

time_counter = {}
time_counter['i_t']=time.time()

#### collects initial structures into INPUT folder

user_at_input = gro.collect_input(g_var.c, g_var.a)

#### saves flags used into INPUT folder
gen.flags_used()

if g_var.at2cg:
    print('\nThis script is now hopefully doing the following (Good luck):\n\nReading in your AT representation \n')
    at_residues, box_vec =read_in.read_initial_at_pdb()
    print('converting your atomistic system to coarse grain (Time for a Becherovka)\n')
    cg_residues = at2cg.convert_AT2CG(at_residues, box_vec)
    at2cg.write_topology(cg_residues)
    at2cg.print_system_info(cg_residues)
    print('To remake the topology of your martini protein, a copy of martinise is in the scripts directory.\n')
    print(g_var.scripts_dir+'martinise.py -f '+g_var.input_directory+'conversion_input.pdb -o '+g_var.final_dir+'/protein.top')
    print('you\'ll need to added the flags you require ')

else:

    print('\nThis script is now hopefully doing the following (Good luck):\n\nReading in your CG representation\n')

    #### reads in CG file and separates into residue types
    cg_residues, box_vec_initial = read_in.read_initial_cg_pdb()

    #### box size update 
    if g_var.box != None:
        box_vec, box_shift = gen.new_box_vec(box_vec_initial, g_var.box)
    else:
        box_vec=box_vec_initial
        box_shift=np.array([0,0,0])

    #### simple pbc fix and residue truncation if required
    cg_residues =read_in.fix_pbc(cg_residues, box_vec_initial, box_vec, box_shift)
    if 'PROTEIN' in cg_residues:
        protein = True
    else:
        protein = False
    #### checks if fragment database is correct
    at_mod.sanity_check(cg_residues)

    time_counter['r_i_t']=time.time()
    system={}
    ### convert protein to atomistic representation
    if protein:
        if user_at_input:
            atomistic_protein_input_raw, chain_count = at_mod_p.read_in_atomistic(g_var.input_directory+'AT_input.pdb')  ## reads in user structure
        p_system, backbone_coords, coordinates_atomistic, sequence=at_mod_p.build_protein_atomistic_system(cg_residues['PROTEIN'])

        if not user_at_input and g_var.v >= 1:
            print('coarse grain protein sequence:\n')
            for index in sequence:
                print('chain:', index,sequence[index], '\n') 
        system['PROTEIN']=p_system['PROTEIN']
        time_counter['p_d_n_t']=time.time()
        if user_at_input:
            seq_user = at_mod_p.check_sequence(atomistic_protein_input_raw, chain_count)
            atomistic_protein_input, group_chain, user_at_input = at_mod_p.align_chains(atomistic_protein_input_raw, seq_user, sequence)
            user_cys_bond = at_mod_p.find_disulphide_bonds_user_sup(atomistic_protein_input)
        else:
            user_cys_bond = {}
        user_cys_bond = at_mod_p.find_disulphide_bonds_de_novo(coordinates_atomistic, user_cys_bond)
        coordinates_atomistic = at_mod_p.correct_disulphide_bonds(coordinates_atomistic, user_cys_bond)
        final_coordinates_atomistic = at_mod_p.finalise_novo_atomistic(coordinates_atomistic, cg_residues['PROTEIN'], box_vec)
        if user_at_input:
            atomistic_protein_centered, cg_com = at_mod_p.center_atomistic(atomistic_protein_input, backbone_coords, group_chain) ## centers each monomer by center of mass
            at_com_group, cg_com_group = at_mod_p.rotate_protein_monomers(atomistic_protein_centered, final_coordinates_atomistic, backbone_coords, cg_com, box_vec, group_chain) ## rigid fits each monomer
            final_user_supplied_coord = at_mod_p.apply_rotations_to_chains(final_coordinates_atomistic, atomistic_protein_centered, at_com_group,cg_com_group,cg_com, box_vec, group_chain)
            final_user_supplied_coord = at_mod_p.correct_disulphide_bonds(final_user_supplied_coord, user_cys_bond)
            at_mod_p.write_user_chains_to_pdb(final_user_supplied_coord, box_vec)
        #### minimise each protein chain
        print('Minimising '+str(p_system['PROTEIN'])+' protein chains')
        gro.minimise_protein(system['PROTEIN'], p_system, user_at_input, box_vec)
        #### read in minimised de novo protein chains and merges chains
        merge_de_novo = at_mod_p.read_in_protein_pdbs(system['PROTEIN'], g_var.working_dir+'PROTEIN/MIN/PROTEIN_de_novo', '.pdb')
        at_mod_p.write_merged_pdb(merge_de_novo, '_de_novo', box_vec)
        #### runs steered MD on user supplied protein chains
        if user_at_input:
            if g_var.o in ['all', 'steer']:
                print('Running steered MD on input atomistic structure')
            #### runs steered MD on atomistic structure on CA and CB atoms
                for chain in range(system['PROTEIN']):
                    gro.steered_md_atomistic_to_cg_coord(chain)
                #### read in minimised user supplied protein chains and merges chains
                merge_at_user = at_mod_p.read_in_protein_pdbs(system['PROTEIN'], g_var.working_dir+'PROTEIN/STEERED_MD/PROTEIN_steered', '.pdb')
                at_mod_p.write_merged_pdb(merge_at_user, '_steered', box_vec)
            if g_var.o in ['all', 'align']:
                merge_at_user_no_steer = at_mod_p.read_in_protein_pdbs(system['PROTEIN'], g_var.working_dir+'PROTEIN/PROTEIN_aligned', '_gmx.pdb')
                at_mod_p.write_merged_pdb(merge_at_user_no_steer, '_aligned', box_vec)

    time_counter['f_p_t']=time.time()

    #### converts non protein residues into atomistic (runs on all cores)
    if len([key for value, key in enumerate(cg_residues) if key not in ['PROTEIN']]) > 0:
        print('\nConverting the following residues concurrently: ')
        np_system={}
        pool = mp.Pool(mp.cpu_count())
        pool_process = pool.starmap_async(at_mod_np.build_atomistic_system, [(cg_residues, residue_type, box_vec) for residue_type in [key for value, key in enumerate(cg_residues) if key not in ['PROTEIN']]]).get()          ## minimisation grompp parallised  
        pool.close()
        for residue_type in pool_process:
            np_system.update(residue_type)
        #### minimises each residue separately
        print('\nThis may take some time....(probably time for a coffee)\n')
        for residue_type in [key for value, key in enumerate(cg_residues) if key not in ['PROTEIN', 'ION']]:
            print('Processing individual residues: '+residue_type)
            gro.non_protein_minimise(np_system[residue_type], residue_type)
            at_mod_np.merge_minimised(residue_type, np_system, box_vec)
            print('Minimising merged: '+residue_type+'\n')
            gro.minimise_merged(residue_type, np_system)
        system.update(np_system)

    time_counter['n_p_t']=time.time()

    #### creates merged folder
    print('\nMerging all residue types to single file. (Or possibly tea)\n')

    if len(system)>0:
    #### make final topology in merged directory
        gro.write_merged_topol(system, '_de_novo')
    #### copies all itp files and topologies from whereever they are stored
        for file_name in os.listdir(g_var.merged_directory):
            if not any(s in file_name for s in ['steered_posre.itp', 'low_posre.itp', 'high_posre.itp']):
                if file_name.endswith('.itp') or file_name.endswith('final.top') :
                   gen.file_copy_and_check(g_var.merged_directory+file_name, g_var.final_dir+file_name)
    #### make minimisation directory
        gro.make_min('merged_cg2at')
    #### merges provided atomistic protein and residues types into a single pdb file into merged directory
        at_mod.merge_system_pdbs(system, '_de_novo', cg_residues, box_vec)
        gro.minimise_merged_pdbs(system, '_de_novo')
        gro.run_nvt(g_var.merged_directory+'MIN/merged_cg2at_de_novo_minimised.pdb', protein)
        ringed_lipids = at_mod.check_ringed_lipids(g_var.merged_directory+'NVT/merged_cg2at_de_novo_nvt.pdb', box_vec)
        if len(system) > 1 and g_var.alchembed and len(ringed_lipids) > 0 and protein:
            gro.alchembed(system['PROTEIN'], 'de_novo')  
            ringed_lipids = at_mod.check_ringed_lipids(g_var.final_dir+'final_cg2at_de_novo.pdb', box_vec)
            if len(ringed_lipids) > 0:
                print('Check final output as alchembed cannot fix ringed lipid: ', ringed_lipids)
        else:
            gro.run_npt(g_var.merged_directory+'NVT/merged_cg2at_de_novo_nvt', protein)
        time_counter['m_t']=time.time()
        if user_at_input and protein:
            print()
            if g_var.o in ['all', 'steer']:
                print('Creating steered system')
                time_counter['s_s']=time.time()
                at_mod.merge_system_pdbs(system, '_steered', cg_residues, box_vec)
                gro.reverse_steer('steered', 'low', g_var.merged_directory+'/NPT/merged_cg2at_de_novo_npt')
                gro.reverse_steer('steered', 'mid', g_var.merged_directory+'REVERSE_STEER/merged_cg2at_steered_reverse_steer_low')
                gro.reverse_steer('steered', 'high', g_var.merged_directory+'REVERSE_STEER/merged_cg2at_steered_reverse_steer_mid')
                gen.file_copy_and_check(g_var.merged_directory+'REVERSE_STEER/merged_cg2at_steered_reverse_steer_high.pdb', g_var.final_dir+'final_cg2at_steered.pdb')
                time_counter['s_e']=time.time()
            if g_var.o in ['all', 'align']:   
                print('Creating aligned system') 
                time_counter['a_s']=time.time()
                at_mod.merge_system_pdbs(system, '_aligned', cg_residues, box_vec)
                gro.reverse_steer('aligned', 'low', g_var.merged_directory+'/NPT/merged_cg2at_de_novo_npt')
                gro.reverse_steer('aligned', 'mid', g_var.merged_directory+'REVERSE_STEER/merged_cg2at_aligned_reverse_steer_low')
                gro.reverse_steer('aligned', 'high', g_var.merged_directory+'REVERSE_STEER/merged_cg2at_aligned_reverse_steer_mid')
                gen.file_copy_and_check(g_var.merged_directory+'REVERSE_STEER/merged_cg2at_aligned_reverse_steer_high.pdb', g_var.final_dir+'final_cg2at_aligned.pdb')
                time_counter['a_e']=time.time()

    if protein:
    #### calculates final RMS
        RMSD={}
        de_novo_atoms, chain_count = at_mod_p.read_in_atomistic(g_var.final_dir+'final_cg2at_de_novo.pdb')
        if chain_count != system['PROTEIN']:
            sys.exit('number of chains in atomistic protein input ('+str(chain_count)+') does not match CG representation ('+str(system['PROTEIN'])+')')
        RMSD['de novo '] = at_mod_p.RMSD_measure(de_novo_atoms, system,backbone_coords)

        if user_at_input and protein:
            if g_var.o in ['all', 'steer']:
                at_input_atoms, chain_count = at_mod_p.read_in_atomistic(g_var.final_dir+'final_cg2at_steered.pdb')
                RMSD['at steered'] = at_mod_p.RMSD_measure(at_input_atoms, system,backbone_coords)   
            if g_var.o in ['all', 'align']: 
                at_input_atoms, chain_count = at_mod_p.read_in_atomistic(g_var.final_dir+'final_cg2at_aligned.pdb')
                RMSD['at aligned'] = at_mod_p.RMSD_measure(at_input_atoms, system,backbone_coords)   
        with open(g_var.final_dir+'structure_quality.dat', 'w') as qual_out:   
            qual_out.write('\n{0:^10}{1:^25}{2:^10}\n'.format('output ','chain','RMSD ('+chr(197)+')'))
            qual_out.write('{0:^10}{1:^25}{2:^10}\n'.format('-------','-----','---------'))
            print('\n{0:^10}{1:^25}{2:^10}'.format('output ','chain','RMSD ('+chr(197)+')'))
            print('{0:^10}{1:^25}{2:^10}'.format('-------','-----','---------'))
            for rmsd in RMSD:
                for chain in RMSD[rmsd]:
                    qual_out.write('{0:^10}{1:^25}{2:^10}\n'.format(rmsd, str(chain), float(RMSD[rmsd][chain])))
                    print('{0:^10}{1:^25}{2:^10}'.format(rmsd, str(chain), float(RMSD[rmsd][chain])))

    #### removes temp file from script, anything with temp in really
    if not g_var.messy:
        gen.clean(cg_residues)

    time_counter['f_t']=time.time()

    #### prints out system information
    print('\n{:-<100}'.format(''))
    print('{0:^100}'.format('Script has completed, time for a beer'))
    print('\n{0:^20}{1:^10}'.format('molecules','number'))
    print('{0:^20}{1:^10}'.format('---------','------'))
    for section in system:
        print('{0:^20}{1:^10}'.format(section, system[section]))

    #### prints out script timings for each section
    if g_var.v >= 1:
        gen.print_script_timings(time_counter, system, user_at_input)



