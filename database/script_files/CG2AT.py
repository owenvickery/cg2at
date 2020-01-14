#!/usr/bin/env python3

import os, sys
import numpy as np
from shutil import copyfile
import time
import multiprocessing as mp
import gen, gro, at_mod, at_mod_p, at_mod_np, read_in, g_var, f_loc

def CG2AT_run(user_at_input):
    gen.flags_used()

    time_counter = {}
    time_counter['i_t']=time.time()

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

    #### checks if fragment database is correct

    at_mod.sanity_check(cg_residues)

    time_counter['r_i_t']=time.time()

    system={}

    ### convert protein to atomistic representation
    if 'PROTEIN' in cg_residues:
        p_system, backbone_coords, final_coordinates_atomistic, sequence=at_mod_p.build_protein_atomistic_system(cg_residues['PROTEIN'], box_vec)
        system['PROTEIN']=p_system['PROTEIN']
        time_counter['p_d_n_t']=time.time()
        #### reads in user supplied atomistic structure 
        if user_at_input and 'PROTEIN' in system:
            atomistic_protein_input = at_mod_p.read_in_atomistic(g_var.input_directory+'AT_input.pdb', system['PROTEIN'], sequence, True)  ## reads in user structure
            atomistic_protein_centered, cg_com = at_mod_p.center_atomistic(atomistic_protein_input, backbone_coords) ## centers each monomer by center of mass
            at_mod_p.rotate_protein_monomers(atomistic_protein_centered, final_coordinates_atomistic, backbone_coords, cg_com, box_vec) ## rigid fits each monomer
        #### minimise each protein chain
        gro.minimise_protein(system['PROTEIN'], p_system, user_at_input)
        #### read in minimised de novo protein chains and merges chains
        merge_de_novo = at_mod_p.read_in_protein_pdbs(system['PROTEIN'], g_var.working_dir+'PROTEIN/min/PROTEIN_novo', '.pdb')
        at_mod_p.write_merged_pdb(merge_de_novo, '_novo', box_vec)
        #### runs steered MD on user supplied protein chains
        if user_at_input and 'PROTEIN' in system:
            print('\tRunning steered MD on input atomistic structure\n')
        #### runs steered MD on atomistic structure on CA and CB atoms
            for chain in range(system['PROTEIN']):
                gro.steered_md_atomistic_to_cg_coord(chain)
            #### read in minimised user supplied protein chains and merges chains
            merge_at_user = at_mod_p.read_in_protein_pdbs(system['PROTEIN'], g_var.working_dir+'PROTEIN/steered_md/PROTEIN_at_rep_user_supplied', '.pdb')
            at_mod_p.write_merged_pdb(merge_at_user, '_at_rep_user_supplied', box_vec)
            merge_at_user_no_steer = at_mod_p.read_in_protein_pdbs(system['PROTEIN'], g_var.working_dir+'PROTEIN/PROTEIN_at_rep_user_supplied', '_gmx.pdb')
            at_mod_p.write_merged_pdb(merge_at_user_no_steer, '_no_steered', box_vec)

    time_counter['f_p_t']=time.time()

    #### converts non protein residues into atomistic (runs on all cores)
    if len([key for value, key in enumerate(cg_residues) if key not in ['PROTEIN']]) > 0:
        np_system={}
        pool = mp.Pool(mp.cpu_count())
        pool_process = pool.starmap_async(at_mod_np.build_atomistic_system, [(cg_residues, residue_type, box_vec) for residue_type in [key for value, key in enumerate(cg_residues) if key not in ['PROTEIN']]]).get()          ## minimisation grompp parallised  
        pool.close()
        for residue_type in pool_process:
            np_system.update(residue_type)
        #### minimises each residue separately
        print('\nThis may take some time....(probably time for a coffee)\n')
        for residue_type in [key for value, key in enumerate(cg_residues) if key not in ['PROTEIN', 'ION']]:
            print('Minimising individual residues: '+residue_type)
            gro.non_protein_minimise(np_system[residue_type], residue_type)
            at_mod_np.merge_minimised(residue_type, np_system, box_vec)
            print('Minimising merged: '+residue_type)
            gro.minimise_merged(residue_type, np_system)
        system.update(np_system)
        time_counter['b_n_p_t']=time.time()

    time_counter['n_p_t']=time.time()

    #### creates merged folder
    print('\nMerging all residue types to single file. (Or possibly tea)\n')

    if len(system)>0:
    #### make final topology in merged directory
        gro.write_merged_topol(system, '_novo')
    #### make minimisation directory
        gro.make_min('merged_cg2at')
    #### merges provided atomistic protein and residues types into a single pdb file into merged directory
        if user_at_input and 'PROTEIN' in system:
            at_mod.merge_system_pdbs(system, '_no_steered', cg_residues, box_vec)
            at_mod.merge_system_pdbs(system, '_at_rep_user_supplied', cg_residues, box_vec)
            gro.minimise_merged_pdbs(system, '_at_rep_user_supplied')
            if len(system) > 1 and g_var.alchembed:
                gro.alchembed(system['PROTEIN'])
            else:
                gen.file_copy_and_check(g_var.working_dir+'MERGED/min/merged_cg2at_at_rep_user_supplied_minimised.pdb', g_var.final_dir+'final_cg2at_at_rep_user_supplied.pdb')
                gen.file_copy_and_check(g_var.working_dir+'MERGED/merged_cg2at_no_steered.pdb', g_var.final_dir+'final_cg2at_no_steered.pdb')
    #### merges de novo protein and residues types into a single pdb file into merged directory
        at_mod.merge_system_pdbs(system, '_novo', cg_residues, box_vec)
        gro.minimise_merged_pdbs(system, '_novo')
        gen.file_copy_and_check('merged_cg2at_novo_minimised.pdb', g_var.final_dir+'final_cg2at_de_novo.pdb')
        time_counter['m_t']=time.time()

    #### copies all itp files and topologies from whereever they are stored
        for file_name in os.listdir(g_var.working_dir+'MERGED'):
            if file_name.endswith('.itp') or file_name.endswith('final.top'):
                gen.file_copy_and_check(g_var.working_dir+'MERGED/'+file_name, g_var.final_dir+file_name)


    if 'PROTEIN' in cg_residues:
    #### creates mdp file if user wants to pull the structure to initial input
        if not os.path.exists(g_var.final_dir+'steered_md.mdp'):
            with open(g_var.final_dir+'steered_md.mdp', 'w') as steered_md:
                steered_md.write('define = -DPOSRES\nintegrator = md\nnsteps = 3000\ndt = 0.001\ncontinuation   = no\nconstraint_algorithm = lincs\n')
                steered_md.write('constraints = h-bonds\nns_type = grid\nnstlist = 25\nrlist = 1\nrcoulomb = 1\nrvdw = 1\ncoulombtype  = PME\n')
                steered_md.write('pme_order = 4\nfourierspacing = 0.16\ntcoupl = V-rescale\ntc-grps = system\ntau_t = 0.1\nref_t = 310\npcoupl = no\n')
                steered_md.write('pbc = xyz\nDispCorr = no\ngen_vel = yes\ngen_temp = 310\ngen_seed = -1')    
    #### calculates final RMS
        RMSD={}
        de_novo_atoms = at_mod_p.read_in_atomistic(g_var.final_dir+'final_cg2at_de_novo.pdb', system['PROTEIN'], sequence, False)
        RMSD['de novo '] = at_mod_p.RMSD_measure(de_novo_atoms, system,backbone_coords)

        if user_at_input and 'PROTEIN' in system:
            at_input_atoms = at_mod_p.read_in_atomistic(g_var.final_dir+'final_cg2at_at_rep_user_supplied.pdb', system['PROTEIN'], sequence, False)
            RMSD['at input'] = at_mod_p.RMSD_measure(at_input_atoms, system,backbone_coords)         
        print('\n{0:^10}{1:^25}{2:^10}'.format('output ','chain','RMSD ('+chr(197)+')'))
        print('{0:^10}{1:^25}{2:^10}'.format('-------','-----','---------'))
        for rmsd in RMSD:
            for chain in RMSD[rmsd]:
                print('{0:^10}{1:^25}{2:^10}'.format(rmsd, str(chain), float(RMSD[rmsd][chain])))

    #### removes temp file from script, anything with temp in really
    if g_var.clean:
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



