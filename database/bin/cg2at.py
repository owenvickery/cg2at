#!/usr/bin/env python3

import os, sys
import numpy as np
import time
import multiprocessing as mp
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/database/bin')
import gen, gro, at_mod, at_mod_p, at_mod_np, read_in, g_var, check_library


if __name__ == '__main__':
    mp.freeze_support()
    ## hardcoded varibles used by the script
    ## I've tried to make them as comprehensive as possible but they may need updating occasionally
    g_var.version = 0.9

    g_var.script_update = '08-07-2021'

    g_var.other = {'DA':'A', 'DG':'G', 'DC':'C', 'DT':'T', 'DAX':'A', 'DGX':'G', 'DCX':'C', 'DTX':'T'}

    g_var.termini_selections = {'charmm':{'N_TERMINAL':{'PRO':{'NH2+':0,'NH':1,'NH3+':2, '5TER':3, 'NONE':4},
                                                        'NORM':{'NH3+':0,'NH2':1,'5TER':2, 'NONE':3},}, 
                                          'C_TERMINAL':{'NORM':{'COO-':0,'COOH':1,'CT2':2, '3TER':3, 'NONE':4},
                                                        'PRO':{'COO-':0,'COOH':1,'CT2':2, '3TER':3, 'NONE':4}}},
                                'opls':{'N_TERMINAL':{'PRO':{'NH':2, 'NH3+':3, 'NONE':5},
                                                      'NORM':{'NH3+':0,'NH2':2, 'NONE':3}}, 
                                        'C_TERMINAL':{'NORM':{'COO-':0,'COOH':2, 'NONE':3},
                                                     'PRO':{'COO-':0,'COOH':2, 'NONE':3}}}
                                }

# Nothing in the script below here should need changing by the user
    
    g_var.tc['i_t']=time.time()
    ### initialise script 
    gen.cg2at_header()
    gen.fetch_forcefield_water_info()
    gen.check_input_flag() #### if missing structure file print help and quit
    gen.correct_number_cpus()
    gen.find_gromacs()
    gen.read_database_directories()
    gen.forcefield_selection()
    gen.fragment_selection()
    gen.check_water_molecules()
    if g_var.args.posre != None and len(g_var.np_directories) > 0:
        check_library.add_posres_file()
    if g_var.args.compare != None and len(g_var.np_directories) > 0:
        check_library.compare_forcefield_to_database()
    if g_var.args.info:
        gen.database_information()
    if g_var.args.v >= 1:
        print(gen.fragments_in_use())

    gen.fetch_fragment_multi() 
    gen.fetch_fragment_single()   
    gen.fetch_chain_groups()
    gen.sort_swap_group()
    print(gen.print_swap_residues())
    ###
    #### collects initial structures into INPUT folder
    gro.collect_input()
    #### saves flags used into INPUT folder
    gen.flags_used()
    g_var.tc['i_t_e']=time.time()
    #### reads in CG file and separates into residue types
    box_vec_initial = read_in.read_initial_cg_pdb()
    #### box size update 
    if g_var.args.box != None:
        print('box cutting only works for cubic boxes currently')
        g_var.box_vec, box_shift = gen.new_box_vec(box_vec_initial, g_var.args.box)
    else:
        g_var.box_vec=box_vec_initial
        box_shift=np.array([0,0,0])
    read_in.real_box_vectors(g_var.box_vec)
    #### pbc fix and residue truncation if required
    read_in.fix_pbc(box_vec_initial, g_var.box_vec, box_shift)
    #### checks if fragment database and input files match  
    at_mod.sanity_check()
    ### convert protein to atomistic representation
    g_var.tc['r_i_t']=time.time()
    if 'PROTEIN' in g_var.cg_residues:          
        g_var.coord_atomistic = at_mod_p.build_multi_residue_atomistic_system(g_var.cg_residues, 'PROTEIN') ## converts protein to atomistic
        if not g_var.user_at_input and g_var.args.v >= 1:  ## prints protein sequences 
            print(gen.print_sequnce_info('PROTEIN'))
        ## reads in user chain, runs a sequence alignment and finds existing disulphide bonds
        g_var.tc['p_d_n_t']=time.time()
        if g_var.user_at_input:
            for file_num, file_name in enumerate(g_var.args.a):
                atomistic_protein_input_raw, g_var.chain_count = read_in.read_in_atomistic(g_var.input_directory+'AT_INPUT_'+str(file_num)+'.pdb')  ## reads in user structure
                g_var.atomistic_protein_input_raw.update(atomistic_protein_input_raw)
            read_in.duplicate_chain()  ## duplicates user chcains
            at_mod_p.check_sequence()  ## gets user sequence
            at_mod_p.align_chain_sequence('PROTEIN') ## aligns chains 
            at_mod_p.find_disulphide_bonds_user_sup() ## finds user disulphide bonds
        at_mod_p.find_disulphide_bonds_de_novo() ## finds CG disulphide bonds 
        g_var.coord_atomistic = at_mod_p.correct_disulphide_bonds(g_var.coord_atomistic) ## fixes sulphur distances
        final_coordinates_atomistic_de_novo = at_mod_p.finalise_novo_atomistic(g_var.coord_atomistic, 'PROTEIN') ## fixes carbonyl oxygens, hydrogens and writes pdb 
        ## aligns user chains to the CG system
        if g_var.user_at_input:
            at_mod_p.align_user_chains(final_coordinates_atomistic_de_novo)
        #### read in minimised de novo protein chains and merges chains
        if not os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN_de_novo_merged.pdb'):
            gro.run_parallel_pdb2gmx_min('PROTEIN', g_var.ter_res['PROTEIN'])### runs pdb2gmx and minimises each protein chain
            print('Merging de_novo protein chains')
            at_mod.merge_indivdual_chain_pdbs(g_var.working_dir+'PROTEIN/MIN/PROTEIN_de_novo', '.pdb', 'PROTEIN') ## merge protein chains

        #### read in aligned protein chains and merges chains
        if g_var.user_at_input and not os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN_aligned_merged.pdb'):
            print('Merging aligned protein chains')
            if g_var.args.o not in ['none', 'align']:
                at_mod.merge_indivdual_chain_pdbs(g_var.working_dir+'PROTEIN/MIN/PROTEIN_aligned', '.pdb', 'PROTEIN') ## merge aligned chains
            else:
                at_mod.merge_indivdual_chain_pdbs(g_var.working_dir+'PROTEIN/PROTEIN_aligned', '_gmx_checked.pdb', 'PROTEIN')

    ### converts other linked residues  
    g_var.tc['f_p_t']=time.time()
    if 'OTHER' in g_var.cg_residues:  
        g_var.other_atomistic = at_mod_p.build_multi_residue_atomistic_system(g_var.cg_residues, 'OTHER')   
        if g_var.args.v >= 1:  ## prints protein sequences 
            print(gen.print_sequnce_info('OTHER'))   
        fin_at_NP_linked_de_novo = at_mod_p.finalise_novo_atomistic(g_var.other_atomistic, 'OTHER')
        gro.run_parallel_pdb2gmx_min('OTHER', g_var.ter_res['OTHER'])
        if not os.path.exists(g_var.working_dir+'OTHER/OTHER_de_novo_merged.pdb'):
            at_mod.merge_indivdual_chain_pdbs(g_var.working_dir+'OTHER/MIN/OTHER_de_novo', '.pdb', 'OTHER') ## merge  chains
    g_var.tc['f_o_t']=time.time()

    #### converts non protein residues into atomistic (runs on all cores)
    if len([key for value, key in enumerate(g_var.cg_residues) if key not in ['PROTEIN', 'OTHER']]) > 0:
        print('\nConverting the following residues: \n')
        # os.chdir(g_var.start_dir)
        for residue_type in g_var.cg_residues.keys():
            if residue_type not in ['PROTEIN', 'OTHER']:
                number = at_mod_np.build_atomistic_system(residue_type) 
                g_var.system.update(number)
        # with mp.Pool(g_var.args.ncpus) as pool:
        #     pool_process = pool.starmap_async(at_mod_np.build_atomistic_system, [(residue_type, g_var.working_dir, g_var.cg_residues) 
        #                                 for residue_type in [key for key in g_var.cg_residues if key not in ['PROTEIN', 'OTHER']]]).get() ## fragment fitting done in parrallel  
        # for residue_type in pool_process:
        #     g_var.system.update(residue_type) ## updates residue counts 
        #### attempts to minimise all residues at once 
        print('\nThis may take some time....(probably time for a coffee)\n')
        for residue_type in [key for key in g_var.cg_residues if key not in ['PROTEIN', 'OTHER']]:
            if not os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_merged.pdb'):
                print('Minimising: '+residue_type) 
                error = gro.minimise_merged(residue_type, g_var.working_dir+residue_type+'/'+residue_type+'_all.pdb')
                if error == True and residue_type not in ['SOL']:
                    print('Failed to minimise as a group: '+residue_type)
                    print('please check your input file as there is likely something wrong')
 
    ### MERGES system
    g_var.tc['n_p_t']=time.time()

    print('Merging all residue types to single file. (Or possibly tea)\n')

    gro.write_merged_topol() ## make final topology in merged directory

    #### copies all itp files and topologies from wherever they are stored into the FINAL folder
    for file_name in os.listdir(g_var.merged_directory):
        if not any(f in file_name for f in ['steered_posre.itp', 'low_posre.itp','mid_posre.itp', 'high_posre.itp']):
            if file_name.endswith('.itp') or file_name.endswith('final.top') :
               gen.file_copy_and_check(g_var.merged_directory+file_name, g_var.final_dir+file_name)

    #### merges provided atomistic protein and residues types into a single pdb file into merged directory
    if not os.path.exists(g_var.merged_directory+'merged_cg2at_de_novo.pdb'):
        at_mod.merge_system_pdbs('_de_novo') ## merge all minimised residues into a complete system 

    ## minimise merged system
    if not os.path.exists(g_var.merged_directory+'MIN/merged_cg2at_de_novo_minimised.pdb'):
        gro.make_min('merged_cg2at') 
        gro.minimise_merged_pdbs( '_de_novo') ## minimise system pdb
    g_var.tc['m_t']=time.time()
    ## checks for threaded lipids, e.g. abnormal bonds lengths (not had a issue for a long time might delete) 
    if not os.path.exists(g_var.merged_directory+'checked_ringed_lipid_de_novo.pdb'):
        at_mod.check_ringed_lipids(g_var.merged_directory+'MIN/merged_cg2at_de_novo_minimised.pdb') ## check for abnormal bond lengths 
    ## runs short NVT on de_novo system with disres on if available 
    if g_var.args.o not in ['none', 'align']:
        gro.run_nvt(g_var.merged_directory+'checked_ringed_lipid_de_novo') ## run npt on system 
    else:
        print('Completed initial minimisation, please find final de_novo system: \n'+g_var.final_dir+'final_cg2at_de_novo.pdb')
        gen.file_copy_and_check(g_var.merged_directory+'checked_ringed_lipid_de_novo.pdb', g_var.final_dir+'final_cg2at_de_novo.pdb')
    g_var.tc['eq_t']=time.time()
    ## creates aligned system 

    if g_var.user_at_input and 'PROTEIN' in g_var.cg_residues and g_var.args.o in ['all', 'align']:   
        g_var.tc['a_s']=time.time()
        gro.create_aligned()
        g_var.tc['a_e']=time.time()

    ## removes temp file from script
    if not g_var.args.messy:
        gen.clean() 

    ## prints out system information
    print(gen.write_system_components())

    ## prints out RMSD of converted proteins
    if 'PROTEIN' in g_var.cg_residues:
        at_mod_p.write_RMSD()

    g_var.tc['f_t']=time.time()

    ## prints out script timings for each section
    gen.print_script_timings()



