#!/usr/bin/env python3

import os, sys
import numpy as np
from shutil import copyfile
import time
import re
import multiprocessing as mp
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/database/script_files')
import gen, gro, at_mod, at_mod_p, at_mod_np, read_in, g_var

# Nothing in the script should need changing by the user

g_var.tc['i_t']=time.time()
### initialise script 
gen.correct_number_cpus()
gen.find_gromacs()
gen.read_database()
gen.forcefield_selection()
gen.fragment_selection()
gen.check_water_molecules()
if g_var.info:
    database_information()
gen.fetch_fragment()    
gen.fetch_chain_groups()
gen.sort_swap_group()
###

#### collects initial structures into INPUT folder
gro.collect_input()

#### saves flags used into INPUT folder
gen.flags_used()

#### reads in CG file and separates into residue types
box_vec_initial = read_in.read_initial_cg_pdb()
#### box size update 
if g_var.box != None:
    g_var.box_vec, box_shift = gen.new_box_vec(box_vec_initial, g_var.box)
else:
    g_var.box_vec=box_vec_initial
    box_shift=np.array([0,0,0])
#### pbc fix and residue truncation if required
read_in.fix_pbc(box_vec_initial, g_var.box_vec, box_shift)
#### checks if fragment database and input files match  
at_mod.sanity_check()
### convert protein to atomistic representation
g_var.tc['r_i_t']=time.time()
if 'PROTEIN' in g_var.cg_residues:          
    g_var.coord_atomistic = at_mod_p.build_multi_residue_atomistic_system(g_var.cg_residues, 'PROTEIN') ## converts protein to atomistic
    if not g_var.user_at_input and g_var.v >= 1:  ## prints protein sequences 
        print('coarse grain protein sequence:\n')
        for index in g_var.seq_cg:
            print('chain:', index,g_var.seq_cg[index], '\n') 
    ## reads in user chain, runs a sequence alignment and finds existing disulphide bonds
    g_var.tc['p_d_n_t']=time.time()
    if g_var.user_at_input:
        for file_num, file_name in enumerate(g_var.a):
            atomistic_protein_input_raw, g_var.chain_count = read_in.read_in_atomistic(g_var.input_directory+'AT_INPUT_'+str(file_num)+'.pdb')  ## reads in user structure
            g_var.atomistic_protein_input_raw.update(atomistic_protein_input_raw)
        read_in.duplicate_chain()
        at_mod_p.check_sequence()  ## gets user sequence
        at_mod_p.align_chain_sequence('PROTEIN') ## aligns chains 
        at_mod_p.find_disulphide_bonds_user_sup() ## finds user disulphide bonds
    at_mod_p.find_disulphide_bonds_de_novo() ## finds CG disulphide bonds 
    g_var.coord_atomistic = at_mod_p.correct_disulphide_bonds(g_var.coord_atomistic) ## fixes sulphur distances
    final_coordinates_atomistic_de_novo = at_mod_p.finalise_novo_atomistic(g_var.coord_atomistic, 'PROTEIN') ## fixes carbonyl oxygens, hydrogens and writes pdb 

    ## aligns user chains to the CG system
    if g_var.user_at_input:
        at_mod_p.align_user_chains(final_coordinates_atomistic_de_novo)
    ### runs pdb2gmx and minimises each protein chain
    gro.run_parallel_pdb2gmx_min('PROTEIN', g_var.ter_res['PROTEIN'])

    #### read in minimised de novo protein chains and merges chains
    if not os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN_de_novo_merged.pdb'):
        at_mod.merge_indivdual_chain_pdbs(g_var.working_dir+'PROTEIN/MIN/PROTEIN_de_novo', '.pdb', 'PROTEIN') ## merge protein chains

    #### read in aligned protein chains and merges chains
    if g_var.user_at_input and not os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN_aligned_merged.pdb'):
            at_mod.merge_indivdual_chain_pdbs(g_var.working_dir+'PROTEIN/MIN/PROTEIN_aligned', '.pdb', 'PROTEIN') ## merge aligned chains

### converts other linked residues  
g_var.tc['f_p_t']=time.time()
if 'OTHER' in g_var.cg_residues:  
    g_var.other_atomistic = at_mod_p.build_multi_residue_atomistic_system(g_var.cg_residues, 'OTHER')      
    fin_at_NP_linked_de_novo = at_mod_p.finalise_novo_atomistic(g_var.other_atomistic, 'OTHER')
    gro.run_parallel_pdb2gmx_min('OTHER', g_var.ter_res['OTHER'])
    if not os.path.exists(g_var.working_dir+'OTHER/OTHER_de_novo_merged.pdb'):
        at_mod.merge_indivdual_chain_pdbs(g_var.working_dir+'OTHER/MIN/OTHER_de_novo', '.pdb', 'OTHER') ## merge  chains

#### converts non protein residues into atomistic (runs on all cores)
if len([key for value, key in enumerate(g_var.cg_residues) if key not in ['PROTEIN', 'OTHER']]) > 0:
    print('\nConverting the following residues concurrently: ')
    pool = mp.Pool(g_var.ncpus)
    pool_process = pool.starmap_async(at_mod_np.build_atomistic_system, [(residue_type, 1) 
                                    for residue_type in [key for key in g_var.cg_residues if key not in ['PROTEIN', 'OTHER']]]).get() ## fragment fitting done in parrallel  
    pool.close()
    for residue_type in pool_process:
        g_var.system.update(residue_type) ## updates residue counts 

    #### attempts to minimise all residues at once else falls back to doing individually
    print('\nThis may take some time....(probably time for a coffee)\n')
    for residue_type in [key for key in g_var.cg_residues if key not in ['PROTEIN', 'ION', 'OTHER']]:
        if not os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_merged.pdb'):
            print('Minimising merged: '+residue_type+'\n') 
            error = gro.minimise_merged(residue_type, g_var.working_dir+residue_type+'/'+residue_type+'_all.pdb')
            if error == True:
                print('Failed to minimise as a group now processing individual residues: '+residue_type)
                gro.non_protein_minimise_ind(residue_type) ## runs grompp and minimises each residue
                at_mod_np.merge_minimised(residue_type) ## merges minimised residues
                gro.minimise_merged(residue_type, g_var.working_dir+residue_type+'/MIN/'+residue_type+'_merged.pdb') ## minimises merged residues

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

## checks for threaded lipids, e.g. abnormal bonds lengths (not had a issue for a long time might delete) 
if not os.path.exists(g_var.merged_directory+'checked_ringed_lipid_de_novo.pdb'):
    ringed_lipids = at_mod.check_ringed_lipids(g_var.merged_directory+'MIN/merged_cg2at_de_novo_minimised.pdb') ## check for abnormal bond lengths 
    if g_var.alchembed and len(ringed_lipids) > 0 and 'PROTEIN' in g_var.cg_residues:
        gro.alchembed() ## runs alchembed on protein chains 
    else:
        gen.file_copy_and_check(g_var.merged_directory+'MIN/merged_cg2at_de_novo_minimised.pdb', g_var.merged_directory+'checked_ringed_lipid_de_novo.pdb')


## runs short NVT on de_novo system with disres on if available 
if g_var.o in ['all', 'de_novo']:
    gro.run_nvt(g_var.merged_directory+'checked_ringed_lipid_de_novo') ## run npt on system 
else:
    gen.file_copy_and_check(g_var.merged_directory+'checked_ringed_lipid_de_novo.pdb', g_var.final_dir+'final_cg2at_de_novo.pdb')

## creates aligned system 
g_var.tc['m_t']=time.time()
if g_var.user_at_input and 'PROTEIN' in g_var.cg_residues and g_var.o in ['all', 'align']:   
    g_var.tc['a_s']=time.time()
    gro.create_aligned()
    g_var.tc['a_e']=time.time()

## removes temp file from script
if not g_var.messy:
    gen.clean() 

## prints out system information
gen.write_system_components()

## prints out RMSD of converted proteins
if 'PROTEIN' in g_var.cg_residues:
    at_mod_p.write_RMSD()

g_var.tc['f_t']=time.time()

## prints out script timings for each section
gen.print_script_timings()



