#!/usr/bin/env python3

import os, sys
import numpy as np
from shutil import copyfile
import time
import re
import multiprocessing as mp
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/database/script_files')
import gen, gro, at_mod, at_mod_p, at_mod_np, read_in, g_var, f_loc

time_counter = {}
time_counter['i_t']=time.time()

#### collects initial structures into INPUT folder
user_at_input = gro.collect_input(g_var.c, g_var.a)

#### saves flags used into INPUT folder
gen.flags_used()

print('\nThis script is now hopefully doing the following (Good luck):\n\nReading in your CG representation\n')

#### reads in CG file and separates into residue types
cg_residues, box_vec_initial = read_in.read_initial_cg_pdb()

#### box size update 
if g_var.box != None:
    print('box cutting only works for cubic boxes currently')
    g_var.box_vec, box_shift = gen.new_box_vec(box_vec_initial, g_var.box)
else:
    g_var.box_vec=box_vec_initial
    box_shift=np.array([0,0,0])

#### pbc fix and residue truncation if required
cg_residues = read_in.fix_pbc(cg_residues, box_vec_initial, g_var.box_vec, box_shift)

#### checks if fragment database is correct
at_mod.sanity_check(cg_residues)

### convert protein to atomistic representation
time_counter['r_i_t']=time.time()
system={}    
if 'PROTEIN' in cg_residues:          
    p_system, backbone_coords, coordinates_atomistic, sequence=at_mod_p.build_protein_atomistic_system(cg_residues['PROTEIN']) ## converts protein to atomistic
    system['PROTEIN']=p_system['PROTEIN']

    if not user_at_input and g_var.v >= 1:  ## prints protein sequences 
        print('coarse grain protein sequence:\n')
        for index in sequence:
            print('chain:', index,sequence[index], '\n') 

    ## reads in user chain, runs a sequence alignment and finds existing disulphide bonds
    time_counter['p_d_n_t']=time.time()
    if user_at_input:
        atomistic_protein_input_raw, chain_count = read_in.read_in_atomistic(g_var.input_directory+'AT_input.pdb', True)  ## reads in user structure
        seq_user = at_mod_p.check_sequence(atomistic_protein_input_raw, chain_count)  ## gets user sequence
        atomistic_protein_input, group_chain, user_at_input = at_mod_p.align_chains(atomistic_protein_input_raw, seq_user, sequence) ## aligns chains 
        user_cys_bond = at_mod_p.find_disulphide_bonds_user_sup(atomistic_protein_input) ## finds user disulphide bonds
    else:
        user_cys_bond = {}

    user_cys_bond = at_mod_p.find_disulphide_bonds_de_novo(coordinates_atomistic, user_cys_bond) ## finds CG disulphide bonds 
    coordinates_atomistic = at_mod_p.correct_disulphide_bonds(coordinates_atomistic, user_cys_bond) ## fixes sulphur distances
    final_coordinates_atomistic = at_mod_p.finalise_novo_atomistic(coordinates_atomistic, cg_residues['PROTEIN']) ## fixes carbonyl oxygens, hydrogens and writes pdb 

    ## aligns user chains to the CG system
    if user_at_input:
        at_mod_p.align_user_chains(atomistic_protein_input, backbone_coords, group_chain, final_coordinates_atomistic, user_cys_bond)

    ### runs pdb2gmx and minimises each protein chain
    pool = mp.Pool(g_var.ncpus)
    m = mp.Manager()
    q = m.Queue()
    pool_process = pool.starmap_async(gro.pdb2gmx_minimise, [(chain, p_system, user_at_input, q) for chain in range(0, system['PROTEIN'])])
    while not pool_process.ready():
        gro.report_complete('pdb2gmx/minimisation', q.qsize(), system['PROTEIN'])
    print('{:<130}'.format(''), end='\r')
    print('pdb2gmx/minimisation completed on residue type: PROTEIN')     
    pool.close()
    #### read in minimised de novo protein chains and merges chains
    if not os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN_de_novo_merged.pdb'):
        merge_de_novo = at_mod_p.read_in_protein_pdbs(system['PROTEIN'], g_var.working_dir+'PROTEIN/MIN/PROTEIN_de_novo', '.pdb') ## merge protein chains
        at_mod_p.write_merged_pdb(merge_de_novo, '_de_novo') ## write merged chain to pdb 
    #### read in aligned protein chains and merges chains
    if user_at_input:
        if not os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN_aligned_merged.pdb'):
            merge_at_user_no_steer = at_mod_p.read_in_protein_pdbs(system['PROTEIN'], g_var.working_dir+'PROTEIN/PROTEIN_aligned', '_gmx_checked.pdb') ## merge aligned chains
            at_mod_p.write_merged_pdb(merge_at_user_no_steer, '_aligned') ## write merged chain to pdb 

time_counter['f_p_t']=time.time()

#### converts non protein residues into atomistic (runs on all cores)
if len([key for value, key in enumerate(cg_residues) if key not in ['PROTEIN']]) > 0:
    print('\nConverting the following residues concurrently: ')
    np_system={}
    pool = mp.Pool(g_var.ncpus)
    pool_process = pool.starmap_async(at_mod_np.build_atomistic_system, [(cg_residues, residue_type) 
                                    for residue_type in [key for value, key in enumerate(cg_residues) if key not in ['PROTEIN']]]).get() ## fragment fitting done in parrallel  
    pool.close()
    for residue_type in pool_process:
        np_system.update(residue_type) ## updates residue counts
    #### attempts to minimise all residues at once else falls back to doing individually
    print('\nThis may take some time....(probably time for a coffee)\n')
    for residue_type in [key for value, key in enumerate(cg_residues) if key not in ['PROTEIN', 'ION']]:
        if not os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_merged.pdb'):
            print('Minimising merged: '+residue_type+'\n') 
            error = gro.minimise_merged(residue_type, np_system, g_var.working_dir+residue_type+'/'+residue_type+'_all.pdb')
            if error == True:
                print('Failed to minimise as a group now processing individual residues: '+residue_type)
                gro.non_protein_minimise_ind(np_system[residue_type], residue_type) ## runs grompp and minimises each residue
                at_mod_np.merge_minimised(residue_type, np_system) ## merges minimised residues
                print('Minimising merged: '+residue_type+'\n') 
                gro.minimise_merged(residue_type, np_system, g_var.working_dir+residue_type+'/MIN/'+residue_type+'_merged.pdb') ## minimises merged residues
    system.update(np_system)

time_counter['n_p_t']=time.time()

print('Merging all residue types to single file. (Or possibly tea)\n')

gro.write_merged_topol(system, '_de_novo') ## make final topology in merged directory

#### copies all itp files and topologies from wherever they are stored into the FINAL folder
for file_name in os.listdir(g_var.merged_directory):
    if not any(f in file_name for f in ['steered_posre.itp', 'low_posre.itp','mid_posre.itp', 'high_posre.itp']):
        if file_name.endswith('.itp') or file_name.endswith('final.top') :
           gen.file_copy_and_check(g_var.merged_directory+file_name, g_var.final_dir+file_name)
gro.make_min('merged_cg2at')
#### merges provided atomistic protein and residues types into a single pdb file into merged directory
if not os.path.exists(g_var.merged_directory+'merged_cg2at_de_novo.pdb'):
    at_mod.merge_system_pdbs(system, '_de_novo', cg_residues) ## merge all minimised residues into a complete system 
## minimise merged system
if not os.path.exists(g_var.merged_directory+'MIN/merged_cg2at_de_novo_minimised.pdb'):
    gro.minimise_merged_pdbs(system, '_de_novo') ## minimise system pdb
## checks for threaded lipids, e.g. abnormal bonds lengths (not had a issue for a long time might delete) 
if not os.path.exists(g_var.merged_directory+'checked_ringed_lipid_de_novo.pdb'):
    ringed_lipids = at_mod.check_ringed_lipids(g_var.merged_directory+'MIN/merged_cg2at_de_novo_minimised.pdb') ## check for abnormal bond lengths 
    if len(system) > 1 and g_var.alchembed and len(ringed_lipids) > 0 and 'PROTEIN' in cg_residues:
        gro.alchembed(system['PROTEIN'], 'de_novo') ## runs alchembed on protein chains 
        ringed_lipids = at_mod.check_ringed_lipids(g_var.merged_directory+'checked_ringed_lipid_de_novo.pdb') ## rechecks for abnormal bond lengths
        if len(ringed_lipids) > 0:
            print('Check final output as alchembed cannot fix ringed lipid: ', ringed_lipids) ## warning that the script failed to fix bonds
    else:
        gen.file_copy_and_check(g_var.merged_directory+'MIN/merged_cg2at_de_novo_minimised.pdb', g_var.merged_directory+'checked_ringed_lipid_de_novo.pdb')
## runs short NPT on de_novo system with disres on if available 
if g_var.o in ['all', 'de_novo']:
    gro.run_npt(g_var.merged_directory+'checked_ringed_lipid_de_novo', user_at_input) ## run npt on system 
else:
    gen.file_copy_and_check(g_var.merged_directory+'checked_ringed_lipid_de_novo.pdb', g_var.final_dir+'final_cg2at_de_novo.pdb')

## creates aligned system 
time_counter['m_t']=time.time()
if user_at_input and 'PROTEIN' in cg_residues and g_var.o in ['all', 'align']:   
    time_counter['a_s']=time.time()
    gro.create_aligned(system, cg_residues)
    time_counter['a_e']=time.time()

## removes temp file from script, anything with temp in really
if not g_var.messy:
    gen.clean(cg_residues) 

## prints out system information
gen.write_system_components(system)

## prints out RMSD of converted proteins
if 'PROTEIN' in cg_residues:
    at_mod.write_RMSD(system)

time_counter['f_t']=time.time()

## prints out script timings for each section
gen.print_script_timings(time_counter, system, user_at_input)



