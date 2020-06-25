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
    box_vec, box_shift = gen.new_box_vec(box_vec_initial, g_var.box)
else:
    box_vec=box_vec_initial
    box_shift=np.array([0,0,0])

#### simple pbc fix and residue truncation if required
cg_residues = read_in.fix_pbc(cg_residues, box_vec_initial, box_vec, box_shift)
if 'PROTEIN' in cg_residues:
    protein = True
else:
    protein = False

#### checks if fragment database is correct
at_mod.sanity_check(cg_residues)

### convert protein to atomistic representation
time_counter['r_i_t']=time.time()
system={}    
if protein:          
    p_system, backbone_coords, coordinates_atomistic, sequence=at_mod_p.build_protein_atomistic_system(cg_residues['PROTEIN']) ## converts protein to atomistic

    if not user_at_input and g_var.v >= 1:  ## prints protein sequences 
        print('coarse grain protein sequence:\n')
        for index in sequence:
            print('chain:', index,sequence[index], '\n') 

    system['PROTEIN']=p_system['PROTEIN']
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
    final_coordinates_atomistic = at_mod_p.finalise_novo_atomistic(coordinates_atomistic, cg_residues['PROTEIN'], box_vec) ## fixes carbonyl oxygens, hydrogens and writes pdb 
    if user_at_input:
        atomistic_protein_centered, cg_com = at_mod_p.center_atomistic(atomistic_protein_input, backbone_coords, group_chain) ## centers each chain by center of mass
        at_com_group, cg_com_group = at_mod_p.rotate_protein_monomers(atomistic_protein_centered, final_coordinates_atomistic, 
                                                                    backbone_coords, cg_com, box_vec, group_chain) ## finds rotation matrix for rigid fit
        final_user_supplied_coord = at_mod_p.apply_rotations_to_chains(final_coordinates_atomistic, atomistic_protein_centered, 
                                                                    at_com_group,cg_com_group,cg_com, box_vec, group_chain) ## apply rotation matrix to atoms and build in missing residues
        final_user_supplied_coord = at_mod_p.correct_disulphide_bonds(final_user_supplied_coord, user_cys_bond) ## fixes sulphur distances in user structure
        pool = mp.Pool(mp.cpu_count())
        pool_process = pool.starmap_async(at_mod_p.write_user_chains_to_pdb, [(final_user_supplied_coord[chain], box_vec, chain) ## write structure to pdb
                                        for chain in final_user_supplied_coord]).get()
        pool.close()
    #### minimise each protein chain
    print('Minimising '+str(p_system['PROTEIN'])+' protein chains')
    gro.minimise_protein(system['PROTEIN'], p_system, user_at_input, box_vec) ## minimise user and de novo structures
    #### read in minimised de novo protein chains and merges chains
    if not os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN_de_novo_merged.pdb'):
        merge_de_novo = at_mod_p.read_in_protein_pdbs(system['PROTEIN'], g_var.working_dir+'PROTEIN/MIN/PROTEIN_de_novo', '.pdb') ## merge protein chains
        at_mod_p.write_merged_pdb(merge_de_novo, '_de_novo', box_vec) ## write merged chain to pdb 
    #### runs steered MD on user supplied protein chains
    if user_at_input:
        if g_var.o in ['all', 'align'] and not os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN_aligned_merged.pdb'):
            merge_at_user_no_steer = at_mod_p.read_in_protein_pdbs(system['PROTEIN'], g_var.working_dir+'PROTEIN/PROTEIN_aligned', '_gmx.pdb') ## merge aligned chains
            at_mod_p.write_merged_pdb(merge_at_user_no_steer, '_aligned', box_vec) ## write merged chain to pdb 

time_counter['f_p_t']=time.time()

#### converts non protein residues into atomistic (runs on all cores)
if len([key for value, key in enumerate(cg_residues) if key not in ['PROTEIN']]) > 0:
    print('\nConverting the following residues concurrently: ')
    np_system={}
    pool = mp.Pool(mp.cpu_count())
    pool_process = pool.starmap_async(at_mod_np.build_atomistic_system, [(cg_residues, residue_type, box_vec) 
                                    for residue_type in [key for value, key in enumerate(cg_residues) if key not in ['PROTEIN']]]).get() ## fragment fitting done in parrallel  
    pool.close()
    for residue_type in pool_process:
        np_system.update(residue_type) ## updates residue counts
    #### minimises each residue separately
    print('\nThis may take some time....(probably time for a coffee)\n')
    for residue_type in [key for value, key in enumerate(cg_residues) if key not in ['PROTEIN', 'ION']]:
        print('Processing individual residues: '+residue_type)
        gro.non_protein_minimise(np_system[residue_type], residue_type) ## runs grompp and minimises each residue
        at_mod_np.merge_minimised(residue_type, np_system, box_vec) ## merges minimised residues
        print('Minimising merged: '+residue_type+'\n') 
        gro.minimise_merged(residue_type, np_system) ## minimises merged residues
    system.update(np_system)

time_counter['n_p_t']=time.time()

print('\nMerging all residue types to single file. (Or possibly tea)\n')

gro.write_merged_topol(system, '_de_novo') ## make final topology in merged directory
#### copies all itp files and topologies from wherever they are stored
for file_name in os.listdir(g_var.merged_directory):
    if not any(f in file_name for f in ['steered_posre.itp', 'low_posre.itp', 'high_posre.itp']):
        if file_name.endswith('.itp') or file_name.endswith('final.top') :
           gen.file_copy_and_check(g_var.merged_directory+file_name, g_var.final_dir+file_name)
gro.make_min('merged_cg2at') ## make minimisation directory
#### merges provided atomistic protein and residues types into a single pdb file into merged directory
if not os.path.exists(g_var.merged_directory+'merged_cg2at_de_novo.pdb'):
    at_mod.merge_system_pdbs(system, '_de_novo', cg_residues, box_vec) ## merge all minimised residues into a complete system 
if not os.path.exists(g_var.merged_directory+'MIN/merged_cg2at_de_novo_minimised.pdb'):
    gro.minimise_merged_pdbs(system, '_de_novo') ## minimise system pdb
if not os.path.exists(g_var.merged_directory+'NVT/merged_cg2at_de_novo_nvt.pdb'):
    gro.run_nvt(g_var.merged_directory+'MIN/merged_cg2at_de_novo_minimised.pdb', protein) ## run nvt on system
if not os.path.exists(g_var.final_dir+'final_cg2at_de_novo.pdb'):
    ringed_lipids = at_mod.check_ringed_lipids(g_var.merged_directory+'NVT/merged_cg2at_de_novo_nvt.pdb', box_vec) ## check for abnormal bond lengths 
    if len(system) > 1 and g_var.alchembed and len(ringed_lipids) > 0 and protein:
        gro.alchembed(system['PROTEIN'], 'de_novo') ## runs alchembed on protein chains 
        ringed_lipids = at_mod.check_ringed_lipids(g_var.final_dir+'final_cg2at_de_novo.pdb', box_vec) ## rechecks for abnormal bond lengths
        if len(ringed_lipids) > 0:
            print('Check final output as alchembed cannot fix ringed lipid: ', ringed_lipids) ## warning that the script failed to fix bonds
    else:
        gro.run_npt(g_var.merged_directory+'NVT/merged_cg2at_de_novo_nvt', protein) ## run npt on system if no abnormal bonds found
time_counter['m_t']=time.time()
if user_at_input and protein:

    if g_var.o in ['all', 'align']:   
        print('\nCreating aligned system') 
        time_counter['a_s']=time.time()
        at_mod.merge_system_pdbs(system, '_aligned', cg_residues, box_vec) ## create restraint positions for aligned system
        gro.steer_to_aligned('aligned', 'low', g_var.final_dir+'final_cg2at_de_novo') ## run steered md with low restraints
        gro.steer_to_aligned('aligned', 'mid', g_var.merged_directory+'STEER/merged_cg2at_aligned_steer_low') ## run steered md with medium restraints
        gro.steer_to_aligned('aligned', 'high', g_var.merged_directory+'STEER/merged_cg2at_aligned_steer_mid') ## run steered md with high restraints
        gen.file_copy_and_check(g_var.merged_directory+'STEER/merged_cg2at_aligned_steer_high.pdb', g_var.final_dir+'final_cg2at_aligned.pdb') ## copy to final folder
        time_counter['a_e']=time.time()
    if g_var.o in ['all', 'steer']:
        print('\nCreating steered system')
        time_counter['s_s']=time.time()
        # at_mod.merge_system_pdbs(system, '_steered', cg_residues, box_vec) ## create restraint positions for steered system
        gro.steer_to_de_novo(g_var.final_dir+'final_cg2at_de_novo', g_var.final_dir+'final_cg2at_aligned') ## run steered md with low restraints
        # gro.steer('steered', 'mid', g_var.merged_directory+'STEER/merged_cg2at_steered_steer_low') ## run steered md with medium restraints
        # gro.steer('steered', 'high', g_var.merged_directory+'STEER/merged_cg2at_steered_steer_mid') ## run steered md with high restraints
        gen.file_copy_and_check(g_var.merged_directory+'STEER/merged_cg2at_steered.pdb', g_var.final_dir+'final_cg2at_steered.pdb') ## copy to final folder
        time_counter['s_e']=time.time()

#### removes temp file from script, anything with temp in really
if not g_var.messy:
    gen.clean(cg_residues) ## removes files labelled with temp

#### prints out system information
print('\n{:-<100}'.format(''))
print('{0:^100}'.format('Script has completed, time for a beer'))
print('\n{0:^10}{1:^25}'.format('molecules','number'))
print('{0:^10}{1:^25}'.format('---------','------'))
for section in system:
    print('{0:^10}{1:^25}'.format(section, system[section]))

if protein:
#### calculates final RMS
    RMSD={}
    de_novo_atoms, chain_count = read_in.read_in_atomistic(g_var.final_dir+'final_cg2at_de_novo.pdb', False) ## reads in final pdb
    if chain_count != system['PROTEIN']:
        sys.exit('number of chains in atomistic protein input ('+str(chain_count)+') does not match CG representation ('+str(system['PROTEIN'])+')')
    RMSD['de novo '] = at_mod_p.RMSD_measure(de_novo_atoms, system,backbone_coords) ## gets rmsd of de novo

    if user_at_input and protein:
        if g_var.o in ['all', 'steer']:
            at_input_atoms, chain_count = read_in.read_in_atomistic(g_var.final_dir+'final_cg2at_steered.pdb', False)
            RMSD['at steered'] = at_mod_p.RMSD_measure(at_input_atoms, system,backbone_coords)   
        if g_var.o in ['all', 'align']: 
            at_input_atoms, chain_count = read_in.read_in_atomistic(g_var.final_dir+'final_cg2at_aligned.pdb', False)
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
        print()



time_counter['f_t']=time.time()

#### prints out script timings for each section

gen.print_script_timings(time_counter, system, user_at_input)



