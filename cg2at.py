#!/usr/bin/env python3

import os, sys
import numpy as np
from shutil import copyfile
from time import gmtime
import datetime
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/database/script_files')
import gen, gro, at_mod, at_mod_p, at_mod_np, cg_mod, g_var, f_loc


initialisation_time=np.array(gmtime()[3:6])

user_at_input = gro.collect_input(g_var.c, g_var.a)

print('\nThis script is now hopefully doing the following (Good luck):\n')

#### read in CG file
print('Reading in your CG representation\n')
cg_residues, box_vec_initial = cg_mod.read_initial_pdb()

if g_var.box != None:
    box_vec = gen.new_box_vec(box_vec_initial, g_var.box)
else:
    box_vec=box_vec_initial

cg_residues =cg_mod.fix_pbc(cg_residues, box_vec_initial, box_vec)
at_mod.sanity_check(cg_residues)

read_in_time=np.array(gmtime()[3:6])
system={}
### convert protein to atomistic
if 'PROTEIN' in cg_residues:
    p_system, backbone_coords, final_coordinates_atomistic, sequence=at_mod_p.build_protein_atomistic_system(cg_residues['PROTEIN'], box_vec)
    system['PROTEIN']=p_system['PROTEIN']
    protein_de_novo_time=np.array(gmtime()[3:6])
    if user_at_input and 'PROTEIN' in system:
    #### reads in atomistic structure   
        atomistic_protein_input = at_mod_p.read_in_atomistic(g_var.input_directory+'AT_input.pdb', system['PROTEIN'], sequence, True)  
        atomistic_protein_centered, cg_com = at_mod_p.center_atomistic(atomistic_protein_input, backbone_coords)
        at_mod_p.rotate_protein_monomers(atomistic_protein_centered, final_coordinates_atomistic, backbone_coords, cg_com, box_vec)
    gro.minimise_protein(system['PROTEIN'], p_system, user_at_input)

    merge_de_novo = at_mod_p.read_in_protein_pdbs(system['PROTEIN'], g_var.working_dir+'PROTEIN/min/PROTEIN_novo', '.pdb')
    at_mod_p.write_merged_pdb(merge_de_novo, '_novo', box_vec)

    if user_at_input and 'PROTEIN' in system:
        print('\tRunning steered MD on input atomistic structure\n')
    #### runs steered MD on atomistic structure on CA and CB atoms
        for chain in range(system['PROTEIN']):
            gro.steered_md_atomistic_to_cg_coord(chain)
        merge_at_user = at_mod_p.read_in_protein_pdbs(system['PROTEIN'], g_var.working_dir+'PROTEIN/steered_md/PROTEIN_at_rep_user_supplied', '.pdb')
        at_mod_p.write_merged_pdb(merge_at_user, '_at_rep_user_supplied', box_vec)
        merge_at_user_no_steer = at_mod_p.read_in_protein_pdbs(system['PROTEIN'], g_var.working_dir+'PROTEIN/PROTEIN_at_rep_user_supplied', '_gmx.pdb')
        at_mod_p.write_merged_pdb(merge_at_user_no_steer, '_no_steered', box_vec)


final_protein_time=np.array(gmtime()[3:6])

#### converts non protein residues into atomistic and minimises 
if len([key for value, key in enumerate(cg_residues) if key not in ['PROTEIN']]) > 0:
    np_system=at_mod_np.build_atomistic_system(cg_residues, box_vec)
    print('\nThis may take some time....(probably time for a coffee)\n')
    for residue_type in cg_residues:
        if residue_type not in ['PROTEIN', 'ION']:
            print('Minimising individual residues: '+residue_type)
            gro.non_protein_minimise(np_system[residue_type], residue_type)
            at_mod_np.merge_minimised(residue_type, np_system, box_vec)
            print('Minimising merged: '+residue_type)
            gro.minimise_merged(residue_type, np_system)
    system.update(np_system)
    build_non_protein_time=np.array(gmtime()[3:6])

non_protein_time=np.array(gmtime()[3:6])

#### creates merged folder
print('\nMerging all residue types to single file. (Or a possibly tea)\n')

if len(system)>0:
    gen.mkdir_directory(g_var.working_dir+'MERGED')
#### make final topology in merged directory
    gro.write_merged_topol(system, '_novo')
#### make minimisation directory
    gro.make_min('merged_cg2at')
#### merges provided atomistic protein and residues types into a single pdb file into merged directory
    if user_at_input and 'PROTEIN' in system:
        at_mod.merge_system_pdbs(system, '_no_steered', cg_residues, box_vec)
        at_mod.merge_system_pdbs(system, '_at_rep_user_supplied', cg_residues, box_vec)
        gro.minimise_merged_pdbs(system, '_at_rep_user_supplied')
        if len(system) > 1:
            gro.alchembed(system['PROTEIN'])
        else:
            copyfile(g_var.working_dir+'MERGED/min/merged_cg2at_at_rep_user_supplied_minimised.pdb', g_var.final_dir+'final_cg2at_at_rep_user_supplied.pdb')
            copyfile(g_var.working_dir+'MERGED/merged_cg2at_no_steered.pdb', g_var.final_dir+'final_cg2at_no_steered.pdb')
#### merges de novo protein and residues types into a single pdb file into merged directory
    at_mod.merge_system_pdbs(system, '_novo', cg_residues, box_vec)
    gro.minimise_merged_pdbs(system, '_novo')
    copyfile('merged_cg2at_novo_minimised.pdb', g_var.final_dir+'final_cg2at_de_novo.pdb')
    merge_time=np.array(gmtime()[3:6])

#### copies all itp files and topologies from whereever they are stored
    for file_name in os.listdir(g_var.working_dir+'MERGED'):
        if file_name.endswith('.itp') or file_name.endswith('final.top'):
            copyfile(g_var.working_dir+'MERGED/'+file_name, g_var.final_dir+file_name)


if 'PROTEIN' in cg_residues:
#### creates mdp file if user wants to pull the structure to initial input
    with open(g_var.final_dir+'steered_md.mdp', 'w') as steered_md:
        steered_md.write('define = -DPOSRES\nintegrator = md\nnsteps = 3000\ndt = 0.001\ncontinuation   = no\nconstraint_algorithm = lincs')
        steered_md.write('constraints = h-bonds\nns_type = grid\nnstlist = 25\nrlist = 1\nrcoulomb = 1\nrvdw = 1\ncoulombtype  = PME')
        steered_md.write('pme_order = 4\nfourierspacing = 0.16\ntcoupl = V-rescale\ntc-grps = system\ntau_t = 0.1\nref_t = 310\npcoupl = no')
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


final_time=np.array(gmtime()[3:6])

#### prints out system information

print('\n{:-<100}'.format(''))
print('{0:^100}'.format('Script has completed, time for a beer'))
print('\n{0:^20}{1:^10}'.format('molecules','number'))
print('{0:^20}{1:^10}'.format('---------','------'))
for section in system:
    print('{0:^20}{1:^10}'.format(section, system[section]))


#### prints out script timings for each section
if g_var.v >= 1:
    print('\n{0:^47}{1:^22}'.format('Job','Time'))
    print('{0:^47}{1:^22}'.format('---','----'))
    t1 = np.sqrt((read_in_time-initialisation_time)**2)
    print('\n{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Read in CG system: ',t1[0],'hours',t1[1],'min',t1[2],'sec')) 
    if user_at_input and 'PROTEIN' in system:
        t2=np.sqrt((protein_de_novo_time-read_in_time)**2)
        t3=np.sqrt((final_protein_time-protein_de_novo_time)**2)
        t4=np.sqrt((final_protein_time-read_in_time)**2)
        print('{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Build de novo protein system: ',t2[0],'hours',t2[1],'min',t2[2],'sec'))        
        print('{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Build protein system from provided structure: ',t3[0],'hours',t3[1],'min',t3[2],'sec'))
        print('{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Total protein system build: ',t4[0],'hours',t4[1],'min',t4[2],'sec'))
    else:
        t5=np.sqrt((final_protein_time-read_in_time)**2)
        print('{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Build de novo protein system: ',t5[0],'hours',t5[1],'min',t5[2],'sec'))
    t6=np.sqrt((non_protein_time-final_protein_time)**2)
    t7=np.sqrt((merge_time-non_protein_time)**2)
    t8=np.sqrt((final_time-initialisation_time)**2)
    print('{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Build non protein system: ',t6[0],'hours',t6[1],'min',t6[2],'sec'))
    print('{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Merge protein and non protein system: ', t7[0],'hours',t7[1],'min',t7[2],'sec'))
    print('{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Total run time: ',t8[0],'hours',t8[1],'min',t8[2],'sec'))
