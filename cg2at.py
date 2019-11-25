#!/usr/bin/env python3
import os, sys
import numpy as np
from subprocess import Popen, PIPE
import subprocess, shlex
from time import gmtime, strftime
import math
import multiprocessing as mp
import argparse
import copy
from shutil import copyfile
from distutils.dir_util import copy_tree
import time
from string import ascii_uppercase
from pathlib import Path
import re
import datetime
import glob
from scipy.spatial import KDTree
import difflib
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/database/script_files')
import gen, gro, at_mod, at_mod_p, at_mod_np, cg_mod, g_var, f_loc


initialisation_time=time.time()

user_at_input = gro.collect_input(g_var.c, g_var.a)

print('\nThis script is now hopefully doing the following (Good luck):\n')

#### read in CG file
print('Reading in your CG representation\n')
cg_residues, box_vec = cg_mod.read_initial_pdb()
cg_residues=cg_mod.fix_pbc(cg_residues, box_vec)
read_in_time=time.time()
system={}
### convert protein to atomistic
if 'PROTEIN' in cg_residues:
    p_system, backbone_coords, final_coordinates_atomistic, sequence=at_mod_p.build_protein_atomistic_system(cg_residues['PROTEIN'], box_vec)
    system['PROTEIN']=p_system['PROTEIN']
    protein_de_novo_time=time.time()
    if user_at_input and 'PROTEIN' in system:
    #### reads in atomistic structure   
        atomistic_protein_input = at_mod_p.read_in_atomistic(g_var.input_directory+'AT_input.pdb', system['PROTEIN'], sequence, True)  
        atomistic_protein_centered, cg_com = at_mod_p.center_atomistic(atomistic_protein_input, backbone_coords)
        at_mod_p.rotate_protein_monomers(atomistic_protein_centered, final_coordinates_atomistic, backbone_coords, cg_com, box_vec)
    gro.minimise_protein(system['PROTEIN'], p_system, user_at_input)
    at_mod_p.merge_protein(system['PROTEIN'], '_novo', box_vec)
    if user_at_input and 'PROTEIN' in system:
        print('\tRunning steered MD on input atomistic structure\n')
    #### runs steered MD on atomistic structure on CA and CB atoms
        for chain in range(system['PROTEIN']):
            gro.steered_md_atomistic_to_cg_coord(chain)
        at_mod_p.merge_protein(system['PROTEIN'], '_at_rep_user_supplied', box_vec)
        
final_protein_time=time.time()

#### converts non protein residues into atomistic and minimises 
if len([key for value, key in enumerate(cg_residues) if key not in ['PROTEIN']]) > 0:
    np_system=at_mod_np.build_atomistic_system(cg_residues, box_vec)
    print('\nThis may take some time....(probably time for a coffee)\n')
    for residue_type in cg_residues:
        if residue_type not in ['PROTEIN', 'ION']:
            print('Minimising individual residues: '+residue_type)
            gro.non_protein_minimise(np_system[residue_type], residue_type)
            gro.merge_minimised(residue_type, np_system, box_vec)
            print('Minimising merged: '+residue_type)
            gro.minimise_merged(residue_type, np_system)
    system.update(np_system)
    build_non_protein_time=time.time()

non_protein_time=time.time()

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
    merge_time=time.time()

#### copies all itp files and topologies from whereever they are stored
    for file_name in os.listdir(g_var.working_dir+'MERGED'):
        if file_name.endswith('.itp') or file_name.endswith('final.top'):
            copyfile(g_var.working_dir+'MERGED/'+file_name, g_var.final_dir+file_name)

#### calculates final RMS
if 'PROTEIN' in cg_residues:
    with open(g_var.final_dir+'steered_md.mdp', 'w') as steered_md:
        steered_md.write('define = -DPOSRES\nintegrator = md\nnsteps = 3000\ndt = 0.001\ncontinuation   = no\nconstraint_algorithm = lincs\nconstraints = h-bonds\nns_type = grid\nnstlist = 25\n\
rlist = 1\nrcoulomb = 1\nrvdw = 1\ncoulombtype  = PME\npme_order = 4\nfourierspacing = 0.16\ntcoupl = V-rescale\ntc-grps = system\ntau_t = 0.1\nref_t = 310\npcoupl = no\n\
pbc = xyz\nDispCorr = no\ngen_vel = yes\ngen_temp = 310\ngen_seed = -1')    
    RMSD={}
    if len(cg_residues['PROTEIN'])>0:
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

if g_var.clean:
    gen.clean(cg_residues)

final_time=time.time()

#### prints out system information

print('\n{:-<100}'.format(''))
print('{0:^100}'.format('Script has completed, time for a beer'))
print('\n{0:^20}{1:^10}'.format('molecules','number'))
print('{0:^20}{1:^10}'.format('---------','------'))
for section in system:
    print('{0:^20}{1:^10}'.format(section, system[section]))
print()


#### prints out script timings for each section
if g_var.v == 1:
    print('\nRead in CG system: ', str(datetime.timedelta(minutes=np.round(read_in_time-initialisation_time, 2))).rsplit(':', 1)[0], ' min')
    if user_at_input and 'PROTEIN' in system:
        print('Build de novo protein system: ', str(datetime.timedelta(minutes=np.round(protein_de_novo_time-read_in_time, 2))).rsplit(':', 1)[0], ' min',\
        '\nBuild protein system from provided structure: ', str(datetime.timedelta(minutes=np.round(final_protein_time-protein_de_novo_time, 2))).rsplit(':', 1)[0], ' min', \
        '\nTotal protein system build: ', str(datetime.timedelta(minutes=np.round(final_protein_time-read_in_time, 2))).rsplit(':', 1)[0], ' min')
    else:
        print('Build de novo protein system: ', str(datetime.timedelta(minutes=np.round(final_protein_time-read_in_time, 2))).rsplit(':', 1)[0], ' min')
    print('Build non protein system: ', str(datetime.timedelta(minutes=np.round(non_protein_time-final_protein_time, 2))).rsplit(':', 1)[0], ' min', \
    '\nMerge protein and non protein system: ', str(datetime.timedelta(minutes=np.round(merge_time-non_protein_time, 2))).rsplit(':', 1)[0], ' min')
    print('Total run time: ', str(datetime.timedelta(minutes=np.round(final_time-initialisation_time, 2))).rsplit(':', 1)[0], ' min')