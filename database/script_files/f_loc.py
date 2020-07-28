#!/usr/bin/env python3
import os,sys
import numpy as np
import ntpath
import gen, g_var


gen.correct_number_cpus()
gen.find_gromacs()
forcefield_available, fragments_available = gen.read_database_directories()

if g_var.info:
    gen.database_information(forcefield_available, fragments_available)
elif not os.path.exists(g_var.c):
    sys.exit('Cannot find CG input file: '+g_var.c)



##### forcefield selection
if os.path.exists(g_var.ff):
    forcefield_location, forcefield = gen.path_leaf(g_var.ff)
    gen.folder_copy_and_check(g_var.ff, g_var.final_dir+forcefield)
    print('\nYou have chosen to use your own forcefield: '+forcefield+' in '+forcefield_location)
else:
    try: 
        forcefield_number = forcefield_available.index(g_var.ff.split('.')[0]+'.ff')
    except:
        if g_var.ff != None: 
            print('Cannot find forcefield: '+g_var.ff+'  please select one from below\n')
        forcefield_number = gen.database_selection(forcefield_available, 'forcefields')    
    print('\nYou have selected the forcefield: '+forcefield_available[forcefield_number].split('.')[0])
    gen.folder_copy_and_check(g_var.database_dir+'/forcefields/'+forcefield_available[forcefield_number], g_var.final_dir+forcefield_available[forcefield_number])
    forcefield_location, forcefield=g_var.database_dir+'forcefields/', forcefield_available[forcefield_number]
    if g_var.ff == None:
        g_var.opt['ff'] = forcefield_available[forcefield_number]

##### fragment selection
frag_location, fragment_number, fragments_available = [],[],[]
if g_var.fg != None:
    for frag_val, frag_path in enumerate(g_var.fg):
        if os.path.exists(frag_path):
            frag_loc, fragments = gen.path_leaf(frag_path)
            frag_location.append(frag_loc)
            fragment_number.append(frag_val)
            fragments_available.append(fragments)
if len(fragment_number) == 0:
    fragment_number = gen.fetch_frag_number(fragments_available)
    frag_location = [g_var.database_dir+'fragments/']*len(fragments_available)
    if g_var.fg == None:
        g_var.opt['fg'] = ''
        for database in fragment_number:
            g_var.opt['fg'] += fragments_available[database]+' '

p_directories_unsorted, mod_directories_unsorted, np_directories_unsorted = gen.fetch_residues(frag_location, fragments_available, fragment_number)

np_residues, p_residues, mod_residues, np_directories, p_directories, mod_directories = gen.sort_directories(p_directories_unsorted, 
																						mod_directories_unsorted, np_directories_unsorted)
### reads in water molecules
water_dir, water = gen.check_water_molecules(g_var.w, np_directories)
if g_var.w == None:
    g_var.opt['w'] = water
    ### return backbone information
res_top, sorted_connect, hydrogen, heavy_bond, ions, at_mass = gen.fetch_fragment(p_residues, p_directories, mod_directories,  
                                                                    np_directories, forcefield_location+forcefield, mod_residues)

swap_dict=gen.sort_swap_group()

if g_var.group != None:
    group_chains = gen.fetch_chain_groups()
else:
    group_chains = None
