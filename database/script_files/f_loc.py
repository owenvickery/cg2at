#!/usr/bin/env python3
import os,sys
import numpy as np
import gen, g_var

forcefield_available, fragments_available = gen.read_database_directories()




if g_var.info:
    gen.database_information(forcefield_available, fragments_available)

##### select forcefield
try: 
    forcefield_number = forcefield_available.index(g_var.ff.split('.')[0]+'.ff')
except:
    if g_var.ff != None: 
        print('Cannot find forcefield: '+g_var.ff.split('.')[0]+'.ff  please select one from below\n')
    forcefield_number = gen.database_selection(forcefield_available, 'forcefields')
print('\nYou have selected: '+forcefield_available[forcefield_number].split('.')[0])
gen.folder_copy_and_check(g_var.database_dir+'/forcefields/'+forcefield_available[forcefield_number], g_var.final_dir+forcefield_available[forcefield_number])
forcefield_location, forcefield=g_var.database_dir+'forcefields/', forcefield_available[forcefield_number]

### reads in and sorts fragment information

fragment_number = gen.fetch_frag_number(fragments_available)
p_directories_unsorted, mod_directories_unsorted, np_directories_unsorted = gen.fetch_residues(fragments_available, fragment_number)

np_residues, p_residues, mod_residues, np_directories, p_directories, mod_directories = gen.sort_directories(p_directories_unsorted, 
																						mod_directories_unsorted, np_directories_unsorted)
chiral = gen.fetch_chiral(np_directories, p_directories)

if not g_var.at2cg:
    ### reads in water molecules
    water_dir, water = gen.check_water_molecules(g_var.w, np_directories)

    ### return backbone information
    backbone, sorted_connect, hydrogen, heavy_bond = gen.fetch_fragment(p_residues, p_directories, mod_directories,  
                                                                    np_directories, forcefield_location+forcefield, mod_residues)
swap_dict=gen.sort_swap_group()

if g_var.group != None:
    group_chains = gen.fetch_chain_groups()
else:
    group_chains = None


### finds initial rotation matrices
x_rot, y_rot, z_rot=[],[],[]
for angle in range(0,360, 5):
    angle=np.radians(angle)
    x_rot.append(gen.eulerAnglesToRotationMatrix([angle,0,0]))
    y_rot.append(gen.eulerAnglesToRotationMatrix([0,angle,0]))
    z_rot.append(gen.eulerAnglesToRotationMatrix([0,0,angle]))


