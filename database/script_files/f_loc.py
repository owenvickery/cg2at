#!/usr/bin/env python3
import os,sys
import numpy as np
import gen, g_var

forcefield_available, fragments_available = gen.read_database_directories()

if g_var.info:
    print('{0:30}'.format('\nThis script is a fragment based conversion of the coarsegrain representation to atomistic.\n'))
    print('{0:^90}'.format('Current version number: 0.00'))
    print('{0:^90}'.format('Written by Owen Vickery'))
    print('{0:^90}'.format('Project leader Phillip Stansfeld'))
    print('\n{0:^90}\n{1:^90}'.format('Contact email address:','owen.vickery@warwick.ox.ac.uk'))
    print('\n{0:^90}\n{1:^90}\n{2:^90}\n{3:-<90}'.format('Address:','School of Life Sciences, University of Warwick,','Gibbet Hill Road, Coventry, CV4 7AL, UK', ''))
    print('\n{0:^90}\n{1:-<90}\n'.format('The available forcefields within your database are (flag -ff):', ''))
    for forcefields in forcefield_available:
        print('{0:^90}'.format(forcefields))
    print('\n\n{0:^90}\n{1:-<90}\n'.format('The available fragment libraries within your database are (flag -fg):', ''))
    for fragments in fragments_available:
        print('{0:^90}'.format(fragments))    
    sys.exit('\n\"If all else fails, immortality can always be assured by spectacular error.\" (John Kenneth Galbraith)\n')


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

try: 
    fragment_number = []
    for frag in g_var.fg:
        fragment_number.append(fragments_available.index(frag))
except:
    if g_var.fg != None: 
        print('Cannot find fragment library: '+frag+' please select library from below\n')
    fragment_number = gen.database_selection(fragments_available, 'fragments')


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




### finds initial rotation matrices
x_rot, y_rot, z_rot=[],[],[]
for angle in range(0,360, 5):
    angle=np.radians(angle)
    x_rot.append(gen.eulerAnglesToRotationMatrix([angle,0,0]))
    y_rot.append(gen.eulerAnglesToRotationMatrix([0,angle,0]))
    z_rot.append(gen.eulerAnglesToRotationMatrix([0,0,angle]))


