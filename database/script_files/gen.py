#!/usr/bin/env python3

import os, sys
import numpy as np
import math
from distutils.dir_util import copy_tree
from shutil import copyfile
import glob
import re
import shlex
import g_var


def calculate_distance(p1, p2):
    return np.sqrt(((p1[0]-p2[0])**2)+((p1[1]-p2[1])**2)+((p1[2]-p2[2])**2))

def file_copy_and_check(file_in,file_out):
    if not os.path.exists(file_out):
        copyfile(file_in, file_out)

def folder_copy_and_check(folder_in,folder_out):
    if not os.path.exists(folder_out):
        copy_tree(folder_in, folder_out)

def flags_used():
    os.chdir(g_var.input_directory)
    with open('script_inputs.dat', 'w') as scr_input:
        for var in g_var.variables_to_save:
            line='{0:15}{1:15}\n'.format(var,str(g_var.variables_to_save[var]))
            scr_input.write(line)

def is_hydrogen(atom):
    if str.isdigit(atom[0]) and atom[1] != 'H':
        return False
    elif not str.isdigit(atom[0]) and not atom.startswith('H'):
        return False
    else:
        return True

def fetch_chain_groups():
    group_chains = {}
    if g_var.group[0] not in ['all','chain']:
        for group_val, group in enumerate(g_var.group):
            for chain in group.split(','):
                group_chains[int(chain)]=group_val 
        return group_chains
    else:
        return g_var.group[0]
    

def split_swap(swap):
    try:
        res_range = re.split(':', swap)[2].split(',')
        res_id = []
        for resid_section in res_range:
            if '-' in resid_section:
                spt = resid_section.split('-')
                for res in range(int(spt[0]), int(spt[1])+1):
                    res_id.append(res)
            else:
                res_id.append(resid_section)
        return res_range, res_id
    except:
        return 'ALL', 'ALL'

def sort_swap_group():
    s_res_d = {}
    if g_var.swap != None:
        s_res_d = {}
        for swap in g_var.swap:
            res_s = re.split(':', swap)[0].split(',')
            if re.split(':', swap)[1].split(',') is not type(int):
                res_e = re.split(':', swap)[1].split(',')
            else:
                sys.exit('swap layout is not correct')
            if len(res_s) == len(res_e):
                res_range, res_id = split_swap(swap)
                if res_s[0] not in s_res_d:
                    s_res_d[res_s[0]]={}
                if len(res_s) == 1:
                    s_res_d[res_s[0]][res_s[0]+':'+res_e[0]]={'ALL':'ALL'}
                else:
                    s_res_d[res_s[0]][res_s[0]+':'+res_e[0]]={}
                    for bead in range(1,len(res_s)):
                        s_res_d[res_s[0]][res_s[0]+':'+res_e[0]][res_s[bead]]=res_e[bead]
                s_res_d[res_s[0]][res_s[0]+':'+res_e[0]]['resid']=res_id
            else:
                sys.exit('The length of your swap groups do not match')
        print_swap_residues(s_res_d, res_range)
    return s_res_d

def print_swap_residues(s_res_d, res_range):
    print('\nYou have chosen to swap the following residues\n')
    print('{0:^10}{1:^5}{2:^11}{3:^11}{4:^11}{5:^11}'.format('residue', 'bead', '     ', 'residue', 'bead', 'range'))
    print('{0:^10}{1:^5}{2:^11}{3:^11}{4:^11}{5:^11}'.format('-------', '----', '     ', '-------', '----', '-----'))
    for residue in s_res_d:
        for swap in s_res_d[residue]:
            bead_s, bead_e='', ''
            for bead in s_res_d[residue][swap]:
                if bead != 'resid':
                    bead_s+=bead+' '
                    bead_e+=s_res_d[residue][swap][bead]+' '
                else:
                    if res_range != 'ALL':
                        ran=''
                        for resid_section in res_range:
                            ran += resid_section+', '
                        ran = ran[:-2]
                    else:
                        ran = res_range
            print('{0:^10}{1:^5}{2:^11}{3:^11}{4:^11}{5:^11}'.format(swap.split(':')[0], bead_s, ' --> ', swap.split(':')[1], bead_e, ran))
    

def new_box_vec(box_vec, box):
    box_vec_split = box_vec.split()[1:4]
    box_vec_values, box_shift = [], []
    for xyz_val, xyz in enumerate(box):
        if xyz != 0:
            box_shift.append((float(box_vec_split[xyz_val])/2) - (float(xyz)/2))
            box_vec_values.append(np.round(float(xyz), 3))
        else:
            box_shift.append(0)
            box_vec_values.append(float(box_vec_split[xyz_val]))
    box_vec = g_var.box_line%(box_vec_values[0], box_vec_values[1], box_vec_values[2])
    return box_vec, np.array(box_shift)

def fetch_chiral(np_directories,p_directories):
    processing={}     
    for directory_type in [np_directories,p_directories]:
        for directory in range(len(directory_type)):
            for residue in directory_type[directory][1:]:   
                if os.path.exists(directory_type[directory][0]+residue+'/chiral.dat') and residue not in processing:
                    processing[residue]={'atoms':[]}
                    with open(directory_type[directory][0]+residue+'/chiral.dat', 'r') as chir_input:            
                        for line_nr, line in enumerate(chir_input.readlines()):
                            if line[0] != '#':
                                line_sep =line.split()
                                if len(line_sep) == 5:
                                    processing[residue]['atoms']+=line_sep
                                    processing[residue][line_sep[0]]={'m':line_sep[1], 'c1':line_sep[2], 'c2':line_sep[3], 'c3':line_sep[4]}
                                else:
                                    sys.exit('The following chiral group file is incorrect: \n'+directory_type[directory][0]+residue+'/chiral.dat')
                else:
                    pass
    return processing

def sep_fragments_header(line, residue_name):
    line = line.replace('[','')
    line = line.replace(']','')
    line_sep = shlex.split(line)
    residue = {}
    residue['ter']=False
    for top in line_sep:
        try:
            if top in g_var.topology:
                t_header = top
                residue[top]={}
            elif top == 'ter':
                residue['ter']=True
            elif top.count(':') >= 1:
                top_split_grouped = top.split()
                for group in top_split_grouped: 
                    group_split=group.split(':')
                    residue[t_header][group_split[0]]=group_split[1].split(',')
            elif top.count(',') >= 1:
                residue[t_header]=top.split(',')
            elif 't_header' in locals() and len(top) > 0 :
                residue[t_header]=top  
        except:
            print('Something is wrong in the residue: ',residue_name,'\n',line)
            sys.exit()
    return residue

def switch_num_name(dictionary, input_val, num_to_letter):
    if num_to_letter:
        dictionary = {v: k for k, v in dictionary.items()}
    if input_val in dictionary:
        return dictionary[input_val]
    else:
        return input_val        

def sort_connectivity(atom_dict, heavy_bond, connect):
    cut_group = {}
    if len(atom_dict) > 1:
        for group in atom_dict:
            cut_group[group]={}
            for frag in atom_dict[group]:
                for atom in atom_dict[group][frag]:
                    if atom in heavy_bond:
                        for bond in heavy_bond[atom]:
                            for group_2 in atom_dict:
                                if group_2 != group:
                                    for frag in atom_dict[group_2]:                          
                                        if bond in atom_dict[group_2][frag]:
                                            cut_group[group][atom] = [frag]

    return cut_group

def fetch_fragment(p_residues, p_directories, mod_directories, np_directories, forcefield_location, mod_residues):
#### fetches the Backbone heavy atoms and the connectivity with pre/proceeding residues 
    amino_acid_itp = fetch_amino_rtp_file_location(forcefield_location) 
    processing={}     ### dictionary of backbone heavy atoms and connecting atoms eg backbone['ASP'][atoms/b_connect]
    sorted_connect={}
    hydrogen = {}
    heavy_bond = {}
    atoms_dict={}
    for directory in range(len(p_directories)):
        for residue in p_directories[directory][1:]:    
            if residue not in processing:
                atoms_dict={}
                location = fragment_location(residue,p_residues, p_directories, mod_directories, np_directories)
                hydrogen[residue], heavy_bond[residue] = fetch_bond_info(residue, amino_acid_itp, mod_residues, p_residues)
                processing, grouped_atoms, heavy_bond[residue], connect = get_fragment_topology(residue, location, processing, heavy_bond)
                sorted_connect[residue]  = sort_connectivity(grouped_atoms, heavy_bond[residue], connect)
                
    for directory in range(len(np_directories)):
        for residue in np_directories[directory][1:]:    
            if residue not in processing:
                atoms_dict={}
                location = fragment_location(residue,p_residues, p_directories, mod_directories, np_directories)
                if residue in ['SOL','ION']: 
                    hydrogen[residue], heavy_bond[residue], atoms_dict = {},{},{}
                else:
                    hydrogen[residue], heavy_bond[residue] = fetch_bond_info(residue, location[:-4]+'.itp', mod_residues, p_residues)
                processing, grouped_atoms, heavy_bond[residue], connect = get_fragment_topology(residue, location, processing, heavy_bond)
                if residue in ['SOL', 'ION']: 
                    sorted_connect[residue]={}
                else:
                    sorted_connect[residue]  = sort_connectivity(grouped_atoms, heavy_bond[residue], connect)
    
    return processing, sorted_connect, hydrogen, heavy_bond 

def atom_bond_check(line_sep):
    if line_sep[1] == 'atoms':
        return True, False
    elif line_sep[1] == 'bonds':
        return False,True
    else:
        return False, False

def fetch_amino_rtp_file_location(forcefield_loc):
    for file in os.listdir(forcefield_loc):
        if file in ['aminoacids.rtp', 'merged.rtp']:
            return forcefield_loc+'/'+file

def fetch_bond_info(residue, rtp, mod_residues, p_residues):
    bond_dict=[]
    heavy_dict, H_dict=[],[]
    residue_present = False
    atom_conversion = {}
    with open(rtp, 'r') as itp_input:
        for line in itp_input.readlines():
            if len(line.split()) >= 2 and not line.startswith(';'):
                line_sep = line.split()
                if line_sep[1] == residue:
                    residue_present = True
                elif line_sep[1] in ['HSE', 'HIE'] and residue == 'HIS': 
                    residue_present = True
                elif residue_present or residue not in p_residues:
                    if line_sep[0] == '[':
                        atoms, bonds = atom_bond_check(line_sep)
                    elif atoms:
                        if residue in p_residues:
                            atom_conversion[line_sep[0]]=int(line_sep[3])+1
                            if is_hydrogen(line_sep[0]):
                                H_dict.append(line_sep[0])
                            else:
                                heavy_dict.append(line_sep[0])
                        else:
                            if is_hydrogen(line_sep[4]):
                                H_dict.append(int(line_sep[0]))
                            else:
                                heavy_dict.append(int(line_sep[0]))
                    elif bonds:
                        try:
                            bond_dict.append([int(line_sep[0]),int(line_sep[1])])
                        except:
                            bond_dict.append([line_sep[0],line_sep[1]])
                    elif not atoms and not bonds and residue in p_residues:
                        break
    bond_dict=np.array(bond_dict)
    hydrogen = {}
    heavy_bond = {}
    if residue in p_residues and residue not in mod_residues:
        at_conv = {}
        for key_val, key in enumerate(heavy_dict):
            atom_conversion[key] = key_val+1

    for bond in bond_dict:
        hydrogen = add_to_topology_list(bond[0], bond[1], hydrogen, heavy_dict, H_dict, atom_conversion, residue, p_residues)
        heavy_bond = add_to_topology_list(bond[0], bond[1], heavy_bond, heavy_dict, heavy_dict, atom_conversion, residue, p_residues)

    return hydrogen, heavy_bond

def add_to_topology_list(bond_1, bond_2, top_list, dict1, dict2, conversion, residue, p_residues):
    for bond in [[bond_1, bond_2], [bond_2, bond_1]]:
        if bond[0] in dict1 and bond[1] in dict2:
            if residue in p_residues:
                bond[0], bond[1] = conversion[bond[0]],conversion[bond[1]] 
            if bond[0] not in top_list:
                top_list[bond[0]]=[]
            top_list[bond[0]].append(bond[1])
    return top_list

def get_fragment_topology(residue, location, processing, heavy_bond):
    with open(location, 'r') as pdb_input:
        processing[residue] = {'C_ter':'C', 'N_ter':'N', 'posres':[], 'ter':False, 'sul':False}
        group=1
        atom_list=[]
        connect={}
        grouped_atoms={}
        for line_nr, line in enumerate(pdb_input.readlines()):
            if line.startswith('['):
                header_line = sep_fragments_header(line, residue)
                if g_var.mod:
                    header_line['group']=group
                    group+=1
                if 'frag' not in header_line or 'group' not in header_line:
                    sys.exit('\nThere is a issue with the fragment header: '+line+'found in: '+location)
                if int(header_line['group']) not in connect:
                    connect[int(header_line['group'])] = {'g_frag':[header_line['frag']]}
                    grouped_atoms[int(header_line['group'])]={header_line['frag']:[]}
                else:
                    connect[int(header_line['group'])]['g_frag'].append(header_line['frag'])
                    grouped_atoms[int(header_line['group'])][header_line['frag']]=[]
                processing = get_posres(residue, processing, header_line)
            if line.startswith('ATOM'):
                line_sep = pdbatom(line)
                grouped_atoms[int(header_line['group'])][header_line['frag']].append(line_sep['atom_number'])
            ### return backbone info for each aminoacid residue
            try:
                if header_line['frag'] == 'BB':
                    if line.startswith('ATOM'):
                        line_sep = pdbatom(line)
                        if not is_hydrogen(line_sep['atom_name']):
                            atom_list.append(line_sep['atom_name'])    ### list of backbone heavy atoms
                    processing[residue]['atoms']=atom_list
            except:
                sys.exit('The residue: '+residue+' is missing fragment information')
    return processing, grouped_atoms, heavy_bond[residue], connect

def get_posres(residue, processing, header_line):
    for top in processing[residue]:
        if top in header_line and not processing[residue][top]:
            if top == 'posres':
                if type(header_line[top]) == list:
                    for atom in header_line[top]:
                        processing[residue][top].append(atom)
                else:
                    processing[residue][top].append(header_line[top])
            else:
                processing[residue][top] = header_line[top]
        elif top == 'posres' and 'posres' in header_line:
            if type(header_line[top]) == list:
                for atom in header_line[top]:
                    processing[residue][top].append(atom)
            else:
                processing[residue][top].append(header_line[top])
    return processing

def fragment_location(residue, p_residues,  p_directories, mod_directories, np_directories):  
#### runs through dirctories looking for the atomistic fragments returns the correct location
    if residue in p_residues:
        for directory in range(len(p_directories)):
            if os.path.exists(p_directories[directory][0]+residue+'/'+residue+'.pdb'):
                return p_directories[directory][0]+residue+'/'+residue+'.pdb'
        for directory in range(len(mod_directories)):
            if os.path.exists(mod_directories[directory][0]+residue+'/'+residue+'.pdb'):
                return mod_directories[directory][0]+residue+'/'+residue+'.pdb'
    else:
        for directory in range(len(np_directories)):
            if os.path.exists(np_directories[directory][0]+residue+'/'+residue+'.pdb'):
                return np_directories[directory][0]+residue+'/'+residue+'.pdb'
    sys.exit('cannot find fragment: '+residue+'/'+residue+'.pdb')


# def get_forcefields_from_gromacs():
# #### maybe implement fetching forcefields from gromacs install
#     gmxdat = os.environ.get("GMXDATA")
#     if gmxdat:
#         if os.path.exists(os.path.join(gmxdat,'gromacs','top')):
#             return os.path.join(gmxdat,'gromacs','top')

#         elif os.path.exists(os.path.join(gmxdat,'top')):
#             return os.path.join(gmxdat,'top')
#         else:
#             print('cannot find gromacs forcefields')
#             return False
#     else:
#         return False

def read_database_directories():
#### maybe implement fetching forcefields from gromacs install  
    # gromacs_directory=get_forcefields_from_gromacs()
    # gromacs_provided=[]
    # if gromacs_directory: 
    #     for root, dirs, files in os.walk(gromacs_directory):
    #         gromacs_provided=dirs
    #         break
#### Read in forcefields provided
    available_provided_database=[]
    for directory_type in ['forcefields', 'fragments']:
        if os.path.exists(g_var.database_dir+directory_type):
            for root, dirs, files in os.walk(g_var.database_dir+directory_type):
                # available_provided=sorted(dirs)
                available_provided = [x for x in sorted(dirs) if not x.startswith('_')]
                break
        else:
            sys.exit('no '+directory_type+' found')
            available_provided=[]
        available_provided_database.append(available_provided)

    return  available_provided_database[0], available_provided_database[1]


def database_selection(provided, selection_type):
#### print out selection of forcefields
    print('\n\n{0:^45}\n'.format('Provided '+selection_type))
    print('{0:^20}{1:^30}'.format('Selection',selection_type))
    print('{0:^20}{1:^30}'.format('---------','----------'))
    for force_num_prov, line in enumerate(provided):
        print('{0:^20}{1:^30}'.format(force_num_prov,line.split('.')[0]))
    return ask_database(provided,  selection_type)

def ask_database(provided, selection_type):
#### ask which database to use
    while True:
        try:
            if len(provided)==1:
                print('\nOnly 1 '+selection_type[:-1]+' database is currently available, therefore you have no choice but to accept the following choice.')
                return 0
        #### if asking about fragments accept a list of libraries 
            if selection_type=='fragments': 
                number = np.array(input('\nplease select fragment libraries (in order of importance: eg. "1 0" then ENTER): ').split())
                number=number.astype(int)
                if len(number[np.where(number >= len(provided))]) == 0:
                    return number
        #### if forcefield only accept one selection
            else:
                number = int(input('\nplease select a forcefield: '))
                if number < len(provided):
                    return number
        except KeyboardInterrupt:
            sys.exit('\nInterrupted')
        except:
            print("Oops!  That was a invalid choice")

def add_to_list(root, dirs, list_to_add):
    list_to_add.append([])
    list_to_add[-1].append(root)
    list_to_add[-1]+=dirs
    list_to_add[-1] = [x for x in list_to_add[-1] if not x.startswith('_')]    
    return list_to_add

def fetch_frag_number(fragments_available):
    try: 
        fragment_number = []
        for frag in g_var.fg:
            fragment_number.append(fragments_available.index(frag))
    except:
        if g_var.fg != None or g_var.info: 
            if g_var.info:
                sys.exit('Cannot find find database: '+frag)
            print('Cannot find fragment library: '+frag+' please select library from below\n')
        fragment_number = database_selection(fragments_available, 'fragments')
    return fragment_number

def fetch_residues(fragments_available_prov, fragment_number):
#### list of directories and water types  [[root, folders...],[root, folders...]]
    np_directories, p_directories,mod_directories=[], [],[]

    if type(fragment_number) == int:
        fragment_number=[0]
#### run through selected fragments
    for database in fragment_number:
        if not g_var.info:
            print('\nYou have selected: '+fragments_available_prov[database])
    #### separate selection between provided and user
        location = g_var.database_dir+'fragments/'+ fragments_available_prov[database]
    #### runs through protein and non protein
        for directory_type in ['/non_protein/', '/protein/']:
    #### adds non protein residues locations to np_directories
            if os.path.exists(location+directory_type):
                for root, dirs, files in os.walk(location+directory_type):
                    if directory_type =='/non_protein/':
                        np_directories = add_to_list(root, dirs, np_directories)
        #### adds protein residues locations to p_directories
                    else:
                        p_directories = add_to_list(root, dirs, p_directories)
                    #### adds modified residues to mod directories and removes MOD from p_directories
                        if os.path.exists(location+directory_type+'MOD/'):
                            p_directories[-1].remove('MOD')
                            for root, dirs, files in os.walk(location+directory_type+'MOD/'):
                                p_directories = add_to_list(root, dirs, p_directories)
                                mod_directories = add_to_list(root, dirs, mod_directories)
                                break
                    break
    return p_directories, mod_directories, np_directories


def sort_directories(p_directories, mod_directories, np_directories):
#### sorts directories alphabetically and creates residue database
    p_residues, np_residues, mod_residues = [],[],[]
    for directory in range(len(mod_directories)):
        mod_directories[directory].sort()
        mod_residues+=mod_directories[directory][1:]
    for directory in range(len(p_directories)):
        p_directories[directory].sort()
        p_residues+=p_directories[directory][1:]
    for directory in range(len(np_directories)):
        np_directories[directory].sort()
        np_residues+=np_directories[directory][1:]
#### if verbose prints all fragments found
    if g_var.v >= 2:
        for directory in range(len(np_directories)):
            print('\nnon protein residues fragment directories found: \n\nroot file system\n')
            print(np_directories[directory][0],'\n\nresidues\n\n',np_directories[directory][1:], '\n')
        for directory in range(len(p_directories)):
            print('\nprotein residues fragment directories found: \n\nroot file system\n')
            print(p_directories[directory][0],'\n\nresidues\n\n',p_directories[directory][1:], '\n')
    return np_residues, p_residues, mod_residues, np_directories, p_directories, mod_directories

def print_water_selection(water_input, water, directory):
    if water_input != None:
        print('\nThe water type '+water_input+' doesn\'t exist')
    if len(water) == 0:
        sys.exit('\nCannot find any water models in: \n\n'+directory[0]+'SOL/'+'\n')
    print('\nPlease select a water molecule from below:\n')
    print('{0:^20}{1:^30}'.format('Selection','water_molecule'))
    print('{0:^20}{1:^30}'.format('---------','----------'))
    for selection, water_model in enumerate(water):
        print('{0:^20}{1:^30}'.format(selection,water_model))

def ask_for_water_model(directory, water):
    while True:
        try:
            number = int(input('\nplease select a water model: '))
            if number < len(water):
                return directory[0]+'SOL/', water[number]
        except KeyboardInterrupt:
            sys.exit('\nInterrupted')
        except:
            print("Oops!  That was a invalid choice")

def check_water_molecules(water_input, np_directories):
    water=[]
    for directory in np_directories:
        if os.path.exists(directory[0]+'SOL/SOL.pdb'):
            with open(directory[0]+'SOL/SOL.pdb', 'r') as sol_input:
                for line_nr, line in enumerate(sol_input.readlines()):
                    if line.startswith('['):
                        frag_header = sep_fragments_header(line, 'SOL')
                        water.append(frag_header['frag'])
    if water_input in water:
        return directory[0]+'SOL/', water_input
    else:
        print_water_selection(water_input, water, directory)
    return ask_for_water_model(directory, water)                         

############################################################################################## fragment rotation #################################################################################

def eulerAnglesToRotationMatrix(theta) :
#### rotation matrices for the rotation of fragments. theta is [x,y,z] in radians     
    R_x = np.array([[1,         0,                  0                   ],
                    [0,         math.cos(theta[0]), -math.sin(theta[0]) ],
                    [0,         math.sin(theta[0]), math.cos(theta[0])  ]
                    ])
         
    R_y = np.array([[math.cos(theta[1]),    0,      math.sin(theta[1])  ],
                    [0,                     1,      0                   ],
                    [-math.sin(theta[1]),   0,      math.cos(theta[1])  ]
                    ])
                 
    R_z = np.array([[math.cos(theta[2]),    -math.sin(theta[2]),    0],
                    [math.sin(theta[2]),    math.cos(theta[2]),     0],
                    [0,                     0,                      1]
                    ])
                                        
    R = np.dot(R_z, np.dot( R_y, R_x ))


    return R

def angle_clockwise(A, B):
#### find angle between vectors
    AB = np.linalg.norm(A)*np.linalg.norm(B)
    A_dot_B = A.dot(B)
    angle = np.degrees(np.arccos(A_dot_B/AB))
#### determinant of A, B
    determinant = np.linalg.det([A,B])

    if determinant < 0: 
        return angle
    else: 
        return 360-angle

############################################################################################## fragment rotation done #################################################################################

def connectivity(bead_number, cg_bead, connect, at_residues, cg_residues, resid):
    at_connections,cg_connections=[],[]
#### finds all beads that the cg_bead is connected to
    try:
        run=np.where(connect[:,0]==cg_bead)
    except:
        if cg_bead == 'BB':
            return [],[], cg_residues[resid][cg_bead]['coord']
        sys.exit('cannot find connectivity for :'+str(cg_bead))
#### center of mass of cg_bead
    center=cg_residues[resid][cg_bead]['coord']
#### loop through bead connections from bead of interest
    for con_test in connect[run]:
        cg_temp=[]
    #### fetch connections which have more than one bead 1 to 2 beads and not self      
        cg=connect[np.where(np.logical_and(connect[:,2]==con_test[2],connect[:,0]!=cg_bead))]
    #### for each connecting bead 
        for con_bead in cg[:,0]:
            cg_temp.append(cg_residues[resid][con_bead]['coord']-center)
    #### average position of connecting bead
        cg_connections.append(np.mean(cg_temp, axis=0))
    #### all atoms with bead connections and self. should only ever be one. 
        at = int(connect[np.where(np.logical_and(connect[:,2]==con_test[2],connect[:,0]==cg_bead))][:,1])
        at_connections.append(at_residues[cg_bead][at]['coord']-center)
    return at_connections, cg_connections, center



def pdbatom(line):
### get information from pdb file
### atom number, atom name, residue name,chain, resid,  x, y, z, backbone (for fragment), connect(for fragment)
    try:
        return dict([('atom_number',int(line[7:11].replace(" ", ""))),('atom_name',str(line[12:16]).replace(" ", "")),('residue_name',str(line[17:21]).replace(" ", "")),\
            ('chain',line[21]),('residue_id',int(line[22:26])), ('x',float(line[30:38])),('y',float(line[38:46])),('z',float(line[46:54]))])
    except:
        sys.exit('\npdb line is wrong:\t'+line) 

def create_pdb(file_name, box_vec):
    pdb_output = open(file_name, 'w')
    pdb_output.write('REMARK    GENERATED BY sys_setup_script\nTITLE     SELF-ASSEMBLY-MAYBE\nREMARK    Good luck\n\
'+box_vec+'MODEL        1\n')
    return pdb_output


def mkdir_directory(directory):
#### checks if folder exists, if not makes folder
    if not os.path.exists(directory):
        os.mkdir(directory)


def clean(cg_residues):
#### cleans temp files from residue_types
    for residue_type in cg_residues:
        if residue_type not in ['SOL', 'ION']:
            print('\ncleaning temp files from : '+residue_type)
            os.chdir(g_var.working_dir+residue_type)
            file_list = glob.glob('*temp*', recursive=True)
            for file in file_list:
                os.remove(file)
            os.chdir(g_var.working_dir+residue_type+'/min')
            file_list = glob.glob('*temp*', recursive=True)
            for file in file_list:
                os.remove(file) 

def fix_time(t1, t2):
    minutes, seconds= divmod(t1-t2, 60)
    if minutes > 60:
        hours, minutes = divmod(minutes, 60)
    else:
        hours = 0
    return int(np.round(hours)), int(np.round(minutes)), int(np.round(seconds,0))

def print_script_timings(tc, system, user_at_input):
    print('\n{:-<100}'.format(''))
    print('\n{0:^47}{1:^22}'.format('Job','Time'))
    print('{0:^47}{1:^22}'.format('---','----'))
    t1 = fix_time(tc['r_i_t'], tc['i_t'])
    print('\n{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Read in CG system: ',t1[0],'hours',t1[1],'min',t1[2],'sec')) 
    print('{:-<69}'.format(''))

    if user_at_input and 'PROTEIN' in system:
        t2=fix_time(tc['p_d_n_t'],tc['r_i_t'])
        t3=fix_time(tc['f_p_t'],tc['p_d_n_t'])
        t4=fix_time(tc['f_p_t'],tc['r_i_t'])
        print('{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Build de novo protein: ',t2[0],'hours',t2[1],'min',t2[2],'sec'))        
        print('{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Build protein from provided structure: ',t3[0],'hours',t3[1],'min',t3[2],'sec'))
        # print('{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Total protein build: ',t4[0],'hours',t4[1],'min',t4[2],'sec'))
    # else:
    #     t5=fix_time(tc['f_p_t'],tc['r_i_t'])
        # print('{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Total protein build: ',t5[0],'hours',t5[1],'min',t5[2],'sec'))
    # print('{:-<69}'.format(''))
    t6=fix_time(tc['n_p_t'],tc['f_p_t'])
    t7=fix_time(tc['m_t'],tc['n_p_t'])
    
    print('{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Build non protein system: ',t6[0],'hours',t6[1],'min',t6[2],'sec'))
    print('{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Equilibrate de novo: ', t7[0],'hours',t7[1],'min',t7[2],'sec'))
    if g_var.o in ['all', 'steer']:
        t8=fix_time(tc['s_e'],tc['s_s'])
        print('{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Creating steered system: ', t8[0],'hours',t8[1],'min',t8[2],'sec'))
    if g_var.o in ['all', 'align']:
        t9=fix_time(tc['a_e'],tc['a_s'])
        print('{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Creating aligned system: ', t9[0],'hours',t9[1],'min',t9[2],'sec'))
    print('{:-<69}'.format(''))
    t10=fix_time(tc['f_t'],tc['i_t'])
    print('{0:47}{1:^3}{2:^6}{3:^3}{4:^4}{5:^3}{6:^4}'.format('Total run time: ',t10[0],'hours',t10[1],'min',t10[2],'sec'))

def database_information(forcefield_available, fragments_available):
    print('{0:30}'.format('\nThis script is a fragment based conversion of the coarsegrain representation to atomistic.\n'))
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
    if g_var.fg != None:
        
        fragment_number = fetch_frag_number(fragments_available)
        p_directories_unsorted, mod_directories_unsorted, np_directories_unsorted = fetch_residues(fragments_available, fragment_number)
        for database_val, database in enumerate(sorted(g_var.fg)):
            print('\n\n{0:^90}\n{1:-<90}\n'.format('The following residues are available in the database: '+database, ''))
            res_type_name = ['protein residues', 'modified protein residues', 'non protein residues']
            for res_val, residue in enumerate([p_directories_unsorted, mod_directories_unsorted, np_directories_unsorted]):
                try:
                    res_type = sorted(residue[database_val][1:])
                    print('\n{0:^90}\n{1:^90}'.format(res_type_name[res_val], '-'*len(res_type_name[res_val])))
                    if len(', '.join(map(str, res_type))) <= 80:
                        print('{0:^90}'.format(', '.join(map(str, res_type))))
                    else:
                        start, end = 0, 1                       
                        while end < len(res_type):
                            line = ', '.join(map(str, res_type[start:end]))
                            while len(line) <= 80:
                                if end < len(res_type):
                                    end+=1
                                    line = ', '.join(map(str, res_type[start:end]))
                                    if len(line) > 80:
                                        end-=1
                                        line = ', '.join(map(str, res_type[start:end]))
                                        break
                                else:
                                    break
                            print('{0:^90}'.format(line))
                            start = end
                except:
                    pass
    sys.exit('\n\"If all else fails, immortality can always be assured by spectacular error.\" (John Kenneth Galbraith)\n')
