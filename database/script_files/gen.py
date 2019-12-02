#!/usr/bin/env python3

import os, sys
import numpy as np
import math
from distutils.dir_util import copy_tree
import glob
import g_var


def fetch_chiral(np_directories,p_directories):
    processing={}     
    for directory_type in [np_directories,p_directories]:
        for directory in range(len(directory_type)):
            for residue in directory_type[directory][1:]:   
                if residue not in processing:
                    if os.path.exists(directory_type[directory][0]+residue+'/chiral.dat'):
                        processing[residue]={'atoms':[]}
                        with open(directory_type[directory][0]+residue+'/chiral.dat', 'r') as chir_input:            #     for file_name in os.listdir(p_directories[directory][0]+residue):
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




def fetch_fragment(p_directories):
#### fetches the Backbone heavy atoms and the connectivity with pre/proceeding residues 
    
    processing={}     ### dictionary of backbone heavy atoms and connecting atoms eg backbone['ASP'][atoms/b_connect]
    for directory in range(len(p_directories)):
        for residue in p_directories[directory][1:]:    
            if residue not in processing:
                atom_list, bb_list, restraint, disulphide, terminal, dihedral=[], [], [], '', False, []
                for file_name in os.listdir(p_directories[directory][0]+residue):
                    if file_name.endswith('.pdb'):
                        with open(p_directories[directory][0]+residue+'/'+file_name, 'r') as pdb_input:
                            for line_nr, line in enumerate(pdb_input.readlines()):
                                if line.startswith('ATOM'):
                                    line_sep = pdbatom(line)
                                    if 'H' not in line_sep['atom_name']:
                                        if file_name.startswith('B'):
                                            atom_list.append(line_sep['atom_name'])    ### list of backbone heavy atoms
                                        if line_sep['backbone'] == 2:
                                            bb_list.append(line_sep['atom_name'])  ### connecting atoms
                                        if line_sep['backbone'] in [3,4,5]:
                                            restraint.append(line_sep['atom_name'])  ### position restrained atoms
                                        if line_sep['backbone'] in [5]:
                                            disulphide = line_sep['atom_name'] ### position restrained atoms
                                        if line_sep['backbone'] == 4:
                                            terminal=True
                                        if line_sep['backbone'] in [6] and file_name.startswith('B'):
                                            dihedral.append(line_sep['atom_name'])
                processing[residue]={'atoms':atom_list,'b_connect':bb_list,'restraint':restraint, 'disulphide':disulphide, 'ter':terminal, 'dihedral':dihedral}  ### adds heavy atoms and connecting atoms to backbone dictionary 
                atom_list, bb_list, restraint=[], [], []  ### resets residue lists of heavy atoms, connecting atoms and restraint 
#### if verbose prints out all heavy atoms and connecting atoms for each backbone

    if g_var.v >= 2:
        print('\n{:-<75}'.format('>  Verbose level 2 start'))
        print('backbone atoms for each residue and connecting atoms:\n')
        for residue in processing:
            print(residue, '\tbackbone atoms:', processing[residue]['atoms'], '\n\tbackbone connecting atoms:',
                  processing[residue]['b_connect'], '\n\trestrained atoms:', processing[residue]['restraint'],
                  '\n\tTerminal residue:', processing[residue]['ter'],'\n\tdihedral atoms:', processing[residue]['dihedral'],'\n')
        print('\n{:-<75}'.format('>  Verbose level 2 end\n'))
    return processing

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
                available_provided=dirs
                break
        else:
            available_provided=[]
        available_provided_database.append(available_provided)

    if len(available_provided_database[0]) != 0:
        if len(available_provided_database[1]) != 0:
            return  available_provided_database[0], available_provided_database[1]
        else:
            sys.exit('no fragments found')
    else:
        sys.exit('no forcefields found')

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


def sort_forcefield(forcefield_available_prov, f_number):
#### returns forcefield location and forcefield name
#### if forcefield selection is in provided copy forcefield to FORCEFIELD and FINAL directories
    print('\nYou have selected: '+forcefield_available_prov[f_number].split('.')[0])
    copy_tree(g_var.database_dir+'/forcefields/'+forcefield_available_prov[f_number], g_var.working_dir+'FORCEFIELD/'+forcefield_available_prov[f_number]+'/.')
    copy_tree(g_var.database_dir+'/forcefields/'+forcefield_available_prov[f_number], g_var.final_dir+forcefield_available_prov[f_number]+'/.')
    return g_var.database_dir+'/forcefields/', forcefield_available_prov[f_number].split('.')[0]

def fetch_residues(fragments_available_prov, fragment_number):
#### list of directories and water types  [[root, folders...],[root, folders...]]
    np_directories, p_directories,mod_directories=[], [],[]

    if type(fragment_number) == int:
        fragment_number=[0]

#### run through selected fragments
    for database in fragment_number:
        print('\nYou have selected: '+fragments_available_prov[database])
    #### separate selection between provided and user
        location = g_var.database_dir+'fragments/'+ fragments_available_prov[database]
    #### runs through protein and non protein
        for directory_type in ['/non_protein/', '/protein/']:
    #### adds non protein residues locations to np_directories
            if os.path.exists(location+directory_type):
                for root, dirs, files in os.walk(location+directory_type):
                    if directory_type =='/non_protein/':
                        np_directories.append([])
                        np_directories[-1].append(root)
                        np_directories[-1]+=dirs
        #### adds protein residues locations to p_directories
                    else:
                        p_directories.append([])
                        p_directories[-1].append(root)
                        p_directories[-1]+=dirs 
                    #### adds modified residues to mod directories and removes MOD from p_directories
                        if os.path.exists(location+directory_type+'MOD/'):
                            p_directories[-1].remove('MOD')
                            p_directories.append([])
                            for root, dirs, files in os.walk(location+directory_type+'MOD/'):
                                p_directories[-1].append(root)
                                p_directories[-1]+=dirs
                                mod_directories.append([])
                                mod_directories[-1].append(root)
                                mod_directories[-1]+=dirs
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
    if g_var.v >= 1:
        print('\n{:-<75}'.format('>  Verbose level 1 start'))
        for directory in range(len(np_directories)):
            print('\nnon protein residues fragment directories found: \n\nroot file system\n')
            print(np_directories[directory][0],'\n\nresidues\n\n',np_directories[directory][1:], '\n')
        for directory in range(len(p_directories)):
            print('\nprotein residues fragment directories found: \n\nroot file system\n')
            print(p_directories[directory][0],'\n\nresidues\n\n',p_directories[directory][1:], '\n')
        print('\n{:-<75}'.format('>  Verbose level 1 end\n'))

    return np_residues, p_residues, mod_residues, np_directories, p_directories, mod_directories

def check_water_molecules(water_input, np_directories):
    exists = False
    if water_input != None:
        for directory in np_directories:
            if os.path.exists(directory[0]+'SOL/'+water_input+'.pdb'):
                return directory[0]+'SOL/', water_input
    if not exists:
        if water_input != None:
            print('\nThe water type '+water_input+' doesn\'t exist')
        water=[]
        for directory in np_directories:
            water.append([directory[0]+'SOL/'])
            for root, dirs, files in os.walk(directory[0]+'SOL/'):
                for file in files:
                    if file.endswith('.pdb'):
                        water[-1].append(file)
                break
        print('\nPlease select a water molecule from below:\n')
        print('{0:^20}{1:^30}'.format('Selection','water_molecule'))
        print('{0:^20}{1:^30}'.format('---------','----------'))
        offset=0
        for d_val, directory in enumerate(np_directories):
            print('the following water models are found in: \n\n'+water[d_val][0]+'\n')
            for selection, water_file in enumerate(water[d_val][1:]):
                print('{0:^20}{1:^30}'.format(selection+offset,water_file.split('.')[0]))
            offset+=len(water[d_val][1:])
    while True:
        try:
            number = int(input('\nplease select a water molecule: '))
            minim = 0
            for water_val, water_selection in enumerate(water):
                if minim <= number < minim + len(water_selection):
                    return water_selection[0], water_selection[number-minim+water_val+1].split('.')[0]
                minim += len(water_selection)
        except KeyboardInterrupt:
            sys.exit('\nInterrupted')
        except:
            print("Oops!  That was a invalid choice")    

############################################################################################## fragment rotation #################################################################################

def eulerAnglesToRotationMatrix(theta) :
#### rotaion matrices for the rotation of fragments. theta is [x,y,z] in radians     
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
            ('chain',line[21]),('residue_id',int(line[22:26])), ('x',float(line[30:38])),('y',float(line[38:46])),('z',float(line[46:54])), ('backbone',int(float(line[55:61]))),\
            ('connect',int(float(line[62:67])))])
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


