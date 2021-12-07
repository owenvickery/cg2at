#!/usr/bin/env python3

import os, sys
import numpy as np
import math
from distutils.dir_util import copy_tree
import multiprocessing as mp
import distutils.spawn
from shutil import copyfile
import glob
import re
import copy
import ntpath
import g_var


def check_alternate_resname(resname):
    if resname in g_var.alt_res_name:
        return  g_var.alt_res_name[resname]
    elif resname not in g_var.res_top:
        sys.exit('The residue '+resname+' cannot be found in the topology or alternate resnames')
    else:
        return resname

def fetch_forcefield_water_info():
    if g_var.args.info or g_var.args.posre != None or g_var.args.compare != None:
        g_var.get_forcefield = False


def forcefield_selection(test=False):
    ##### forcefield selection
    if g_var.get_forcefield:
        if g_var.args.ff != None:
            if os.path.exists(g_var.args.ff):
                g_var.forcefield_location, g_var.forcefield = path_leaf(g_var.args.ff)
                folder_copy_and_check(g_var.args.ff, g_var.final_dir+g_var.forcefield)
                g_var.opt['ff'] = g_var.forcefield
                if not test:
                    print('\nYou have chosen to use your own forcefield: '+g_var.forcefield+' in '+g_var.forcefield_location)
            elif path_leaf(g_var.args.ff)[1]+'.ff' in g_var.forcefield_available:
                forcefield_number = g_var.forcefield_available.index(path_leaf(g_var.args.ff)[1]+'.ff')
            elif path_leaf(g_var.args.ff)[1] in g_var.forcefield_available:
                forcefield_number = g_var.forcefield_available.index(path_leaf(g_var.args.ff)[1])
            else:
                print('Cannot find forcefield: '+g_var.args.ff+'  please select one from below\n')   
        if 'forcefield_number' not in locals() and g_var.forcefield == '':
            forcefield_number = database_selection(g_var.forcefield_available, 'forcefields', test)    
        if g_var.forcefield == '':
            if not test:
                print('\nYou have selected the forcefield: '+g_var.forcefield_available[forcefield_number].split('.')[0])
            folder_copy_and_check(g_var.database_dir+'/forcefields/'+g_var.forcefield_available[forcefield_number], g_var.final_dir+g_var.forcefield_available[forcefield_number])
            g_var.forcefield_location, g_var.forcefield=g_var.database_dir+'forcefields/', g_var.forcefield_available[forcefield_number]
            g_var.opt['ff'] = g_var.forcefield_available[forcefield_number]

def fragment_selection(test=False):
    ##### fragment selection
    frag_location, fragment_number, fragments_available_other = [],[],[]
    if g_var.args.fg != None:
        for frag_val, frag_path in enumerate(g_var.args.fg):
            if os.path.exists(frag_path):
                frag_loc, fragments = path_leaf(frag_path)
                frag_location.append(os.path.abspath(frag_loc)+'/')
                fragment_number.append(frag_val)
                fragments_available_other.append(fragments)
    if len(fragment_number) == 0:
        fragment_number = fetch_frag_number(g_var.fragments_available, test)
        frag_location = [g_var.database_dir+'fragments/']*len(g_var.fragments_available)
        if g_var.args.fg is None:
            g_var.opt['fg'] = []
            for database in fragment_number:
                g_var.opt['fg'].append(g_var.fragments_available[database])
    else:
        g_var.fragments_available = fragments_available_other
    fetch_residues(frag_location, g_var.fragments_available, fragment_number, test)

def correct_number_cpus():
    if g_var.args.ncpus != None:
        if g_var.args.ncpus > mp.cpu_count():
            print('you have selected to use more CPU cores than are available: '+str(g_var.args.ncpus))
            print('defaulting to the maximum number of cores: '+str(mp.cpu_count()))
            g_var.args.ncpus = mp.cpu_count()
    else:
        if mp.cpu_count() >= 8:
            g_var.args.ncpus = 8
        else:
            g_var.args.ncpus = mp.cpu_count()
    g_var.opt['ncpus'] = g_var.args.ncpus

def check_input_flag():
    if g_var.get_forcefield and g_var.args.c is None:
        g_var.parser.print_help(sys.stderr)
        sys.exit('\nError: the following arguments are required: -c\n')

def path_leaf(path):
    head, tail = ntpath.split(path)
    if not tail:
        return path.replace(ntpath.basename(head)+'/', ''), ntpath.basename(head)
    else:
        return path.replace(tail, ''), tail

### finds gromacs installation
def find_gromacs():
    if g_var.args.gmx != None:
        g_var.args.gmx=distutils.spawn.find_executable(g_var.args.gmx)
    else:
        g_var.args.gmx=distutils.spawn.find_executable('gmx')
    if g_var.args.gmx is None or type(g_var.args.gmx) != str:
        if os.environ.get("GMXBIN") != None:
            for root, dirs, files in os.walk(os.environ.get("GMXBIN")):
                for file_name in files:
                    if file_name.startswith('gmx') and file_name.islower() and '.' not in file_name:
                        g_var.args.gmx=distutils.spawn.find_executable(file_name)
                        if type(g_var.args.gmx) == str and g_var.args.gmx != None :
                            break
                        else:
                            g_var.args.gmx=None
                break
        if g_var.args.gmx is None:
            sys.exit('Cannot find gromacs installation')
    g_var.opt['gmx'] = g_var.args.gmx

def trunc_coord(xyz):
    xyz_new = []
    for coord in xyz:
        if len(str(coord)) > 8:
            if '.' in str(coord):
                xyz_new.append(np.round(coord, 7-len(str(int(coord)))))
            else:
                # print(np.round(coord, 1), coord, 9-len(str(int(coord))))
                xyz_new.append(np.round(coord, 8-len(str(int(coord)))))
        else:
            xyz_new.append(coord)
    return xyz_new[0],xyz_new[1],xyz_new[2]

def calculate_distance(p1, p2):
    return np.sqrt(((p1[0]-p2[0])**2)+((p1[1]-p2[1])**2)+((p1[2]-p2[2])**2))

def file_copy_and_check(file_in,file_out):
    if not os.path.exists(file_out) and os.path.exists(file_in):
        copyfile(file_in, file_out)

def folder_copy_and_check(folder_in,folder_out):
    if not os.path.exists(folder_out):
        copy_tree(folder_in, folder_out)

def flags_used():
    print('\nAll variables supplied have been saved in : \n'+g_var.input_directory+'script_inputs.dat')
    os.chdir(g_var.input_directory)
    with open('script_inputs.dat', 'w') as scr_input:
        scr_input.write('\n'+g_var.opt['input']+'\n')
        for var in g_var.opt:
            if var != 'input':
                scr_input.write('{0:9}{1:15}\n'.format(var,str(g_var.opt[var])))

def is_hydrogen(atom):
    if str.isdigit(atom[0]) and atom[1] != 'H':
        return False
    elif not str.isdigit(atom[0]) and not atom.startswith('H'):
        return False
    else:
        return True

def fetch_chain_groups():
    if g_var.args.group != None:
        if g_var.args.group[0] not in ['all','chain']:
            g_var.group_chains = {}
            for group_val, group in enumerate(g_var.args.group):
                for chain in group.split(','):
                    g_var.group_chains[int(chain)]=group_val 
        else:
            g_var.group_chains =  g_var.args.group[0]   

def split_swap(swap):

    if len(re.split(':', swap)) > 2:
        res_range = re.split(':', swap)[2].split(',')
        res_id = []
        for resid_section in res_range:
            if 'ALL' in resid_section.upper():
                return 'ALL', 'ALL'
            else:
                if '-' in resid_section:
                    spt = resid_section.split('-')
                    for res in range(int(spt[0]), int(spt[1])+1):
                        res_id.append(res)
                else:
                    res_id.append(int(resid_section))
        return res_range, res_id
    else:
        return 'ALL', 'ALL'

def sort_swap_group():
    if g_var.args.swap != None:
        for swap in g_var.args.swap:
            res_s = re.split(':', swap)[0].split(',')
            if re.split(':', swap)[1].split(',') is not type(int):
                res_e = re.split(':', swap)[1].split(',')
            else:
                sys.exit('swap layout is not correct')
            
            if len(res_s) == len(res_e):
                res_range, res_id = split_swap(swap)
                if res_s[0] not in g_var.swap_dict:
                    g_var.swap_dict[res_s[0]]={}
                if len(res_s) == 1:
                    g_var.swap_dict[res_s[0]][res_s[0]+':'+res_e[0]]={'ALL':'ALL'}
                else:
                    g_var.swap_dict[res_s[0]][res_s[0]+':'+res_e[0]]={}
                    for bead in range(1,len(res_s)):
                        g_var.swap_dict[res_s[0]][res_s[0]+':'+res_e[0]][res_s[bead]]=res_e[bead]
                g_var.swap_dict[res_s[0]][res_s[0]+':'+res_e[0]]['resid']=res_id
                g_var.swap_dict[res_s[0]][res_s[0]+':'+res_e[0]]['range']=res_range
            else:
                sys.exit('The length of your swap groups do not match')
        
def print_swap_residues():
    if g_var.args.swap != None:
        to_print = '\nYou have chosen to swap the following residues\n\n'
        to_print += '{0:^10}{1:^5}{2:^11}{3:^11}{4:^11}{5:^11}\n'.format('residue', 'bead', '     ', 'residue', 'bead', 'range')
        to_print += '{0:^10}{1:^5}{2:^11}{3:^11}{4:^11}{5:^11}\n'.format('-------', '----', '     ', '-------', '----', '-----')
        for residue in g_var.swap_dict:
            for swap in g_var.swap_dict[residue]:
                bead_s, bead_e='', ''
                for bead in g_var.swap_dict[residue][swap]:
                    if bead not in ['resid', 'range']:
                        bead_s+=bead+' '
                        bead_e+=g_var.swap_dict[residue][swap][bead]+' '
                    elif bead == 'range':
                        if g_var.swap_dict[residue][swap]['range'] != 'ALL':
                            ran=''
                            for resid_section in g_var.swap_dict[residue][swap]['range']:
                                ran += resid_section+', '
                            ran = ran[:-2]
                        else:
                            ran = g_var.swap_dict[residue][swap]['range']
                to_print += '{0:^10}{1:^5}{2:^11}{3:^11}{4:^11}{5:^11}\n'.format(swap.split(':')[0], bead_s, ' --> ', swap.split(':')[1], bead_e, ran)
        return to_print
    else:
        return ''
    

def new_box_vec(box_vec, box):
    box_vec_split = box_vec.split()[1:]
    box_vec_values, box_shift = [], []
    for xyz_val, xyz in enumerate(box):
        if xyz != 0:
            box_shift.append((float(box_vec_split[xyz_val])/2) - (float(xyz)/2))
            box_vec_values.append(np.round(float(xyz), 3))
        else:
            box_shift.append(0)
            box_vec_values.append(float(box_vec_split[xyz_val]))
    box_vec = g_var.box_line%(float(box_vec_values[0]), float(box_vec_values[1]), float(box_vec_values[2]), float(box_vec_split[3]), float(box_vec_split[4]),float(box_vec_split[5]))
    return box_vec, np.array(box_shift)

def strip_header(line):
    line_new = line.replace('[','').split(']', 1)[0]
    if len(line_new.split())>1 or len(line_new.split())==0:
        sys.exit('There is a issue in one of the fragment headers: \n'+line)
    return line_new.strip()

def topology_header(line, topology, location):
    top = strip_header(line).upper()
    if top in topology:
        return top
    else:
        print('The topology header line is incorrect, therefore ignoring: \n'+location+'.top')
        return ''  

def sep_fragments_topology(residue, location):
    topology={}
    group = 1
    topology = copy.deepcopy(g_var.topology)
    if os.path.exists(location+'.top'):
        with open(location+'.top', 'r') as top_input:
            for line_nr, line in enumerate(top_input.readlines()):
                if not line.startswith('#') and len(line) > 0:
                    if line.startswith('['):
                        topology_function = topology_header(line, topology, location)
                    else:
                        line_sep = line.split()
                        if len(line_sep) > 0:
                            if topology_function == 'GROUPS':
                                topology, group = add_groups(topology, line_sep, group)
                            elif topology_function in ['N_TERMINAL', 'C_TERMINAL']:
                                topology[topology_function] = ''.join(line_sep)
                            elif topology_function == 'CHIRAL':
                                topology = add_chiral(topology, line_sep)
                            elif topology_function == 'CONNECT':
                                topology = add_connections(topology, line_sep)
                            elif topology_function == 'ALT_RES':
                                sort_alternate_residues(line_sep, residue)
                            elif topology_function == 'HYDRATION':
                                sort_hydration(line_sep, residue,)
    return topology

def add_groups(topology, line_sep, group):
    if not g_var.args.mod:
        for bead in line_sep:
            topology['GROUPS'][bead] = group
        group += 1
    topology['GROUPS']['group_max'] = group
    return topology, group

def add_chiral(topology, line_sep):
    if len(line_sep) == 5:
        topology['CHIRAL']['atoms']+=line_sep
        topology['CHIRAL'][line_sep[0]]={'m':line_sep[1], 'c1':line_sep[2], 'c2':line_sep[3], 'c3':line_sep[4]}
    else:
        print('The chiral group section is incorrect: \n'+location+'.top')
    return topology

def add_connections(topology, line_sep):
    if len(line_sep) == 4:
        topology['CONNECT']['atoms'][line_sep[1]] = int(line_sep[3])
        if line_sep[0] in topology['CONNECT']:
            topology['CONNECT'][line_sep[0]]['atom']+=[line_sep[1]]
            topology['CONNECT'][line_sep[0]]['Con_Bd']+=[line_sep[2]]
            topology['CONNECT'][line_sep[0]]['dir']+=[int(line_sep[3])]
        else:
            topology['CONNECT'][line_sep[0]]={'atom':[line_sep[1]], 'Con_Bd':[line_sep[2]], 'dir':[int(line_sep[3])]}
    else:
        print('The bead connection group section is incorrect: \n'+location+'.top')
    return topology

def sort_alternate_residues(line_sep, residue):
    for alt_res in line_sep:
        if alt_res not in g_var.alt_res_name:
            g_var.alt_res_name[alt_res] = residue
        else:
            sys.exit('The alternate residue name: \"'+alt_res+'\" already corresponds to: \"'+
                     g_var.alt_res_name[alt_res]+'\"')

def sort_hydration(line_sep, residue):
    if len(line_sep) == 1 and residue not in g_var.hydration:
        g_var.hydration[residue] = line_sep[0]
    else:
        print('There is a issue with the hydration section of: ', residue)

def get_fragment_topology(residue, location):
    topology = sep_fragments_topology(residue, location[:-4])
    g_var.res_top[residue] = {'ATOMS': [], 'C_TERMINAL':topology['C_TERMINAL'], 'N_TERMINAL':topology['N_TERMINAL'], \
                             'CHIRAL':topology['CHIRAL'], 'GROUPS':{}, 'CONNECT':topology['CONNECT']}
    with open(location, 'r') as pdb_input:
        group=topology['GROUPS']['group_max']
        atom_list=[]
        grouped_atoms={}
        for line_nr, line in enumerate(pdb_input.readlines()):
            if line.startswith('['):
                bead = strip_header(line)
                if bead in topology['GROUPS']:
                    group_temp = topology['GROUPS'][bead]
                else:
                    group_temp = group
                    group+=1
                if group_temp not in grouped_atoms:
                    grouped_atoms[group_temp]={bead:[]}
                    g_var.res_top[residue]['GROUPS'][bead]=group_temp
                else:
                    grouped_atoms[group_temp][bead]=[]
                    g_var.res_top[residue]['GROUPS'][bead]=group_temp
            if line.startswith('ATOM'):
                line_sep = pdbatom(line)
                grouped_atoms[group_temp][bead].append(line_sep['atom_number'])
            if 'bead' not in locals():
                sys.exit('error reading:\n'+location)
            ### return backbone info for each aminoacid residue
            if bead in g_var.res_top[residue]['CONNECT']:
                if line.startswith('ATOM'):
                    line_sep = pdbatom(line)
                    if not is_hydrogen(line_sep['atom_name']):
                        atom_list.append(line_sep['atom_name'])    ### list of backbone heavy atoms
                g_var.res_top[residue]['ATOMS']=atom_list
    return grouped_atoms   

def sort_connectivity(atom_dict, heavy_bond):
    cut_group = {}
    if len(atom_dict) > 1:
        for group in atom_dict:
            cut_group[group]={}
            group_atoms = [atom for frag in atom_dict[group] for atom in atom_dict[group][frag] if atom in heavy_bond]
            non_self_group = [x for x in atom_dict.keys() if x != group]
            for atom in group_atoms:
                for bond in heavy_bond[atom]:
                    for group_2 in non_self_group:
                        for frag in atom_dict[group_2]:                          
                            if bond in atom_dict[group_2][frag]:
                                cut_group[group][atom] = [frag]
    return cut_group

def fetch_fragment_multi():
#### fetches the Backbone heavy atoms and the connectivity with pre/proceeding residues 
    amino_acid_itp = fetch_amino_rtp_file_location(g_var.forcefield_location+g_var.forcefield) 
    g_var.at_mass = fetch_atom_masses(g_var.forcefield_location+g_var.forcefield)
    for residue_type in [g_var.p_directories, g_var.o_directories, g_var.mod_directories]:
        for directory in range(len(residue_type)):
            for residue in residue_type[directory][1:]:    
                if residue not in g_var.res_top:
                    location = fragment_location(residue)
                    grouped_atoms = get_fragment_topology(residue, location)
                    g_var.hydrogen[residue], g_var.heavy_bond[residue], residue_list, at_mass, amide_h = fetch_bond_info(residue, amino_acid_itp, g_var.at_mass, location)
                    
                    g_var.sorted_connect[residue]  = sort_connectivity(grouped_atoms, g_var.heavy_bond[residue])
                    g_var.res_top[residue]['RESIDUE'] = residue_list
                    g_var.res_top[residue]['atom_masses'] = at_mass
                    if residue in g_var.p_residues:
                        g_var.res_top[residue]['amide_h'] = amide_h
                    if residue in g_var.alt_res_name:
                        g_var.res_top[g_var.alt_res_name[residue]]=g_var.res_top[residue]

def fetch_fragment_single():
    for frag_val, frag_dir in enumerate([g_var.sol_directories, g_var.ion_directories, g_var.np_directories]):
        for directory in range(len(frag_dir)):
            for residue in frag_dir[directory][1:]:    
                if residue not in g_var.res_top:
                    location = fragment_location(residue)
                    grouped_atoms = get_fragment_topology(residue, location)
                    if frag_val <=1: 
                        at_mass = fetch_atoms_water_ion(frag_dir[directory][0]+residue+'/', g_var.at_mass)
                        g_var.hydrogen[residue], g_var.heavy_bond[residue], g_var.sorted_connect[residue] = {},{},{}
                        residue_list = [residue]
                    else:
                        g_var.hydrogen[residue], g_var.heavy_bond[residue], residue_list, at_mass, amide_h  = fetch_bond_info(residue, [location[:-4]+'.itp'], g_var.at_mass, location)
                        g_var.sorted_connect[residue]  = sort_connectivity(grouped_atoms, g_var.heavy_bond[residue])
                     
                    g_var.res_top[residue]['RESIDUE'] = residue_list
                    g_var.res_top[residue]['atom_masses'] = at_mass

def atom_bond_check(line_sep):
    if line_sep[1] == 'atoms':
        return True, False
    elif line_sep[1] == 'bonds':
        return False, True
    else:
        return False, False

def fetch_amino_rtp_file_location(forcefield_loc):
    rtp=[]
    for file in os.listdir(forcefield_loc):
        if file.endswith('.rtp'):
            rtp.append(forcefield_loc+'/'+file)
    return rtp

def fetch_atom_masses(forcefield_loc):
    at_mass = {}
    if os.path.exists(forcefield_loc+'/atomtypes.atp'):
        with open(forcefield_loc+'/atomtypes.atp', 'r') as itp_input:
            for line in itp_input.readlines():
                if not (line.isspace() or len(line) == 0) and line[0] not in [';', '#']:
                    line_sep = line.split()       
                    at_mass[line_sep[0]]=line_sep[1] 
    else:
        sys.exit('cannot find atomtypes.dat in the forcefield: '+forcefield_loc)
    return at_mass

def fetch_atoms_water_ion(forcefield_loc, at_mass_p):
    at_mass = {}
    for file in os.listdir(forcefield_loc):
        if file.endswith('itp'):
            with open(forcefield_loc+file, 'r') as itp_input:
                strip_atoms = False
                for line in itp_input.readlines():
                    line_sep = line.split()
                    if not (line.isspace() or len(line) == 0) and line[0] not in [';', '#']:
                        if line.strip().startswith('[') and strip_header(line) == 'atoms':
                            strip_atoms = True
                        elif line.strip().startswith('['):
                            strip_atoms = False
                        elif strip_atoms:
                            at_mass[line_sep[4]] = float(at_mass_p[line_sep[1]])
    return at_mass


def check_res_name(residue, residue_itp):
    if residue_itp in g_var.alt_res_name or residue_itp == residue:
        if residue_itp == residue:
            return True
        elif g_var.alt_res_name[residue_itp] == residue:
            return True
        else:
            return False
    else:
        return False


def fetch_bond_info(residue, rtp, at_mass, location):
    bond_dict=[]
    heavy_dict, H_dict=[],[]
    residue_present, mol_type = False,False
    atom_conversion = {}
    residue_list=[]
    res_at_mass = {}
    for rtp_file in rtp:
        with open(rtp_file, 'r') as itp_input:
            for line in itp_input.readlines():
                if not (line.isspace() or len(line) == 0) and line[0] not in [';', '#']:
                    line_sep = line.split()
                    if line.strip().startswith('[') and not residue_present: 
                        residue_present = check_res_name(residue, strip_header(line))  
                        if strip_header(line) == 'moleculetype':
                            mol_type = True
                    elif mol_type:
                        residue_present = check_res_name(residue, line.split()[0])  
                        mol_type = False
                    elif residue_present or residue not in g_var.p_residues+g_var.o_residues:
                        if line.strip().startswith('[') :
                            atoms, bonds = atom_bond_check(line_sep)
                        elif 'atoms' in locals():
                            if atoms and residue in g_var.p_residues + g_var.o_residues:
                                residue_list, atom_conversion, H_dict, res_at_mass, heavy_dict = fetch_bond_info_atoms_linked(residue, line_sep, residue_list, atom_conversion, H_dict, res_at_mass, heavy_dict)
                            elif atoms:
                                residue_list, H_dict, res_at_mass, heavy_dict = fetch_bond_info_atoms_NP(residue, line_sep, residue_list, H_dict, res_at_mass, heavy_dict)
                            elif bonds:
                                try:
                                    bond_dict.append([int(line_sep[0]),int(line_sep[1])])
                                except:
                                    bond_dict.append([line_sep[0],line_sep[1]])
                            elif not atoms and not bonds and residue in g_var.p_residues+g_var.o_residues:
                                break
            if 'atoms' not in locals():
                print('Issue finding information for residue: ',residue)
                sys.exit('There is a issue with: \n'+rtp_file)
        if len(heavy_dict) > 0:
            break
    if not residue_present:
        sys.exit('cannot find topology information for: '+residue)
    bond_dict=np.array(bond_dict)
    hydrogen = {}
    heavy_bond = {}
    if residue not in g_var.mod_residues:
        atom_conversion = get_atomistic(location)
    for bond in bond_dict:
        hydrogen, amide_h = add_to_topology_list(bond[0], bond[1], hydrogen, heavy_dict, H_dict, atom_conversion, residue, g_var.p_residues+g_var.o_residues)
        if residue in g_var.p_residues and amide_h != None:
            amide_hydrogen = amide_h
        heavy_bond, amide_h= add_to_topology_list(bond[0], bond[1], heavy_bond, heavy_dict, heavy_dict, atom_conversion, residue, g_var.p_residues+g_var.o_residues)
    if 'amide_hydrogen' in locals():
        return hydrogen, heavy_bond, residue_list, res_at_mass, amide_hydrogen
    else:
        return hydrogen, heavy_bond, residue_list, res_at_mass, None

def fetch_bond_info_atoms_linked(residue, line_sep, residue_list, atom_conversion, H_dict, res_at_mass, heavy_dict):
    if residue not in residue_list:
        residue_list.append(residue)
    atom_conversion[line_sep[0]]=int(line_sep[3])+1
    if is_hydrogen(line_sep[0]):
        H_dict.append(line_sep[0])
    else:
        res_at_mass[line_sep[0]] = float(g_var.at_mass[line_sep[1]])
        heavy_dict.append(line_sep[0])
    return residue_list, atom_conversion, H_dict, res_at_mass, heavy_dict

def fetch_bond_info_atoms_NP(residue, line_sep, residue_list, H_dict, res_at_mass, heavy_dict):
    if line_sep[3] not in residue_list:
        residue_list.append(line_sep[3])
    if is_hydrogen(line_sep[4]):
        H_dict.append(int(line_sep[0]))
    else:
        res_at_mass[line_sep[4]] = float(line_sep[7])
        heavy_dict.append(int(line_sep[0]))
    return residue_list, H_dict, res_at_mass, heavy_dict


def get_atomistic(frag_location):
#### read in atomistic fragments into dictionary    
    residue = {} 
    with open(frag_location, 'r') as pdb_input:
        for line_nr, line in enumerate(pdb_input.readlines()):
            if line.startswith('ATOM'):
                line_sep = pdbatom(line) ## splits up pdb line
                residue[line_sep['atom_name']] = line_sep['atom_number']
    return residue


def add_to_topology_list(bond_1, bond_2, top_list, dict1, dict2, conversion, residue, linked_residues):
    amide_hydrogen = None
    for bond in [[bond_1, bond_2], [bond_2, bond_1]]:
        if bond[0] in dict1 and bond[1] in dict2:
            if residue in linked_residues:
                if bond[0] == 'N' and is_hydrogen(bond[1]):
                    amide_hydrogen = bond[1]
                if bond[0] in conversion and bond[1] in conversion:
                    bond[0], bond[1] = conversion[bond[0]],conversion[bond[1]] 
            if bond[0] not in top_list:
                top_list[bond[0]]=[]
            top_list[bond[0]].append(bond[1])
    return top_list, amide_hydrogen

def fragment_location(residue):  
#### runs through dirctories looking for the atomistic fragments returns the correct location
    for res_type in [g_var.np_directories, g_var.p_directories, g_var.mod_directories, g_var.o_directories, g_var.sol_directories, g_var.ion_directories]:
        for directory in range(len(res_type)):
            if os.path.exists(res_type[directory][0]+residue+'/'+swap_to_solvent(residue)+'.pdb'):
                return res_type[directory][0]+residue+'/'+swap_to_solvent(residue)+'.pdb'            
    sys.exit('Cannot find fragment: '+residue+'/'+swap_to_solvent(residue)+'.pdb')


def read_database_directories():
#### Read in forcefields provided
    available_provided_database=[]
    for directory_type in ['forcefields', 'fragments']:
        if os.path.exists(g_var.database_dir+directory_type):
            for root, dirs, files in os.walk(g_var.database_dir+directory_type):
                available_provided = [x for x in sorted(dirs) if not x.startswith('_')]
                break
        else:
            sys.exit('no '+directory_type+' found')
            available_provided=[]
        available_provided_database.append(available_provided)
    g_var.forcefield_available, g_var.fragments_available = available_provided_database[0], available_provided_database[1]


def database_selection(provided, selection_type, test=False):
#### print out selection of forcefields
    if not test:
        print('\n\n{0:^79}\n'.format('Provided '+selection_type))
        print('{0:^90}'.format('{0:^20}{1:^41}'.format('Selection',selection_type)))
        print('{0:^90}'.format('{0:^20}{1:^41}'.format('---------','-'*len(selection_type))))
        for force_num_prov, line in enumerate(provided):
            print('{0:^90}'.format('{0:^20}{1:^41}'.format(force_num_prov,line.split('.')[0])))
    return ask_database(provided,  selection_type)

def ask_database(provided, selection_type, test=False):
#### ask which database to use
    attempt=0
    while True:
        try:
            if len(provided)==1:
                print('\nOnly 1 '+selection_type[:-1]+' database is currently available, therefore you have no choice but to accept the following choice.')
                return number
        #### if asking about fragments accept a list of libraries 
            if selection_type=='fragments': 
                number = np.array(input('\nplease select fragment libraries (in order of importance: eg. "1 0" then ENTER): ').split())
                number=number.astype(int)
                if len(number[np.where(number >= len(provided))]) == 0:
                    return number
                elif test:
                    return True
        #### if forcefield only accept one selection
            else:
                number = int(input('\nplease select a forcefield: '))
                if number < len(provided):
                    return number
                elif test:
                    return True
        except KeyboardInterrupt:
            sys.exit('\nInterrupted')
        except BaseException:
            if test:
                return True
            print("Oops!  That was a invalid choice")
            attempt+=1
            if attempt > 3:
                sys.exit('Too many invalid choices')

def fetch_frag_number(fragments_available, test=False):
    fragment_number = []
    if g_var.args.fg != None and len(g_var.args.fg) > 0 :
        for frag in g_var.args.fg:
            if frag in fragments_available:
                fragment_number.append(fragments_available.index(frag))
            elif not g_var.args.info:
                print('Cannot find fragment library: '+frag+' please select library from below\n')
                fragment_number += database_selection(fragments_available, 'fragments', test).tolist()
    else:
        fragment_number = database_selection(fragments_available, 'fragments', test)
    if len(fragment_number) > 0:
        return fragment_number 
    else:
        if g_var.args.info:
            return []
        sys.exit('no fragment databases selected')

def add_to_list(root, dirs, list_to_add, residues):
    list_to_add.append([])
    list_to_add[-1].append(root)
    list_to_add[-1] += dirs
    list1 = [x for x in list_to_add[-1] if not x.startswith('_')]
    list1.sort()
    list_to_add[-1] = list1
    residues += list_to_add[-1][1:]
    residues.sort()    

def fetch_residues(frag_dir, fragments_available_prov, fragment_number, test=False):
#### list of directories and water types  [[root, folders...],[root, folders...]]
#### run through selected fragments
    for database in fragment_number:
        if not g_var.args.info and not test:
            if g_var.database_dir in frag_dir[database]:
                print('\nYou have selected the fragment library: '+fragments_available_prov[database])
            else:
                print('\nYou have chosen to use your own fragment library: '+fragments_available_prov[database]+' in '+frag_dir[database])
    #### separate selection between provided and user
        location = frag_dir[database]+ fragments_available_prov[database]
    #### runs through protein and non protein
        for directory_type in ['/non_protein/', '/protein/', '/other/', '/protein_modified/', '/solvent/', '/ions/']:
            if os.path.exists(location+directory_type):
                for root, dirs, files in os.walk(location+directory_type):
                    #### adds non protein residues locations to np_directories
                    if directory_type =='/non_protein/':
                        add_to_list(root, dirs, g_var.np_directories, g_var.np_residues)
                    #### adds other mutli residue locations to o_directories
                    elif directory_type =='/other/':
                        add_to_list(root, dirs, g_var.o_directories, g_var.o_residues)
                    #### adds protein residues locations to p_directories
                    elif directory_type =='/protein/':
                        add_to_list(root, dirs, g_var.p_directories, g_var.p_residues)
                    elif directory_type =='/protein_modified/':
                        add_to_list(root, dirs, g_var.mod_directories, g_var.mod_residues)
                        g_var.p_residues+=g_var.mod_residues
                    elif directory_type =='/ions/':
                        add_to_list(root, dirs, g_var.ion_directories, g_var.ion_residues)
                    elif directory_type =='/solvent/':
                        add_to_list(root, dirs, g_var.sol_directories, g_var.sol_residues)
                    break
            else:
                print('Cannot find fragments for: ', directory_type[1:-1] )
    

def print_water_selection(water):
    to_print =''
    if g_var.args.w != None:
        to_print +='\nThe water type '+g_var.args.w+' doesn\'t exist\n'
    to_print +='\nPlease select a water molecule from below:\n\n'
    to_print +='{0:^20}{1:^30}\n'.format('Selection','water_molecule')
    to_print +='{0:^20}{1:^30}\n'.format('---------','----------')
    for selection, water_model in enumerate(water):
        to_print +='{0:^20}{1:^30}\n'.format(selection,water_model)
    return to_print

def ask_for_water_model(water):
    attempt=0
    while True:
        try:
            number = int(input('\nplease select a water model: '))
            if number < len(water):
                return water[number]
        except KeyboardInterrupt:
            sys.exit('\nInterrupted')
        except BaseException:
            print("Oops!  That was a invalid choice")
            attempt+=1
            if attempt > 3:
                sys.exit('Too many invalid choices')

def check_water_molecules(test=False):
    water_info=[]        
    water=np.array([])
    for directory in g_var.sol_directories:
        water_info.append([directory[0]])
        for fragment in directory[1:]:
            for filename in os.listdir(directory[0]+fragment+'/'):
                if filename.endswith('.pdb') and not filename.startswith('_'):
                    water = np.append(water, filename[:-4])
                    water_info[-1].append(filename[:-4])
        water_info[-1] = np.sort(np.unique(np.array(water_info[-1])))
    g_var.water_info = water_info  
    g_var.water = np.sort(np.unique(water))
    if g_var.get_forcefield:
        if g_var.args.w != None:
            g_var.args.w=g_var.args.w.upper()
        if len(water) == 0:
            print('WARNING cannot find any solvent fragments')
        else:
            if g_var.args.w in water:
                print('\nYou have selected the water model: '+g_var.args.w)
            else:
                print(print_water_selection(g_var.water))
                g_var.args.w = ask_for_water_model(g_var.water)
            g_var.opt['w'] = g_var.args.w 

def swap_to_solvent(residue_type):
    if residue_type in g_var.sol_residues:
        return g_var.args.w
    else:
        return residue_type

############################################################################################## fragment rotation #################################################################################

def AnglesToRotMat(theta) :
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

def pdbatom(line):
### get information from pdb file
### atom number, atom name, residue name,chain, resid,  x, y, z
    try:
        return dict([('atom_number',int(line[7:11].replace(" ", ""))),('atom_name',str(line[12:16]).replace(" ", "")),('residue_name',str(line[16:21]).replace(" ", "")),\
            ('chain',line[21]),('residue_id',int(line[22:26])), ('x',float(line[30:38])),('y',float(line[38:46])),('z',float(line[46:54]))])
    except BaseException:
        sys.exit('\npdb line is wrong:\t'+line) 

def create_pdb(file_name):
    pdb_output = open(file_name, 'w')
    pdb_output.write('TITLE     GENERATED BY CG2AT\nREMARK    Please don\'t explode\nREMARK    Good luck\n\
'+g_var.box_vec+'MODEL        1\n')
    return pdb_output

def mkdir_directory(directory):
#### checks if folder exists, if not makes folder
    if not os.path.exists(directory):
        os.mkdir(directory)

def clean(test=False):
#### cleans temp files from residue_types
    if not test:
        print()
    for residue_type in g_var.cg_residues:

        if not test:
            print('Cleaning temp files from : '+residue_type)
        os.chdir(g_var.working_dir+residue_type)
        file_list = glob.glob('*temp*', recursive=True)
        file_list += glob.glob(residue_type+'*pdb', recursive=True)
        for file in file_list:
            if not file.endswith('.tpr') and not file.endswith('_merged.pdb') and 'gmx' not in file:
                os.remove(file)
        os.chdir(g_var.working_dir+residue_type+'/MIN')
        file_list = glob.glob('*temp*', recursive=True)
        file_list += glob.glob('*trr', recursive=True)
        for file in file_list:
            if not file.endswith('.tpr'):
                os.remove(file) 

def fix_time(t1, t2):
    minutes, seconds= divmod(t1-t2, 60)
    if minutes > 60:
        hours, minutes = divmod(minutes, 60)
    else:
        hours = 0
    return '{0:^3}{1:^6}{2:^3}{3:^4}{4:^3}{5:^4}'.format(int(np.round(hours)),'hours',int(np.round(minutes)),'min',int(np.round(seconds,0)),'sec') 

def print_script_timings():
    to_print=[]
    to_print.append('\n{:-<100}'.format(''))
    to_print.append('\n{0:^47}{1:^22}'.format('Job','Time'))
    to_print.append('{0:^47}{1:^22}'.format('---','----'))
    to_print.append('\n{0:47}{1}'.format('Initialisation: ', fix_time(g_var.tc['i_t_e'],g_var.tc['i_t'])))
    to_print.append('{0:47}{1}'.format('Read in CG system: ', fix_time(g_var.tc['r_i_t'],g_var.tc['i_t_e']))) 
    if 'PROTEIN' in g_var.system:
        to_print.append('{0:47}{1}'.format('Build protein systems: ',fix_time(g_var.tc['f_p_t'],g_var.tc['r_i_t'])))
    if 'OTHER' in g_var.system:
        to_print.append('{0:47}{1}'.format('Build other systems: ',fix_time(g_var.tc['f_o_t'],g_var.tc['f_p_t'])))        
    to_print.append('{0:47}{1}'.format('Build non protein system: ',fix_time(g_var.tc['n_p_t'],g_var.tc['f_o_t'])))
    to_print.append('{0:47}{1}'.format('Merge and minimise de novo: ',fix_time(g_var.tc['m_t'],g_var.tc['n_p_t'])))
    to_print.append('{0:47}{1}'.format('NVT on de novo: ',fix_time(g_var.tc['eq_t'],g_var.tc['m_t'])))
    if g_var.args.o in ['all', 'align'] and g_var.user_at_input:
        to_print.append('{0:47}{1}'.format('Creating aligned system: ',fix_time(g_var.tc['a_e'],g_var.tc['a_s'])))
    to_print.append('{:-<69}'.format(''))
    to_print.append('{0:47}{1}'.format('Total run time: ',fix_time(g_var.tc['f_t'],g_var.tc['i_t'])))
    with open(g_var.final_dir+'script_timings.dat', 'w') as time_out:  
        for line in to_print:
            time_out.write(line+'\n')
            if g_var.args.v >= 1:
                print(line)
        print('\nAll script timings have been saved in: \n'+g_var.final_dir+'script_timings.dat\n')

def cg2at_header():
    print('\n{0:^90}\n'.format('CG2AT2 is a fragment based conversion of coarsegrain to atomistic.'))
    print('{0:^90}\n'.format('CG2AT2 version: '+str(g_var.version)))
    print('{0:^90}\n'.format('Last updated : '+str(g_var.script_update)))
    print('{0:^90}'.format('CG2AT2 is written by Owen Vickery'))
    print('{0:^90}'.format('Project leader Phillip Stansfeld'))
    print('\n{0:^90}\n{1:^90}'.format('Contact email address:','cg2at2@gmail.com'))
    print('\n{0:^90}\n{1:^90}\n{2:^90}\n{3:-<90}'.format('Address:','School of Life Sciences, University of Warwick,','Gibbet Hill Road, Coventry, CV4 7AL, UK', ''))
    print('{0:^90}'.format('Please email me any new residues for the database!'))
    print('\n{0:^90}\n{1:^90}'.format('If you are using this script please acknowledge me (Dr Owen Vickery)','and cite the following DOI: 10.5281/zenodo.3890163'))    
    print('\n{0:-<90}\n{1:^90}'.format('', 'File locations'))
    print('\n{0:^90}'.format('Executable: '+g_var.opt['input'].split()[0]))
    print('{0:^90}'.format('Database locations: '+g_var.database_dir))
    print('{0:^90}\n\n{1:-<90}'.format('Script locations: '+g_var.scripts_dir, ''))

def database_information():
    
    to_print = '\n{0:^90}\n{1:-<90}\n\n'.format('The available forcefields within your database are (flag -ff):', '')
    for forcefields in g_var.forcefield_available:
        to_print += '{0:^90}\n'.format(forcefields)
    to_print += '\n\n{0:^90}\n{1:-<90}\n\n'.format('The available fragment libraries within your database are (flag -fg):', '')
    for fragments in g_var.fragments_available:
        to_print += '{0:^90}\n'.format(fragments)   
    if g_var.args.fg != None :
        to_print = fragments_in_use(to_print)
    sys.exit(to_print+'\n\"If all else fails, immortality can always be assured by spectacular error.\" (John Kenneth Galbraith)\n')

def fragments_in_use(to_print=''):
    protein_directories=[]
    if np.any(np.array([g_var.np_directories, protein_directories, g_var.mod_directories, g_var.o_directories, g_var.water_info, g_var.ion_directories], dtype=object)):
        for database_val, database in enumerate(sorted(g_var.args.fg)):
            to_print += '\n\n{0:^90}\n{1:-<90}\n\n'.format('The following residues are available in the database: '+database,'')
            res_type_name = ['Non protein residues', 'Protein residues', 'Modified protein residues', 'Other linked residues', 'Solvent residues', 'Solvent models', 'Ions']
            for res_val, residue in enumerate([g_var.np_directories, g_var.p_directories, g_var.mod_directories, g_var.o_directories, g_var.sol_directories, g_var.water_info, g_var.ion_directories]):
                try:
                    res_type = sorted(residue[database_val][1:])
                    to_print += '\n{0:^90}\n{1:^90}\n'.format(res_type_name[res_val], '-'*len(res_type_name[res_val]))
                    if len(', '.join(map(str, res_type))) <= 80:
                        to_print += '{0:^90}\n'.format(', '.join(map(str, res_type)))
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
                            to_print += '{0:^90}\n'.format(line)
                            start = end
                except:
                    pass
        to_print += '\n{0:-<90}\n\n'.format('')
    return to_print

def write_system_components():
    to_write = '\n{:-<100}\n'.format('')
    to_write += '{0:^100}\n'.format('Script has completed, time for a beer')
    to_write += '\n{0:^10}{1:^25}\n'.format('molecules','number')
    to_write += '{0:^10}{1:^25}\n'.format('---------','------')
    for section in g_var.system:
        to_write += '{0:^10}{1:^25}\n'.format(section, g_var.system[section])
    return to_write

def print_sequnce_info(sys_type):
    sequence_info = [g_var.seq_cg[sys_type], g_var.seq_at[sys_type]] if g_var.user_at_input and sys_type == 'PROTEIN' else [g_var.seq_cg[sys_type]]
    to_print = ''
    for rep_val, rep in enumerate(sequence_info):
        if rep_val == 0 and len(rep) != 0:
            to_print += 'Summary of coarsegrain '+sys_type+' chains\n'
        elif rep_val > 0 and len(rep) != 0:
            to_print += '\nSummary of atomistic '+sys_type+' chains\n'
        else:
            break
        to_print += '\n{0:^15}{1:^12}\n'.format('chain number', 'length of chain') #   \nchain number\tDelta A\t\tno in pdb\tlength of chain')
        to_print += '\n{0:^15}{1:^12}\n'.format('------------', '---------------')
        for chain in rep:
            to_print += '{0:^15}{1:^12}\n'.format(chain, len(rep[chain]))
        to_print += '\nSequences:\n'
        counter=0
        for index in rep:
            rep, to_print, counter = print_sequnce_info_header(rep_val, rep, to_print, counter, index)
            to_print += '{0:9}{1:10}{2:10}{3:10}{4:10}{5:10}{6:10}{7:10}{8:10}{9:10}\n'.format('1','10','20','30','40','50','60','70','80','90')
            to_print = print_to_100_char(rep[index], to_print)
    return to_print

def print_to_100_char(list_to_print, to_print):
    if len(''.join(map(str, list_to_print))) <= 100:
        to_print += '{0:100}\n'.format(''.join(map(str, list_to_print)))
    else:
        start, end = 0, 1                       
        while end < len(list_to_print):
            line = ''.join(map(str, list_to_print[start:end]))
            while len(line) <= 100:
                if end < len(list_to_print):
                    end+=1
                    line = ''.join(map(str, list_to_print[start:end]))
                    if len(line) > 100:
                        end-=1
                        line = ''.join(map(str, list_to_print[start:end]))
                        break
                else:
                    break
            to_print += '{0:100}\n'.format(line)
            start = end
    return to_print

def print_sequnce_info_header(rep_val, rep, to_print, counter, index):
    if rep_val == 0:
        to_print += '\nCG chain: '+str(index)+'\n'
    else:
        to_print += '\nAT chain: '+str(index)+' -> Group '+str(g_var.group_chains[index])+' -> CG chain '+str(g_var.cg_chain_group[index])+'\n'
        if len(g_var.atomistic_protein_input_aligned[g_var.cg_chain_group[index]].keys()) > 1:
            for seq in list(g_var.atomistic_protein_input_aligned[g_var.cg_chain_group[index]].keys()):
                chain_max = int(seq.split(':')[1])
            chain_max = np.max(chain_max)
            seq_range = list(g_var.atomistic_protein_input_aligned[g_var.cg_chain_group[index]].keys())[counter].split(':')
            rep[index] = ['-']*int(seq_range[0])+rep[index]+['-']*(chain_max-int(seq_range[1]))
            counter += 1
        else:
            seq_range = list(g_var.atomistic_protein_input_aligned[g_var.cg_chain_group[index]].keys())[counter].split(':')
            if int(seq_range[0]) > 0:
                rep[index] = ['-']*int(seq_range[0])+rep[index]
            counter = 0
    return rep, to_print, counter
