#!/usr/bin/env python3

import sys, os
import numpy as np
import math
import copy
import gen, g_var



def read_initial_cg_pdb(test=False):
    if not test:
        print('\nThis script is now hopefully doing the following (Good luck):\n')
    residue_list={} ## a dictionary of bead in each residue eg residue_list[bead name(BB)][residue_name(PO4)/coordinates(coord)]
    count=0  ### residue counter initialisation
    with open(g_var.input_directory+'CG_INPUT.pdb', 'r') as pdb_input:
        pdb_lines_atoms, box_vec = filter_input(pdb_input.readlines())
        for line_val, line_sep in enumerate(pdb_lines_atoms):
            if np.round((line_val/len(pdb_lines_atoms))*100,2).is_integer() and not test:
                print('Reading in your CG representation: ',np.round((line_val/len(pdb_lines_atoms))*100,2),'%', end='\r')
            line_sep = update_residue_names(line_sep)
            line_sep['atom_name'], line_sep['residue_name'] = swap(line_sep['atom_name'], line_sep['residue_name'], line_sep['residue_id']) ## implements swap group
            if 'SKIP' not in [line_sep['atom_name'].upper(), line_sep['residue_name'].upper()]:
#### set up resnames in dictionaries
                add_residue_to_dictionary(line_sep)
#### sets up previous resid id 
                if 'residue_prev' not in locals(): 
                    residue_prev=line_sep.copy() 
#### if resid the same as previous line
                if residue_prev['residue_id'] == line_sep['residue_id'] and line_sep['residue_name'] == residue_prev['residue_name'] \
                                                                         and line_sep['atom_name'] not in residue_list:   ### if resid is the same as the previous line, it adds resname and coordinates to the atom name key in residue_list 
                    residue_list[line_sep['atom_name']]={'residue_name':line_sep['residue_name'],
                                                        'coord':np.array([line_sep['x'],line_sep['y'],line_sep['z']])}
                    line_sep_prev=line_sep.copy()
                else: 
                    add_to_cg_database(line_sep_prev, count, residue_list)
#### updates dictionaries and counters
                    residue_list={}  ### resets residue list
                    count+=1 ### moves counter along to next residue
                    residue_list[line_sep['atom_name']]={'residue_name':line_sep['residue_name'],
                                                        'coord':np.array([line_sep['x'],line_sep['y'],line_sep['z']])} ### it adds resname and coordinates to the atom name key in residue_list
                    residue_prev=line_sep.copy()    ### updates residue_prev with new resid
                    line_sep_prev=line_sep.copy()
                if line_val+1 == len(pdb_lines_atoms):
                        residue_list[line_sep['atom_name']]={'residue_name':line_sep['residue_name'],
                                                            'coord':np.array([line_sep['x'],line_sep['y'],line_sep['z']])}
                        add_to_cg_database(line_sep, count, residue_list)
    for key in g_var.cg_residues:
        if len(g_var.cg_residues[key]) == 0:
            sys.exit('there is a issue with the residue type: '+key)
    return box_vec[0]

def filter_input(pdb_lines_raw, CG=True):
    pdb_lines_atoms = [gen.pdbatom(j) for j in pdb_lines_raw if j.startswith('ATOM ')] 
    if len(pdb_lines_atoms) == 0:
        sys.exit('input coarsegrain structure seems to contain no beads')
    if CG:
        box_vec =  [j for j in pdb_lines_raw if j.startswith('CRYST')]
        if len(box_vec) == 0:
            sys.exit('The input file is missing the Box vectors')
        return pdb_lines_atoms, box_vec
    return pdb_lines_atoms


def add_to_cg_database(line_sep_prev, count, residue_list):
    if line_sep_prev['residue_name'] in g_var.sol_residues+g_var.ion_residues+g_var.np_residues :
        g_var.cg_residues[line_sep_prev['residue_name']][count]={} ### then create sub dictionary cg_residues[resname][count]
        g_var.cg_residues[line_sep_prev['residue_name']][count]=residue_list ### adds residue list to dictionary key cg_residues[resname][count]
        if line_sep_prev['residue_name'] in g_var.hydration:
            sol_res_list={}
            sol_res_list[g_var.hydration[line_sep_prev['residue_name']]]=residue_list[line_sep_prev['atom_name']].copy()
            sol_res_list[g_var.hydration[line_sep_prev['residue_name']]]['residue_name']=g_var.hydration[line_sep_prev['residue_name']]
            g_var.cg_residues[g_var.hydration[line_sep_prev['residue_name']]][count]=sol_res_list
    elif line_sep_prev['residue_name'] in g_var.o_residues:
        g_var.cg_residues['OTHER'][count]={} ### then create sub dictionary cg_residues['PROTEIN'][count]
        g_var.cg_residues['OTHER'][count]=residue_list ### adds residue list to dictionary key cg_residues['PROTEIN'][count]
    else:
        g_var.cg_residues['PROTEIN'][count]={} ### then create sub dictionary cg_residues['PROTEIN'][count]
        g_var.cg_residues['PROTEIN'][count]=residue_list ### adds residue list to dictionary key cg_residues['PROTEIN'][count]


def update_residue_names(line_sep):
    if line_sep['residue_name'] =='ION':
        line_sep['residue_name'] = line_sep['atom_name']
    if line_sep['residue_name'] in g_var.alt_res_name:
        line_sep['residue_name'] = g_var.alt_res_name[line_sep['residue_name']]
    if line_sep['residue_name'] in g_var.sol_residues+g_var.ion_residues and line_sep['atom_name'] in g_var.alt_res_name:
        line_sep['atom_name'] = g_var.alt_res_name[line_sep['atom_name']]
    return line_sep

def add_residue_to_dictionary(line_sep):
    if line_sep['residue_name'] in g_var.p_residues: ## if in protein database 
        if 'PROTEIN' not in g_var.cg_residues:  ## if protein does not exist add to dict
            g_var.cg_residues['PROTEIN']={}
    elif line_sep['residue_name'] in g_var.o_residues: ## if in protein database 
        if 'OTHER' not in g_var.cg_residues:  ## if protein does not exist add to dict
            g_var.cg_residues['OTHER']={}
    elif line_sep['residue_name'] in g_var.sol_residues: 
        if line_sep['residue_name'] not in g_var.cg_residues: ## if residue type does not exist add to dict
            g_var.cg_residues[line_sep['residue_name']]={}
    elif line_sep['residue_name'] in g_var.ion_residues: 
        if line_sep['residue_name'] not in g_var.cg_residues: ## if residue type does not exist add to dict
            g_var.cg_residues[line_sep['residue_name']]={}
        if line_sep['residue_name'] in g_var.hydration:
            if g_var.hydration[line_sep['residue_name']] not in g_var.cg_residues:
                g_var.cg_residues[g_var.hydration[line_sep['residue_name']]] = {}
    elif line_sep['residue_name'] in g_var.np_residues:
        if line_sep['residue_name'] not in g_var.cg_residues:
            g_var.cg_residues[line_sep['residue_name']]={}
    elif line_sep['residue_name'] == 'SKIP':
        pass
    else:
        sys.exit('\n'+line_sep['residue_name']+' is not in the fragment database!') 


def check_new_box(coord,box, new_box):
    for xyz_val, xyz in enumerate(new_box):
        if g_var.args.box[xyz_val] != 0:
            lower, upper = (float(box[xyz_val])/2)-(float(xyz)/2), (float(box[xyz_val])/2)+(float(xyz)/2)
            if  lower >= coord[xyz_val] or coord[xyz_val] >= upper: 
                return True
    return False

def real_box_vectors(box_vec):
    x, y, z, yz, xz, xy = [float(i) for i in box_vec.split()[1:7]]
    # return real-space vector
    cyz, cxz, cxy, sxy = math.cos(np.radians(yz)), math.cos(np.radians(xz)), math.cos(np.radians(xy)) , math.sin(np.radians(xy))
    wx, wy             = z*cxz, z*(cyz-cxz*cxy)/sxy
    wz                 = math.sqrt(z**2 - wx**2 - wy**2)
    g_var.r_b_vec = np.array([[x, 0, 0], [y*cxy, y*sxy, 0], [wx, wy, wz]]).T
    g_var.r_b_inv = np.linalg.inv(g_var.r_b_vec).T

def brute_mic(p1, p2):
    result = None
    n = 2
    if gen.calculate_distance(p1, p2) > 10:
        for x in range(-n, n+1):
            for y in range(-n, n+1):
                for z in range(-n, n+1):
                    rp = p2+np.dot(g_var.r_b_vec, [x,y,z])
                    d = gen.calculate_distance(p1, rp)
                    if (result is None) or (result[1] > d):
                        result = (rp, d)
        if result[1] < 10:
            return result[0]
        else:
            return p2
    else:
        return p2

def fix_pbc(box_vec, new_box, box_shift):
#### fixes box PBC
    new_box = new_box.split()[1:4]
    BB_pre_resid = 0
    for residue_type in g_var.cg_residues:
        cut_keys=[]
        print('{:<100}'.format(''), end='\r')
        for res_val, residue in enumerate(g_var.cg_residues[residue_type]):
            if np.round((res_val/len(g_var.cg_residues[residue_type]))*100,2).is_integer():
                print('Fixing PBC of residue type '+residue_type+': ',np.round((res_val/len(g_var.cg_residues[residue_type]))*100,2),'%', end='\r')
            for bead_val, bead in enumerate(g_var.cg_residues[residue_type][residue]):
                bead_info = g_var.cg_residues[residue_type][residue][bead]
                if g_var.args.box != None and residue_type not in ['PROTEIN', 'OTHER']:
                    cut = check_new_box(g_var.cg_residues[residue_type][residue][bead]['coord'],box_vec.split()[1:4], new_box)
                    if cut:
                        cut_keys.append(residue)
                        break
                if residue_type in ['PROTEIN', 'OTHER'] and bead in g_var.res_top[bead_info['residue_name']]['CONNECT']:
                    if residue != BB_pre_resid and residue != 0:
                        BB_pre_resid = residue
                        con_info = g_var.res_top[bead_info['residue_name']]['CONNECT'][bead]
                        for con_val, connection in enumerate(con_info['dir']):
                            if connection < 0:
                                BB_cur = g_var.cg_residues[residue_type][residue][bead]['coord']
                                if residue+connection in g_var.cg_residues[residue_type]:
                                    BB_pre = g_var.cg_residues[residue_type][residue+connection][con_info['Con_Bd'][con_val]]['coord']
                                    g_var.cg_residues[residue_type][residue][bead]['coord'] = brute_mic(BB_pre, BB_cur)

                if bead_val != 0:
                    g_var.cg_residues[residue_type][residue][bead]['coord'] = brute_mic(g_var.cg_residues[residue_type][residue][bead_prev]['coord'],
                                                                                    g_var.cg_residues[residue_type][residue][bead]['coord'])
                bead_prev=bead
                if g_var.args.box != None:
                    g_var.cg_residues[residue_type][residue][bead]['coord'] = g_var.cg_residues[residue_type][residue][bead]['coord']-box_shift
        for key in cut_keys:
            g_var.cg_residues[residue_type].pop(key)
    print('{:<100}'.format(''), end='\r')

def swap(atom, residue, resid):
    if residue in g_var.swap_dict:
        for key, value in g_var.swap_dict[residue].items():
            break
        if 'ALL' in g_var.swap_dict[residue][key]['resid'] or resid in g_var.swap_dict[residue][key]['resid']:
            if atom in g_var.swap_dict[residue][key]:
                atom = g_var.swap_dict[residue][key][atom].upper()
            residue = key.split(':')[1].upper()
    return atom, residue


##################################################################  User supplied protein ##############

def read_in_atomistic(protein):
#### reset location and check if pdb exists  
    os.chdir(g_var.start_dir)
    if not os.path.exists(protein):
        sys.exit('cannot find atomistic protein : '+protein)
#### read in atomistic fragments into dictionary residue_list[0]=x,y,z,atom_name    
    atomistic_protein_input={}
    if g_var.input_directory in protein:
        chain_count=g_var.chain_count
    else:
        chain_count = 0
#### read in pdb
    new_chain= False
    with open(protein, 'r') as pdb_input:
        pdb_lines_atoms = filter_input(pdb_input.readlines(), False)
        atomistic_protein_input[chain_count]={}
        for line_nr, line_sep in enumerate(pdb_lines_atoms):
            if line_sep['residue_name'] in g_var.alt_res_name:
                line_sep['residue_name'] = g_var.alt_res_name[line_sep['residue_name']]
            if not gen.is_hydrogen(line_sep['atom_name']) or line_sep['residue_name'] in g_var.mod_residues:
                if line_sep['residue_name'] in g_var.p_residues:
                    if 'line_sep_prev' not in locals():
                        line_sep_prev = line_sep.copy()
                        line_sep_prev['residue_id'] = 'X'
                    #### sorts out wrong atoms in terminal residues
                    if line_sep['atom_name'] in ['OT', 'O1', 'O2']:
                        line_sep['atom_name']='O'
                #### makes C_terminal connecting atom variable  
                    if 'prev_atom_coord' in locals():
                        line_sep['x'],line_sep['y'],line_sep['z'] = brute_mic(prev_atom_coord, [line_sep['x'],line_sep['y'],line_sep['z']])

                    if line_sep['residue_id'] !=  line_sep_prev['residue_id']:
                        dir_list = []
                        for directions in g_var.res_top[line_sep['residue_name']]['CONNECT']['atoms'].values():
                            dir_list.append(directions)
                        line_sep_prev = line_sep.copy()
                        if not np.any(np.array(dir_list) < 0 ) or new_chain:
                            if len(atomistic_protein_input[chain_count]) != 0 :
                                chain_count+=1
                                atomistic_protein_input[chain_count]={}
                    if line_sep['atom_name'] in g_var.res_top[line_sep['residue_name']]['CONNECT']['atoms']:
                        if g_var.res_top[line_sep['residue_name']]['CONNECT']['atoms'][line_sep['atom_name']] > 0:
                            C_ter=[line_sep['x'],line_sep['y'],line_sep['z']]
                            C_resid=line_sep['residue_id']
                        elif 'C_ter' in locals():
                            N_resid=line_sep['residue_id']
                            N_ter=[line_sep['x'],line_sep['y'],line_sep['z']]
                            dist=gen.calculate_distance(N_ter, C_ter)
                            if C_resid != N_resid and dist > 3.5:
                                del N_ter, C_ter
                                chain_count+=1
                                atomistic_protein_input[chain_count]={} ### new chain key

                    
                    prev_atom_coord = [line_sep['x'],line_sep['y'],line_sep['z']]
                    if line_sep['residue_id'] not in atomistic_protein_input[chain_count]:  ## if protein does not exist add to dict
                        atomistic_protein_input[chain_count][line_sep['residue_id']]={}
                #### adds atom to dictionary, every atom is given a initial mass of zero 
                    atomistic_protein_input[chain_count][line_sep['residue_id']][line_sep['atom_number']]={'coord':np.array([line_sep['x'],line_sep['y'],line_sep['z']]),'atom':line_sep['atom_name'], 'res_type':line_sep['residue_name'],'frag_mass':0, 'resid':line_sep['residue_id']}
                #### if atom is in the backbone list then its mass is updated to the correct one
                    if line_sep['atom_name'] in g_var.res_top[line_sep['residue_name']]['ATOMS']:
                        if line_sep['atom_name'] in line_sep['atom_name'] in g_var.res_top[line_sep['residue_name']]['atom_masses']:
                            atomistic_protein_input[chain_count][line_sep['residue_id']][line_sep['atom_number']]['frag_mass']=g_var.res_top[line_sep['residue_name']]['atom_masses'][line_sep['atom_name']]    
    return atomistic_protein_input, chain_count+1    

def duplicate_chain(test=False):
    if len(g_var.args.d) != 0:
        if not test:
            print('{:<100}'.format(''), end='\r')
            print('Now duplicating the supplied chains')
        for ch_d in g_var.args.d:
            duplicate = [int(x) for x in ch_d.split(':')]
            if len(duplicate) == 2 and duplicate[0] in g_var.atomistic_protein_input_raw:
                if not test:
                    print('Using '+str(duplicate[1])+' copies of the atomistic chain '+str(duplicate[0]))
                for chain_duplication in range(duplicate[1]-1):
                    g_var.atomistic_protein_input_raw[g_var.chain_count]=copy.deepcopy(g_var.atomistic_protein_input_raw[duplicate[0]])
                    g_var.chain_count+=1
            else:
                sys.exit('your atomistic chain duplication input is incorrrect')   
