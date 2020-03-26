#!/usr/bin/env python3

import sys
import numpy as np
import copy
import gen, g_var, f_loc

def read_initial_cg_pdb():
#### initialisation of dictionaries etc
    cg_residues={}  ## dictionary of CG beads eg cg_residues[residue type(POPE)][resid(1)][bead name(BB)][residue_name(PO4)/coordinates(coord)]
    residue_list={} ## a dictionary of bead in each residue eg residue_list[bead name(BB)][residue_name(PO4)/coordinates(coord)]
    count=0  ### residue counter initialisation
    with open(g_var.input_directory+'conversion_input.pdb', 'r') as pdb_input:
        for line in pdb_input.readlines():
#### separates lines
            if line.startswith('ATOM'):
                line_sep = gen.pdbatom(line)
                line_sep['atom_name'], line_sep['residue_name'] = swap(line_sep['atom_name'], line_sep['residue_name'], line_sep['residue_id'])
                if line_sep['atom_name'].upper() != 'SKIP' or line_sep['residue_name'].upper() != 'SKIP':
#### set up resnames in dictionaries
                    cg_residues = add_residue_to_dictionary(cg_residues, line_sep)
    #### sets up previous resid id 
                    if 'residue_prev' not in locals(): 
                        residue_prev=line_sep['residue_id'] 
    #### if resid the same as previous line
                    if residue_prev == line_sep['residue_id']:   ### if resid is the same as the previous line, it adds resname and coordinates to the atom name key in residue_list 
                        residue_list[line_sep['atom_name']]={'residue_name':line_sep['residue_name'],'coord':np.array([line_sep['x'],line_sep['y'],line_sep['z']])}
                        line_sep_prev=line_sep.copy()
    #### if resids are different then the residue list is added to cg_residues
                    else: 
                        if line_sep_prev['residue_name'] not in f_loc.p_residues:
                            cg_residues[line_sep_prev['residue_name']][count]={} ### then create sub dictionary cg_residues[resname][count]
                            cg_residues[line_sep_prev['residue_name']][count]=residue_list ### adds residue list to dictionary key cg_residues[resname][count]
                            if line_sep_prev['residue_name'] == 'ION':
                                cg_residues['SOL'][count]={}
                                sol_res_list={}

                                sol_res_list[f_loc.water]=residue_list[line_sep_prev['atom_name']].copy()
                                sol_res_list[f_loc.water]['residue_name']='SOL'
                                cg_residues['SOL'][count]=sol_res_list
                        else:
                            cg_residues['PROTEIN'][count]={} ### then create sub dictionary cg_residues['PROTEIN'][count]
                            cg_residues['PROTEIN'][count]=residue_list ### adds residue list to dictionary key cg_residues['PROTEIN'][count]
    #### updates dictionaries and counters
                        residue_list={}  ### resets residue list
                        count+=1 ### moves counter along to next residue
                        residue_list[line_sep['atom_name']]={'residue_name':line_sep['residue_name'],'coord':np.array([line_sep['x'],line_sep['y'],line_sep['z']])} ### it adds resname and coordinates to the atom name key in residue_list
                        residue_prev=line_sep['residue_id']   ### updates residue_prev with new resid
                        line_sep_prev=line_sep.copy()
#### finds box vectors
            if line.startswith('CRYST'): ### collects box vectors from pdb
                box_vec=line
#### adds final residue to cg_residues in the same manner as above
    if line_sep['residue_name'] in f_loc.p_residues: 
        if count not in cg_residues['PROTEIN']:
            cg_residues['PROTEIN'][count]={}
        cg_residues['PROTEIN'][count]=residue_list
    else:
        if count not in cg_residues[line_sep['residue_name']]:
            cg_residues[line_sep['residue_name']][count]={}
        cg_residues[line_sep['residue_name']][count]=residue_list
        if line_sep['residue_name'] == 'ION':
            cg_residues['SOL'][count]={}
            sol_res_list={}
            sol_res_list[f_loc.water]=residue_list[line_sep['atom_name']].copy()
            sol_res_list[f_loc.water]['residue_name']='SOL'
            cg_residues['SOL'][count]=sol_res_list
#### checks if box vectors exist
    if 'box_vec' not in locals():### stops script if it cannot find box vectors
        sys.exit('missing box vectors')

    return cg_residues, box_vec

def add_residue_to_dictionary(cg_residues, line_sep):
    if line_sep['residue_name'] in f_loc.p_residues: ## if in protein database 
        if 'PROTEIN' not in cg_residues:  ## if protein does not exist add to dict
            cg_residues['PROTEIN']={}
    elif line_sep['residue_name'].startswith('W') and line_sep['atom_name'].startswith('W'):
        line_sep['residue_name']='SOL'
        if line_sep['residue_name'] not in cg_residues: ## if residue type does not exist add to dict
            cg_residues[line_sep['residue_name']]={}
        line_sep['atom_name']=f_loc.water
    elif line_sep['residue_name'] in f_loc.np_residues:
        if line_sep['residue_name'] not in cg_residues:
            cg_residues[line_sep['residue_name']]={}
        if line_sep['residue_name'] == 'ION' and 'SOL' not in cg_residues:
            cg_residues['SOL']={}
    else:
        sys.exit('\n'+line_sep['residue_name']+' is not in the fragment database!') 
    return cg_residues


def check_new_box(coord,box, new_box):
    for xyz_val, xyz in enumerate(new_box):
        if g_var.box[xyz_val] != 0:
            lower, upper = (float(box[xyz_val])/2)-(float(xyz)/2), (float(box[xyz_val])/2)+(float(xyz)/2)
            if  lower >= coord[xyz_val] or coord[xyz_val] >= upper: 
                return True
    return False

def connect_pbc(temp_coord, prev_coord, box, protein):
#### for x, y, z if the distance between bead is more than half the box length
    if np.sqrt((temp_coord-prev_coord)**2) > float(box)/2:
    #### if the bead if in the opposite 1/3 of the box the position the box length is add/subtracted
        if temp_coord <= float(box)/2:
            temp = temp_coord+float(box)
        elif temp_coord > float(box)/2:
            temp = temp_coord-float(box)
    #### if distance between corrected coordinate is still > 1/2 the box length then counts as a new chain
        if np.sqrt((temp-prev_coord)**2) < 8 or not protein:
            return temp
        else:
            return temp_coord
    else:
        return temp_coord

def fix_pbc(cg_residues, box_vec, new_box, box_shift):
#### fixes box PBC
    box = box_vec.split()[1:4]
    new_box = new_box.split()[1:4]
    for residue_type in cg_residues:
        cut_keys=[]
        for res_val, residue in enumerate(cg_residues[residue_type]):
            for bead_val, bead in enumerate(cg_residues[residue_type][residue]):
                if g_var.box != None and residue_type not in ['PROTEIN']:
                    cut = check_new_box(cg_residues[residue_type][residue][bead]['coord'],box, new_box)
                    if cut:
                        cut_keys.append(residue)
                        break
                if residue_type in ['PROTEIN'] and bead == f_loc.backbone[cg_residues[residue_type][residue][bead]['residue_name']]['BB']:
                    BB_bead = f_loc.backbone[cg_residues[residue_type][residue][bead]['residue_name']]['BB']
                    if res_val != 0 :
                        BB_cur = cg_residues[residue_type][residue][BB_bead]['coord']
                        for xyz in range(3):
                            cg_residues[residue_type][residue][BB_bead]['coord'][xyz] = connect_pbc(BB_cur[xyz], BB_pre[xyz], box[xyz], True)
                    BB_pre = cg_residues[residue_type][residue][BB_bead]['coord'].copy()
                if bead_val != 0 and residue_type not in ['ION','SOL']:
                    for xyz in range(3):
                        cg_residues[residue_type][residue][bead]['coord'][xyz] = connect_pbc(cg_residues[residue_type][residue][bead]['coord'][xyz], 
                                                                                            cg_residues[residue_type][residue][bead_prev]['coord'][xyz], box[xyz], False)
                bead_prev=bead
                if g_var.box != None:
                    cg_residues[residue_type][residue][bead]['coord'] = cg_residues[residue_type][residue][bead]['coord']-box_shift
        for key in cut_keys:
            cg_residues[residue_type].pop(key)
    return cg_residues

def swap(atom, residue, resid):
    if residue in f_loc.swap_dict:
        for key, value in f_loc.swap_dict[residue].items():
            break
        if 'ALL' in f_loc.swap_dict[residue][key]['resid'] or resid in f_loc.swap_dict[residue][key]['resid']:
            if atom in f_loc.swap_dict[residue][key]:
                atom = f_loc.swap_dict[residue][key][atom]
            residue = key.split(':')[1]
    if atom in f_loc.ions and residue != 'ION':
        residue = 'ION'
    return atom, residue

def read_initial_at_pdb():
    at_residues={}  ## dictionary of CG beads eg cg_residues[residue type(POPE)][resid(1)][bead name(BB)][residue_name(PO4)/coordinates(coord)]
    residue_list={} ## a dictionary of bead in each residue eg residue_list[bead name(BB)][residue_name(PO4)/coordinates(coord)]
    count=0  ### residue counter initialisation
    with open(g_var.input_directory+'conversion_input.pdb', 'r') as pdb_input:
        for line in pdb_input.readlines():
            if line.startswith('ATOM'):
                line_sep = gen.pdbatom(line)
                line_sep['atom_name'], line_sep['residue_name'] = swap(line_sep['atom_name'], line_sep['residue_name'], line_sep['residue_id'])
                if line_sep['atom_name'].upper() != 'SKIP' or line_sep['residue_name'].upper() != 'SKIP':
                    if line_sep['residue_name'] in f_loc.p_residues:
                        if 'PROTEIN' not in at_residues:  ## if protein does not exist add to dict
                            at_residues['PROTEIN']={}
                    else:
                        if line_sep['residue_name'] not in at_residues:
                            at_residues[line_sep['residue_name']]={}
                    if 'residue_prev' not in locals(): 
                        residue_prev=line_sep['residue_id'] 
        #### if resid the same as previous line
                    if residue_prev == line_sep['residue_id']:   ### if resid is the same as the previous line, it adds resname and coordinates to the atom name key in residue_list 
                        residue_list[line_sep['atom_name']]={'residue_name':line_sep['residue_name'],'coord':np.array([line_sep['x'],line_sep['y'],line_sep['z']])}
                        line_sep_prev=line_sep.copy()
                    else:
                        if line_sep_prev['residue_name'] not in f_loc.p_residues:
                                at_residues[line_sep_prev['residue_name']][count]={} ### then create sub dictionary cg_residues[resname][count]
                                at_residues[line_sep_prev['residue_name']][count]=residue_list
                        else:
                            at_residues['PROTEIN'][count]={} ### then create sub dictionary cg_residues['PROTEIN'][count]
                            at_residues['PROTEIN'][count]=residue_list ### adds residue list to dictionary key cg_residues['PROTEIN'][count]
    #### updates dictionaries and counters
                        residue_list={}  ### resets residue list
                        count+=1 ### moves counter along to next residue
                        residue_list[line_sep['atom_name']]={'residue_name':line_sep['residue_name'],'coord':np.array([line_sep['x'],line_sep['y'],line_sep['z']])} ### it adds resname and coordinates to the atom name key in residue_list
                        residue_prev=line_sep['residue_id']   ### updates residue_prev with new resid
                        line_sep_prev=line_sep.copy()
#### finds box vectors
            if line.startswith('CRYST'): ### collects box vectors from pdb
                box_vec=line
#### adds final residue to cg_residues in the same manner as above
    if line_sep['residue_name'] in f_loc.p_residues: 
        if count not in at_residues['PROTEIN']:
            at_residues['PROTEIN'][count]={}
        at_residues['PROTEIN'][count]=residue_list
    else:
        if count not in at_residues[line_sep['residue_name']]:
            at_residues[line_sep['residue_name']][count]={}
        at_residues[line_sep['residue_name']][count]=residue_list
    if 'box_vec' not in locals():### stops script if it cannot find box vectors
        sys.exit('missing box vectors')
    return at_residues, box_vec
    