#!/usr/bin/env python3

import sys
import numpy as np
import copy
import gen, g_var, f_loc

def read_initial_pdb():
#### initialisation of dictionaries etc
    cg_residues={}  ## dictionary of CG beads eg cg_residues[residue type(POPE)][resid(1)][bead name(BB)][residue_name(PO4)/coordinates(coord)]
    residue_list={} ## a dictionary of bead in each residue eg residue_list[bead name(BB)][residue_name(PO4)/coordinates(coord)]
    box_line="CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00 P 1           1\n"  ## box vectors format for pdbs
    count=0  ### residue counter initialisation
    with open(g_var.input_directory+'CG_input.pdb', 'r') as pdb_input:
        for line in pdb_input.readlines():
#### separates lines
            if line.startswith('ATOM'):
                line_sep = gen.pdbatom(line)
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
                            sol_res_list[f_loc.water]=residue_list[line_sep_prev['atom_name']]
                            sol_res_list[f_loc.water]['residue_name']='SOL'
                            cg_residues['SOL'][count]=sol_res_list
                    else:
                        for bead in residue_list:
                            if bead.startswith('B'):
                                residue_list['BB']= residue_list.pop(bead)
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
            sol_res_list[f_loc.water]=residue_list[line_sep['atom_name']]
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
        # if line_sep['residue_name'] in water: ## renames waters to SOL saves headache later on
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



def fix_pbc(cg_residues, box_vec):
#### fixes box PBC
    box = box_vec.split()[1:4]
    # print(residues_all)
    for residue_type in cg_residues:
        if residue_type not in ['ION','SOL']:
            for res_val, residue in enumerate(cg_residues[residue_type]):
                for bead_val, bead in enumerate(cg_residues[residue_type][residue]):
                    if bead_val != 0:
                        for xyz in range(3):
                        #### for x, y, z if the distance between bead is more than half the box length
                            if np.sqrt((cg_residues[residue_type][residue][bead]['coord'][xyz]-cg_residues[residue_type][residue][bead_prev]['coord'][xyz])**2) > float(box[xyz])/2:
                            #### if the bead if in the opposite 1/3 of the box the position the box length is add/subtracted
                                if cg_residues[residue_type][residue][bead]['coord'][xyz] <= float(box[xyz])/2:
                                    temp = cg_residues[residue_type][residue][bead]['coord'][xyz]+float(box[xyz])
                                elif cg_residues[residue_type][residue][bead]['coord'][xyz] > float(box[xyz])/2:
                                    temp = cg_residues[residue_type][residue][bead]['coord'][xyz]-float(box[xyz])
                            #### if distance between corrected coordinate is still > 1/2 the box length then counts as a new chain
                                if np.sqrt((temp-cg_residues[residue_type][residue][bead_prev]['coord'][xyz])**2) > 8:
                                    bead_prev=bead
                                else:
                                    cg_residues[residue_type][residue][bead]['coord'][xyz] = temp
                    if residue_type != 'PROTEIN':
                        bead_prev=bead
                    elif res_val == 0 and residue_type == 'PROTEIN':
                        bead_prev=bead
    return cg_residues
