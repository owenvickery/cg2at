#!/usr/bin/env python3

import os, sys
import numpy as np
import math
import itertools
from scipy.spatial import cKDTree
import re
import gen, g_var, f_loc, at_mod_p, read_in


### sanity checking

def sanity_check_fragments(res, cg, sin_bead):
#### fetches bead and atom info from fragment database
    location = fragment_location(res)
    residue, fragment_mass = get_atomistic(location)
    atom_list = []
    bead_list = []
    for group in residue:
        for bead in residue[group]:
            bead_list.append(bead)
            if not sin_bead:
                for atom in residue[group][bead]:
                    atom_list.append(atom)
            else:
                if bead == sin_bead:
                    for atom in residue[group][bead]:
                        atom_list.append(atom)
    return bead_list, atom_list

def sanity_check_atoms(atom_list, res):
#### checks atom order
    for at_num in range(1, len(atom_list)+1):
        if at_num not in atom_list:
            sys.exit('atom number '+str(at_num)+' is missing from fragment library: '+res+'\n')

def sanity_check_beads(bead_list, cg, res):
#### checks if bead is in fragment library
    bead_list=[]
    for bead in cg:
        bead_list.append(bead)
        if bead not in bead_list:
            sys.exit('The bead '+bead+' is missing from the fragment library: '+res+'\n')   
    return bead_list

def sanity_check():
#### runs through every bead and checks whether it exists
    for res_type in g_var.cg_residues:
        if res_type == 'PROTEIN':
            for residue in g_var.cg_residues['PROTEIN']:
                for bead in g_var.cg_residues['PROTEIN'][residue]:
                    resname = g_var.cg_residues['PROTEIN'][residue][bead]['residue_name']
                    break
                bead_list, atom_list = sanity_check_fragments(resname, g_var.cg_residues['PROTEIN'][residue], False)
                bead_list_cg = sanity_check_beads(bead_list, g_var.cg_residues['PROTEIN'][residue], resname)  
                sanity_check_atoms(atom_list, resname)
                if sorted(bead_list) != sorted(bead_list_cg):
                    if len(bead_list) == len(bead_list_cg):
                        fix_atom_wrap(bead_list, bead_list_cg, 'PROTEIN', residue)
                    else:
                        print('There is a issue with residue: '+resname+' '+str(residue+1))
                        sys.exit('number of atomistic fragments: '+str(len(bead_list))+' does not equal number of CG beads: '+str(len(bead_list_cg)))
        elif res_type in ['SOL', 'ION']:
            for residue in g_var.cg_residues[res_type]:
                for bead in g_var.cg_residues[res_type][residue]:
                    sin_bead=bead
                    break
                bead_list, atom_list = sanity_check_fragments(res_type, g_var.cg_residues[res_type][residue], sin_bead)
                sanity_check_beads(bead_list, g_var.cg_residues[res_type][residue], res_type) 
                sanity_check_atoms(atom_list, res_type)
        else:
            bead_list, atom_list = sanity_check_fragments(res_type, g_var.cg_residues[res_type], False)
            sanity_check_atoms(atom_list, res_type)
            for residue in g_var.cg_residues[res_type]:
                bead_list_cg = sanity_check_beads(bead_list, g_var.cg_residues[res_type][residue], res_type)
                if sorted(bead_list) != sorted(bead_list_cg):
                    if len(bead_list) == len(bead_list_cg):
                        fix_atom_wrap(bead_list, bead_list_cg, res_type, residue)
                    else:
                        print('There is a issue with residue: '+res_type+' '+str(residue+1))
                        sys.exit('number of atomistic fragments: '+str(len(bead_list))+' does not equal number of CG beads: '+str(len(bead_list_cg)))
                        
def fix_atom_wrap(bead_list_frag, bead_list_cg, section, resid):
    for bead in bead_list_cg:
        if bead not in bead_list_frag:
            new_bead = bead[1:]+bead[0]
            if new_bead in bead_list_frag and new_bead not in bead_list_cg:
                g_var.cg_residues[section][resid][new_bead] = g_var.cg_residues[section][resid][bead]
                del g_var.cg_residues[section][resid][bead]
            else:
                print('There is a issue with residue: '+section+' '+str(resid+1))
                print('cannot find: '+bead+' or '+new_bead+' in fragment list:')
                sys.exit(bead_list_frag)

#####  Sanity check end


def rotate_atom(coord, center,xyz_rot_apply):
#### rotates atom around center
    coord =  coord-center  #### centers COM coordinates to 0,0,0
    if np.any(xyz_rot_apply): 
        coord=np.array(coord).dot(xyz_rot_apply)
    coord =  coord+center #### translates coord back by original offset
    return coord

def kabsch_rotate(at_connections,cg_connections):
    ## A solution for the best rotation to relate two sets of vectors. By WOLFGANG KABSCH
    ## http://en.wikipedia.org/wiki/Kabsch_algorithm
    Cov_mat = np.dot(np.transpose(at_connections), cg_connections)
    V, S, W = np.linalg.svd(Cov_mat)
    reflection = np.linalg.det(V) * np.linalg.det(W)
    if reflection < 0:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    rot_mat = np.dot(V, W)
    return rot_mat

def add_to_sequence(sequence, residue, chain_count):
    if residue  not in f_loc.mod_residues:
        sequence[chain_count]+=g_var.aas[residue]
    else:
        sequence[chain_count]+='X'    
    return sequence

def overlapping_atoms(tree):
    overlapped_ndx = tree.query_ball_tree(tree, r=0.3)
    overlapped_cut = [ndx for ndx in overlapped_ndx if len(ndx) >1]
    overlapped_cut.sort()
    overlapped=list(overlapped_cut for overlapped_cut,_ in itertools.groupby(overlapped_cut))
    return overlapped

def check_atom_overlap(coordinates):
#### creates tree of atom coordinates
    tree = cKDTree(coordinates)
    overlapped = overlapping_atoms(tree)
#### runs through overlapping atoms and moves atom in a random diection until it is no longer overlapping
    while len(overlapped) > 0:
        for ndx_val, ndx in enumerate(overlapped):
            if len(overlapped) > 30:
                print('{:<130}'.format(''), end='\r')
                print('fixing '+str(len(overlapped))+' overlapped atoms: '+str(np.round((ndx_val/len(overlapped))*100,1))+' %', end='\r')
            xyz_check = np.array([coordinates[ndx[0]][0]+np.random.uniform(-0.2, 0.2), coordinates[ndx[0]][1]+np.random.uniform(-0.2, 0.2),coordinates[ndx[0]][2]+np.random.uniform(-0.2, 0.2)])
            while len(tree.query_ball_point(xyz_check, r=0.3)) > 1:
                xyz_check = np.array([coordinates[ndx[0]][0]+np.random.uniform(-0.2, 0.2), coordinates[ndx[0]][1]+np.random.uniform(-0.2, 0.2),coordinates[ndx[0]][2]+np.random.uniform(-0.2, 0.2)])
            coordinates[ndx[0]]=xyz_check
            tree = cKDTree(coordinates)
        overlapped = overlapping_atoms(tree)
    return coordinates


def fragment_location(residue):  
#### runs through dirctories looking for the atomistic fragments returns the correct location
    if residue in f_loc.p_residues:
        for directory in range(len(f_loc.p_directories)):
            if os.path.exists(f_loc.p_directories[directory][0]+residue+'/'+residue+'.pdb'):
                return f_loc.p_directories[directory][0]+residue+'/'+residue+'.pdb'
        for directory in range(len(f_loc.mod_directories)):
            if os.path.exists(f_loc.mod_directories[directory][0]+residue+'/'+residue+'.pdb'):
                return f_loc.mod_directories[directory][0]+residue+'/'+residue+'.pdb'
    else:
        for directory in range(len(f_loc.np_directories)):
            if os.path.exists(f_loc.np_directories[directory][0]+residue+'/'+residue+'.pdb'):
                return f_loc.np_directories[directory][0]+residue+'/'+residue+'.pdb'
    sys.exit('cannot find fragment: '+residue+'/'+residue+'.pdb')


def split_fragment_names(line, residue, resname):
    bead = gen.strip_header(line)
    group = f_loc.res_top[resname]['GROUPS'][bead]
    if group not in residue:
        residue[group] = {}
    if bead not in residue[group]:
        residue[group][bead]={} 
    return residue, group, bead   

def get_atomistic(frag_location):
#### read in atomistic fragments into dictionary    
    residue = {} ## a dictionary of bead in each residue eg residue[group][bead][atom number(1)][residue_name(ASP)/coordinates(coord)/atom name(C)/connectivity(2)/atom_mass(12)]
    fragment_mass = {}
    group=0
    resname = frag_location.split('/')[-1][:-4]
    with open(frag_location, 'r') as pdb_input:
        for line_nr, line in enumerate(pdb_input.readlines()):
            if line.startswith('['):
                residue, group, bead = split_fragment_names(line, residue, resname)
                fragment_mass[bead]=[]
            if line.startswith('ATOM'):
                line_sep = gen.pdbatom(line) ## splits up pdb line
                residue[group][bead][line_sep['atom_number']]={'coord':np.array([line_sep['x']*g_var.sf,line_sep['y']*g_var.sf,line_sep['z']*g_var.sf]),
                                                                'atom':line_sep['atom_name'],'resid':1, 'res_type':line_sep['residue_name'],
                                                                'frag_mass':1}    
#### updates fragment mass   
                if not gen.is_hydrogen(line_sep['atom_name']):
                    if line_sep['atom_name'] in f_loc.res_top[resname]['atom_masses']:
                        residue[group][bead][line_sep['atom_number']]['frag_mass']=f_loc.res_top[resname]['atom_masses'][line_sep['atom_name']]  ### updates atom masses with crude approximations
                        fragment_mass[bead].append([line_sep['x']*g_var.sf,line_sep['y']*g_var.sf,line_sep['z']*g_var.sf,f_loc.res_top[resname]['atom_masses'][line_sep['atom_name']]])               
                else:
                    fragment_mass[bead].append([line_sep['x']*g_var.sf,line_sep['y']*g_var.sf,line_sep['z']*g_var.sf,1])
    return residue, fragment_mass

def COM(mass, fragment):
#### returns center of mass of fragment
    
    try:
        if np.any(np.array(mass)[:,3]):      
            return np.average(np.array(mass)[:,:3], axis=0, weights=np.array(mass)[:,3])
        else:
            print('bead has no mass: \n')
            sys.exit(fragment)
    except:
        print(fragment, mass)
        if len(fragment) == 1:
            for key in fragment:
                if len(fragment[key]) == 1:
                    return fragment[key][1]['coord']
                else:
                    print('missing the mass one of the atoms\n')
                    print(mass)
                    sys.exit(fragment)
        else:
            print('missing the mass one of the atoms\n')
            print(mass)
            sys.exit(fragment)

def rigid_fit(group, frag_mass, resid, cg):
#### rigid fits group to CG beads
    rigid_mass_at = []
    rigid_mass_cg = []
    for bead in group:
        rigid_mass_at+=frag_mass[bead]
        rigid_mass_cg.append(cg[bead]['coord'])
    rigid_mass_at = COM(rigid_mass_at, group)
    rigid_mass_cg = np.mean(rigid_mass_cg, axis=0)
    
    group, COM_vector = align_at_frag_to_CG_frag(rigid_mass_at, rigid_mass_cg, group)
    at_frag_centers = {}
    cg_frag_centers = {}
    for bead in group:
        cg_frag_centers[bead] = cg[bead]['coord']
        at_frag_centers[bead] = COM(frag_mass[bead], bead)
        at_frag_centers[bead] -= COM_vector
    return rigid_mass_cg, at_frag_centers, cg_frag_centers, group


def align_at_frag_to_CG_frag(at_com, cg_com, group):
#### aligns atomistic fragment to cg bead
    COM_vector=at_com-cg_com ### gets vector between COM of atoms in fragment and cg bead 
    for bead in group:
        for atom in group[bead]: ### runs through atoms in fragments and centers on the cg bead 
            group[bead][atom]['coord']=group[bead][atom]['coord']-COM_vector
    return group, COM_vector


def connectivity(cg, at_frag_centers, cg_frag_centers, group, group_number):
#### returns the connections between the atomistic and coarse grain beads
    at_connection, cg_connection=[],[]
    for bead in cg:
        resname = cg[bead]['residue_name']
        break
    if resname in f_loc.sorted_connect:
        if len(f_loc.sorted_connect[resname]) > 0:
            for group_bead in group: 
                for bead_atom in group[group_bead]:
                    if bead_atom in f_loc.sorted_connect[resname][int(group_number)]:
                        for bead_connect in f_loc.sorted_connect[resname][int(group_number)][bead_atom]:
                            if 'cg_connect' not in locals():
                                cg_connect = [cg[bead_connect]['coord']]
                            else:
                                cg_connect.append(cg[bead_connect]['coord'])
                        at_connection.append(group[group_bead][bead_atom]['coord'])
                        cg_connection.append(np.mean(np.array(cg_connect), axis=0))    
    if len(at_frag_centers)  > 1:         
        for bead in at_frag_centers:
            cg_connection.append(cg_frag_centers[bead])
            at_connection.append(at_frag_centers[bead])     
    return at_connection, cg_connection    

def find_cross_vector(ca):
#### finds cross vector of the atoms CA, CA+1, CA+2 
    AB = ca[0]-ca[1]
    AC = ca[0]-ca[2] 
    cross_vector = np.cross(AB, AC)
    cross_vector = cross_vector/np.linalg.norm(cross_vector)
    return cross_vector

def noramlised_vector(c1, c2):
    initial_vector= c1-c2
    initial_vector = initial_vector/np.linalg.norm(initial_vector)
    return initial_vector


def align_to_vector(v1, v2):
#### returns the rotation matrix to rotate v1 to v2
    v = np.cross(v1,v2)
    c = np.dot(v1,v2)
    s = np.linalg.norm(v)

    rotation=np.array([[   0, -v[2],  v[1]],
                       [v[2],     0, -v[0]],
                       [-v[1], v[0],    0], 
                       ])

    r = np.identity(3) - rotation + np.matmul(rotation,rotation) * ((1 - c)/(s**2))
    return r

################################################################### Merged system

def merge_system_pdbs(protein):
    os.chdir(g_var.merged_directory)
#### create merged pdb 
    if not os.path.exists(g_var.merged_directory+'merged_cg2at'+protein+'.pdb'):
        pdb_output=gen.create_pdb(g_var.merged_directory+'merged_cg2at'+protein+'.pdb') 
        merge=[]
        merge_coords=[]
    #### run through every residue type in cg_residues
        for segment, residue_type in enumerate(g_var.cg_residues):
        #### if file contains user input identifier 
            if residue_type != 'PROTEIN':
                input_type=''
            else:
                input_type=protein
            merge, merge_coords = read_in_merged_pdbs(merge, merge_coords, g_var.working_dir+residue_type+'/'+residue_type+input_type+'_merged.pdb')
        if 'novo' in protein:
            print('checking for atom overlap in : '+protein[1:])
            merge_coords = check_atom_overlap(merge_coords)
        for line_val, line in enumerate(merge):
            x, y, z = gen.trunc_coord(merge_coords[line_val])
            pdb_output.write(g_var.pdbline%((int(line['atom_number']), line['atom_name'], line['residue_name'],' ',line['residue_id'],\
                x,y,z,1,0))+'\n')
        pdb_output.write('TER\nENDMDL')
        pdb_output.close()

def read_in_merged_pdbs(merge, merge_coords, location):
    if os.path.exists(location):
    #### opens pdb files and writes straight to merged_cg2at pdb
        with open(location, 'r') as pdb_input:
            for line in pdb_input.readlines():
                if line.startswith('ATOM'):
                    line_sep = gen.pdbatom(line)
                    line_sep['x'],line_sep['y'],line_sep['z'] = gen.trunc_coord([line_sep['x'],line_sep['y'],line_sep['z']])
                    merge.append(line_sep)
                    merge_coords.append([line_sep['x'],line_sep['y'],line_sep['z']])
        return merge, merge_coords
    else:
        sys.exit('cannot find minimised residue: \n'+ location) 

def check_overlap_chain(chain, input, collate=False):
    checked=[]
    if not os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN_'+input+str(chain)+'_gmx_checked.pdb'):
        lines, coords = read_in_merged_pdbs([], [], 'PROTEIN_'+input+str(chain)+'_gmx.pdb')
        coords = at_mod_p.correct_amide_h(lines, coords)
        updated_coords = check_atom_overlap(coords)
        pdb_output=gen.create_pdb(g_var.working_dir+'PROTEIN/PROTEIN_'+input+str(chain)+'_gmx_checked.pdb') 

        for line_val, line in enumerate(lines):
            x, y, z = gen.trunc_coord(updated_coords[line_val])
            pdb_output.write(g_var.pdbline%((int(line['atom_number']), line['atom_name'], line['residue_name'],' ',line['residue_id'],\
                            x,y,z,1,0))+'\n')
            if collate:
                checked.append(g_var.pdbline%((int(line['atom_number']), line['atom_name'], line['residue_name'],' ',line['residue_id'],\
                            x,y,z,1,0))+'\n')
    if collate:
        return checked

def fetch_chiral_coord(merge_temp, residue_type):
    chiral_atoms={}
    coord=[]
    for atom in range(len(merge_temp)):
        if residue_type == 'PROTEIN':
            resname = merge_temp[atom]['residue_name']
        else:
            resname=residue_type
        if len(f_loc.res_top[resname]['CHIRAL']['atoms']) > 0:
            if merge_temp[atom]['atom_name'] in f_loc.res_top[resname]['CHIRAL']['atoms']:
                if merge_temp[atom]['residue_id'] not in chiral_atoms:
                    chiral_atoms[merge_temp[atom]['residue_id']]={}
                chiral_atoms[merge_temp[atom]['residue_id']][merge_temp[atom]['atom_name']]=atom
        coord.append(np.array([merge_temp[atom]['x'],merge_temp[atom]['y'],merge_temp[atom]['z']]))
    return chiral_atoms, coord

def fix_chirality(merge, merge_temp, merged_coords, residue_type):
#### fixes chiral groups
    r_b_vec, r_b_inv = read_in.real_box_vectors(g_var.box_vec)
    chiral_atoms, coord= fetch_chiral_coord(merge_temp, residue_type)
    for residue in chiral_atoms:
        if residue_type == 'PROTEIN':
            for atom in chiral_atoms[residue]:
                resname = merge_temp[chiral_atoms[residue][atom]]['residue_name']
                resid = merge_temp[chiral_atoms[residue][atom]]['residue_id']
                break
        else:
            resname=residue_type
        for chiral_group in f_loc.res_top[resname]['CHIRAL']:

            if chiral_group != 'atoms':
                stat = merge_temp[chiral_atoms[residue][chiral_group]].copy()
                atom_move = {'stat':np.array([stat['x'],stat['y'],stat['z']]), 'm':'', 'c1':'', 'c2':'', 'c3':''}
                for chir_atom in atom_move:
                    if chir_atom != 'stat':
                        test = merge_temp[chiral_atoms[residue][f_loc.res_top[resname]['CHIRAL'][chiral_group][chir_atom]]].copy()
                        atom_move[chir_atom]= np.array([test['x'],test['y'],test['z']])
                        if gen.calculate_distance(atom_move['stat'], atom_move[chir_atom]) > 10:
                            atom_move[chir_atom] = np.array(read_in.brute_mic(atom_move['stat'],atom_move[chir_atom], r_b_vec))
                S_M = atom_move['m'] - atom_move['stat']
                rotation = align_to_vector(S_M, [0,0,1])
                c1_coord = (atom_move['c1'] - atom_move['stat']).dot(rotation)
                c2_coord = (atom_move['c2'] - atom_move['stat']).dot(rotation)
                c3_coord = (atom_move['c3'] - atom_move['stat']).dot(rotation)

                C1_C2_a = gen.angle_clockwise(c1_coord[0:2], c2_coord[0:2])
                C1_C3_a = gen.angle_clockwise(c1_coord[0:2], c3_coord[0:2])

                if C1_C2_a > C1_C3_a:
                    for ax_val, axis in enumerate(['x', 'y', 'z']):
                        merge_temp[chiral_atoms[residue][f_loc.res_top[resname]['CHIRAL'][chiral_group]['m']]][axis] = merge_temp[chiral_atoms[residue][f_loc.res_top[resname]['CHIRAL'][chiral_group]['m']]][axis] - (3*S_M[ax_val])   
                        merge_temp[chiral_atoms[residue][chiral_group]][axis] = merge_temp[chiral_atoms[residue][chiral_group]][axis] - (S_M[ax_val])        
                    coord[chiral_atoms[residue][f_loc.res_top[resname]['CHIRAL'][chiral_group]['m']]] -=  (2*S_M) #move_coord -
                    coord[chiral_atoms[residue][chiral_group]] -=  (0.25*S_M) #stat_coord -
                
    merge+=merge_temp
    merged_coords+=coord
    return merge, merged_coords

def check_hydrogens(residue):
#### finds the connecting carbons and their associated carbons [carbon atom, hydrogen ref number, connecting ref number]    
    for atom_num, atom in enumerate(residue):
        resname=residue[atom]['res_type']
        break

    for atom in f_loc.hydrogen[resname]:
        h_coord = []
        for group in f_loc.sorted_connect[resname]:
            if atom in f_loc.sorted_connect[resname][group]:
                for hydrogen in f_loc.hydrogen[resname][atom]:
                    h_coord.append(residue[hydrogen]['coord'])
                h_com=np.mean(np.array(h_coord), axis=0)
                for heavy_bond in f_loc.heavy_bond[resname][atom]:
                    for group_check in f_loc.sorted_connect[resname]:
                        if heavy_bond in f_loc.sorted_connect[resname][group_check] and group_check != group:
                            con_heavy_atom = heavy_bond
                            con_heavy_atom_co =  residue[con_heavy_atom]['coord']
            #### vector between H COM and bonded carbon 
                vector=np.array([h_com[0]-residue[atom]['coord'][0],h_com[1]-residue[atom]['coord'][1],h_com[2]-residue[atom]['coord'][2]])
                h_com_f=h_com+vector*2
                d1 = gen.calculate_distance(h_com, con_heavy_atom_co)    
                d2 = gen.calculate_distance(h_com_f, con_heavy_atom_co)   
                if d2 < d1:
                    for h_at in f_loc.hydrogen[resname][atom]:
                        residue[h_at]['coord']=residue[h_at]['coord']-vector*2
    return residue

def check_ringed_lipids(protein):
    box_vec = g_var.box_vec.split()[1:4]
    os.chdir(g_var.merged_directory)
    merge, merge_coords = read_in_merged_pdbs([], [], protein)
    resid_prev=0
    ringed=[]
    # offset =0
    # and atom['atom_number']-offset > max(f_loc.hydrogen[atom['residue_name']], key=f_loc.hydrogen[atom['residue_name']].get):# and atom['residue_name'] n:

    for at_val, atom in enumerate(merge):         
        if atom['residue_name'] in f_loc.np_residues:
            if atom['residue_id'] != resid_prev:
                if 'offset' in locals() and len(f_loc.heavy_bond[resname])>0: 
                    if at_val-offset >= max(f_loc.heavy_bond[resname]): 
                        offset=at_val
                else:
                    offset=at_val
                    resname = atom['residue_name']
            resid_prev=atom['residue_id']
            if atom['atom_number']-offset in f_loc.heavy_bond[atom['residue_name']]:
                resname = atom['residue_name']
                at_coord = [atom['x'], atom['y'], atom['z']]
                for at_bond in f_loc.heavy_bond[atom['residue_name']][atom['atom_number']-offset]:
                    at_bond_coord = [merge[at_bond+offset-1]['x'], merge[at_bond+offset-1]['y'], merge[at_bond+offset-1]['z']]
                    for xyz in range(3):
                    #### for x, y, z if the distance between bead is more than half the box length
                        if gen.calculate_distance(at_coord, at_bond_coord) > float(box_vec[xyz])/2:
                        #### if the bead if in the opposite 1/3 of the box the position the box length is add/subtracted
                            if at_bond_coord[xyz] <= float(box_vec[xyz])/2:
                                temp = at_bond_coord[xyz]+float(box_vec[xyz])
                            elif at_bond_coord[xyz] > float(box_vec[xyz])/2:
                                temp = at_bond_coord[xyz]-float(box_vec[xyz])
                        #### if distance between corrected coordinate is still > 1/2 the box length then counts as a new chain
                            if np.sqrt((temp-at_coord[xyz])**2) < 2:
                                at_bond_coord[xyz] = temp
                    dist = gen.calculate_distance(at_coord, at_bond_coord)
                    if dist > 2:
                         ringed.append([at_val, at_bond+offset-1])

    return ringed

