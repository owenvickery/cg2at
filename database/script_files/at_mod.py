#!/usr/bin/env python3

import os, sys
import numpy as np
from scipy.spatial import KDTree
import gen, g_var, f_loc

def rotate_atom(coord, center,xyz_rot_apply):
    coord =  coord-center  #### centers COM coordinates to 0,0,0
    coord =  coord.dot(gen.eulerAnglesToRotationMatrix([xyz_rot_apply[0],0,0]))  #### rotates coord around x
    coord =  coord.dot(gen.eulerAnglesToRotationMatrix([0,xyz_rot_apply[1],0]))  #### rotates coord around y
    coord =  coord.dot(gen.eulerAnglesToRotationMatrix([0,0,xyz_rot_apply[2]]))  #### rotates coord around z
    coord =  coord+center #### translates coord back by original offset
    return coord

def rotate(at_connections, cg_connections, same):
    xyz_rot_apply=[]
#### iterates through rotation matrices
    for xyz_rot in [f_loc.x_rot,f_loc.y_rot,f_loc.z_rot]:
        dist=[]
    #### iterates through rotation matrices 
        for rot_val, rotation in enumerate(xyz_rot):
        #### applies matrix to coordinates saved as check
            check = at_connections.dot(rotation)
        #### for each connection the distance is calculated and added to list
            individual_connections=[]
            for connect in range(len(cg_connections)):
                individual_connections.append(np.sqrt(((check[connect][0]-cg_connections[connect][0])**2)+((check[connect][1]-cg_connections[connect][1])**2)+((check[connect][2]-cg_connections[connect][2])**2)))
        #### for each rotation the connection distances are added to dist list 
            dist.append(individual_connections)
    #### the RMS is calculated for each rotation    
        dist=np.array(dist)
        inter= np.sqrt(np.mean(dist**2,axis=1))
        if len(dist[0])==2 and same:
            ratio=dist/np.min(dist, axis=1)[:,np.newaxis]
            rotation_index=np.argmin(inter[np.where(np.sum(ratio, axis=1)<np.min(np.sum(ratio, axis=1))*1.02)])
            for i in range(len(ratio)):
                if np.all(ratio[i]<1.05):
                    if 'rotation_RMS' in locals():
                        if inter[i] < rotation_RMS:
                            rotation_RMS = inter[i]
                            rotation_index=i
                    else:
                        rotation_RMS=inter[i]
                        rotation_index=i
        else:
            rotation_index=np.argmin(inter)

    #### the rotation with the lowest RMS applied to the at_connections
        at_connections = at_connections.dot(xyz_rot[rotation_index])
    #### the optimal rotation is added to xyz_rot_apply list as radians
        xyz_rot_apply.append(np.radians(rotation_index*5))
    return xyz_rot_apply

def add_to_sequence(sequence, residue, chain_count):
    if residue  not in f_loc.mod_residues:
        sequence[chain_count]+=g_var.aas[residue]
    else:
        sequence[chain_count]+='X'    
    return sequence

def check_atom_overlap(coordinates):
#### creates tree of atom coordinates
    tree = KDTree(coordinates)
#### provides index of any atoms that are within 0.3A of each other
    overlapped_ndx = tree.query_ball_tree(tree, r=0.3)  ### takes a while 
    done=[]
    moved_coord=[]
    dist=0.35
#### runs through overlapping atoms and moves atom in a random diection until it is no longer overlapping
    for ndx in overlapped_ndx:
        if len(ndx) == 2 and ndx[0] not in done:
            xyz_check = np.array([coordinates[ndx[0]][0]+np.random.uniform(-0.2, 0.2), coordinates[ndx[0]][1]+np.random.uniform(-0.2, 0.2),coordinates[ndx[0]][2]+np.random.uniform(-0.2, 0.2)])
            if len(moved_coord)>0:
                dist=np.min(np.sqrt(((xyz_check[0]-np.array(moved_coord)[:,0])**2)+((xyz_check[1]-np.array(moved_coord)[:,1])**2)+((xyz_check[2]-np.array(moved_coord)[:,2])**2)))
            while len(tree.query_ball_point(xyz_check, r=0.3)) == 2 or dist < 0.3:
                if len(moved_coord)>0:
                    dist=np.min(np.sqrt(((xyz_check[0]-np.array(moved_coord)[:,0])**2)+((xyz_check[1]-np.array(moved_coord)[:,1])**2)+((xyz_check[2]-np.array(moved_coord)[:,2])**2)))
                xyz_check = np.array([coordinates[ndx[0]][0]+np.random.uniform(-0.2, 0.2), coordinates[ndx[0]][1]+np.random.uniform(-0.2, 0.2),coordinates[ndx[0]][2]+np.random.uniform(-0.2, 0.2)])
            coordinates[ndx[0]]=xyz_check
            moved_coord.append(xyz_check)
            done.append(ndx[0])
    return coordinates


def fragment_location(residue, fragment):  
#### runs through dirctories looking for the atomistic fragments returns the correct location
    if residue in f_loc.p_residues:
        for directory in range(len(f_loc.p_directories)):
            if os.path.exists(f_loc.p_directories[directory][0]+residue+'/'+fragment):
                return f_loc.p_directories[directory][0]+residue+'/'+fragment
        for directory in range(len(f_loc.mod_directories)):
            if os.path.exists(f_loc.mod_directories[directory][0]+residue+'/'+fragment):
                return f_loc.mod_directories[directory][0]+residue+'/'+fragment
    else:
        for directory in range(len(f_loc.np_directories)):
            if os.path.exists(f_loc.np_directories[directory][0]+residue+'/'+fragment):
                return f_loc.np_directories[directory][0]+residue+'/'+fragment
    sys.exit('cannot find fragment: '+residue+'/'+fragment)

def get_atomistic(residue,cg_fragment, cg_coord,resid):
    SF=1   #### scaling factor for fragments
#### find atomistic residues
    residue_list={} ## a dictionary of bead in each residue eg residue_list[atom number(1)][residue_name(ASP)/coordinates(coord)/atom name(C)/connectivity(2)/atom_mass(12)]
    frag_location=fragment_location(residue, cg_fragment+'.pdb') ### get fragment location from database
    fragment_masses=[] ### list [[coord, mass],[coord, mass]]
#### read in atomistic fragments into dictionary    
    with open(frag_location, 'r') as pdb_input:
        for line_nr, line in enumerate(pdb_input.readlines()):
            if line.startswith('ATOM'):
                line_sep = gen.pdbatom(line) ## splits up pdb line
                residue_list[line_sep['atom_number']]={'coord':np.array([line_sep['x']*SF,line_sep['y']*SF,line_sep['z']*SF]),'atom':line_sep['atom_name'], 'res_type':line_sep['residue_name'],'extra':line_sep['backbone'], 'connect':line_sep['connect'], 'frag_mass':1}
#### updates fragment mass   
                if 'H' not in line_sep['atom_name']:
                    for atom in line_sep['atom_name']:
                        if atom in g_var.mass:
                            residue_list[line_sep['atom_number']]['frag_mass']=g_var.mass[atom]  ### updates atom masses with crude approximations
                            fragment_masses.append([line_sep['x']*SF,line_sep['y']*SF,line_sep['z']*SF,g_var.mass[atom]])
                else:
                    fragment_masses.append([line_sep['x']*SF,line_sep['y']*SF,line_sep['z']*SF,1])
    return align_at_frag_to_CG_frag(fragment_masses, cg_coord, residue_list)

def align_at_frag_to_CG_frag(fragment_masses, cg_coord, residue_list):
#### aligns atomistic fragment to cg bead
    COM_vector=np.average(np.array(fragment_masses)[:,:3], axis=0, weights=np.array(fragment_masses)[:,3])-np.array(cg_coord['coord']) ### gets vector between COM of atoms in fragment and cg bead 
    for at_id, residue in enumerate(residue_list): ### runs through atoms in fragments and centers on the cg bead 
        residue_list[residue]['coord']=residue_list[residue]['coord']-COM_vector
    return residue_list


def get_atomistic_fragments(cg_residue_type,cg_residue, cg_resid):
    at_residues={}
    connect=[]
#### runs through every in bead in residue 
    for cg_bead in cg_residue:
    #### gets atoms from database for each bead 
        at_residues[cg_bead]=get_atomistic(cg_residue_type,cg_bead, cg_residue[cg_bead], cg_resid+1)
    #### if not SOL/ION the connectivity is read from the fragment dictionary key (connect)
        if cg_residue_type not in ['SOL', 'ION']:
            for atom_num, atom in enumerate(at_residues[cg_bead]):
            #### if atom has a connection which is not zero (0 = does not connect)
                if at_residues[cg_bead][atom]['connect'] > 0:
                    connect.append([cg_bead,atom, at_residues[cg_bead][atom]['connect']]) 
    connect=np.array(connect)   
    return at_residues, connect


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


def find_cross_vector(ca, C, O):
    initial_vector= O-C
    initial_vector = initial_vector/np.linalg.norm(initial_vector)
    AB = ca[0]-ca[1]
    AC = ca[0]-ca[2] 
    cross_vector = np.cross(AB, AC)
    cross_vector = cross_vector/np.linalg.norm(cross_vector)
    
    return initial_vector, cross_vector


def align_to_vector(v1, v2):
    v = np.cross(v1,v2)
    c = np.dot(v1,v2)
    s = np.linalg.norm(v)

    rotation=np.array([[0,    -v[2],  v[1]],
                       [v[2],     0, -v[0]],
                       [-v[1], v[0],    0], 
                       ])

    r = np.identity(3) - rotation + np.matmul(rotation,rotation) * ((1 - c)/(s**2))
    return r

################################################################### Merged system

def merge_system_pdbs(system, protein, cg_residues, box_vec):
    os.chdir(g_var.working_dir+'MERGED')
#### create merged pdb 
    pdb_output=gen.create_pdb(g_var.working_dir+'MERGED/merged_cg2at'+protein+'.pdb', box_vec) 
    merge=[]
    merge_coords=[]
#### run through every residue type in cg_residues
    for segment, residue_type in enumerate(cg_residues):
    #### if file contains user input identifier 
        if residue_type != 'PROTEIN':
            input_type=''
        else:
            input_type=protein
        if os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+input_type+'_merged.pdb'):
        #### opens pdb files and writes straight to merged_cg2at pdb
            with open(g_var.working_dir+residue_type+'/'+residue_type+input_type+'_merged.pdb', 'r') as pdb_input:
                for line in pdb_input.readlines():
                    if line.startswith('ATOM'):
                        line_sep = gen.pdbatom(line)
                        merge.append(line_sep)
                        merge_coords.append([line_sep['x'],line_sep['y'],line_sep['z']])
        else:
            sys.exit('cannot find minimised residue: \n'+ g_var.working_dir+residue_type+'/'+residue_type+input_type+'_merged.pdb')
    if protein in ['_at_rep_user_supplied','_novo']:
        print('checking for atom overlap in : '+protein[1:])
        merge_coords = check_atom_overlap(merge_coords)
    for line_val, line in enumerate(merge):
        pdb_output.write(g_var.pdbline%((int(line['atom_number']), line['atom_name'], line['residue_name'],' ',line['residue_id'],\
            merge_coords[line_val][0],merge_coords[line_val][1],merge_coords[line_val][2],1,0))+'\n')
    pdb_output.write('TER\nENDMDL')
    pdb_output.close()

def fix_chirality(merge, merge_temp, merged_coords):
    chiral_atoms={}
    coord=[]
    for atom in range(len(merge_temp)):
        if merge_temp[atom]['residue_name'] in f_loc.chiral:
            if merge_temp[atom]['residue_id'] not in chiral_atoms:
                chiral_atoms[merge_temp[atom]['residue_id']]={}
            if merge_temp[atom]['atom_name'] in f_loc.chiral[merge_temp[atom]['residue_name']]['atoms']:
                chiral_atoms[merge_temp[atom]['residue_id']][merge_temp[atom]['atom_name']]=atom
        coord.append(np.array([merge_temp[atom]['x'],merge_temp[atom]['y'],merge_temp[atom]['z']]))

    for residue in chiral_atoms:
        for atom in chiral_atoms[residue]:
            resname = merge_temp[chiral_atoms[residue][atom]]['residue_name']
            break
        for chiral_group in f_loc.chiral[resname]:
            if chiral_group != 'atoms':
                stat = merge_temp[chiral_atoms[residue][chiral_group]]
                move = merge_temp[chiral_atoms[residue][f_loc.chiral[resname][chiral_group]['m']]]
                c1   = merge_temp[chiral_atoms[residue][f_loc.chiral[resname][chiral_group]['c1']]]
                c2   = merge_temp[chiral_atoms[residue][f_loc.chiral[resname][chiral_group]['c2']]]
                c3   = merge_temp[chiral_atoms[residue][f_loc.chiral[resname][chiral_group]['c3']]]

                stat_coord = np.array([stat['x'], stat['y'], stat['z']])
                move_coord = np.array([move['x'], move['y'], move['z']])
                c1_coord   = np.array([c1['x'], c1['y'], c1['z']])
                c2_coord   = np.array([c2['x'], c2['y'], c2['z']])
                c3_coord   = np.array([c3['x'], c3['y'], c3['z']])

                S_M=move_coord-stat_coord
                rotation = align_to_vector(S_M, [0,0,1])
                center = stat_coord

                c1_coord = (c1_coord - center).dot(rotation)
                c2_coord = (c2_coord - center).dot(rotation)
                c3_coord = (c3_coord - center).dot(rotation)

                C1_C2_a = gen.angle_clockwise(c1_coord[0:2], c2_coord[0:2])
                C1_C3_a = gen.angle_clockwise(c1_coord[0:2], c3_coord[0:2])
                if C1_C2_a > C1_C3_a:
                    for ax_val, axis in enumerate(['x', 'y', 'z']):
                        merge_temp[chiral_atoms[residue][f_loc.chiral[resname][chiral_group]['m']]][axis] = merge_temp[chiral_atoms[residue][f_loc.chiral[resname][chiral_group]['m']]][axis] - (3*S_M[ax_val])   
                        merge_temp[chiral_atoms[residue][chiral_group]][axis] = merge_temp[chiral_atoms[residue][chiral_group]][axis] - (S_M[ax_val])        
                    coord[chiral_atoms[residue][f_loc.chiral[resname][chiral_group]['m']]] = move_coord - (3*S_M) 
                    coord[chiral_atoms[residue][chiral_group]] = stat_coord - (1.5*S_M) 
    merge+=merge_temp
    merged_coords+=coord

    return merge, merged_coords