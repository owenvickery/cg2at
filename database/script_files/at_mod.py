#!/usr/bin/env python3

import os, sys
import numpy as np
from scipy.spatial import KDTree
import gen, g_var, f_loc

def sanity_check_fragments(res, cg, sin_bead):
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
                if bead ==sin_bead:
                    for atom in residue[group][bead]:
                        atom_list.append(atom)
    return bead_list, atom_list

def sanity_check_atoms(atom_list, res):
    for at_num in range(1, len(atom_list)+1):
        if at_num not in atom_list:
            sys.exit('atom number '+str(at_num)+' is missing from fragment library: '+res)
def sanity_check_beads(bead_list, cg, res):
    for bead in cg:
        if bead not in bead_list:
            sys.exit('The bead '+bead+' is missing from the fragment library: '+res)   

def sanity_check(cg_residues):
    for res_type in cg_residues:
        if res_type == 'PROTEIN':
            for p_res in cg_residues['PROTEIN']:
                for bead in cg_residues['PROTEIN'][p_res]:
                    resname = cg_residues['PROTEIN'][p_res][bead]['residue_name']
                    break
                bead_list, atom_list = sanity_check_fragments(resname, cg_residues['PROTEIN'][p_res], False)
                sanity_check_beads(bead_list, cg_residues['PROTEIN'][p_res], resname)  
                sanity_check_atoms(atom_list, resname)
        elif res_type in ['SOL', 'ION']:
            for p_res in cg_residues[res_type]:
                for bead in cg_residues[res_type][p_res]:
                    sin_bead=bead
                    break
                bead_list, atom_list = sanity_check_fragments(res_type, cg_residues[res_type][p_res], sin_bead )
                sanity_check_beads(bead_list, cg_residues[res_type][p_res], res_type) 
                sanity_check_atoms(atom_list, res_type)
        else:
            bead_list, atom_list = sanity_check_fragments(res_type, cg_residues[res_type], False)
            for residue in cg_residues[res_type]:
                sanity_check_beads(bead_list, cg_residues[res_type][residue], res_type) 
                sanity_check_atoms(atom_list, res_type)


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

def get_atomistic(frag_location):
    SF=1   #### scaling factor for fragments
#### find atomistic residues
    
#### read in atomistic fragments into dictionary    
    residue = {} ## a dictionary of bead in each residue eg residue[group][bead][atom number(1)][residue_name(ASP)/coordinates(coord)/atom name(C)/connectivity(2)/atom_mass(12)]
    fragment_mass = {}
    group=0
    with open(frag_location, 'r') as pdb_input:
        for line_nr, line in enumerate(pdb_input.readlines()):
            if line.startswith('['):
                line_split=line.split()
                if not g_var.ind:
                    bead, group = line_split[1], line_split[2]
                else:
                    bead = line_split[1]
                    group+=1
                fragment_mass[bead]=[]
                if group != ']' and group not in residue:
                    residue[group] = {}
                if bead not in residue[group]:
                    residue[group][bead]={}
            if line.startswith('ATOM'):
                line_sep = gen.pdbatom(line) ## splits up pdb line
                residue[group][bead][line_sep['atom_number']]={'coord':np.array([line_sep['x']*SF,line_sep['y']*SF,line_sep['z']*SF]),'atom':line_sep['atom_name'],'resid':line_sep['residue_id'], 'res_type':line_sep['residue_name'],'extra':line_sep['backbone'], 'connect':line_sep['connect'], 'frag_mass':1}
#### updates fragment mass   
                if 'H' not in line_sep['atom_name']:
                    for atom in line_sep['atom_name']:
                        if atom in g_var.mass:
                            residue[group][bead][line_sep['atom_number']]['frag_mass']=g_var.mass[atom]  ### updates atom masses with crude approximations
                            fragment_mass[bead].append([line_sep['x']*SF,line_sep['y']*SF,line_sep['z']*SF,g_var.mass[atom]])
                else:
                    fragment_mass[bead].append([line_sep['x']*SF,line_sep['y']*SF,line_sep['z']*SF,1])
    return residue, fragment_mass

def COM(mass, fragment):
    if np.any(np.array(mass)[:,3]):      
        mass = np.average(np.array(mass)[:,:3], axis=0, weights=np.array(mass)[:,3])
    else:
        print('bead has no mass: \n')
        sys.exit(fragment)
    return mass

def rigid_fit(group, frag_mass, resid, cg):
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


def connection(residue):
    connect=[]
#### runs through every in bead in residue 
    for group in residue:
        for bead in residue[group]:
            for atom_num, atom in enumerate(residue[group][bead]):
            #### if atom has a connection which is not zero (0 = does not connect)
                if residue[group][bead][atom]['connect'] > 0:
                    connect.append([group, bead, atom, residue[group][bead][atom]['connect']]) 
    connect=np.array(connect)   
    return connect

def connectivity(cg, connect, at_frag_centers, cg_frag_centers, group, group_number):
    at_connection, cg_connection=[],[]
    if len(connect) > 0:
        for info in connect[np.where(str(group_number) == connect[:,0])]:
            for ind in connect:
                if info[3] == ind[3] and group_number != ind[0]: 
                    try:
                        cg_connection.append(cg[ind[1]]['coord'])
                    except:
                        sys.exit('cannot find CG bead: '+ind[1])
                    try:
                        at_connection.append(group[info[1]][int(info[2])]['coord'])  
                    except:
                        sys.exit('cannot find atomistic bead: '+info[1])
        if len(at_frag_centers)  > 1:         
            for bead in at_frag_centers:
                cg_connection.append(cg_frag_centers[bead])
                at_connection.append(at_frag_centers[bead])          
    return at_connection, cg_connection    

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