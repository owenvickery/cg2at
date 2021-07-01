#!/usr/bin/env python3

import os, sys
import numpy as np
import itertools
import math
from scipy.spatial import cKDTree
import gen, g_var, at_mod_p, read_in, gro


### sanity checking

def sanity_check_fragments(res, cg, sin_bead):
#### fetches bead and atom info from fragment database
    location = gen.fragment_location(res)
    residue, fragment_mass = get_atomistic(location, res)
    atom_list = []
    bead_list = []
    for group in residue.values():
        for bead in group:
            bead_list.append(bead)
            if not sin_bead:
                for atom in group[bead]:
                    atom_list.append(atom)
            else:
                if bead == sin_bead:
                    for atom in group[bead]:
                        atom_list.append(atom)
    return sorted(bead_list), atom_list

def sanity_check_atoms(atom_list, res):
#### checks atom order
    for at_num in range(1, len(atom_list)+1):
        if at_num not in atom_list:
            sys.exit('atom number '+str(at_num)+' is missing from fragment library: '+res+'\n')

def sanity_check_beads(bead_list, cg, res):
#### checks if bead is in fragment library
    bead_list_new=[]
    for bead in cg:
        if bead not in bead_list:
            new_bead = bead[1:]+bead[0]
            if new_bead in bead_list and new_bead not in cg:
                pass
            else:
                sys.exit('The bead '+bead+' is missing from the fragment library: '+res+'\n')   
        bead_list_new.append(bead)
    return bead_list_new

def sanity_check_protein_other(res_type, test=False):
    bead_list, atom_list = {},{}
    for residue in g_var.cg_residues[res_type]:
        for bead in g_var.cg_residues[res_type][residue]:
            resname = g_var.cg_residues[res_type][residue][bead]['residue_name']
            break
        if resname not in bead_list or resname not in atom_list:
            bead_list[resname], atom_list[resname] = sanity_check_fragments(resname, g_var.cg_residues[res_type][residue], False)
            sanity_check_atoms(atom_list[resname], resname)
        bead_list_cg = sanity_check_beads(bead_list[resname], g_var.cg_residues[res_type][residue], resname)  
        if bead_list[resname] != sorted(bead_list_cg):
            if len(bead_list[resname]) == len(bead_list_cg):
                fix_atom_wrap(bead_list[resname], bead_list_cg, res_type, residue)
            else:
                if not test:
                    print('There is a issue with residue: '+resname+' '+str(residue+1)+'. If expected ignore this message.' )
                if res_type != 'OTHER':
                    sys.exit('number of atomistic fragments: '+str(len(bead_list[resname]))+' does not equal number of CG beads: '+str(len(bead_list_cg)))

def sanity_check_solvent(res_type):
    bead_list, atom_list = {},{}
    for residue in g_var.cg_residues[res_type]:
        for bead in g_var.cg_residues[res_type][residue]:
            sin_bead=bead
            break
        if sin_bead not in bead_list or sin_bead not in atom_list:
            bead_list[sin_bead], atom_list[sin_bead] = sanity_check_fragments(res_type, g_var.cg_residues[res_type][residue], sin_bead)
            sanity_check_atoms(atom_list[sin_bead], res_type)
        sanity_check_beads(bead_list[sin_bead], g_var.cg_residues[res_type][residue], res_type)

def sanity_check_non_protein(res_type):
    bead_list, atom_list = {},{}
    if res_type not in bead_list or res_type not in atom_list:
        bead_list[res_type], atom_list[res_type] = sanity_check_fragments(res_type, g_var.cg_residues[res_type], False)
    sanity_check_atoms(atom_list[res_type], res_type)
    for residue in g_var.cg_residues[res_type]:
        bead_list_cg = sanity_check_beads(bead_list[res_type], g_var.cg_residues[res_type][residue], res_type)
        if bead_list[res_type] != sorted(bead_list_cg):
            if len(bead_list[res_type]) == len(bead_list_cg):
                fix_atom_wrap(bead_list[res_type], bead_list_cg, res_type, residue)
            else:
                print('There is a issue with residue: '+res_type+' '+str(residue+1))
                sys.exit('number of atomistic fragments: '+str(len(bead_list[res_type]))+' does not equal number of CG beads: '+str(len(bead_list_cg)))

def sanity_check():
#### runs through every bead and checks whether it exists
    for res_type in g_var.cg_residues:
        if res_type in ['PROTEIN', 'OTHER']:
            sanity_check_protein_other(res_type)
        elif res_type in g_var.sol_residues or res_type in g_var.ion_residues:
            sanity_check_solvent(res_type)
        else:
            sanity_check_non_protein(res_type)

def fix_atom_wrap(bead_list_frag, bead_list_cg, section, resid):
    for bead in bead_list_cg:
        if bead not in bead_list_frag:
            new_bead = bead[1:]+bead[0]
            if new_bead in bead_list_frag and new_bead not in bead_list_cg:
                g_var.cg_residues[section][resid][new_bead] = g_var.cg_residues[section][resid][bead]
                del g_var.cg_residues[section][resid][bead]
            else:
                print('There is a issue with residue: '+section+' '+str(resid+1))
                print('input file list:\n',sorted(bead_list_cg))
                print('\ncannot find: '+bead+' or '+new_bead+' in fragment list:')
                sys.exit(sorted(bead_list_frag))

#####  Sanity check end

#### rotations section

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

## rotations end
## alignments 

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

def align_at_frag_to_CG_frag(at_com, cg_com, group):
#### aligns atomistic fragment to cg bead
    COM_vector=at_com-cg_com ### gets vector between COM of atoms in fragment and cg bead 
    for bead in group:
        for atom in group[bead]: ### runs through atoms in fragments and centers on the cg bead 
            group[bead][atom]['coord']=group[bead][atom]['coord']-COM_vector
    return group, COM_vector

def COM(mass, fragment):
#### returns center of mass of fragment
    try:
        if np.any(np.array(mass)[:,3]):   
            return np.average(np.array(mass)[:,:3], axis=0, weights=np.array(mass)[:,3])
    except:
        for bead in fragment:
            print(bead, fragment[bead], '\n')
            print(mass)
            sys.exit('missing the mass one of the atoms in '+fragment[bead][1]['res_type'])      

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
    rigid_mass_cg, at_frag_centers, cg_frag_centers, group
    return rigid_mass_cg, at_frag_centers, cg_frag_centers, group

## alignment end

## overlap checker

def overlapping_atoms(tree):
    overlapped_ndx = tree.query_ball_tree(tree, r=g_var.args.ov)
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
            if np.round((ndx_val/len(overlapped))*100,2).is_integer() and len(overlapped) > 30:
                print('fixing '+str(len(overlapped))+' overlapped atoms: '+str(np.round((ndx_val/len(overlapped))*100,2))+' %', end='\r')
            xyz_check = np.array([coordinates[ndx[0]][0]+np.random.uniform(-0.2, 0.2), coordinates[ndx[0]][1]+np.random.uniform(-0.2, 0.2),coordinates[ndx[0]][2]+np.random.uniform(-0.2, 0.2)])
            while len(tree.query_ball_point(xyz_check, r=g_var.args.ov)) > 1:
                xyz_check = np.array([coordinates[ndx[0]][0]+np.random.uniform(-0.2, 0.2), coordinates[ndx[0]][1]+np.random.uniform(-0.2, 0.2),coordinates[ndx[0]][2]+np.random.uniform(-0.2, 0.2)])
            coordinates[ndx[0]]=xyz_check
            tree = cKDTree(coordinates)
        if len(overlapped) > 30:
            print('{:<100}'.format(''), end='\r')
        overlapped = overlapping_atoms(tree)
    return coordinates

## overlap end
## get fragment information

def split_fragment_names(line, residue, resname):
    bead = gen.strip_header(line)
    group = g_var.res_top[resname]['GROUPS'][bead]
    if group not in residue:
        residue[group] = {}
    if bead not in residue[group]:
        residue[group][bead]={} 
    return residue, group, bead   

def get_atomistic(frag_location, resname=False):
    if not resname:
        resname = frag_location.split('/')[-1][:-4]
#### read in atomistic fragments into dictionary    
    residue = {} ## a dictionary of bead in each residue eg residue[group][bead][atom number(1)][residue_name(ASP)/coordinates(coord)/atom name(C)/connectivity(2)/atom_mass(12)]
    fragment_mass = {}
    with open(frag_location, 'r') as pdb_input:
        for line_nr, line in enumerate(pdb_input.readlines()):
            if line.startswith('['):
                residue, group, bead = split_fragment_names(line, residue, resname)
                fragment_mass[bead]=[]
            if line.startswith('ATOM'):
                line_sep = gen.pdbatom(line) ## splits up pdb line
                residue[group][bead][line_sep['atom_number']]={'coord':np.array([line_sep['x']*g_var.sf,line_sep['y']*g_var.sf,line_sep['z']*g_var.sf]),
                                                                'atom':line_sep['atom_name'],'resid':1, 'resid_ori':line_sep['residue_id'],'res_type':line_sep['residue_name'],
                                                                'frag_mass':1}    
#### updates fragment mass   
                if not gen.is_hydrogen(line_sep['atom_name']):
                    if line_sep['atom_name'] in g_var.res_top[resname]['atom_masses']:
                        residue[group][bead][line_sep['atom_number']]['frag_mass']=g_var.res_top[resname]['atom_masses'][line_sep['atom_name']]  ### updates atom masses with crude approximations
                        fragment_mass[bead].append([line_sep['x']*g_var.sf,line_sep['y']*g_var.sf,line_sep['z']*g_var.sf,g_var.res_top[resname]['atom_masses'][line_sep['atom_name']]])               
                else:
                    fragment_mass[bead].append([line_sep['x']*g_var.sf,line_sep['y']*g_var.sf,line_sep['z']*g_var.sf,1])
    return residue, fragment_mass

def connectivity(cg, at_frag_centers, cg_frag_centers, group, group_number):
#### returns the connections between the atomistic and coarse grain beads
    at_connection, cg_connection=[],[]
    for bead in cg:
        resname = cg[bead]['residue_name']
        break
    if resname in g_var.sorted_connect:
        if len(g_var.sorted_connect[resname]) > 0:
            for group_bead in group: 
                for bead_atom in group[group_bead]:
                    if bead_atom in g_var.sorted_connect[resname][int(group_number)]:
                        skip = False
                        for bead_connect in g_var.sorted_connect[resname][int(group_number)][bead_atom]:
                            if 'cg_connect' not in locals():
                                if bead_connect in cg:
                                    cg_connect = [cg[bead_connect]['coord']]
                                else:
                                    skip = True
                            else:
                                if bead_connect in cg:
                                    cg_connect.append(cg[bead_connect]['coord'])
                        if not skip:
                            at_connection.append(group[group_bead][bead_atom]['coord'])
                            cg_connection.append(np.mean(np.array(cg_connect), axis=0))    
    if len(at_frag_centers)  > 1:         
        for bead in at_frag_centers:
            cg_connection.append(cg_frag_centers[bead])
            at_connection.append(at_frag_centers[bead])   
    return at_connection, cg_connection    

def get_rotation(cg_connect, at_connect, center, resname, group, cg_resid):
    if len(at_connect) == len(cg_connect):
        if len(cg_connect) == 0:
            return gen.AnglesToRotMat([np.random.uniform(0, math.pi*2), np.random.uniform(0, math.pi*2), np.random.uniform(0, math.pi*2)])
        else:
            return kabsch_rotate(np.array(at_connect)-center, np.array(cg_connect)-center)
    else:
        print('atom connections: '+str(len(at_connect))+' does not match CG connections: '+str(len(cg_connect)))
        sys.exit('residue number: '+str(cg_resid)+', residue type: '+str(resname)+', group: '+group)

def apply_rotations(atomistic_fragments,cg_resid, group_fit, center, xyz_rot_apply):
    for bead in group_fit:
        for atom in group_fit[bead]:
            group_fit[bead][atom]['coord'] = rotate_atom(group_fit[bead][atom]['coord'], center, xyz_rot_apply)   
            atomistic_fragments[cg_resid][atom] = group_fit[bead][atom].copy()
    return atomistic_fragments

### end fragment information
### linked residue connectivity start

def BB_connectivity(at_connections,cg_connections, cg_residues, at_residues, residue_number, BB_bead, resname):
    con_atoms = {}
    for atom in at_residues:
        # resname = at_residues[atom]['res_type']
        if at_residues[atom]['atom'] in g_var.res_top[resname]['CONNECT'][BB_bead]['atom']:
            con_atoms[g_var.res_top[resname]['CONNECT'][BB_bead]['atom'].index(at_residues[atom]['atom'])]=atom
    new_chain=False
    for con in con_atoms:
        con_resid = residue_number+g_var.res_top[resname]['CONNECT'][BB_bead]['dir'][con]
        if con_resid in cg_residues:
            xyz_cur = cg_residues[residue_number][BB_bead]['coord']
            if g_var.res_top[resname]['CONNECT'][BB_bead]['Con_Bd'][con] in cg_residues[con_resid]:
                xyz_con = cg_residues[con_resid][g_var.res_top[resname]['CONNECT'][BB_bead]['Con_Bd'][con]]['coord']
                if gen.calculate_distance(xyz_con, xyz_cur) < 6:
                    cg_connections.append(xyz_con)
                    at_connections.append(at_residues[con_atoms[con]]['coord'])
                else:
                    if g_var.res_top[resname]['CONNECT'][BB_bead]['dir'][con] > 0:
                        new_chain=True
            else:
                    if g_var.res_top[resname]['CONNECT'][BB_bead]['dir'][con] > 0:
                        new_chain=True
        else:
            if g_var.res_top[resname]['CONNECT'][BB_bead]['dir'][con] > 0:
                new_chain=True
    return at_connections,cg_connections, new_chain

################################################################### Merged system

def merge_indivdual_chain_pdbs(file, end, res_type):
#### reads in each chain into merge list
    merge, merged_coords = [],[]
    count = 0
    restraint_count = -1
    for chain in range(0,g_var.system[res_type]):
        merge_temp = []
        if os.path.exists(file+'_'+str(chain)+end):
            with open(file+'_'+str(chain)+end, 'r') as pdb_input:
                merge_temp += read_in.filter_input(pdb_input.readlines(), False)
        else:
            sys.exit('cannot find chain: '+file+'_'+str(chain)+end)
        if res_type+'_aligned' in file:  
            count, restraint_count = at_mod_p.create_disres(merge_temp, chain, file, count, restraint_count)
        merge, merge_coords = fix_chirality(merge,merge_temp,merged_coords, res_type)   
    if res_type+'_aligned' not in file:
        coords, index_conversion = index_conversion_generate(merge, merge_coords)
    if 'aligned' in file.split('/')[-1]:
        write_pdb(merge, merge_coords, {}, g_var.working_dir+'PROTEIN/PROTEIN_aligned_merged.pdb')
    else:
        write_pdb(merge, coords,index_conversion, g_var.working_dir+res_type+'/'+res_type+'_de_novo_merged.pdb')

def index_conversion_generate(merge, merge_coords):
    index_conversion = {}
    coords = []
    for at_val, atom in enumerate(merge):
        if not atom['atom_name'].startswith('M'):
            index_conversion[at_val] = len(coords)
            coords.append(merge_coords[at_val])
    coords = check_atom_overlap(coords)
    return coords, index_conversion

def write_pdb(merge, merge_coords, index_conversion, write_file):
#### creates merged pdb and writes chains to it
    if not os.path.exists(write_file):
        pdb_output=gen.create_pdb(write_file)
        atom_counter=0
        for line_val, line in enumerate(merge):
            if line_val in index_conversion:
                x, y, z = gen.trunc_coord(merge_coords[index_conversion[line_val]])
            else:
                if 'coord' in line:
                    x, y, z = gen.trunc_coord(line['coord'])
                else:
                    x, y, z = gen.trunc_coord([line['x'],line['y'],line['z']])
            if atom_counter >= 99_999:
                atom_counter=1
            else:
                atom_counter+=1      
            pdb_output.write(g_var.pdbline%((atom_counter, line['atom_name'], line['residue_name'],' ',line['residue_id'],\
            x,y,z,1,0))+'\n')
        pdb_output.close()

def merge_system_pdbs(protein):
    os.chdir(g_var.merged_directory)
#### create merged pdb 
    if not os.path.exists(g_var.merged_directory+'merged_cg2at'+protein+'.pdb'):
        merge=[]
        merge_coords=[]
    #### run through every residue type in cg_residues
        done = []
        if 'novo' in protein:
            print('checking for atom overlap in : '+protein[1:])
        for segment, residue_type in enumerate(g_var.system):
            if residue_type not in done:
                if residue_type not in ['PROTEIN', 'OTHER']:
                        start_num = len(merge_coords)
                        merge, merge_coords = read_in_merged_pdbs(merge, merge_coords, 
                                                                g_var.working_dir+residue_type+'/'+residue_type+'_merged.pdb')
                        g_var.np_blocks[residue_type]=[start_num,len(merge_coords)]
                elif residue_type == 'OTHER':
                    merge, merge_coords = read_in_merged_pdbs(merge, merge_coords, 
                                                            g_var.working_dir+residue_type+'/'+residue_type+'_de_novo_merged.pdb')
                elif residue_type == 'PROTEIN':
                    merge, merge_coords = read_in_merged_pdbs(merge, merge_coords, 
                                                            g_var.working_dir+residue_type+'/'+residue_type+protein+'_merged.pdb')     
                done.append(residue_type)
        index_conversion = {}
        if 'novo' in protein:
            coords, index_conversion = index_conversion_generate(merge, merge_coords)
            write_pdb(merge, coords, index_conversion, g_var.merged_directory+'merged_cg2at'+protein+'.pdb')
        else:
            write_pdb(merge, merge_coords, {}, g_var.merged_directory+'merged_cg2at'+protein+'.pdb')

def read_in_merged_pdbs(merge, merge_coords, location):
    if os.path.exists(location):
    #### opens pdb files and writes straight to merged_cg2at pdb
        with open(location, 'r') as pdb_input:
            read_in_atoms = read_in.filter_input(pdb_input.readlines(), False)
            merge_coords += [[line_sep['x'],line_sep['y'],line_sep['z']] for line_sep in read_in_atoms ]
            return merge+read_in_atoms, merge_coords
    else:
        sys.exit('cannot find minimised residue: \n'+ location) 

def check_overlap_chain(chain, input,res_type):
    if not os.path.exists(g_var.working_dir+res_type+'/'+res_type+'_'+input+str(chain)+'_gmx_checked.pdb'):
        lines, coords = read_in_merged_pdbs([], [], res_type+'_'+input+str(chain)+'_gmx.pdb')
        if res_type == 'PROTEIN':
            coords = at_mod_p.correct_amide_h(lines, coords)

        coords, index_conversion = index_conversion_generate(lines, coords)
        write_pdb(lines, coords, index_conversion, g_var.working_dir+res_type+'/'+res_type+'_'+input+str(chain)+'_gmx_checked.pdb')

## fix chirality errors

def fetch_chiral_coord(merge_temp, residue_type):
    chiral_atoms={}
    coord=[]
    for atom in range(len(merge_temp)):
        if residue_type in ['PROTEIN', 'OTHER']:
            resname = gen.check_alternate_resname(merge_temp[atom]['residue_name'])
        else:
            resname = gen.check_alternate_resname(residue_type)
        if len(g_var.res_top[resname]['CHIRAL']['atoms']) > 0:
            if merge_temp[atom]['atom_name'] in g_var.res_top[resname]['CHIRAL']['atoms']:
                if merge_temp[atom]['residue_id'] not in chiral_atoms:
                    chiral_atoms[merge_temp[atom]['residue_id']]={}
                chiral_atoms[merge_temp[atom]['residue_id']][merge_temp[atom]['atom_name']]=atom
        coord.append(np.array([merge_temp[atom]['x'],merge_temp[atom]['y'],merge_temp[atom]['z']]))
    return chiral_atoms, coord


def get_atom_move(merge_temp, resname, residue, chiral_group, chiral_atoms):
    stat = merge_temp[chiral_atoms[residue][chiral_group]].copy()
    atom_move = {'stat':np.array([stat['x'],stat['y'],stat['z']]), 'm':'', 'c1':'', 'c2':'', 'c3':''}
    for chir_atom in atom_move:
        if chir_atom != 'stat':
            test = merge_temp[chiral_atoms[residue][g_var.res_top[resname]['CHIRAL'][chiral_group][chir_atom]]].copy()
            atom_move[chir_atom]= np.array([test['x'],test['y'],test['z']])
            if gen.calculate_distance(atom_move['stat'], atom_move[chir_atom]) > 10:
                atom_move[chir_atom] = np.array(read_in.brute_mic(atom_move['stat'],atom_move[chir_atom]))
    return atom_move


def fix_chirality(merge, merge_temp, merged_coords, residue_type):
#### fixes chiral groups
    chiral_atoms, coord = fetch_chiral_coord(merge_temp, residue_type)
    for residue in chiral_atoms:
        if residue_type in ['PROTEIN', 'OTHER']:
            for atom in chiral_atoms[residue]:
                resname = merge_temp[chiral_atoms[residue][atom]]['residue_name']
                resname = gen.check_alternate_resname(resname)
                break
        else:
            resname = gen.check_alternate_resname(residue_type)

        for chiral_group in g_var.res_top[resname]['CHIRAL']:
            if chiral_group != 'atoms':
                atom_move = get_atom_move(merge_temp, resname, residue, chiral_group, chiral_atoms)

                S_M = atom_move['m'] - atom_move['stat']
                rotation = align_to_vector(S_M, [0,0,1])
                c1_coord = (atom_move['c1'] - atom_move['stat']).dot(rotation)
                c2_coord = (atom_move['c2'] - atom_move['stat']).dot(rotation)
                c3_coord = (atom_move['c3'] - atom_move['stat']).dot(rotation)

                if gen.angle_clockwise(c1_coord[0:2], c2_coord[0:2]) > gen.angle_clockwise(c1_coord[0:2], c3_coord[0:2]):
                    for ax_val, axis in enumerate(['x', 'y', 'z']):
                        merge_temp[chiral_atoms[residue][g_var.res_top[resname]['CHIRAL'][chiral_group]['m']]][axis] = merge_temp[chiral_atoms[residue][g_var.res_top[resname]['CHIRAL'][chiral_group]['m']]][axis] - (3*S_M[ax_val])   
                        merge_temp[chiral_atoms[residue][chiral_group]][axis] = merge_temp[chiral_atoms[residue][chiral_group]][axis] - (S_M[ax_val])        
                    coord[chiral_atoms[residue][g_var.res_top[resname]['CHIRAL'][chiral_group]['m']]] -=  (2*S_M) #move_coord -
                    coord[chiral_atoms[residue][chiral_group]] -=  (0.25*S_M) 
    merge+=merge_temp
    merged_coords+=coord
    return merge, merged_coords

## chirality end
## hydrogen orientaion checker

def check_hydrogens(residue):
#### finds the connecting carbons and their associated carbons [carbon atom, hydrogen ref number, connecting ref number]    
    for atom_num, atom in enumerate(residue):
        resname=residue[atom]['res_type']
        break
    for atom in g_var.hydrogen[resname]:
        h_coord = []
        for group in g_var.sorted_connect[resname]:
            if atom in g_var.sorted_connect[resname][group]:
                for hydrogen in g_var.hydrogen[resname][atom]:
                    h_coord.append(residue[hydrogen]['coord'])
                h_com=np.mean(np.array(h_coord), axis=0)
                for heavy_bond in g_var.heavy_bond[resname][atom]:
                    for group_check in g_var.sorted_connect[resname]:
                        if heavy_bond in g_var.sorted_connect[resname][group_check] and group_check != group:
                            skip = False
                            if heavy_bond in residue:
                                con_heavy_atom = heavy_bond
                                con_heavy_atom_co =  residue[con_heavy_atom]['coord']
                            else:
                                skip = True
                if not skip:
                #### vector between H COM and bonded carbon 
                    vector=np.array([h_com[0]-residue[atom]['coord'][0],h_com[1]-residue[atom]['coord'][1],h_com[2]-residue[atom]['coord'][2]])
                    h_com_f=h_com+vector*2
                    d1 = gen.calculate_distance(h_com, con_heavy_atom_co)    
                    d2 = gen.calculate_distance(h_com_f, con_heavy_atom_co)   
                    if d2 < d1:
                        for h_at in g_var.hydrogen[resname][atom]:
                            residue[h_at]['coord']=residue[h_at]['coord']-vector*2
    return residue

## check for threaded lipids
def check_ringed_lipids(protein):
    print('Checking for ringed lipids')
    if not os.path.exists(g_var.merged_directory+'checked_ringed_lipid_de_novo.pdb'):
        if not os.path.exists(g_var.merged_directory+'merged_cg2at_threaded.pdb'):
            os.chdir(g_var.merged_directory)
            merge, merge_coords = read_in_merged_pdbs([], [], protein)
            ringed=False
            lipid_atoms = []
            with open(g_var.merged_directory+'threaded_lipids.dat', 'w') as ring_ouput:
                for at_val, atom in enumerate(merge): 
                    resname = get_np_resname(at_val)
                    if resname in g_var.np_residues:
                        offset = fetch_start_of_residue_np(at_val, resname)
                        if atom['atom_number']-offset in g_var.heavy_bond[resname]:
                            for at_bond in g_var.heavy_bond[resname][atom['atom_number']-offset]:
                                at_bond -=1
                                if merge[at_bond+offset]['atom_number'] > merge[at_val]['atom_number']:
                                    merge[at_bond+offset]['x'], merge[at_bond+offset]['y'], merge[at_bond+offset]['z'] = np.array(read_in.brute_mic(merge_coords[at_val],merge_coords[at_bond+offset]))
                                    merge_coords[at_bond+offset] = merge[at_bond+offset]['x'], merge[at_bond+offset]['y'], merge[at_bond+offset]['z']
                                    dist = gen.calculate_distance(merge_coords[at_val], merge_coords[at_bond+offset])
                                    if 2 < dist < 6:
                                        lipid_atoms.append([at_val, at_bond+offset, (np.array(merge_coords[at_val])+np.array(merge_coords[at_bond+offset]))/2])
                                        ring_ouput.write('{0:6}{1:6}{2:2}{3:4}{4:2}{5:7}{6:5}{7:5}{8:5}{9:5}{10:5}{11:5}\n'.format(
                                                            'distance: ',str(np.round(dist,2)),'residue: ', merge[at_val]['residue_name'], merge[at_val]['residue_id'],
                                                            ' atom_1: ', merge[at_val]['atom_name'],
                                                            'atom_2: ', merge[at_bond+offset]['atom_name'], 'rough line num: ', at_val, at_bond+offset))
                                        ringed = True
        if ringed or os.path.exists(g_var.merged_directory+'merged_cg2at_threaded.pdb'):
            print('Found '+str(len(lipid_atoms))+' abnormal bonds, now attempting to fix.')
            print('See this file for a complete list: '+g_var.merged_directory+'threaded_lipids.dat')
            fix_threaded_lipids(lipid_atoms, merge, merge_coords)
        else:
            gen.file_copy_and_check(g_var.merged_directory+'MIN/merged_cg2at_de_novo_minimised.pdb', g_var.merged_directory+'checked_ringed_lipid_de_novo.pdb')

def fetch_start_of_residue(at, merge):
    count = at
    while True:
        if merge[count]['residue_id'] != merge[at]['residue_id']:
            count+=1
            return count
        else:
            count-=1

def fetch_start_of_residue_np(atom, resname):
            prev=g_var.np_blocks[resname][0]
            range_np_block = g_var.np_blocks[resname][1]-g_var.np_blocks[resname][0]
            for range_np_block_ind in range(g_var.np_blocks[resname][0],g_var.np_blocks[resname][1]+1, int(range_np_block/g_var.system[resname])):
                if prev <= atom < range_np_block_ind:
                    return prev
                prev = range_np_block_ind

def get_np_resname(atom):
    for res_check in g_var.np_blocks:
        if g_var.np_blocks[res_check][0] <= atom < g_var.np_blocks[res_check][1]:
            return res_check 

def fix_threaded_lipids(lipid_atoms, merge, merge_coords):
    if not os.path.exists(g_var.merged_directory+'merged_cg2at_threaded.pdb'):
        tree = cKDTree(merge_coords)         
        for threaded in lipid_atoms:  
            resname = get_np_resname(threaded[0]) 
            
            atoms = tree.query_ball_point(threaded[2], r=3)
            for at in atoms:
                if merge[at]['residue_id'] != merge[threaded[0]]['residue_id']:
                    P_count = fetch_start_of_residue(at, merge)
                    break

            NP_count = fetch_start_of_residue_np(threaded[0], resname)
            bb = []
            if 'P_count' not in locals():
                sys.exit('There is an issue with the bond length detection')
            for at in merge[P_count:]:
                if at['residue_id'] != merge[P_count]['residue_id']:
                    break
                if at['atom_name'] in g_var.res_top[at['residue_name']]['ATOMS']:
                    bb.append([at['x'], at['y'], at['z']])
            bb = np.mean(np.array(bb), axis=0)
            BB_M3 = (threaded[2]-bb)/np.linalg.norm((threaded[2]-bb))      
            
            for heavy_atom in threaded[:2]:
                merge_coords[heavy_atom] += BB_M3*3 
                merge[heavy_atom]['x'], merge[heavy_atom]['y'], merge[heavy_atom]['z'] = merge_coords[heavy_atom]
                for hydrogen in g_var.hydrogen[resname][heavy_atom-NP_count+1]:
                    merge_coords[NP_count+hydrogen-1] += BB_M3*3 
                    merge[NP_count+hydrogen-1]['x'], merge[NP_count+hydrogen-1]['y'], merge[NP_count+hydrogen-1]['z'] = merge_coords[NP_count+hydrogen-1]

        coords, index_conversion = index_conversion_generate(merge, merge_coords)
        write_pdb(merge, coords, index_conversion, g_var.merged_directory+'merged_cg2at_threaded.pdb')

    if not os.path.exists(g_var.merged_directory+'MIN/merged_cg2at_threaded_minimised.pdb'):
        gro.minimise_merged_pdbs('_threaded')
    gen.file_copy_and_check(g_var.merged_directory+'MIN/merged_cg2at_threaded_minimised.pdb', g_var.merged_directory+'checked_ringed_lipid_de_novo.pdb')
