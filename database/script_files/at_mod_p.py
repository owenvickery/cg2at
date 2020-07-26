#!/usr/bin/env python3

import os, sys
import numpy as np
import copy
from string import ascii_uppercase
import difflib
from scipy.spatial import cKDTree
import multiprocessing as mp
import gen, g_var, f_loc, at_mod, read_in
import math

def build_protein_atomistic_system(cg_residues):

    
#### initisation of counters
    chain_count=0
    g_var.backbone_coords[chain_count]=[]
    g_var.p_system[chain_count]=[]
    g_var.coord_atomistic[chain_count]={}
    g_var.seq_cg[chain_count]=[]
    print('Converting Protein')
    gen.mkdir_directory(g_var.working_dir+'PROTEIN')  ### make and change to protein directory
#### for each residue in protein
    initial=True
    residue_type={}
    residue_type_mass={}
    for cg_residue_id, residue_number in enumerate(cg_residues):
        resname = cg_residues[residue_number][next(iter(cg_residues[residue_number]))]['residue_name']
        BB_bead = f_loc.res_top[resname]['BACKBONE']
        if cg_residue_id == 0:
            g_var.p_system[chain_count].append(terminal_residue(resname)[0]) 
        g_var.coord_atomistic[chain_count][residue_number]={}
        frag_location=at_mod.fragment_location(resname) ### get fragment location from database
        residue_type[resname], residue_type_mass[resname] = at_mod.get_atomistic(frag_location)
        for group in residue_type[resname]:
            center, at_frag_centers, cg_frag_centers, group_fit = at_mod.rigid_fit(residue_type[resname][group], residue_type_mass[resname], residue_number, cg_residues[residue_number])
            at_connect, cg_connect = at_mod.connectivity(cg_residues[residue_number], at_frag_centers, cg_frag_centers, group_fit, group)
            if BB_bead in group_fit:
                BB_connect = []
                at_connect, cg_connect, new_chain = BB_connectivity(at_connect,cg_connect, cg_residues, group_fit[BB_bead], residue_number, BB_bead)
                g_var.seq_cg = at_mod.add_to_sequence(g_var.seq_cg, resname, chain_count)
                g_var.backbone_coords[chain_count].append(np.append(cg_residues[residue_number][BB_bead]['coord'], 1))     
            if len(at_connect) == len(cg_connect) and len(at_connect) != 0:
                xyz_rot_apply=at_mod.kabsch_rotate(np.array(at_connect)-center, np.array(cg_connect)-center)
            elif len(at_connect) == 0:
                xyz_rot_apply = False
                print('Cannot find any connectivity for residue number: '+str(residue_number)+', residue type: '+str(resname)+', group: '+str(group))
            else:
                print('atom connections: '+str(len(at_connect))+' does not equal CG connections: '+str(len(cg_connect)))
                sys.exit('residue number: '+str(residue_number)+', residue type: '+str(resname)+', group: '+group)

            for bead in group_fit:
                for atom in group_fit[bead]:
                    group_fit[bead][atom]['coord'] = at_mod.rotate_atom(group_fit[bead][atom]['coord'], center, xyz_rot_apply)   
                    atom_new = group_fit[bead][atom].copy()
                    g_var.coord_atomistic[chain_count][residue_number][atom] = atom_new
        if new_chain:
            g_var.p_system[chain_count].append(terminal_residue(resname)[1])
            chain_count+=1
            if cg_residue_id+1 != len(cg_residues):
                g_var.backbone_coords[chain_count]=[]
                g_var.coord_atomistic[chain_count]={}
                g_var.p_system[chain_count]=[]
                BB_bead = f_loc.res_top[resname]['BACKBONE']
                g_var.p_system[chain_count].append(terminal_residue(cg_residues[residue_number+1][next(iter(cg_residues[residue_number+1]))]['residue_name'])[0])
                g_var.seq_cg[chain_count]=[]
    if g_var.v >=1:
        print('\n{0:^15}{1:^12}'.format('chain number', 'length of chain')) #   \nchain number\tDelta A\t\tno in pdb\tlength of chain')
        print('\n{0:^15}{1:^12}'.format('------------', '---------------'))
        for chain in g_var.seq_cg:
            print('{0:^15}{1:^12}'.format(chain, len(g_var.seq_cg[chain])))
        print()
    g_var.system['PROTEIN']=chain_count

def terminal_residue(resname):
    ter = [False, False]
    if f_loc.res_top[resname]['N_TERMINAL'] == 'TER':
        ter[0] = True
    if f_loc.res_top[resname]['C_TERMINAL'] == 'TER':
        ter[1] = True
    return ter


def BB_connectivity(at_connections,cg_connections, cg_residues, at_residues, residue_number, BB_bead):
    for atom in at_residues:
        resname = at_residues[atom]['res_type']
        # print(at_residues)
        if at_residues[atom]['atom'] == f_loc.res_top[at_residues[atom]['res_type']]['N_TERMINAL']:
            N_ter=atom
        if at_residues[atom]['atom'] == f_loc.res_top[at_residues[atom]['res_type']]['C_TERMINAL']:
            C_ter=atom
#### connect to preceeding backbone bead in chain
    BB_cur = f_loc.res_top[cg_residues[residue_number][next(iter(cg_residues[residue_number]))]['residue_name']]['BACKBONE']
    new_chain=False
    if residue_number-1 in cg_residues and 'N_ter' in locals(): 
        prev_resname = cg_residues[residue_number-1][next(iter(cg_residues[residue_number-1]))]['residue_name']
        BB_prev = f_loc.res_top[prev_resname]['BACKBONE']
        xyz_cur = cg_residues[residue_number][BB_cur]['coord']
        xyz_prev = cg_residues[residue_number-1][BB_prev]['coord']
        if gen.calculate_distance(xyz_prev, xyz_cur) < 6 and f_loc.res_top[resname]['C_TERMINAL'] != 'TER':
            cg_connections.append(cg_residues[residue_number-1][BB_prev]['coord'])
            at_connections.append(at_residues[N_ter]['coord'])
#### connect to next backbone bead in chain
    if residue_number+1 in cg_residues and 'C_ter' in locals(): 
        next_resname = cg_residues[residue_number+1][next(iter(cg_residues[residue_number+1]))]['residue_name']
        BB_next = f_loc.res_top[next_resname]['BACKBONE']
        xyz_cur = cg_residues[residue_number][BB_cur]['coord']
        xyz_next = cg_residues[residue_number+1][BB_next]['coord']
        if gen.calculate_distance(xyz_next, xyz_cur) < 6 and f_loc.res_top[next_resname]['N_TERMINAL'] != 'TER':
            cg_connections.append(cg_residues[residue_number+1][BB_next]['coord'])
            at_connections.append(at_residues[C_ter]['coord'])
        else:
            new_chain=True
    else:
        new_chain=True
    return at_connections,cg_connections, new_chain

################# Fixes disulphide bond, martini cysteine bone is too far apart to be picked up by pdb2gmx. 
#### 

def ask_if_disulphide(chain, res_1, res_2):
    while True:
        try:
            answer = str(input('\nAre these residues connected by a disulphide bond in chain '+str(chain)+' (Y/N): '+str(res_1)+'--'+str(res_2)+': '))
            # print(answer)
            if answer.lower() in ['yes','y']:
                return True
            elif answer.lower() in ['no','n']:
                return False
            else:
                print("Oops!  That was a invalid choice")
        except KeyboardInterrupt:
            sys.exit('\nInterrupted')
        except:
            print("Oops!  That was a invalid choice")

def find_disulphide_bonds_user_sup(atomistic_protein_input):
    for chain in atomistic_protein_input:
        g_var.user_cys_bond[chain]=[]
        cysteines = []
        cys_resid = []
        for part in atomistic_protein_input[chain]:
            for resid in atomistic_protein_input[chain][part]:
                atom = next(iter(atomistic_protein_input[chain][part][resid]))
                if atomistic_protein_input[chain][part][resid][atom]['res_type'] == 'CYS':
                    for atom in atomistic_protein_input[chain][part][resid]:
                        if 'S' in atomistic_protein_input[chain][part][resid][atom]['atom'].upper() :
                            cys_resid.append(resid)
                            cysteines.append(atomistic_protein_input[chain][part][resid][atom]['coord'])
        if len(cysteines)>=2:
            tree = cKDTree(cysteines)
            done_query=[]
            for cys_index, cys in enumerate(cysteines):
                query = tree.query_ball_point(cys, r=2.1)
                if len(query) == 2 and query not in done_query:
                    g_var.user_cys_bond[chain].append([cys_resid[query[0]],cys_resid[query[1]]])
                    done_query.append(query)

def find_disulphide_bonds_de_novo():
    for chain in g_var.coord_atomistic:
        offset = next(iter(g_var.coord_atomistic[chain]))
        if chain not in g_var.user_cys_bond:
            g_var.user_cys_bond[chain]=[]
        cysteines = []
        cys_resid = []
        for resid in g_var.coord_atomistic[chain]:
            resid_corrected=np.array(resid)-offset
            run = True
            if len(g_var.user_cys_bond) > 0:
                for bonds in g_var.user_cys_bond[chain]:
                    if resid_corrected in bonds:
                        run=False
            if run:
                atom = next(iter(g_var.coord_atomistic[chain][resid]))
                if g_var.coord_atomistic[chain][resid][atom]['res_type'] == 'CYS':
                    for atom in g_var.coord_atomistic[chain][resid]:
                        if 'S' in g_var.coord_atomistic[chain][resid][atom]['atom'].upper():
                            cys_resid.append(resid_corrected)
                            cysteines.append(g_var.coord_atomistic[chain][resid][atom]['coord'])
        if len(cys_resid)>=2:
            tree = cKDTree(cysteines)
            done_query=[]
            for cys_index, cys in enumerate(cysteines):
                query = tree.query_ball_point(cys, r=g_var.cys)
                if len(query) == 2 and query not in done_query:
                    if cys_resid[query[1]] - cys_resid[query[0]] > 4:
                        if g_var.silent==True:
                            disul = True
                        else:
                            disul = ask_if_disulphide(chain, cys_resid[query[0]],cys_resid[query[1]])
                        if disul:
                            g_var.user_cys_bond[chain].append([cys_resid[query[0]],cys_resid[query[1]]])
                    done_query.append(query)

def correct_disulphide_bonds(coordinates_atomistic):
    for chain in coordinates_atomistic:
        offset = next(iter(coordinates_atomistic[chain]))
        for bonds in g_var.user_cys_bond[chain]:
            bonds_corrected=np.array(bonds)+offset
            at_num=[]
            at_coord=[]
            for resid in bonds_corrected: 
                for atom in coordinates_atomistic[chain][resid]:
                    if 'S' in coordinates_atomistic[chain][resid][atom]['atom'].upper():
                        at_num.append(atom)
                        at_coord.append(coordinates_atomistic[chain][resid][atom]['coord'])
            new_c1, new_c2 = shrink_coordinates(at_coord[0],at_coord[1])
            coordinates_atomistic[chain][bonds_corrected[0]][at_num[0]]['coord']=new_c1
            coordinates_atomistic[chain][bonds_corrected[1]][at_num[1]]['coord']=new_c2
    return coordinates_atomistic

def shrink_coordinates(c1,c2):
    dist = gen.calculate_distance(c1,c2)
    vector = c1-c2
    scale=0.05
    while True:
        new_c1 = c1 - (vector*scale)
        new_c2 = c2 + (vector*scale)
        dist = gen.calculate_distance(new_c1,new_c2)
        if dist > 2.1:
            scale+=0.0025
        elif dist < 1.9:
            scale-=0.001
        else:
            return new_c1, new_c2
    return c1, c2

def finalise_novo_atomistic(cg_residues):
    final_at_residues={}
    final_at = {}
    for chain in g_var.coord_atomistic: 
        at_number=0    
        final_at_residues[chain]={}  
        final_at[chain]={}
        coords=[]
        skip = os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN_de_novo_'+str(chain)+'.pdb')
        if not skip:
            pdb_output = gen.create_pdb(g_var.working_dir+'PROTEIN/PROTEIN_de_novo_'+str(chain)+'.pdb')
        for res_index, residue_id in enumerate(g_var.coord_atomistic[chain]):
            if g_var.coord_atomistic[chain][residue_id][1]['res_type'] in f_loc.mod_residues:
                g_var.coord_atomistic[chain][residue_id] = at_mod.check_hydrogens(g_var.coord_atomistic[chain][residue_id])
            if res_index <= len(g_var.coord_atomistic[chain])-3:
                g_var.coord_atomistic[chain][residue_id], cross_vector = fix_carbonyl(residue_id, cg_residues, g_var.coord_atomistic[chain][residue_id], False)
            elif res_index < len(g_var.coord_atomistic[chain]) and 'cross_vector' in locals():
                g_var.coord_atomistic[chain][residue_id], cross_vector = fix_carbonyl(residue_id, cg_residues, g_var.coord_atomistic[chain][residue_id], cross_vector)
            for at_val, atom in enumerate(g_var.coord_atomistic[chain][residue_id]):
                g_var.coord_atomistic[chain][residue_id][at_val+1]['resid'] = res_index
                coords.append(g_var.coord_atomistic[chain][residue_id][at_val+1]['coord'])
                final_at[chain][at_number]=g_var.coord_atomistic[chain][residue_id][at_val+1]
                at_number+=1
            final_at_residues[chain][res_index]=g_var.coord_atomistic[chain][residue_id]

        coords = at_mod.check_atom_overlap(coords)
        for atom in final_at[chain]:
            final_at[chain][atom]['coord']=coords[atom]
            final_at_residues[chain][final_at[chain][atom]['resid']][atom]=final_at[chain][atom]
            if not skip:
                x, y, z = gen.trunc_coord(final_at[chain][atom]['coord'])
                pdb_output.write(g_var.pdbline%((atom,final_at[chain][atom]['atom'],final_at[chain][atom]['res_type'],' ',\
                                final_at[chain][atom]['resid'],x,y,z,1,0))+'\n')
    if 'pdb_output' in locals():
        pdb_output.close()
    return final_at_residues

########################################### fix carbonyl section 

def fix_carbonyl(residue_id, cg, at, cross_vector):
    chiral = {}
    carbonyl = {}
    for atom in at:
        if at[atom]['atom'] in ['N','CA', 'C', 'O']:
            carbonyl[at[atom]['atom']] = atom
        if at[atom]['atom'] in f_loc.res_top[at[atom]['res_type']]['CHIRAL']['atoms']:
            chiral[at[atom]['atom']] = atom
    ca=[]
    prev_BB, next_BB=False, False
    cur_BB =  cg[residue_id][f_loc.res_top[cg[residue_id][next(iter(cg[residue_id]))]['residue_name']]['BACKBONE']]['coord']
    if residue_id-1 in cg:
        BB_bead = f_loc.res_top[cg[residue_id-1][next(iter(cg[residue_id-1]))]['residue_name']]['BACKBONE']
        prev_BB = cg[residue_id-1][BB_bead]['coord']
    if residue_id+1 in cg:
        BB_bead = f_loc.res_top[cg[residue_id+1][next(iter(cg[residue_id+1]))]['residue_name']]['BACKBONE']
        next_BB = cg[residue_id+1][BB_bead]['coord']
    if not np.any(cross_vector):
        for index in range(3):
            BB_bead = f_loc.res_top[cg[residue_id+index][next(iter(cg[residue_id+index]))]['residue_name']]['BACKBONE']
            ca.append(cg[residue_id+index][BB_bead]['coord'])
        cross_vector = at_mod.find_cross_vector( ca )
    rotation = at_mod.align_to_vector(at_mod.noramlised_vector(at[carbonyl['O']]['coord'],at[carbonyl['C']]['coord']), cross_vector)
    at[carbonyl['O']]['coord'] = (at[carbonyl['O']]['coord']-at[carbonyl['C']]['coord']).dot(rotation)+at[carbonyl['C']]['coord']
    at[carbonyl['C']]['coord'] = at[carbonyl['C']]['coord'] + cross_vector*0.2
    at[carbonyl['N']]['coord'] = at[carbonyl['N']]['coord'] - cross_vector*0.5
    if len(f_loc.res_top[at[1]['res_type']]['CHIRAL']) >= 2:
        for chiral_group in f_loc.res_top[at[1]['res_type']]['CHIRAL']:
            if chiral_group != 'atoms':
                p1 = chiral[chiral_group]
                c1 = chiral[f_loc.res_top[at[1]['res_type']]['CHIRAL'][chiral_group]['c1']]
                c2 = chiral[f_loc.res_top[at[1]['res_type']]['CHIRAL'][chiral_group]['c2']]
                c3 = chiral[f_loc.res_top[at[1]['res_type']]['CHIRAL'][chiral_group]['c3']]
                cross_vector_chiral = at_mod.find_cross_vector( [at[c3]['coord'], at[c2]['coord'], at[c1]['coord']])
                at[p1]['coord'] = at[p1]['coord'] + cross_vector_chiral*0.5
                if f_loc.res_top[at[1]['res_type']]['CHIRAL'][chiral_group]['m'] in chiral:
                    m = chiral[f_loc.res_top[at[1]['res_type']]['CHIRAL'][chiral_group]['m']]
                    at[m]['coord'] = at[m]['coord'] + cross_vector_chiral*1
    return at, cross_vector

  

def check_sequence(atomistic_protein_input, chain_count):
    for chain in range(len(atomistic_protein_input)):
        g_var.seq_at[chain]=[]
        for resid in atomistic_protein_input[chain]:
            for atom in atomistic_protein_input[chain][resid]:
                g_var.seq_at = at_mod.add_to_sequence(g_var.seq_at, atomistic_protein_input[chain][resid][atom]['res_type'], chain)
                break

def align_chains(atomistic_protein_input):
    if g_var.v >= 2:
        print('coarse grain protein sequence:\n')
        for index in g_var.seq_cg:
            print('chain:', index,g_var.seq_cg[index], '\n')
        print('\nuser supplied structure:\n')
        for index in g_var.seq_at:
            print('chain:', index, g_var.seq_at[index], '\n')        
    at={}
    sequence_temp = g_var.seq_cg.copy()
    test_chain={}
    for chain_at in range(len(atomistic_protein_input)):
        skip_sequence=False
        chain_cg=0
        s = difflib.SequenceMatcher(None, g_var.seq_at[chain_at], g_var.seq_cg[chain_cg], autojunk=False)
        seq_info = s.get_matching_blocks()
        while seq_info[0][2] != len(g_var.seq_at[chain_at]):
            if chain_cg >= len(g_var.seq_cg)-1:
                print('\nCannot find a match for user supplied chain: '+str(chain_at))#+'\n\nAtomistic chain:\n'+str(seq_user[chain_at]),'\n\nIn CG:\n'+str(sequence))
                skip_sequence = True
                break
            if not skip_sequence:
                chain_cg+=1
                s = difflib.SequenceMatcher(None, g_var.seq_at[chain_at], g_var.seq_cg[chain_cg], autojunk=False)
                seq_info = s.get_matching_blocks()
        if not skip_sequence:
            temp={}
            if chain_cg not in at:
                at[chain_cg]={}
            if seq_info[0][2] == len(g_var.seq_at[chain_at]):
                for resid,  residue in enumerate(atomistic_protein_input[chain_at]):
                    temp[resid + seq_info[0][1]] = atomistic_protein_input[chain_at][residue]
                at[chain_cg][str(seq_info[0][1])+':'+str(seq_info[0][1]+seq_info[0][2])]=temp  
            g_var.seq_cg[chain_cg] = mask_sequence(g_var.seq_cg[chain_cg], seq_info[0][1], seq_info[0][1]+seq_info[0][2])
            test_chain[chain_at]=chain_cg

    if len(at) > 0:
        g_var.user_at_input = True
    else: 
        g_var.user_at_input = False
    if f_loc.group_chains == 'chain':
        return at, test_chain
    else:
        return at, f_loc.group_chains 

def mask_sequence(sequence, st, end):
    for index, residue in enumerate(sequence):
        if st <= index < end:
            sequence[index]='-'
    return sequence

def center_atomistic(atomistic_protein_input, group_chain):
    cg_com={}
#### for each protein chain center on cg representation 
    for chain in atomistic_protein_input:
        cg_com[chain]=[]
        for part_val, part in enumerate(atomistic_protein_input[chain]):
            sls, sle= int(part.split(':')[0]),int(part.split(':')[1])
            if group_chain==None:
                protein_mass=fetch_backbone_mass(atomistic_protein_input[chain][part], [])
                atomistic_protein_mass = at_mod.COM(protein_mass, 'protein at: '+str(chain)+' '+part)#### add center of mass of CG_proteins
                cg_com[chain].append(at_mod.COM(g_var.backbone_coords[chain][sls:sle], 'protein cg: '+str(chain)+' '+part))
                atomistic_protein_input[chain][part] = update_part_coordinate(atomistic_protein_input[chain][part], atomistic_protein_mass, cg_com[chain][part_val])
            elif group_chain=='all':
                if 'protein_mass' not in locals():
                    protein_mass=[]
                protein_mass=fetch_backbone_mass(atomistic_protein_input[chain][part], protein_mass)
                if 'cg_backbone_masses' not in locals():
                    cg_backbone_masses = {}
                    cg_backbone_masses['all'] = g_var.backbone_coords[chain][sls:sle]
                else:
                    cg_backbone_masses['all'] += g_var.backbone_coords[chain][sls:sle]
            else:
                if 'protein_mass' not in locals():
                    protein_mass={}
                if chain in group_chain:
                    if group_chain[chain] not in protein_mass:
                        protein_mass[group_chain[chain]]=fetch_backbone_mass(atomistic_protein_input[chain][part], [])
                    else:
                        protein_mass[group_chain[chain]]=fetch_backbone_mass(atomistic_protein_input[chain][part], protein_mass[group_chain[chain]])
                    if 'cg_backbone_masses' not in locals():
                        cg_backbone_masses = {}
                    if group_chain[chain] not in cg_backbone_masses:
                        cg_backbone_masses[group_chain[chain]]= g_var.backbone_coords[chain][sls:sle]
                    else:
                        cg_backbone_masses[group_chain[chain]] += g_var.backbone_coords[chain][sls:sle]                  
    if group_chain=='all':
        atomistic_protein_input, cg_com = center_at_protein_all_chains(atomistic_protein_input, cg_com, protein_mass, cg_backbone_masses, group_chain)
    if group_chain not in ['all',None]:
        atomistic_protein_input, cg_com = center_at_protein_chain_groups(atomistic_protein_input, cg_com, protein_mass, cg_backbone_masses, group_chain)
    return atomistic_protein_input, cg_com

def center_at_protein_chain_groups(atomistic_protein_input, cg_com, protein_mass, cg_backbone_masses,  group_chain):
    for chain in atomistic_protein_input:
        for part_val, part in enumerate(atomistic_protein_input[chain]):
            sls, sle= int(part.split(':')[0]),int(part.split(':')[1])
            if chain in group_chain:
                atomistic_protein_mass = at_mod.COM(protein_mass[group_chain[chain]], 'AT protein chain: '+str(chain))
                cg_com[chain].append(at_mod.COM(cg_backbone_masses[group_chain[chain]], 'CG protein chain: '+str(chain)))
            else:
                protein_mass=fetch_backbone_mass(atomistic_protein_input[chain][part], [])
                atomistic_protein_mass = at_mod.COM(protein_mass, 'protein at: '+str(chain)+' '+part)#### add center of mass of CG_proteins
                cg_com[chain].append(at_mod.COM(g_var.backbone_coords[chain][sls:sle], 'protein cg: '+str(chain)+' '+part))
            atomistic_protein_input[chain][part] = update_part_coordinate(atomistic_protein_input[chain][part], atomistic_protein_mass, cg_com[chain][part_val])
    return atomistic_protein_input, cg_com


def center_at_protein_all_chains(atomistic_protein_input, cg_com,  protein_mass, cg_backbone_masses, group_chain):
    atomistic_protein_mass = at_mod.COM(protein_mass, 'All AT protein chains')
    for chain in atomistic_protein_input:
        for part_val, part in enumerate(atomistic_protein_input[chain]):
            cg_com[chain].append(at_mod.COM(cg_backbone_masses['all'], 'All CG protein chains'))
            sls, sle= int(part.split(':')[0]),int(part.split(':')[1])
            atomistic_protein_input[chain][part] = update_part_coordinate(atomistic_protein_input[chain][part], atomistic_protein_mass, cg_com[chain][part_val])
    return atomistic_protein_input, cg_com

def update_part_coordinate(part, at_com, cg_com):
    #### each atoms coord is updated so the AT COM is the same as the CG
    for residue in part:
        for atom in part[residue]:  
            part[residue][atom]['coord']=part[residue][atom]['coord']-(at_com-cg_com)
    return part

def fetch_backbone_mass(part, protein_mass):
    for residue in part:
    #### creates a list of all coordinates and masses [[coord, mass],[coord, mass]]
        for atom in part[residue]:
            if part[residue][atom]['atom'] in f_loc.res_top[part[residue][atom]['res_type']]['ATOMS']:
                protein_mass.append([part[residue][atom]['coord'][0],part[residue][atom]['coord'][1],part[residue][atom]['coord'][2],part[residue][atom]['frag_mass']])    
    return protein_mass

def rotate_protein_monomers(atomistic_protein_centered, final_coordinates_atomistic, cg_com, group_chain):
#### run through each chain in proteins
    at_com_group, cg_com_group={},{}
    for chain in range(len(final_coordinates_atomistic)):
    #### creates atomistic pdb
        if chain in atomistic_protein_centered:
            for part_val, part in enumerate(atomistic_protein_centered[chain]):
                sls, sle= int(part.split(':')[0]),int(part.split(':')[1])        
                at_centers=[]
            #### runs through every residue and atom  
                for residue in atomistic_protein_centered[chain][part]:
                #### gets center of mass of each residue (note only backbone heavy atoms have a mass)
                    at_centers_iter=[]

                    for atom in atomistic_protein_centered[chain][part][residue]:

                        if atomistic_protein_centered[chain][part][residue][atom]['atom'] in f_loc.res_top[atomistic_protein_centered[chain][part][residue][atom]['res_type']]['ATOMS']:
                            at_centers_iter.append(np.append(atomistic_protein_centered[chain][part][residue][atom]['coord'],atomistic_protein_centered[chain][part][residue][atom]['frag_mass']))
                    try:
                        at_centers.append(np.average(np.array(at_centers_iter)[:,:3], axis=0, weights=np.array(at_centers_iter)[:,3]))
                    except:
                        for atom in atomistic_protein_centered[chain][part][residue]:
                            print(atomistic_protein_centered[chain][part][residue][atom])
                        sys.exit()
                if len(at_centers) == len(np.array(g_var.backbone_coords[chain])[sls:sle,:3]):
                    if group_chain=='all':
                        at_com_group, cg_com_group = return_all_rotations(chain, at_centers, at_com_group, cg_com_group, sls, sle)
                    elif group_chain==None:
                        at_com_group = return_indivdual_rotations(chain, part_val, at_centers, cg_com, at_com_group, cg_com_group, sls, sle)
                    else:
                        at_com_group, cg_com_group = return_grouped_rotations(chain, part_val, at_centers, cg_com, at_com_group, cg_com_group, sls, sle,group_chain)
                else:
                    sys.exit('In chain '+str(chain)+' the atomistic input does not match the CG. \n\
                            number of CG residues '+str(len(g_var.backbone_coords[chain]))+'\nnumber of AT residues '+str(len(at_centers)))
    return at_com_group, cg_com_group

def apply_rotations_to_chains(final_coordinates_atomistic, atomistic_protein_centered, at_com_group,cg_com_group,cg_com, group_chain):
    if group_chain=='all':
        rotate_all = return_all_rotations_final(at_com_group,cg_com_group,cg_com)
    final_rotations = []
    final_user_supplied_coord={}
    for chain in range(len(final_coordinates_atomistic)):
        rotations=[]
        if chain in atomistic_protein_centered:
            if group_chain not in ['all',None]: 
                if chain in group_chain:
                    rotate_chain=at_mod.kabsch_rotate(at_com_group[group_chain[chain]]-cg_com[group_chain[chain]][0], cg_com_group[group_chain[chain]]-cg_com[group_chain[chain]][0])
            for part_val, part in enumerate(atomistic_protein_centered[chain]):
                if 'rotate_all' in locals():
                    rotations.append(rotate_all)
                elif group_chain==None:
                    rotations = at_com_group[chain]
                else:
                    if chain in group_chain:
                        rotations.append(rotate_chain)
                    else:
                        rotations = at_com_group[chain]

            final_user_supplied_coord[chain] = hybridise_protein_inputs(final_coordinates_atomistic[chain], atomistic_protein_centered[chain], cg_com[chain], rotations, chain)
        else:
            final_user_supplied_coord[chain] = hybridise_protein_inputs(final_coordinates_atomistic[chain], [], [], [], chain)
    return final_user_supplied_coord

def return_all_rotations_final(at_com_group,cg_com_group,cg_com):
    if len(at_com_group['all']) == len(cg_com_group['all']):
        return at_mod.kabsch_rotate(at_com_group['all']-cg_com[0][0], cg_com_group['all']-cg_com[0][0])
    else:
        sys.exit('The atomistic input does not match the CG. \n\
                number of CG residues '+str(len(cg_com_group['all']))+'\nnumber of AT residues '+str(len(at_com_group['all'])))

def return_indivdual_rotations(chain, part_val, at_centers, cg_com, at_com_group, cg_com_group, sls, sle):
    if chain not in at_com_group:
        at_com_group[chain]=[]
## finds optimal rotation of each monomer  
    if len(at_centers) == len(np.array(g_var.backbone_coords[chain])[sls:sle,:3]):
        at_com_group[chain].append(at_mod.kabsch_rotate(np.array(at_centers)-cg_com[chain][part_val], 
                             np.array(g_var.backbone_coords[chain])[sls:sle,:3]-cg_com[chain][part_val]))
    else:
        sys.exit('In chain '+str(chain)+' the atomistic input does not match the CG. \n\
                    number of CG residues '+str(len(g_var.backbone_coords[chain]))+'\nnumber of AT residues '+str(len(at_centers)))    
    return at_com_group

def return_grouped_rotations(chain, part_val, at_centers, cg_com, at_com_group, cg_com_group, sls, sle, group_chain):
    if chain in group_chain:
        if group_chain[chain] not in at_com_group:
            at_com_group[group_chain[chain]]=[]   
            at_com_group[group_chain[chain]]=np.array(at_centers)
            cg_com_group[group_chain[chain]]=np.array(g_var.backbone_coords[chain])[sls:sle,:3]  
        else:
            at_com_group[group_chain[chain]]=np.append(at_com_group[group_chain[chain]], np.array(at_centers), axis=0)
            cg_com_group[group_chain[chain]]=np.append(cg_com_group[group_chain[chain]], np.array(g_var.backbone_coords[chain])[sls:sle,:3], axis=0)                                                                               
    else:
        at_com_group = return_indivdual_rotations(chain, part_val, at_centers, cg_com, at_com_group, cg_com_group, sls, sle)
    return at_com_group, cg_com_group

def return_all_rotations(chain, at_centers, at_com_group, cg_com_group, sls, sle):
    if 'all' not in at_com_group:
        at_com_group['all']=np.array(at_centers)
        cg_com_group['all']=np.array(g_var.backbone_coords[chain])[sls:sle,:3]
    else:
        at_com_group['all']=np.append(at_com_group['all'], np.array(at_centers), axis=0)
        cg_com_group['all']=np.append(cg_com_group['all'], np.array(g_var.backbone_coords[chain])[sls:sle,:3], axis=0) 
    return at_com_group, cg_com_group   

def hybridise_protein_inputs(final_coordinates_atomistic, atomistic_protein_centered, cg_com, xyz_rot_apply, chain):

    complete_user_at = {}
    for residue in final_coordinates_atomistic:
        exists=False
        resname = final_coordinates_atomistic[residue][next(iter(final_coordinates_atomistic[residue]))]['res_type']
        if resname in f_loc.mod_residues:
            complete_user_at[residue]=final_coordinates_atomistic[residue]
        elif resname not in f_loc.mod_residues:
            for part_val, part in enumerate(atomistic_protein_centered):
                if residue in atomistic_protein_centered[part]:
                    exists=True
                    for atom in atomistic_protein_centered[part][residue]:   
                        if atomistic_protein_centered[part][residue][atom]['res_type'] != resname:
                            print('de_novo' , resname,'at_user', atomistic_protein_centered[part][residue][atom]['res_type'])
                            sys.exit('de novo and at user supplied don\'t match')
                        atomistic_protein_centered[part][residue][atom]['coord'] = at_mod.rotate_atom(atomistic_protein_centered[part][residue][atom]['coord'], cg_com[part_val], xyz_rot_apply[part_val])
                    complete_user_at[residue]=atomistic_protein_centered[part][residue]
        if not exists:
            complete_user_at[residue]=final_coordinates_atomistic[residue]
    return complete_user_at

def write_user_chains_to_pdb(atomistic_user_supplied, chain):
    if not os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN_aligned_'+str(chain)+'.pdb'):
        pdb_output = gen.create_pdb(g_var.working_dir+'PROTEIN/PROTEIN_aligned_'+str(chain)+'.pdb')
        final_atom={}
        at_id=0
        coord=[]
        for resid in atomistic_user_supplied:
            for atom in atomistic_user_supplied[resid]:
                coord.append(atomistic_user_supplied[resid][atom]['coord'])
                short_line=atomistic_user_supplied[resid][atom]
                final_atom[at_id]={'atom':short_line['atom'], 'res_type':short_line['res_type'], 'chain':chain, 'residue':resid,\
                                    'x':short_line['coord'][0],'y':short_line['coord'][1],'z':short_line['coord'][2]}
                at_id+=1
        merge_coords=at_mod.check_atom_overlap(coord)

        for at_id, coord in enumerate(merge_coords):
            x, y, z = gen.trunc_coord(coord)
            pdb_output.write(g_var.pdbline%((at_id+1,final_atom[at_id]['atom'],final_atom[at_id]['res_type'],'A',final_atom[at_id]['residue'],
                 x,y,z,1,0))+'\n') 

################################################################## Merge chains ####################

def merge_protein_pdbs(file, end):
#### reads in each chain into merge list
    merge, merged_coords = [],[]
    count = 0
    for chain in range(0,g_var.system['PROTEIN']):
        merge_temp = []
        if os.path.exists(file+'_'+str(chain)+end):
            with open(file+'_'+str(chain)+end, 'r') as pdb_input:
                for line in pdb_input.readlines():
                    if line.startswith('ATOM'):
                        line_sep=gen.pdbatom(line)
                        merge_temp.append(line_sep)
        else:
            sys.exit('cannot find minimised protein chain: '+file+'_'+str(chain)+end)
        if 'PROTEIN_aligned' in file and '_gmx_checked.pdb' in end:  
            count += write_disres(merge_temp, chain, file, count)
        merge, merge_coords = at_mod.fix_chirality(merge,merge_temp,merged_coords, 'PROTEIN')   
    merged_coords = at_mod.check_atom_overlap(merge_coords)
    merged=[]
    for line_val, line in enumerate(merge):
        x, y, z = gen.trunc_coord(merged_coords[line_val])
        merged.append(g_var.pdbline%((int(line['atom_number']), line['atom_name'], line['residue_name'],' ',line['residue_id'],\
            x,y,z,1,0))+'\n')
    if 'aligned' in file.split('/')[-1]:
        write_merged_pdb(merged, '_aligned')
    else:
        write_merged_pdb(merged,'_de_novo')
    

def correct_amide_h(lines, coords):
    for at_val, atom in enumerate(lines):
        if atom['atom_name'] == 'C': 
            C = np.array(coords[at_val])
        if atom['atom_name'] == 'O':   
            O = np.array(coords[at_val]) 
        if atom['atom_name'] == f_loc.res_top[atom['residue_name']]['amide_h'] and atom['residue_id'] != 0:  
            HN = np.array(coords[at_val]) 
            HN_index = at_val
        if atom['atom_name'] == 'N' and atom['residue_id'] != 0:
            N = np.array(coords[at_val])
        if 'C' in locals() and 'N' in locals() and 'O' in locals() and 'HN' in locals():
            O_C = O-C
            coords[at_val] = N - O_C
            del N, C, O, HN 
    return coords

def write_disres(coord, chain, file, at_start):
    if not os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN_disres.itp'):
        with open(g_var.working_dir+'PROTEIN/PROTEIN_disres.itp', 'a') as disres_out:
            disres_out.write(';backbone hydrogen bonding distance restraints\n\n')
            disres_out.write('[ intermolecular_interactions ]\n[ distance_restraints ]\n')
            disres_out.write(';   i     j type label      funct         lo        up1        up2     weight')
            HN, O = [],[]
            for atom in coord:
                if atom['residue_name'] in f_loc.p_residues and atom['atom_name'] == f_loc.res_top[atom['residue_name']]['amide_h']:
                    HN.append([int(atom['atom_number']), atom['x'], atom['y'], atom['z']])
                if atom['residue_name'] in f_loc.p_residues and atom['atom_name'] == 'O':
                    O.append([int(atom['atom_number']), atom['x'], atom['y'], atom['z']])
            HN, O = np.array(HN), np.array(O)
            tree = cKDTree(HN[:,1:])
            count=0
            for carbonyl in O:
                ndx = tree.query_ball_point(carbonyl[1:], r=3)
                if len(ndx) > 0:
                    for at in ndx:
                        if coord[int(HN[ndx[0]][0])-1]['residue_id'] < coord[int(carbonyl[0])-1]['residue_id']-1 or coord[int(HN[ndx[0]][0])-1]['residue_id'] > coord[int(carbonyl[0])-1]['residue_id']+1:
                            count+=1
                            xyz1 = [coord[int(carbonyl[0])-1]['x'], coord[int(carbonyl[0])-1]['y'], coord[int(carbonyl[0])-1]['z']]
                            xyz2 = [coord[int(HN[at][0])-1]['x'], coord[int(HN[at][0])-1]['y'], coord[int(HN[at][0])-1]['z']]
                            dist = np.round((gen.calculate_distance(xyz1, xyz2)/10)-0.05, 4)
                            disres_out.write('\n{0:10}{1:10}{2:3}{3:12}{4:12}{5:^12}{6:14}{7:14}{8:5}'.format(str(at_start+int(HN[at][0])), str(at_start+int(carbonyl[0])), 
                                                                                                    '1', str(count),'1', '0', str(dist), str(dist), '5'))
    return len(coord) 

def write_merged_pdb(merge, protein):
#### creates merged pdb and writes chains to it
    if not os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN'+protein+'_merged.pdb'):
        pdb_output=gen.create_pdb(g_var.working_dir+'PROTEIN/PROTEIN'+protein+'_merged.pdb')
        for line in merge:
            pdb_output.write(line)
        pdb_output.close()

def align_user_chains(atomistic_protein_input, group_chain, final_coordinates_atomistic):
    atomistic_protein_centered, cg_com = center_atomistic(atomistic_protein_input, group_chain) ## centers each chain by center of mass
    at_com_group, cg_com_group = rotate_protein_monomers(atomistic_protein_centered, final_coordinates_atomistic, cg_com, group_chain) 
    atomistic_protein_rotated = apply_rotations_to_chains(final_coordinates_atomistic, atomistic_protein_centered, 
                                                                at_com_group,cg_com_group,cg_com, group_chain) ## apply rotation matrix to atoms and build in missing residues
    final_user_supplied_coord = correct_disulphide_bonds(atomistic_protein_rotated) ## fixes sulphur distances in user structure

    pool = mp.Pool(g_var.ncpus)
    pool_process = pool.starmap_async(write_user_chains_to_pdb, [(final_user_supplied_coord[chain], chain) ## write structure to pdb
                                        for chain in final_user_supplied_coord]).get()
    pool.close()


########################################################### RMSD

def write_RMSD():
    RMSD={}
    de_novo_atoms, chain_count = read_in.read_in_atomistic(g_var.final_dir+'final_cg2at_de_novo.pdb', False) ## reads in final pdb
    if chain_count != g_var.system['PROTEIN']:
        sys.exit('number of chains in atomistic protein input ('+str(chain_count)+') does not match CG representation ('+str(g_var.system['PROTEIN'])+')')
    RMSD['de novo '] = RMSD_measure(de_novo_atoms) ## gets rmsd of de novo

    if g_var.user_at_input and 'PROTEIN' in g_var.cg_residues: 
        if g_var.o in ['all', 'align']: 
            at_input_atoms, chain_count = read_in.read_in_atomistic(g_var.final_dir+'final_cg2at_aligned.pdb', False)
            RMSD['at aligned'] = RMSD_measure(at_input_atoms)   
    with open(g_var.final_dir+'structure_quality.dat', 'w') as qual_out:   
        qual_out.write('\n{0:^10}{1:^25}{2:^10}\n'.format('output ','chain','RMSD ('+chr(197)+')'))
        qual_out.write('{0:^10}{1:^25}{2:^10}\n'.format('-------','-----','---------'))
        print('\n{0:^10}{1:^25}{2:^10}'.format('output ','chain','RMSD ('+chr(197)+')'))
        print('{0:^10}{1:^25}{2:^10}'.format('-------','-----','---------'))
        for rmsd in RMSD:
            for chain in RMSD[rmsd]:
                qual_out.write('{0:^10}{1:^25}{2:^10}\n'.format(rmsd, str(chain), float(RMSD[rmsd][chain])))
                print('{0:^10}{1:^25}{2:^10}'.format(rmsd, str(chain), float(RMSD[rmsd][chain])))
        print()

def RMSD_measure(structure_atoms):
    RMSD_dict = {}
    for chain in range(g_var.system['PROTEIN']):
        at_centers=[]
    #### runs through every residue and atom  
        for residue in structure_atoms[chain]:
        #### gets center of mass of each residue (note only backbone heavy atoms have a mass)
            at_centers_iter=[]
            for atom in structure_atoms[chain][residue]:
                at_centers_iter.append(np.append(structure_atoms[chain][residue][atom]['coord'],structure_atoms[chain][residue][atom]['frag_mass']))
            try:
                at_centers.append(np.average(np.array(at_centers_iter)[:,:3], axis=0, weights=np.array(at_centers_iter)[:,3]))
            except:
                print('The fragment probably has no mass\n')
                for atom in structure_atoms[chain][residue]:
                    print(structure_atoms[chain][residue][atom])
                sys.exit()
    #### checks that the number of residues in the chain are the same between CG and AT
        if len(at_centers) != len(g_var.backbone_coords[chain]):
            sys.exit('In chain '+str(chain)+' the atomistic input does not match the CG. \n\
    number of CG residues '+str(len(g_var.backbone_coords[chain]))+'\nnumber of AT residues '+str(len(at_centers)))
        cg_center = np.mean(np.array(g_var.backbone_coords[chain])[:,:3], axis=0)
        at_align = np.array(at_centers) - (np.mean(np.array(at_centers), axis=0) - cg_center)
        xyz_rot_apply=at_mod.kabsch_rotate(np.array(at_align)-cg_center, np.array(np.array(g_var.backbone_coords[chain])[:,:3])-cg_center)
        for at_val, atom in enumerate(at_align):
            at_align[at_val] = at_mod.rotate_atom(atom, cg_center, xyz_rot_apply) 
        #### finds distance between backbone COM and cg backbone beads
        dist=np.sqrt((np.array(at_align) - np.array(g_var.backbone_coords[chain])[:,:3])**2)
        RMSD_val = np.sqrt(np.mean(dist**2)) #### RMSD calculation
        RMSD_dict[chain]=np.round(RMSD_val, 3)  #### stores RMSD in dictionary
    return RMSD_dict

