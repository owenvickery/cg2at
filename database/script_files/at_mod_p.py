#!/usr/bin/env python3

import os, sys
import numpy as np
import copy
from string import ascii_uppercase
import difflib
from scipy.spatial import cKDTree
import gen, g_var, f_loc, at_mod
import math

def build_protein_atomistic_system(cg_residues, box_vec, user_supplied):
#### initisation of counters
    chain_count=0
    system={}
    backbone_coords={}
    backbone_coords[chain_count]=[]
    terminal={}
    terminal[chain_count]=[]
    coordinates_atomistic={}
    coordinates_atomistic[chain_count]={}
    sequence={}
    sequence[chain_count]=[]
    print('Converting Protein')
    gen.mkdir_directory(g_var.working_dir+'PROTEIN')  ### make and change to protein directory
#### for each residue in protein
    initial=True
    residue_type={}
    residue_type_mass={}
    for cg_residue_id, residue_number in enumerate(cg_residues):
        resname = cg_residues[residue_number][next(iter(cg_residues[residue_number]))]['residue_name']
        if cg_residue_id == 0:
            terminal[chain_count].append(f_loc.backbone[resname]['ter'])
        coordinates_atomistic[chain_count][residue_number]={}
        frag_location=at_mod.fragment_location(resname) ### get fragment location from database
        residue_type[resname], residue_type_mass[resname] = at_mod.get_atomistic(frag_location)
        for group in residue_type[resname]:
            center, at_frag_centers, cg_frag_centers, group_fit = at_mod.rigid_fit(residue_type[resname][group], residue_type_mass[resname], residue_number, cg_residues[residue_number])

            at_connect, cg_connect = at_mod.connectivity(cg_residues[residue_number], at_frag_centers, cg_frag_centers, group_fit, group)
            if 'BB' in group_fit:
                BB_connect = []
                for atom in group_fit['BB']:
                    if group_fit['BB'][atom]['atom'] == f_loc.backbone[resname]['N_ter']:
                        N_ter=atom
                    if group_fit['BB'][atom]['atom'] == f_loc.backbone[resname]['C_ter']:
                        C_ter=atom
                at_connect, cg_connect, new_chain = BB_connectivity(at_connect,cg_connect, cg_residues, group_fit['BB'], residue_number, N_ter, C_ter)
                sequence = at_mod.add_to_sequence(sequence, resname, chain_count)
                backbone_coords[chain_count].append(np.append(cg_residues[residue_number]['BB']['coord'], 1))

            if len(at_connect) == len(cg_connect):
                xyz_rot_apply=at_mod.rotate(np.array(at_connect)-center, np.array(cg_connect)-center, False)
            else:
                print('atom connections: '+str(len(at_connect))+' does not equal CG connections: '+str(len(cg_connect)))
                sys.exit('residue number: '+str(residue_number)+', residue type: '+str(resname)+', group: '+group)

            for bead in group_fit:
                for atom in group_fit[bead]:
                    group_fit[bead][atom]['coord'] = at_mod.rotate_atom(group_fit[bead][atom]['coord'], center, xyz_rot_apply)   
                    atom_new = group_fit[bead][atom].copy()
                    coordinates_atomistic[chain_count][residue_number][atom] = atom_new
        if new_chain:
            terminal[chain_count].append(f_loc.backbone[resname]['ter'])
            chain_count+=1
            if cg_residue_id+1 != len(cg_residues):
                backbone_coords[chain_count]=[]
                coordinates_atomistic[chain_count]={}
                terminal[chain_count]=[]
                terminal[chain_count].append(f_loc.backbone[cg_residues[residue_number+1]['BB']['residue_name']]['ter'])
                sequence[chain_count]=[]
    if g_var.v >=1:
        print('\n{0:^15}{1:^12}'.format('chain number', 'length of chain')) #   \nchain number\tDelta A\t\tno in pdb\tlength of chain')
        print('\n{0:^15}{1:^12}'.format('------------', '---------------'))
        for chain in sequence:
            print('{0:^15}{1:^12}'.format(chain, len(sequence[chain])))
        print()
    
    system['terminal_residue']=terminal
    system['PROTEIN']=chain_count
    return system, backbone_coords, coordinates_atomistic, sequence    


def BB_connectivity(at_connections,cg_connections, cg_residues, at_residues, residue_number, N_ter, C_ter):
#### connect to preceeding backbone bead in chain
    new_chain=False
    try:
        xyz_cur = cg_residues[residue_number]['BB']['coord']
        xyz_prev = cg_residues[residue_number-1]['BB']['coord']
        dist=gen.calculate_distance(xyz_prev, xyz_cur)
        if dist < 5:
            cg_n = cg_residues[residue_number-1]['BB']['coord']
            at_n = at_residues[N_ter]['coord']
            cg_connections.append(cg_n)
            at_connections.append(at_n)
    except:
        pass
#### connect to next backbone bead in chain
    try:
        xyz_cur = cg_residues[residue_number]['BB']['coord']
        xyz_next = cg_residues[residue_number+1]['BB']['coord']
        dist=gen.calculate_distance(xyz_next, xyz_cur)
        if dist < 5:
            cg_c = cg_residues[residue_number+1]['BB']['coord']
            at_c = at_residues[C_ter]['coord']
            cg_connections.append(cg_c)
            at_connections.append(at_c)
        else:
            new_chain=True
    except:
        new_chain=True
        pass
    return at_connections,cg_connections, new_chain

################# Fixes disulphide bond, martini cysteine bone is too far apart to be picked up by pdb2gmx. 
#### 

def ask_if_disulphide(chain, res_1, res_2):
    while True:
        try:
            answer = str(input('\nAre these residues connected by a disulphide bond in chain '+str(chain)+' (Y/N): '+str(res_1)+'--'+str(res_2)+': '))
            print(answer)
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
    cys_seq={}
    for chain in atomistic_protein_input:
        cys_seq[chain]=[]
        cysteines = []
        cys_resid = []
        for part in atomistic_protein_input[chain]:
            for resid in atomistic_protein_input[chain][part]:
                atom = next(iter(atomistic_protein_input[chain][part][resid]))
                if atomistic_protein_input[chain][part][resid][atom]['res_type'] == 'CYS':
                    for atom in atomistic_protein_input[chain][part][resid]:
                        if atomistic_protein_input[chain][part][resid][atom]['atom'] == f_loc.backbone['CYS']['sul']:
                            cys_resid.append(resid)
                            cysteines.append(atomistic_protein_input[chain][part][resid][atom]['coord'])
        if len(cysteines)>=2:
            tree = cKDTree(cysteines)
            done_query=[]
            for cys_index, cys in enumerate(cysteines):
                query = tree.query_ball_point(cys, r=2.1)
                if len(query) == 2 and query not in done_query:
                    cys_seq[chain].append([cys_resid[query[0]],cys_resid[query[1]]])
                    done_query.append(query)
    return cys_seq

def find_disulphide_bonds_de_novo(atomistic_protein_input, user_cys_bond):
    for chain in atomistic_protein_input:
        offset = next(iter(atomistic_protein_input[chain]))
        if chain not in user_cys_bond:
            user_cys_bond[chain]=[]
        cysteines = []
        cys_resid = []
        for resid in atomistic_protein_input[chain]:
            resid_corrected=np.array(resid)-offset
            run = True
            if len(user_cys_bond) > 0:
                for bonds in user_cys_bond[chain]:
                    if resid_corrected in bonds:
                        run=False
            if run:
                atom = next(iter(atomistic_protein_input[chain][resid]))
                if atomistic_protein_input[chain][resid][atom]['res_type'] == 'CYS':
                    for atom in atomistic_protein_input[chain][resid]:
                        if atomistic_protein_input[chain][resid][atom]['atom'] == f_loc.backbone['CYS']['sul']:
                            cys_resid.append(resid_corrected)
                            cysteines.append(atomistic_protein_input[chain][resid][atom]['coord'])
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
                            user_cys_bond[chain].append([cys_resid[query[0]],cys_resid[query[1]]])
                    done_query.append(query)
    return user_cys_bond

def correct_disulphide_bonds(coordinates_atomistic, user_cys_bond):
    for chain in coordinates_atomistic:
        offset = next(iter(coordinates_atomistic[chain]))
        for bonds in user_cys_bond[chain]:
            bonds_corrected=np.array(bonds)+offset
            at_num=[]
            at_coord=[]
            for resid in bonds_corrected: 
                for atom in coordinates_atomistic[chain][resid]:
                    if coordinates_atomistic[chain][resid][atom]['atom'] == f_loc.backbone['CYS']['sul']:
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

def finalise_novo_atomistic(atomistic, cg_residues, box_vec):
    final_at_residues={}
    final_at = {}
    for chain in atomistic: 
        at_number=0    
        final_at_residues[chain]={}  
        final_at[chain]={}
        coords=[]
        skip = os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN_de_novo_'+str(chain)+'.pdb')
        if not skip:
            pdb_output = gen.create_pdb(g_var.working_dir+'PROTEIN/PROTEIN_de_novo_'+str(chain)+'.pdb', box_vec)
        for res_index, residue_id in enumerate(atomistic[chain]):
            if atomistic[chain][residue_id][1]['res_type'] in f_loc.mod_residues:
                atomistic[chain][residue_id] = at_mod.check_hydrogens(atomistic[chain][residue_id])
            if res_index < len(atomistic[chain])-2:
                atomistic[chain][residue_id] = fix_carbonyl(residue_id, cg_residues, atomistic[chain][residue_id])
            else:
                atomistic[chain][residue_id] = atomistic[chain][residue_id]
            for at_val, atom in enumerate(atomistic[chain][residue_id]):
                atomistic[chain][residue_id][at_val+1]['resid'] = res_index
                coords.append(atomistic[chain][residue_id][at_val+1]['coord'])
                final_at[chain][at_number]=atomistic[chain][residue_id][at_val+1]
                at_number+=1
            final_at_residues[chain][res_index]=atomistic[chain][residue_id]

        coords = at_mod.check_atom_overlap(coords)
        for atom in final_at[chain]:
            final_at[chain][atom]['coord']=coords[atom]
            final_at_residues[chain][final_at[chain][atom]['resid']][atom]=final_at[chain][atom]
            if not skip:
                pdb_output.write(g_var.pdbline%((atom,final_at[chain][atom]['atom'],final_at[chain][atom]['res_type'],' ',\
                                final_at[chain][atom]['resid'],final_at[chain][atom]['coord'][0],final_at[chain][atom]['coord'][1],final_at[chain][atom]['coord'][2],1,0))+'\n')
    if 'pdb_output' in locals():
        pdb_output.close()
    return final_at_residues

########################################### fix carbonyl section 

def fix_carbonyl(residue_id, cg, at):
    ca=[]
    for index in range(3):
        ca.append(cg[residue_id+index]['BB']['coord'])

    for atom in at:
        if at[atom]['atom'] == 'C': 
            C = atom
        if at[atom]['atom'] in 'O':   
            O = atom                 

    initial_vector, cross_vector = at_mod.find_cross_vector( ca, at[C]['coord'], at[O]['coord'])
    rotation = at_mod.align_to_vector(initial_vector, cross_vector)
    center = ca[0]+(ca[1]-ca[0])/3
    at[C]['coord'] = (at[C]['coord']-center).dot(rotation)+center
    at[O]['coord'] = (at[O]['coord']-center).dot(rotation)+center
    return at

##################################################################  User supplied protein ##############

def read_in_atomistic(protein):
#### reset location and check if pdb exists  
    os.chdir(g_var.start_dir)
    if not os.path.exists(protein):
        sys.exit('cannot find atomistic protein : '+protein)
#### read in atomistic fragments into dictionary residue_list[0]=x,y,z,atom_name    
    atomistic_protein_input={}
    chain_count=0
#### read in pdb
    ter_residues=[]
    with open(protein, 'r') as pdb_input:
        atomistic_protein_input[chain_count]={}
        for line_nr, line in enumerate(pdb_input.readlines()):
            #### separate line 
            run=False ## turns to true is line is a bead/atom
            if line.startswith('ATOM'):
                line_sep = gen.pdbatom(line)
                # print(line_sep['atom_name'])
                if line_sep['residue_name'] in f_loc.mod_residues:
                    run=True
                if not gen.is_hydrogen(line_sep['atom_name']):
                    run=True
            #### if line is correct
            if run:
                if line_sep['residue_name'] in f_loc.p_residues:
                    if not gen.is_hydrogen(line_sep['atom_name']) or line_sep['residue_name'] in f_loc.mod_residues:  
                    #### sorts out wrong atoms in terminal residues
                        if line_sep['atom_name'] in ['OT', 'O1', 'O2']:
                            line_sep['atom_name']='O'
                    #### makes C_terminal connecting atom variable  
                        if line_sep['atom_name'] == f_loc.backbone[line_sep['residue_name']]['C_ter']:
                            C_ter=[line_sep['x'],line_sep['y'],line_sep['z']]
                            C_resid=line_sep['residue_id']
                            C=True
                        try:
                        #### tries to make a N_terminal connecting atom variable
                            if line_sep['atom_name'] == f_loc.backbone[line_sep['residue_name']]['N_ter']:
                                N_resid=line_sep['residue_id']
                                N_ter=[line_sep['x'],line_sep['y'],line_sep['z']]
                                N=True
                        #### measures distance between N and C atoms. if the bond is over 3 A it counts as a new protein
                            dist=gen.calculate_distance(N_ter, C_ter)
                            if N and C and C_resid != N_resid and dist > 3.5:# and aas[line_sep['residue_name']] != sequence[chain_count][line_sep['residue_id']]:
                                N_ter, C_ter=False, False
                                ter_residues.append(line_sep['residue_id'])
                                chain_count+=1
                                atomistic_protein_input[chain_count]={} ### new chain key
                        except:
                            pass
                        if line_sep['residue_id'] not in atomistic_protein_input[chain_count]:  ## if protein does not exist add to dict
                            atomistic_protein_input[chain_count][line_sep['residue_id']]={}
                    #### adds atom to dictionary, every atom is given a initial mass of zero 
                        atomistic_protein_input[chain_count][line_sep['residue_id']][line_sep['atom_number']]={'coord':np.array([line_sep['x'],line_sep['y'],line_sep['z']]),'atom':line_sep['atom_name'], 'res_type':line_sep['residue_name'],'frag_mass':0, 'resid':line_sep['residue_id']}
                    #### if atom is in the backbone list then its mass is updated to the correct one
                        if line_sep['atom_name'] in f_loc.backbone[line_sep['residue_name']]['atoms']:
                            for atom in line_sep['atom_name']:
                                if atom in g_var.mass:
                                    atomistic_protein_input[chain_count][line_sep['residue_id']][line_sep['atom_number']]['frag_mass']=g_var.mass[atom]
                
                    

    return atomistic_protein_input, chain_count+1      

def check_sequence(atomistic_protein_input, chain_count):
    at={}
    seq_user={}
    for chain in range(len(atomistic_protein_input)):
        seq_user[chain]=[]
        for resid in atomistic_protein_input[chain]:
            for atom in atomistic_protein_input[chain][resid]:
                seq_user = at_mod.add_to_sequence(seq_user, atomistic_protein_input[chain][resid][atom]['res_type'], chain)
                break
    return seq_user

def align_chains(atomistic_protein_input, seq_user, sequence):
    if g_var.v >= 1:
        print('coarse grain protein sequence:\n')
        for index in sequence:
            print('chain:', index,sequence[index], '\n')
        print('\nuser supplied structure:\n')
        for index in seq_user:
            print('chain:', index, seq_user[index], '\n')        
    at={}
    sequence_temp = sequence.copy()
    test_chain={}
    for chain_at in range(len(atomistic_protein_input)):
        skip_sequence=False
        chain_cg=0
        s = difflib.SequenceMatcher(None, seq_user[chain_at], sequence[chain_cg])
        seq_info = s.get_matching_blocks()
        while seq_info[0][2] != len(seq_user[chain_at]):
            if chain_cg >= len(sequence)-1:
                print('\nCannot find a match for user supplied chain: '+str(chain_at)+'\n'+str(seq_user[chain_at]),'\nIn\n'+str(sequence))
                print('Using de novo instead')
                skip_sequence = True
                break
            if not skip_sequence:
                chain_cg+=1
                s = difflib.SequenceMatcher(None, seq_user[chain_at], sequence[chain_cg])
                seq_info = s.get_matching_blocks()
        if not skip_sequence:
            temp={}
            if chain_cg not in at:
                at[chain_cg]={}
            if seq_info[0][2] == len(seq_user[chain_at]):
                for resid,  residue in enumerate(atomistic_protein_input[chain_at]):
                    temp[resid + seq_info[0][1]] = atomistic_protein_input[chain_at][residue]
                at[chain_cg][str(seq_info[0][1])+':'+str(seq_info[0][1]+seq_info[0][2])]=temp  
            sequence[chain_cg] = mask_sequence(sequence[chain_cg], seq_info[0][1], seq_info[0][1]+seq_info[0][2])
            test_chain[chain_at]=chain_cg

    if f_loc.group_chains == 'chain':
        return at, test_chain
    else:
        return at, f_loc.group_chains

def mask_sequence(sequence, st, end):
    for index, residue in enumerate(sequence):
        if st <= index < end:
            sequence[index]='-'
    return sequence

def center_atomistic(atomistic_protein_input, backbone_coords, group_chain):
    cg_com={}
#### for each protein chain center on cg representation 
    for chain in atomistic_protein_input:
        cg_com[chain]=[]
        for part_val, part in enumerate(atomistic_protein_input[chain]):
            sls, sle= int(part.split(':')[0]),int(part.split(':')[1])
            if group_chain==None:
                protein_mass=fetch_backbone_mass(atomistic_protein_input[chain][part], [])
                atomistic_protein_mass = at_mod.COM(protein_mass, 'protein at: '+str(chain)+' '+part)#### add center of mass of CG_proteins
                cg_com[chain].append(at_mod.COM(backbone_coords[chain][sls:sle], 'protein cg: '+str(chain)+' '+part))
                atomistic_protein_input[chain][part] = update_part_coordinate(atomistic_protein_input[chain][part], atomistic_protein_mass, cg_com[chain][part_val])
            elif group_chain=='all':
                if 'protein_mass' not in locals():
                    protein_mass=[]
                protein_mass=fetch_backbone_mass(atomistic_protein_input[chain][part], protein_mass)
                if 'cg_backbone_masses' not in locals():
                    cg_backbone_masses = {}
                    cg_backbone_masses['all'] = backbone_coords[chain][sls:sle]
                else:
                    cg_backbone_masses['all'] += backbone_coords[chain][sls:sle]
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
                        cg_backbone_masses[group_chain[chain]]= backbone_coords[chain][sls:sle]
                    else:
                        cg_backbone_masses[group_chain[chain]] += backbone_coords[chain][sls:sle]                  
    if group_chain=='all':
        atomistic_protein_input, cg_com = center_at_protein_all_chains(atomistic_protein_input, cg_com, protein_mass, cg_backbone_masses, group_chain)
    if group_chain not in ['all',None]:
        atomistic_protein_input, cg_com = center_at_protein_chain_groups(atomistic_protein_input, cg_com, protein_mass, cg_backbone_masses, backbone_coords, group_chain)
    return atomistic_protein_input, cg_com

def center_at_protein_chain_groups(atomistic_protein_input, cg_com, protein_mass, cg_backbone_masses, backbone_coords, group_chain):
    for chain in atomistic_protein_input:
        for part_val, part in enumerate(atomistic_protein_input[chain]):
            sls, sle= int(part.split(':')[0]),int(part.split(':')[1])
            if chain in group_chain:
                atomistic_protein_mass = at_mod.COM(protein_mass[group_chain[chain]], 'AT protein chain: '+str(chain))
                cg_com[chain].append(at_mod.COM(cg_backbone_masses[group_chain[chain]], 'CG protein chain: '+str(chain)))
            else:
                protein_mass=fetch_backbone_mass(atomistic_protein_input[chain][part], [])
                atomistic_protein_mass = at_mod.COM(protein_mass, 'protein at: '+str(chain)+' '+part)#### add center of mass of CG_proteins
                cg_com[chain].append(at_mod.COM(backbone_coords[chain][sls:sle], 'protein cg: '+str(chain)+' '+part))
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
            if part[residue][atom]['atom'] in f_loc.backbone[part[residue][atom]['res_type']]['atoms']:
                protein_mass.append([part[residue][atom]['coord'][0],part[residue][atom]['coord'][1],part[residue][atom]['coord'][2],part[residue][atom]['frag_mass']])    
    return protein_mass

def rotate_protein_monomers(atomistic_protein_centered, final_coordinates_atomistic, backbone_coords, cg_com,  box_vec, group_chain):
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
                        if atomistic_protein_centered[chain][part][residue][atom]['atom'] in f_loc.backbone[atomistic_protein_centered[chain][part][residue][atom]['res_type']]['atoms']:
                            at_centers_iter.append(np.append(atomistic_protein_centered[chain][part][residue][atom]['coord'],atomistic_protein_centered[chain][part][residue][atom]['frag_mass']))
                    try:
                        at_centers.append(np.average(np.array(at_centers_iter)[:,:3], axis=0, weights=np.array(at_centers_iter)[:,3]))
                    except:
                        for atom in atomistic_protein_centered[chain][part][residue]:
                            print(atomistic_protein_centered[chain][part][residue][atom])
                        sys.exit()
                if len(at_centers) == len(np.array(backbone_coords[chain])[sls:sle,:3]):
                    if group_chain=='all':
                        at_com_group, cg_com_group = return_all_rotations(chain, at_centers, backbone_coords, at_com_group, cg_com_group, sls, sle)
                    elif group_chain==None:
                        at_com_group = return_indivdual_rotations(chain, part_val, at_centers, backbone_coords, cg_com, at_com_group, cg_com_group, sls, sle)
                    else:
                        at_com_group, cg_com_group = return_grouped_rotations(chain, part_val, at_centers, backbone_coords, cg_com, at_com_group, cg_com_group, sls, sle,group_chain)
                else:
                    sys.exit('In chain '+str(chain)+' the atomistic input does not match the CG. \n\
                            number of CG residues '+str(len(backbone_coords[chain]))+'\nnumber of AT residues '+str(len(at_centers)))
    return at_com_group, cg_com_group

def apply_rotations_to_chains(final_coordinates_atomistic, atomistic_protein_centered, at_com_group,cg_com_group,cg_com, box_vec,group_chain):
    if group_chain=='all':
        rotate_all = return_all_rotations_final(at_com_group,cg_com_group,cg_com)
    final_rotations = []
    final_user_supplied_coord={}
    for chain in range(len(final_coordinates_atomistic)):
        rotations=[]
        if chain in atomistic_protein_centered:
            if group_chain not in ['all',None]: 
                if chain in group_chain:
                    rotate_chain=at_mod.rotate(at_com_group[group_chain[chain]]-cg_com[group_chain[chain]][0], cg_com_group[group_chain[chain]]-cg_com[group_chain[chain]][0], False)
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
                if g_var.v >= 1:
                    if 'seg' not in locals():
                        seg=0
                    else:
                        seg+=1
                    print_rotations(cg_com, chain, part_val, rotations,group_chain, seg)

            final_user_supplied_coord[chain] = hybridise_protein_inputs(final_coordinates_atomistic[chain], atomistic_protein_centered[chain], cg_com[chain], rotations, chain, box_vec)
        else:
            final_user_supplied_coord[chain] = hybridise_protein_inputs(final_coordinates_atomistic[chain], [], [], [], chain, box_vec, user_cys_bond)
    return final_user_supplied_coord

def return_all_rotations_final(at_com_group,cg_com_group,cg_com):
    if len(at_com_group['all']) == len(cg_com_group['all']):
        return at_mod.rotate(at_com_group['all']-cg_com[0][0], cg_com_group['all']-cg_com[0][0], False)
    else:
        sys.exit('The atomistic input does not match the CG. \n\
                number of CG residues '+str(len(cg_com_group['all']))+'\nnumber of AT residues '+str(len(at_com_group['all'])))

def return_indivdual_rotations(chain, part_val, at_centers, backbone_coords, cg_com, at_com_group, cg_com_group, sls, sle):
    if chain not in at_com_group:
        at_com_group[chain]=[]
## finds optimal rotation of each monomer  
    if len(at_centers) == len(np.array(backbone_coords[chain])[sls:sle,:3]):
        at_com_group[chain].append(at_mod.rotate(np.array(at_centers)-cg_com[chain][part_val], 
                             np.array(backbone_coords[chain])[sls:sle,:3]-cg_com[chain][part_val], False))
    else:
        sys.exit('In chain '+str(chain)+' the atomistic input does not match the CG. \n\
                    number of CG residues '+str(len(backbone_coords[chain]))+'\nnumber of AT residues '+str(len(at_centers)))    
    return at_com_group

def return_grouped_rotations(chain, part_val, at_centers, backbone_coords, cg_com, at_com_group, cg_com_group, sls, sle, group_chain):
    if chain in group_chain:
        if chain not in at_com_group:
            at_com_group[group_chain[chain]]=[]   
            at_com_group[group_chain[chain]]=np.array(at_centers)
            cg_com_group[group_chain[chain]]=np.array(backbone_coords[chain])[sls:sle,:3]  
        else:
            at_com_group[group_chain[chain]]=np.append(at_com_group[group_chain[chain]], np.array(at_centers), axis=0)
            cg_com_group[group_chain[chain]]=np.append(cg_com_group[group_chain[chain]], np.array(backbone_coords[chain])[sls:sle,:3], axis=0)                                                                               
    else:
        at_com_group = return_indivdual_rotations(chain, part_val, at_centers, backbone_coords, cg_com, at_com_group, cg_com_group, sls, sle)
    return at_com_group, cg_com_group

def return_all_rotations(chain, at_centers, backbone_coords, at_com_group, cg_com_group, sls, sle):
    if 'all' not in at_com_group:
        at_com_group['all']=np.array(at_centers)
        cg_com_group['all']=np.array(backbone_coords[chain])[sls:sle,:3]
    else:
        at_com_group['all']=np.append(at_com_group['all'], np.array(at_centers), axis=0)
        cg_com_group['all']=np.append(cg_com_group['all'], np.array(backbone_coords[chain])[sls:sle,:3], axis=0) 
    return at_com_group, cg_com_group   

def print_rotations(cg_com, chain, part_val, xyz_rot_apply,group_chain, seg):
    if seg == 0:
        print('\nThe proteins chains are rotated around the COM of backbone heavy atoms.\n')
    if group_chain != None:
        if group_chain == 'all':
            group = [g for g in cg_com]
            print('User supplied AT chain '+str(chain)+' is grouped together with chains: '+', '.join(map(str, group)))
        elif seg in group_chain:
            group = [g for g in group_chain if group_chain[g] == group_chain[seg]]
            print('User supplied AT chain '+str(seg)+' is grouped together with chains: '+', '.join(map(str, group))+'\n')
        else:
            print('User supplied AT chain '+str(seg)+' is treated separately')
        
    print('The COM of chain', chain,'is :', np.round(cg_com[chain][part_val][0], 2),',', np.round(cg_com[chain][part_val][1], 2),',', 
          np.round(cg_com[chain][part_val][2], 2))
    print('rotating chain ', chain, 'by :',np.round(np.degrees(xyz_rot_apply[part_val][0]),2),',',np.round(np.degrees(xyz_rot_apply[part_val][1]),2),
          ',',np.round(np.degrees(xyz_rot_apply[part_val][2]),2))
    print()    

def hybridise_protein_inputs(final_coordinates_atomistic, atomistic_protein_centered, cg_com, xyz_rot_apply, chain, box_vec):

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

def write_user_chains_to_pdb(atomistic_user_supplied, box_vec):
    for chain in atomistic_user_supplied:
        pdb_output = gen.create_pdb(g_var.working_dir+'PROTEIN/PROTEIN_aligned_'+str(chain)+'.pdb', box_vec)
        final_atom={}
        at_id=0
        coord=[]
        for resid in atomistic_user_supplied[chain]:
            for atom in atomistic_user_supplied[chain][resid]:
                coord.append(atomistic_user_supplied[chain][resid][atom]['coord'])
                short_line=atomistic_user_supplied[chain][resid][atom]
                final_atom[at_id]={'atom':short_line['atom'], 'res_type':short_line['res_type'], 'chain':ascii_uppercase[chain], 'residue':resid,\
                                    'x':short_line['coord'][0],'y':short_line['coord'][1],'z':short_line['coord'][2]}
                at_id+=1
        merge_coords=at_mod.check_atom_overlap(coord)

        for at_id, coord in enumerate(merge_coords):
            pdb_output.write(g_var.pdbline%((at_id+1,final_atom[at_id]['atom'],final_atom[at_id]['res_type'],final_atom[at_id]['chain'],final_atom[at_id]['residue'],
                 coord[0],coord[1],coord[2],1,0))+'\n') 
################################################################## Merge chains ####################

def read_in_protein_pdbs(no_chains, file, end):
#### reads in each chain into merge list
    merge, merged_coords = [],[]
    for chain in range(0,no_chains):
        merge_temp = []
        if os.path.exists(file+'_'+str(chain)+end):
            with open(file+'_'+str(chain)+end, 'r') as pdb_input:
                for line in pdb_input.readlines():
                    if line.startswith('ATOM'):
                        line_sep=gen.pdbatom(line)
                        merge_temp.append(line_sep)
        else:
            sys.exit('cannot find minimised protein chain: '+str(chain)) 
        merge, merge_coords = at_mod.fix_chirality(merge,merge_temp,merged_coords)    
    merged_coords = at_mod.check_atom_overlap(merge_coords)
    merged=[]
    for line_val, line in enumerate(merge):
        merged.append(g_var.pdbline%((int(line['atom_number']), line['atom_name'], line['residue_name'],' ',line['residue_id'],\
            merged_coords[line_val][0],merged_coords[line_val][1],merged_coords[line_val][2],1,0))+'\n')
    return merged



def write_merged_pdb(merge, protein, box_vec):
#### creates merged pdb and writes chains to it
    if not os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN'+protein+'_merged.pdb'):
        pdb_output=gen.create_pdb(g_var.working_dir+'PROTEIN/PROTEIN'+protein+'_merged.pdb', box_vec)
        for line in merge:
            pdb_output.write(line)
        pdb_output.close()


########################################################### RMSD

def RMSD_measure(structure_atoms, system, backbone_coords):
    RMSD_dict = {}
    for chain in range(system['PROTEIN']):
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
        if len(at_centers) != len(backbone_coords[chain]):
            sys.exit('In chain '+str(chain)+' the atommistic input does not match the CG. \n\
    number of CG residues '+str(len(backbone_coords[chain]))+'\nnumber of AT residues '+str(len(at_centers)))
        #### finds distance between backbone COM and cg backbone beads
        dist=np.sqrt((np.array(at_centers) - np.array(backbone_coords[chain])[:,:3])**2)
        RMSD_val = np.sqrt(np.mean(dist**2)) #### RMSD calculation
        RMSD_dict[chain]=np.round(RMSD_val, 3)  #### stores RMSD in dictionary
    return RMSD_dict


