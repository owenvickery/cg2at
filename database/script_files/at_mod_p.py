#!/usr/bin/env python3

import os, sys
import numpy as np
import copy
from string import ascii_uppercase
import difflib
import gen, g_var, f_loc, at_mod
import math




def build_protein_atomistic_system(cg_residues, box_vec):
#### initisation of counters
    chain_information=[]
    chain_count=0
    system={}
    at_residues={}
    backbone_coords={}
    backbone_coords[chain_count]=[]
    terminal={}
    terminal[chain_count]=[]
    temporary_coordinates_atomistic={}
    temporary_coordinates_atomistic[chain_count]={}
    sequence={}
    sequence[chain_count]=[]
    res=0
    print('Converting Protein')
    gen.mkdir_directory(g_var.working_dir+'PROTEIN')  ### make and change to protein directory
#### create new pdb file for chain 0 
    pdb_output = gen.create_pdb(g_var.working_dir+'PROTEIN/PROTEIN_novo_'+str(chain_count)+'.pdb', box_vec)
#### for each residue in protein
    initial=True

    for cg_residue_id, residue_number in enumerate(cg_residues):
    #### temporary index/dictionaries       
        at_residues[residue_number]={}
    #### fetch fragments in residue and connectivity
        at_residues[residue_number], connect=at_mod.get_atomistic_fragments(cg_residues[residue_number][next(iter(cg_residues[residue_number]))]['residue_name'],cg_residues[residue_number], residue_number)
    #### if residue contains BB bead a index of the BB connectivity is collected
        if 'BB' in at_residues[residue_number]:
            BB_connect=[] ### backbone connectivity
            for atom_num, atom in enumerate(at_residues[residue_number]['BB']):
                if at_residues[residue_number]['BB'][atom]['atom'] in f_loc.backbone[cg_residues[residue_number]['BB']['residue_name']]['b_connect']:
                    BB_connect.append(atom)
        
    #### for each bead in residue
        for frag_val,cg_fragments in enumerate(cg_residues[residue_number]):
        #### gets connectivity between fragents 
            at_connections, cg_connections, center=at_mod.connectivity(frag_val, cg_fragments, connect, at_residues[residue_number], cg_residues, residue_number)
        #### if BB bead adds the N and C terminal atoms to connectivity
            if cg_fragments=='BB':
                at_connections, cg_connections, res = BB_connectivity(at_connections,cg_connections, cg_residues, at_residues, residue_number, BB_connect, res, center)
            #### measures the distance between BB beads. 
                if not initial:
                    residue_number, terminal, chain_count, temporary_coordinates_atomistic, backbone_coords, chain_information, res, sequence = check_chain(
                        cg_residues, residue_number, terminal, chain_count, temporary_coordinates_atomistic, backbone_coords, chain_information, res, sequence)
                #### if not prev residue the xyz coord of the cg_bead are added to the backbone_coords dictionary
                else:
                    xyz_cur=[cg_residues[residue_number]['BB']['coord'][0],cg_residues[residue_number]['BB']['coord'][1],cg_residues[residue_number]['BB']['coord'][2]]
                    backbone_coords[chain_count].append(xyz_cur+[1])
                    initial=False
                    terminal[chain_count].append(f_loc.backbone[cg_residues[residue_number]['BB']['residue_name']]['ter'])
                    sequence = at_mod.add_to_sequence(sequence, cg_residues[residue_number]['BB']['residue_name'], chain_count)
            if cg_residues[residue_number][cg_fragments]['residue_name'] == 'CYS' and cg_fragments != 'BB':
                at_connections, cg_connections, disulphide, disul_at_info, disul_cg_info= find_closest_cysteine(at_connections, cg_connections, cg_residues, at_residues, residue_number, BB_connect, res, center)
        #### finds optimum rotation of fragment
            if len(at_connections) == len(cg_connections):
                if cg_residues[residue_number][cg_fragments]['residue_name'] == 'PRO' and cg_fragments != 'BB':
                    xyz_rot_apply=at_mod.rotate(np.array(at_connections), np.array(cg_connections), True)
                else:
                    xyz_rot_apply=at_mod.rotate(np.array(at_connections), np.array(cg_connections), False)
            else:
                print('atom connections: '+str(len(at_connections))+' does not equal CG connections: '+str(len(cg_connections)))
                sys.exit('residue number: '+str(residue_number)+', residue type: '+str(cg_residues[residue_number][cg_fragments]['residue_name'])+', bead: '+cg_fragments)


        #### applies rotation to each atom
            for atom in at_residues[residue_number][cg_fragments]:
                at_residues[residue_number][cg_fragments][atom]['coord'] = at_mod.rotate_atom(at_residues[residue_number][cg_fragments][atom]['coord'], center, xyz_rot_apply)     
    #### if disulphide bond found move the S atoms to within 2 A of each other
        if 'disulphide' in locals():
            if disulphide:
                at_residues[residue_number] = shift_sulphur(residue_number, disul_at_info, disul_cg_info, at_residues, cg_residues) 
                disulphide = False
        temporary_coordinates_atomistic[chain_count][residue_number]=at_residues[residue_number]

    final_coordinates_atomistic = finalise_novo_atomistic(temporary_coordinates_atomistic, cg_residues, box_vec)

    terminal[chain_count].append(f_loc.backbone[cg_residues[residue_number]['BB']['residue_name']]['ter'])

    system['terminal_residue']=terminal
    system['PROTEIN']=chain_count+1
    if g_var.v >=1:
        print('\n{:-<75}'.format('>  Verbose level 1 start'))
        print('\nchain number\tDelta A\t\tno in pdb\tlength of chain')
        print('------------\t-------\t\t---------\t---------------')
        for chain in range(chain_count):
            print('\t',chain,'\t',np.round(chain_information[chain][0], 1),'\t\t',chain_information[chain][1]-chain_information[chain][2]+1,'-',chain_information[chain][1],'\t\t',chain_information[chain][2])
        if chain_count >1:
            print('\t', chain_count,'\tN/A\t\t',chain_information[chain][1]+2,'-',residue_number+1,'\t\t',residue_number-res)
        else:
            print('\t', chain_count,'\tN/A\t\t',res+2,'-',residue_number+1,'\t\t',residue_number-res)

        print('\n{:-<75}'.format('>  Verbose level 1 end\n'))
    return system, backbone_coords, final_coordinates_atomistic, sequence

def check_chain(cg_residues, residue_number, terminal, chain_count, temporary_coordinates_atomistic, backbone_coords, chain_information, res, sequence):
    xyz_prev=[cg_residues[residue_number-1]['BB']['coord'][0],cg_residues[residue_number-1]['BB']['coord'][1],cg_residues[residue_number-1]['BB']['coord'][2]]              
    xyz_cur=[cg_residues[residue_number]['BB']['coord'][0],cg_residues[residue_number]['BB']['coord'][1],cg_residues[residue_number]['BB']['coord'][2]]
    dist=np.sqrt(((xyz_prev[0]-xyz_cur[0])**2)+((xyz_prev[1]-xyz_cur[1])**2)+((xyz_prev[2]-xyz_cur[2])**2))
#### if distance between BB beads is more than 5 A then it is considered a new chain.
    if dist > 5:
        terminal[chain_count].append(f_loc.backbone[cg_residues[residue_number-1]['BB']['residue_name']]['ter'])
        chain_count+=1  ### adds to to the protein count
        temporary_coordinates_atomistic[chain_count]={}
        backbone_coords[chain_count]=[]   #### creates another dictionary key for bb fragments 
        backbone_coords[chain_count].append(xyz_cur+[1])  #### adds xyz coord and mass of 1 to list
        chain_information.append([dist, residue_number, residue_number-res])  ### info for verbose flag
        res=residue_number-1 #### updates residue
        terminal[chain_count]=[]
        terminal[chain_count].append(f_loc.backbone[cg_residues[residue_number]['BB']['residue_name']]['ter'])
        sequence[chain_count]=[]
        sequence = at_mod.add_to_sequence(sequence, cg_residues[residue_number]['BB']['residue_name'], chain_count)

    else:
    #### the xyz coord of the BB bead are added to the backbone_coords dictionary
        backbone_coords[chain_count].append(xyz_cur+[1])
        sequence = at_mod.add_to_sequence(sequence, cg_residues[residue_number]['BB']['residue_name'], chain_count)
    return residue_number, terminal, chain_count, temporary_coordinates_atomistic, backbone_coords, chain_information, res, sequence

################# Fixes disulphide bond, martini cysteine bone is too far apart to be picked up by pdb2gmx. 
#### 

def find_closest_cysteine(at_connections, cg_connections, cg_residues, at_residues, residue_number, BB_connect, res, center):
    disulphide, atom_number=False, 0
    for res_id in cg_residues: 
        try:
            #### checks distance between cysteines if closer than 7A then adds another connection between the S and the sidechain of the other cysteine CG bead
            if cg_residues[res_id]['SC1']['residue_name'] == 'CYS' and res_id not in [residue_number-1, residue_number, residue_number+1]:
                xyz_cur=[cg_residues[residue_number]['SC1']['coord'][0],cg_residues[residue_number]['SC1']['coord'][1],cg_residues[residue_number]['SC1']['coord'][2]] ### cysteine of interest
                xyz_check=[cg_residues[res_id]['SC1']['coord'][0],cg_residues[res_id]['SC1']['coord'][1],cg_residues[res_id]['SC1']['coord'][2]] ### cysteine to check distance
                dist=np.sqrt(((xyz_check[0]-xyz_cur[0])**2)+((xyz_check[1]-xyz_cur[1])**2)+((xyz_check[2]-xyz_cur[2])**2)) ### distance
                if dist < 6:
                    for atom_number in at_residues[residue_number]['SC1']:
                        if at_residues[residue_number]['SC1'][atom_number]['atom']==f_loc.backbone[cg_residues[residue_number]['SC1']['residue_name']]['disulphide']: ### if sulphur
                            at_connections.append(at_residues[residue_number]['SC1'][atom_number]['coord']-center) ### add at centered coordinates
                    cg_connections.append(cg_residues[res_id]['SC1']['coord']-center) ### add cg centered coordinates
                    disulphide=True
                    break
        except:
            pass
    return at_connections, cg_connections, disulphide, atom_number, res_id

def shift_sulphur(residue_number, disul_at_info, disul_cg_info, at_residues, cg_residues ):
#### to shift sidechains roughly equally the 1st sidechain is shifted to within 3.2A and the second to 2A (pdb2gmx cutoff)
    if disul_cg_info > residue_number:
        xyz_cur=np.array([cg_residues[disul_cg_info]['SC1']['coord'][0],cg_residues[disul_cg_info]['SC1']['coord'][1],cg_residues[disul_cg_info]['SC1']['coord'][2]])
        cutoff=3.2
    else:
        xyz_cur=np.array([at_residues[disul_cg_info]['SC1'][disul_at_info]['coord'][0],at_residues[disul_cg_info]['SC1'][disul_at_info]['coord'][1],at_residues[disul_cg_info]['SC1'][disul_at_info]['coord'][2]])
        cutoff=2

    xyz_check=np.array([at_residues[residue_number]['SC1'][disul_at_info]['coord'][0],at_residues[residue_number]['SC1'][disul_at_info]['coord'][1],at_residues[residue_number]['SC1'][disul_at_info]['coord'][2]])
    dist=np.sqrt(((xyz_check[0]-xyz_cur[0])**2)+((xyz_check[1]-xyz_cur[1])**2)+((xyz_check[2]-xyz_cur[2])**2))
#### moves sidechains closer together in increments of 5% of the length of the vector 
    offset=0
    if dist >= cutoff:
        vector=xyz_cur-xyz_check
        x=0.05
        while True:
            xyz_check_new = xyz_check + ( vector * x )
            dist=np.sqrt(((xyz_check_new[0]-xyz_cur[0])**2)+((xyz_check_new[1]-xyz_cur[1])**2)+((xyz_check_new[2]-xyz_cur[2])**2))
            if dist >= cutoff:
                x+=0.05
            else:
                offset = vector * x 
                break
#### applies final shift to the rest of the atoms in the sidechain
    for atom in at_residues[residue_number]['SC1']:
        at_residues[residue_number]['SC1'][atom]['coord']=at_residues[residue_number]['SC1'][atom]['coord']+offset
    return at_residues[residue_number]

def BB_connectivity(at_connections,cg_connections, cg_residues, at_residues, residue_number, BB_connect, res, center):
#### connect to preceeding backbone bead in chain
    try:
        cg_connections.append(cg_residues[residue_number-1]['BB']['coord']-center)
        at_connections.append(at_residues[residue_number]['BB'][BB_connect[0]]['coord']-center)
    except:
        res=residue_number
        pass
#### connect to next backbone bead in chain
    try:
        cg_connections.append(cg_residues[residue_number+1]['BB']['coord']-center)
        at_connections.append(at_residues[residue_number]['BB'][BB_connect[1]]['coord']-center)
    except:
        pass
    return at_connections,cg_connections, res

def finalise_novo_atomistic(temporary_coordinates_atomistic, cg_residues, box_vec):
    final_at_residues={}
    final_at_resids = {}
    for chain in temporary_coordinates_atomistic:
        pdb_output = gen.create_pdb(g_var.working_dir+'PROTEIN/PROTEIN_novo_'+str(chain)+'.pdb', box_vec)
        final_at_residues[chain], final_at_resids[chain] = fix_dihedrals(temporary_coordinates_atomistic[chain], cg_residues)
    #### check if any atom overlap if so give the atom a kick
        coords=[]
        for atom in final_at_residues[chain]:
            coords.append(final_at_residues[chain][atom]['coord'])
        coords = at_mod.check_atom_overlap(coords)
        for line_index, atom in enumerate(final_at_residues[chain]):
            final_at_residues[chain][atom]['coord']=coords[line_index]
            pdb_output.write(g_var.pdbline%((line_index,final_at_residues[chain][atom]['atom'],final_at_residues[chain][atom]['res_type'],ascii_uppercase[chain],\
                final_at_residues[chain][atom]['resid'],final_at_residues[chain][atom]['coord'][0],final_at_residues[chain][atom]['coord'][1],final_at_residues[chain][atom]['coord'][2],1,0))+'\n')
        pdb_output.close()
    return final_at_resids


########################################### fix carbonyl section 


def fix_carbonyl(residue_id, cg, at):
    ca=[]
    for index in range(3):
        for atom in at[residue_id+index]['BB']:
            if at[residue_id+index]['BB'][atom]['atom'] in f_loc.backbone[cg[residue_id+index]['BB']['residue_name']]['restraint']: 
                ca.append(at[residue_id+index]['BB'][atom]['coord'])
            if index == 0 :
                if at[residue_id+index]['BB'][atom]['atom'] == f_loc.backbone[cg[residue_id+index]['BB']['residue_name']]['b_connect'][1]: 
                    C = atom
                if at[residue_id+index]['BB'][atom]['atom'] in f_loc.backbone[cg[residue_id+index]['BB']['residue_name']]['dihedral']:   
                    O = atom                 

    initial_vector, cross_vector = at_mod.find_cross_vector( ca, at[residue_id]['BB'][C]['coord'], at[residue_id]['BB'][O]['coord'])
    rotation = at_mod.align_to_vector(initial_vector, cross_vector)
    center = ca[0]+(ca[1]-ca[0])/3
    at[residue_id]['BB'][C]['coord'] = (at[residue_id]['BB'][C]['coord']-center).dot(rotation)+center
    at[residue_id]['BB'][O]['coord'] = (at[residue_id]['BB'][O]['coord']-center).dot(rotation)+center
    return at



def fix_dihedrals(atomistic_residues, cg_residues):
    for res_index, residue_id in enumerate(atomistic_residues):
        if res_index < len(atomistic_residues)-2:
            atomistic_residues = fix_carbonyl(residue_id, cg_residues, atomistic_residues)
    final = {}
    atom_count=0
    residue_number=0
    final_resid={}
    for residue_id in atomistic_residues:
        final_resid[residue_id]={}
        residue_length=0
        residue_temp={}
        for cg_fragments in atomistic_residues[residue_id]:
            residue_temp.update(atomistic_residues[residue_id][cg_fragments])
        for atom in range(1, len(residue_temp)+1):
            final[atom_count]=residue_temp[atom]
            final[atom_count].update({'resid':residue_number})
            final_resid[residue_id][atom]=residue_temp[atom]
            atom_count+=1
        residue_number+=1
    return final, final_resid

##################################################################  User supplied protein ##############

def read_in_atomistic(protein, cg_chain_count, sequence, check_alignment):
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
                run=True
            #### if line is correct
            if run:
                if line_sep['residue_name'] in f_loc.p_residues or line_sep['residue_name'] in f_loc.mod_residues:
                    
                    if 'H' not in line_sep['atom_name'][0] or line_sep['residue_name'] in f_loc.mod_residues:  
                    #### sorts out wrong atoms in terminal residues
                        if line_sep['atom_name'] in ['OT', 'O1', 'O2']:
                            line_sep['atom_name']='O'
                    #### makes C_terminal connecting atom variable  
                        if line_sep['atom_name'] == f_loc.backbone[line_sep['residue_name']]['b_connect'][1]:
                            C_ter=[line_sep['x'],line_sep['y'],line_sep['z']]
                            C_resid=line_sep['residue_id']
                            C=True
                        try:
                        #### tries to make a N_terminal connecting atom variable
                            if line_sep['atom_name'] == f_loc.backbone[line_sep['residue_name']]['b_connect'][0]:
                                N_resid=line_sep['residue_id']
                                N_ter=[line_sep['x'],line_sep['y'],line_sep['z']]
                                N=True
                        #### measures distance between N and C atoms. if the bond is over 3 A it counts as a new protein
                            dist=np.sqrt(((N_ter[0]-C_ter[0])**2)+((N_ter[1]-C_ter[1])**2)+((N_ter[2]-C_ter[2])**2))
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
                else:
                    if check_alignment:
                        sys.exit('The residue '+line_sep['residue_name']+' does not exist in the fragment database')
    if check_alignment:
        seq_user = check_sequence(atomistic_protein_input, chain_count+1)
        atomistic_protein_input = align_chains(atomistic_protein_input, seq_user, sequence)
#### check if number of monomers is the same
    elif chain_count+1 != cg_chain_count:
        sys.exit('number of chains in atomistic protein input ('+str(chain_count+1)+') does not match CG representation ('+str(cg_chain_count)+')')
    return atomistic_protein_input

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
            print(sequence[index], '\n')
        print('\nuser supplied structure:\n')
        for index in seq_user:
            print(seq_user[index], '\n')        



    at={}
    sequence_temp = sequence.copy()
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
    return at

def mask_sequence(sequence, st, end):
    for index, residue in enumerate(sequence):
        if st <= index < end:
            sequence[index]='-'
    return sequence

def center_atomistic(atomistic_protein_input, backbone_coords): 
    cg_com=[]
#### for each protein chain center on cg representation 
    for chain in range(len(atomistic_protein_input)):
        cg_com.append([])
        for part_val, part in enumerate(atomistic_protein_input[chain]):
            sls, sle= int(part.split(':')[0]),int(part.split(':')[1])
            protein_mass=[]
            for residue in atomistic_protein_input[chain][part]:
            #### creates a list of all coordinates and masses [[coord, mass],[coord, mass]]
                for atom in atomistic_protein_input[chain][part][residue]:
                    short_line=atomistic_protein_input[chain][part][residue][atom]
                    protein_mass.append([short_line['coord'][0],short_line['coord'][1],short_line['coord'][2],short_line['frag_mass']])
        #### returns the COM of the atomistic protein
            atomistic_protein_mass=np.average(np.array(protein_mass)[:,:3], axis=0, weights=np.array(protein_mass)[:,3])#### add center of mass of CG_proteins
        #### for each chain the COM of the CG representation is stored (only cg is needed)
            cg_com[chain].append(np.average(np.array(backbone_coords[chain])[sls:sle,:3], axis=0, weights=np.array(backbone_coords[chain])[sls:sle,3]))
        #### each atoms coord is updated so the monomer COM is the same as the CG
            for residue in atomistic_protein_input[chain][part]:
                for atom in atomistic_protein_input[chain][part][residue]:
                    atomistic_protein_input[chain][part][residue][atom]['coord']=atomistic_protein_input[chain][part][residue][atom]['coord']-(atomistic_protein_mass-cg_com[chain][part_val])
    return atomistic_protein_input, cg_com

def rotate_protein_monomers(atomistic_protein_centered, final_coordinates_atomistic, backbone_coords, cg_com,  box_vec):
#### run through each chain in proteins
    for chain in range(len(atomistic_protein_centered)):
    #### creates atomistic pdb
        xyz_rot_apply=[]
        for part_val, part in enumerate(atomistic_protein_centered[chain]):
            sls, sle= int(part.split(':')[0]),int(part.split(':')[1])        
            at_centers=[]
        #### runs through every residue and atom  
            for residue in atomistic_protein_centered[chain][part]:
            #### gets center of mass of each residue (note only backbone heavy atoms have a mass)
                at_centers_iter=[]
                for atom in atomistic_protein_centered[chain][part][residue]:
                    at_centers_iter.append(np.append(atomistic_protein_centered[chain][part][residue][atom]['coord'],atomistic_protein_centered[chain][part][residue][atom]['frag_mass']))
                try:
                    at_centers.append(np.average(np.array(at_centers_iter)[:,:3], axis=0, weights=np.array(at_centers_iter)[:,3]))
                except:
                    for atom in atomistic_protein_centered[chain][part][residue]:
                        print(atomistic_protein_centered[chain][part][residue][atom])
                    sys.exit()
        #### finds optimal rotation of each monomer  
            if len(at_centers) == len(np.array(backbone_coords[chain])[sls:sle,:3]):
                xyz_rot_apply.append(at_mod.rotate(np.array(at_centers)-cg_com[chain][part_val], 
                                     np.array(backbone_coords[chain])[sls:sle,:3]-cg_com[chain][part_val], False))
            else:
                sys.exit('In chain '+str(chain)+' the atomistic input does not match the CG. \n\
    number of CG residues '+str(len(backbone_coords[chain]))+'\nnumber of AT residues '+str(len(at_centers)))

            if g_var.v >= 1:
                print('\nThe proteins chains are rotated around the COM of all the backbone heavy atoms.')
                print('The COM of chain', chain,'is :', np.round(cg_com[chain][part_val][0], 2),',', np.round(cg_com[chain][part_val][1], 2),',', 
                      np.round(cg_com[chain][part_val][2], 2))
                print('rotating chain ', chain, 'by :',np.round(np.degrees(xyz_rot_apply[part_val][0]),2),',',np.round(np.degrees(xyz_rot_apply[part_val][1]),2),
                      ',',np.round(np.degrees(xyz_rot_apply[part_val][2]),2))
                print()
        hybridise_protein_inputs(final_coordinates_atomistic[chain], atomistic_protein_centered[chain], cg_com[chain], xyz_rot_apply, chain, box_vec)

def hybridise_protein_inputs(final_coordinates_atomistic, atomistic_protein_centered, cg_com, xyz_rot_apply, chain, box_vec):
    pdb_output = gen.create_pdb(g_var.working_dir+'PROTEIN/PROTEIN_at_rep_user_supplied_'+str(chain)+'.pdb', box_vec)
    final_atom={}
    coord=[]
    at_id=0
    for residue in final_coordinates_atomistic:
        exists=False
        for initial_index in final_coordinates_atomistic[residue]:
            if final_coordinates_atomistic[residue][initial_index]['res_type'] in f_loc.mod_residues:
                for atom in final_coordinates_atomistic[residue]:
                    short_line=final_coordinates_atomistic[residue][atom]
                    final_atom[at_id]={'atom':short_line['atom'], 'res_type':short_line['res_type'], 'chain':ascii_uppercase[chain], 'residue':residue,\
                                 'x':short_line['coord'][0],'y':short_line['coord'][1],'z':short_line['coord'][2]}
                    at_id+=1
                    coord.append(short_line['coord'])
            elif final_coordinates_atomistic[residue][initial_index]['res_type'] not in f_loc.mod_residues:
                for part_val, part in enumerate(atomistic_protein_centered):
                    if residue in atomistic_protein_centered[part]:
                        exists=True
                        for atom in atomistic_protein_centered[part][residue]:                               
                            atomistic_protein_centered[part][residue][atom]['coord'] = at_mod.rotate_atom(atomistic_protein_centered[part][residue][atom]['coord'], cg_com[part_val], xyz_rot_apply[part_val])
                            short_line = atomistic_protein_centered[part][residue][atom]
                            final_atom[at_id]={'atom':short_line['atom'], 'res_type':short_line['res_type'], 'chain':ascii_uppercase[chain], 'residue':residue,\
                                        'x':short_line['coord'][0],'y':short_line['coord'][1],'z':short_line['coord'][2]}
                            at_id+=1
                            coord.append(short_line['coord'])
            if not exists:
                for atom in final_coordinates_atomistic[residue]:
                    short_line=final_coordinates_atomistic[residue][atom]  
                    final_atom[at_id]={'atom':short_line['atom'], 'res_type':short_line['res_type'], 'chain':ascii_uppercase[chain], 'residue':residue,\
                                'x':short_line['coord'][0],'y':short_line['coord'][1],'z':short_line['coord'][2]}
                    at_id+=1
                    coord.append(short_line['coord'])
            break
    merge_coords = at_mod.check_atom_overlap(coord)
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


