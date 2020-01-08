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

            if resname == 'CYS' and 'BB' not in group_fit:
                at_connect, cg_connect, disulphide, disul_at_info, disul_cg_info= find_closest_cysteine(at_connect, cg_connect, cg_residues, group_fit, residue_number)

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
        #### if disulphide bond found move the S atoms to within 2 A of each other
            if 'disulphide' in locals():
                if disulphide:
                    coordinates_atomistic[chain_count][residue_number] = shift_sulphur(residue_number, disul_at_info, disul_cg_info, coordinates_atomistic[chain_count], cg_residues) 
                    disulphide = False

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
        print('\n{:-<75}'.format('>  Verbose level 1 start'))
        print('\n{0:^15}{1:^12}'.format('chain number', 'length of chain')) #   \nchain number\tDelta A\t\tno in pdb\tlength of chain')
        print('\n{0:^15}{1:^12}'.format('------------', '---------------'))
        for chain in sequence:
            print('{0:^15}{1:^12}'.format(chain, len(sequence[chain])))
        print()
    final_coordinates_atomistic = finalise_novo_atomistic(coordinates_atomistic, cg_residues, box_vec)
    system['terminal_residue']=terminal
    system['PROTEIN']=chain_count
    return system, backbone_coords, final_coordinates_atomistic, sequence    


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

def find_closest_cysteine(at_connections, cg_connections, cg_residues, group, residue_number):
    disulphide, atom_number=False, 0
    for res_id in cg_residues: 
        try:
            #### checks distance between cysteines if closer than 7A then adds another connection between the S and the sidechain of the other cysteine CG bead
            if cg_residues[res_id]['SC1']['residue_name'] == 'CYS' and res_id not in [residue_number-1, residue_number, residue_number+1]:
                xyz_cur=[cg_residues[residue_number]['SC1']['coord'][0],cg_residues[residue_number]['SC1']['coord'][1],cg_residues[residue_number]['SC1']['coord'][2]] ### cysteine of interest
                xyz_check=[cg_residues[res_id]['SC1']['coord'][0],cg_residues[res_id]['SC1']['coord'][1],cg_residues[res_id]['SC1']['coord'][2]] ### cysteine to check distance
                dist = gen.calculate_distance(xyz_check, xyz_cur)
                if dist < g_var.cys:
                    for atom_number in group['SC1']:
                        if group['SC1'][atom_number]['atom']==f_loc.backbone[cg_residues[residue_number]['SC1']['residue_name']]['sul']: ### if sulphur
                            at_connections.append(group['SC1'][atom_number]['coord']) ### add at centered coordinates
                    cg_connections.append(cg_residues[res_id]['SC1']['coord']) ### add cg centered coordinates
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
        xyz_cur=np.array([at_residues[disul_cg_info][disul_at_info]['coord'][0],at_residues[disul_cg_info][disul_at_info]['coord'][1],at_residues[disul_cg_info][disul_at_info]['coord'][2]])
        cutoff=2
    xyz_check=np.array([at_residues[residue_number][disul_at_info]['coord'][0],at_residues[residue_number][disul_at_info]['coord'][1],at_residues[residue_number][disul_at_info]['coord'][2]])
    dist=gen.calculate_distance(xyz_check,xyz_cur) 
#### moves sidechains closer together in increments of 5% of the length of the vector 
    offset=0
    if dist >= cutoff:
        vector=xyz_cur-xyz_check
        x=0.05
        while True:
            xyz_check_new = xyz_check + ( vector * x )
            dist=gen.calculate_distance(xyz_check_new,xyz_cur) 
            if dist >= cutoff:
                x+=0.05
            else:
                offset = vector * x 
                break
#### applies final shift to the rest of the atoms in the sidechain
    for atom in at_residues[residue_number]:
        at_residues[residue_number][atom]['coord']=at_residues[residue_number][atom]['coord']+offset
    return at_residues[residue_number]

def finalise_novo_atomistic(atomistic, cg_residues, box_vec):
    final_at_residues={}
    final_at = {}
    for chain in atomistic: 
        at_number=0    
        final_at_residues[chain]={}  
        final_at[chain]={}
        coords=[]
        skip = os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN_novo_'+str(chain)+'.pdb')
        if not skip:
            pdb_output = gen.create_pdb(g_var.working_dir+'PROTEIN/PROTEIN_novo_'+str(chain)+'.pdb', box_vec)
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
                pdb_output.write(g_var.pdbline%((atom,final_at[chain][atom]['atom'],final_at[chain][atom]['res_type'],ascii_uppercase[chain],\
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
            print('chain:', index,sequence[index], '\n')
        print('\nuser supplied structure:\n')
        for index in seq_user:
            print('chain:', index, seq_user[index], '\n')        
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
            atomistic_protein_mass = at_mod.COM(protein_mass, 'protein at: '+str(chain)+' '+part)#### add center of mass of CG_proteins
            cg_com[chain].append(at_mod.COM(backbone_coords[chain][sls:sle], 'protein cg: '+str(chain)+' '+part))
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
                            if atomistic_protein_centered[part][residue][atom]['res_type'] != final_coordinates_atomistic[residue][initial_index]['res_type']:
                                print('de_novo' , final_coordinates_atomistic[residue][initial_index]['res_type'],'at_user', atomistic_protein_centered[part][residue][atom]['res_type'])
                                sys.exit('de novo and at user supplied don\'t match')
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


