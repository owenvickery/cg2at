#!/usr/bin/env python3

import os, sys
import numpy as np
import difflib
from scipy.spatial import cKDTree
import copy
import gen, g_var, at_mod, read_in

def build_multi_residue_atomistic_system(cg_residues, sys_type):   
#### initisation of counters
    chain_count=0
    coord_atomistic={}
    g_var.seq_cg={sys_type:{}}
    g_var.ter_res={sys_type:{}}
    gen.mkdir_directory(g_var.working_dir+sys_type)  ### make and change to protein directory
#### for each residue in protein
    residue_type={}
    residue_type_mass={}
    new_chain = True
    for cg_residue_id, residue_number in enumerate(cg_residues[sys_type]):

        if np.round((cg_residue_id/len(cg_residues[sys_type]))*100,2).is_integer():
            print('Converting de_novo '+sys_type+': ',np.round((cg_residue_id/len(cg_residues[sys_type]))*100,2),'%', end='\r')
        resname = cg_residues[sys_type][residue_number][next(iter(cg_residues[sys_type][residue_number]))]['residue_name']
        if new_chain: 
            if chain_count not in coord_atomistic:
                if sys_type == 'PROTEIN':
                    g_var.backbone_coords[chain_count]=[]
                coord_atomistic[chain_count]={}
                g_var.seq_cg[sys_type][chain_count]=[]
                g_var.ter_res[sys_type][chain_count]=[resname, False]
            new_chain = False
        coord_atomistic[chain_count][residue_number]={}
        frag_location=gen.fragment_location(resname) ### get fragment location from database
        residue_type[resname], residue_type_mass[resname] = at_mod.get_atomistic(frag_location)
        g_var.seq_cg[sys_type] = add_to_sequence(g_var.seq_cg[sys_type], resname, chain_count)
        new_chain = False
        for group in residue_type[resname]:
            for key in list(residue_type[resname][group].keys()):
                if key not in cg_residues[sys_type][residue_number]:
                    del residue_type[resname][group][key]
            if len(residue_type[resname][group]) > 0:
                center, at_frag_centers, cg_frag_centers, group_fit = at_mod.rigid_fit(residue_type[resname][group], residue_type_mass[resname], residue_number, cg_residues[sys_type][residue_number])
                at_connect, cg_connect = at_mod.connectivity(cg_residues[sys_type][residue_number], at_frag_centers, cg_frag_centers, group_fit, group)
                for group_bead in group_fit:
                    if group_bead in g_var.res_top[resname]['CONNECT']:
                        at_connect, cg_connect, new_chain = at_mod.BB_connectivity(at_connect,cg_connect, cg_residues[sys_type], group_fit[group_bead], residue_number, group_bead, resname)
                        if sys_type == 'PROTEIN':
                            g_var.backbone_coords[chain_count].append(np.append(cg_residues[sys_type][residue_number][group_bead]['coord'], 1)) 
                xyz_rot_apply = at_mod.get_rotation(cg_connect, at_connect, center, resname, group, residue_number)
                coord_atomistic[chain_count] = at_mod.apply_rotations(coord_atomistic[chain_count],residue_number, group_fit, center, xyz_rot_apply)
        if new_chain:
            g_var.ter_res[sys_type][chain_count][1] = resname
            chain_count+=1

    print('Completed initial conversion of '+sys_type+'\n')        
    g_var.system[sys_type]=chain_count
    if sys_type == 'PROTEIN':
        for chain in range(chain_count):
            g_var.skip_disul[chain]=False
    return coord_atomistic

################# Fixes disulphide bond, martini cysteine bone is too far apart to be picked up by pdb2gmx. 
#### 

def ask_if_disulphide(chain, res_1, res_2):
    while True:
        try:
            answer = str(input('\nAre these residues connected by a disulphide bond in chain '+str(chain)+' (Y/N): '+str(res_1)+'--'+str(res_2)+': '))
            if answer.lower() in ['yes','y']:
                return True
            elif answer.lower() in ['no','n']:
                return False
            else:
                print("Oops!  That was a invalid choice")
        except KeyboardInterrupt:
            sys.exit('\nInterrupted')
        except BaseException:
            print("Oops!  That was a invalid choice")

def find_disulphide_bonds_user_sup():
    for chain in g_var.atomistic_protein_input_aligned:
        g_var.user_cys_bond[chain]=[]
        cysteines = []
        cys_resid = []
        for part in g_var.atomistic_protein_input_aligned[chain]:
            for resid in g_var.atomistic_protein_input_aligned[chain][part]:
                atom = next(iter(g_var.atomistic_protein_input_aligned[chain][part][resid]))
                if g_var.atomistic_protein_input_aligned[chain][part][resid][atom]['res_type'] == 'CYS':
                    for atom in g_var.atomistic_protein_input_aligned[chain][part][resid]:
                        if 'S' in g_var.atomistic_protein_input_aligned[chain][part][resid][atom]['atom'].upper() :
                            cys_resid.append(resid)
                            cysteines.append(g_var.atomistic_protein_input_aligned[chain][part][resid][atom]['coord'])
        if len(cysteines)>=2:
            tree = cKDTree(cysteines)
            done_query=[]
            for cys_index, cys in enumerate(cysteines):
                query = tree.query_ball_point(cys, r=2.3)
                if len(query) == 2 and query not in done_query:
                    g_var.user_cys_bond[chain].append([cys_resid[query[0]],cys_resid[query[1]]])
                    done_query.append(query)

def find_disulphide_bonds_de_novo():
    for chain in g_var.coord_atomistic:
        if not g_var.skip_disul[chain]:
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
                    query = tree.query_ball_point(cys, r=g_var.args.cys)
                    if len(query) == 2 and query not in done_query:
                        if cys_resid[query[1]] - cys_resid[query[0]] > 4:
                            disul = True if g_var.args.silent else ask_if_disulphide(chain, cys_resid[query[0]],cys_resid[query[1]])
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

#### disulphide checker end

### writes out final de_novo atomistic coordinates

def finalise_novo_atomistic(coord_atomistic, sys_type):
    final_at_residues={}
    final_at = {}
    for chain in coord_atomistic: 
        at_number=0    
        final_at_residues[chain]={}  
        final_at[chain]={}
        coords=[]
        skip = os.path.exists(g_var.working_dir+sys_type+'/'+sys_type+'_de_novo_'+str(chain)+'.pdb')
        if not skip:
            pdb_output = gen.create_pdb(g_var.working_dir+sys_type+'/'+sys_type+'_de_novo_'+str(chain)+'.pdb')
        for res_index, residue_id in enumerate(coord_atomistic[chain]):
            if coord_atomistic[chain][residue_id][next(iter(coord_atomistic[chain][residue_id]))]['res_type'] in g_var.mod_residues+g_var.o_residues:
                coord_atomistic[chain][residue_id] = at_mod.check_hydrogens(coord_atomistic[chain][residue_id])
            if sys_type == 'PROTEIN':
                if res_index <= len(coord_atomistic[chain])-3:
                    coord_atomistic[chain][residue_id], cross_vector = fix_carbonyl_chiral(residue_id, g_var.cg_residues['PROTEIN'], coord_atomistic[chain][residue_id], False)
                elif res_index < len(coord_atomistic[chain]) and 'cross_vector' in locals():
                    coord_atomistic[chain][residue_id], cross_vector = fix_carbonyl_chiral(residue_id, g_var.cg_residues['PROTEIN'], coord_atomistic[chain][residue_id], cross_vector)
            order = np.sort(np.array(list(coord_atomistic[chain][residue_id].keys())))
            for at_val, atom in enumerate(order):
                coord_atomistic[chain][residue_id][atom]['resid'] = res_index
                coords.append(coord_atomistic[chain][residue_id][atom]['coord'])
                final_at[chain][at_number]=coord_atomistic[chain][residue_id][atom]
                at_number+=1               
            final_at_residues[chain][res_index]=coord_atomistic[chain][residue_id]
        coords = at_mod.check_atom_overlap(coords)
        for atom in final_at[chain]:
            final_at[chain][atom]['coord']=coords[atom]
            if not skip:
                x, y, z = gen.trunc_coord(final_at[chain][atom]['coord'])
                pdb_output.write(g_var.pdbline%((atom,final_at[chain][atom]['atom'],final_at[chain][atom]['res_type'],' ',\
                                final_at[chain][atom]['resid'],x,y,z,1,0))+'\n')
    if 'pdb_output' in locals():
        pdb_output.close()
    return final_at_residues

########################################### fix carbonyl section 

def get_crossvector(cg, residue_id):
    ca=[]
    res_off = 0
    for x in range(2):
        for bead in cg[residue_id+res_off]:
            resname = cg[residue_id+res_off][next(iter(cg[residue_id+res_off]))]['residue_name']
            if bead in g_var.res_top[resname]['CONNECT']:
                if x == 0 :
                    ca.append(cg[residue_id+res_off][bead]['coord'])
                for con_val, con_dir in enumerate(g_var.res_top[resname]['CONNECT'][bead]['dir']):
                    if con_dir > 0:
                        ca.append(cg[residue_id+con_dir+res_off][g_var.res_top[resname]['CONNECT'][bead]['Con_Bd'][con_val]]['coord'])
                        res_off = con_dir      
    return ca

def get_chiral_carbonyl(at):
    chiral = {}
    carbonyl = {}
    for atom in at:
        if at[atom]['atom'] in ['N','CA', 'C', 'O']:
            carbonyl[at[atom]['atom']] = atom
        if at[atom]['atom'] in g_var.res_top[at[atom]['res_type']]['CHIRAL']['atoms']:
            chiral[at[atom]['atom']] = atom
    return chiral, carbonyl

def correct_carbonyl_alignment(at, cross_vector, carbonyl):
    rotation = at_mod.align_to_vector(at_mod.noramlised_vector(at[carbonyl['O']]['coord'],at[carbonyl['C']]['coord']), cross_vector)
    at[carbonyl['O']]['coord'] = (at[carbonyl['O']]['coord']-at[carbonyl['C']]['coord']).dot(rotation)+at[carbonyl['C']]['coord']
    at[carbonyl['C']]['coord'] = at[carbonyl['C']]['coord'] + cross_vector*0.2
    at[carbonyl['N']]['coord'] = at[carbonyl['N']]['coord'] - cross_vector*0.5
    return at

def correct_protein_chiral(at, chiral):
    for chiral_group in g_var.res_top[at[1]['res_type']]['CHIRAL']:
        if chiral_group != 'atoms':            
            atom_list = [chiral[chiral_group]]
            for chiral_atoms in ['c1','c2','c3']:
                atom_list.append(chiral[g_var.res_top[at[1]['res_type']]['CHIRAL'][chiral_group][chiral_atoms]])
            cross_vector_chiral = at_mod.find_cross_vector( [at[atom_list[3]]['coord'], at[atom_list[2]]['coord'], at[atom_list[1]]['coord']])
            at[atom_list[0]]['coord'] = at[atom_list[0]]['coord'] + cross_vector_chiral*0.5
            if g_var.res_top[at[1]['res_type']]['CHIRAL'][chiral_group]['m'] in chiral:
                m = chiral[g_var.res_top[at[1]['res_type']]['CHIRAL'][chiral_group]['m']]
                at[m]['coord'] = at[m]['coord'] + cross_vector_chiral*1
    return at

def fix_carbonyl_chiral(residue_id, cg, at, cross_vector):
    chiral, carbonyl = get_chiral_carbonyl(at)
    if not np.any(cross_vector):
        ca = get_crossvector(cg, residue_id)
        cross_vector = at_mod.find_cross_vector( ca )
    at = correct_carbonyl_alignment(at, cross_vector, carbonyl)
    if len(g_var.res_top[at[1]['res_type']]['CHIRAL']) >= 2:
        at = correct_protein_chiral(at, chiral)
    return at, cross_vector

### align sequences of user AT structure and CG input

def add_to_sequence(sequence, residue, chain_count):
    if residue in g_var.o_residues:
        if residue in g_var.other:
            sequence[chain_count]+=g_var.other[residue]
        else:
            sequence[chain_count]+=residue
    elif residue in g_var.aas:
        sequence[chain_count]+=g_var.aas[residue]
    else:
        sequence[chain_count]+='X'    
    return sequence

def check_sequence():
    g_var.seq_at['PROTEIN'] = {}
    for chain in range(len(g_var.atomistic_protein_input_raw)):
        g_var.seq_at['PROTEIN'][chain]=[]
        for resid in g_var.atomistic_protein_input_raw[chain]:
            for atom in g_var.atomistic_protein_input_raw[chain][resid]:
                g_var.seq_at['PROTEIN'] = add_to_sequence(g_var.seq_at['PROTEIN'], g_var.atomistic_protein_input_raw[chain][resid][atom]['res_type'], chain)
                break

def align_chain_sequence(sys_type):
    cg_sequence = copy.deepcopy(g_var.seq_cg) 
    for chain_at in range(len(g_var.atomistic_protein_input_raw)):
        skip_sequence=False
        chain_cg=0
        s = difflib.SequenceMatcher(None, g_var.seq_at[sys_type][chain_at], cg_sequence[sys_type][chain_cg], autojunk=False)
        seq_info = s.get_matching_blocks()
        while seq_info[0][2] != len(g_var.seq_at[sys_type][chain_at]):
            if chain_cg >= len(cg_sequence[sys_type])-1:
                print('\nCannot find a match for user supplied chain: '+str(chain_at))#+'\n\nAtomistic chain:\n'+str(seq_user[chain_at]),'\n\nIn CG:\n'+str(sequence))
                skip_sequence = True
                break
            else:
                chain_cg+=1
                s = difflib.SequenceMatcher(None, g_var.seq_at[sys_type][chain_at], cg_sequence[sys_type][chain_cg], autojunk=False)
                seq_info = s.get_matching_blocks()
        if not skip_sequence:
            temp={}
            if chain_cg not in g_var.atomistic_protein_input_aligned:
                g_var.atomistic_protein_input_aligned[chain_cg]={}
            if seq_info[0][2] == len(g_var.seq_at[sys_type][chain_at]):
                for resid,  residue in enumerate(g_var.atomistic_protein_input_raw[chain_at]):
                    temp[resid + seq_info[0][1]] = g_var.atomistic_protein_input_raw[chain_at][residue]
                g_var.atomistic_protein_input_aligned[chain_cg][str(seq_info[0][1])+':'+str(seq_info[0][1]+seq_info[0][2])]=temp  
            cg_sequence[sys_type][chain_cg] = mask_sequence(cg_sequence[sys_type][chain_cg], seq_info[0][1], seq_info[0][1]+seq_info[0][2])
            g_var.cg_chain_group[chain_at]=chain_cg

    check_chain_alignment_coverage(sys_type)

    if len(g_var.atomistic_protein_input_raw) < len(g_var.seq_cg[sys_type]):
        print('### WARNING you have supplied fewer chains than exist in the CG system ###\n')
    if len(g_var.atomistic_protein_input_aligned) > 0:
        g_var.user_at_input = True
    else: 
        g_var.user_at_input = False
    sort_chains(sys_type)
    if g_var.args.v >= 2:
        print(gen.print_sequnce_info('PROTEIN'))

def sort_chains(sys_type):
    if g_var.group_chains is None:
        g_var.group_chains = {}
        for i in range(len(g_var.seq_at[sys_type])):
            g_var.group_chains[i] = i
    elif g_var.group_chains == 'chain':
        g_var.group_chains = g_var.cg_chain_group
    elif g_var.group_chains == 'all':
        g_var.group_chains = {}
        for i in range(len(g_var.cg_chain_group)):
            g_var.group_chains[i] = 0
    else:
        sys.exit('Failed to parse chain sorting')


def check_chain_alignment_coverage(sys_type):
    for key, chain in g_var.atomistic_protein_input_aligned.items():
        total = 0
        for fragment in chain:
            resid_range = fragment.split(':')
            total+=(int(resid_range[1])-int(resid_range[0]))
        if total == len(g_var.seq_cg[sys_type][key]):
            g_var.skip_disul[key]=True
        else:
            g_var.skip_disul[key]=False

def mask_sequence(sequence, st, end):
    for index, residue in enumerate(sequence):
        if st <= index < end:
            sequence[index]='-'
    return sequence

#####   aligning user protein chains

def align_user_chains(final_coordinates_atomistic):
    atomistic_protein_centered, cg_com = center_atomistic() ## centers each chain by center of mass
    at_com_group, cg_com_group = rotate_protein_monomers(atomistic_protein_centered, final_coordinates_atomistic, cg_com) 
    atomistic_protein_rotated = apply_rotations_to_chains(final_coordinates_atomistic, atomistic_protein_centered, 
                                                                at_com_group,cg_com_group,cg_com) ## apply rotation matrix to atoms and build in missing residues
    final_user_supplied_coord = correct_disulphide_bonds(atomistic_protein_rotated) ## fixes sulphur distances in user structure
    for chain, protein in final_user_supplied_coord.items():
        at_id=0
        coord=[]
        final_atom = []
        for resid, residue in protein.items():
            for atom in residue.values():
                coord.append(atom['coord'])
                final_atom.append({'atom_name':atom['atom'], 'residue_name':atom['res_type'], 'chain':chain, 'residue_id':resid,\
                                    'x':atom['coord'][0],'y':atom['coord'][1],'z':atom['coord'][2]})
                at_id+=1
        corrected_coords, index_conversion = at_mod.index_conversion_generate(final_atom, coord)
        at_mod.write_pdb(final_atom, corrected_coords, index_conversion, g_var.working_dir+'PROTEIN/PROTEIN_aligned_'+str(chain)+'.pdb')

def center_atomistic():
    cg_com={}
#### for each protein chain center on cg representation 
    for chain in g_var.atomistic_protein_input_aligned:
        cg_com[chain]=[]
        for part_val, part in enumerate(g_var.atomistic_protein_input_aligned[chain]):
            sls, sle= int(part.split(':')[0]),int(part.split(':')[1])
            if 'protein_mass' not in locals():
                protein_mass={}
            if chain in g_var.group_chains:
                if g_var.group_chains[chain] not in protein_mass:
                    protein_mass[g_var.group_chains[chain]]=fetch_backbone_mass(g_var.atomistic_protein_input_aligned[chain][part], [])
                else:
                    protein_mass[g_var.group_chains[chain]]=fetch_backbone_mass(g_var.atomistic_protein_input_aligned[chain][part], protein_mass[g_var.group_chains[chain]])
                if 'cg_backbone_masses' not in locals():
                    cg_backbone_masses = {}
                if g_var.group_chains[chain] not in cg_backbone_masses:
                    cg_backbone_masses[g_var.group_chains[chain]]= g_var.backbone_coords[chain][sls:sle]
                else:
                    cg_backbone_masses[g_var.group_chains[chain]] += g_var.backbone_coords[chain][sls:sle]                  
    g_var.atomistic_protein_input_aligned, cg_com = center_at_protein_chain_groups(g_var.atomistic_protein_input_aligned, cg_com, protein_mass, cg_backbone_masses)
    return g_var.atomistic_protein_input_aligned, cg_com

def center_at_protein_chain_groups(atomistic_protein_input, cg_com, protein_mass, cg_backbone_masses):
    for chain in atomistic_protein_input:
        for part_val, part in enumerate(atomistic_protein_input[chain]):
            sls, sle= int(part.split(':')[0]),int(part.split(':')[1])
            if chain in g_var.group_chains:
                atomistic_protein_mass = at_mod.COM(protein_mass[g_var.group_chains[chain]], 'AT protein chain: '+str(chain))
                cg_com[chain].append(at_mod.COM(cg_backbone_masses[g_var.group_chains[chain]], 'CG protein chain: '+str(chain)))
            else:
                protein_mass=fetch_backbone_mass(atomistic_protein_input[chain][part], [])
                atomistic_protein_mass = at_mod.COM(protein_mass, 'protein at: '+str(chain)+' '+part)#### add center of mass of CG_proteins
                cg_com[chain].append(at_mod.COM(g_var.backbone_coords[chain][sls:sle], 'protein cg: '+str(chain)+' '+part))
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
            if part[residue][atom]['atom'] in g_var.res_top[part[residue][atom]['res_type']]['ATOMS']:
                protein_mass.append([part[residue][atom]['coord'][0],part[residue][atom]['coord'][1],part[residue][atom]['coord'][2],part[residue][atom]['frag_mass']])    
    return protein_mass

def rotate_protein_monomers(atomistic_protein_centered, final_coordinates_atomistic, cg_com):
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

                        if atomistic_protein_centered[chain][part][residue][atom]['atom'] in g_var.res_top[atomistic_protein_centered[chain][part][residue][atom]['res_type']]['ATOMS']:
                            at_centers_iter.append(np.append(atomistic_protein_centered[chain][part][residue][atom]['coord'],atomistic_protein_centered[chain][part][residue][atom]['frag_mass']))
                    try:
                        at_centers.append(np.average(np.array(at_centers_iter)[:,:3], axis=0, weights=np.array(at_centers_iter)[:,3]))
                    except BaseException:
                        for atom in atomistic_protein_centered[chain][part][residue]:
                            print(atomistic_protein_centered[chain][part][residue][atom])
                        sys.exit()
                if len(at_centers) == len(np.array(g_var.backbone_coords[chain])[sls:sle,:3]):
                    at_com_group, cg_com_group = return_grouped_rotations(chain, part_val, at_centers, cg_com, at_com_group, cg_com_group, sls, sle)
                else:
                    sys.exit('In chain '+str(chain)+' the atomistic input does not match the CG. \n\
                            number of CG residues '+str(len(g_var.backbone_coords[chain]))+'\nnumber of AT residues '+str(len(at_centers)))
    return at_com_group, cg_com_group

def apply_rotations_to_chains(final_coordinates_atomistic, atomistic_protein_centered, at_com_group,cg_com_group,cg_com):
    final_user_supplied_coord={}
    for chain in range(len(final_coordinates_atomistic)):
        rotations=[]
        if chain in atomistic_protein_centered:
            rotate_chain=at_mod.kabsch_rotate(at_com_group[g_var.group_chains[chain]]-cg_com[g_var.group_chains[chain]][0], cg_com_group[g_var.group_chains[chain]]-cg_com[g_var.group_chains[chain]][0])
            for part_val, part in enumerate(atomistic_protein_centered[chain]):
                if chain in g_var.group_chains:
                    rotations.append(rotate_chain)
                else:
                    rotations = at_com_group[chain]
            final_user_supplied_coord[chain] = hybridise_protein_inputs(final_coordinates_atomistic[chain], atomistic_protein_centered[chain], cg_com[chain], rotations, chain)
        else:
            final_user_supplied_coord[chain] = hybridise_protein_inputs(final_coordinates_atomistic[chain], [], [], [], chain)
    return final_user_supplied_coord

def return_indivdual_rotations(chain, part_val, at_centers, cg_com, at_com_group, cg_com_group, sls, sle):
    if chain not in at_com_group:
        at_com_group[chain]=[]
    if len(at_centers) == len(np.array(g_var.backbone_coords[chain])[sls:sle,:3]):
        at_com_group[chain].append(at_mod.kabsch_rotate(np.array(at_centers)-cg_com[chain][part_val], 
                             np.array(g_var.backbone_coords[chain])[sls:sle,:3]-cg_com[chain][part_val]))
    else:
        sys.exit('In chain '+str(chain)+' the atomistic input does not match the CG. \n\
                    number of CG residues '+str(len(g_var.backbone_coords[chain]))+'\nnumber of AT residues '+str(len(at_centers)))    
    return at_com_group

def return_grouped_rotations(chain, part_val, at_centers, cg_com, at_com_group, cg_com_group, sls, sle):
    if chain in g_var.group_chains:
        if g_var.group_chains[chain] not in at_com_group:
            at_com_group[g_var.group_chains[chain]]=[]   
            at_com_group[g_var.group_chains[chain]]=np.array(at_centers)
            cg_com_group[g_var.group_chains[chain]]=np.array(g_var.backbone_coords[chain])[sls:sle,:3]  
        else:
            at_com_group[g_var.group_chains[chain]]=np.append(at_com_group[g_var.group_chains[chain]], np.array(at_centers), axis=0)
            cg_com_group[g_var.group_chains[chain]]=np.append(cg_com_group[g_var.group_chains[chain]], np.array(g_var.backbone_coords[chain])[sls:sle,:3], axis=0)                                                                               
    else:
        at_com_group = return_indivdual_rotations(chain, part_val, at_centers, cg_com, at_com_group, cg_com_group, sls, sle)
    return at_com_group, cg_com_group

def hybridise_protein_inputs(final_coordinates_atomistic, atomistic_protein_centered, cg_com, xyz_rot_apply, chain):
    complete_user_at = {}
    for residue in final_coordinates_atomistic:
        exists=False
        resname = final_coordinates_atomistic[residue][next(iter(final_coordinates_atomistic[residue]))]['res_type']
        if resname in g_var.mod_residues:
            complete_user_at[residue]=final_coordinates_atomistic[residue]
        elif resname not in g_var.mod_residues:
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

################################################################## Merge chains ####################    

def correct_amide_h(lines, coords):
    for at_val, atom in enumerate(lines):
        if atom['atom_name'] == 'C': 
            C = np.array(coords[at_val])
        if atom['atom_name'] == 'O':   
            O = np.array(coords[at_val]) 

        resname = gen.check_alternate_resname(atom['residue_name'])
        if atom['atom_name'] == g_var.res_top[resname]['amide_h'] and atom['residue_id'] != 0:  
            HN = np.array(coords[at_val]) 
            HN_index = at_val
        if atom['atom_name'] == 'N' and atom['residue_id'] != 0:
            N = np.array(coords[at_val])
        if 'C' in locals() and 'N' in locals() and 'O' in locals() and 'HN' in locals():
            O_C = O-C
            coords[at_val] = N - O_C
            del N, C, O, HN 
    return coords

def create_disres(coord, chain, file, at_start, count):
    P_R = np.array([])
    if chain in g_var.atomistic_protein_input_aligned:
        for key in g_var.atomistic_protein_input_aligned[chain].keys():
            P_R = np.append(P_R, np.arange(int(key.split(':')[0]), int(key.split(':')[1])+1))
        HN, O = get_backbone(coord)
        count, disres = find_connect(coord, P_R, HN, O, at_start, count)
        write_disres(disres)
    return len(coord)+at_start, count

def write_disres(disres):        
    header = True if not os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN_disres.itp') else False
    with open(g_var.working_dir+'PROTEIN/PROTEIN_disres.itp', 'a') as disres_out:
        if header:
            disres_out.write(';backbone hydrogen bonding distance restraints\n\n')
            disres_out.write('[ intermolecular_interactions ]\n[ distance_restraints ]\n')
            disres_out.write(';   i     j type label      funct         lo        up1        up2     weight\n')
        for restraint in disres:
            disres_out.write(restraint)

def get_backbone(coord):
    HN, O = [],[]
    for atom in coord:
        if atom['residue_name'] in g_var.p_residues and atom['atom_name'] == g_var.res_top[atom['residue_name']]['amide_h']:
            HN.append([int(atom['atom_number']), atom['x'], atom['y'], atom['z']])
        if atom['residue_name'] in g_var.p_residues and atom['atom_name'] == 'O':
            O.append([int(atom['atom_number']), atom['x'], atom['y'], atom['z']])
    return np.array(HN), np.array(O)     
            
def find_connect(coord, P_R, HN, O, at_start, count):
    disres = []
    tree = cKDTree(HN[:,1:])
    for carbonyl in O:
        ndx = tree.query_ball_point(carbonyl[1:], r=3)
        if len(ndx) > 0:
            for at in ndx:
                HN_resid = coord[int(HN[at][0])-1]['residue_id']
                O_resid = coord[int(carbonyl[0])-1]['residue_id']
                if HN_resid in P_R and O_resid in P_R and HN_resid not in [O_resid-1 ,O_resid, O_resid+1]:
                    count+=1
                    xyz1 = [coord[int(carbonyl[0])-1]['x'], coord[int(carbonyl[0])-1]['y'], coord[int(carbonyl[0])-1]['z']]
                    xyz2 = [coord[int(HN[at][0])-1]['x'], coord[int(HN[at][0])-1]['y'], coord[int(HN[at][0])-1]['z']]
                    dist = (gen.calculate_distance(xyz1, xyz2)/10)-0.05
                    disres.append('{0:10}{1:10}{2:3}{3:12}{4:12}{5:^12}{6:14}{7:14}{8:5}\n'.format(str(at_start+int(HN[at][0])), str(at_start+int(carbonyl[0])), 
                                                                                                         '1', str(count),'2', '0', str(np.round(dist,4)), str(np.round(dist+0.01,4)), '1'))
    return count, disres 

########################################################### RMSD

def write_RMSD():
    print('\nCalculating backbone RMSDs\n')
    RMSD={}
    de_novo_atoms, chain_count = read_in.read_in_atomistic(g_var.final_dir+'final_cg2at_de_novo.pdb') ## reads in final pdb
    if chain_count != g_var.system['PROTEIN']:
        sys.exit('number of chains in atomistic protein input ('+str(chain_count)+') does not match CG representation ('+str(g_var.system['PROTEIN'])+')')
    RMSD_de_novo = RMSD_measure_de_novo(de_novo_atoms) ## gets rmsd of de novo

    if g_var.user_at_input and 'PROTEIN' in g_var.cg_residues: 
        if g_var.args.o in ['all', 'align']: 
            at_input_atoms, chain_count = read_in.read_in_atomistic(g_var.final_dir+'final_cg2at_aligned.pdb')
            RMSD_aligned = RMSD_measure_de_novo(at_input_atoms)   
            seg_rmsd = RMSD_measure_aligned(at_input_atoms)   

    with open(g_var.final_dir+'structure_quality.dat', 'w') as qual_out:   
        line_1   = ' chain    De novo BB RMSD ('+chr(197)+')'
        line_2   = ' -----    -------------------'
        if 'seg_rmsd' in locals():
            line_1+= '    Aligned BB RMSD ('+chr(197)+')    Seg backbone RMSD ('+chr(197)+')'
            line_2+= '    -------------------    ---------------------'
        qual_out.write(line_1+'\n'+line_2+'\n')
        print(line_1+'\n'+line_2)
        write_warning = False
        for chain in RMSD_de_novo:
            line = ' {0:^5}{1:^28}'.format(str(chain), float(RMSD_de_novo[chain]))
            if 'seg_rmsd' in locals():
                line += '{0:^18}'.format(float(RMSD_aligned[chain]))
                if chain in seg_rmsd:
                    line+= '{0:^28}'.format(', '.join(seg_rmsd[chain])) 
                    if '*' in ', '.join(seg_rmsd[chain]):
                        write_warning = True
            print(line)
            qual_out.write(line+'\n')
        if write_warning:
            print('\n * Segment alignments may have minor deviations due to either clashes or structure hybridisation.')
            qual_out.write('\n * Segment alignments may have minor deviations due to either clashes or structure hybridisation.')
        print('\n\nAll RMSDs have been saved in: \n'+g_var.final_dir+'structure_quality.dat\n')

def RMSD_measure_de_novo(structure_atoms):
    RMSD_dict = {}
    for chain in range(g_var.system['PROTEIN']):
        at_centers=[]
    #### runs through every residue and atom  
        for residue in structure_atoms[chain]:
        #### gets center of mass of each residue (note only backbone heavy atoms have a mass)
            at_centers_iter=[]
            for atom in structure_atoms[chain][residue].values():
                if atom['atom'] in g_var.res_top[atom['res_type']]['ATOMS']:
                    at_centers_iter.append(np.append(atom['coord'],atom['frag_mass']))
            if np.any(np.array(at_centers_iter)[:,3] > 0):
                at_centers.append(np.average(np.array(at_centers_iter)[:,:3], axis=0, weights=np.array(at_centers_iter)[:,3]))
            else:
                sys.exit('Missing masses for RMSD')

    #### checks that the number of residues in the chain are the same between CG and AT
        if len(at_centers) != len(g_var.backbone_coords[chain]):
            sys.exit('In chain '+str(chain)+' the atomistic input does not match the CG. \n\
    number of CG residues '+str(len(g_var.backbone_coords[chain]))+'\nnumber of AT residues '+str(len(at_centers)))

        cg_center = np.mean(np.array(g_var.backbone_coords[chain])[:,:3], axis=0)
        at_align = np.array(at_centers) - (np.mean(np.array(at_centers), axis=0) - cg_center)
        at_align = RMSD_align(at_align, np.array(g_var.backbone_coords[chain])[:,:3])
        RMSD_dict[chain]=Calculate_RMSD(np.array(at_align), np.array(g_var.backbone_coords[chain])[:,:3])  #### stores RMSD in dictionary
    return RMSD_dict

def RMSD_align(coord_set_1, coord_set_2):
    center = np.mean(coord_set_2, axis=0)
    xyz_rot_apply=at_mod.kabsch_rotate(coord_set_1-center, coord_set_2-center)
    ali= []
    for at_val, atom in enumerate(coord_set_1):
        ali.append( at_mod.rotate_atom(atom, center, xyz_rot_apply) )
    return np.array(ali)

def Calculate_RMSD(C1, C2):
    dist=np.sqrt((C1 - C2)**2)
    return np.round(np.sqrt(np.mean(dist**2)),3) #### RMSD calculation

def get_coordinates(input_coord, P_R, chain):
    coord_dict = {}
    for residue_val, residue in input_coord[chain].items():
        for seg_val, seg in enumerate(P_R):
            if seg_val not in coord_dict:
                coord_dict[seg_val]=[]
            if residue_val in seg:
                for atom in residue.values():
                    if atom['atom'] in g_var.res_top[atom['res_type']]['ATOMS'] and atom['res_type'] in g_var.p_residues:
                        coord_dict[seg_val].append(np.append(atom['coord'],atom['frag_mass']))
                break
    for segment in range(len(coord_dict)):
        coord_dict[segment] = np.array(coord_dict[segment])
    return coord_dict

def RMSD_measure_aligned(Final_structure):
    seg_rmsd={}
    total_initial, total_final = np.array([]),np.array([])
    for chain in g_var.atomistic_protein_input_aligned:
        initial_structure, chain_count = read_in.read_in_atomistic(g_var.working_dir+'PROTEIN/MIN/PROTEIN_aligned_'+str(chain)+'.pdb')
        P_R = []
        for key in g_var.atomistic_protein_input_aligned[chain].keys():
            P_R.append(np.arange(int(key.split(':')[0])-1, int(key.split(':')[1])))
        final_backbone = get_coordinates(Final_structure, P_R, chain)
        initial_backbone = get_coordinates(initial_structure, P_R, 0)
        seg_rmsd[chain]=[]
        for segment in range(len(final_backbone)):
            initial_backbone_fitted = RMSD_align(initial_backbone[segment][:,:3],final_backbone[segment][:,:3])
            total_initial = np.append(total_initial, initial_backbone_fitted)
            total_final = np.append(total_final, final_backbone[segment][:,:3])
            RMSD_val = Calculate_RMSD(np.array(initial_backbone_fitted), final_backbone[segment][:,:3])
            if RMSD_val > 0.15:
                seg_rmsd[chain].append(str(RMSD_val)+' *')
            else:
                seg_rmsd[chain].append(str(RMSD_val))
    # print(Calculate_RMSD(total_initial, total_final))
    return seg_rmsd
