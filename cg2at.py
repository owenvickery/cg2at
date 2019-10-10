import os, sys
import numpy as np
from subprocess import Popen, PIPE
import subprocess, shlex
from time import gmtime, strftime
import math
import multiprocessing as mp
import argparse
import copy
from shutil import copyfile
from distutils.dir_util import copy_tree
import time
# from numba import njit, jit
from string import ascii_uppercase
from pathlib import Path
import re
import datetime

parser = argparse.ArgumentParser(description='Converts CG representation into atomistic', epilog='Enjoy the program and best of luck!', allow_abbrev=True)
parser.add_argument('-c', help='coarse grain coordinates',metavar='pdb/gro',type=str, required=True)
parser.add_argument('-a', help='atomistic coordinates',metavar='pdb/gro',type=str)
parser.add_argument('-l', help='additional fragment library location',type=str)
parser.add_argument('-v', help='verbose', action='store_true')
parser.add_argument('-gromacs', help='prints gromacs commands as well', action='store_true')
parser.add_argument('-f', help='additional forcefield location',type=str)
parser.add_argument('-steer', help='do not run steered MD on atomistic structure (requires atomistic structure)', action='store_false')
args = parser.parse_args()
options = vars(args)




def read_forcefield():
#### Read in forcefields provided
	forcefield_available_user=[] ## needed to calculate length of list of 0
	for root, dirs, files in os.walk(script_dir+'database/forcefields'):
			forcefield_available_prov=dirs
			break
#### Read in forcefields user provided
	if args.f!= None:
		for root, dirs, files in os.walk(args.f+'/database/forcefields'):
			forcefield_available_user=dirs
			break
	return 	forcefield_available_prov, forcefield_available_user

def forcefield_selection(forcefield_available_prov, forcefield_available_user):
#### print out selection of forcefields
	print('\n\n\tselection\t\tforcefield\n')
	print('\t---------\t\t----------\n')
	print('\t\tprovided forcefields\n')
	for force_num_prov, line in enumerate(forcefield_available_prov):
		print('\t',force_num_prov,'\t\t',line[:-3])
	if args.f != None:
		print('\nUser provided forcefields\n')
		for force_num_user, line in enumerate(forcefield_available_user):
			print('\t',force_num_user+force_num_prov,'\t\t',line[:-3])	
#### ask which forcefield to use
	while True:
		try:
			if len(forcefield_available_prov)+len(forcefield_available_user)==1:
				print('\nOnly 1 forcefield currently available, therefore you have no choice but to accept the following choice.')
				return 0 
			number = int(input("\nplease select a forcefield: "))
			if number < len(forcefield_available_prov)+len(forcefield_available_user):
				return number
		except KeyboardInterrupt:
			sys.exit('\nInterrupted')
		except:
			print("Oops!  That was a invalid choice")

def sort_forcefield(forcefield_available_prov, forcefield_available_user, f_number):
#### returns forcefield location and forcefield name
#### if forcefield selection is in provided copy forcefield to FORCEFIELD and FINAL directories
	if f_number < len(forcefield_available_prov):
		print('\nYou have selected: '+forcefield_available_prov[f_number][:-3]+'\n\nGood luck\n')
		copy_tree(script_dir+'/database/forcefields/'+forcefield_available_prov[f_number], working_dir+'FORCEFIELD/'+forcefield_available_prov[f_number]+'/.')
		copy_tree(script_dir+'/database/forcefields/'+forcefield_available_prov[f_number], final_dir+forcefield_available_prov[f_number]+'/.')
		return script_dir+'/database/forcefields/', forcefield_available_prov[f_number][:-3]
#### if forcefield selection is in user provided copy forcefield to FORCEFIELD and FINAL directories
	else:
		print('\nYou have selected: '+forcefield_available_user[f_number-len(forcefield_available_prov)][:-3]+'\n\nGood luck\n')
		copy_tree(args.f+'/database/forcefields/'+forcefield_available_user[f_number-len(forcefield_available_prov)], working_dir+'FORCEFIELD/'+forcefield_available_user[f_number-len(forcefield_available_prov)]+'/.')
		copy_tree(args.f+'/database/forcefields/'+forcefield_available_user[f_number-len(forcefield_available_prov)], final_dir+forcefield_available_user[f_number-len(forcefield_available_prov)]+'/.')		
		return args.f+'/database/forcefields/', forcefield_available_user[f_number-len(forcefield_available_prov)][:-3]

def fetch_residues():
#### list of directories and water types  [[root, folders...],[root, folders...]]
	np_directories, p_directories, water=[[]], [[]],[]
	for directory_type in ['/non_protein/', '/protein/']:
#### adds non protein residues locations to np_directories
		for root, dirs, files in os.walk(script_dir+'database/'+forcefield+directory_type):
			if directory_type =='/non_protein/':
				np_directories[-1].append(root)
				np_directories[-1]+=dirs
#### adds water residues to water list
				if os.path.exists(script_dir+'database/'+forcefield+directory_type+'SOL'):
					for root, dirs, files in os.walk(script_dir+'database/'+forcefield+directory_type+'SOL'):
						for water_type in files:
							water.append(water_type[:-4])
						break
#### adds protein residues locations to p_directories
			else:
				p_directories[-1].append(root)
				p_directories[-1]+=dirs	
				if os.path.exists(script_dir+'database/'+forcefield+directory_type+'MOD'):
					p_directories.append([])
					for root, dirs, files in os.walk(script_dir+'database/'+forcefield+directory_type+'MOD'):
						p_directories[-1].append(root)
						p_directories[-1]+=dirs
						break
			break
#### adds non protein residues locations to np_directories from user defined
		if args.l != None:
			if os.path.exists(args.l+'/database/'+forcefield):
				for root, dirs, files in os.walk(args.l+'/database/'+forcefield+directory_type):
					if directory_type =='/non_protein/':
						np_directories[-1].append(root)
						np_directories[-1]+=dirs
#### adds water residues to water list
						if os.path.exists(args.l+'/database/'+forcefield+directory_type+'SOL'):
							for root, dirs, files in os.walk(script_dir+'/database/'+forcefield+directory_type+'SOL'):
								for water_type in files:
									water.append(water_type[:-4])
								break
#### adds protein residues locations to p_directories from user defined
					else:
						p_directories[-1].append(root)
						p_directories[-1]+=dirs	
						if os.path.exists(args.l+'database/'+forcefield+directory_type+'MOD'):
							p_directories.append([])
							for root, dirs, files in os.walk(args.l+'database/'+forcefield+directory_type+'MOD'):
								p_directories[-1].append(root)
								p_directories[-1]+=dirs
								break
					break
			else:
				os.exit('Cannot find :'+args.l+'/database/forcefield')
#### sorts directories alphabetically and creates residue database
	p_residues, np_residues = [],[]
	for directory in range(len(p_directories)):
		p_directories[directory].sort()
		p_residues+=p_directories[directory][1:]
	for directory in range(len(np_directories)):
		np_directories[directory].sort()
		np_residues+=np_directories[directory][1:]
#### if verbose prints all fragments found
	if args.v:
		print('\nnon protein residues fragment directories found: \n\nroot file system\n')
		for directory in range(len(np_directories)):
			print(np_directories[directory][0],'\n\nresidues\n\n',np_directories[directory][1:], '\n')
		print('\nprotein residues fragment directories found: \n\nroot file system\n')
		for directory in range(len(p_directories)):
			print(p_directories[directory][0],'\n\nresidues\n\n',p_directories[directory][1:], '\n')
	return np_residues, p_residues, np_directories, p_directories, water

def fetch_fragment(p_directories):
#### fetches the Backbone heavy atoms and the connectivity with pre/proceeding residues 
	atom_list, bb_list=[], []
	backbone={}     ### dictionary of backbone heavy atoms and connecting atoms eg backbone['ASP'][atoms/b_connect]
	for directory in range(len(p_directories)):
		for residue in p_directories[directory][1:]:	
			with open(p_directories[directory][0]+residue+'/BB.pdb', 'r') as pdb_input:
				for line_nr, line in enumerate(pdb_input.readlines()):
					if line.startswith('ATOM'):
						line_sep = pdbatom(line)
						if 'H' not in line_sep['atom_name']:
							atom_list.append(line_sep['atom_name'])    ### list of heavy atoms
							if line_sep['backbone'] != 0:
								bb_list.append(line_sep['atom_name'])  ### connecting atoms
			backbone[residue]={'atoms':atom_list,'b_connect':bb_list}  ### adds heavy atoms and connecting atoms to backbone dictionary 
			atom_list, bb_list=[], []  ### resets residue lists of heavy atoms and connecting atoms 
#### if verbose prints out all heavy atoms and connecting atoms for each backbone
	if args.v:
		print('backbone atoms for each residue and connecting atoms:\n')
		for residue in backbone:
			print(residue, '\tbackbone atoms:', backbone[residue]['atoms'], '\n\tbackbone connecting atoms:', backbone[residue]['b_connect'],'\n')
	return backbone

def make_min(residue, fragments):
#### makes minimisation folder
	if not os.path.exists('min'):
		os.mkdir('min')
#### makes em.mdp file for each residue
	if not os.path.exists('em_'+residue+'.mdp'):
		with open('em_'+residue+'.mdp','w') as em:
			em.write('define = -DPOSRE_CA\n integrator = steep\nnsteps = 10000\nemtol = 1500\nemstep = 0.001\ncutoff-scheme = Verlet\n')
#### if the residue is not in ['SOL', 'PROTEIN', 'ION'] at translational center of mass removal is used on the fragments
			if residue not in  ['SOL', 'PROTEIN', 'ION']:
				em.write('comm-mode = Linear\nnstcomm = 1\ncomm-grps = ')
				for frag in fragments:
					em.write(frag+' ')

### runs gromacs commands
def gromacs(cmd):
#### if the flag gromacs is used every gromacs command will be printed to the terminal 
	if args.gromacs:
		print('\nrunning gromacs: \n '+cmd+'\n')
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	err, out = output.communicate()
	exitcode = output.returncode
	out=out.decode("utf-8")
#### all gromacs outputs will be saved into gromacs_outputs within the folder it is run
	with open('gromacs_outputs', 'a') as checks:
		checks.write(out)
#### standard catch for failed gromacs commands
		if 'File input/output error:' in out:
			sys.exit('\n'+out)
		elif 'Error in user input:' in out:
			sys.exit('\n'+out)
		elif 'did not converge to Fmax ' in out:
			sys.exit('\n'+out)
		elif 'Segmentation fault (core dumped):' in out:
			sys.exit('\n'+out)
		elif 'Fatal error:' in out:
			sys.exit('\n'+out)


def pdbatom(line):
### get information from pdb file
### atom number, atom name, residue name,chain, resid,  x, y, z, backbone (for fragment), connect(for fragment)
	try:
		return dict([('atom_number',str(line[7:11]).replace(" ", "")),('atom_name',str(line[12:16]).replace(" ", "")),('residue_name',str(line[17:21]).replace(" ", "")),('chain',line[21]),('residue_id',int(line[22:26])), ('x',float(line[30:38])),('y',float(line[38:46])),('z',float(line[46:54])), ('backbone',int(float(line[56:62]))),('connect',int(float(line[62:67])))])
	except:
		sys.exit('\npdb line is wrong:\t'+line) 

def groatom(line):
### get information from gro file
### atom name, res name, chain, res id, x, y, z     
	try:
		return dict([('atom_name',str(line[10:15]).replace(" ", "")),('residue_name',str(line[5:10]).replace(" ", "")),('chain',line[21]),('residue_id',int(line[:5])), ('x',float(line[20:28])*10),('y',float(line[28:36])*10),('z',float(line[36:44])*10)])
	except:
		sys.exit('\ngro line is wrong:\t'+line) 

def read_initial_pdb():
#### initialisation of dictionaries etc
	cg_residues={}  ## dictionary of CG beads eg cg_residues[residue type(POPE)][resid(1)][bead name(BB)][residue_name(PO4)/coordinates(coord)]
	residue_list={} ## a dictionary of bead in each residue eg residue_list[bead name(BB)][residue_name(PO4)/coordinates(coord)]
	initial=True  ## if 1st atom
	box_line="CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00 P 1           1\n"  ## box vectors format for pdbs
	count=0  ### residue counter initialisation
	with open(args.c, 'r') as pdb_input:
		for line in pdb_input.readlines():
#### separate line dependant on gro or pdb
			run=False ## turns to true is line is a bead/atom
			if args.c.endswith('gro') and len(line.split())==6:
				line_sep = groatom(line)
				run=True
			elif args.c.endswith('pdb') and line.startswith('ATOM'):
				line_sep = pdbatom(line)
				run=True
#### if a bead it will be read in			
			if run: 
#### set up resnames in dictionaries
				if line_sep['residue_name'] in p_residues: ## if in protein database 
					if 'protein' not in cg_residues:  ## if protein does not exist add to dict
						cg_residues['protein']={}
				elif line_sep['residue_name'] in np_residues or line_sep['residue_name'] in water:
					if line_sep['residue_name'] in water: ## renames waters to SOL saves headache later on
						line_sep['residue_name']='SOL'
					if line_sep['residue_name'] not in cg_residues: ## if residue type does not exist add to dict
						cg_residues[line_sep['residue_name']]={}
#### sets up previous resid id 
				if initial: 
					residue_prev=line_sep['residue_id'] 
					initial=False
#### if resid the same as previous line
				if residue_prev == line_sep['residue_id']:   ### if resid is the same as the previous line, it adds resname and coordinates to the atom name key in residue_list 
					residue_list[line_sep['atom_name']]={'residue_name':line_sep['residue_name'],'coord':np.array([line_sep['x'],line_sep['y'],line_sep['z']])}
					line_sep_prev=line_sep.copy()
#### if resids are different then the residue list is added to cg_residues
				else: 
					if line_sep_prev['residue_name'] not in p_residues:
						cg_residues[line_sep_prev['residue_name']][count]={} ### then create sub dictionary cg_residues[resname][count]
						cg_residues[line_sep_prev['residue_name']][count]=residue_list ### adds residue list to dictionary key cg_residues[resname][count]
					else:
						cg_residues['protein'][count]={} ### then create sub dictionary cg_residues['protein'][count]
						cg_residues['protein'][count]=residue_list ### adds residue list to dictionary key cg_residues['protein'][count]
#### updates dictionaries and counters
					residue_list={}  ### resets residue list
					count+=1 ### moves counter along to next residue
					residue_list[line_sep['atom_name']]={'residue_name':line_sep['residue_name'],'coord':np.array([line_sep['x'],line_sep['y'],line_sep['z']])} ### it adds resname and coordinates to the atom name key in residue_list
					residue_prev=line_sep['residue_id']   ### updates residue_prev with new resid
					line_sep_prev=line_sep.copy()
				run=False ### resets check if atom
#### finds box vectors
			if line.startswith('CRYST') and args.c.endswith('pdb') : ### collects box vectors from pdb
				box_vec=line
			if args.c.endswith('gro') and len(line.split())==3 and not initial: ### collects box vectors from gro
				box_vec=box_line%(float(line.split()[0])*10,float(line.split()[1])*10,float(line.split()[2])*10)
#### adds final residue to cg_residues in the same manner as above
	if line_sep['residue_name'] in p_residues: 
		if count not in cg_residues['protein']:
			cg_residues['protein'][count]={}
		cg_residues['protein'][count]=residue_list
	else:
		if count not in cg_residues[line_sep['residue_name']]:
			cg_residues[line_sep['residue_name']][count]={}
		cg_residues[line_sep['residue_name']][count]=residue_list
#### checks if box vectors exist
	if 'box_vec' not in locals():### stops script if it cannot find box vectors
		sys.exit('missing box vectors')

	return cg_residues, box_vec

def fragment_location(residue, fragment):  
#### runs through dirctories looking for the atomistic fragments returns the correct location
	if residue in p_residues:
		for directory in range(len(p_directories)):
			if os.path.exists(p_directories[directory][0]+residue+'/'+fragment):
				return p_directories[directory][0]+residue+'/'+fragment
	else:
		for directory in range(len(np_directories)):
			if os.path.exists(np_directories[directory][0]+residue+'/'+fragment):
				return np_directories[directory][0]+residue+'/'+fragment
	sys.exit('cannot find fragment: '+residue+'/'+fragment)

def get_atomistic(residue,cg_fragment, cg_coord,resid):
#### find atomistic residues
	residue_list={} ## a dictionary of bead in each residue eg residue_list[atom number(1)][residue_name(ASP)/coordinates(coord)/atom name(C)/connectivity(2)/atom_mass(12)]
	frag_location=fragment_location(residue, cg_fragment+'.pdb') ### get fragment location from database
	fragment_masses=[] ### list [[coord, mass],[coord, mass]]
#### read in atomistic fragments into dictionary	
	with open(frag_location, 'r') as pdb_input:
		for line_nr, line in enumerate(pdb_input.readlines()):
			if line.startswith('ATOM'):
				line_sep = pdbatom(line) ## splits up pdb line
				residue_list[line_sep['atom_number']]={'coord':np.array([line_sep['x'],line_sep['y'],line_sep['z']]),'atom':line_sep['atom_name'], 'res_type':line_sep['residue_name'], 'connect':line_sep['connect'], 'frag_mass':1}
#### updates fragment mass   
				if 'H' not in line_sep['atom_name']:
					for atom in line_sep['atom_name']:
						if atom in mass:
							residue_list[line_sep['atom_number']]['frag_mass']=mass[atom]  ### updates atom masses with crude approximations
							fragment_masses.append([line_sep['x'],line_sep['y'],line_sep['z'],mass[atom]])
				else:
					fragment_masses.append([line_sep['x'],line_sep['y'],line_sep['z'],1])
#### aligns atomistic fragment to cg bead
	COM_vector=np.average(np.array(fragment_masses)[:,:3], axis=0, weights=np.array(fragment_masses)[:,3])-np.array(cg_coord['coord']) ### gets vector between COM of atoms in fragment and cg bead 
	for at_id, residue in enumerate(residue_list): ### runs through atoms in fragments and centers on the cg bead 
		residue_list[residue]['coord']=residue_list[residue]['coord']-COM_vector
	return residue_list


def build_atomistic_system(cg_residues, box_vec):
	# pdbline = "ATOM  %5d %4s %4s %1s%3d    %8.3f%8.3f%8.3f%6.2f%6.2f"
	system=[]
	ion_sep={}
	atomistic_fragments={}
	fragment, fragment_names=[],[]
#### for each residue type covert to atomistic except protein
	for residue_type in [key for value, key in enumerate(cg_residues) if key not in ['protein']]:
	#### reset counters for each residue type
		print('Converting residue type: ' +residue_type)
	#### adds a mulitpier for number of residues in CG bead (water=4)
		if residue_type == 'SOL':
			multiplier=water_mol_number
		else:
			multiplier=1
	#### creates folder for residue type
		if not os.path.exists(working_dir+residue_type):
			os.mkdir(working_dir+residue_type)
	#### fetches atoms for the residue type and centers then on cg bead
		atomistic_fragments[residue_type] = reorientate(residue_type, cg_residues[residue_type])
	#### if the residue type is in ['SOL', 'ION'] a single pdb is created
		if residue_type in ['SOL', 'ION']:
			pdb_output = open(working_dir+residue_type+'/mapping_'+residue_type+'.pdb', 'w')
			pdb_output.write('REMARK    GENERATED BY sys_setup_script\nTITLE     SELF-ASSEMBLY-MAYBE\nREMARK    Good luck\n\
'+box_vec+'MODEL        1\n')
	#### loop through all resids of that residue type 
		for resid in atomistic_fragments[residue_type]:

		#### if the residue type is not in ['SOL', 'ION'] a individual pdbs are created for each resid
			if residue_type not in ['SOL', 'ION']:
				pdb_output = open(working_dir+residue_type+'/'+residue_type+'_'+str(resid)+'.pdb', 'w')
				pdb_output.write('REMARK    GENERATED BY sys_setup_script\nTITLE     SELF-ASSEMBLY-MAYBE\nREMARK    Good luck\n\
'+box_vec+'MODEL        1\n')
		#### for every fragment in that resid
			for at_id, atom in enumerate(atomistic_fragments[residue_type][resid]):
				if resid==0 and residue_type not in ['ION', 'SOL']:
					if at_id==0:
						fragment.append('[ '+str(atomistic_fragments[residue_type][resid][at_id+1]['cg_bead'])+' ]')
						fragment_names.append(str(atomistic_fragments[residue_type][resid][at_id+1]['cg_bead']))
					if 'previous_bead' in locals():
						if previous_bead != atomistic_fragments[residue_type][resid][at_id+1]['cg_bead']:
							fragment.append('[ '+str(atomistic_fragments[residue_type][resid][at_id+1]['cg_bead'])+' ]')
							fragment_names.append(str(atomistic_fragments[residue_type][resid][at_id+1]['cg_bead']))
					fragment.append(str(at_id+1)+' ')
					previous_bead=atomistic_fragments[residue_type][resid][atom]['cg_bead']
				if residue_type=='ION':
					if atomistic_fragments[residue_type][resid][at_id+1]['res_type'] != 'SOL' and atomistic_fragments[residue_type][resid][at_id+1]['atom'] not in ion_sep:
						ion_sep[atomistic_fragments[residue_type][resid][at_id+1]['atom']]=0
					if atomistic_fragments[residue_type][resid][at_id+1]['res_type'] != 'SOL':
						ion_sep[atomistic_fragments[residue_type][resid][at_id+1]['atom']]
					pdb_output.write(pdbline%((int(at_id+1),atomistic_fragments[residue_type][resid][at_id+1]['atom'],atomistic_fragments[residue_type][resid][at_id+1]['res_type'],' ',1,\
					atomistic_fragments[residue_type][resid][at_id+1]['coord'][0],atomistic_fragments[residue_type][resid][at_id+1]['coord'][1],atomistic_fragments[residue_type][resid][at_id+1]['coord'][2],1,0))+'\n')
				else:
					pdb_output.write(pdbline%((int(atom),atomistic_fragments[residue_type][resid][at_id+1]['atom'],atomistic_fragments[residue_type][resid][at_id+1]['res_type'],' ',1,\
					atomistic_fragments[residue_type][resid][at_id+1]['coord'][0],atomistic_fragments[residue_type][resid][at_id+1]['coord'][1],atomistic_fragments[residue_type][resid][at_id+1]['coord'][2],1,0))+'\n')

#### works to here

		pdb_output.close()
		system.append([residue_type,(int(resid)+1)*multiplier])
		if residue_type == 'ION':
			os.chdir(working_dir+'ION')
			make_min(residue_type, [])
			os.chdir('..')
			copyfile(working_dir+residue_type+'/mapping_'+residue_type+'.pdb', working_dir+residue_type+'/min/'+residue_type+'_merged.pdb')
			system[-1][1]=[]
			for res_val, residue_type in enumerate(system):
				if residue_type[0]=='SOL':
					sol_index=res_val
			for ion in ion_sep:
				system[-1][1].append([ion,ion_sep[ion]])
				try:
					system[sol_index][1]=system[sol_index][1]+(ion_sep[ion]*4)
				except:
					system.append(['SOL', ion_sep[ion]*4])
		else:
			simulation_setup(resid, residue_type, fragment, fragment_names, multiplier)
	return system

def get_atomistic_non_protein(cg_residue_type,cg_residue, cg_resid):
	at_residues={}
	connect=[]
	for cg_bead in cg_residue:
		at_residues[cg_bead]=get_atomistic(cg_residue_type,cg_bead, cg_residue[cg_bead], cg_resid+1)
		if cg_residue_type not in ['SOL', 'ION']:
			for atom_num, atom in enumerate(at_residues[cg_bead]):
				if at_residues[cg_bead][atom]['connect'] > 0:
					connect.append([cg_bead,atom, at_residues[cg_bead][atom]['connect']]) 
	connect=np.array(connect)	

	return at_residues, connect


def reorientate(cg_residue_type,cg_residues):
	pdbline = "ATOM  %5d %4s %4s %1s%3d    %8.3f%8.3f%8.3f%6.2f%6.2f"
	atomistic_fragments={}
	for cg_resid, cg_residue in enumerate(cg_residues):
		initial, initial_second=True, False
		at_residues, connect = get_atomistic_non_protein(cg_residue_type,cg_residues[cg_residue], cg_resid)	
		atomistic_fragments[cg_resid]={}
		if cg_residue_type not in ['SOL', 'ION']:
			for bead_number, cg_bead in enumerate(cg_residues[cg_residue]):
				at_connections,cg_connections=[],[]
				run=np.where(connect[:,0]==str(cg_bead))
				center=cg_residues[cg_residue][cg_bead]['coord']
				repeat=[]
				cg_temp=[]
				for val,con in enumerate(connect):
					repeat = len(np.where(connect[:,2]==con[2])[0])
					if con[2] in connect[run][:,2] and connect[val][0] != cg_bead:
						cg_temp.append(cg_residues[cg_residue][connect[val][0]]['coord']-center)
						if len(cg_temp) == repeat-1:
							cg_temp=np.array(cg_temp)
							cg_connections.append(np.mean(cg_temp, axis=0))
							cg_temp=[]
					if con[0] == cg_bead and con[2] in connect[run][:,2]:
						at_connections.append(at_residues[con[0]][con[1]]['coord']-center)

				xyz_rot_apply=rotate(np.array(at_connections), np.array(cg_connections))
				for atom in at_residues[cg_bead]:
					at_residues[cg_bead][atom]['coord']=rotate_atom(at_residues[cg_bead][atom]['coord'], center, xyz_rot_apply)
					at_residues[cg_bead][atom].update({'cg_bead':bead_number+1})
					atomistic_fragments[cg_resid][int(atom)]=at_residues[cg_bead][atom]
		else:
			for cg_bead in at_residues:
				for atom in at_residues[cg_bead]:
					atomistic_fragments[cg_resid][int(atom)]=at_residues[cg_bead][atom]
	return atomistic_fragments


def build_protein_atomistic_system(cg_residues, box_vec):
	chain_count=0
	backbone_coords={} 
	at_counter=1
	system=[]
	ion_sep={}
	atomistic_fragments={}
	print('Converting Protein')
	if not os.path.exists(working_dir+'PROTEIN'):
		os.mkdir(working_dir+'PROTEIN')		
	pdb_output = open(working_dir+'PROTEIN/PROTEIN_novo_'+str(chain_count)+'.pdb', 'w')
	pdb_output.write('REMARK    GENERATED BY sys_setup_script\nTITLE     SELF-ASSEMBLY-MAYBE\nREMARK    Good luck\n\
'+box_vec+'MODEL        1\n')
	at_residues={}
	backbone_coords[chain_count]=[]
	for residue_number in cg_residues:
		connect=[]
		BB_connect=[]
		final_at_residues={}
		at_residues[residue_number]={}
		for cg_fragments in cg_residues[residue_number]:
			at_residues[residue_number][cg_fragments]=get_atomistic(cg_residues[residue_number][cg_fragments]['residue_name'],cg_fragments, cg_residues[residue_number][cg_fragments], residue_number)
			for atom_num, atom in enumerate(at_residues[residue_number][cg_fragments]):
				if at_residues[residue_number][cg_fragments][atom]['connect'] > 0:
					connect.append([cg_fragments,atom, at_residues[residue_number][cg_fragments][atom]['connect']]) 
				if cg_fragments=='BB' and at_residues[residue_number][cg_fragments][atom]['atom'] in backbone[cg_residues[residue_number][cg_fragments]['residue_name']]['b_connect']:
					BB_connect.append(str(atom))
		terminal_residue=False
		for frag_val,cg_fragments in enumerate(cg_residues[residue_number]):
			center=cg_residues[residue_number][cg_fragments]['coord']
			at_connections,cg_connections=[],[]
			if len(connect)>0:
				connect=np.array(connect)
				connect=connect[connect[:,2].argsort()]
				run=np.where(connect[:,0]==str(cg_fragments))
				count=0
				cg_con={}
				repeat=[]
				cg_temp=[]
				for val,con in enumerate(connect):
					repeat = len(np.where(connect[:,2]==con[2])[0])
					if con[2] in connect[run][:,2] and connect[val][0] != cg_fragments:
						cg_temp.append(cg_residues[residue_number][connect[val][0]]['coord']-center)
						if len(cg_temp) == repeat-1:
							cg_temp=np.array(cg_temp)
							cg_connections.append(np.mean(cg_temp, axis=0))
							cg_temp=[]
					if con[0] == cg_fragments and con[2] in connect[run][:,2]:
						at_connections.append(at_residues[residue_number][con[0]][con[1]]['coord']-center)
			if cg_fragments=='BB':
				try:
					cg_connections.append(cg_residues[residue_number-1]['BB']['coord']-center)
					at_connections.append(at_residues[residue_number]['BB'][BB_connect[0]]['coord']-center)
				except:
					res=residue_number
					pass
				try:
					cg_connections.append(cg_residues[residue_number+1]['BB']['coord']-center)
					at_connections.append(at_residues[residue_number]['BB'][BB_connect[1]]['coord']-center)
				except:
					pass
				try:
					xyz_prev=[cg_residues[residue_number-1]['BB']['coord'][0],cg_residues[residue_number-1]['BB']['coord'][1],cg_residues[residue_number-1]['BB']['coord'][2]]				
					xyz_cur=[cg_residues[residue_number]['BB']['coord'][0],cg_residues[residue_number]['BB']['coord'][1],cg_residues[residue_number]['BB']['coord'][2]]
					dist=np.sqrt(((xyz_prev[0]-xyz_cur[0])**2)+((xyz_prev[1]-xyz_cur[1])**2)+((xyz_prev[2]-xyz_cur[2])**2))
					if dist > 5:
						terminal_residue=True
						chain_count+=1
						backbone_coords[chain_count]=[]
						backbone_coords[chain_count].append(xyz_cur+[1])
						if args.v:
							print('\nchain number\tresidue no in pdb\tlength of chain')
							print('\n',chain_count,'\t',dist,'\t',residue_number,'\t',residue_number-res)
						res=residue_number-1
						pdb_output.write('TER\n')
						pdb_output = open(working_dir+'PROTEIN/PROTEIN_novo_'+str(chain_count)+'.pdb', 'w')
						pdb_output.write('REMARK    GENERATED BY sys_setup_script\nTITLE     SELF-ASSEMBLY-MAYBE\nREMARK    Good luck\n\
'+box_vec+'MODEL        1\n')
						at_counter=1
					else:
						backbone_coords[chain_count].append(xyz_cur+[1])
				except:
					xyz_cur=[cg_residues[residue_number]['BB']['coord'][0],cg_residues[residue_number]['BB']['coord'][1],cg_residues[residue_number]['BB']['coord'][2]]
					backbone_coords[chain_count].append(xyz_cur+[1])
					pass
			xyz_rot_apply=rotate(np.array(at_connections), np.array(cg_connections))
			for atom in at_residues[residue_number][cg_fragments]:
				at_residues[residue_number][cg_fragments][atom]['coord'] = rotate_atom(at_residues[residue_number][cg_fragments][atom]['coord'], center, xyz_rot_apply)
				final_at_residues[atom]=at_residues[residue_number][cg_fragments][atom]
		for at_val, atom in enumerate(final_at_residues):
			pdb_output.write(pdbline%((int(at_counter),final_at_residues[str(at_val+1)]['atom'],final_at_residues[str(at_val+1)]['res_type'],ascii_uppercase[chain_count],int(residue_number),\
final_at_residues[str(at_val+1)]['coord'][0],final_at_residues[str(at_val+1)]['coord'][1],final_at_residues[str(at_val+1)]['coord'][2],1,0))+'\n')
			at_counter+=1
		terminal_residue=False
	
	pdb_output.close()
	if args.a != None:
		read_in_atomistic(backbone_coords, chain_count+1)

	return [['PROTEIN', chain_count+1]]



def read_in_atomistic(cg_coords, chain_count):
	os.chdir(start_dir)
	if not os.path.exists(args.a):
		sys.exit('cannot find atomistic protein : '+args.a)
#### read in atomistic fragments into dictionary residue_list[0]=x,y,z,atom_name	
	atomistic_protein_input={}
	chain_count=0

	with open(args.a, 'r') as pdb_input:
		atomistic_protein_input[chain_count]={}
		for line_nr, line in enumerate(pdb_input.readlines()):
			if line.startswith('ATOM'):
				line_sep = pdbatom(line)
				if line_sep['atom_name'] in ['OT', 'O1', 'O2']:
					line_sep['atom_name']='O'
				if line_sep['atom_name'] == backbone[line_sep['residue_name']]['b_connect'][1]:
					C_ter=[line_sep['x'],line_sep['y'],line_sep['z']]
					C_resid=line_sep['residue_id']
					C=True
				try:
					if line_sep['atom_name'] == backbone[line_sep['residue_name']]['b_connect'][0]:
						N_resid=line_sep['residue_id']
						N_ter=[line_sep['x'],line_sep['y'],line_sep['z']]
						N=True
					dist=np.sqrt(((N_ter[0]-C_ter[0])**2)+((N_ter[1]-C_ter[1])**2)+((N_ter[2]-C_ter[2])**2))
					if N and C and C_resid != N_resid and dist > 3:
						N_ter, C_ter=False, False
						chain_count+=1
						atomistic_protein_input[chain_count]={}
				except:
					pass
#### get COM of fragment 
				atomistic_protein_input[chain_count][line_sep['atom_number']]={'coord':np.array([line_sep['x'],line_sep['y'],line_sep['z']]),'atom':line_sep['atom_name'], 'res_type':line_sep['residue_name'],'frag_mass':0, 'resid':line_sep['residue_id']}
				if 'H' not in line_sep['atom_name'] and line_sep['atom_name'] in backbone[line_sep['residue_name']]['atoms']:
					for atom in line_sep['atom_name']:
						if atom in mass:
							atomistic_protein_input[chain_count][line_sep['atom_number']]['frag_mass']=mass[atom]
	atomistic_protein_centered, cg_coords_mass = center_atomistic(atomistic_protein_input, cg_coords, chain_count+1)

	for chain in range(chain_count+1):
		pdb_output = open(working_dir+'PROTEIN/PROTEIN_at-input_'+str(chain)+'.pdb', 'w')
		pdb_output.write('REMARK    GENERATED BY sys_setup_script\nTITLE     SELF-ASSEMBLY-MAYBE\nREMARK    Good luck\n\
'+box_vec+'MODEL        1\n')
		at_centers, at_centers_iter=[],[]
		initial=True
		for atoms in atomistic_protein_input[chain]:

			if initial:
				resid_prev=atomistic_protein_input[chain][atoms]['resid']
				at_centers_iter.append(np.append(atomistic_protein_input[chain][atoms]['coord'],atomistic_protein_input[chain][atoms]['frag_mass']))
				initial=False				
			elif atomistic_protein_input[chain][atoms]['resid'] != resid_prev:
				at_centers.append(np.average(np.array(at_centers_iter)[:,:3], axis=0, weights=np.array(at_centers_iter)[:,3]))

				at_centers_iter=[]
				at_centers_iter.append(np.append(atomistic_protein_input[chain][atoms]['coord'],atomistic_protein_input[chain][atoms]['frag_mass']))
				resid_prev=atomistic_protein_input[chain][atoms]['resid']
			else:
				at_centers_iter.append(np.append(atomistic_protein_input[chain][atoms]['coord'],atomistic_protein_input[chain][atoms]['frag_mass']))
				resid_prev=atomistic_protein_input[chain][atoms]['resid']
		at_centers.append(np.average(np.array(at_centers_iter)[:,:3], axis=0, weights=np.array(at_centers_iter)[:,3]))

		if len(at_centers) != len(cg_coords[chain]):
			os.exit('In chain '+str(chain)+' the atommistic input does not match the CG. \n\
number of CG residues '+str(len(cg_coords[chain]))+'\nnumber of AT residues '+str(len(at_centers)))
		xyz_rot_apply = rotate(np.array(at_centers)-cg_coords_mass[chain], np.array(cg_coords[chain])[:,:3]-cg_coords_mass[chain])
		if args.v:
			print('\nThe proteins chains are rotated around the COM of all the backbone heavy atoms.\n')
			print('rotating chain ', chain, 'by ',np.round(np.degrees(xyz_rot_apply[0]),2),', ',np.round(np.degrees(xyz_rot_apply[1]),2),', ',np.round(np.degrees(xyz_rot_apply[2]),2))
		for atom in atomistic_protein_input[chain]:
			atomistic_protein_input[chain][atom]['coord'] = rotate_atom(atomistic_protein_input[chain][atom]['coord'], cg_coords_mass[chain], xyz_rot_apply)
		for at_val, atom in enumerate(atomistic_protein_input[chain]):
			pdb_output.write(pdbline%((int(atom),atomistic_protein_input[chain][atom]['atom'],atomistic_protein_input[chain][atom]['res_type'],ascii_uppercase[chain],atomistic_protein_input[chain][atom]['resid'],\
atomistic_protein_input[chain][atom]['coord'][0],atomistic_protein_input[chain][atom]['coord'][1],atomistic_protein_input[chain][atom]['coord'][2],1,0))+'\n')


def center_atomistic(atomistic_protein_input, cg_coords, chain_count): 
	cg_masses=[]
	for chain in range(chain_count):
		protein_mass=[]
		for at in atomistic_protein_input[chain]:
			protein_mass.append([atomistic_protein_input[chain][at]['coord'][0],atomistic_protein_input[chain][at]['coord'][1],atomistic_protein_input[chain][at]['coord'][2],atomistic_protein_input[chain][at]['frag_mass']])
		atomistic_protein_mass=np.average(np.array(protein_mass)[:,:3], axis=0, weights=np.array(protein_mass)[:,3])#### add center of mass of CG_proteins
		cg_masses.append(np.average(np.array(cg_coords[chain])[:,:3], axis=0, weights=np.array(cg_coords[chain])[:,3]))
		for at_id, atom in enumerate(atomistic_protein_input[chain]):
			atomistic_protein_input[chain][atom]['coord']=atomistic_protein_input[chain][atom]['coord']-(atomistic_protein_mass-cg_masses[chain])
	return atomistic_protein_input, cg_masses

def minimise_protein(chain_count):
	os.chdir(working_dir+'/PROTEIN')
	if not os.path.exists(working_dir+'FORCEFIELD'):
		os.mkdir(working_dir+'FORCEFIELD')
	copy_tree(forcefield_location+forcefield+'.ff', working_dir+'PROTEIN/'+forcefield+'.ff/.')
	if not os.path.exists('min'):
		os.mkdir('min')
	for chain in range(chain_count):
		with open('em_'+str(chain)+'.mdp','w') as em:
			em.write('define = -DPROTEIN_'+str(chain)+'_CA_posre.itp\nintegrator = steep\nnsteps = 10000\nemtol = 1000\nemstep = 0.001\ncutoff-scheme = Verlet\n')
		if args.a != None:
			minimise_protein_chain(chain, 'at-input_')
		minimise_protein_chain(chain, 'novo_')
	os.chdir('..')
	if args.a != None:
		merge_residues('PROTEIN', chain_count, '_at-input')
	merge_residues('PROTEIN', chain_count, '_novo')

def minimise_protein_chain(chain, input):
	gromacs(gmx+' pdb2gmx -f PROTEIN_'+input+str(chain)+'.pdb -o PROTEIN_'+input+str(chain)+'_gmx.pdb -ignh -water none \
-p PROTEIN_'+input+str(chain)+'.top  -i PROTEIN_'+str(chain)+'_posre.itp << EOF \n1\nEOF') ### single chains
	convert_topology('PROTEIN_'+input, chain)
	write_topol('PROTEIN_'+input, chain, 1, str(chain))
	gromacs(gmx+' grompp -f em_'+str(chain)+'.mdp -p PROTEIN_'+input+str(chain)+'.top -c PROTEIN_'+input+str(chain)+'_gmx.pdb -o min/PROTEIN_'+input+str(chain)+' -maxwarn 1')
	os.chdir('min')
	gromacs(gmx+' mdrun -v -nt 10 -deffnm PROTEIN_'+input+str(chain))
	gromacs(gmx+' editconf -f PROTEIN_'+input+str(chain)+'.gro -o PROTEIN_'+input+str(chain)+'.pdb -pbc')
	os.chdir('..')	

def simulation_setup(resid, residue_type, fragment, fragment_names, multiplier):
	resid=int(resid)
	if residue_type not in ['ION','SOL']:
		with open(residue_type+'/mapping_'+residue_type+'.ndx', 'w') as ndx_output:
			for line in fragment:
				ndx_output.write(line+'\n')
		os.chdir(working_dir+residue_type)
		write_topol(residue_type, 1, multiplier,'')
		make_min(residue_type, fragment_names)
		pool = mp.Pool(mp.cpu_count())
		pool_process = pool.map_async(gromacs, [(gmx+' grompp \
-po md_out-'+residue_type+'_'+str(rid)+' \
-f em_'+residue_type+'.mdp \
-p topol_'+residue_type+'.top \
-n mapping_'+residue_type+'.ndx \
-r '+residue_type+'_'+str(rid)+'.pdb \
-c '+residue_type+'_'+str(rid)+'.pdb \
-o min/'+residue_type+'_'+str(rid)+' -maxwarn 1') \
	for rid in range(0, resid+1)]).get()			## minimisation grompp parallised
		pool.close()
		os.chdir('min')
		pool = mp.Pool(mp.cpu_count())
		pool.map_async(gromacs, [(gmx+' mdrun -v -nt 1 -s '+residue_type+'_'+str(rid)+' -deffnm '+residue_type+'_'+str(rid)) \
for rid in range(0, resid+1)]).get()
		pool.close()
		pool = mp.Pool(mp.cpu_count())
		pool.map_async(gromacs, [(gmx+' editconf -f '+residue_type+'_'+str(rid)+'.gro \
-o '+residue_type+'_'+str(rid)+'.pdb -pbc') for rid in range(0, resid+1)]).get()
		pool.close()
		os.chdir('../..')
		merge_residues(residue_type, resid+1,'')
	elif residue_type == 'SOL':
		os.chdir(working_dir+residue_type)
		write_topol(residue_type, resid+1, multiplier,'')
		make_min(residue_type, [])
		gromacs(gmx+' grompp \
-po md_out-'+residue_type+' \
-f em_'+residue_type+'.mdp \
-p topol_'+residue_type+'.top \
-r mapping_'+residue_type+'.pdb \
-c mapping_'+residue_type+'.pdb \
-o min/'+residue_type)
		os.chdir('min')
		gromacs(gmx+' mdrun -v -nt 10 -deffnm '+residue_type)
		gromacs(gmx+' editconf -f '+residue_type+'.gro -o '+residue_type+'_merged.pdb -pbc')
		os.chdir('../..')

def rotate_atom(coord, center,xyz_rot_apply):
	coord =  coord-center
	coord =  coord.dot(eulerAnglesToRotationMatrix([xyz_rot_apply[0],0,0])) 
	coord =  coord.dot(eulerAnglesToRotationMatrix([0,xyz_rot_apply[1],0])) 
	coord =  coord.dot(eulerAnglesToRotationMatrix([0,0,xyz_rot_apply[2]])) 
	coord =  coord+center
	return coord



def rotate(at_connections, cg_connections):
	xyz_rot_apply=[]
	for xyz_rot in [x_rot,y_rot,z_rot]:
		min_dist=[]
		for rot_val, rotation in enumerate(xyz_rot):
			check = at_connections.dot(rotation)
			individual_connections=[]
			for connect in range(len(cg_connections)):
				individual_connections.append(np.sqrt(((check[connect][0]-cg_connections[connect][0])**2)+((check[connect][1]-cg_connections[connect][1])**2)+((check[connect][2]-cg_connections[connect][2])**2)))
			min_dist.append(individual_connections)
		inter=np.array([])
		for mdist in np.array(min_dist):
			inter = np.append(inter, np.sqrt(np.mean(mdist**2)))
		at_connections = at_connections.dot(xyz_rot[int(np.where(inter==np.min(inter))[0])])
		xyz_rot_apply.append(np.radians(int(np.where(inter==np.min(inter))[0])*5))
	return xyz_rot_apply


def write_topol(residue_type, residue_number, multiplier, chain):
	found=False
	with open('topol_'+residue_type+chain+'.top', 'w') as topol_write:
		topol_write.write('; Include forcefield parameters\n#include \"'+working_dir+'FORCEFIELD/charmm36-jul2017-updated.ff/forcefield.itp\"\n')
		topol_write.write('#include \"'+working_dir+'/FORCEFIELD/'+forcefield+'.ff/tip3p.itp\"\n\n#include \"'+working_dir+'/FORCEFIELD/'+forcefield+'.ff/ions.itp\"\n\n')
		if residue_type not in ['ION','SOL']:
			for directory in range(len(np_directories)):
				if os.path.exists(np_directories[directory][0]+residue_type+'/'+residue_type+'.itp'):
					topol_write.write('#include \"'+np_directories[directory][0]+residue_type+'/'+residue_type+'.itp\"\n')
					found=True
			if os.path.exists(working_dir+'/PROTEIN/'+residue_type+chain+'.itp'):
				topol_write.write('#include \"'+residue_type+chain+'.itp\"\n')
				found=True
			if not found:
				sys.exit('cannot find itp : '+residue_type+'/'+residue_type+chain)
		topol_write.write('\n\n[ system ]\n; Name\nSomething clever....\n\n[ molecules ]\n; Compound        #mols\n')
		topol_write.write(residue_type+chain+'    '+str(residue_number*multiplier))



def merge_residues(residue_type, resid, protein):
	merge=[]
	for res_val in range(0,resid):
		if os.path.exists(working_dir+residue_type+'/min/'+residue_type+protein+'_'+str(res_val)+'.pdb'):
			if residue_type == 'PROTEIN':
				posres_output = open(working_dir+'PROTEIN/PROTEIN_'+str(res_val)+'_CA_posre.itp', 'w')
				posres_output.write('[ position_restraints ]\n; atom  type      fx      fy      fz\n')
			with open(working_dir+residue_type+'/min/'+residue_type+protein+'_'+str(res_val)+'.pdb', 'r') as pdb_input:
				at_counter=0
				for line in pdb_input.readlines():
					if line.startswith('ATOM'):
						merge.append(line)
						line_sep = pdbatom(line)
						at_counter+=1
						if line_sep['atom_name'] in ['CA','CB']:
							posres_output.write(str(at_counter)+'     1  1000  1000  1000\n')
			if residue_type == 'PROTEIN':
				posres_output.close()
		else:
			sys.exit('cannot find minimised residue: \n'+residue_type+'/min/'+residue_type+protein+'_'+str(res_val)+'.pdb')	
	with open(working_dir+residue_type+'/min/'+residue_type+protein+'_merged.pdb', 'w') as pdb_output:
		pdb_output.write('REMARK    GENERATED BY sys_setup_script\nTITLE     SELF-ASSEMBLY-MAYBE\nREMARK    Good luck\n\
'+box_vec+'MODEL        1\n')
		for line in merge:
			pdb_output.write(line)

def convert_topology(topol, protein_number):
	if Path(topol+str(protein_number)+'.top').exists():
		read=False
		with open(topol+str(protein_number)+'.itp', 'w') as itp_write:
			for line in open(topol+str(protein_number)+'.top', 'r').readlines():
				if len(line.split()) > 1: 
					if read == False and line.split()[1] == 'moleculetype':
						read = True
					if read == True and line.split()[1] == 'system':
						read = False
					if line.split()[0] == 'Protein' or line.split()[0][:-1] == 'Protein_chain_':
						line= re.sub('Protein_chain_.?', 'Protein',line)
						line= re.sub('Protein', 'protein_'+str(protein_number),line)
				if read ==True:
					itp_write.write(line)
			itp_write.write('; Include CA Position restraint file\n#ifdef POSRES_CA\n#include \"PROTEIN_'+str(protein_number)+'_CA_posre.itp\"\n#endif')
	else:
		sys.exit('cannot find : '+topol+'_'+str(protein_number)+'.top')


def merge_system_pdbs(system, box_vec, protein):
	print('Merging atomistic files')
	system_input = copy.deepcopy(system)
	for res_val, residue_type in enumerate(system):
		system_input[res_val][1]=[]
		if residue_type[0]=='SOL':
			sol_index=res_val
	for res_val, residue_type in enumerate(system):
		if residue_type[0] == 'PROTEIN':
			prot_input=protein
		else:
			prot_input=''
		if os.path.exists(working_dir+residue_type[0]+'/min/'+residue_type[0]+prot_input+'_merged.pdb'):
			with open(working_dir+residue_type[0]+'/min/'+residue_type[0]+prot_input+'_merged.pdb', 'r') as pdb_input:
				for line in pdb_input.readlines():
					if line.startswith('ATOM'):
						line_sep = pdbatom(line)
						if line_sep['residue_name']=='SOL':
							system_input[sol_index][1].append(line)
						else:
							system_input[res_val][1].append(line)
		else:
			sys.exit('cannot find minimised residue: \n'+ working_dir+residue_type[0]+'/min/'+residue_type[0]+prot_input+'_merged.pdb')	
	with open(working_dir+'MERGED/merged_cg2at'+protein+'.pdb', 'w') as pdb_output:
		pdb_output.write('REMARK    GENERATED BY sys_setup_script\nTITLE     SELF-ASSEMBLY-MAYBE\nREMARK    Good luck\n\
'+box_vec+'MODEL        1\n')
		for residue_type in system_input:
			for line in residue_type[1]:
				pdb_output.write(line)
		pdb_output.write('TER\nENDMDL')

def write_merged_topol(system, box_vec, protein):
	for topology in ['topol_final', 'topol_alchembed']:
		with open(topology+'.top', 'w') as topol_write:
			topol_write.write('; Include forcefield parameters\n#include \"'+working_dir+'FORCEFIELD/charmm36-jul2017-updated.ff/forcefield.itp\"\n')
			topol_write.write('#include \"'+working_dir+'/FORCEFIELD/'+forcefield+'.ff/tip3p.itp\"\n\n#include \"'+working_dir+'/FORCEFIELD/'+forcefield+'.ff/ions.itp\"\n\n')
			for residue_type in system:
				found=False
				if residue_type[0] not in ['ION','SOL']:
					for directory in range(len(np_directories)):
						if os.path.exists(np_directories[directory][0]+residue_type[0]+'/'+residue_type[0]+'.itp'):		
							topol_write.write('#include \"'+residue_type[0]+'.itp\"\n')
							copyfile(np_directories[directory][0]+residue_type[0]+'/'+residue_type[0]+'.itp', residue_type[0]+'.itp')
							found=True
					if residue_type[0] == 'PROTEIN':
						if topology == 'topol_alchembed':
							topol_write.write('#include \"PROTEIN.itp\"\n')
						else:
							for protein_unit in range(residue_type[1]): 
								topol_write.write('#include \"PROTEIN_'+str(protein_unit)+'.itp\"\n')
								copyfile('../PROTEIN/PROTEIN'+protein+'_'+str(protein_unit)+'.itp', 'PROTEIN_'+str(protein_unit)+'.itp')
								copyfile('../PROTEIN/PROTEIN_'+str(protein_unit)+'_CA_posre.itp', 'PROTEIN_'+str(protein_unit)+'_CA_posre.itp')
								copyfile('../PROTEIN/PROTEIN_'+str(protein_unit)+'_posre.itp', 'PROTEIN_'+str(protein_unit)+'_posre.itp')
			topol_write.write('[ system ]\n; Name\nSomething clever....\n\n[ molecules ]\n; Compound        #mols\n')


			for residue_type in system:
				if residue_type[0] not in  ['ION', 'PROTEIN']:
					topol_write.write(residue_type[0]+'    '+str(residue_type[1])+'\n')	
				elif residue_type[0] == 'ION':
					for ion_type in residue_type[1]:
						 topol_write.write(ion_type[0]+'    '+str(ion_type[1])+'\n')	
				elif residue_type[0] == 'PROTEIN':
					if topology == 'topol_alchembed':
						topol_write.write('PROTEIN    1\n')	
					else:
						for protein_unit in range(residue_type[1]):
							 topol_write.write('PROTEIN_'+str(protein_unit)+'    1\n')	

def minimise_merged_pdbs(system, box_vec, protein):
	print('Minising merged atomistic files')
	write_merged_topol(system, box_vec, protein)
	make_min('merged_cg2at', [])
	gromacs(gmx+' grompp -po md_out-merged_cg2at -f em_merged_cg2at.mdp -p topol_final.top \
-r merged_cg2at'+protein+'.pdb -c merged_cg2at'+protein+'.pdb -o min/merged_cg2at'+protein+'_minimised')
	os.chdir('min')
	gromacs(gmx+' mdrun -v -nt 10 -deffnm merged_cg2at'+protein+'_minimised')
	gromacs(gmx+' editconf -f merged_cg2at'+protein+'_minimised.gro -o merged_cg2at'+protein+'_minimised.pdb -pbc')

def alchembed(system, box_vec, protein):
	print('Running alchembed')
	os.chdir(working_dir+'MERGED')
	if not os.path.exists('alchembed'):
		os.mkdir('alchembed')
	for chain in range(system[0][1]):
		with open('alchembed_'+str(chain)+'.mdp', 'w') as alchembed:
			alchembed.write('define = -DPOSRES\nintegrator = sd\nnsteps = 1000\ndt = 0.0005\ncontinuation = no\nconstraint_algorithm = lincs\nconstraints	= h-bonds\nns_type = grid\nnstlist = 25\n\
rlist = 1\nrcoulomb	= 1\nrvdw = 1\ncoulombtype	= PME\npme_order = 4\nfourierspacing = 0.16\ntcoupl	= V-rescale\ntc-grps = system\ntau_t = 0.1\nref_t = 310\npcoupl	= no\n\
pbc = xyz\nDispCorr	= no\ngen_vel = yes\ngen_temp = 310\ngen_seed = -1\nfree_energy = yes\ninit_lambda = 0.00\ndelta_lambda = 1e-3\nsc-alpha = 0.1000\nsc-power = 1\nsc-r-power = 6\n\
couple-moltype = protein_'+str(chain)+'\ncouple-lambda0 = none\ncouple-lambda1 = vdw')
		if chain == 0:
			gromacs(gmx+' grompp -po md_out-merged_cg2at_alchembed_'+str(chain)+' -f alchembed_'+str(chain)+'.mdp -p topol_final.top -r min/merged_cg2at_at-input_minimised.gro \
-c min/merged_cg2at_at-input_minimised.gro -o alchembed/merged_cg2at_at-input_alchembed_'+str(chain)+' -maxwarn 1')
		else:
			gromacs(gmx+' grompp -po md_out-merged_cg2at_alchembed_'+str(chain)+' -f alchembed_'+str(chain)+'.mdp -p topol_final.top -r min/merged_cg2at_at-input_minimised.gro \
-c alchembed/merged_cg2at_at-input_alchembed_'+str(chain-1)+'.gro -o alchembed/merged_cg2at_at-input_alchembed_'+str(chain)+' -maxwarn 1')			

		os.chdir('alchembed')
		gromacs(gmx+' mdrun -v -nt 10 -deffnm merged_cg2at_at-input_alchembed_'+str(chain))
		os.chdir('..')
	gromacs(gmx+' editconf -f alchembed/merged_cg2at_at-input_alchembed_'+str(chain)+'.gro -o alchembed/merged_cg2at_at-input_alchembed.pdb -pbc')
	copyfile('alchembed/merged_cg2at_at-input_alchembed.pdb', final_dir+'merged_cg2at_at-input_alchembed.pdb')

		

def steered_md_atomistic_to_cg_coord():
	print('Creating tpr file for steered MD on atoms CA and CB to fit CG backbone')

	os.chdir(working_dir+'MERGED')
	if not os.path.exists('steered_md'):
		os.mkdir('steered_md')
	with open('steered_md.mdp', 'w') as steered_md:
		steered_md.write('define = -DPOSRES_CA\nintegrator = md\nnsteps = 10000\ndt = 0.001\ncontinuation	= no\nconstraint_algorithm = lincs\nconstraints	= h-bonds\nns_type = grid\nnstlist = 25\n\
rlist = 1\nrcoulomb	= 1\nrvdw = 1\ncoulombtype	= PME\npme_order = 4\nfourierspacing = 0.16\ntcoupl	= V-rescale\ntc-grps = system\ntau_t = 0.1\nref_t = 310\npcoupl	= no\n\
pbc = xyz\nDispCorr	= no\ngen_vel = yes\ngen_temp = 310\ngen_seed = -1')	
	gromacs(gmx+' grompp -po md_out-merged_cg2at_pull -f steered_md.mdp -p topol_final.top \
-r min/merged_cg2at_novo_minimised.pdb -c alchembed/merged_cg2at_at-input_alchembed.pdb -o steered_md/merged_cg2at_steered_md -maxwarn 1')
	if args.steer:
		print('Running steered MD')
		os.chdir('steered_md')
		gromacs(gmx+' mdrun -v -nt 10 -deffnm merged_cg2at_steered_md')
		gromacs(gmx+' editconf -f merged_cg2at_steered_md.gro -o merged_cg2at_steered_md.pdb -pbc')
		copyfile('merged_cg2at_steered_md.pdb', final_dir+'final_cg2at_CA_CB_steered_md.pdb')



def eulerAnglesToRotationMatrix(theta) :
     
    R_x = np.array([[1,         0,                  0                   ],
                    [0,         math.cos(theta[0]), -math.sin(theta[0]) ],
                    [0,         math.sin(theta[0]), math.cos(theta[0])  ]
                    ])
         
    R_y = np.array([[math.cos(theta[1]),    0,      math.sin(theta[1])  ],
                    [0,                     1,      0                   ],
                    [-math.sin(theta[1]),   0,      math.cos(theta[1])  ]
                    ])
                 
    R_z = np.array([[math.cos(theta[2]),    -math.sin(theta[2]),    0],
                    [math.sin(theta[2]),    math.cos(theta[2]),     0],
                    [0,                     0,                      1]
                    ])
                                        
    R = np.dot(R_z, np.dot( R_y, R_x ))


    return R

start_time=time.time()

print('\nInitialisation of CG2AT v2')

### finds initial rotation matrices
x_rot, y_rot, z_rot=[],[],[]
for angle in range(0,360, 5):
	angle=np.radians(angle)
	x_rot.append(eulerAnglesToRotationMatrix([angle,0,0]))
	y_rot.append(eulerAnglesToRotationMatrix([0,angle,0]))
	z_rot.append(eulerAnglesToRotationMatrix([0,0,angle]))

np_system,p_system=[],[]

water_mol_number=4
gmx='gmx'
pdbline = "ATOM  %5d %4s %4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f"
mass = {'H': 1,'C': 12,'N': 14,'O': 16,'S': 32,'P': 31,'M': 0, 'B': 32}

### sets up file locations
start_dir=os.getcwd()
working_dir=os.getcwd()+'/CG2AT/'
final_dir=os.getcwd()+'/CG2AT/FINAL/'
if not os.path.exists(working_dir):
	os.mkdir(working_dir)
	os.mkdir(final_dir)
script_dir=os.path.dirname(os.path.realpath(__file__))+'/'
initialisation_time=time.time()

### read in and sort forcefield info
forcefield_available_prov, forcefield_available_user = read_forcefield()
f_number = forcefield_selection(forcefield_available_prov, forcefield_available_user)
forcefield_location, forcefield=sort_forcefield(forcefield_available_prov, forcefield_available_user, f_number)

### reads in and sorts fragment information
np_residues, p_residues,np_directories, p_directories, water=fetch_residues()
backbone=fetch_fragment(p_directories)


print('This script is now hopefully doing the following:\n')

#### read in CG file
print('Reading in your CG representation')
cg_residues, box_vec = read_initial_pdb()
os.chdir(working_dir)
read_in_time=time.time()


#### converts protein into atomistic and minimises
print('Converting CG into a atomistic representation')
if 'protein' in cg_residues:
	if len(cg_residues['protein'])>0:
		p_system=build_protein_atomistic_system(cg_residues['protein'], box_vec)
		minimise_protein(p_system[0][1])
sort_protein_time=time.time()

#### converts non protein residues into atomistic and minimises 
if len([key for value, key in enumerate(cg_residues) if key not in ['protein']]) > 0:
	np_system=build_atomistic_system(cg_residues, box_vec)
	build_non_protein_time=time.time()
mid_point_time=time.time()


#### creates merged folder
print('Merging all residue types to single files')
system=p_system+np_system
if len(system)>0:
	if not os.path.exists(working_dir+'MERGED'):
		os.mkdir(working_dir+'MERGED')
	os.chdir(working_dir+'MERGED') 	
#### copies all itp files and topologies from whereever they are stored
	for file_name in os.listdir('.'):
		if file_name.endswith('.itp') or file_name.endswith('.top'):
			copyfile(file_name, final_dir+file_name)
#### merges provided atomistic protein and residues types into a single pdb file into merged directory
	if args.a != None:
		merge_system_pdbs(system, box_vec,'_at-input' )
		minimise_merged_pdbs(system, box_vec, '_at-input')
		alchembed(p_system, box_vec, '_at-input')
#### merges de novo protein and residues types into a single pdb file into merged directory
	merge_system_pdbs(system, box_vec,'_novo' )
	minimise_merged_pdbs(system, box_vec, '_novo')
	copyfile('merged_cg2at_novo_minimised.pdb', final_dir+'final_cg2at_novo_minimised.pdb')

#### runs steered MD on atomistic structure on CA and CB atoms
merge_time=time.time()
if args.a != None:
	steered_md_atomistic_to_cg_coord()
pull_time=time.time()

#### prints out system information
print('\nmolecules\tnumber')
print('---------\t------')
for section in system:
	if section[0] == 'ION':
		for solvent in section[1]:
			print(solvent[0],'\t', solvent[1])
	else:
		print(section[0],'\t', section[1])

#### prints out script timings for each section
print('\nInitialisation: ', str(datetime.timedelta(minutes=np.round(initialisation_time-start_time, 2))).rsplit(':', 1)[0], ' min',\
'\nRead in CG system: ', str(datetime.timedelta(minutes=np.round(read_in_time-initialisation_time, 2))).rsplit(':', 1)[0], ' min' ,\
'\nBuild atomistic protein system: ', str(datetime.timedelta(minutes=np.round(sort_protein_time-read_in_time, 2))).rsplit(':', 1)[0], ' min',\
'\nBuild non protein system: ', str(datetime.timedelta(minutes=np.round(mid_point_time-sort_protein_time, 2))).rsplit(':', 1)[0], ' min', \
'\nMerge protein and non protein system: ', str(datetime.timedelta(minutes=np.round(merge_time-mid_point_time, 2))).rsplit(':', 1)[0], ' min', \
'\nSteered MD atomistic input to match CG: ', str(datetime.timedelta(minutes=np.round(pull_time-merge_time, 2))).rsplit(':', 1)[0], ' min', \
'\nTotal run time: ', str(datetime.timedelta(minutes=np.round(pull_time-start_time, 2))).rsplit(':', 1)[0], ' min')