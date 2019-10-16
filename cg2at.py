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
import random

parser = argparse.ArgumentParser(description='Converts CG representation into atomistic', epilog='Enjoy the program and best of luck!', allow_abbrev=True)
parser.add_argument('-c', help='coarse grain coordinates',metavar='pdb/gro',type=str, required=True)
parser.add_argument('-a', help='atomistic coordinates',metavar='pdb/gro',type=str)
parser.add_argument('-l', help='additional fragment library location',metavar='fragments folder',type=str)
parser.add_argument('-v', action="count", default=0, help="increase output verbosity (eg -vv, 3 levels)")
parser.add_argument('-f', help='additional forcefield location',metavar='forcefield folder', type=str)
parser.add_argument('-steer', help='do not run steered MD on atomistic structure (requires atomistic structure)', action='store_false')
args = parser.parse_args()
options = vars(args)

def mkdir_directory(directory):
	if not os.path.exists(directory):
		os.mkdir(directory)

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
	print('\n\n{0:^45}\n'.format('provided forcefields'))
	print('{0:^20}{1:^30}'.format('Selection','Forcefield'))
	print('{0:^20}{1:^30}'.format('---------','----------'))

	for force_num_prov, line in enumerate(forcefield_available_prov):
		# print('\t',force_num_prov,'\t\t',line[:-3])
		print('{0:^20}{1:^30}'.format(force_num_prov,line[:-3]))
	if args.f != None:
		print('\n{0:^45}\n'.format('User provided forcefields'))
		for force_num_user, line in enumerate(forcefield_available_user):
			print('{0:^20}{1:^30}'.format(force_num_user+force_num_prov,line[:-3]))	
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
		copy_tree(args.f+'/'+forcefield_available_user[f_number-len(forcefield_available_prov)], working_dir+'FORCEFIELD/'+forcefield_available_user[f_number-len(forcefield_available_prov)]+'/.')
		copy_tree(args.f+'/'+forcefield_available_user[f_number-len(forcefield_available_prov)], final_dir+forcefield_available_user[f_number-len(forcefield_available_prov)]+'/.')		
		return args.f+'/', forcefield_available_user[f_number-len(forcefield_available_prov)][:-3]

def fetch_residues():
#### list of directories and water types  [[root, folders...],[root, folders...]]
	np_directories, p_directories,mod_directories=[[]], [[]],[[]]
	water=[]
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
				if os.path.exists(script_dir+'database/'+forcefield+directory_type+'MOD/'):
					p_directories[-1].remove('MOD')
					p_directories.append([])
					for root, dirs, files in os.walk(script_dir+'database/'+forcefield+directory_type+'MOD/'):
						p_directories[-1].append(root)
						p_directories[-1]+=dirs
						mod_directories[-1].append(root)
						mod_directories[-1]+=dirs
						break
			break
#### adds non protein residues locations to np_directories from user defined
		if args.l != None:
			if os.path.exists(args.l+'/'+forcefield):
				for root, dirs, files in os.walk(args.l+'/'+forcefield+directory_type):
					if directory_type =='/non_protein/':
						np_directories.append([])
						np_directories[-1].append(root)
						np_directories[-1]+=dirs
#### adds water residues to water list
						if os.path.exists(args.l+'/'+forcefield+directory_type+'SOL'):
							for root, dirs, files in os.walk(script_dir+'/database/'+forcefield+directory_type+'SOL'):
								for water_type in files:
									water.append(water_type[:-4])
								break
#### adds protein residues locations to p_directories from user defined
					else:
						p_directories.append([])
						p_directories[-1].append(root)
						p_directories[-1]+=dirs	
						if os.path.exists(args.l+'/'+forcefield+directory_type+'MOD'):
							p_directories[-1].remove('MOD')
							p_directories.append([])
							mod_directories.append([])
							for root, dirs, files in os.walk(args.l+'database/'+forcefield+directory_type+'MOD'):
								p_directories[-1].append(root)
								p_directories[-1]+=dirs
								mod_directories[-1].append(root)
								mod_directories[-1]+=dirs
								break
					break
			else:
				sys.exit('Cannot find :'+args.l+'/database/forcefield')
#### sorts directories alphabetically and creates residue database
	p_residues, np_residues, mod_residues = [],[],[]
	for directory in range(len(mod_directories)):
		mod_directories[directory].sort()
		mod_residues+=mod_directories[directory][1:]
	for directory in range(len(p_directories)):
		p_directories[directory].sort()
		p_residues+=p_directories[directory][1:]
	for directory in range(len(np_directories)):
		np_directories[directory].sort()
		np_residues+=np_directories[directory][1:]
#### if verbose prints all fragments found
	if args.v >= 1:
		print('\n---------------------------------------------->  Verbose level 1 start')
		for directory in range(len(np_directories)):
			print('\nnon protein residues fragment directories found: \n\nroot file system\n')
			print(np_directories[directory][0],'\n\nresidues\n\n',np_directories[directory][1:], '\n')
		for directory in range(len(p_directories)):
			print('\nprotein residues fragment directories found: \n\nroot file system\n')
			print(p_directories[directory][0],'\n\nresidues\n\n',p_directories[directory][1:], '\n')
		print('\n---------------------------------------------->  Verbose level 1 end\n')

	return np_residues, p_residues, mod_residues, np_directories, p_directories, mod_directories, water

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
	if args.v >= 2:
		print('\n---------------------------------------------->  Verbose level 2 start')
		print('backbone atoms for each residue and connecting atoms:\n')
		for residue in backbone:
			print(residue, '\tbackbone atoms:', backbone[residue]['atoms'], '\n\tbackbone connecting atoms:', backbone[residue]['b_connect'],'\n')
		print('\n---------------------------------------------->  Verbose level 2 end\n')

	return backbone

def make_min(residue, fragments):
#### makes minimisation folder
	mkdir_directory('min')
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
	if args.v >= 3:
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

def rotate_atom(coord, center,xyz_rot_apply):
	coord =  coord-center  #### centers COM coordinates to 0,0,0
	coord =  coord.dot(eulerAnglesToRotationMatrix([xyz_rot_apply[0],0,0]))  #### rotates coord around x
	coord =  coord.dot(eulerAnglesToRotationMatrix([0,xyz_rot_apply[1],0]))  #### rotates coord around y
	coord =  coord.dot(eulerAnglesToRotationMatrix([0,0,xyz_rot_apply[2]]))  #### rotates coord around z
	coord =  coord+center #### translates coord back by original offset
	return coord



def rotate(at_connections, cg_connections):
	xyz_rot_apply=[]
#### iterates through rotation matrices
	for xyz_rot in [x_rot,y_rot,z_rot]:
		min_dist=[]
	#### iterates through rotation matrices 
		for rot_val, rotation in enumerate(xyz_rot):
		#### applies matrix to coordinates saved as check
			check = at_connections.dot(rotation)
		#### for each connection the distance is calculated and added to list
			individual_connections=[]
			for connect in range(len(cg_connections)):
				individual_connections.append(np.sqrt(((check[connect][0]-cg_connections[connect][0])**2)+((check[connect][1]-cg_connections[connect][1])**2)+((check[connect][2]-cg_connections[connect][2])**2)))
		#### for each rotation the connection distances are added to min_dist list 
			min_dist.append(individual_connections)
		inter=np.array([])
	#### the RMS is calculated for each rotation	
		for mdist in np.array(min_dist):
			inter = np.append(inter, np.sqrt(np.mean(mdist**2)))
	#### the rotation with the lowest RMS applied to the at_connections
		at_connections = at_connections.dot(xyz_rot[int(np.where(inter==np.min(inter))[0])])
	#### the optimal rotation is added to xyz_rot_apply list as radians
		xyz_rot_apply.append(np.radians(int(np.where(inter==np.min(inter))[0])*5))
	return xyz_rot_apply

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

def fragment_location(residue, fragment):  
#### runs through dirctories looking for the atomistic fragments returns the correct location
	if residue in p_residues:
		for directory in range(len(p_directories)):
			if os.path.exists(p_directories[directory][0]+residue+'/'+fragment):
				return p_directories[directory][0]+residue+'/'+fragment
		for directory in range(len(mod_directories)):
			if os.path.exists(mod_directories[directory][0]+residue+'/'+fragment):
				return mod_directories[directory][0]+residue+'/'+fragment
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
		at_connections.append(at_residues[cg_bead][str(at)]['coord']-center)
	return at_connections, cg_connections, center

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

def create_pdb(file_name):
	pdb_output = open(file_name, 'w')
	pdb_output.write('REMARK    GENERATED BY sys_setup_script\nTITLE     SELF-ASSEMBLY-MAYBE\nREMARK    Good luck\n\
'+box_vec+'MODEL        1\n')
	return pdb_output

############################################################ Read in CG file Section ################################################################

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
					if 'PROTEIN' not in cg_residues:  ## if protein does not exist add to dict
						cg_residues['PROTEIN']={}
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
						cg_residues['PROTEIN'][count]={} ### then create sub dictionary cg_residues['PROTEIN'][count]
						cg_residues['PROTEIN'][count]=residue_list ### adds residue list to dictionary key cg_residues['PROTEIN'][count]
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
		if count not in cg_residues['PROTEIN']:
			cg_residues['PROTEIN'][count]={}
		cg_residues['PROTEIN'][count]=residue_list
	else:
		if count not in cg_residues[line_sep['residue_name']]:
			cg_residues[line_sep['residue_name']][count]={}
		cg_residues[line_sep['residue_name']][count]=residue_list
#### checks if box vectors exist
	if 'box_vec' not in locals():### stops script if it cannot find box vectors
		sys.exit('missing box vectors')

	return cg_residues, box_vec


############################################################ Build Non Protein Section ################################################################

def build_atomistic_system(cg_residues, box_vec):
	system={}
	ion_sep={}
	atomistic_fragments={}
	fragment, fragment_names={},{}
#### for each residue type covert to atomistic except protein
	for residue_type in [key for value, key in enumerate(cg_residues) if key not in ['PROTEIN']]:
		fragment_names[residue_type], fragment[residue_type]=[], []
	#### reset counters for each residue type
		print('Converting residue type: ' +residue_type)
	#### creates folder for residue type
		mkdir_directory(working_dir+residue_type)
	#### fetches atoms for the residue type and centers then on cg bead
		atomistic_fragments[residue_type] = atomistic_non_protein(residue_type, cg_residues[residue_type])
	#### if the residue type is in ['SOL', 'ION'] a single pdb is created
		if residue_type in ['SOL', 'ION']:
			water_count=0     #### resets water counter to 0 
			mkdir_directory(working_dir+'SOL')	
		#### creates solvent directory and SOL key in system dictionay otherwise it appends solvent molecules to sol pdb
			if not os.path.exists(working_dir+'SOL'+'/SOL_0.pdb'):   
				system['SOL']=0
				pdb_sol = open(working_dir+'SOL'+'/SOL_0.pdb', 'w')
				pdb_sol.write('REMARK    GENERATED BY sys_setup_script\nTITLE     SELF-ASSEMBLY-MAYBE\nREMARK    Good luck\n\
'+box_vec+'MODEL        1\n')
			else:
				pdb_sol = open(working_dir+'SOL'+'/SOL_0.pdb', 'a')
		#### creates ion pdb with header
			if residue_type == 'ION':
			#### make minimisation directory and makes homogeneous diectory structure
				mkdir_directory(working_dir+'ION/min')
				pdb_ion = create_pdb(working_dir+residue_type+'/min/'+residue_type+'_merged.pdb')
	#### loop through all resids of that residue type 
		for resid in atomistic_fragments[residue_type]:
		#### if the residue type is not in ['SOL', 'ION'] a individual pdbs are created for each resid
			if residue_type not in ['SOL', 'ION']:
				pdb_output = create_pdb(working_dir+residue_type+'/'+residue_type+'_'+str(resid)+'.pdb')
		#### for every fragment in that resid
			for at_id, atom in enumerate(atomistic_fragments[residue_type][resid]):
			#### creates fragment index of the residue for the removal of COM translation of the fragment
				if resid==0 and residue_type not in ['ION', 'SOL']:
					if at_id==0 and 'previous_bead' not in locals():
						fragment[residue_type].append('[ '+str(atomistic_fragments[residue_type][resid][at_id+1]['cg_bead'])+' ]')
						fragment_names[residue_type].append(str(atomistic_fragments[residue_type][resid][at_id+1]['cg_bead']))
					if 'previous_bead' in locals():
						if previous_bead != atomistic_fragments[residue_type][resid][at_id+1]['cg_bead']:
							fragment[residue_type].append('[ '+str(atomistic_fragments[residue_type][resid][at_id+1]['cg_bead'])+' ]')
							fragment_names[residue_type].append(str(atomistic_fragments[residue_type][resid][at_id+1]['cg_bead']))
					fragment[residue_type].append(str(at_id+1)+' ')
					previous_bead=atomistic_fragments[residue_type][resid][atom]['cg_bead']
			#### separates out the water molecules from the ion in the fragment
				if residue_type in ['ION', 'SOL']:
				#### if not restype SOL it creates iontype variable and writes ion to pdb
					if atomistic_fragments[residue_type][resid][at_id+1]['res_type'] != 'SOL':
						if atomistic_fragments[residue_type][resid][at_id+1]['res_type'] not in system:
							system[atomistic_fragments[residue_type][resid][at_id+1]['res_type']]=1
						else:
							system[atomistic_fragments[residue_type][resid][at_id+1]['res_type']]+=1
						pdb_ion.write(pdbline%((int(at_id+1),atomistic_fragments[residue_type][resid][at_id+1]['atom'],atomistic_fragments[residue_type][resid][at_id+1]['res_type'],' ',1,\
					atomistic_fragments[residue_type][resid][at_id+1]['coord'][0],atomistic_fragments[residue_type][resid][at_id+1]['coord'][1],atomistic_fragments[residue_type][resid][at_id+1]['coord'][2],1,0))+'\n')
					else:
					#### if restype SOL and has a mass over 1 (eg H) adds a count to the water_count also writes to solvent pdb
						if atomistic_fragments[residue_type][resid][at_id+1]['frag_mass'] > 1:					
							water_count+=1
						pdb_sol.write(pdbline%((int(at_id+1),atomistic_fragments[residue_type][resid][at_id+1]['atom'],atomistic_fragments[residue_type][resid][at_id+1]['res_type'],' ',1,\
					atomistic_fragments[residue_type][resid][at_id+1]['coord'][0],atomistic_fragments[residue_type][resid][at_id+1]['coord'][1],atomistic_fragments[residue_type][resid][at_id+1]['coord'][2],1,0))+'\n')
			#### if residue_type not in ['ION', 'SOL'] write out to separate pdb
				else:
					pdb_output.write(pdbline%((int(atom),atomistic_fragments[residue_type][resid][at_id+1]['atom'],atomistic_fragments[residue_type][resid][at_id+1]['res_type'],' ',1,\
					atomistic_fragments[residue_type][resid][at_id+1]['coord'][0],atomistic_fragments[residue_type][resid][at_id+1]['coord'][1],atomistic_fragments[residue_type][resid][at_id+1]['coord'][2],1,0))+'\n')
	#### if restype ION add ion to system dictionary and add solvent molecules
		if residue_type == 'ION':
		#### creates a index group for each residue type containing lists of fragment
			if not os.path.exists(working_dir+'SOL/mapping_SOL.ndx'):
				ndx_output = open(working_dir+'SOL/mapping_SOL.ndx', 'w')
				ndx_output.close()
			system['SOL']+=water_count
		else:
		#### if restype is solvent updates the system dictionary
			if residue_type == 'SOL':
				system['SOL']+=water_count
			else:
			#### adds retype to system dictionary
				system[residue_type]=int(resid)+1
		#### creates a index group for each residue type containing lists of fragment
			if not os.path.exists(working_dir+residue_type+'/mapping_'+residue_type+'.ndx'):
				ndx_output = open(working_dir+residue_type+'/mapping_'+residue_type+'.ndx', 'w')
		#### writes fragments index for COM removal
			for line in fragment[residue_type]:
				ndx_output.write(line+'\n')
	return system, fragment_names

def atomistic_non_protein(cg_residue_type,cg_residues):
	atomistic_fragments={}  #### residue dictionary
#### run through every residue in a particular residue type

	for cg_resid, cg_residue in enumerate(cg_residues):
	#### get atomistic fragments for each bead and connectivity with other beads
		at_residues, connect = get_atomistic_fragments(cg_residue_type,cg_residues[cg_residue], cg_resid)	
		atomistic_fragments[cg_resid]={}  #### creats key in atomistic_fragments for each residue eg atomistic_fragments[1]
	#### runs through all beads in each resid
		for bead_number, cg_bead in enumerate(cg_residues[cg_residue]):
		#### if cg_residue_type not in ['SOL', 'ION'] as they have no connectivity 
			if cg_residue_type not in ['SOL', 'ION']:
			#### finds all beads that the cg_bead is connected to, and returns atom coord, cg coord and center of cg_bead
				at_connections,cg_connections, center=connectivity(bead_number, cg_bead, connect, at_residues, cg_residues,cg_residue)
			#### rotates at_connections to finds minimum RMS distance with cg_connections 
				xyz_rot_apply=rotate(np.array(at_connections), np.array(cg_connections))
		#### if ION/SOL a random rotation is applied to the cluster 
			else:
				center=cg_residues[cg_residue][cg_bead]['coord']
				xyz_rot_apply=[random.uniform(0, math.pi*2), random.uniform(0, math.pi*2), random.uniform(0, math.pi*2)]
			#### applies optimum rotation to each atom in the fragment 
			for atom in at_residues[cg_bead]:
				at_residues[cg_bead][atom]['coord']=rotate_atom(at_residues[cg_bead][atom]['coord'], center, xyz_rot_apply)  ### applies rotation
			#### adds bead number to at_residues dictionary
				at_residues[cg_bead][atom].update({'cg_bead':bead_number+1})
			#### adds atomistic fragment to new dictionary atomistic_fragments[resid][atom number] allows reordering of atoms by atom number in fragment database
				atomistic_fragments[cg_resid][int(atom)]=at_residues[cg_bead][atom]
	return atomistic_fragments

def non_protein_minimise(resid, residue_type, fragment_names):
#### in the case of SOL all residues are minimised, whilst in all other cases individual residues are minimised separately
	if residue_type != 'SOL':
		individual = 1
		resid=resid
	else:
		individual=resid
		resid=1
	os.chdir(working_dir+residue_type)
### write topology and minimisation parts (min folder and em.mdp)
	write_topol(residue_type, individual, '')
	make_min(residue_type, fragment_names)

#### spin up multiprocessing for grompp 
	pool = mp.Pool(mp.cpu_count())
	pool_process = pool.map_async(gromacs, [(gmx+' grompp \
-po md_out-'+residue_type+'_'+str(rid)+' \
-f em_'+residue_type+'.mdp \
-p topol_'+residue_type+'.top \
-n mapping_'+residue_type+'.ndx \
-r '+residue_type+'_'+str(rid)+'.pdb \
-c '+residue_type+'_'+str(rid)+'.pdb \
-o min/'+residue_type+'_'+str(rid)+' -maxwarn 1') \
for rid in range(0, resid)]).get()			## minimisation grompp parallised
	pool.close()
#### close grompp multiprocessing and change to min directory and spin up mdrun multiprocessing
	os.chdir('min')
	pool = mp.Pool(mp.cpu_count())
	pool.map_async(gromacs, [(gmx+' mdrun -v -nt 1 -s '+residue_type+'_'+str(rid)+' -deffnm '+residue_type+'_'+str(rid)+' -c '+residue_type+'_'+str(rid)+'.pdb') \
for rid in range(0, resid)]).get()
	pool.close()
	os.chdir(working_dir)

def merge_minimised(residue_type):
	os.chdir(working_dir+residue_type+'/min')
	print('\nMerging non protein atomistic files')
#### create merged pdb in min folder
	pdb_output=create_pdb(working_dir+residue_type+'/min/'+residue_type+'_merged.pdb')	
	if residue_type =='SOL':
		resid_range=1
	else:
		resid_range=np_system[residue_type]
#### run through every resid 
	for resid in range(resid_range):
	#### check if it exists
		if os.path.exists(working_dir+residue_type+'/min/'+residue_type+'_'+str(resid)+'.pdb'):
		#### read in resid and write straight to merged pdb
			with open(working_dir+residue_type+'/min/'+residue_type+'_'+str(resid)+'.pdb', 'r') as pdb_input:
				for line in pdb_input.readlines():
					if line.startswith('ATOM'):
						pdb_output.write(line)
		else:
			sys.exit('cannot find minimised residue: \n'+ working_dir+residue_type+'/min/'+residue_type+'_'+str(resid)+'.pdb')		
	pdb_output.write('TER\nENDMDL')
	pdb_output.close()



############################################################ Build Protein Section ################################################################


def BB_connectivity(at_connections,cg_connections, cg_residues, at_residues, residue_number, BB_connect, res, center):
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
	return at_connections,cg_connections, res

def build_protein_atomistic_system(cg_residues, box_vec):
#### initisation of counters
	chain_information=[]
	chain_count=0
	at_counter=1
	system={}
	at_residues={}
	backbone_coords={}
	backbone_coords[chain_count]=[]
	res=0
	print('Converting Protein')
	mkdir_directory(working_dir+'PROTEIN')	### make and change to protein directory
#### create new pdb file for chain 0 
	pdb_output = create_pdb(working_dir+'PROTEIN/PROTEIN_novo_'+str(chain_count)+'.pdb')
#### for each residue in protein
	initial=True
	for residue_number in cg_residues:
	#### temporary index/dictionaries	
		
		final_at_residues={}  
		at_residues[residue_number]={}
	#### fetch fragments in residue and connectivity
		at_residues[residue_number], connect=get_atomistic_fragments(cg_residues[residue_number][next(iter(cg_residues[residue_number]))]['residue_name'],cg_residues[residue_number], residue_number)
	#### if residue contains BB bead a index of the BB connectivity is collected
		if 'BB' in at_residues[residue_number]:
			BB_connect=[] ### backbone connectivity
			for atom_num, atom in enumerate(at_residues[residue_number]['BB']):
				if at_residues[residue_number]['BB'][atom]['atom'] in backbone[cg_residues[residue_number]['BB']['residue_name']]['b_connect']:
					BB_connect.append(str(atom))
		
	#### for each bead in residue
		for frag_val,cg_fragments in enumerate(cg_residues[residue_number]):
		#### gets connectivity between fragents	
			at_connections, cg_connections, center=connectivity(frag_val, cg_fragments, connect, at_residues[residue_number], cg_residues, residue_number)
		#### if BB bead adds the N and C terminal atoms to connectivity
			if cg_fragments=='BB':
#art#################################
				# if len(at_connections) == 0:
					##### write a function to a add a artifial sidechain to fix randon rotation of GLY/ALA
#art###############################
				at_connections, cg_connections, res = BB_connectivity(at_connections,cg_connections, cg_residues, at_residues, residue_number, BB_connect, res, center)
			#### measures the distance between BB beads. 
				if not initial:
					xyz_prev=[cg_residues[residue_number-1]['BB']['coord'][0],cg_residues[residue_number-1]['BB']['coord'][1],cg_residues[residue_number-1]['BB']['coord'][2]]				
					xyz_cur=[cg_residues[residue_number]['BB']['coord'][0],cg_residues[residue_number]['BB']['coord'][1],cg_residues[residue_number]['BB']['coord'][2]]
					dist=np.sqrt(((xyz_prev[0]-xyz_cur[0])**2)+((xyz_prev[1]-xyz_cur[1])**2)+((xyz_prev[2]-xyz_cur[2])**2))
				#### if distance between BB beads is more than 5 A then it is considered a new chain.
					if dist > 5:
						chain_count+=1  ### adds to to the protein count
						backbone_coords[chain_count]=[]   #### creates another dictionary key for bb fragments 
						backbone_coords[chain_count].append(xyz_cur+[1])  #### adds xyz coord and mass of 1 to list
						chain_information.append([dist, residue_number, residue_number-res])  ### info for verbose flag
						res=residue_number-1 #### updates residue
					#### creates a new protein pdb 
						pdb_output = create_pdb(working_dir+'PROTEIN/PROTEIN_novo_'+str(chain_count)+'.pdb')
						at_counter=1  ### resets protein_chain atom count
					else:
					#### the xyz coord of the BB bead are added to the backbone_coords dictionary
						backbone_coords[chain_count].append(xyz_cur+[1])
				#### if not prev residue the xyz coord of the cg_bead are added to the backbone_coords dictionary
				else:
					xyz_cur=[cg_residues[residue_number]['BB']['coord'][0],cg_residues[residue_number]['BB']['coord'][1],cg_residues[residue_number]['BB']['coord'][2]]
					backbone_coords[chain_count].append(xyz_cur+[1])
					initial=False
		#### finds optimum rotation of fragement
			xyz_rot_apply=rotate(np.array(at_connections), np.array(cg_connections))
		#### applies rotation to each atom
			for atom in at_residues[residue_number][cg_fragments]:
				at_residues[residue_number][cg_fragments][atom]['coord'] = rotate_atom(at_residues[residue_number][cg_fragments][atom]['coord'], center, xyz_rot_apply)
				final_at_residues[atom]=at_residues[residue_number][cg_fragments][atom]
	#### writes fragment to pdb
		for at_val, atom in enumerate(final_at_residues):
			pdb_output.write(pdbline%((int(at_counter),final_at_residues[str(at_val+1)]['atom'],final_at_residues[str(at_val+1)]['res_type'],ascii_uppercase[chain_count],int(residue_number),\
final_at_residues[str(at_val+1)]['coord'][0],final_at_residues[str(at_val+1)]['coord'][1],final_at_residues[str(at_val+1)]['coord'][2],1,0))+'\n')
			at_counter+=1  #### adds 1 to atom counter
	pdb_output.close()   #### close file write
	if args.v >=1:
		print('\n---------------------------------------------->  Verbose level 1 start')

		print('\nchain number\tDelta A\tno in pdb\tlength of chain')
		print('------------\t-------\t---------\t---------------')
		for chain in range(chain_count):
			print('\t',chain,'\t',np.round(chain_information[chain][0], 1),'\t',chain_information[chain][1]-chain_information[chain][2]+1,'-',chain_information[chain][1],'\t\t',chain_information[chain][2])
		if chain_count >1:
			print('\t', chain_count,'\tN/A\t',chain_information[chain][1]+1,'-',residue_number+1,'\t\t',residue_number-res)
		else:
			print('\t', chain_count,'\tN/A\t',res+1,'-',residue_number+1,'\t\t',residue_number-res+1)

		print('\n---------------------------------------------->  Verbose level 1 end\n')

	system['PROTEIN']=chain_count+1
	return system, backbone_coords


############################################################ Processes atomistic protein input ################################################################


def read_in_atomistic(chain_count):
#### reset location and check if pdb exists  
	os.chdir(start_dir)
	if not os.path.exists(args.a):
		sys.exit('cannot find atomistic protein : '+args.a)
#### read in atomistic fragments into dictionary residue_list[0]=x,y,z,atom_name	
	atomistic_protein_input={}
	chain_count=0
#### read in pdb
	with open(args.a, 'r') as pdb_input:
		atomistic_protein_input[chain_count]={}
		for line_nr, line in enumerate(pdb_input.readlines()):
			#### separate line dependant on gro or pdb
			run=False ## turns to true is line is a bead/atom
			if args.c.endswith('gro') and len(line.split())==6:
				line_sep = groatom(line)
				run=True
			elif args.c.endswith('pdb') and line.startswith('ATOM'):
				line_sep = pdbatom(line)
				run=True
			#### if line is correct
			if run:
				if 'H' not in line_sep['atom_name'][0] or line_sep['residue_name'] in mod_residues:  
				#### sorts out wrong atoms in terminal residues
					if line_sep['atom_name'] in ['OT', 'O1', 'O2']:
						line_sep['atom_name']='O'
				#### makes C_terminal connecting atom variable  
					if line_sep['atom_name'] == backbone[line_sep['residue_name']]['b_connect'][1]:
						C_ter=[line_sep['x'],line_sep['y'],line_sep['z']]
						C_resid=line_sep['residue_id']
						C=True
					try:
					#### tries to make a N_terminal connecting atom variable
						if line_sep['atom_name'] == backbone[line_sep['residue_name']]['b_connect'][0]:
							N_resid=line_sep['residue_id']
							N_ter=[line_sep['x'],line_sep['y'],line_sep['z']]
							N=True
					#### measures distance between N and C atoms. if the bond is over 3 A it counts as a new protein
						dist=np.sqrt(((N_ter[0]-C_ter[0])**2)+((N_ter[1]-C_ter[1])**2)+((N_ter[2]-C_ter[2])**2))
						if N and C and C_resid != N_resid and dist > 3:
							N_ter, C_ter=False, False
							chain_count+=1
							atomistic_protein_input[chain_count]={} ### new chain key
					except:
						pass
					if line_sep['residue_id'] not in atomistic_protein_input[chain_count]:  ## if protein does not exist add to dict
						atomistic_protein_input[chain_count][line_sep['residue_id']]={}
				#### adds atom to dictionary, every atom is given a initial mass of zero 
					atomistic_protein_input[chain_count][line_sep['residue_id']][line_sep['atom_number']]={'coord':np.array([line_sep['x'],line_sep['y'],line_sep['z']]),'atom':line_sep['atom_name'], 'res_type':line_sep['residue_name'],'frag_mass':0, 'resid':line_sep['residue_id']}
				#### if atom is in the backbone list then its mass is updated to the correct one
					if line_sep['atom_name'] in backbone[line_sep['residue_name']]['atoms']:
						for atom in line_sep['atom_name']:
							if atom in mass:
								atomistic_protein_input[chain_count][line_sep['residue_id']][line_sep['atom_number']]['frag_mass']=mass[atom]
#### check if number of monomers is the same
	if chain_count+1 != system['PROTEIN']:
		sys.exit('number of chains in atomistic protein input ('+str(chain_count+1)+') does not match CG representation ('+str(system['PROTEIN']))
	return atomistic_protein_input

def center_atomistic(atomistic_protein_input): 
	cg_com=[]
#### for each protein chain center on cg representation 
	for chain in range(system['PROTEIN']):
		protein_mass=[]
		for residue in atomistic_protein_input[chain]:
		#### creates a list of all coordinates and masses [[coord, mass],[coord, mass]]
			for atom in atomistic_protein_input[chain][residue]:
				protein_mass.append([atomistic_protein_input[chain][residue][atom]['coord'][0],atomistic_protein_input[chain][residue][atom]['coord'][1],\
					atomistic_protein_input[chain][residue][atom]['coord'][2],atomistic_protein_input[chain][residue][atom]['frag_mass']])
	#### returns the COM of the atomistic protein
		atomistic_protein_mass=np.average(np.array(protein_mass)[:,:3], axis=0, weights=np.array(protein_mass)[:,3])#### add center of mass of CG_proteins
	#### for each chain the COM of the CG representation is stored (only cg is needed)
		cg_com.append(np.average(np.array(backbone_coords[chain])[:,:3], axis=0, weights=np.array(backbone_coords[chain])[:,3]))
	#### each atoms coord is updated so the monomer COM is the same as the CG
		for residue in atomistic_protein_input[chain]:
			for atom in atomistic_protein_input[chain][residue]:
				atomistic_protein_input[chain][residue][atom]['coord']=atomistic_protein_input[chain][residue][atom]['coord']-(atomistic_protein_mass-cg_com[chain])
	return atomistic_protein_input, cg_com

def rotate_protein_monomers():
#### run through each chain in proteins
	for chain in range(system['PROTEIN']):
	#### creates atomistic pdb
		pdb_output = create_pdb(working_dir+'PROTEIN/PROTEIN_at-input_'+str(chain)+'.pdb')
		at_centers=[]
	#### runs through every residue and atom  
		for residue in atomistic_protein_input[chain]:
		#### gets center of mass of each residue (note only backbone heavy atoms have a mass)
			at_centers_iter=[]
			for atom in atomistic_protein_input[chain][residue]:
				at_centers_iter.append(np.append(atomistic_protein_centered[chain][residue][atom]['coord'],atomistic_protein_centered[chain][residue][atom]['frag_mass']))
			
			try:
				at_centers.append(np.average(np.array(at_centers_iter)[:,:3], axis=0, weights=np.array(at_centers_iter)[:,3]))
			except:
				for atom in atomistic_protein_input[chain][residue]:
					print(atomistic_protein_input[chain][residue][atom])
				sys.exit()
	#### checks that the number of residues in the chain are the same between CG and AT
		if len(at_centers) != len(backbone_coords[chain]):
			sys.exit('In chain '+str(chain)+' the atommistic input does not match the CG. \n\
number of CG residues '+str(len(backbone_coords[chain]))+'\nnumber of AT residues '+str(len(at_centers)))
	#### finds optimal rotation of each monomer  
		xyz_rot_apply = rotate(np.array(at_centers)-cg_com[chain], np.array(backbone_coords[chain])[:,:3]-cg_com[chain])
		if args.v >= 1:
			print('\nThe proteins chains are rotated around the COM of all the backbone heavy atoms.')
			print('The COM of chain', chain,'is :', np.round(cg_com[chain][0], 2),',', np.round(cg_com[chain][1], 2),',', np.round(cg_com[chain][2], 2))
			print('rotating chain ', chain, 'by :',np.round(np.degrees(xyz_rot_apply[0]),2),',',np.round(np.degrees(xyz_rot_apply[1]),2),',',np.round(np.degrees(xyz_rot_apply[2]),2))
			print()
	#### applies optimal rotation to each atom 
		for residue in atomistic_protein_centered[chain]:
			for atom in atomistic_protein_centered[chain][residue]:
				atomistic_protein_centered[chain][residue][atom]['coord'] = rotate_atom(atomistic_protein_centered[chain][residue][atom]['coord'], cg_com[chain], xyz_rot_apply)
			#### writes out new pdb for each optimised chain
				pdb_output.write(pdbline%((int(atom),atomistic_protein_centered[chain][residue][atom]['atom'],atomistic_protein_centered[chain][residue][atom]['res_type'],\
					ascii_uppercase[chain],atomistic_protein_centered[chain][residue][atom]['resid'],atomistic_protein_centered[chain][residue][atom]['coord'][0],\
					atomistic_protein_centered[chain][residue][atom]['coord'][1],atomistic_protein_centered[chain][residue][atom]['coord'][2],1,0))+'\n')
	return atomistic_protein_centered

######################################################################## GROMACS protein ###################################################################

def minimise_protein():
	os.chdir(working_dir+'/PROTEIN')
	mkdir_directory(working_dir+'FORCEFIELD')
	copy_tree(forcefield_location+forcefield+'.ff', working_dir+'PROTEIN/'+forcefield+'.ff/.')
	mkdir_directory('min')
	for chain in range(system['PROTEIN']):
		with open('em_'+str(chain)+'.mdp','w') as em:
			em.write('define = -DPROTEIN_'+str(chain)+'_CA_posre.itp\nintegrator = steep\nnsteps = 10000\nemtol = 1000\nemstep = 0.001\ncutoff-scheme = Verlet\n')
		if args.a != None:
			minimise_protein_chain(chain, 'at-input_')
		minimise_protein_chain(chain, 'novo_')
	os.chdir('..')
	merge_protein(system['PROTEIN'], '_novo')
	if args.a != None:
		merge_protein(system['PROTEIN'], '_at-input')

def minimise_protein_chain(chain, input):
	copyfile(working_dir+'PROTEIN/'+forcefield+'.ff/residuetypes.dat', 'residuetypes.dat')
	gromacs(gmx+' pdb2gmx -f PROTEIN_'+input+str(chain)+'.pdb -o PROTEIN_'+input+str(chain)+'_gmx.pdb -water none \
-p PROTEIN_'+input+str(chain)+'.top  -i PROTEIN_'+input+str(chain)+'_posre.itp << EOF \n1\nEOF') ### single chains
	convert_topology('PROTEIN_'+input, chain)
	write_topol('PROTEIN_'+input, 1, str(chain))
	gromacs(gmx+' grompp -f em_'+str(chain)+'.mdp -p PROTEIN_'+input+str(chain)+'.top -c PROTEIN_'+input+str(chain)+'_gmx.pdb -o min/PROTEIN_'+input+str(chain)+' -maxwarn 1')
	os.chdir('min')
	gromacs(gmx+' mdrun -v -nt 10 -deffnm PROTEIN_'+input+str(chain))
	gromacs(gmx+' editconf -f PROTEIN_'+input+str(chain)+'.gro -o PROTEIN_'+input+str(chain)+'.pdb -pbc')
	os.chdir('..')	

def write_topol(residue_type, residue_number, chain):
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
		topol_write.write(residue_type+chain+'    '+str(residue_number))

def merge_protein(resid, protein):
	merge=[]
	for res_val in range(0,resid):
		if os.path.exists(working_dir+'PROTEIN/min/PROTEIN'+protein+'_'+str(res_val)+'.pdb'):
			posres_output = open(working_dir+'PROTEIN/PROTEIN_'+str(res_val)+'_CA_posre.itp', 'w')
			posres_output.write('[ position_restraints ]\n; atom  type      fx      fy      fz\n')
			with open(working_dir+'PROTEIN/min/PROTEIN'+protein+'_'+str(res_val)+'.pdb', 'r') as pdb_input:
				at_counter=0
				for line in pdb_input.readlines():
					if line.startswith('ATOM'):
						merge.append(line)
						line_sep = pdbatom(line)
						at_counter+=1
						if line_sep['atom_name'] in ['CA', 'CB']:
							posres_output.write(str(at_counter)+'     1  500  500  500\n')

						if line_sep['residue_name'] in mod_residues:
							if 'H' not in line_sep['atom_name'] and line_sep['atom_name'] not in backbone[line_sep['residue_name']]['atoms']:
								posres_output.write(str(at_counter)+'     1  500  500  500\n')
			posres_output.close()
		else:
			sys.exit('cannot find minimised residue: \n'+'PROTEIN/min/PROTEIN'+protein+'_'+str(res_val)+'.pdb')	
	pdb_output=create_pdb(working_dir+'PROTEIN/min/PROTEIN'+protein+'_merged.pdb')
	for line in merge:
		pdb_output.write(line)
	pdb_output.close()

def convert_topology(topol, protein_number):
	if Path(topol+str(protein_number)+'.top').exists():
		read=False
		with open(topol+str(protein_number)+'.itp', 'w') as itp_write:
			for line in open(topol+str(protein_number)+'.top', 'r').readlines():
				if len(line.split()) > 1: 
					if read == False and line.split()[1] == 'moleculetype':
						read = True
					if read == True and line.split()[1] == 'POSRES':
						read = False
					if line.split()[0] == 'Protein' or line.split()[0][:-1] == 'Protein_chain_':
						line= re.sub('Protein_chain_.?', 'Protein',line)
						line= re.sub('Protein', 'protein_'+str(protein_number),line)
				if read:
					itp_write.write(line)
			itp_write.write('#ifdef POSRES\n#include \"PROTEIN_'+str(protein_number)+'_posre.itp\"\n#endif\n') 
			itp_write.write('\n; Include CA Position restraint file\n#ifdef POSRES_CA\n#include \"PROTEIN_'+str(protein_number)+'_CA_posre.itp\"\n#endif')
	else:
		sys.exit('cannot find : '+topol+'_'+str(protein_number)+'.top')


def merge_system_pdbs(system, box_vec, protein):
	os.chdir(working_dir+'MERGED')
	pdb_output=create_pdb(working_dir+'MERGED/merged_cg2at'+protein+'.pdb')	
	for residue_type in cg_residues:
		if residue_type =='ION' and 'SOL' not in cg_residues:
			cg_residues['SOL']=[]
		if residue_type != 'PROTEIN':
			protein=''
		if os.path.exists(working_dir+residue_type+'/min/'+residue_type+protein+'_merged.pdb'):
			with open(working_dir+residue_type+'/min/'+residue_type+protein+'_merged.pdb', 'r') as pdb_input:
				for line in pdb_input.readlines():
					if line.startswith('ATOM'):
						pdb_output.write(line)
		else:
			sys.exit('cannot find minimised residue: \n'+ working_dir+residue_type+'/min/'+residue_type+protein+'_merged.pdb')		
	pdb_output.write('TER\nENDMDL')
	pdb_output.close()


#################################################################### GROMACS non protein ###################################################################



def write_merged_topol(system, box_vec, protein):
	os.chdir(working_dir+'MERGED')
	for topology in ['topol_final', 'topol_alchembed']:
		with open(topology+'.top', 'w') as topol_write:
			topol_write.write('; Include forcefield parameters\n#include \"'+working_dir+'FORCEFIELD/charmm36-jul2017-updated.ff/forcefield.itp\"\n')
			topol_write.write('#include \"'+working_dir+'/FORCEFIELD/'+forcefield+'.ff/tip3p.itp\"\n\n#include \"'+working_dir+'/FORCEFIELD/'+forcefield+'.ff/ions.itp\"\n\n')
			for residue_type in system:
				found=False
				if residue_type not in ['ION','SOL']:
					for directory in np_directories:
						if os.path.exists(directory[0]+residue_type+'/'+residue_type+'.itp'):		
							topol_write.write('#include \"'+residue_type+'.itp\"\n')
							copyfile(directory[0]+residue_type+'/'+residue_type+'.itp', residue_type+'.itp')
							found=True
					if residue_type == 'PROTEIN':
						if topology == 'topol_alchembed':
							topol_write.write('#include \"PROTEIN.itp\"\n')
						else:
							for protein_unit in range(system[residue_type]): 
								topol_write.write('#include \"PROTEIN_'+str(protein_unit)+'.itp\"\n')
								copyfile(working_dir+'PROTEIN/PROTEIN'+protein+'_'+str(protein_unit)+'.itp', 'PROTEIN_'+str(protein_unit)+'.itp')
								copyfile(working_dir+'PROTEIN/PROTEIN_'+str(protein_unit)+'_CA_posre.itp', 'PROTEIN_'+str(protein_unit)+'_CA_posre.itp')
								copyfile(working_dir+'PROTEIN/PROTEIN'+protein+'_'+str(protein_unit)+'_posre.itp', 'PROTEIN_'+str(protein_unit)+'_posre.itp')
			topol_write.write('[ system ]\n; Name\nSomething clever....\n\n[ molecules ]\n; Compound        #mols\n')


			for residue_type in system:
				if residue_type not in  ['PROTEIN']:
					topol_write.write(residue_type+'    '+str(system[residue_type])+'\n')	
				if residue_type == 'PROTEIN':
					if topology == 'topol_alchembed':
						topol_write.write('PROTEIN    1\n')	
					else:
						for protein_unit in range(system[residue_type]):
							 topol_write.write('PROTEIN_'+str(protein_unit)+'    1\n')	


############################################################################  GROMACS system ############################################################							 

def minimise_merged_pdbs(system, box_vec, protein):
	print('Minimising merged atomistic files : '+protein[1:])
	os.chdir(working_dir+'MERGED')
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
	mkdir_directory('alchembed')
	for chain in range(system['PROTEIN']):
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
	mkdir_directory('steered_md')
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
		if os.path.exists('merged_cg2at_steered_md.gro'):
			gromacs(gmx+' editconf -f merged_cg2at_steered_md.gro -o merged_cg2at_steered_md.pdb -pbc')
			copyfile('merged_cg2at_steered_md.pdb', final_dir+'final_cg2at_CA_CB_steered_md.pdb')
		else:
			print('steered MD failed! Starting atomistic input may be too far from CG structure')


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

gmx='gmx'
pdbline = "ATOM  %5d %4s %4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f"
mass = {'H': 1,'C': 12,'N': 14,'O': 16,'S': 32,'P': 31,'M': 0, 'B': 32}

### sets up file locations
timestamp =  strftime("%Y-%m-%d_%H-%M-%S", gmtime())
start_dir=os.getcwd()
working_dir=os.getcwd()+'/CG2AT_'+timestamp+'/'
final_dir=os.getcwd()+'/CG2AT_'+timestamp+'/FINAL/'
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
np_residues, p_residues, mod_residues,np_directories, p_directories, mod_directories, water=fetch_residues()
backbone=fetch_fragment(p_directories)

fragment_selection_time=time.time()
print('This script is now hopefully doing the following:\n')

#### read in CG file
print('Reading in your CG representation')
cg_residues, box_vec = read_initial_pdb()
os.chdir(working_dir)
read_in_time=time.time()

system={}
#### converts protein into atomistic and minimises
print('\nConverting CG into a atomistic representation')
if 'PROTEIN' in cg_residues:
	if len(cg_residues['PROTEIN'])>0:
		p_system, backbone_coords=build_protein_atomistic_system(cg_residues['PROTEIN'], box_vec)
		system.update(p_system)
		protein_de_novo_time=time.time()
		if args.a != None:
		#### reads in atomistic structure	
			atomistic_protein_input = read_in_atomistic(system['PROTEIN'])	
			atomistic_protein_centered, cg_com = center_atomistic(atomistic_protein_input)
			atomistic_protein_rotated = rotate_protein_monomers()
		minimise_protein()
		
final_protein_time=time.time()


#### converts non protein residues into atomistic and minimises 
if len([key for value, key in enumerate(cg_residues) if key not in ['PROTEIN']]) > 0:
	np_system, fragment_names=build_atomistic_system(cg_residues, box_vec)
	for residue_type in cg_residues:
		if residue_type =='ION' and 'SOL' not in cg_residues:
			non_protein_minimise(np_system['SOL'], 'SOL', '')
		elif residue_type in np_system:
			non_protein_minimise(np_system[residue_type], residue_type, fragment_names[residue_type])
		if residue_type not in ['PROTEIN', 'ION']:
			merge_minimised(residue_type)
	system.update(np_system)
	build_non_protein_time=time.time()

non_protein_time=time.time()

#### creates merged folder
print('Merging all residue types to single file')

if len(system)>0:
	mkdir_directory(working_dir+'MERGED')
#### merges provided atomistic protein and residues types into a single pdb file into merged directory
	if args.a != None:
		merge_system_pdbs(system, box_vec,'_at-input' )
		minimise_merged_pdbs(system, box_vec, '_at-input')
		alchembed(p_system, box_vec, '_at-input')
#### merges de novo protein and residues types into a single pdb file into merged directory
	merge_system_pdbs(system, box_vec,'_novo' )
	minimise_merged_pdbs(system, box_vec, '_novo')
	copyfile('merged_cg2at_novo_minimised.pdb', final_dir+'final_cg2at_novo_minimised.pdb')
	merge_time=time.time()
	if args.a != None:
	#### runs steered MD on atomistic structure on CA and CB atoms
		steered_md_atomistic_to_cg_coord()
		steered_md_time=time.time()
#### copies all itp files and topologies from whereever they are stored
	for file_name in os.listdir(working_dir+'MERGED'):
		if file_name.endswith('.itp') or file_name.endswith('.top'):
			copyfile(working_dir+'MERGED/'+file_name, final_dir+file_name)
#### if atomistic structure has been supplies 

final_time=time.time()




#### prints out system information

print('\n{:-<100}'.format(''))
print('{0:^100}'.format('Script has completed'))
print('\n{0:^20}{1:^10}'.format('molecules','number'))
print('{0:^20}{1:^10}'.format('---------','------'))
for section in system:
	# print('\t',str(section),'\t', +str(system[section])).format())
	print('{0:^20}{1:^10}'.format(section, system[section]))

#### prints out script timings for each section
print('\nInitialisation: ', str(datetime.timedelta(minutes=np.round(initialisation_time-start_time, 2))).rsplit(':', 1)[0], ' min',\
'\nFragment selection: ', str(datetime.timedelta(minutes=np.round(fragment_selection_time-initialisation_time, 2))).rsplit(':', 1)[0], ' min' ,\
'\nRead in CG system: ', str(datetime.timedelta(minutes=np.round(read_in_time-fragment_selection_time, 2))).rsplit(':', 1)[0], ' min')
if args.a != None:
	print('Build de novo protein system: ', str(datetime.timedelta(minutes=np.round(protein_de_novo_time-read_in_time, 2))).rsplit(':', 1)[0], ' min',\
	'\nBuild protein system from provided structure: ', str(datetime.timedelta(minutes=np.round(final_protein_time-protein_de_novo_time, 2))).rsplit(':', 1)[0], ' min', \
	'\nTotal protein system build: ', str(datetime.timedelta(minutes=np.round(final_protein_time-read_in_time, 2))).rsplit(':', 1)[0], ' min')
else:
	print('Build de novo protein system: ', str(datetime.timedelta(minutes=np.round(final_protein_time-read_in_time, 2))).rsplit(':', 1)[0], ' min')
print('Build non protein system: ', str(datetime.timedelta(minutes=np.round(non_protein_time-final_protein_time, 2))).rsplit(':', 1)[0], ' min', \
'\nMerge protein and non protein system: ', str(datetime.timedelta(minutes=np.round(merge_time-non_protein_time, 2))).rsplit(':', 1)[0], ' min')
if args.a != None:
	print('Steered MD atomistic input to match CG: ', str(datetime.timedelta(minutes=np.round(steered_md_time-merge_time, 2))).rsplit(':', 1)[0], ' min')
print('Total run time: ', str(datetime.timedelta(minutes=np.round(final_time-start_time, 2))).rsplit(':', 1)[0], ' min')