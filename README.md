                                   **CG2AT v2 README**

If you are using this script please acknowledge me (Dr Owen Vickery) and cite the following DOI.

DOI: 10.5281/zenodo.3592318

                                    **SCRIPT OVERVIEW**

the script is a fragment based conversion of CG representation into atomistic.


                                     **REQUIREMENTS**

- Python v3 or higher
- Numpy
- Scipy
- GROMACS > v5

                                        **FLAGS**

REQUIRED
-c          coarse grain input (pdb/gro/tpr)

OPTIONAL
-a          atomistic_input.pdb (pdb/gro/tpr)
-v          verbose with 3 levels  of verbosity (-vvv) 
-mod        treats every fragment individually 
-sf         scale factor for fragments, shrinks fragments before minimisation (default=0.9)
-clean      removes temp file after cg2at 
-ter        interactive terminal residue selection
-nt         charged N terminal residue
-ct         charged C terminal residue
-swap       allows the swapping of residues/beads (eg PIP2,D3A:PIP2,C3A GLU,SC2:ASP,skip)
-box        new box vectors in Angstrom, (0 is ignored and uses input eg 100 150 100 )
-w          choose your solvent, common choices are: tip3p, tip4p, spc and spce. 
-ff         choose your forcefield. 
-fg         choose your fragment libraries. 
-gromacs    provide gromacs executable name if none standard. 






you can supply 2 input files:
coarse grain input as a pdb/gro/tpr
atomistic representation as a pdb/gro/tpr

The script has 3 outputs:

The 1st requires no information other than the fragment corresponding to the bead
The script will convert any protein using a basic de novo method of fragment alignment to the CG beads.

The 2nd maps a users supplied atomistic structure onto the the CG representation. 
The mapped structure is then steered for a very short period to fit the CG coordinates.

Note:
the mapping of the user structure now includes a hybrid approach.
If any residues are missing from the supplied atomistic structure they are built in using the de novo structure.

for example if you have a signal peptide linked via a flexible linker to the main protein.
You can just supply the coordinates for the signal peptide and main protein and the script will build in the linker from the de novo.

In martini this is a single protein. (eg a membrane signal helix linked to to a cytosolic protein)

You only have atomistic structures for the helix and the cytoslic protein.
the script will rigid body fit the helix and cytosolic protein separately and fill in the linker from the de novo construct.

                              CG                        AT from user
    input            helix --------- protein           helix  protein 
                               |                             |
                               V                             |
    conversion       helix --------- protein                 |
     (de novo)                 |                             |
                               V                             V
           de novo linker  ---------             user   helix, protein
                                    \                  /
                                             |
                                             V
    final hybrid system            helix --------- protein



The 3rd just does a rigid body fit on the backbone.
There is no alchembed run on this system as you can use it as a restraint file and pull the steered coordinates into the rigid fit.

This can be done by:
add "define = -DPOSRES" to the MDP file applies restraints to all heavy atoms
gmx grompp -f md.mdp -c final_cg2at_at_rep_user_supplied.pdb -r final_cg2at_no_steered.pdb -p topol_final.top -o rigid_fit.tpr
gmx mdrun -v -deffnm rigid_fit


#############################################################################




example input.

python cg2at.py -c cg_input.pdb -a atomistic_input.pdb -v -clean -w tip3p -fg martini_2-2_charmm36 -ff charmm36-jul2017-update -swap GLU,SC2:ASP,skip PIP2,D3B:PIP2,C3B -box 0 100 100 -nt -ct

############################################################################

Available fragments

Charmm36
________

Non protein
-----------
SOL: 'tip3p', 'tip4p', 'spc', 'spce'
PC:  'POPC', 'DGPC', 'DLPC', 'DNPC', 'DOPC', 'DPPC', 'DRPC', 'DYPC'
PE:  'POPE', 'DGPE', 'DLPE', 'DNPE', 'DOPE', 'DPPE', 'DYPE'
PG:  'POPG', 'DLPG', 'DNPA', 'DOPA', 'DPPA', 'DYPA',
PA:  'POPA', 'DGPA', 'DLPA', 'DNPA', 'DOPA', 'DPPA', 'DYPA' 
PS:  'POPS', 'DGPS', 'DLPS', 'DMPS', 'DNPS', 'DOPS', 'DPPS', 'DYPS'
PI:  'PIP2', 'POPI2C', 'POPI2D','DMPI',  
Sterols: 'CHOL'
Lipopolysaccharides: 'LIPA'
Cardiolipin: 'CARD' (PVCL2) 
misc: 'SPA', 'SPU'


Protein
-------
all amino acids

diacyl (CYSD) and triacyl (CYST) cysteine

-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

Amber99sb
_________

Protein
-------
All amino acids

-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

lipid-14
________

Non protein
-----------
Currently empty


-----------------------------------------------------------------------------

Forcefields

amber14sb.ff
amber99sb-ildn.ff
charmm36-mar2019.ff
charmm36-jul2017-updated.ff
oplsaam.ff


#############################################################################

Script requirements

Gromacs >= v5

Non standard modules

numpy
scipy

standard modules included in base python install

argparse
copy
datetime
difflib
distutils
glob
math
multiprocessing
os
pathlib
re
shutil
string
sys
subprocess
time


#############################################################################

Script overview

 Read in databases 
        |-> choose forcefield and fragments
                      |-> fetch available fragements and bonded information

 collate input files as pdb into INPUT directory 
                      |-> Read in CG input
 
 Build non protein atomistic system
                |-> align COM of fragments to bead
                                | If not SOL or ION
                                |-> rotate fragment to minimise distance between connecting atoms the connected beads 
                                |                                       |-> energy minimise each residue separately
                                |                                                         |-> merge all minimised residues 
                                |                                                                         |-> minimise merged residue type incase of any clashes 
                                | If ION or SOL
                                |-> Randomly rotate water around COM of water/ion bead
                                                         |-> concatonate water molecules in one pdb
                                                         |              |-> minimise all water molecules
                                                         |if ion
                                                         |-> separate ions into separate file (no need to do anything to them) 

Build de novo protein atomistic system
                |-> align COM of fragments to bead
                                |-> rotate fragment to minimise distance between connecting atoms the connected beads 
                                                                           |-> each chain is minimised separately
                                                                                               |-> merge all minimised chains
                                                                                                            |-> minimise merged chains incase of any clashes                             
Build user supplied protein atomistic system
                |-> align COM of user supplied protein chain to CG chain
                                         |-> fit backbone COMs of user supplied protein chain to CG backbone beads
                                                                        |-> each chain is minimised separately
                                                                                        |-> Run steered MD on selected atoms to fit the user supplied protein chain more accurately to the CG protein_mass

 Merge protein and non protein 
           |if de novo
           |-> minimise merged system
                        |-> copy to final directory 
           |if user supplied protein
           |-> run alchembed on each chain iteratively to prevent any stuck acyl tails 
                                            |-> copy to final directory
