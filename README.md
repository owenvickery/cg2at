[![License: GPL v3](https://img.shields.io/badge/License-GPL_v3-blue)](https://www.gnu.org/licenses/gpl-3.0.html)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/owenvickery/cg2at.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/owenvickery/cg2at/alerts/)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/owenvickery/cg2at.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/owenvickery/cg2at/context:python)
[![Github Releases](https://img.shields.io/github/v/release/owenvickery/cg2at)](https://github.com/owenvickery/cg2at/releases)
<p align="center">
                                   <b>**CG2AT v2 README**</b>
</p>
If you are using this script please acknowledge me (Dr Owen Vickery) and cite the following DOI.

DOI: 10.5281/zenodo.3890163

<p align="center">
                                   <b>**CG2AT SCRIPT OVERVIEW**</b>
</p>

The script is a fragment based conversion of the CG representation into atomistic. 

This script roughly follows the following workflow.

- Center fragments based on the center of mass (COM) of heavy atoms.
- Rotate fragments to find minimum distance between the atoms connecting to other beads.
- Minimise residue
- Merge all residues and minimise
- check for abnormal bonds (marker for lipid through rings)
- Run NVT and NPT
- Morph protein to either rigid body alignment or steered alignment


For a detailed description of the script workflow please read the file in 
<b>
database/script_files/workflow.pdf
</b>

![Alt text](./database/script_files/images/Fragments.png)


<pre>
           CG beads         Fragments           COM aligned fragments           Aligned fragments                 Atomistic
                                                                              
           --------                                    --------                      --------      
          (        )                                  (        )                    (        )     
         (          )             O1                 (      O1  )                  (  O1  O2  )                     
        (            )           /                  (       /    )                (    \  /    )                      
       (     SC1      )    [CB]-CG                 (  [CB]-CG     )              (      CG      )                     
        (            )           \                  (       \    )                (     |      )                    O1  O2
         (          )             O2                 (      O2  )                  (   [CB]   )                      \  /
          (        )                                  (        )                    (        )                        CG
           --------                       COM          --------       rotation       --------      Minimisation       |
               |                      ---------->          |         ---------->         |         ----------->       CB  
           --------                    Alignment       --------       Alignment      --------                         | 
          (        )                                  (        )                    (        )                    X1  CA   X2
         (          )          [C]-O                 (   [C]-O  )                  (    CA    )                    \ /  \ /
        (            )         /                    (    /       )                (    /  \    )                    N    C
    X1-(      BB      )-X2   [CA]               X1-(   [CA]       )-X2        X1-(   (N)  [C]   )-X2                     |
        (            )         \                    (    \       )                (        |   )                         O
         (          )          [N]                   (   [N]    )                  (       O  )    
          (        )                                  (        )                    (        )     
           --------                                    --------                      --------       

The connecting atoms are highlighted by square brackets

</pre> 

This workflow allows each fragment to be treated individually, with no knowledge of what any other bead contains.

<p align="center">
                                   <b>**REQUIREMENTS**</b>
</p>
                                     
Software:

- Python v3 or higher
- GROMACS > v5

Non standard python modules:

- Numpy
- Scipy

Standard modules included in base python install (ubuntu 18):

- argparse
- copy
- datetime
- difflib
- distutils
- glob
- math
- multiprocessing
- os
- pathlib
- re
- shutil
- string
- sys
- subprocess
- time

<p align="center">
                                   <b>**FLAGS**</b>
</p>
                                        

REQUIRED
- -c          (pdb/gro/tpr)

OPTIONAL
- -a          (pdb/gro/tpr)
- -o          ['all', 'steer', 'align', 'none']
- -group      (str)(e.g. 0,1 2,3 or all or chain)
- -npcus      (int)
- -mod        (True/False)
- -sf         (float)(default=0.9)
- -swap       (str)
- -cys        (float)
- -ter        (True/False)
- -nt         (True/False)
- -ct         (True/False)
- -al         (True/False)
- -vs         (True/False)
- -swap       (str)(e.g. PIP2,D3A:PIP2,C3A GLU,SC2:ASP,skip)
- -box        (int) (Angstrom)(0 is ignored and uses input eg 100 150 100 ) 
- -w          [tip3p, tip4p, spc and spce] 
- -ff         (str) 
- -fg         (str) 
- -gromacs    (str) 
- -loc        (str)  
- -messy      (True/False) 
- -info       (True/False)
- -version    (True/False)
- -v          (-vvv) 

<p align="center">
                                   <b>**INPUT**</b>
</p>

The script contains two methods to rebuild the system:

- de novo method, following the protocol described in the script overview section.
- flexible fitting of a user supplied atomistic structure (prefereably the one used to make the CG representation).

To convert your system, all you need to do is supply the coarsegrain system.

However, the quality of the conversion is improved if you can provide the starting atomistic structure used to create the CG model.

All user supplied structures should be in a single file. The atomistic segments are then aligned to the coarsegrain initially by sequence.

Both input files must be in the following file formats:

- pdb
- gro
- tpr

If only partial structures are supplied, then the script will build in the missing residues from the de novo build. 

For example if you have a signal peptide linked via a flexible linker to the main protein.

You can just supply the atomistic coordinates for the signal peptide and main protein and the script will build in the linker from the de novo method.

the script will rigid body fit the helix and cytosolic protein separately and fill in the linker from the de novo construct. 

<pre>
                              CG                        AT from user
    input            helix --------- protein           helix, protein 
                               |                             |
                               V                             |
    conversion       helix --------- protein                 |
     (de novo)                 |                             |
                               V                             V
          de novo: linker  ---------            user:   helix, protein
                                    \                  /
                                     \________________/
                                             |
                                             V
    final hybrid system            helix --------- protein
</pre>

e.g.
<pre>
    python cg2at.py -c cg_input.gro -a atomistic_input.gro
</pre>

<p align="center">
                                   <b>**Advanced Usage**</b>
</p>

<p align="center">
                                   <b>Number of CPUs</b>
</p>

This script will try to use 8 cores preferentially, if your computer has less than 8 it will use all the cores. 

To override this and limit the number of cores you can use the flag:

- -ncpus (int)

<p align="center">
                                   <b>Extra information</b>
</p>

To print out all available forcefields and fragments use the flag:

- -info

if a fragment library is also specified then the available fragments will also been shown (-fg)
<p align="center">
                                   <b>Version</b>
</p>
                                        
To print out the version of the script use the flag:

- -version

<p align="center">
                                   <b>Verbosity</b>
</p>                                      

There are 3 levels of verbosity within this script, each includes the previous levels :

no verbosity:

- Basic information of current process.
- RMSD of final backbone to original CG representation (Angstrom)
- Number of residues converted

1st level (-v)

- Chain number and length
- Sequence for chains
- COM of each chain
- Rotation of each chain to find optimum fitting
- Timings for each section of the script

2nd level (-vv)

- location of each fragment library 
- all residues available in the fragment library  

3rd level (-vvv)

- prints out all gromacs commands called by the script

<p align="center">
                                   <b>Cleaning up folders</b>
</p>
                                        
To keep all the part files created by the script use the flag:

- -messy

<p align="center">
                                   <b>Virtual sites</b>
</p>

To apply virtual sites to you protein use the flag:

- -vs
<p align="center">
                                   <b>Disulphide Bonds</b>
</p>
                                        

The script currently finds disulphide bonds in the user supplied atomistic structure (S-S < 2.1 A).
As well as searching the CG representation for disulphide bonds (SC1-SC1 < 7 A and more than 4 residues apart).
If the disulphide only exists in the CG, then the script will ask if it is a disuphide. 
To silence the question and automatically select disulphides use the flag: 

- -silent

In most cases the script catches the extra long disulphide bonds in martini simulations, however in some cases the disulphide can extend up to 10A.
If you recieve a error that the pdb and topology don't match and the atom number is out by 2. It is most likely a disulpide bond not being treated correctly.
You may be able to fix it by increasing the disulphide bond search radius catch using the flag:

- -cys (default = 7 A)

<p align="center">
                                   <b>Rigid fitting</b>
</p>
                                        

If you are converting a multimeric protein, such as a potassium channel. You can fit the atomistic structure in several ways:

fit by coarse grain chain.

- -group chain

fit entire atomistic structure rigidly

- -group all

fit chains in specific groups (treat chains 0, 2 and 1, 3 as individual groups).
each group is separated by a space.

- -group 0,2 1,3   

<p align="center">
                                   <b>Swap residues and beads</b>
</p>
                                        

Due to the modular nature of CG representation, you can switch residues during the conversion if you wish to make simple mutations.

If the residue you are swapping to has the same number or fewer CG beads you can use the following flag to change the residue.

- -swap  

Usage

             From        :              To           :    range
      resname,bead,bead  :     resname,bead,bead     :  0-10,30-40

Examples: 

If the beads are the same:

swap all ASP to ASN:
- -swap ASP:ASN

swap ASP to ASN in the resid range 0-10 and 30-40:
- -swap ASP:ASN:0-10,30-40

If the beads are different:

switch residues with different beads:
- -swap POPC,NC3:POPG,GL0

if you wish to switch within the same residue:
- -swap POPG,D2B:POPG,C2B

To skip residues or beads:

to skip a bead.
- -swap GLU,SC2:ASP,skip

to skip a residue.
- -swap POPG:skip

You are not limited to a single residue

The following switches all POPE to POPG and all POPG to POPE
- -swap POPE,NH3:POPG,GL0 POPG,GL0:POPE,NH3

The following will skip all NA+ between resid 4000 and 4100.
- -swap NA+:skip:4000-4100

<p align="center">
                                   <b>**OUTPUT**</b>
</p>
                                        

The script will create a output file system as below.

    | --    CG2AT_(timestamp)
                | --    INPUT
                                - cg_input.gro, conversion_input.pdb, atomistic_input.gro, AT_input.pdb, script_inputs.dat
                | --    RESIDUE_TYPE (PROTEIN, POPE)
                                - converted indivudal residues
                | --    MERGED
                                -  merged residue types
                | --    FINAL
                                - Forcefield selected, all topology files, final conversions

Directories

- INPUT
  - supplied cg file (pdb,gro,tpr)
  - CG converted to pdb (CG_input.pdb)
  - supplied atomistic file (pdb,gro,tpr)
  - supplied atomistic file converted to pdb (AT_input.pdb)
  - script inputs, all flags used in the conversion saved for future reference
- RESIDUE_TYPE (e.g POPE, PROTEIN)
  - individual residues after initial conversion
  - topology for residues
  - mdp for minisation
  - gromacs outputs saved
  - MIN folder containing all individual minimised residues
  - merged pdb containing all minimised residues in a single pdb
- MERGED
  - topology for all residues 
  - mdp for minisation
  - all residue types merged into single pdb 
  - gromacs outputs saved
  - MIN folder containing all merged misiation files
  - ALCHEMBED folder containing all alchembed steps if protein is present
- FINAL
  - FORCEFIELD folder 
  - topology for all residues
  - final atomistic structures in pdb format

<p align="center">
                                   <b>Output conversions</b>
</p>
                                        

The script provides 3 types of coarsegrain conversions. 

You can select which of these is supplied using the flag:
- -o ['all', 'steer', 'align', 'none']

Standard output:

- De novo method. (final_cg2at_de_novo.pdb)
  - fragments are fitted individually to beads and rotated to minimise bond length to connecting beads.
  - if any lipid tails are detected to be trapped within aromatic rings, a short alchembed step is run on each monomer.
  - if no alchembed step is run, a short NVT run is done instead.

Steered output:

- flexible user supplied atomistic fitting. (final_cg2at_steered.pdb)  
  - Atomistic segments are aligned via sequence to the CG protein.
  - The supplied atomistic structures undertake 2 fitting steps.
    - Each segment is rigid body fitted to the CG backbone beads. 
      - Any missing residues or residue types found in the MOD directory are added from the de novo method.
    - Short steered MD step on selected atoms align the AT protein to the CG representation, preserving most interactions (e.g. backbone H-bonds).   
  - The de novo system is then steered to the coordinates from the AT steered protein.  

Aligned output:

- Rigid body user supplied atomistic fitting. (final_cg2at_aligned.pdb)
  - Atomistic segments are aligned via sequence to the CG protein.
  - Each segment is rigid body fitted to the CG backbone beads. 
    - Any missing residues or residue types found in the MOD directory are added from the de novo method.
  - The de novo system is then steered to the coordinates from the rigidly fit protein.  

<p align="center">
                                   <b>**Automation**</b>
</p>
                                        
If you know in advance which settings you wish to use, these can be supplied by the following flags. 

- w   (str) water model         e.g. tip3p
- fg  (str) fragment databases  e.g. martini_2-2_charmm36
- ff  (str) forcefield          e.g. charmm36-jul2017
- box (int) box size (Angstrom) e.g. 100 100 100
- nt  (T/F) neutral N-terminus
- ct  (T/F) neutral C-terminus

example input.

<pre>
    python cg2at.py -c cg_input.pdb -a atomistic_input.pdb -w tip3p -fg martini_2-2_charmm36 -ff charmm36-jul2017-update  -box 100 100 100 -nt -ct
</pre>

<p align="center">
                                   <b>**Database**</b>
</p>                                        

The database has the following file structure and is checked everytime the script is run.

New forcefields and fragments can be added very easily by creating a new folder within the directory structure below. 

<pre>
    | --    database
                | --    scripts_files
                                - run files
                | --    forcefields
                                -  forcefield directories for gromacs (eg. charmm36.ff)
                | --    fragments
                              | -- forcefield type (eg. charmm36)
                                            | -- protein
                                                    | -- Aminoacids (eg. cg residue name ASP)
                                                                 - fragment pdb called the same as bead names (eg. ASP.pdb)
                                                                 - topology file (eg. ASP.top) (optional)
                                                    | -- MOD
                                                          | -- modified residues (eg. cg residue name CYSD) 
                                                                    - fragment pdb called the same as bead names (eg. CYSD.pdb)
                                                                    - topology file (eg. CYSD.top) (optional)            
                                            | -- non_protein                                             
                                                      | --non protein molecules (eg. lipids POPC)
                                                                    - fragment pdb called the same as residue names (eg. POPC.pdb)
                                                                    - itp file of residue called the same as residue names (eg. POPC.itp)
                                                                    - topology file (eg. POPC.top) (optional)
</pre>

You can prevent the script from reading any file or folder by prefixing the name with a underscore.

<p align="center">
                                   <b>**Fragments**</b>
</p>

The fragment database is separated into two parts (protein and non-protein).

Each fragment file contains the following flags:

<pre>
[ BB ]
atom 1
atom 2
...
[ SC1 ]
atom 8
atom 9
...

</pre>

The protein section contains all the normal amino acids and post-translationally modified residues.

The normal amino acids fragments do not contain any hydrogens, these are incorporated by pdb2gmx.

Whilst the modified protein and non protein fragments retain their hydrogens. This is due to problematic adding of every residue to the hydrogen database. 

For ring like structures, I have implemented a grouping system.

In this particular case, the fragments SC2 and SC3 are treated as a rigid body in the fitting and rotation steps.

This default grouping can be switched off using the flag:

- -mod

The grouping is especially useful for converting sugar groups in which the hydrogen geometry should be retained as much as possible.

An example of a normal amino acid fragment files:

<pre>
Phenylalanine

[ BB ]
ATOM      1  N   PHE     1      42.030  16.760  10.920  2.00  0.00           N
ATOM      2  CA  PHE     1      42.770  17.920  11.410  3.00  1.00           C
ATOM     10  C   PHE     1      44.240  17.600  11.550  2.00  0.00           C
ATOM     11  O   PHE     1      44.640  16.530  12.080  6.00  0.00           O
[ SC1 ]
ATOM      3  CB  PHE     1      42.220  18.360  12.800  1.00  1.00           C
ATOM      4  CG  PHE     1      40.730  18.730  12.860  3.00  2.00           C
ATOM      5  CD1 PHE     1      39.780  17.730  13.110  1.00  4.00           C
[ SC2 ]
ATOM      8  CD2 PHE     1      40.300  20.030  12.600  1.00  2.00           C
ATOM      9  CE2 PHE     1      38.940  20.340  12.590  1.00  3.00           C
[ SC3 ]
ATOM      6  CE1 PHE     1      38.420  18.030  13.090  1.00  4.00           C
ATOM      7  CZ  PHE     1      38.000  19.330  12.830  1.00  3.00           C
</pre>

The optional topology file contains the information about grouping and connectivity

<pre>

[ STEER ]
# atoms to steer into CG structure
CA CZ CD1 CD2
[ N_TERMINAL ]
# N_terminal connecting atom, if terminal residue write TER
N
[ C_TERMINAL ]
# C_terminal connecting atom, if terminal residue write TER
C
[ BACKBONE ]
# name of backbone bead, will connect to terminal atoms
BB
[ GROUPS ]
# atoms to put into groups, each line is a separate group
SC2 SC3
[ CHIRAL ]
# column must be in the order:
# central atom, atom to move, atom_1, atom_2, atom_3
# eg: CA HA CB N C
CA HA CB N C

</pre>


For non protein residues you can create a position restraint file which is applied during the creation of the aligned and steered systems.
A script exists within the scripts directory called make_fragments_posre.py this can either create the correct posre files for every residue in the system or for a single residue.

To apply to every folder in directory in non_protein:
<pre>
python make_fragments_posre.py -dir
</pre>
To apply to a single residue:
<pre>
python make_fragments_posre.py -f POPE.itp
</pre>

In the case of solvent. All ions and water molecules are in single repective pdb files with separate groups.

In martini water, 4 atomistic water molecules are condensed into a single bead, therefore the fragment has 4 water molecules.

<pre>
[ tip3p ]
ATOM      1  OW  SOL     1      20.910  21.130  75.300  1.00  0.00
ATOM      2  HW1 SOL     1      20.580  21.660  76.020  1.00  0.00
ATOM      3  HW2 SOL     1      21.640  21.640  74.940  1.00  0.00
ATOM      4  OW  SOL     2      21.000  22.960  77.110  1.00  0.00
ATOM      5  HW1 SOL     2      21.390  23.220  76.270  1.00  0.00
ATOM      6  HW2 SOL     2      21.450  22.160  77.360  1.00  0.00
ATOM      7  OW  SOL     3      22.650  23.020  75.080  1.00  0.00
ATOM      8  HW1 SOL     3      22.830  22.410  75.790  1.00  0.00
ATOM      9  HW2 SOL     3      23.510  23.340  74.810  1.00  0.00
ATOM     10  OW  SOL     4      22.890  21.190  76.980  1.00  0.00
ATOM     11  HW1 SOL     4      22.130  20.970  76.440  1.00  0.00
ATOM     12  HW2 SOL     4      23.090  20.380  77.460  1.00  0.00
[ tip4p ]
ATOM      1  OW  SOL     1      38.680  58.360  49.620  1.00  0.00
ATOM      2  HW1 SOL     1      37.800  58.250  49.280  1.00  0.00
ATOM      3  HW2 SOL     1      38.580  58.980  50.350  1.00  0.00
ATOM      4  MW  SOL     1      38.560  58.430  49.670  1.00  0.00
ATOM      5  OW  SOL     2      36.110  58.430  49.680  1.00  0.00
ATOM      6  HW1 SOL     2      36.530  59.190  50.080  1.00  0.00
ATOM      7  HW2 SOL     2      36.480  57.680  50.140  1.00  0.00
ATOM      8  MW  SOL     2      36.210  58.430  49.790  1.00  0.00
ATOM      9  OW  SOL     3      37.520  59.880  51.350  1.00  0.00
ATOM     10  HW1 SOL     3      37.430  58.950  51.530  1.00  0.00
ATOM     11  HW2 SOL     3      37.480  60.300  52.210  1.00  0.00
ATOM     12  MW  SOL     3      37.500  59.810  51.480  1.00  0.00
ATOM     13  OW  SOL     4      37.410  57.240  51.590  1.00  0.00
ATOM     14  HW1 SOL     4      37.960  57.460  50.840  1.00  0.00
ATOM     15  HW2 SOL     4      37.630  56.330  51.790  1.00  0.00
ATOM     16  MW  SOL     4      37.510  57.150  51.520  1.00  0.00
</pre>

In the case of ions only the ion is stored as a fragment. During the conversion a water bead is superimposed over the ion, to replace the hydration shell of the ion.

<pre>
[ NA+ ]
ATOM      1  NA   NA     1      21.863  22.075  76.118  1.00  0.00
[ NA ]
ATOM      1  NA   NA     1      21.863  22.075  76.118  1.00  0.00
</pre>

<p align="center">
                                   <b>**AT2CG SCRIPT OVERVIEW**</b>
</p>
                                    
Within this script I have included a the reverse of CG2AT. Here the fragment library is used to convert a atomistic system into a coarsegrain representation.

This script roughly follows the following workflow.

- Finds the center of mass of each fragment's heavy atoms as described by the database.
- Writes out a pdb file containing the beads placed at the COM of each fragment.

<pre>
   Atomistic                     CG fragments                            CG beads
                                   --------                              --------    
                                  (        )                            (        )    
                                 (  O1  O2  )                          (          )      
                                (    \  /    )                        (            )   
                               (      CG      )                      (     SC1      )        
    O1  O2                      (     |      )                        (            )    
     \  /                        (    CB    )                          (          )    
      CG                          (        )                            (        )       
      |        Fragments           --------            Coarse            --------  
      CB      ----------->            |             ----------->             |        
      |                            --------             grain            --------   
  X1  CA   X2                     (        )                            (        )  
   \ /  \ /                      (    CA    )                          (          )
    N    C                      (    /  \    )                        (            ) 
         |                  X1-(    N    C    )-X2                X1-(      BB      )-X2
         O                      (        |   )                        (            )     
                                 (       O  )                          (          )
                                  (        )                            (        ) 
                                   --------                              -------- 
</pre> 

This workflow allows each fragment to be treated individually, with no knowledge of what any other bead contains.

<p align="center">
                                   <b>**FLAGS**</b>
</p>
           
The conversion follows a similar syntax to CG2AT.

The following flags are required:

- -c (str)
- -at2cg (True/False)

The following flags are optional
- -fg (str)
- -ff (str)
- -swap (eg POPI:POP2)


example input.

<pre>
python cg2at.py -c at_input.gro -fg martini_2-2_charmm36 -ff martini_2-2 -at2cg 
</pre>

Due the truncation of PI lipid types residue names to 4 characters. The swap group is required to specify to correct lipid type.

example input.

resname POPI (which in reality is POPI2_3-5 not POPI1_3) in martini is called POP2.

<pre>
python cg2at.py -c at_input.gro -fg martini_2-2_charmm36 -ff martini_2-2 -at2cg -swap POPI:POP2
</pre>

<p align="center">
                                   <b>**AT2CG OUTPUT DIRECTORIES**</b>
</p>
                                        
The script will create a output file system as below.

    | --    AT2CG_(timestamp)
                | --    INPUT
                                - at_input.gro, conversion_input.pdb, script_inputs.dat
                | --    FINAL
                                - Forcefield selected, final conversion

Directories

- INPUT
  - supplied AT file (pdb,gro,tpr)
  - AT converted to pdb (conversion_input.pdb)
  - script inputs, all flags used in the conversion saved for future reference
- FINAL
  - FORCEFIELD folder 
  - final CG structures in pdb format
  - topology file
  
A topology will be provided, however it will require updating with the correct protein information. Places to update are highlighted as "XXX" 

To remake your martini topology of the protein. You will need to rerun martinise.py.

This step was not incorporated as it requires the package DSSP to correctly make the topology.

However the martinise.py script is available with the scripts directory. The initial section of the martinise command to run will be printed at the end of the conversion. 
