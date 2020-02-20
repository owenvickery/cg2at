                                        **Database**

Once upon a time in the database. There were 2 folders, one liked to know how the world worked, therefore collected Molecular Dynamics forcefields, whilst the other was more material minded and hoarded all the atoms and kept them locked away in fragments. Only to see the world when the folder had a martini in hand.


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
                                                                 - chiral file called chiral.dat (optional)
                                                    | -- MOD
                                                          | -- modified residues (eg. cg residue name CYSD) 
                                                                    - fragment pdb called the same as bead names (eg. CYSD.pdb)
                                                                    - chiral file called chiral.dat (optional)            
                                            | -- non_protein                                             
                                                      | --non protein molecules (eg. lipids POPC)
                                                                    - fragment pdb called the same as residue names (eg. POPC.pdb)
                                                                    - itp file of residue called the same as residue names (eg. POPC.itp)
                                                                    - chiral file called chiral.dat (optional)

You can prevent the script from reading any file or folder by prefixing the name with a underscore.

                                        **Fragments**

The fragment database is separated into two parts (protein and non-protein).

Each fragment file contains the following flags:

<pre>
REQUIRED ALL

- frag 'fragment name'
- group 'group number'

REQUIRED PROTEIN ONLY

 - posres 'atoms to position restrain (separated by a comma)'
 - N_ter  'atom connecting to preceeding aminoacid'
 - C_ter  'atom connecting to proceeding aminoacid'

REQUIRED MODIFIED TERMINAL AMINOACID
if the modified residue has a non standard terminus eg CYST
 - ter  true/false 
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

[ frag 'BB'  group '1' posres 'CA,CZ' N_ter 'N' C_ter 'C' ]
ATOM      1  N   PHE     1      42.030  16.760  10.920  2.00  0.00           N
ATOM      2  CA  PHE     1      42.770  17.920  11.410  3.00  1.00           C
ATOM     10  C   PHE     1      44.240  17.600  11.550  2.00  0.00           C
ATOM     11  O   PHE     1      44.640  16.530  12.080  6.00  0.00           O
[ frag 'SC1'  group '2' ]
ATOM      3  CB  PHE     1      42.220  18.360  12.800  1.00  1.00           C
ATOM      4  CG  PHE     1      40.730  18.730  12.860  3.00  2.00           C
ATOM      5  CD1 PHE     1      39.780  17.730  13.110  1.00  4.00           C
[ frag 'SC2'  group '3' ]
ATOM      8  CD2 PHE     1      40.300  20.030  12.600  1.00  2.00           C
ATOM      9  CE2 PHE     1      38.940  20.340  12.590  1.00  3.00           C
[ frag 'SC3'  group '3' ]
ATOM      6  CE1 PHE     1      38.420  18.030  13.090  1.00  4.00           C
ATOM      7  CZ  PHE     1      38.000  19.330  12.830  1.00  3.00           C
</pre>

In the case of solvent. All ions and water molecules are in single repective pdb files with separate groups.

In martini water, 4 atomistic water molecules are condensed into a single bead, therefore the fragment has 4 water molecules.

<pre>
[ frag 'tip3p' group '1' ]
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
[ frag 'tip4p' group '2' ]
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
[ frag 'NA+'  group '1' ]
ATOM      1  NA   NA     1      21.863  22.075  76.118  1.00  0.00
[ frag 'NA'  group '2' ]
ATOM      1  NA   NA     1      21.863  22.075  76.118  1.00  0.00
</pre>