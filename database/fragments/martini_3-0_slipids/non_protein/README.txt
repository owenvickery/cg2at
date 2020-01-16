Slipids force field for various lipids
---------------------------------------

Authors:  Joakim Jämbeck, Inna Ermilova, Alexander Lyubartsev
	  Department of Materials and Environmental Chemistry,
	  Stockholm University,  Stockholm   10691  Sweden
	  e-mail:  alexander.lyubartsev@mmk.su.se


See also: http://www.fos.su.se/~sasha/SLipids/

Content:

SLipids_FF:  directory containing the force field. Included into the Gromacs
	     topology file by:  
	     #include "SLipids_FF/forcefield.itp"

itp_files:   itp files for various lipids



The force field can be used togther with the AMBER99SB/AMBER99SB-ILDN/AMBER03
FF for proteins


!!!! MAKE SURE YOU CITE THE FOLLOWING REFERENCES WHEN USING THIS FORCE FIELD !


Saturated PC lipids:

Joakim P. M. Jämbeck and Alexander P. Lyubartsev, "Derivation and Systematic
Validation of a Refined All-Atom Force Field for Phosphatidylcholine Lipids",
J. Phys. Chem. B, 2012, 116, 3164-3179 (2012)  DOI: 10.1021/jp212503e

POPC, DOPC, SOPC, DOPE, POPE and similar:
 
Joakim P. M. Jämbeck and Alexander P. Lyubartsev, "An Extension and Further
Validation of an All-Atomistic Force Field for Biological Membranes"
J. Chem. Theory Comput., 8, 2938-2948, (2012) DOI: 10.1021/ct300342n

PS, PG, SM lipids and Cholesterol:

Joakim P. M. Jämbeck and Alexander P. Lyubartsev, "Another Piece of the
Membrane Puzzle: Extending Slipids Further"
J. Chem. Theory Comput.,  9 (1), 774-784 (2013) DOI: 10.1021/ct300777p

Polyinsaturated lipids:

Inna Ermilova and Alexander Lyubartsev:, "Extension of the Slipids Force
Field for Polyunsaturated Lipids",
J. Phys. Chem. B, 120 (50), 12826–12842 (2016)


Cite also this paper on Charmm36 force field since bond and angle parameters,
as well as a part of Lennard-Jones parameters and torsion angles in the lipid
headgroups in SLipids are taken from the Charmm36 force field:

Jeffery B. Klauda, Richard M. Venable, Alfredo Freites, Joseph W. O’Connor,
Douglas J. Tobias, Carlos Mondragon-Ramirez, Igor Vorobyov, Alexander D.
MacKerell, Jr. and Richard W. Pastor, "Update of the CHARMM all-atom additive
force field for lipids: Validation on six lipid types", J.Phys.Chem B, 114,
7830–78 (2010)
