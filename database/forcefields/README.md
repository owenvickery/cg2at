Note:
For the amber99sb-ildn forcefield the ion paramters have been updated to those from Joung et al 2008.
This is prevent the ions aggregation at higher salt concentrations.
A copy of this paper can be found within amber99sb-ildn.ff and amber99sb-ildn_slipids.ff folders. 

- amber99sb-ildn.ff         = original parameters
- amber99sb-ildn-tip3p.ff   = updated for tip3p water model
- amber99sb-ildn-spc.ff     = updated for spc/e water model

Updated with stockholm lipids

- amber99sb-ildn_slipids.ff         = original parameters with slipids parameters added
- amber99sb-ildn_slipids_tip3p.ff   = updated for tip3p water model
- amber99sb-ildn_slipids_spc.ff     = updated for spc/e water model


Amber14_OL15.ff downloaded from https://www.gromacs.org/Downloads/User_contributions/Force_fields
Amber OL15 force field for DNA and OL3 (chiOL3) force field for RNA, combined with ff14SB protein force field. More details on ffol.upol.cz. Contributed by Petr Jurecka.
Zgarbova, M.; Sponer, J.; Otyepka, M.; Cheatham, T. E.; Galindo-Murillo, R.; Jurecka, P., Refinement of the Sugar-Phosphate Backbone Torsion Beta for the AMBER Force Fields Improves the Description of Z-DNA and B-DNA. J. Chem. Theory Comput. 2015, 11 (12), 5723-5736


Amber14_OL21.ff downloaded from https://fch.upol.cz/ff_ol/
OL21 is a 2021 DNA force field, so far tested on B-DNA and Z-DNA. 
Zgarbova, M.; Sponer, J.; Jurecka, P., Z-DNA as a Touchstone for Additive Empirical Force Fields and a Refinement of the Alpha/Gamma DNA Torsions for AMBER J. Chem. Theory Comput. 2021, https://doi.org/10.1021/acs.jctc.1c00697


charmm36 forecefields are available:

Base level forcefield with no updates (cannot use with lipids):
charmm36-mar2019.ff
charmm36-jul2020.ff

Base level forcefield with OW,HW added to atomtypes (cannot use with modified side chains):
charmm36-jul2021.ff
charmm36_ljpme-jul2021.ff.  # updated LJ long-range dispersion correction: ljpme, Yu et al., doi: 10.1021/acs.jctc.0c01326

Updated forcefields to contain lipid atom types, non-bonded and bonded parameters.
charmm36-mar2019-updated.ff
charmm36-jul2020-updated.ff

Charmm36 references:
Vanommeslaeghe, K. Hatcher, E. Acharya, C. Kundu, S. Zhong, S. Shim, J. E. Darian, E. Guvench, O. Lopes, P. Vorobyov, I. and MacKerell, Jr. A.D. "CHARMM General Force Field (CGenFF): A force field for drug-like molecules compatible with the CHARMM all-atom additive biological force fields," Journal of Computational Chemistry 31: 671-90, 2010, PMC2888302 
Vanommeslaeghe, K., and MacKerell Jr., A.D., "Automation of the CHARMM General Force Field (CGenFF) I: bond perception and atom typing," Journal of Chemical Informationa and Modeling, 52: 3144-3154, 2012, PMC3528824 
Vanommeslaeghe, K., Raman, E.P., and MacKerell Jr., A.D., "Automation of the CHARMM General Force Field (CGenFF) II: Assignment of bonded parameters and partial atomic charges, Journal of Chemical Informationa and Modeling, 52: 3155-3168, 2012, PMC3528813 
Yu, W., He, X., Vanommeslaeghe, K. and MacKerell, A.D., Jr., "Extension of the CHARMM General Force Field to Sulfonyl-Containing Compounds and Its Utility in Biomolecular Simulations," Journal of Computational Chemistry, 33: 2451-2468, 2012, PMC3477297 
Soteras Gutierrez, I., Lin, F.-Y., Vanommeslaeghe, K., Lemkul, J.A., Armacost, K.A., Brooks, Cl., III, and MacKerell, A.D., Jr., "Parametrization of Halogen Bonds in the CHARMM General Force Field: Improved treatment of ligand-protein interactions," Bioorganic & Medicinal Chemistry, In Press, 2016,
