# Generate rotamer libraries
# Input: PDB of the query amine
# Output: Rotamer libraries of the external aldimine intermediate form of the query amine
yasara -txt generate_rotamers.mcr 

# Notes
-  ListMacroTargets() = '01R','01S','06R','06S' # List of input amines. PDB file should be present: 01R.pdb, 01S.pdb, etc...
-  Cofactor ='PMP' # Requires PDB file, PMP.pdb
-  Prune, MaxEaboveMin, and maxNumberRotamers settings can be altered by the user.
-  DihedralAngleCH() = '-90.00' # Dihedral angle chi_1.  -90 means the atom H99 is pointing directly towards Lys285-NH2
-  molfiletoparams # Points to the molfile_to_params.py ROSETTA file to generate the SUB.params file (to be used for docking)
