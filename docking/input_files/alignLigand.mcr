# This YASARA macro aligns the external aldimine intermediate (ligand) to the PMP atoms of the crystal structure (4E3Q) of Vf-TA
LOADPDB 4e3q_cleaned.pdb, center=no # Crystal structure without water molecules and with optimized hydrogen bonding network
LOADPDB SUB_A_0001.pdb, center=no
SUPORDEREDATOM OBJ 2 ATOM N50, OBJ 1 MOL A RES PMP ATOM N1, OBJ 2 ATOM C3, OBJ 1 MOL A RES PMP ATOM C3, OBJ 2 ATOM C50, OBJ 1 MOL A RES PMP ATOM C4A, OBJ 2 ATOM C6, OBJ 1 MOL A RES PMP ATOM C5
SAVEPDB 2, SUB_A_0001_alignedWithA.pdb
SHELL cp Scaffold.pdb 4e3q_forRosetta.pdb
SHELL cat SUB_A_0001_alignedWithA.pdb | grep -e "END" -e "HETATM" | grep -v "REMARK" | sed 's/1.00  0.00/1.00 20.00/g' >> 4e3q_forRosetta.pdb
SHELL rm SUB_A_0001_alignedWithA.pdb
EXIT
