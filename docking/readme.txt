# Step 1: Put the ligand (01S, 02S, 03S.... etc) into the binding site of Vf-TA. You can do it manually or by running the following YASARA macro:
yasara -txt alignLigand.mcr

# Step 2: Rosetta docking
# This command will do docking of the external aldimine intermediate into the binding site of the enzyme
/home/carlos/soft/rosetta_bin_linux_2015.25.57927_bundle/main/source/bin/enzyme_design.static.linuxgccrelease @flags -resfile resfile -database /home/carlos/soft/rosetta_bin_linux_2015.25.57927_bundle/main/database/ -enzdes::cstfile enzdes.cst -nstruct 10 -s 4e3q_forRosetta.pdb > log

# Step 3: Select the best docked structures (decoys) for MD simulations
# The best decoys are selected based on the Rosetta Interface Energy, which can be found in the enz_score.out file
awk '{print $29}' enz_score.out

# Note: Don't forget to add the following: "PDB_ROTAMERS SUB_rotamers.pdb" as the last line of the SUB.params file. This line should have been added automatically by the create_rotamers.mcr yasara script already. But, in case you did not provide a valid Rosetta path, you'll have to add it manually.
