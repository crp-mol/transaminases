# Run a YASARA MD simulation of the docked complex (VfTA + external aldimine intermediate)
# The MD simulations can be 5x20ps (5 replicas of 20 ps simulation) or 1x100ps or any other value. 
# We recommend you keep the simulations short because Rosetta will have already made the larger jumps in conformational space

# Input: *.yob file
# Output: *tab files containing the NAC count

# Step 1: Run the MD simulation using a *yob file as input. The *yob file contains ligand + enzyme + crystallographic waters
yasara -txt MDSimulation.mcr "MacroTarget ='example_input/4e3q_forRosetta__DE_1_B__DE_1_WOW'" "CurrentSetting ='MultiShort'" >> logMD

# Step 2: Count NACs
  ## Option A: If NACs definitions were provided for MDSimulation.mcr to use, then the %NACs of all replicas is already calculated and can be found in example_output/*NAC_Results_ZZZ_Combinations_All5Seeds_Summary.tab
  ## Option B: If you changed your mind about the NAC definitions you provided to MDSimulation.mcr (from NACGeometricCriteria.mcr), then you can easily find all measurements in the output *tab files. Re-calculating NACS can be done very easily. For example:
  >> cat 4e3q_forRosetta__DE_1_B__DE_1_WOW_F_??_A_LSOn00001.tab | grep -v "Mixed_Data" | sed '/^$/d' | awk '($3<3.5)&&($5<20)&&($5>-20) {print $3 " " $5}' | wc -l
     ^ This command line will count the number of frames in which the distance d_1 is closer than 3.5 Agstrom and the angle chi_2 is between 20 and -20 deg.
  ## Option C: If you did not provide any NAC definitions to MDSimulation.mcr to use, then there should be trajectory output frames available (*sim files). You have to create your own YASARA script to read the trajectory and count the frames that have a NAC.

# Step 3: If you're following the protocol provide in the manuscript, then you should average the number of NACs you obtain from this run (ligand 01S, library with chi_1 = -90) with the NACs obtained from all the other libraries, e.g. chi_1 = {-168.75 ,-157.50 ,-146.25 ,-135.00 ,-123.75 ,-112.50 ,-101.25 ,-90.00 ,-78.75 ,-67.50 ,-56.25 ,-45.00 ,-33.75 ,-22.50 ,-11.25 , +0.00 , +11.25 , +22.50 , +33.75 , +45.00 , +56.25 , +67.50 , +78.75 , +90.00 , +101.25 , +112.50 , +123.75 , +135.00 , +146.25 , +157.50 , +168.75 , +180.00}
This will give you the number of NACs that the ligand 01S produced.

# Step 4: Calculate ee%
The ee% can be calculated by comparing the number of NACs that the ligand 01S produced with the number of NACs that the ligand 01R produced.
The formula is: ee% = (NACs_S - NACs_R) / (NACs_S + NACs_R) * 100%

