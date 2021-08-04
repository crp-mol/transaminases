# Convert the docked complex (enzyme + ligand) to a *yob file (enzyme + ligand + crystal waters)
# MacroTargetA = enzyme + ligand docked in binding site A
# MacroTargetB = enzyme + ligand docked in binding site B
# Command:
>> yasara -txt ConvertRosettaPdbB2WOW.mcr "Scaffold = 'example_input/4e3q_cleaned'" "MacroTargetA = 'example_input/4e3q_forRosetta__DE_6'" "MacroTargetB = 'example_input/4e3q_forRosetta__DE_6_B__DE_1'"
  ^ This command produces enzyme + ligand_subA + ligand_subB + waters as s *yob file
The *yob file can be used as starting frame to run the MD simulations


