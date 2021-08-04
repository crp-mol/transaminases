## Predict the ee% for the asymmetric synthesis of amines from ketones catalyzed by Vf-TA
This is a computational protocol to predict the enantioselectivity (ee%) of ω-transaminases using a combination of molecular docking and molecular dynamics simulations. The protocol has been tested in the ω-transaminase from *Vibrio fluvialis* (Vf-TA).

### Requirements
*	[YASARA](http://yasara.org/)
*	[Rosetta suite](https://www.rosettacommons.org/)
*	[Gaussian09](https://gaussian.com/glossary/g09/) (optional)
*	[Avogadro](https://avogadro.cc/) (optional)
*	Linux (preferred) or MacOS (untested, but should work)

## STEP 1: Generate the rotamer library
![rotamer library](imgs/rotamer_library.png)
First, we need to generate rotamer libraries for the external aldimine intermediate form of the query compound (ligand). You need the PDB structure of the query amine (R- or S-amine).
To do so, execute the following YASARA macro:
```bash
$PATH_YASARA/yasara -txt generate_rotamers.mcr
```

The generate_rotamers.mcr script takes the amine form of the query compound as input (PDB or MOL format) and automatically generates rotamer libraries of the external aldimine form. The script should have produce YASARA scene (SCE) files with the rotamer libraries before and after prunning by RMSD and energy criteria. The rotamer library is also contained inside the file named `SUB_rotamers.pdb` to be used in the next stage by Rosetta. If a path for the Rosetta script `molfile_to_params.py` was provided, the `generate_rotamers.mcr` macro should have also generated a Rosetta parameters file (`SUB.params`).

:ledger: **Notes**
*	The user-defined parameters in the `generate_rotamers.mcr` macro are the RMSD cutoff for considering two rotamers as unique (`Prune`), the number of rotamers generated before prunning (`maxNumberRotamers`), and the dihedral angle chi_1 of the external aldimine intermediate library (`DihedralAngleCH`).
*	A chi_1 angle of -90 degrees places the hydrogen-to-be-abstracted (H99) pointing directly toward the catalytic lysine.
*	The `generate_rotamers.mcr` script can also be run in visual mode, allowing the user to see each step taking place. Running the YASARA macro in visual mode can be useful for visual inspection, but for large runs the text-only mode is recommended.
*	To run the YASARA macro in visual mode simply delete the `-txt` flag: `$PATH_YASARA/yasara generate_rotamers.mcr`

**Screenshots of the **`generate_rotamers.mcr`** YASARA macro in visual mode**
![screenshots generate rotamers](imgs/screenshots.png)


## STEP 2: Generate the enzyme-ligand complex (docking)
![Docking](imgs/docking.png)
Using the rotamer library (`SUB_rotamers.pdb`) and the Rosetta params file (`SUB.params`) already generated in the previous step, we can run the [Rosetta Enzyme Design](https://new.rosettacommons.org/docs/latest/application_documentation/design/enzyme-design) application to dock the ligand into the binding site of the enzyme.
First, we align the ligand into the binding site of the enzyme by RMSD fit of the PMP ligand atoms (the atoms of the ligand that originally come from the pyridoxal ring) to the PMP atoms of the original crystal structure (the PMP atoms originally found in e.g., [4E3Q](https://files.rcsb.org/download/4e3q.pdb)). The command is as follows:
```bash
$PATH_YASARA/yasara -txt alignLigand.mcr
```
The YASARA macro should generate `4e3q_forRosetta.pdb`. This is simply the initial positioning of the ligand. You can place the ligand anywhere else you want, but Rosetta will have a harder time finding the binding site than if you simply put the ligand in the binding site.

Once the ligand has been placed in the binding site of the enzyme, you can run Rosetta docking using the following command:
```bash
$ROSETTA_PATH/main/source/bin/enzyme_design.static.linuxgccrelease @flags -resfile resfile -database $ROSETTA_PATH/main/database/ -enzdes::cstfile enzdes.cst -nstruct 10 -s 4e3q_forRosetta.pdb
```
Because the ligand was initially placed aligned to the PMP cofactor (`4e3q_forRosetta.pdb`), the docking algorithm will probably only select among the rotamers provided the one that best fits the binding site, and re-arrange the residues around it. The settings can be found in the `flags` file. You can modify `flags` to adjust the intensity of the search. The `resfile` flags is for mutants. And the `enzdes.cst` contains the user-defined constraints. In `enzdes.cst` we define which is the catalytic lysine, and what distances / angles / dihedrals need to be constrained. The `-nstruct` flag tells Rosetta how many decoys to generate. The `-s` flag is the input file. The `SUB_rotamers.pdb` and `SUB.params` are referenced in the `flags` file.

After running the previous command, Rosetta should have generated 10 decoys (`4e3q_forRosetta__DE_1.pdb`, `4e3q_forRosetta__DE_2.pdb`, ...), and a file called `enz_score.out` that contains a list of Rosetta measurements about the decoys. The most useful one is `SR_2_interf_E_1_3`, which contains the Rosetta Interface Energy between the ligand and the enzyme. Not all the decoys generated are good enough. Sometimes, the ligand ends up outside the binding site (easily seen by the poor Interface Energy). Visual inspection of the genenrated decoys should be done previous to MD simulation (recommended).

:ledger: **Notes**
*	If your query transaminase does not have a crystallized PMP cofactor, then you cannot use the `alignLigand.mcr` YASARA macro. 
*	But, you can align the crystal structure of your transaminase to the crystal structure of a similar transaminase that does contain PMP, and delete the structure of the second transaminase but leave the PMP cofactor. This way, you'll have your original transaminase with PMP (again this step is for the initial placement of the ligand, and can be done in many ways. It's up to the user to decide what's more appropiate for their own specific case). Beware that some ω-transaminases have distinct conformations when bound/unbound to the cofactor ([Sirin et al., 2014](dx.doi.org/10.1021/ci5002185)).
*	If you also want to introduce mutations at the same time, it can easily be done by editing the file called `resfile`. For example:
```
85 B PIKAA ACD
151 A PIKAA ACDEFGHIKLMNOPQRSTVWY
118 A PIKAA NATRO
```

Tells Rosetta to mutate (`PIKAA`) position 85 of subunit B into either Alanine, Cysteine, or Aspartate. `AND` to mutate position 151 of subunit A into either of the 20 cannonical amino acids (including the wild-type amino acid). `AND` to keep the residue 118 of subunit A in its NATural ROtamer form (`NATRO`). For more information visit: https://www.rosettacommons.org/docs/latest/application_documentation/design/enzyme-design 

## STEP 3: Prepare the docked complex for MD simulations.
![Prepare MD](imgs/prepare_MD.png)
The previous step should have generated a docked complex that we can use as starting point for MD simulations. The structure of the docked complex does not contain water molecules, since they were deleted because water molecules are usually not handled well by docking algorithms. Before MD simulations, it is a good idea to put the crystallographic water molecules back into their original position. 
```bash
$PATH_YASARA/yasara -txt ConvertRosettaPdbB2WOW.mcr "Scaffold = 'example_input/4e3q_cleaned'" "MacroTargetA = 'example_input/4e3q_forRosetta__DE_6'" "MacroTargetB = 'example_input/4e3q_forRosetta__DE_6_B__DE_1'"
```
The script produces a YASARA object file (YOB) containing the water molecules from the first file (`example_input/4e3q_cleaned`), the ligand from the second file (`example_input/4e3q_forRosetta__DE_6`), and the enzyme structure from the third file (`example_input/4e3q_forRosetta__DE_6_B__DE_1`). What we are defining is the `Scaffold` (original protein + water), the `MacroTargetA` (enzyme + ligand docked in the binding site A), and the `MacroTargetB` (enzyme + ligand docked in the binding site B)

:ledger: **Notes**
*	Adding back the crystallographic water molecules is a necessity if the molecular dynamics simulations are short (<< 1 ns). Else, you can skip this step.
*	Note that `4e3q_forRosetta__DE_6` was the sixth decoy obtained from **step 2**. This decoy contains the ligand (SUB) docked in the binding site formed mainly by residues of subunit A (binding site A), while the binding site B is empty. To perform MD simulations we have three options:   1) leave binding site B empty, 2) add PMP to binding site B, 3) dock a second ligand (external aldimine) to binding site B. The first alternative is OK, but transaminases are less stable without a PMP cofactor ([Börner et al., 2017](https://doi.org/10.1002/cbic.201700236)). The second alternative is better in this regard, but we might want to make use of the binding site B for actual calculations with the ligand. The third alternative requires you do to do a second docking step (**step 2**) on the binding site B. In the example input files shown, `4e3q_forRosetta__DE_6_B__DE_1` is the result of performing a second docking step on `4e3q_forRosetta__DE_6` this time in the binding site B.

## STEP 4: Perform the MD simulations on the docked structures. 
![MD](imgs/MD.png)
We can perform MD simulations on the docked structures to count the number of reactive poses and calculate the ee% of the enzyme toward the query compound. We just need the YOB file generated in STEP 3 as initial frame.
```bash
$PATH_YASARA/yasara -txt MDSimulation.mcr "MacroTarget ='example_input/4e3q_forRosetta__DE_1_B__DE_1_WOW'" "CurrentSetting ='MultiShort'"
```

The YASARA macro should take care of everything. It will do 5 replicas of 20 ps each, and will count NACs on-the-fly.

:ledger: **Notes**

*	The MD simulation parameters can be controlled with the `CurrentSetting` flag (see `MDSimulation.mcr`). 
*	The file listing the geometric criteria is `NACGeometricCriteria.mcr`. If a file containing the geometric criteria defining reactive poses (NACs) was provided, then the script should also have counted the NACs on-the-fly. 
*	Example output files are provided in `MD/example_output`, for the 5x20 ps setup: 5 replicas of 20 ps each. 
*	In the example provided, the number of NACs counted for subunit A is 3.78 (average across 5 replicas). This number can be found in the file `4e3q_forRosetta__DE_1_B__DE_1_WOW_LSOn_F_2000fs_20000fs_NAC_Results_ZZZ_Combinations_All5Seeds_Summary.tab`.
*	NACs% can be recalculated by using the TAB files that contain the geometric criteria measurements taken across the 20 ps trajectory, e.g. `4e3q_forRosetta__DE_1_B__DE_1_WOW_F_01_A_LSOn00001.tab` is from replica 01 binding site A. The first column is the *time*, the second is whether this frame is a NAC or not (`0.0000` means `no`, `1.0000` means `yes`), the third column is the *criteria_1*, the fourth column is whether *criteria_1* passed the NAC criteria (`0.000` means `no`, `1.0000` means `yes`), the fifth colum is the *criteria_2*, the sixth column is whether *criteria_2* was passed or not, etc...

## STEP 5: Calculate the enantiomeric excess (ee%) 
Once the MD simulations have been run, you can calculate the ee% of the query compound by comparing the %NACs generated by the (R)- and (S)-enantiomer (`NACs_R` and `NACs_S`, respectively).
The formula is: 
```math
ee% = (NACs_S - NACs_R) / (NACs_S + NACs_R) * 100%
```

