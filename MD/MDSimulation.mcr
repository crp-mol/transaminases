###########################################################################################################################################
#  Name Script: 	MDSimulation.mcr
#               	------------------
#
#  Author:      	Hein J. Wijma, University of Groningen, The Netherlands, Biochemical Laboratory, Biotransformation and Biocatalysis group
#                       H.J.Wijma@rug.nl
#                       Parts of this script were inspired by Elmar Krieger his macros. 
# 
#  Purpose:             - Can be used to create a series of differently initialized trajectories of which the relevant angles and distance are sampled on the fly, during 
#                         the MD simulation.
#                       - script written such that it can be easily adapted for a new enzyme
#
#  Requirements: 	- YASARA-Structure 
#                	- yob file with correct protonation states of the entire protein and of the ligand, this is the MacroTarget. 
#                          
#  Results:     	- MD trajectory(ies) for the (differently seeded) simulations(s) that can be inspected with Elmar Krieger's md_play macro
#               	- Several Table Files with Near Attack Conformation percentages and other statistical data, such as for which NAC criteria the requirements are met
#                       - overview flexibility of the protein versus flexibility of the X-ray structure (_RMSF.tab)
#                       - PDB file with average structure from the MD simulation
#                       - table files with energy and RMSD from starting structure versus simulation time
#  
#  User modify: 	- if no definitions are available for your enzyme, you will have to insert suitable definition for the Near Attack Conformation (NAC) in your enzyme.
#                       - suitable definitions consist of (see below for more details)
#                            	- name of the criterion (12 characters, e.g. '  AngletoNH1', spaces are only possible at the start)
#                               - how to measure the criterion
#                               - minimal requirement, should always be defined even if the smaller the better (e.g. define then -0.01 Angstrom for the minimal distance)
#                               - maximally allowed, should also always be defined, even if the larger the better (e.g. define an angle of 1000 degrees as maximum).
# 
#  Important:   	- The produced average PDB files standard have RMSF instead of b-factors. This is done since with RMSF, with its unit in Angstrom, the flexibility is 
#                         more easy to grasp than with b-factors. This can be changed back to b-factors if desired as described below under UseRMSFInPdbYOBFile
# 	
#  Further:     	Under linux this script is easiest printed with the following or an equivalent command: mpage -l2 -W180 -H -m50 -r DMDAnalysisNacs.mcr | lpr
#

ONERROR EXIT

# --------------------------------------------------------------------------------------------------------------------------- #
# ----------------------- The four settings below have to altered by the users regulary ------------------------------------- # 

# Define whether to automatically exit at end, YES if you are running macros under -txt mode
# Set here whether we want Yasara to exit after finishing the script by choosing for 'YES' or 'NO'
AutomaticExitAtEnd = 'YES'

# Set the number of processors that Yasara is allowed to use. 
PROCESSORS 1

# Define Current Target of the script, the script needs to know which enzyme your are working with
# What it means is defined in: NACGeometricCriteria.mcr
CT = 'AT_vibriofluvialis'

# Set the mode (see below for what these settings of ForDebuggingDefinitions, MultiShort (good sampling of conformations at low cpu cost)
# You can also define your own settings
#CurrentSetting = 'MultiShort'

# At the end of the script, delete all *tab files (except *ZZZ* and criteria_vs_time)
DeleteUnnecessaryFiles='YES'

# --------------------------------------------------------------------------------------------------------------------------- #
# -- The part below does not need to be adapted unless a new enzyme (CT) is adopted or unless a new MD protocol is desired -- #

# Parameters 1: Simulation time periods for the phases of warming, equilibration, and production 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if (CurrentSetting == 'ForDebuggingDefinitions')
  WarmUpRegime           = 'F'
  EquilibrationTime      = 1000
  ProductionTime         = 3000
  MS                     = 2
  InterValAnalysis       = 'On'
  TrajectorySubAnalysisIntervals() =1000,2000
  UseRMSFInPdbYOBFile    = 'TRUE'
  MoreLOGFile            = 'TRUE'
  SnapShotInterval       = 1000
  NACSamplingInterval    = 20
  NACShotEveryXSnapshots = 1
  AssayTemperature       = 298
  EmployedForcefield     =  'Yamber3'
  LincsSettle            = 'Off'
  CalculateAverageTrajectory = 'Off'
  SkipDifferentMinimization = 'TRUE'
  ApplyConstraintsOrCharges = 'Yes'
  ConstraintsChargesName    = 'CorrectCharges'

if (CurrentSetting == 'MultiShort')
  WarmUpRegime           = 'F'         # normal (N) or fast (F) warming
  EquilibrationTime      = 2000        # in fs, 2 ps
  ProductionTime         = 20000       # in fs, 20 ps
  MS                     = 5           # Number of replicas
  InterValAnalysis       = 'Off'       # distribution of NACs over time intervals (not necessary).
  UseRMSFInPdbYOBFile    = 'TRUE'
  MoreLOGFile            = 'UnTRUE'
  SnapShotInterval       = 20000       # take snapshots every dt (in fs). Saved as *sim file
  NACSamplingInterval    = 20          # How often to count the NACs (in fs). It's output to a *tab file
  NACShotEveryXSnapshots = 1           # How often to update the output file.
  AssayTemperature       = 298         # in Kelvin
  EmployedForcefield     =  'Yamber3'  
  LincsSettle            = 'On'        # Constraints
  CalculateAverageTrajectory = 'Off'   # whether to print out a PDB file containing an average trajectory
  SkipDifferentMinimization = 'TRUE'
  ApplyConstraintsOrCharges = 'No'     # It was yes, but commented out.
  ConstraintsChargesName    = 'CorrectCharges'

if (CurrentSetting == 'slightlyLonger')
  WarmUpRegime           = 'N'
  EquilibrationTime      = 20000
  ProductionTime         = 1000000
  MS                     = 1
  InterValAnalysis       = 'Off'
  UseRMSFInPdbYOBFile    = 'TRUE'
  MoreLOGFile            = 'UnTRUE'
  SnapShotInterval       = 2000
  NACSamplingInterval    = 100
  NACShotEveryXSnapshots = 1
  AssayTemperature       = 298
  EmployedForcefield     =  'Yamber3'
  LincsSettle            = 'On'
  CalculateAverageTrajectory = 'Off'
  SkipDifferentMinimization = 'TRUE'
  ApplyConstraintsOrCharges = 'Yes'
  ConstraintsChargesName    = 'CorrectCharges'


# This is the easiest way to get the ligand name distributed. Just need to make sure a suitable file is present. 
LigandName = '(MacroTarget)'

# Count the number of NACs on-the-fly with predefined geometric criteria
include NACGeometricCriteria.mcr

# These 20 lines are commented out, uncomment if you want RESP charges (recommended but not necessary). 
#lkslkdef CorrectCharges LigandName
#  # get the charges
#  # =============================================
#  NameLigand = 'SUB'
# for LetterSubunit in 'A','B'
#    LOADTAB (LigandName)_RESP_names(LetterSubunit).tab
#  
#    # convert them to proper lists
#    RawNewCharges() = TAB 1
#    DELTAB 1
#    NumberOfEntries = (count RawNewCharges) / 4
#    for z = 1 to NumberOfEntries
#      AtomNumbersRESP(z) = (RawNewCharges((((z) -1)*4) + 1))
#      AtomNamefrmRESP(z) = '(RawNewCharges((((z) -1)*4) + 2))'
#      AtomElementRESP(z) = '(RawNewCharges((((z) -1)*4) + 3))'
#      AtomChargefRESP(z) = (RawNewCharges((((z) -1)*4) + 4))
#      #PRINT just read in memory (0+(AtomNumbersRESP(z))) (AtomNamefrmRESP(z)) (AtomElementRESP(z)) (0.0000 + (AtomChargefRESP(z)))
#  
#    # The following assumes the right atom names are given in the structure. The code below which is commented out was for the case of single ligand with different atom names in simulation.
#    for z = 1 to NumberOfEntries
#      CHARGEATOM MOL (LetterSubunit) res (NameLigand) atom (AtomNamefrmRESP(z)), (AtomChargefRESP(z))
#  
#    PRINT just corrected charges, may not work when doing an neutralization experiment
  
# This verifies that some needed files are present and no setting are wrong.Some other files and settings are checked elsewhere
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 

if MacroTarget == ''
  RaiseError The Enzyme Substrate complex Target is missing
if CT == ''
  RaiseError The Current Target, which enzyme whe are working with, is missing
if WarmUpRegime != 'F'
  if WarmUpRegime != 'N'
    RaiseError Illegal Warm up regime defined in script, (WarmUpRegime) should be either 'N' or 'F' for normal or fast warming


# ========================================================================================================================================================================================
# =========================  Part B of the Script: Neutralization with salt and energy minimisation ======================================================================================
# ========================================================================================================================================================================================
# This removes anything left in memory that can hinder the current calculations.
CLEAR
#  Assumes the protonation states of the His/Glu/Asp residues have been set manually. 
#  What the procedure does is that it uses the existing Yasara Neutralization procedure to place salt ions. Then it reverts to the original structure
#  (with the manually set ionization states), and places waters with the command FillCellWater (does not alter ionization states) and puts all salt ions that are 
#  acceptable inside. After that the cell is unlikely to be neutral, and it positions additional salt ions to neutralize the cell again. 
#

if MoreLOGFile == 'TRUE'
  CONSOLE On
else
  CONSOLE Off


UnMinimizedStructureExists = FILESIZE (MacroTarget)_SaltWaterBoxNotMinimized.sce
if !(UnMinimizedStructureExists)
  FileThere = FILESIZE (MacroTarget).yob
  if (FileThere) 
    LOADYOB (MacroTarget)
  else 
    FileThere = FILESIZE (MacroTarget).pdb
    if (FileThere)
      LOADPDB (MacroTarget)
    else
      RAISEERROR No MacroTarget Found
  # some settings
  NICEORIOBJ 1
  FORCEFIELD YAMBER3, SETPAR = YES
  #CELL AUTO, extension = 7.5, Scale = Yes
  CELL Auto, Extension=7.5,Shape=Dodecahedron,res protein 
  
  # make a copy of the protein and remove it temporarily
  DUPLICATEOBJ 1
  RENUMBEROBJ 3,4
  REMOVEOBJ 4
  
  # Get salt ions vio the normal procedure
  FORCEFIELD (EmployedForcefield)
  EXPERIMENT Neutralization
    WaterDensity 0.997
    pH (pH)
    NaCl 0.5
    pKaFile (MacroTarget).pka
    Speed Fast
  FIXATOM water 
  EXPERIMENT On
  # stuff below seems not to work because it is a Neutralization experiment
  if ApplyConstraintsOrCharges == 'Yes'
    (ConstraintsChargesName) '(LigandName)'
  WAIT ExpEnd
  
  # now keep the salt ions, but delete the protein since the ionizations states have altered now. The macro requires the user to provide a pdb file with the correct ionization states.
  DELOBJ 1
  DELRES element O
  ADDOBJ 4
  RENUMBEROBJ 4,1
  RENUMBEROBJ 3,4
  REMOVEOBJ 4
  
  # use command FillCellWater since with original neutralization too many waters placed inside the protein. Then a protocol based on that used by Elmar Krieger
  # -----------------------------------------------------------------------------------------------------------------------------------------------------------  
  LONGRANGE none
  FILLCELLWATER Density=0.997,Probe=1.4,BumpSum=1.0,DisMax=0
  # Fix all the heavy atoms of object 1, the original
  FIXATOM OBJ 1
  FREERES HOH

  
  # some settings of pressure and fast electrostatics
  PRESSURECTRL Off 
  SIMSPEED Fast
  TEMPERATURE (AssayTemperature)
  # switch electrostatics off to prevent local minima
  INTERACTIONS Bond,Angle,Dihedral,Planarity,VdW
  TIMESTEP 2,2.00
  # now 25 timesteps of steepes descent minimization of OBj 3 water and hydrogen atoms
  TempCtrl SteepDes
  FIXBOND water, water
  SIM ON
  if ApplyConstraintsOrCharges == 'Yes'
    (ConstraintsChargesName) '(LigandName)'
  WAIT 25
  # switch electrostatics back on and do a short simulated annealing, too long will create vacuum bubbles
  INTERACTIONS Bond,Angle,Dihedral,Planarity,Coulomb,VdW
  TEMPCTRL Anneal
  if ApplyConstraintsOrCharges == 'Yes'
    (ConstraintsChargesName) '(LigandName)'
  WAIT 100
  TEMPERATURE (AssayTemperature)  
  # Followed by 100 steps normal dynamics to ensure the water is OK.   
  Wait 100
  FREERES ALL
  Timestep 2,1.25
  SIM OFF
  
  
  # Now add the previously added salt ions and remove waters that are too close to the salt or salt that is too close to the active site. 
  ADDOBJ 4
  RENUMBEROBJ 3,5
  JOINOBJ 5,4
  RENUMBEROBJ 4,3
  DELRES OBJ 3 res element Na Cl with distance < (StayAwayDistance) from res (ActiveSiteResidues)
  DELRES OBJ 3 res element O with distance < 1.75 from res OBJ 3 atom element Cl
  DELRES OBJ 3 res element O with distance < 1.02 from res OBJ 3 atom element Na
  
  
  # some settings of pressure and fast electrostatics
  PRESSURECTRL Type=Combined,Pressure=1.000,Name=HOH,Density=0.997,Axis=XYZ
  TIMESTEP 2,1.25
  LONGRANGE Coulomb
  BOUNDARY periodic
  ENERGYUNIT kJ/mol
  
  # now need an algorithms to check every water from OBj 3 that lives more than 4 angrstrom from the protein, 
  # find the one at the most positive/negative place place, change it, relist all waters, find again the most pos/negative, untill neutral. 
  # present them by Ballres 
  
  # Check what the charge of the system is before neutralization
  CurrentChargeList() = CHARGEOBJ all 
  CurrentCharge = sum CurrentChargeList
  CurrentCharge = 0 + (CurrentCharge)
  
  # determine what kind of ions, and how many, to introduce for neutralization
  ReplaceWithSodiumIons = 0
  ReplaceWithChlorideIons = 0
  if (CurrentCharge < 0.5)
    ReplaceWithSodiumIons = - (CurrentCharge)
  if (CurrentCharge > 0.5)  
    ReplaceWithChlorideIons = (CurrentCharge)
  
  # now do the replacements, start with the ion for neutralization only after those are finished add the ones for salt.
  If ((ReplaceWithSodiumIons) > (ReplaceWithChlorideIons))
    FirstIon = 'SOD'
    SecondIon = 'CHL'
    FirstChargeCount  = 0 + (ReplaceWithSodiumIons)
    SecondChargeCount =  0 + (ReplaceWithChlorideIons)
  else
    FirstIon = 'CHL'
    SecondIon = 'SOD'  
    FirstChargeCount  = 0 + (ReplaceWithChlorideIons)
    SecondChargeCount =  0 + (ReplaceWithSodiumIons)
  
  if (FirstChargeCount)
    HIDERES HOH
    REMOVEENVRES HOH
    for CurrentStage in 'FirstRound', 'SecondRound'
    if CurrentStage == 'FirstRound'
      CurrentIon = 'FirstIon'
      ChargeCount = (FirstChargeCount)
    if CurrentStage == 'SecondRound'
      CurrentIon = 'SecondIon'
      ChargeCount = (SecondChargeCount)
    for i = 1 to ChargeCount
      # reset the lists 
      EligibleH2OList() = 0
      ESPSURFLIST()     = 0  
      # get the list of waters that are part of Object 3 and more than 6 angstrom away from the , the list consists of atom numbers
      SIM on
      SIM pause
      EligibleH2OList() = LISTRES OBJ 3 res HOH with distance > 6 from res SOD CHL CIM CIP
      SIM Off
      # determine the Electrostatic potential at their surface (the list consists of surface and potential energy assuming the atom charge = +1)
      ESPSURFLIST() = SURFESPRES OBJ 3 res HOH with distance > 6 from res SOD CHL, type = accessible, method = PME, unit = RES
      # seems bug in yasara in which the first value for ESPATOM is always nan
      ESPSURFLIST(2) = 0
      for j = 1 to count EligibleH2OList
        CurrentCoulomb = ESPSURFLIST(j*2)
        if j==1
          MostNegativeSpot  = (CurrentCoulomb)
          MostPositiveSpot  = (CurrentCoulomb)
          IndexMostNegative = 1
          IndexMostPositive = 1
        if CurrentCoulomb<MostNegativeSpot
          CriticalDistance = GroupDistance res (ActiveSiteResidues),  atom (EligibleH2OList(j))
          if (CriticalDistance > StayAwayDistance)
            MostNegativeSpot  = (CurrentCoulomb)
            IndexMostNegative = (j)
        if CurrentCoulomb>MostPositiveSpot
          CriticalDistance = GroupDistance res (ActiveSiteResidues),  atom (EligibleH2OList(j))
          if (CriticalDistance > StayAwayDistance)
            MostPositiveSpot  = (CurrentCoulomb)
            IndexMostPositive = (j)
      # do the replacements, do not forget to delete the hydrogen atoms that belong to it (be carefull the atom numbers might shift)
      if (ReplaceWithSodiumIons) > 0.5
        if ((CurrentIon) == 'SOD')
          RENAMERES atom (EligibleH2OList(IndexMostNegative)), SOD
          ReplaceWithSodiumIons = (ReplaceWithSodiumIons) - 1
      if (ReplaceWithChlorideIons) > 0.5
        if ((CurrentIon) == 'CHL')
          RENAMERES atom (EligibleH2OList(IndexMostPositive)), CHL
          ReplaceWithChlorideIons = (ReplaceWithChlorideIons) - 1
      DELATOM RES CHL SOD element H
      SWAPATOM RES SOD, Na
      SWAPATOM RES CHL, Cl
      BALLRES RES SOD CHL
      ADDENVRES SOD CHL
      SIM init

  SIM OFF
  DELRES OBJ 3 Res HOH with distance < 2.8 from Res Sub

  SAVESCE (MacroTarget)_SaltWaterBoxNotMinimized

# also check the conformations in the Unminimized structure, if that has not yet been done

LOADSCE (MacroTarget)_SaltWaterBoxNotMinimized
if ApplyConstraintsOrCharges == 'Yes'
  (ConstraintsChargesName) '(LigandName)'

  
for i = 1 to count L
  b = '(L(i))'
  for c = 1 to count Label_(b)_
    (Label_(b)_(c))_Q_       = ((Label_(b)_(c))_Measurement) 
  # Magic number, i.e. NACs, and all criteria with their NAC
  MagicNumber = 1 + (2*(count Label_(b)_))
  # determine whether NAC is achieved 
  NAC = 1
  for c = 1 to count Label_(b)_
    (Label_(b)_(c))_NAC = 0.00000
    if (  ((Label_(b)_(c))_Q_)  > ((Label_(b)_(c))_Q_Min))
      if (((Label_(b)_(c))_Q_)  < ((Label_(b)_(c))_Q_Max))
        (Label_(b)_(c))_NAC = 1.00000
    NAC = (NAC) * ((Label_(b)_(c))_NAC)  
  TABULATE (100.000 * (NAC))
  for c = 1 to count Label_(b)_
    TABULATE ((Label_(b)_(c))_Q_)
    TABULATE ((Label_(b)_(c))_NAC)
  
  LabelCollection = ' ________NAC?'
  for c = 1 to count Label_(b)_
    LabelCollection = '(LabelCollection) (Label_(b)_(c)) CriteriaMet?'
  SaveTab 1,(MacroTarget)_NAC_Results_(L(i))_UnMinimized,Format=Text,Columns=(MagicNumber),NumFormat=%12.2f,(LabelCollection)
  DELTAB 1



  
# ========================================================================================================================================================================================
# =========================  Part C of the Script: The MD simulation(s) ==================================================================================================================
# ========================================================================================================================================================================================


#  **********************************  HERE STARTS THE LOOP THAT DOES THE ENTIRE MD SIMULATION   *****************************************************************************************
#  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 

# This removes anything left in memory that can hinder the current calculations.
CLEAR

MaximumSeeds = 00 + (MS)
for CSN = 01 to (MaximumSeeds)
  # for multi trajectories
  SeedingNumber = (0.001234567*(CSN))
  
  # The warmuptime is 
  #  do not go under 3000 fs since otherwise it goes so fast that the simulation is out of control by the time room temperature is reached
  if WarmUpRegime == 'F'
    WarmUpTime = 3000
  if WarmUpRegime == 'N'
    WarmUpTime = 30000
  

  # Settings of the forcefield and temperature pressure control, saving snapshots, etcetera
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  
  # Use the above defined ForceField
  FORCEFIELD (EmployedForcefield ),SETPAR = Yes
  
  # Ensure no fixed atoms are present present but quickly allow for constraints
  FREE
  
  
  # use periodic boundary conditions
  BOUNDARY Periodic

  # Use the earlier made saltwater box not minimized structure, and the minimized _water.sce if it exists, else do an independently seeded minimization
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  
  # ensure some files exists, to prevent problems later
  UnMinimizedStructureExists = FILESIZE (MacroTarget)_SaltWaterBoxNotMinimized.sce
  if (UnMinimizedStructureExists)
    MinimizedStructureExists = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)water.sce
    if (MinimizedStructureExists)
      LOADSCE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)water.sce
      if ApplyConstraintsOrCharges == 'Yes'
        (ConstraintsChargesName) '(LigandName)'
    else
      if (SkipDifferentMinimization != 'TRUE')
        LOADSCE (MacroTarget)_SaltWaterBoxNotMinimized
        if ApplyConstraintsOrCharges == 'Yes'
          (ConstraintsChargesName) '(LigandName)'
        # now a trick to get differently initilized minimizations, by putting a water in front
        JOINRES OBJ 3
        for i = 1 to (CSN) 
          AllWaters()    = LISTRES OBJ 3 res HOH element O
          NumberOfWaters = COUNTRES OBJ 3 res HOH element O
          SPLITRES OBJ 3 res (AllWaters(NumberOfWaters)) 
          SPLITOBJ 3
          JOINOBJ 1,4
          RENUMBEROBJ 4,1     
        EXPERIMENT Minimization
        EXPERIMENT On
        if ApplyConstraintsOrCharges == 'Yes'
          (ConstraintsChargesName) '(LigandName)'
        WAIT ExpEnd
        # restore the pre-existing state
        JOINRES HOH
        SPLITOBJ 1
        TotalObjects = COUNTOBJ All
        JOINOBJ 1,3
        for i = 4 to ((CSN)+2) 
        RENUMBEROBJ 4,1
        for i = 5 to TotalObjects
          JOINOBJ (i), 1
        RENAMEOBJ 1, protein
        RENAMEOBJ 3, solution
        SAVESCE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)water
      else
        MinimizedStructureExists = FILESIZE (MacroTarget)_(WarmUpRegime)_01_LS(LincsSettle)water.sce 
        if (MinimizedStructureExists)
          shell cp (MacroTarget)_(WarmUpRegime)_01_LS(LincsSettle)water.sce (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)water.sce
        else 
          LOADSCE (MacroTarget)_SaltWaterBoxNotMinimized
          # now a fix for water 
          FillCellWater Density=0.997,Probe=1.4,BumpSum=1.0,DisMax=0
          DELRES OBJ 4 res HOH with distance > 10 from res protein
          TESTOBJ4 = COUNTOBJ All
          if TESTOBJ4 > 3.5
            DELRES OBJ 4 res HOH with distance < 2 from OBJ 1 3
          TestOBJ4 = COUNTOBJ All
          if TestOBJ4 > 3.5
            JOINOBJ 4,3 
          SAVESCE (MacroTarget)_SaltWaterBoxNotMinimized
          EXPERIMENT Minimization
          EXPERIMENT On
          if ApplyConstraintsOrCharges == 'Yes'
            (ConstraintsChargesName) '(LigandName)'
          WAIT ExpEnd
          SAVESCE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)water
     
  else
    RAISEEROR Could not find (MacroTarget)_NotMinimized.sce, needed later in procedure
  
  
  
  # determine whether NAC is achieved, and tabulate it, after that reset it again 
  LOADSCE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)water.sce
  
  
  
  for i = 1 to count L
    b = '(L(i))'
    for c = 1 to count Label_(b)_
      (Label_(b)_(c))_Q_Seed(CSN)_M   = ((Label_(b)_(c))_Measurement) 
    
    MagicNumber = 1 + (2*(count Label_(b)_))
    # determine whether NAC is achieved 
    NAC_Seed(CSN)_M = 1
    for c = 1 to count Label_(b)_
      (Label_(b)_(c))_NAC_Seed(CSN)_M = 0.00000
      if ( ((Label_(b)_(c))_Q_Seed(CSN)_M)  > ((Label_(b)_(c))_Q_Min))
        if (((Label_(b)_(c))_Q_Seed(CSN)_M)  < ((Label_(b)_(c))_Q_Max))
          (Label_(b)_(c))_NAC_Seed(CSN)_M  = 1.00000
      NAC_Seed(CSN)_M = (NAC_Seed(CSN)_M) * ((Label_(b)_(c))_NAC_Seed(CSN)_M)  
    TABULATE (100.000 * (NAC_Seed(CSN)_M))
    for c = 1 to count Label_(b)_
      TABULATE ((Label_(b)_(c))_Q_Seed(CSN)_M)
      TABULATE ((Label_(b)_(c))_NAC_Seed(CSN)_M)
    
    LabelCollection = ' ________NAC?'
    for c = 1 to count Label_(b)_
      LabelCollection = '(LabelCollection) (Label_(b)_(c)) CriteriaMet?'
    SaveTab 1,(MacroTarget)_NAC_Results_(L(i))_Seed(CSN)_Minimized,Format=Text,Columns=(MagicNumber),NumFormat=%12.2f,(LabelCollection)
    DELTAB 1
   
  # use a cutoff for the non-bonded of 7.86 Angstrom
  CUTOFF 7.86
  
  # calculate further away than 7.86 Angstrom with Particle Mesh Ewald (PME) algorithm
  LONGRANGE coulomb
  
  # use the 4th degree splines for PME
  SIMSPEED normal
  
  # Correct for the diffusion of the protein
  CORRECTDRIFT On
  
  # settings for timesteps 
  if LincsSettle == 'Off'
    # use a timestep of 1.25 fs, update the non-bonded interactions every 2 timesteps
    CalculationTimeStep = 1.250
    UpdateCycle         = 2
    TIMESTEP (UpdateCycle),(CalculationTimeStep)
    # remove LINCS and SETTLE constraints
    FREEBOND All, All
    FREEANGLE All, All, All
  else
    if  LincsSettle == 'On' 
      # make water rigid with SETTLE
      FixHydAngle all
      # constrain hydrogen atom bond distances with LINCS
      FIXBOND All, element H
      # use a timestep of 1.333333 fs, update the non-bonded interactions every 3 timesteps
      CalculationTimeStep = 2.5
      UpdateCycle         = 2
      TIMESTEP (UpdateCycle),(CalculationTimeStep)
    else    
      RAISERROR Illegal LincsSettle defined in script, (LincsSettle) should be either 'On' or 'Off'
         
  # use standard SI units rather than kcal/mol or pound per square inch or stones or other medieval units that are only used in three countries in the world.
  ENERGYUNIT kj/mol
  
  # exit if a warning occurs
  WARNISERROR On
  
  # Start the simulation 
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  SIM On
  CORRECTDRIFT On
  if ApplyConstraintsOrCharges == 'Yes'
    (ConstraintsChargesName) '(LigandName)'
  # Keep the pressure constant with the solvent density
  PRESSURECTRL Type=SolventProbe,Name=HOH,Density=0.997,Axis=XYZ
  
  # Control the temperature by rescaling the velocities
  TEMPCTRL Rescale
  if ApplyConstraintsOrCharges == 'Yes'
    (ConstraintsChargesName) '(LigandName)'
  FREE
   
  # save snapshots with an appropriate name at appropriate intervals
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  IntervalsSnapShots = (SnapShotInterval) / ( (CalculationTimeStep)*(UpdateCycle ))
  i = 00000
  FileName = '(MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(i).sim'
  SAVEsim (FileName),(IntervalsSnapShots)

  # Check if the simulation was already done?
  FinalTime = 0+ (( (WarmUpTime)+(EquilibrationTime)+(ProductionTime)) / 1000)

  # *********************** do the WARMING UP or load a snapshot of it after it is done ***************************************************************************************************
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # check if the snapshot at which time the simulation is done is there
  WarmUpDone = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)WarmUpDone.sim
  if (WarmUpDone)
    LOADSIM (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)WarmUpDone.sim
  else
    # the 5 and 10 come from the do while loop below
    NumberOfTemperatureJumps = 0 + (( (AssayTemperature)  - 5) / 10)
    TimeForEquilibrationInTimeSteps = (WarmUpTime) / ( (CalculationTimeStep)*(UpdateCycle ))
    CurrentWaitPeriod = (TimeForEquilibrationInTimeSteps) / (NumberOfTemperatureJumps)

    # start at 5 K, not 0 K, to conserve the existing motions, their speeds are rescaled only, seems like a poor idea if the original speed was 0 K. 
    t = 5

    # The following section sets the temperature with different random seeds. 
    # reliable than the following method. 
    TEMP ((SeedingNumber)+(t))
    # This do while loop increases the temperature, in steps of 10 K to allow the temperature to equilibrate before the next step 
    do
      TEMP (t), REASSIGN=NO
      wait (CurrentWaitPeriod)
      t = (t) + 10
    while (t < (AssayTemperature))
    # save the snapshot of the warmup being done
    SAVESIM  (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)WarmUpDone
    # ensure still saved continuously (unclear if this is necessary)
    SAVESIM (FileName),(IntervalsSnapShots)

    
  # ensure the temperature and time is set correctly now
  TEMP (AssayTemperature), REASSIGN=NO  
  Time (WarmUpTime)
  
  # ********************** do the equilibrationphase or load a snapshot of it after it was done *******************************************************************************************
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  SnapShotEquilibrationDone = 00000 + ( ((WarmUpTime) + (EquilibrationTime) )/(SnapShotInterval) )
  SnapShotWarmUpDone = 00000 + 1+ ( ((WarmUpTime)  )/(SnapShotInterval) )
  EquilibrationDone = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(SnapShotEquilibrationDone).sim
  if (EquilibrationDone)
    LOADSIM (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(SnapShotEquilibrationDone).sim
  else
    LastSnapShot = 0
    CorrectionForSavedSnapShots = 0
    # check which snapshot is there, continue with the last saved snapshot during equilbration
    for j = SnapShotWarmUpDone to SnapShotEquilibrationDone
      EquilibrationSnapshotHalfway = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(j).sim
      if (EquilibrationSnapshotHalfway)
        LastSnapShot = (j)
    if (LastSnapShot)
      LOADSIM (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(LastSnapShot).sim    
      # correct for the remaining waiting time
      CorrectionForSavedSnapshots =  (SnapShotInterval) * ((LastSnapShot) - (SnapShotWarmUpDone))
    EquilibrationTime = (EquilibrationTime) - (CorrectionForSavedSnapShots)
    WaitingTime = (EquilibrationTime) / ( (CalculationTimeStep)*(UpdateCycle ))
    Wait (WaitingTime)
  
  
  # **************************** This is the PRODUCTION run part, sample NACs on the fly **************************************************************************************************
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  # some time tracking calculations and checking whether part of the production has already been done before running the production phase loop
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # How long does a NAC take expressed in intervals rather than fs?
  SamplingNACInterval = 0 + (NACSamplingInterval) / ( (CalculationTimeStep)*(UpdateCycle ))
  
  # How often to sample a NAC?
  TotalNACIntervals = 0 + ((ProductionTime) / (NACSamplingInterval))
    
  # some initialization
  LastSnapShot = (SnapShotEquilibrationDone)
  LastTimeCounter = 0
  
  # Magic number, i.e. time, NACs, and all criteria with their NAC
  #MagicNumber = 2 + (2*(count Label_))
  
  # Are there already snapshots and periodically saved tables available for part of the production run??
  # check which snapshot is there, continue with the last saved snapshot during equilbration
  SnapShotProductionDone = 00000 + ( ((WarmUpTime) + (EquilibrationTime)+ (ProductionTime))/(SnapShotInterval) )
  for j = SnapShotEquilibrationDone to SnapShotProductionDone
    EquilibrationSnapshotHalfway = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(j).sim
    if (EquilibrationSnapshotHalfway)
      TimeCounterfs = 0 +  ((j) - (SnapShotEquilibrationDone)) * (SnapShotInterval)
      LastSnapShot = (j)
      LastTimeCounter = TimeCounterfs
    else
      break      
  if (LastSnapShot > SnapShotEquilibrationDone)
    LOADSIM (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(LastSnapShot).sim
 
  # This keeps track of how far we are into production time
  TimeCounterfs = 0 + (LastTimeCounter)
  StartInterval = 000001 + ((TimeCounterfs)/(NACSamplingInterval))
  # the counter is to coincide with the snapshots, the snapshots for the table are saved independently
  AtTimeSnapShot = (SnapShotInterval) * (NACShotEveryXSnapshots)/ (NACSamplingInterval)
  AlreadyObtainedSnapShotsShouldBe = (StartInterval) / (AtTimeSnapShot)
  CounterToTimeSnapShot = 0 + ((StartInterval) - ( (AlreadyObtainedSnapShotsShouldBe) * (AtTimeSnapShot) )) - 1
  # do the actual sampling, this loop runs during the production time, when it is finished it saves a final table but during the snapshots tables are saved as well
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  for j =  StartInterval to TotalNACIntervals 
    PRINT (j) interval out of (TotalNACIntervals)
    if (CounterToTimeSnapShot == AtTimeSnapShot)
      CurrentSnapShot = (SnapShotEquilibrationDone) + ((TimeCounterfs)/(SnapShotInterval))
      CounterToTimeSnapShot = 0
      # also save a table with only the new NACs 
      for a = 1 to count L
        b = '(L(a))'
        MagicNumber = (2*(count Label_(b)_)) + 2
        for q = ((j) - (AtTimeSnapShot)) to ((j) - 1)
          # TABULATE The results
          TABULATE (0.001 *(q)*(NACSamplingInterval) )
          TABULATE (NAC_(b)_(q))
          for c = 1 to count Label_(b)_
            TABULATE ((Label_(b)_(c)_Q_(q)))
            TABULATE ((Label_(b)_(c))_NAC(q))
        SAVETAB 1, (MacroTarget)_(WarmUpRegime)_(CSN)_(b)_LS(LincsSettle)(CurrentSnapShot),Format=Text,Columns=(MagicNumber),NumFormat=%12.6f, Mixed_Data 
        DELTAB 1
      # also tabulate the number of NACs at the same time and all the desired pairs of NACs to test.
      for q = ((j) - (AtTimeSnapShot)) to ((j) - 1)
        TABULATE (0.001 *(q)*(NACSamplingInterval) ) 
        TABULATE  (NAC_count_total(q)) 
        for a = 1 to count S
          TABULATE (Set_complete_(a)_(q))  
      MagicNumber = (count S) + 2
      SAVETAB 1, (MacroTarget)_(WarmUpRegime)_(CSN)_NAC_combinations_LS(LincsSettle)(CurrentSnapShot),Format=Text,Columns=(MagicNumber),NumFormat=%12.6f, Mixed_Data 
      DELTAB 1 
    CounterToTimeSnapShot = (CounterToTimeSnapShot) + 1
    wait (SamplingNACInterval)
    TimeCounterfs  = (TimeCounterfs) + (NACSamplingInterval)
    # do the measurements and determine whether NAC is achieved. 
    for a = 1 to count L
      b = '(L(a))'
      # do the measurements
      for c = 1 to count Label_(b)_
        Label_(b)_(c)_Q_(j)  = ((Label_(b)_(c))_Measurement)
       # determine whether NAC is achieved 
      NAC_(b)_(j) = 1
      for c = 1 to count Label_(b)_
        (Label_(b)_(c))_NAC(j) = 0.00000
        if (  (Label_(b)_(c)_Q_(j))  > ((Label_(b)_(c))_Q_Min))
          if ((Label_(b)_(c)_Q_(j))  < ((Label_(b)_(c))_Q_Max))
            (Label_(b)_(c))_NAC(j) = 1.000000
        NAC_(b)_(j) = (NAC_(b)_(j)) * (Label_(b)_(c))_NAC(j) 
    # now count how many NACs are present at the same time
    # ----------------------------------------------------
    NAC_count_total(j) = 0.00000
    for a = 1 to count L
      b = '(L(a))'
      if ( (NAC_(b)_(j)) == 1)
        NAC_count_total(j) = (NAC_count_total(j) ) + 1.00000 
    # and determine the pairs of NACs that are present at the same time.
    # ----------------------------------------------------------------- 
    # for all sets
    for a = 1 to count S
      Set_complete_(a)_(j) = 1.00000
      for b = 1 to count SetsTestCombined_(S(a))()
        c = '(SetsTestCombined_(S(a))(b))'
        Set_complete_(a)_(j) = (Set_complete_(a)_(j)) * (NAC_(c)_(j))
        
  # stop the simulations and save a final scene and intermediate tab file
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
  # ensure last snapshots saved first for the tab file
  CurrentSnapShot = (SnapShotEquilibrationDone) + ((TimeCounterfs)/(SnapShotInterval))
  # Only save it if it is not already there (!) 
  LastTabFileAlreadyExists = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_(L(1))_LS(LincsSettle)(CurrentSnapShot).tab
  if (LastTabFileAlreadyExists == 0)
    for a = 1 to count L
      b = '(L(a))'
      MagicNumber = (2*(count Label_(b)_)) + 2
      for q = ((j) - (AtTimeSnapShot) + 1) to ((j) )
        # TABULATE The results
        TABULATE (0.001 *(q)*(NACSamplingInterval) )
        TABULATE (NAC_(b)_(q))
        for c = 1 to count Label_(b)_
          TABULATE ((Label_(b)_(c)_Q_(q)))
          TABULATE ((Label_(b)_(c))_NAC(q))
      SAVETAB 1, (MacroTarget)_(WarmUpRegime)_(CSN)_(b)_LS(LincsSettle)(CurrentSnapShot),Format=Text,Columns=(MagicNumber),NumFormat=%12.6f, Mixed_Data 
      DELTAB 1
    # also tabulate the number of NACs at the same time and all the desired pairs of NACs to test.
    for q = ((j) - (AtTimeSnapShot) + 1) to ((j) )
      TABULATE (0.001 *(q)*(NACSamplingInterval) ) 
      TABULATE  (NAC_count_total(q)) 
      for a = 1 to count S
        TABULATE (Set_complete_(a)_(q))  
    MagicNumber = (count S) + 2
    SAVETAB 1, (MacroTarget)_(WarmUpRegime)_(CSN)_NAC_combinations_LS(LincsSettle)(CurrentSnapShot),Format=Text,Columns=(MagicNumber),NumFormat=%12.6f, Mixed_Data 
    DELTAB 1 
  
  # and here for the sim snapshot
  wait 2
  
  SIM Off
  # How long has it taken in total)
  FinalTime = 0+ (( (WarmUpTime)+(EquilibrationTime)+(ProductionTime)) / 1000)
  
  

# ========================================================================================================================================================================================
# =========================  Part D of the Script: The Analysis of the MD simulation(s)===================================================================================================
# ========================================================================================================================================================================================


# ***************************** This Loop analyses the NAC data for every seed individually *********************************************************************************************
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
for a = 1 to count L
  b = '(L(a))'
  for CSN = 01 to (MaximumSeeds)  
    # save two sets of tables, one with the original sampled data, one with the averaged data and their standard deviations
    # ------------------------------------------------
    
    # first load the original data from the intermediate tables
    CounterToTimeSnapShot = (AtTimeSnapShot)
    PreviousSnapShot	  = (SnapShotEquilibrationDone)
    PreviousTimeCounter   = 0
    q = 000001
    for j =  000001 to TotalNACIntervals
      if (CounterToTimeSnapShot == AtTimeSnapShot)
        CounterToTimeSnapShot	= 0
        q			= 000001
        LastSnapShot		= (PreviousSnapShot) + (NACShotEveryXSnapshots) 
        LastTimeCounter 	= (PreviousTimeCounter) + (NACShotEveryXSnapshots) * (SnapShotInterval)
        LOADTAB (MacroTarget)_(WarmUpRegime)_(CSN)_(b)_LS(LincsSettle)(LastSnapShot) 
        CounterToTimeSnapShot = 0
        PreviousSnapShot      = (LastSnapShot)    
        PreviousTimeCounter   = (LastTimeCounter) 
        RawData(b)() = TAB 1
        DELTAB 1
      # now convert this loaded table into the lists 
      MagicNumber = (2*(count Label_(b)_)) + 2 
      TotalNacIntervalsAlreadyAnalyzed = 00000 + (count RawData)/ (MagicNumber)
      NACSamplingIntervalList(b)(j) = 1000.000 *(RawData(b)( ((MagicNumber)*((q)-1)) +1))
      NAC(b)(j)			 =	     (RawData(b)( ((MagicNumber)*((q)-1)) +2))
      for c = 1 to count Label_(b)_
        (Label_(b)_(c))_Q_(j)	   =	       (RawData(b)( ((MagicNumber)*((q)-1)) +((2*(c))+1) ))
        (Label_(b)_(c))_NAC(j)	   =	       (RawData(b)( ((MagicNumber)*((q)-1)) +((2*(c))+2) ))
      # to see how far it is in the log file, if it is not simply stock
      if MoreLOGFile == 'TRUE'
        print reading in (q) out of (TotalNacIntervalsAlreadyAnalyzed) in this tab file
        print reading in (j) out of (TotalNacIntervals) in total for this MD run
      # update counter
      CounterToTimeSnapShot = (CounterToTimeSnapShot) + 1   
      q 		    = (q) + 1
    
    # Save a table with all the data. 
    for j = 000001 to TotalNACIntervals
      # TABULATE The results
      TABULATE (0.001 *(j)*(NACSamplingInterval) )
      TABULATE (NAC(b)(j))
      for c = 1 to count Label_(b)_
        TABULATE ((Label_(b)_(c))_Q_(j))
        TABULATE ((Label_(b)_(c))_NAC(j))
    LabelCollection = '____Interval ________NAC?'
    for c = 1 to count Label_(b)_
      LabelCollection = '(LabelCollection) (Label_(b)_(c)) CriteriaMet?'
    SAVETAB 1,(MacroTarget)_(WarmUpRegime)_(CSN)_(b)_LS(LincsSettle)(EquilibrationTime)fs_(ProductionTime)fs_NAC_DATA_ProductionRun,Format=Text,Columns=(MagicNumber),NumFormat=%12.2f,(LabelCollection)
    DELTAB 1
    
    # TABULATE The interval results, see if even distribution of NACs over the intervals or not. 
    if InterValAnalysis == 'On'
      for q = 1 to count TrajectorySubAnalysisIntervals
        CountTrajectoryInterVal = 1
        CurrentTimeInterValEnd = 0
        EndTrajectoryCount =  (TrajectorySubAnalysisIntervals(q))/ (NACSamplingInterval) 
        for j = 000001 to TotalNACIntervals
          #print REPORT (j) (CountTrajectoryInterVal)
          TemporaryIntervalNAC(CountTrajectoryInterVal) = (NAC(b)(j))
          if (CountTrajectoryInterval ==  EndTrajectoryCount)
            CurrentTimeInterValEnd = (CurrentTimeInterValEnd) + (TrajectorySubAnalysisIntervals(q))  
            TABULATE (CurrentTimeInterValEnd)
            CurrentTrajectoryAverage  = (100.0000000000000 * (mean TemporaryIntervalNAC))
            TABULATE (CurrentTrajectoryAverage)
            #print REPORT Current number of intervals is (count TemporaryIntervalNAC)
            CountTrajectoryInterVal = 1
            TemporaryIntervalNAC() = 0
          else 
            CountTrajectoryInterVal = (CountTrajectoryInterVal)+ 1
        SAVETAB 1,(MacroTarget)_(WarmUpRegime)_(CSN)_(b)_LS(LincsSettle)(EquilibrationTime)fs_(ProductionTime)fs_NAC_Interval(TrajectorySubAnalysisIntervals(q))fsProductionRun,Format=Text,Columns=2,NumFormat=%12.2f,____Interval _____NACPerc  
        DELTAB 1
      
    TABULATE 'Total Number'
    TABULATE (TotalNACIntervals)
    TABULATE 'NA'
    TABULATE 'NA'
    
    TABULATE 'Total NACs'
    TABULATE (sum NAC(b))
    TABULATE 'NA'
    TABULATE 'NA'
    
    TABULATE 'TotalTime ps'
    TABULATE ((ProductionTime)/1000)
    TABULATE 'NA'
    TABULATE 'NA'
    
    TABULATE 'NACPercentag'
    TABULATE (100.000*(mean NAC(b)))
    TABULATE 'NA'
    TABULATE 'NA'
    
    for c = 1 to count Label_(b)_
      TABULATE '(Label_(b)_(c))'
      TABULATE (mean (Label_(b)_(c))_Q_)
      TABULATE (stddev (Label_(b)_(c))_Q_)
      TABULATE (100.000*(mean (Label_(b)_(c))_NAC))
    SAVETAB 1,(MacroTarget)_LS(LincsSettle)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_NAC_Results_(b)_Seed(CSN)_Summary,Format=Text,Columns=4,NumFormat=%12.3f,___Parameter Average_Value Standard_Dev PassCritPerc
    DELTAB 1

# ***************************** This Loop analyses the NAC COMBINATIONS data for every seed individually *********************************************************************************
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
for CSN = 01 to (MaximumSeeds)  
  # save two sets of tables, one with the original results per time point, one with the averaged results and their standard deviations
  # ----------------------------------------------------------------------------------------------------------------------------------
  # first load the original data from the intermediate tables
  CounterToTimeSnapShot = (AtTimeSnapShot)
  PreviousSnapShot	= (SnapShotEquilibrationDone)
  PreviousTimeCounter   = 0
  q = 000001
  for j =  000001 to TotalNACIntervals
    if (CounterToTimeSnapShot == AtTimeSnapShot)
      CounterToTimeSnapShot   = 0
      q		              = 000001
      LastSnapShot	              = (PreviousSnapShot) + (NACShotEveryXSnapshots) 
      LastTimeCounter         = (PreviousTimeCounter) + (NACShotEveryXSnapshots) * (SnapShotInterval)
      LOADTAB (MacroTarget)_(WarmUpRegime)_(CSN)_NAC_combinations_LS(LincsSettle)(LastSnapShot) 
      CounterToTimeSnapShot = 0
      PreviousSnapShot      = (LastSnapShot)	
      PreviousTimeCounter   = (LastTimeCounter) 
      RawData(b)() = TAB 1
      DELTAB 1
    # now convert this loaded table into the lists 
    MagicNumber = (count S) + 2
    TotalNacIntervalsAlreadyAnalyzed = 00000 + (count RawData)/ (MagicNumber)
    NACSamplingIntervalList(b)(j) = 1000.000 *(RawData(b)( ((MagicNumber)*((q)-1)) +1))
    NAC_TotalCount(b)(j)		               =	   (RawData(b)( ((MagicNumber)*((q)-1)) +2))
    for a = 1 to count S
      Combination_(a)_(j)	 =	     (RawData(b)( ((MagicNumber)*((q)-1)) +2+(a)))
    # to see how far it is in the log file, if it is not simply stock
    if MoreLOGFile == 'TRUE'
      print reading in (q) out of (TotalNacIntervalsAlreadyAnalyzed) in this tab file
      print reading in (j) out of (TotalNacIntervals) in total for this MD run
    # update counter
    CounterToTimeSnapShot = (CounterToTimeSnapShot) + 1   
    q		  = (q) + 1
  
  # Save a table with all the data. 
  for j = 000001 to TotalNACIntervals
    # TABULATE The results
    TABULATE (0.001 *(j)*(NACSamplingInterval) )
    TABULATE (NAC_TotalCount(b)(j))
    for a = 1 to count S
      TABULATE (Combination_(a)_(j))
  LabelCollection = '____Interval NumberOfNACs'
  for c = 1 to count S
    LabelCollection = '(LabelCollection)            (L(c))'
  SAVETAB 1,(MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(EquilibrationTime)fs_(ProductionTime)fs_Combinations_NAC_DATA_ProductionRun,Format=Text,Columns=(MagicNumber),NumFormat=%12.2f,(LabelCollection)
  DELTAB 1
  
  TABULATE 'Total Number'
  TABULATE (TotalNACIntervals)
  TABULATE 'NA'
  
  TABULATE 'TotalTime ps'
  TABULATE ((ProductionTime)/1000)
  TABULATE 'NA'

  TABULATE 'Average#NACs'
  TABULATE (100.000*(mean NAC_TotalCount(b)))
  TABULATE 'NA'
  
  
  for c = 1 to count S
    TABULATE '(L(c))'
    TABULATE (100.000* (mean (Combination_(c)_)))
    TABULATE (100.000* stddev (Combination_(c)_))
  SAVETAB 1,(MacroTarget)_LS(LincsSettle)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_NAC_ResultsCombined_Seed(CSN)_Summary,Format=Text,Columns=3,NumFormat=%12.3f,___Parameter Average_Value Standard_Dev PassCritPerc
  DELTAB 1




# This calculates the averages and standard deviations of the ensemble of the NACs rather than that of the individual runs
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
for a = 1 to count L
  b = '(L(a))'
  for CSN = 01 to (MaximumSeeds)
    # loat the table and convert it to a temporary list
    LOADTAB (MacroTarget)_LS(LincsSettle)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_NAC_Results_(b)_Seed(CSN)_Summary
    IntermediateList(b)() = TAB 1
    DELTAB 1
    # convert the variables
    #PRINT round (CSN) of CSN = 01 to (MaximumSeeds)
    TotalNumbers(b)(CSN)	 = (IntermediateList(b)(3))
    TotalNacs(b)(CSN)	 = (IntermediateList(b)(8))
    TotalTimes(b)(CSN)	 = (IntermediateList(b)(13))
    NACPercentages(b)(CSN)  = (IntermediateList(b)(17))
    for c = 1 to count Label_(b)_
      
      (Label_(b)_(c))s(CSN)  = (IntermediateList(b)(17 + (c*4)))
      (Label_(b)_(c))SD(CSN) = (IntermediateList(b)(17 + (c*4)+1))
      (Label_(b)_(c))Pass(CSN) = (IntermediateList(b)(17 + (c*4)+2))	

  # Make a table with all the numbers, both mean and standard deviation
  
  TABULATE 'Total Number'
  TABULATE (mean TotalNumbers(b))
  if MS == 1
    TABULATE 'NA'
  else 
    TABULATE (stddev TotalNumbers(b))
  for i = 1 to 4 
    TABULATE 'NA'
  
  TABULATE 'Total NACs'
  TABULATE (mean   TotalNacs(b))
  if MS == 1
    TABULATE 'NA'
  else 
    TABULATE (stddev TotalNacs(b))
  for i = 1 to 4 
    TABULATE 'NA'
  
  TABULATE 'TotalTime ps'
  TABULATE (mean     TotalTimes(b))
  if MS == 1
    TABULATE 'NA'
  else 
    TABULATE (stddev   TotalTimes(b))
  for i = 1 to 4 
    TABULATE 'NA'
  
  TABULATE 'NACPercentag'
  TABULATE (mean       NACPercentages(b))  
  if MS == 1
    TABULATE 'NA'
  else 
    TABULATE (stddev	 NACPercentages(b))  
  for i = 1 to 4 
    TABULATE 'NA'
  
  for c = 1 to count Label_(b)_
    TABULATE '(Label_(b)_(c))'
    TABULATE (mean     (Label_(b)_(c))s)	
    if MS == 1
      TABULATE 'NA'
    else 
      TABULATE (stddev   (Label_(b)_(c))s)   
    TABULATE (mean     (Label_(b)_(c))SD) 
    if MS == 1
      TABULATE 'NA'
    else 
      TABULATE (stddev   (Label_(b)_(c))SD) 
    TABULATE (mean (Label_(b)_(c))Pass)   
    if MS == 1
      TABULATE 'NA'
    else 
      TABULATE (stddev   (Label_(b)_(c))Pass)  
  
  SAVETAB 1, (MacroTarget)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_LS(LincsSettle)_All(MS)Seeds_NAC_Results_(b)_Summary,Format=Text,Columns=7,NumFormat=%12.3f,___Parameter Average_Value SD_from average  mean SD Standard_Dev PercCritPass Standard_Dev
  DELTAB 1

# This calculates the averages and standard deviations of the ensemble of the NAC COMBINATIONS rather than that of the individual runs
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


for CSN = 01 to (MaximumSeeds)
  # loat the table and convert it to a temporary list
  LOADTAB (MacroTarget)_LS(LincsSettle)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_NAC_ResultsCombined_Seed(CSN)_Summary 
  IntermediateListCombinations() = TAB 1
  DELTAB 1
  # convert the variables
  TotalNumbersCombi(CSN)  = (IntermediateListCombinations(3))
  TotalTimesCombi(CSN)      = (IntermediateListCombinations(7))
  NumbersNACsCombi(CSN)     = (IntermediateListCombinations(10))
  for c = 1 to count S 
    Set_(c)mean(CSN)  = (IntermediateListCombinations(10 + (c*3)))
    Set_(c)SD(CSN) = (IntermediateListCombinations(10 + (c*3)+1))
           
# Make a table with all the numbers, both mean and standard deviation

TABULATE 'Total Number'
TABULATE (mean TotalNumbersCombi)
if MS == 1
  TABULATE 'NA'
else 
  TABULATE (stddev TotalNumbersCombi)
for i = 1 to 2 
  TABULATE 'NA'

TABULATE 'Avg NAC sum'
TABULATE (mean   NumbersNACsCombi)
if MS == 1
  TABULATE 'NA'
else 
  TABULATE (stddev NumbersNACsCombi)
for i = 1 to 2 
  TABULATE 'NA'

for c = 1 to count S
  TABULATE '(S(c))'
  TABULATE (mean     Set_(c)mean)        
  if MS == 1
    TABULATE 'NA'
  else 
    TABULATE (stddev   Set_(c)mean)   
  TABULATE (mean     Set_(c)SD) 
  if MS == 1
    TABULATE 'NA'
  else 
    TABULATE (stddev   Set_(c)SD) 

SAVETAB 1, (MacroTarget)_LS(LincsSettle)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_NAC_Results_ZZZ_Combinations_All(MS)Seeds_Summary,Format=Text,Columns=5,NumFormat=%12.3f,___Parameter Average_Value SD_from average  mean SD Standard_Dev 
DELTAB 1





# *************** Now also do the analysis of the snapshots: record energy and RMSD and prepare an averaged structure of the complex with Bfactors (during the equilibration time)****** 
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


LOADSCE (MacroTarget)_SaltWaterBoxNotMinimized.sce
# Analyze the Bfactors of the crystal structure (or at least the earlier structure)
ResidueListProtein() = LISTRES protein, FORMAT=RESNUM
BfactorsOriginal()   = BFACTORATOM Res Protein atom CA
for i = 1 to count BfactorsOriginal
  RMSFOriginal(i) = SQRT (BfactorsOriginal(i)*0.037995443)
# Make a copy to use later for RMSD calculations
OriginalStructure = DUPLICATEOBJ 1
REMOVEOBJ (OriginalStructure)

for CSN = 01 to (MaximumSeeds)  
  i = 00000
  NextSnapShotExists = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(i).sim
  while (NextSnapShotExists)
    LOADSIM (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(i).sim
    CurrentTime = TIME
    TimeSnapShot(i) =(CurrentTime) / 1000
    EnergySnapShot(i) = ENERGY
    SIM OFF
    # Now analyse the rmsd from the crystal structure or designed structure
    ADDOBJ (OriginalStructure)
    RMSDCASnapShot(i)       = SUPATOM OBJ 1 atom CA, OBJ 4 atom CA
    RMSDBackBoneSnapShot(i) = SUPATOM OBJ 1 atom backbone, OBJ 4 atom backbone
    RMSDAllHeavySnapShot(i) = SUPATOM OBJ 1 protein atom element !H, OBJ 4 protein atom element !H
    # if during production phase
    if ((TimeSnapShot(i)*1000) >= ((WarmUpTime)+(EquilibrationTime)))
      ADDPOSATOM OBJ 1
    REMOVEOBJ (OriginalStructure)
    i = i + 1
    NextSnapShotExists = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(i).sim
  # save a table with the energy and RMSD versus time
  for i = 00000 to ((count TimeSnapShot) - 1)
    TABULATE (TimeSnapShot(i))
    TABULATE (EnergySnapShot(i))
    TABULATE (RMSDCASnapShot(i))
    TABULATE (RMSDBackBoneSnapShot(i))
    TABULATE (RMSDAllHeavySnapShot(i))
  # ensure to reset them in case some snapshots are missing
  for i = 00000 to ((count TimeSnapShot) - 1)
    TimeSnapShot(i)         = 0
    EnergySnapShot(i)       = 0
    RMSDCASnapShot(i)       = 0
    RMSDBackBoneSnapShot(i) = 0
    RMSDAllHeavySnapShot(i) = 0
  SAVETAB 1,(MacroTarget)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_(CSN)_LS(LincsSettle)_EnergyRMSDTime,Format=Text,Columns=5,NumFormat=%12.3f,____Timeinps _Energykjmol ______RMSDCA ____RMDSBack ___RMSDHeavy 
  DELTAB 1


# Take the average of the positions to save a PDB and YOB File of it
AVERAGEPOSATOM OBJ 1

# make a list of the CA RMSD values
LongListRMSFMDRun()   = RMSFATOM OBj 1,UNIT=A
LongListAtoms()       = LISTATOM OBJ 1
j = 1
k = 1
for i = 1 to count LongListRMSFMDRun
  CurrentAtomName = NAMEATOM Atom (LongListAtoms(i))
  if (CurrentAtomName == 'CA')
    ShortListRMSFMDRun(j) = (LongListRMSFMDRun(i))
    j = j + 1
  CurrentElement = ElementAtom (LongListAtoms(i))
  if (CurrentElement > 1.5)
    HeavyAtomsRMSFMDRun(k) = (LongListRMSFMDRun(i))
    k = k + 1
        

# Take the average of the positions to save a PDB and YOB File of it
if UseRMSFInPdbYOBFile== 'TRUE'
  for i = 1 to count LongListRMSFMDRun
    BFACTORATOM (i), (LongListRMSFMDRun(i))
  SAVEPDB OBJ 1, (MacroTarget)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_All(MS)Seeds_LS(LincsSettle)_AvgRMSF
else
  RMSFATOM OBj 1,UNIT=bfactor
  for i = 1 to count LongListRMSFMDRun
    Value10timesTooHigh =  BFACTORATOM (i)
    if Value10timesTooHigh < 9999
      BFACTORATOM (i), ((Value10timesTooHigh)/10)
    else
      BFACTORATOM (i),999.9
  SAVEPDB OBJ 1, (MacroTarget)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_All(MS)Seeds_LS(LincsSettle)_Avg




# save a table suitable for a plot of Bfactor versus residue number
for i = 1 to count ResidueListProtein
  TABULATE (ResidueListProtein(i))
  TABULATE (RMSFOriginal(i))
  TABULATE (ShortListRMSFMDRun(i))
SAVETAB 1,(MacroTarget)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_All(MS)Seeds_LS(LincsSettle)_RMSF,Format=Text,Columns=3,NumFormat=%12.3f,__ResnNumber RMSFOriginal ___RMSFMDrun 
DELTAB 1






# if we are doing MD for thermostability, save yob files of averages of trajectories
if (CurrentSetting == '2ndScreening')
  # Now redo the procedure of getting RMSF values for the individual trajectories, store the lists separately
  for CSN = 01 to (MaximumSeeds)  
    # LOAD the Scene again
    CLEAR
    LOADSCE (MacroTarget)_SaltWaterBoxNotMinimized.sce
    # Make a copy to use later for RMSD calculations
    OriginalStructure = DUPLICATEOBJ 1
    REMOVEOBJ (OriginalStructure)
    i = 00000
    NextSnapShotExists = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(i).sim
    # only load the snapshot if not above specified range
    # Criterium stopping
    print test1
    while (NextSnapShotExists)
      print test2
      LOADSIM (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(i).sim
      CurrentTime = TIME
      TimeSnapShot(i) =(CurrentTime) / 1000
      SIM OFF
      # Now analyse the rmsd from the crystal structure or designed structure
      # if during production phase
      if ((TimeSnapShot(i)*1000) >= ((WarmUpTime)+(EquilibrationTime)))
        print test3
        ADDOBJ (OriginalStructure)
        SUPATOM OBJ 1 atom CA, OBJ (OriginalStructure) atom CA
        ADDPOSATOM OBJ 1 
        REMOVEOBJ (OriginalStructure)
      i = i + 1
      if (i) > (SnapShotProductionDone)
        NextSnapShotExists = 0
      else 
        NextSnapShotExists = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(i).sim
    AVERAGEPOSATOM OBJ 1
    LongListRMSFMDRun(CSN)_()	= RMSFATOM OBj 1 ,UNIT=A
    for i = 1 to count LongListRMSFMDRun
      BFACTORATOM (i), (LongListRMSFMDRun(CSN)_(i))
    # Also get rid of all the mobile water atoms and ions, 2 seems to be the right level here
    DELRES HOH CIP CIM bfactor > 2
    SAVEYOB OBJ 1, (MacroTarget)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_(CSN)_LS(LincsSettle)_Avg 
















# ******************** if desired, get the deviation from the original structure of the average structure versus time, average of all the seeds *****************************************
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if CalculateAverageTrajectory == 'On'
  for i = 00000 to ((count TimeSnapShot) - 1)
    CLEAR
    GreenFlag = 1
    LOADSCE (MacroTarget)_SaltWaterBox.sce
    OriginalStructure = DUPLICATEOBJ 1
    REMOVEOBJ (OriginalStructure)
    for CSN = 01 to (MaximumSeeds)
      NextSnaphotExists = FILESIZE (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(i).sim
      if (NextSnaphotExists)
        LOADSIM (MacroTarget)_(WarmUpRegime)_(CSN)_LS(LincsSettle)(i).sim
        CurrentTime = TIME
        TimeSnapShot(i) =(CurrentTime) / 1000
        SIM Off
        ADDOBJ (OriginalStructure)
        SUPATOM OBJ 1 protein atom element !H, OBJ 4 protein atom element !H
        ADDPOSATOM OBJ 1
        REMOVEOBJ (OriginalStructure)
      else 
        GreenFlag = 0
    if (GreenFlag)  
      AVERAGEPOSATOM OBJ 1
      ADDOBJ (OriginalStructure)
      RMSDCASnapShot(i)       = SUPATOM OBJ 1 atom CA, OBJ 4 atom CA
      RMSDBackBoneSnapShot(i) = SUPATOM OBJ 1 atom backbone, OBJ 4 atom backbone
      RMSDAllHeavySnapShot(i) = SUPATOM OBJ 1 protein atom element !H, OBJ 4 protein atom element !H
    else
      TimeSnapShot(i)         = 0
      RMSDCASnapShot(i)       = 0
      RMSDBackBoneSnapShot(i) = 0
      RMSDAllHeavySnapShot(i) = 0
  for i = 00000 to ((count TimeSnapShot) - 1)
    TABULATE (TimeSnapShot(i))
    TABULATE (RMSDCASnapShot(i))
    TABULATE (RMSDBackBoneSnapShot(i))
    TABULATE (RMSDAllHeavySnapShot(i))
  SAVETAB 1,(MacroTarget)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_All(MS)Seeds_LS(LincsSettle)_RMSDversusTime,Format=Text,Columns=4,NumFormat=%12.3f,____Timeinps ______RMSDCA ____RMDSBack ___RMSDHeavy 
  DELTAB 1


# ******************** Now a short file that contains everything that we need for a quick look at the NAC percentages during the stages and the RMSF ************************************
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# LOAD ALL DATA, put the NAC types all in one table as well to get the highest average out
for a = 1 to count L
  LOADTAB (MacroTarget)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_LS(LincsSettle)_All(MS)Seeds_NAC_Results_(L(a))_Summary
  TemporaryListCompleteData() = TAB 1
  DELTAB 1
  FinalNACPercentage(a) = (TemporaryListCompleteData(26))
# Now determine the highest peak
HighestNACPercentage = max FinalNACPercentage

HighestNAC_identifier = 1
if (count L) >= 1.5
  for a = 1 to count L
  if ((FinalNACPercentage(a))) == (HighestNACPercentage)
    HighestNAC_identifier = (a)
    
# calculate product fractions, or enantiomeric excess
SumNACPercentages = 0.00000
for a = 1 to count L
  SumNACPercentages =(SumNACPercentages)+((FinalNACPercentage(a)))
if (SumNACPercentages) == 0
  for a = 1 to count L
    RelativeFractionPercentage(a) = 0.000000
else
  for a = 1 to count L
    RelativeFractionPercentage(a) = ((FinalNACPercentage(a))) /   (SumNACPercentages)   

# Also load The Sets of combinations
LOADTAB (MacroTarget)_LS(LincsSettle)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_NAC_Results_ZZZ_Combinations_All(MS)Seeds_Summary
TemporaryListCompleteData2() = TAB 1
DELTAB 1
AverageNumberOfNACs = TemporaryListCompleteData2(10)
if count S > 0.5
  for i = 1 to count S
    MagicNumber = 5
    AverageCombination(S(i)) = TemporaryListCompleteData2(10+((MagicNumber)*(i)))

 
### NEXT 30 lines have been commented off
# Load the mutations
#LOADTAB (MacroTarget)_listMutations.tab
#ListMutations() = TAB 1
#DELTAB 1
#
#
## First a simple ranking file
#TABULATE '(L(HighestNAC_identifier))'
#TABULATE (100.000* (RelativeFractionPercentage((HighestNAC_identifier))))
#TABULATE (HighestNACPercentage)
#TABULATE (mean HeavyAtomsRMSFMDRun)
#TABULATE (MS)
#TABULATE '(MacroTarget)'
#TABULATE '(ListMutations(2))'
#SAVETAB 1,(MacroTarget)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_LS(LincsSettle)_All(MS)Seeds_Ranking,Format=Text,Columns=7,NumFormat=%12.3f,__PrefAttack RelativeFrac __HighestNAC ___RMSFHeavy _______Seeds MacroTarget____________________
#DELTAB 1
#
## Then a more detailed ranking file 
## Start with the original results
#TABULATE '(L(HighestNAC_identifier))'
#TABULATE (100.000* (RelativeFractionPercentage((HighestNAC_identifier))))
#TABULATE (HighestNACPercentage)
#TABULATE (mean HeavyAtomsRMSFMDRun)
#TABULATE (MS)
#LabelCollection = '__PrefAttack RelativeFrac __HighestNAC ___RMSFHeavy _______Seeds'
## Tabulate the NAC percentage AND the relative fraction
#TABULATE (SumNACPercentages)
#LabelCollection = '(LabelCollection) __SumNacPerc'
#for a = 1 to count L
#  TABULATE (((FinalNACPercentage(a))))
#  TABULATE (RelativeFractionPercentage(a))
#  LabelCollection = '(LabelCollection) ______NACs_(L(a)) _RelFraction'     
#
## TABULATE ALL The combinations (These need a real name)
#if count S > 0.5
#  for i = 1 to count S
#    TABULATE (AverageCombination(S(i)))
#    LabelCollection = '(LabelCollection) Combination(L(i))'  
#
#TABULATE '(MacroTarget)'
#TABULATE '(ListMutations(2))'
#LabelCollection = '(LabelCollection) MacroTarget____________________'
#SAVETAB 1,(MacroTarget)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_LS(LincsSettle)_All(MS)Seeds_RankingDetailed,Format=Text,Columns=700,NumFormat=%12.3f,(LabelCollection) 
#DELTAB 1




# ******************************* The end of the script, leave automatically if set ******************************************************************************************************
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
# Delete unnecessary bulky files.
if DeleteUnnecessaryFiles=='YES'
  PRINT DELETING UNNECESSARY TAB FILES THAT HAVE BEEN PRODUCED DURING THE RUN. # Easy fix, can be turned off
  for CSN = 01 to (MaximumSeeds)
    SHELL rm (MacroTarget)_LS(LincsSettle)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_NAC_ResultsCombined_Seed(CSN)_Summary.tab
    for a = 1 to count L
      b = '(L(a))'
      SHELL rm (MacroTarget)_LS(LincsSettle)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_NAC_Results_(b)_Seed(CSN)_Summary.tab
      SHELL rm (MacroTarget)_(WarmUpRegime)_(EquilibrationTime)fs_(ProductionTime)fs_LS(LincsSettle)_All(MS)Seeds_NAC_Results_(b)_Summary.tab
      SHELL rm (MacroTarget)_(WarmUpRegime)_(CSN)_(b)_LS(LincsSettle)(EquilibrationTime)fs_(ProductionTime)fs_NAC_DATA_ProductionRun.tab
      SHELL rm (MacroTarget)_(WarmUpRegime)_(CSN)_(b)_LS(LincsSettle)(EquilibrationTime)fs_(ProductionTime)fs_Combinations_NAC_DATA_ProductionRun.tab

# EXIT
if AutomaticExitAtEnd == 'YES'
  EXIT
else
  SHOWMESSAGE Finished with analysis trajectories 
  CONSOLE ON





