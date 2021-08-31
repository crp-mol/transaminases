##################################################################################################################################
#                                                                                                                                #
#  Name Script: generate_rotamers.mcr                                                                                            #
#               -----------------------                                                                                          #       
#                                                                                                                                #
#  Purpose:     This scripts converts the poorly formatted ChemDraw files of the amine portion of the substrates and those of    #
#               the PMP cofactor to a polished fused version which has the right atom names for later docking and MD simulations.#
#                                                                                                                                #
#  User modify: See instructions between two lines of pound (#) signs below.                                                     #
#                                                                                                                                #
#  Author:      Hein J. Wijma, H.J.Wijma@rug.nl                                                                                  #  
#               Script has not been tested on Windows                                                                            #
#               Last modified by C.R. Mar-2020                                                                                   #
#                                                                                                                                #
#  Further:     From Mac TextWrangler: set under preferences -> Appearances _> page guide at 130; under printing set custom font #
#               to courier-7, for easy printing. On screen can be 13 (view).                                                     #
#                                                                                                                                #
##################################################################################################################################

#--------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------- User-defined settings -------------------------------------------------------------#

# Whether or not to stop yasara at the end of this script, YES if you want this.
AutomaticExitAtEnd = 'YES'  # End at the end of the MACRO.
OnError exit                # If the script fails, exit YASARA.

# give a list of macrotargets. Obtain it using ls */*pdb| sed 's/.pdb//g'| tr '\n' '?'| sed "s/?/','/g" in input_amines folder
# ls *.pdb| sed 's/.pdb//g'| tr '\n' '?'| sed "s/?/','/g"
ListMacroTargets() = '01R','01S','06R','06S'
#'01R','01S','02R','02S','03R','03S','04R','04S','05R','05S','06R','06S','07R','07S','08R','08S','09R','09S','10R','10S','11R','11S','12R','12S','13R','13S','14R','14S','15R','15S','16R','16S','17R','17S','18R','18S','19R','19S','20R','20S','21R','21S','22R','22S','23R','23S','24R','24S','25R','25S','26R','26S','27R','27S','28R','28S','29R','29S','30R','30S','31R','31S','32R','32S','33R','33S','34R','34S','35R','35S','36R','36S','37R','37S','38R','38S','39R','39S','40R','40S','41R','41S','42R','42S','43R','43S','44R','44S','45R','45S','46R','46S','47R','47S','48R','48S','49R','49S','50R','50S','51R','51S','52R','52S','53R','53S'

# indicate the name of the cofactor file. 
Cofactor ='PMP'  # Requires PDB file, PMP.pdb
AtomsToBeFrozen       = 'C1  C2  C3  C4   C6  C7  C8  H1  H2  H3  H4  N50  H6  H7  H8  O1  O2  O3  O4  O5  P1 N99 H99 C99 C50 H5 H95'
AtomsToBeSuperimposed = 'C1  C2  C3  C4   C6  C7  C8  H1  H2  H3  H4  N50  H6  H7  H8  O1  O2  O3  O4  O5  P1 '

# cut off to use for RMSD (in Angstrom). If less RMSD than this, not unique enough to serve as a rotamer.
Prune = 0.0050

# This is to eliminate rotamers with a very high energy, higher than the lowest energy rotamer.
MaxEaboveMin = 100

# This may have to become 4-fold higher for serious calculations. 
maxNumberRotamers = 100 

# if Yes, more stops for inspection
OKeveryStep = 'Yes'

# set the number of processors to use
PROCESSORS 4

# Define the dihedral angle chi_1 of the generated rotamers. It can be a list of angles
# -90 means the atom H99 is pointing directly towards Lys285-NH2
DihedralAngleCH() = '-90.00'

# A params file is needed for Rosetta docking. The params file contains information about the ligand.
# It can be obtained by running the molfile_to_params.py script (obtained from Rosetta)
GenerateParamsFile = 'NO' # 'YES'
molfiletoparams = '~/rosetta/main/source/src/python/apps/public/molfile_to_params.py'

#--------------------------------------------------------------------------------------------------------------------------------#

##################################################################################################################################

# for getting different seed numbers
RANDOMSEED time
futureNameTmp = rnd (100000)
futureNameTmp = 0 +(futureNameTmp)


# stuff below can be uncommented to make it possible to switch between single targets and tablefiles listing MacroTargets.
# --------------------------------------------------------------------------------------------------------------------------------

# switch the console off
CONSOLE Off

# self explanatory
SHOWMESSAGE A total  of (count ListMacroTargets) pdb files will be loaded and converted to mol2 files suitable for docking and .com files for calculating RESP charges
WAIT ContinueButton

# --------------------------------------------------------------------------------------------------------------------------------

# Make a directory for the new files. At least on a mac this does not give an overwrite of existing directories.
SHELL mkdir EA_library


# This loops goes through all the targeted structures and does all the work
# ================================================================================================================================
for i = 1 to count ListMacroTargets  
  # first some standard stuff, loading the targets, visualization, ensuring realistic conformation amine substrate by AM1 optimiz. 
  # ------------------------------------------------------------------------------------------------------------------------------
  # CleanUp anything in memory that could disturb the current calculation
  CLEAR
  
  # define macrotarget from the list
  MacroTarget = '(ListMacroTargets(i))'

  # get also the file name and directory name separately.  
  SHELL echo '' > (futureNameTmp).tab
  SHELL echo (MacroTarget)| tr '/' '\n'>> (futureNameTmp).tab
  LOADTAB (futureNameTmp).tab
  MacroTargetShort() = TAB 1
  DELTAB 1
  SHELL rm (futureNameTmp).tab
  
  # ensure there is a directory to later save the files in
  SHELL mkdir EA_library/(MacroTargetShort(1))
  
  # load the pdb file of the cofactor
  LOADPDB (Cofactor)
  
  # Load the pdb file of the amine ligand and zoom it in
  LOADPDB input_amines/(MacroTarget)
  
  # for easier viewing
  BALLSTICKRES ALL
  #ZOOMRES All, steps = 4
  
  # The following lines do an AM1 COSMO optimization on only the amine substrate
  REMOVEOBJ 1
  SHOWMESSAGE Loaded (MacroTarget) and will now optimize its geometry
  QuantumMechanics AM1
  OptimizeRes OBJ 2, method = QM
  HIDEMESSAGE
  ADDOBJ 1
  

  
  # do the superposition of amine substrate on the PMP group
  # ------------------------------------------------------------------------------------------------------------------------------
  
  # first identify the atoms that need to be superimposed. 
  


  # A, the nitrogen atom of PMP that needs to be replaced
  atomA = LISTATOM OBJ 1 atom N4A
  
  # B, the nitrogen atom of the substrates that needs to replace it.
  # for substrates with multiple nitrogen atoms it will automatically choose the primary (not nitro) one.
  countNitrogens = COUNTATOM element N
  if countNitrogens > 1.5
    ListNitrogens() = LISTATOM element N
    for j = 1 to countNitrogens
      DetermineCarbons = COUNTATOM element C with bond to (ListNitrogens(j)) 
      DetermineNitrogens = COUNTATOM element O with bond to (ListNitrogens(j))
      DetermineIfPrimary = DetermineCarbons + DetermineNitrogens
      if DetermineIfPrimary == 1 
        atomB = (ListNitrogens(j))
  else
    atomB = LISTATOM OBJ 2 element N
  
  # C, the carbon atom of PMP that binds the nitrogen
  atomC = LISTATOM OBJ 1 atom C4A
  
  # D, one of the hydrogen atoms of atomB
  atomD, atomD2 = LISTATOM element H with bond to (atomB)
  
  # E, the hydrogen atom of the nitrogen of PMP that is furthest away of the PLP ring
  atomE = LISTATOM OBJ 1 atom HN4 with maximum distance from OBJ 1 atom N1
  
  # F, the carbon atom of the substrate that binds the amine
  atomF= LISTATOM OBJ 2 element C with bond to (atomB)
  
  # now do the superimposing
  SUPORDEREDATOM (atomB),(atomA),(atomD), (atomC),(atomF),(atomE)
  SHOWMESSAGE Done the superposition for (MacroTargetShort(1)), check the result 
  if OKEveryStep == 'Yes'
    WAIT CONTINUEBUTTON
  

  
  # get the atom inventory and the bonds correct. Delete some hydrogens to later add them to get them in SP2-like orientation
  # ------------------------------------------------------------------------------------------------------------------------------
  
  # first rename some of the atoms involved. The reason to do this is that their numbers change due to the deleting of atoms
  RENAMEATOM (atomA),N13
  RENAMEATOM (atomB),N99
  RENAMEATOM (atomC),C50
  
  # now we are at it, also rename the nitrogen atom of the PLP ring  and the carbon atom that is attached to the amine.
  RENAMEATOM OBJ 1 atom N1, N50
  RENAMEATOM (atomF), C99 
      
  # first hydrogens on the PMP nitrogen and the PMP nitrogen
  DELATOM element H with bond to atom N13
  DELATOM N13
  
  # then the hydrogens on the N and C atom that will be connected by a double bond
  DELATOM element H with bond to atom N99
  DELATOM element H with bond to atom C50
  
  # now join the objects and form the bond
  JOINOBJ 2,1
  ADDBOND atom N99, atom C50,2, update=NO
  
  # add hydrogens to the nitrogen and carbon atom that now have a double bond and are thus SP2 hybridized
  ADDHYDATOM atom C50,1, update=No
  ADDHYDATOM atom N99,1, update=No
  
  # self explanatory
  SHOWMESSAGE Fused the molecules, please check the results for (MacroTargetShort(1))
  if OKEveryStep == 'Yes'
    WAIT ContinueButton
      
  
  
  # now a long piece of script to get all the names right 
  # ------------------------------------------------------------------------------------------------------------------------------
  
  # get the residue name (Substrate) and residue number right, it needs to have this name and number for later scripts
  RENAMERES all, SUB
  RENAMEMOL all, X
  RENUMBERRES all, 999
  RENUMBERRES 1000, 999
  JOINRES All
 
  # now store the atom number of those atoms that need specific names. Easiest to store them and later rename them. 
  AtomNumberN50 = LISTATOM atom N50
  AtomNumberN99 = LISTATOM atom N99
  AtomNumberC50 = LISTATOM atom C50
  AtomNumberC99 = LISTATOM atom C99  
  
  # give all the S C H O P N I Br Cl F Se atoms in the file a unique atomname-number. 
  
  # initialize numbering
  NumberS  = 1
  NumberC  = 1
  NumberH  = 1
  NumberO  = 1
  NumberP  = 1
  NumberN  = 1
  NumberI  = 1
  NumberBr = 1
  NumberCl = 1
  NumberF  = 1
  NumberSe = 1
  NumberB  = 1
  AtomListSUB() = LISTATOM Res SUB 999
  
  # do the renumbering and renaming
  for i = 1 to count AtomListSUB
    TestCondition = countatom atom (i) element S
    if (TestCondition)
      RENAMEATOM atom (i), S(NumberS)
      NumberS = (NumberS) + 1
    else
      TestCondition = countatom atom (i) element C
      if (TestCondition)
        RENAMEATOM atom (i), C(NumberC)
        NumberC = (NumberC) + 1
      else
        TestCondition = countatom atom (i) element H
        if (TestCondition)
          RENAMEATOM atom (i), H(NumberH)
          NumberH = (NumberH) + 1
        else
          TestCondition = countatom atom (i) element O
          if (TestCondition)
            RENAMEATOM atom (i), O(NumberO)
            NumberO = (NumberO) + 1
          else
            TestCondition = countatom atom (i) element P
            if (TestCondition)
              RENAMEATOM atom (i), P(NumberP)
              NumberP = (NumberP) + 1
            else
              TestCondition = countatom atom (i) element N
              if (TestCondition)
                RENAMEATOM atom (i), N(NumberN)
                NumberN = (NumberN) + 1
              else
                TestCondition = countatom atom (i) element Cl
                if (TestCondition)
                  RENAMEATOM atom (i), Cl(NumberCl)
                  NumberCl = (NumberCl) + 1
                else
                  TestCondition = countatom atom (i) element Br
                  if (TestCondition)
                    RENAMEATOM atom (i), Br(NumberBr)
                    NumberBr = (NumberBr) + 1
                  else
                    TestCondition = countatom atom (i) element F
                    if (TestCondition)
                      RENAMEATOM atom (i), F(NumberF)
                      NumberF = (NumberF) + 1
                    else
                      TestCondition = countatom atom (i) element I
                      if (TestCondition)
                        RENAMEATOM atom (i), I(NumberI)
                        NumberI = (NumberI) + 1
                      else
                        TestCondition = countatom atom (i) element Se
                        if (TestCondition)
                          RENAMEATOM atom (i), Se(NumberSe)
                          NumberSe = (NumberSe) + 1
                        else
                          TestCondition = countatom atom (i) element B
                          if (TestCondition)
                            RENAMEATOM atom (i), B(NumberB)
                            NumberB = (NumberB) + 1
                          
  
  
  # Now give the original names of important atoms back
  RENAMEATOM (atomNumberN50), N50
  RENAMEATOM (atomNumberN99), N99
  RENAMEATOM (atomNumberC50), C50
  RENAMEATOM (atomNumberC99), C99
  
  # now give the hydrogen atom that is connected to atom C99 the name H99. If it are two, name them H99 and H98. 
  NumberOfHydrogenAtomsThatCanBeAbstracted = COUNTATOM element H with bond to atom C99
  if NumberOfHydrogenAtomsThatCanBeAbstracted == 1
    RENAMEATOM element H with bond to atom C99, H99
  else 
    ListHydrogens() = LISTATOM  element H with bond to atom C99
    RENAMEATOM (ListHydrogens(1)), H99
    RENAMEATOM (ListHydrogens(2)), H98
  
  # now give the hydrogen atom that is connected to atom N99 the name H95. 
  RENAMEATOM element H with bond to atom N99, H95
  
  
  
  # For picture, allow the user to inspect the structure and the relevant atom names of the structure
  # ------------------------------------------------------------------------------------------------------------------------------
  
  # a series of instructions for labeling atoms, self explanatory
  BALLSTICKRES ALL
  LabelAtom element !H,Format=ATOMNAME,Height=0.4,Color=Black,X=0.5,Y=0.5,Z=-1
  LabelAtom atom H99 H98,Format=ATOMNAME,Height=0.4,Color=Black,X=0.5,Y=0.5,Z=-1
  SHOWMESSAGE pres continue button to get view with most hydrogens hidden for (MacroTargetShort(1))
  if OKEveryStep == 'Yes'
    WAIT ContinueButton
  HIDEATOM element H
  SHOWATOM H99 H98
  
  # also make the structure flatter, seems to match TS3 in the Himo paper.
  DIHEDRAL O1, C4, C50, N99, set = 0 
  DIHEDRAL C50, N99, C99, H99, set = -90 
  ANGLE N99, C99, H99, set=90
  
  # get a rotamer that is not clashing. Later removed. For another substrate it caused selection of a rotamer that was clashing. 
  #SAVEYOB OBJ 1, (futureNameTmp).yob
  FIXATOM (AtomsToBeFrozen)
  GROUPATOM (AtomsToBeFrozen), frozen
  #numberRotatableBonds = COUNTBOND !frozen, !frozen, type=rotatable  ##counts the number of rotatable atoms that don-t belong to the frozen group
  numberRotatableBonds = COUNTBOND ALL, ALL, type=rotatable
  numberRotatableBonds = numberRotatableBonds - 6  #6 is the number of rotatable bonds in the PMP alone
  DUPLICATEOBJ 1
  CELL auto, extension = 10
  REMOVEOBJ 2
  referenceEnergy = ENERGYRES OBJ 1
  if numberRotatableBonds > 0.5
    for i = 1 to 100
      SAMPLEDIH OBJ 1 atom element C N O P S B , method=staggered,structures=1,bumpsum=10 #bumpsum=0 still generates some molecules with clashes, you can increase it and it shouldn't affect much
      CELL auto, extension = 10
      newEnergy = ENERGYRES OBJ 1
      if newEnergy < referenceEnergy
        ADDOBJ 2
        DELOBJ 2
        DUPLICATEOBJ 1
        REMOVEOBJ 2
        referenceEnergy = (newEnergy)
      PRINT ++++++++++++++++ referenceEnergy is now (0+(referenceEnergy)) and newEnergy is now (0+(newEnergy)) ===========================
    DELOBJ 1 3
    ADDOBJ 2
    RENUMBEROBJ 2,1
    REMOVEOBJ 2
    
  else  
    SHOWMESSAGE There are no extra rotamers found
  FREE

  SHOWMESSAGE made the structure flatter and in an ideal angle for catalysis, about to sav e the files for (MacroTargetShort(1)). Current view will be saved as *.pdb
  WAIT CONTINUEBUTTON
  
  #Added on Aug13-2019
  PRINT changing the dihedral of the molecule to the wanted value
  DIHEDRAL OBJ 1 atom C50, OBJ 1 atom N99, OBJ 1 atom C99, OBJ 1 atom H99 , set = (DihedralAngleCH)

  # This is to save pdb and mol2 files
  SAVEMOL2 Obj 1, EA_library/(MacroTarget)/(MacroTarget).mol2  
  SAVEMOL2 Obj 1, EA_library/(MacroTarget)/(MacroTarget)_RESP.mol2  
  SAVEPDB Obj 1, EA_library/(MacroTarget)/(MacroTarget).pdb  
  SAVECOM OBJ 1, EA_library/(MacroTarget)/(MacroTarget).com

  #Eliminate objects other than Obj 1
  ADDOBJ 2 
  DELOBJ 2 3



  
  # now all this is to set up the right parameters for the RESP calculation of the .com file 
  # ------------------------------------------------------------------------------------------------------------------------------
  
  # get a random number to get a unique name for the chk file. Otherwise calculations might interfere.
  futureNameChkFile = rnd (100000)
  futureNameChkFile = 0 +(futureNameChkFile)
  
  # this is to give the right set up  in the .com file for a RESP calculation, first AM1, then HF/6-31G*
  SHELL echo '%Mem=8192MB' > (futureNameTmp).txt
  SHELL echo '%NProcShared=4' >> (futureNameTmp).txt
  SHELL echo '%chk=(futureNameChkFile).chk' >> (futureNameTmp).txt
  SHELL echo '#P AM1 Opt=XXXMaxCycle=900,InitialHarmonic=50000YYY'| sed "s/XXX/\(/g"| sed 's/YYY/\)/g' >> (futureNameTmp).txt
  SHELL echo '' >> (futureNameTmp).txt
  SHELL echo 'first AM1 optimization' >> (futureNameTmp).txt
  SHELL echo '' >> (futureNameTmp).txt
  SHELL echo '-1  1' >> (futureNameTmp).txt
  SHELL sed '6,1000!d' EA_library/(MacroTarget)/(MacroTarget).com >> (futureNameTmp).txt
  SHELL echo '--link1--'  >> (futureNameTmp).txt
  SHELL echo '%chk=(futureNameChkFile).chk'  >> (futureNameTmp).txt
  SHELL echo '#P HF/6-31G* Opt=XXXMaxCycle=900,InitialHarmonic=50000YYY SCF=XXXTight,MaxCycle=900YYY Geom=AllCheck Guess=Read Pop=MK IOpXXX6/33=2,6/41=10,6/42=17YYY'| sed "s/XXX/\(/g"| sed 's/YYY/\)/g' >> (futureNameTmp).txt
  SHELL mv (futureNameTmp).txt EA_library/(MacroTarget)/(MacroTarget).com 
  FIXATOM (AtomsToBeFrozen)
  GROUPATOM (AtomsToBeFrozen), frozen
  #numberRotatableBonds = COUNTBOND !frozen, !frozen, type=rotatable  ##counts the number of rotatable atoms that don-t belong to the frozen group
  numberRotatableBonds = COUNTBOND ALL, ALL, type=rotatable
  numberRotatableBonds = numberRotatableBonds - 6  #6 is the number of rotatable bonds in the PMP alone
 
 
  # This is to make a batch_script to run the .com file at the peregrine cluster
  # ------------------------------------------------------------------------------------------------------------------------------
  
  SHELL echo '#!/bin/bash' > (futureNameTmp).txt
  SHELL echo '#SBATCH --time=23:59:00' >> (futureNameTmp).txt
  SHELL echo '#SBATCH --ntasks=1' >> (futureNameTmp).txt
  SHELL echo '#SBATCH --cpus-per-task=4' >> (futureNameTmp).txt
  SHELL echo '#SBATCH --mem=8gb' >> (futureNameTmp).txt
  SHELL echo '#SBATCH --job-name=QM-(MacroTargetShort(1))' >> (futureNameTmp).txt
  SHELL echo '#SBATCH --partition=nodes' >> (futureNameTmp).txt
  SHELL echo '\n' >> (futureNameTmp).txt
  SHELL echo 'module load    Gaussian/09-D.01-GaussView-5.0.9' >> (futureNameTmp).txt  
  SHELL echo 'srun g09 (MacroTargetShort(1)).com (MacroTargetShort(1)).out' >> (futureNameTmp).txt
  SHELL mv (futureNameTmp).txt  EA_library/(MacroTargetShort(1))/batch_script
  SHELL chmod +x EA_library/(MacroTargetShort(1))/batch_script

  # self explanatory
  SHOWMESSAGE Saved EA_library/(MacroTarget).mol2, EA_library/(MacroTarget).pdb and EA_library/(MacroTarget).com
  if OKEveryStep == 'Yes'
    WAIT ContinueButton
     
  
  # it turns out that making rotamers only works when using the original structure, if using the earlier generate structure with zero bumps it does not work. 
  #DELOBJ 1
  #LOADYOB (futureNameTmp).yob
  #ZOOMRES all, 0
  
  
  # This is to make a version of substrates in which there is a second potential hydrogen (primary amines!). 
  # ------------------------------------------------------------------------------------------------------------------------------
  
  # determine if it is a primary amine
  PrimaryAmine = COUNTATOM H98
  if PRimaryAMine
    # change the names.
    RENAMEATOM H99, H97
    RENAMEATOM H98, H99
    RENAMEATOM H97, H98
    DIHEDRAL C50, N99, C99, H99, set = -90 
    ANGLE N99, C99, H99, set=90
    SHOWMESSAGE Detected a 2nd potential hydrogen to abstracted. Will now make files with this hydrogen being called H99 but labelled H98. 
    Wait ContinueButton
    SAVEMOL2 Obj 1, EA_library/(MacroTarget)_H98.mol2
    SAVEPDB Obj 1, EA_library/(MacroTarget)_H98.pdb
    SAVEYOB OBJ 1, (futureNameTmp)_H98.yob 
  
  # now start with making the rotamers. 
  # ------------------------------------------------------------------------------------------------------------------------------
  # it seems the sequence of the atoms makes no difference, try to get all the non-plp atoms and then sample dihedrals of those, after doing that, get the unique once using recycled code.
  
  # it turns out that making rotamers only works when using the original structure, if using the earlier generate structure with zero bumps it does not work. 
  ZOOMRES all, 0
  
  FIXATOM (AtomsToBeFrozen)
  WAIT ContinueButton
  if numberRotatableBonds > 0.5 
    SHOWMESSAGE the number of rotatable bonds = (0+(numberRotatableBonds))
    WAIT ContinueButton
    UNLABELATOM all
    
    # now here a BIG fix was necessary. It turned out that often none of the H99 was pointing out of the plane! Looking at some information online:
    # https://chem.libretexts.org/LibreTexts/University_of_California_Davis/UCD_Chem_231A%3A_Methods_of_Organic_Synthesis/Content/Conformational_Analysis/Acyclic_Conformations
    # shows that the preferred dihedrals with the system are actually covered with the methylene group, not anti or gauche. 
    # Decided therefore to set the dihedral such that the ranges 41 to 139 at both sides are sampled. This should be more realistic then what comes out of it without the fix.
    
    # now something to give an equal number of samples for each angle  
    localNumberRotamers =  (maxNumberRotamers) / (count DihedralAngleCH) 
    PRINT (localNumberRotamers) 
    for j = 1 to count DihedralAngleCH       
      # always make rotamers based on object 1, thus set the dihedral angle there appropriately
      DIHEDRAL OBJ 1 atom C50, OBJ 1 atom N99, OBJ 1 atom C99, OBJ 1 atom H99 , set = (DihedralAngleCH(j))  
      # Do get the indexes for the first and last numbers of the current set of rotamers. Needed later
      firstRotamerCurrentSet = 2 + ( ((j) - 1)* (localNumberRotamers) ) 
      firstRotamerCurrentSet_(i)_(j) = (firstRotamerCurrentSet)
      PRINT IMPORTANT (firstRotamerCurrentSet_(i)_(j)) 
      lastRotamerCurrentSet_(i)_(j)  = 2 + ( ((j) - 0)* (localNumberRotamers) ) -1
      PRINT SUMMARIZED (localNumberRotamers) (firstRotamerCurrentSet) (lastRotamerCurrentSet_(i)_(j))
      
      # now need to do an energy minimization, otherwise due to the setting of the rotamers it too often gets trapped at something with very high energy. 
      REMOVEOBJ all
      ADDOBJ 1
      EXPERIMENT Minimization 
      EXPERIMENT On
      WAIT ExpEnd
      DELOBJ (firstRotamerCurrentSet)
      
      # add a 100 structures, that are subsequently removed, because often the first structures (e.g the first 2, or the first 13), have serious bumps. 
      # Furthermore, it turns the algorithm for random conformations works in some kind of random walk way. So 2000 structures really sample different space than the first 200! Therefore adapt to get  times  more rotamers and subsequently get rid of 80 % of them. 10 fold more gave memory problems with a large substrate.
      #numberOfDevBins =0 + (((numberOfDevbins)/3 -1)/2)
      SAMPLEDIH OBJ 1 atom element C N O P S B , method=staggered,structures=(500+(5*(localNumberRotamers))),bumpsum=10,devmax=25,devbins=4 # Devbins 4 seemed to work well. (numberOfDevBins) #check devmax and devbins
      DELOBJ (firstRotamerCurrentSet)-((firstRotamerCurrentSet)+499)
      WAIT ContinueButton
      for newNumberUncorrected = 1 to localNumberRotamers
        ObjectsToBeDeleted() =((firstRotamerCurrentSet)+500) + (5 * (newNumberUncorrected)) -4, ((firstRotamerCurrentSet)+500) + (5 * (newNumberUncorrected)) -3,((firstRotamerCurrentSet)+500) + (5 * (newNumberUncorrected)) -2,((firstRotamerCurrentSet)+500) + (5 * (newNumberUncorrected)) -1 
        DELOBJ (ObjectsToBeDeleted(1)) (ObjectsToBeDeleted(2)) (ObjectsToBeDeleted(3)) (ObjectsToBeDeleted(4))
      RENUMBEROBJ all

      RENUMBEROBJ all
      WAIT CONTINUEBUTTON
    DELOBJ 1
    RENUMBEROBJ all
    ADDOBJ all
    TotalNumberOfObjects = COUNTOBJ all
    SAVESCE EA_library/(MacroTarget)/(MacroTarget)_AllRotamers  
    SHOWMESSAGE all the (maxNumberRotamers) raw rotamers. Pruning for uniqueness done for all substrates at the end. 
    Wait ContinueButton
    # for debugging get a list of energies
  
  # do the same for the alternative hydrogen atom, in case it existis
  if PRimaryAMine
    # get the right structure
    DELOBJ all
    LOADYOB (futureNameTmp)_H98.yob
    ZOOMRES all, 0
    SHELL rm (futureNameTmp)_H98.yob
    
    # 
    FIXATOM (AtomsToBeFrozen)
    SampleDih OBJ 1 atom element C N O P S B , method=staggered,structures=(maxNumberRotamers),bumpsum=1 
    # Now need to verify that indeed generated rotamers. Found that sometimes no rotamers generated. Then bumpsum needs to increase. 
    RMSDTotal= 0
    CurrentBumpSum= 1
    WHILE RMSDtotal == 0
      for i = 2 to (maxNumberRotamers)
        RMSDCurrent=  RMSDATOM OBJ (i) atom !(AtomsToBeFrozen) element !H, OBJ 1 atom !(AtomsToBeFrozen) element !H
        RMSDtotal = (RMSDCurrent) +(RMSDtotal)
      if RMSDtotal == 0
        CurrentBumpSum= 3 * (CurrentBumpSum)
        DELOBJ 2 - ((maxNumberRotamers) + 1)
        SampleDih OBJ 1 atom element C N O P S B , method=staggered,structures=(maxNumberRotamers),bumpsum=(CurrentBumpSum)
        SHOWMESSAGE redid sampling with bumpsum (CurrentBumpsum)
        WAIT CONTINUEBUTTON
      # some molecules simply have no different rotamers. Like aminopropane. If the bumpsum gets too high, stop it. 
      if RMSDtotal == 0
        if CurrentBumpSum > 100
          RMSDtotal = 1000
         
    
      
    DELOBJ 1
    RENUMBEROBJ all

    SAVESCE EA_library/(MacroTarget)/(MacroTarget)_AllRotamersH98  
    SHOWMESSAGE all the (maxNumberRotamers) raw rotamers. Pruning for uniqueness done for all substrates at the end. 
    Wait ContinueButton
    
  # ensure nothing left. 
  # ------------------------------------------------------------------------------------------------------------------------------
  DELOBJ All

# --------------------------------------------------------------------------------------------------------------------------------
# This is to prune the substrates, this is done for all variants at the end
# --------------------------------------------------------------------------------------------------------------------------------

for i = 1 to count ListMacroTargets  
  # define macrotarget from the list
  MacroTarget = '(ListMacroTargets(i))'

  LOADSCE EA_library/(MacroTarget)/(MacroTarget)_AllRotamers.sce

  # turns out that needed to have a registration of everything that still exists
  TotalNumberOfObjects = COUNTOBJ all
  for j = 1 to TotalNumberOfObjects
    StilExistingObjects(j) = 1
  
  
  # do some test if H98 exists, make copies in which atom names are inverted. 
  localNumberRotamers =  (maxNumberRotamers) / (count DihedralAngleCH)
  PRINT SUMMARIZED  (maxNumberRotamers) (count DihedralAngleCH) (localNumberRotamers) 

  for j = 1 to count DihedralAngleCH  
    firstRotamerCurrentSet = 2 + ( ((j) - 1)* (localNumberRotamers) )  - 1
    lastRotamerCurrentSet  = (firstRotamerCurrentSet) + (localNumberRotamers) - 1
    PRINT SUMMARIZED (firstRotamerCurrentSet) (lastRotamerCurrentSet)
          
    startPointCurrent = (firstRotamerCurrentSet)
    lastPointCurrent  = (lastRotamerCurrentSet)
    SHOWMESSAGE will now compare object (startPointCurrent)  to (lastPointCurrent) 
    CurrentPointComparison        = (startPointCurrent)
    CurrentStartingPointComparing = (startPointCurrent)
    CurrentLastPointComparing     = (lastPointCurrent) 
    FutureLastPointComparing      = (lastPointCurrent)
    while (CurrentStartingPointComparing) < (CurrentLastPointComparing)
      #CurrentlyAccepted = (h) - 1 + (startPointCurrent)
      CurrentLastPointComparing = (FutureLastPointComparing)
      CurrentStartingPointComparing = (CurrentStartingPointComparing) + 1
      # only continue if the comparison part exists
      PRINT Currently looking at (CurrentStartingPointComparing)
      if StilExistingObjects(CurrentPointComparison)
        for q = (CurrentStartingPointComparing) to (CurrentLastPointComparing)
          # verify the other compared partner still exist
          
          if StilExistingObjects(q)
            # Note, the following does not compare the hydroxyl hydrogens. Thus, this could eliminate some relevant rotamers
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            TestRMSD = RMSDATOM OBJ (q) atom !(AtomsToBeFrozen) element !H, OBJ (CurrentPointComparison) atom !(AtomsToBeFrozen) element !H
            TestRMSD = 0.000000 + (TestRMSD) 
            # added 0.00000 test because it turned out to be necessary. Weird. Bug in Yasara?
            if (TestRMSD) == 0.000000
              DELOBJ (q)
              PRINT deleted (q)
              StilExistingObjects(q) = 0 
            else
              if  TestRMSD < Prune 
                DELOBJ (q)
                PRINT deleted (q)
                StilExistingObjects(q) = 0 

  TotalNumberOfObjects = COUNTOBJ all
  SHOWMESSAGE after pruning (TotalNumberOfObjects) left of (MacroTarget)

  if OKEveryStep == 'Yes'
    WAIT ContinueButton
  else 
    Wait 3, unit=seconds  
  PruneScaled = 10 * (Prune)
  RENUMBEROBJ all
  SAVESCE EA_library/(MacroTarget)/(MacroTarget)_pruned(PruneScaled)TenthAngstrom_(TotalNumberOfObjects)RotamersLeftAfterPruningForRMSD   

  CELL AUTO, extension = 5
  for j = 1 to TotalNumberOfObjects
    TABULATE (j)
    REMOVEOBJ all
    ADDOBJ (j)
    TABULATE
    ENERGYRES all, Bond , Angle , Dihedral , Planarity , VdW , Coulomb
    EnergyTotal(j) = ENERGYRES all
    TABULATE (EnergyTotal(j))
  SAVETAB 1,  EA_library/(MacroTarget)/(MacroTarget)_AllRotamers_energiesAfterPruningForRMSD, columns = 8 
  DELTAB 1
  # delete the cell again
  DELOBJ ((TotalNumberOfObjects) + 1)
  ADDOBJ all

  MinumumEnergy = (min EnergyTotal)
  for j = 1 to TotalNumberOfObjects
    if (EnergyTotal(j)) > ((MinumumEnergy) +(MaxEaboveMin))
      DELOBJ (j)
      PRINT deleted object (j) because its energy was (0+(EnergyTotal(j))), which is higher than the threshold of  (0.00+((MinumumEnergy) +(MaxEaboveMin)))
  TotalNumberOfObjects = COUNTOBJ all
  DELVAR EnergyTotal # Yasara does not create new lists otherwise
  SHOWMESSAGE after pruning for energy (TotalNumberOfObjects) left of (MacroTarget)
  RENUMBEROBJ all
  SAVESCE EA_library/(MacroTarget)/(MacroTarget)_pruned(PruneScaled)TenthAngstrom_And(MaxEaboveMin)kJ(TotalNumberOfObjects)RotamersLeftAfterPruningForRMSDandForEnergy  
  SAVESCE EA_library/(MacroTarget)/(MacroTarget)_pruned(PruneScaled)TenthAngstrom(TotalNumberOfObjects)_RotamersLeft
  
  # see if it is a primary amine, then also prune the other rotamers.
  PrimaryAmine = COUNTATOM H98
  if PrimaryAmine
    DELOBJ All
    LOADSCE EA_library/(MacroTarget)/(MacroTarget)_AllRotamersH98.sce

    # do some test if H98 exists, make copies in which atom names are inverted. 

    i = 2
    TestIfFinishedPruning = COUNTOBJ All
    while (i) < (TestIfFinishedPruning)
      CurrentlyAccepted = (i) - 1
      for j = i to TestIfFinishedPruning
        # Note, the following does not compare the hydroxyl hydrogens. Thus, this could eliminate some relevant rotamers
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        TestRMSD = RMSDATOM OBJ (j) atom !(AtomsToBeFrozen) element !H, OBJ (CurrentlyAccepted) atom !(AtomsToBeFrozen) element !H
        TestRMSD = 0.000000 +(TestRMSD)
        # PRINT (TestRMSD) (Prune) (j) (i)
        # added 0.00000 test because it turned out to be necessary. Weird. Bug in Yasara?
        if (TestRMSD) == 0.000000
          DELOBJ (j)
          #PRINT deleted it        
        else
          if TestRMSD < (0.00000 + (Prune))
            DELOBJ (j)
            #PRINT deleted it
      RENUMBEROBJ ALL   
      TestIfFinishedPruning = COUNTOBJ All 
      i = i + 1
 
    NumberOfRotamers = COUNTOBJ all
    SHOWMESSAGE after pruning (NumberOfRotamers) left of (MacroTarget) H98

    if OKEveryStep == 'Yes'
      WAIT ContinueButton
    else
      Wait 3, unit=seconds  
    SAVESCE EA_library/(MacroTarget)/(MacroTarget)_pruned(PruneScaled)TenthAngstrom_(TestIfFinishedPruning)RotamersLeftH98   
  
  # ensure nothing left. 
  # ------------------------------------------------------------------------------------------------------------------------------
  DELOBJ All  

  # Convert the generated SCE file into PDBs to be used by Rosetta for docking
  LOADSCE EA_library/(MacroTarget)/(MacroTarget)_pruned(PruneScaled)TenthAngstrom(TotalNumberOfObjects)_RotamersLeft 
  DELATOM ATOM H99 
  numberObj=COUNTOBJ ALL 
  SHELL rm rot_*.pdb
  SAVEMOL2 1, EA_library/(MacroTarget)/SUB_for_params.mol2
  for i=1 to numberObj
    SAVEPDB (i), rot_(i).pdb
  SHELL cat rot_*pdb | grep -e "END" -e "HETATM" | grep -v "REMARK" | sed 's/1.00  0.00/1.00 20.00/g' | sed 's/SUB X 999/SUB X   1/g' > EA_library/(MacroTarget)/SUB_rotamers.pdb
  SHELL cat rot_1.pdb | grep -e "END" -e "HETATM" | grep -v "REMARK" | sed 's/1.00  0.00/1.00 20.00/g' | sed 's/SUB X 999/SUB X   1/g' > EA_library/(MacroTarget)/SUB_A_0001.pdb
  SHELL rm rot_*.pdb
  SHOWMESSAGE Rotamer library has been printed to SUB_rotamers.pdb
  WAIT CONTINUEBUTTON 
  
  # Generate *params file to be used by Rosetta docking
  if GenerateParamsFile == 'YES'
    print "generating SUB.params file for (MacroTarget)"
    SHELL (molfiletoparams) -n SUB --keep-names EA_library/(MacroTarget)/SUB_for_params.mol2 --clobber
    SHELL mv SUB_0001.pdb EA_library/(MacroTarget)/
    SHELL mv SUB.params EA_library/(MacroTarget)/
    SHELL echo "PDB_ROTAMERS SUB_rotamers.pdb" >> EA_library/(MacroTarget)/SUB.params
  

# Now stop the program or mention that it is finished, depending on how the user decided above
# --------------------------------------------------------------------------------------------------------------------------------

if AutomaticExitAtEnd == 'YES'
  SHOWMESSAGE Finished, pressing the continue button will end this Yasara session
  Wait ContinueButton
  EXIT
else 
  SHOWMESSAGE Finished Preparing Substrates

  
