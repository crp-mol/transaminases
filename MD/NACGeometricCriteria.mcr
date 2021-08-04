# This file contains the NAC definitions. It is called on-the-fly by: MDSimulation.mcr
# --------------------------------------------------------------------------------------------------------------------------------------------------------------
# The NAC definitions can be altered accordingly.
# Warning: Some geometric criteria are not turned off in practice
# For example, Asp256A OD2-H1 will count a NAC whenever the distance is lower than 20 Angstrom, which is always true. Hence, Asp256A OD2-H1 is, in practice, turned off
# The NAC definitions of this file are for the binding site A (the one made up mostly of residues from subunit A). But this can be changed easily.
# Choose the enzyme system you're working with.
if CT == 'AT_vibriofluvialis'          # user-defined settings name from MDSimulation.mcr 
  ActiveSiteResidues = 'Lys 285'
  StayAwayDistance = 15
  pH = 7.0
  
  # S means sets
  S() = 'A'
  # This is intended for testing if multiple NACs/productive conformations are present at the same time, should be handy in case of salt-bridge networks. 
  SetsTestCombined_A() = 'A'

  # L means Letter
  L() ='A','B'
  AntiL()='B','A' # For residues of subB that are in the binding site of subA, or viceversa.
  # T means target atom
  TC() = 'C99'
  TH() = 'H99'
  # R means reference atoms
  R1() = 'O1'
  R2() = 'N50'
  R3() = 'C50'
  R4() = 'N99'
    # Atoms from lysine 285
  R5() = 'NZ'
  R6() = 'CE'
  R7() = 'HZ'
  
  #,'C98'
  for i = 1 to count L
    # The labels, measures, and minimal and maximal criteria (if no min or max criteria, use a value that is never reached)
    # i=1,2,3,.... # Only L1, L2, L3, L5, and L6 were used in the manuscript as NAC definitions. But you might find other measurements to be useful

    # first for the distance between the nitrogen and the hydrogen atom that is abstracted
    # Called d_1 in the manuscript
    Label_(L(i))_1		       = '(L(i))_Distan_N_H'
    (Label_(L(i))_1)_Measurement       = 'DISTANCE MOL (L(i)) res SUB atom (TH(1)), MOL (L(i)) res Lys 285 atom NZ with minimum distance from OBJ 1 res SUB atom (TH(1))'
    (Label_(L(i))_1)_Q_Max	       = 2.65
    (Label_(L(i))_1)_Q_Min	       = -0.01
    
    # Then for dihedral to check that the ring system is flat.
    # Called Chi_2 in the manuscript
    Label_(L(i))_2		       = '(L(i))_DihedralAll'
    (Label_(L(i))_2)_Measurement       = 'Dihedral MOL (L(i))  RES SUB atom (R1(1)), MOL (L(i))  RES SUB atom (R2(1)), MOL (L(i)) res SUB atom (R3(1)),MOL (L(i))  Res SUB atom (R4(1))'
    (Label_(L(i))_2)_Q_Max	       = 15.0
    (Label_(L(i))_2)_Q_Min	       = -15.0   
   
    # Then the groupangle between all the atoms of the pi system and the C99 H99 direction. 
    # Called Chi_1 in the manuscript
    Label_(L(i))_3		       = '(L(i))_GroupAngle'
    (Label_(L(i))_3)_Measurement       = 'Dihedral MOL (L(i)) res SUB atom C50, MOL (L(i)) res SUB atom N99, MOL (L(i)) res SUB atom C99, MOL (L(i)) res SUB atom H99'
    (Label_(L(i))_3)_Q_Max	       = -75.0
    (Label_(L(i))_3)_Q_Min	       = -105.0   

# the following geometric criteria are currently not turned on, but it can be changed very easily

    # There needs to be a catalytic orientation of the Lys-NH2 for proton abstraction (nuclephilic/electrophilic attack). Angle ~ 109.5deg
    # Angle H99 - Lys(N) - Lys(CE)
    Label_(L(i))_4		       = '(L(i))_H99_N_C'
    (Label_(L(i))_4)_Measurement       = 'Angle mol (L(i)) RES SUB atom (TH(1)), mol (L(i)) RES 285 atom (R5(1)), mol (L(i)) RES 285 atom (R6(1)) with maximum distance from mol (L(i)) res sub atom (TH(1))'
    (Label_(L(i))_4)_Q_Max             = 180
    (Label_(L(i))_4)_Q_Min             = 0

    # Angle H99 - Lys(N) - Lys(HZ(1)) # Hydrogen from lysine closer to H99
    # Called theta_1 in the manuscript
    Label_(L(i))_5		       = '(L(i))_H99NH1_Angle'
    (Label_(L(i))_5)_Measurement       = 'Angle mol (L(i)) RES SUB atom (TH(1)), mol (L(i)) RES 285 atom (R5(1)), mol (L(i)) RES 285 atom (R7(1)) with minimum distance from mol (L(i)) res sub atom (TH(1))'
    (Label_(L(i))_5)_Q_Max             = 180
    (Label_(L(i))_5)_Q_Min             = 0

    # Angle H99 - Lys(N) - Lys(HZ(2)) # Hydrogren from lysine further away from H99
    # Called theta_2 in the manuscript
    Label_(L(i))_6		       = '(L(i))_H99NH2_Angle'
    (Label_(L(i))_6)_Measurement       = 'Angle mol (L(i)) RES SUB atom (TH(1)), mol (L(i)) RES 285 atom (R5(1)), mol (L(i)) RES 285 atom (R7(1)) with maximum distance from mol (L(i)) res sub atom (TH(1))'
    (Label_(L(i))_6)_Q_Max             = 180
    (Label_(L(i))_6)_Q_Min             = 0
  
    # C99-H99-NZ-(Lys) angle must be around 180 deg for proton abstraction
    Label_(L(i))_7		       = '(L(i))_C99H99NZ_Angle'
    (Label_(L(i))_7)_Measurement       = 'Angle mol (L(i)) RES SUB atom (TC(1)), mol (L(i)) RES SUB atom (TH(1)), mol (L(i)) RES 285 atom (R5(1))'
    (Label_(L(i))_7)_Q_Max             = 180
    (Label_(L(i))_7)_Q_Min             = 0

    # Thr(322)-O --- H-(N-Lys) hydrogen bond, distance less than 2.5
    # Thr322 belongs to the other subunit
    Label_(L(i))_8		       = '(L(i))_Distan_THR_LYS'
    (Label_(L(i))_8)_Measurement       = 'DISTANCE res 322 atom OG1 with minimum distance from MOL (L(i)) res 285 atom HZ, MOL (L(i)) res 285 atom HZ with minimum distance from res 322 atom OG1'
    (Label_(L(i))_8)_Q_Max	       = 20
    (Label_(L(i))_8)_Q_Min	       = -0.01

    # Asp256A OD2-H1
    # Simulations substrates 21-60 had a typo. Atom is now called HN5, not H1.
    Label_(L(i))_9		       = '(L(i))_ASPHN5_Dist'
    (Label_(L(i))_9)_Measurement       = 'DISTANCE MOL (L(i)) res ASP 256 element O with minimum distance from mol (L(i)) res SUB atom HN5, mol (L(i)) res SUB atom HN5'
    (Label_(L(i))_9)_Q_Max	       = 20
    (Label_(L(i))_9)_Q_Min	       = -0.01

    # Alternative dihedral
    # Ideal -90 deg
    Label_(L(i))_10		       = '(L(i))_C7C3N99H99_Dih'
    (Label_(L(i))_10)_Measurement      = 'Dihedral MOL (L(i)) res SUB atom C7, MOL (L(i)) res SUB atom C3, MOL (L(i)) res SUB atom N99, MOL (L(i)) res SUB atom H99'
    (Label_(L(i))_10)_Q_Max	       = 360
    (Label_(L(i))_10)_Q_Min	       = -360

    # Label 11 to 14 are about the water possibly found in between Thr322 and Lys285, we dont' know if it's there.
    # Check if there is any HOH bound in-between THR322 and Lys285
    # HZ from Lys285A would be bound to the oxygen of water
    Label_(L(i))_11		       = '(L(i))_LYSHOH_Dist'
    (Label_(L(i))_11)_Measurement      = 'DISTANCE atom O res HOH with minimum distance from mol (L(i)) res 285 atom HZ, mol (L(i)) res 285 atom HZ with minimum distance from res HOH atom O'
    (Label_(L(i))_11)_Q_Max	       = 20
    (Label_(L(i))_11)_Q_Min	       = -0.01
    
    # OG1 from THR322B would be bound to the hydrogen of the water
    Label_(L(i))_12		       = '(L(i))_THRHOH_Dist'
    (Label_(L(i))_12)_Measurement       = 'DISTANCE atom H res HOH with minimum distance from mol (AntiL(i)) res 322 atom OG1, mol (AntiL(i)) res 322 atom OG1'
    (Label_(L(i))_12)_Q_Max	       = 100
    (Label_(L(i))_12)_Q_Min	       = -0.01

    # Determine whether the water measured in the last two "Label" is the same, that means Label13 = Label14
    # Label13 and Label14 will also give hint on the location of the water in question (angle closer to 180 means that the water is in-between)
    Label_(L(i))_13		       = '(L(i))_LysO_Thr'
    (Label_(L(i))_13)_Measurement      = 'Angle mol (L(i)) res 285 atom NZ, res HOH atom O with minimum distance from mol (L(i)) res 285 atom NZ, mol (AntiL(i)) res 322 atom OG1'
    (Label_(L(i))_13)_Q_Max	       = 180
    (Label_(L(i))_13)_Q_Min	       = 0

    Label_(L(i))_14		       = '(L(i))_Lys_OThr'
    (Label_(L(i))_14)_Measurement      = 'Angle mol (L(i)) res 285 atom NZ, res HOH atom O with minimum distance from mol (AntiL(i)) res 322 atom OG1, mol (AntiL(i)) res 322 atom OG1'
    (Label_(L(i))_14)_Q_Max	       = 180
    (Label_(L(i))_14)_Q_Min	       = 0

