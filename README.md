#####################################################

This is the summary of our implementation of the AMBER99SB-UCB force field for GROMACS.

https://doi.org/10.1021/ct2000183
https://doi.org/10.1021/jp2118373

are the two reference papers for this force field. Please cite them when you are using our implementation.

This particular implementation is done by

Mohammed Khaled, Batuhan Kav, and Birgit Strodel at the Forschungszentrum Juelich, Germany.

For all inquiries, problems, and error reports, please contact b.strodel [at] fz-juelich.de

####################################################



First: add the new atoms HQ,NQ,and OR

edit aminoacids.rtp
1)add HQ (primery amide) for ASN, GLN, and LYN(neutral LYS)

ASN
HD21 H change to HD21 HQ
HD22 H change to HD22 HQ
NASN
HD21 H change to HD21 HQ
HD22 H change to HD22 HQ
CASN
HD21 H change to HD21 HQ
HD22 H change to HD22 HQ

GLN
HE21 H change to HE21  HQ
HE22 H change to HE22  HQ
NGLN
HE21 H change to HE21  HQ
HE22 H change to HE22  HQ 
CGLN
HE21 H change to HE21  HQ
HE22 H change to HE22  HQ

LYN (neutral LYS)
HZ1   H change to HZ1   HQ
HZ2   H change to HZ2   HQ


2)add NQ (ASN, GLN and neutral LYS)  #DONE!

ASN
ND2 N change to ND2 NQ 
NASN
ND2 N change to ND2 NQ
CASN
ND2 N change to ND2 NQ


GLN
NE2 N change to NE2 NQ
NGLN
NE2 N change to NE2 NQ
CGLN
NE2 N change to NE2 NQ

LYN (neutral LYS)
NZ    N change to NZ    NQ

3) add OR (phenols(TYR))   #DONE!

TYR
OH OH change to OH OR
CTYR
OH OH change to OH OR
NTYR
OH OH change to OH OR

#########################################################################
edit atomtype.atp     #DONE!
add:
HQ                 1.00800  ; new amber99sb-ucb,ASN,GLN,LYN(neutral LYS) 
NQ                14.01000  ; new amber99sb-ucb,ASN,GLN,LYN(neutral LYS) 
OR                16.00000  ; new amber99sb-ucb,TYR 

#########################################################################
edit ffbonded.itp

First: edit phi dihedral angle potenatil term n=2: 8.36800 change to 7.53120       
CT  CT  N   C     9       0.0      7.53120     2  ; new amber99sb-ucb(8.36800 for amber99sb)

*************************************************************************
Second:We need to define the new bonds,angles,proper dihederals and impropers dihederals for ASN, GLN, LYN and TYR from aminoacids.rtp and then add them to  ffbonded.itp

######### New bonds ##############
1)ASN
CG   ND2: define new bond C NQ ,same bonded parameters as C N bond
ND2  HD21 and ND2  HD22: define new bond NQ HQ ,same bonded parameters as N H bond

2)GlN 
CD   NE2:define new C NQ, same bonded parameters as C N
NE2  HE21 and NE2  HE22: define new bond NQ HQ,same bonded parameters as N H

3)LYN 
CE NZ: define new bond CT NQ, same bonded parameters as CT N3
NZ HZ1 and NZ HZ2: define new bond NQ HQ, same bonded parameters as N H  

4)TYR 
CZ OH change to C OR, same bonded parameters as C OH bond
OH HH change to OR HO,same bonded parameters as OH HO bond 
*****************************************************************************
add new bonds to ffbonded.itp: #DONE!    

C  OR         1    0.13640   376560.0 ; new amber99sb-ucb,TYR
HO OR         1    0.09600   462750.4 ; new amber99sb-ucb,TYR
C  NQ         1    0.13350   410032.0 ; new amber99sb-ucb,ASN,GLN
CT NQ         1    0.14710   307105.6 ; new amber99sb-ucb,LYN(neutral LYS)
HQ NQ         1    0.10100   363171.2 ; new amber99sb-ucb,ASN,GLN,LYN(neutral LYS) 
*************************************************************************
######### New Angles ##############
new angles: 

1) ASN and GLN
C   NQ  HQ           1   120.000    418.400 ; new amber99sb-ucb,ASN,GLN /same angle parameters as: C   N   H            1   120.000    418.400 ; new99 general, gln, asn,changed based on NMA nmodes
HQ  NQ  HQ           1   120.000    292.880 ; new amber99sb-ucb,ASN,GLN /same angle parameters as: H   N   H            1   120.000    292.880 ;ade,cyt,gua,gln,asn
CT  C   NQ           1   116.600    585.760 ; new amber99sb-ucb,ASN,GLN /same angle parameters as: CT  C   N            1   116.600    585.760 ; AA general
NQ  C   O            1   122.900    669.440 ; new amber99sb-ucb,ASN,GLN /same angle parameters as: N   C   O            1   122.900    669.440 ; AA general

2)LYN (neutral LYS)
HP  CT  NQ           1   109.500    418.400 ; new amber99sb-ucb,LYN(neutral LYS) /same angle parameters as: HP  CT   N           1   109.500    418.400 ; AA lys, ch3nh3+, changed based on NMA nmodes
CT  CT  NQ           1   111.200    669.440 ; new amber99sb-ucb,LYN(neutral LYS) /same angle parameters as: CT  CT   N           1   111.200    669.440 ; AA lys             (JCP 76, 1439)
CT  NQ  HQ           1   109.500    418.400 ; new amber99sb-ucb,LYN(neutral LYS) /same angle parameters as: CT  N    H           1   109.500    418.400 ; AA lys,     changed based on NMA nmodes
HQ  NQ  HQ           1   120.000    292.880 ; new amber99sb-ucb,LYN(neutral LYS) /same angle parameters as: H   N    H           1   109.500    292.880 ; AA lys, AA(end)

3)TYR
CA  C   OR           1   120.000    585.760 ; new amber99sb-ucb,TYR /same angle parameters as: CA  C   OH           1   120.000    585.760 ; AA tyr
C   OR  HO           1   113.000    418.400 ; new amber99sb-ucb,TYR /same angle parameters as: C   OH  HO           1   113.000    418.400 ; new99

***************************************************************************
add new angles to ffbonded.itp: #DONE! 

C   NQ  HQ           1   120.000    418.400 ; new amber99sb-ucb,ASN,GLN 
CT  C   NQ           1   116.600    585.760 ; new amber99sb-ucb,ASN,GLN   
NQ  C   O            1   122.900    669.440 ; new amber99sb-ucb,ASN,GLN
HQ  NQ  HQ           1   120.000    292.880 ; new amber99sb-ucb,ASN,GLN,LYN(neutral LYS) 
HP  CT  NQ           1   109.500    418.400 ; new amber99sb-ucb,LYN(neutral LYS) 
CT  CT  NQ           1   111.200    669.440 ; new amber99sb-ucb,LYN(neutral LYS)
CT  NQ  HQ           1   109.500    418.400 ; new amber99sb-ucb,LYN(neutral LYS)  
CA  C   OR           1   120.000    585.760 ; new amber99sb-ucb,TYR
C   OR  HO           1   113.000    418.400 ; new amber99sb-ucb,TYR  

**********************************************************************************
######### New propers dihederals  ##############
New propers dihederals:

1)ASN and GLN 
CT C NQ HQ define as: X   C   NQ   X     9     180.0     10.46000     2  ; new amber99sb-ucb,ASN,GLN /same propers dih parameters as: X   C   N   X     9     180.0     10.46000     2  ; AA,NMA
O  C NQ HQ define as: X   C   NQ   X     9     180.0     10.46000     2  ; new amber99sb-ucb,ASN,GLN /same propers dih parameters as: X   C   N   X     9     180.0     10.46000     2  ; AA,NMA 

2)LYN (neutral LYS)
CT CT NQ HQ define as: X   CT  NQ  X     9       0.0      0.65084     3  ; new amber99sb-ucb,LYN /same propers dih parameters as: X   CT  N3  X     9       0.0      0.65084     3  ; JCC,7,(1986),230
HP CT NQ HQ define as: X   CT  NQ  X     9       0.0      0.65084     3  ; new amber99sb-ucb,LYN /same propers dih parameters as: X   CT  N3  X     9       0.0      0.65084     3  ; JCC,7,(1986),230

3)TYR
CA C  OR HO define as: X   C   OR  X     9     180.0      9.62320     2 ; new amber99sb-ucb,TYR /same propers dih parameters as: X   C   OH  X     9     180.0      9.62320     2  ; Junmei et al, 1999

**************************************************************************************
add new propers dihederals to ffbonded.itp: #DONE!

X   C   NQ  X     9     180.0     10.46000     2  ; new amber99sb-ucb,ASN,GLN 
X   C   NQ  X     9     180.0     10.46000     2  ; new amber99sb-ucb,ASN,GLN 
X   CT  NQ  X     9       0.0      0.65084     3  ; new amber99sb-ucb,LYN(neutral LYS) 
X   C   OR  X     9     180.0      9.62320     2  ; new amber99sb-ucb,TYR 

************************************************************************************
######### New impropers dihederals  ##############

1)ASN and GLN:(form aminoasids.rtp file)
impropers for ASN form aminoasids.rtp file:  
CB   ND2    CG   OD1: define new improper dih CT NQ C O, same impropers dih parameters as CT N C O
CG  HD21   ND2  HD22: define new improper dih C HQ NQ HQ, same impropers dih parameters as C H N H 
impropers for GLN form aminoasids.rtp file:   
CG   NE2    CD   OE1: define new improper dih CT NQ C O, same impropers dih parameters as CT N C O 
CD  HE21   NE2  HE22: define new improper dih C HQ NQ HQ, same impropers dih parameters as C H N H 

improper dih CT N C O in ffbonded.itp define as: X   X   C   O        4      180.00    43.93200     2	; JCC,7,(1986),230 /## no need to add it ##
improper dih C  H N H in ffbonded.itp define as: X   X   N   H        4      180.00     4.18400     2	; JCC,7,(1986),230

2)TYR
CA  CA  C   OR       4      180.00     4.60240     2    ; new amber99sb-ucb,TYR /same impropers parameters as: CA  CA  C   OH       4      180.00     4.60240     2

*************************************************************************************
add new impropers dih to ffbonded.itp: #DONE! 
	
X   X   NQ  HQ       4      180.00     4.18400     2	; new amber99sb-ucb,ASN,GLN
CA  CA  C   OR       4      180.00     4.60240     2    ; new amber99sb-ucb,TYR 

##################################################################################
edit ffbonded.itp (added to avoid errors in gmx grompp)

HQ           1       1.008   0.0000  A   1.06908e-01  6.56888e-02 ; new amber99sb-ucb;ASN,GLN,LYN(neutral LYS)  /same parameters as H
NQ           7      14.01    0.0000  A   3.25000e-01  7.11280e-01 ; new amber99sb-ucb;ASN,GLN,LYN(neutral LYS)  /same parameters as N 
OR           8      16.00    0.0000  A   3.06647e-01  8.80314e-01 ; new amber99sb-ucb;TYR /same parameters as OH

################################################################################
create new file nbfix.itp

##################################################################################
########## Partial atomic charges #########
change the partial atomic charges for C-terminus non-protonated ALA, GLY, and VAL. The data from SI Table 8 (https://doi.org/10.1021/ct2000183):

aminoacids.rtp

1) CALA 
 [ atoms ]
     N    N           -0.43970     1 ; new amber99sb-ucb
     H    H            0.26110     2 ; new amber99sb-ucb
    CA    CT           0.06060     3 ; new amber99sb-ucb
    HA    H1           0.08790     4 ; new amber99sb-ucb
    CB    CT          -0.10780     5 ; new amber99sb-ucb
   HB1    HC           0.04180     6 ; new amber99sb-ucb
   HB2    HC           0.04180     7 ; new amber99sb-ucb
   HB3    HC           0.04180     8 ; new amber99sb-ucb
     C    C            0.83290     9 ; new amber99sb-ucb
   OC1    O2          -0.60740    10 ; new amber99sb-ucb
   OC2    O2          -0.66320    11 ; new amber99sb-ucb

2) CGLY
 [ atoms ]
     N    N           -0.37180     1 ; new amber99sb-ucb
     H    H            0.25810     2 ; new amber99sb-ucb
    CA    CT          -0.11370     3 ; new amber99sb-ucb
   HA1    H1           0.09270     4 ; new amber99sb-ucb
   HA2    H1           0.09270     5 ; new amber99sb-ucb
     C    C            0.87960     6 ; new amber99sb-ucb
   OC1    O2          -0.61740     7 ; new amber99sb-ucb
   OC2    O2          -0.66550     8 ; new amber99sb-ucb

3) CVAL 
 [ atoms ]
     N    N           -0.33850     1 ; new amber99sb-ucb
     H    H            0.22970     2 ; new amber99sb-ucb
    CA    CT          -0.09100     3 ; new amber99sb-ucb
    HA    H1           0.10340     4 ; new amber99sb-ucb
    CB    CT           0.20170     5 ; new amber99sb-ucb
    HB    HC          -0.00380     6 ; new amber99sb-ucb
   CG1    CT          -0.19790     7 ; new amber99sb-ucb
  HG11    HC           0.04710     8 ; new amber99sb-ucb
  HG12    HC           0.04710     9 ; new amber99sb-ucb
  HG13    HC           0.04710    10 ; new amber99sb-ucb
   CG2    CT          -0.19790    11 ; new amber99sb-ucb
  HG21    HC           0.04710    12 ; new amber99sb-ucb
  HG22    HC           0.04710    13 ; new amber99sb-ucb
  HG23    HC           0.04710    14 ; new amber99sb-ucb
     C    C            0.83150    15 ; new amber99sb-ucb
   OC1    O2          -0.60710    16 ; new amber99sb-ucb
   OC2    O2          -0.66310    17 ; new amber99sb-ucb


Reverting the partial charges of CVAL, CGLY, CALA back to the original A99disp, as they are not protonated and their charges are not supposed to be changed.

Changed O2 atom type in CVAH, CGLH, and CALH to O. Apparently this is the correct version.


#############################################################################################################################################################
