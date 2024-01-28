!===================================================================================================================================
! Single file version of HFBTHO 
! MODULES:
!  Module HFBTHO_VERSION  
!  Module UNEDF               (hidden if LDflags   = -Dhide_dme) 
!  Module HFBTHO
!  Module pairing_HFBTHO 
!  Module THO_MODULE          (hidden if LDflags   = -Dhide_tho)
!  Module qrpa_HFBTHO         (hidden if LDflags   = -Dhide_qrpa)
! MAIN PROGRAM
!  Program hfbthoprog         (hidden if LDflags   = -Dhide_mainprog)
! SOLVER
!  Subroutine HFBTHO_SOLVER 
!===================================================================================================================================
!
!===================================================================================================================================
!#START HFBTHO_VERSION MODULE
!===================================================================================================================================
Module HFBTHO_VERSION
  Implicit None
  Character(6) :: Version='133'  
  !--------------------------------------------------------------------------------------
  ! Version History
  !--------------------------------------------------------------------------------------
  ! ver#133: single-file HFBTHO
  ! ver#130: tho.dat mdifications due to blocking, error indicator introduced
  ! ver#129: Even-even tested and equivalent with ptho101spt15sp.f90 used in ANL fit
  ! ver#128: EQP,U,V and their dimentions NUV,NEQ required for qrpa incoporated
  !          permanently in HFBTHO substituting old arrays eqp and uv  
  ! ver#127: For easy development the module split in different F90 files which are
  !          invoced using INCLUDE statements (remove that when developement is over)
  !          Optimized qrpa_DENSIT_PLUS and qrpa_GAMDEL to be twice faster    
  ! ver#126: Cleaning, optimizing, and separating THO stuff in one place 
  !             parity good:   Time per iteration: 3.841 seconds
  !             parity broken: Time per iteration: 9.933 seconds
  !          HFB+THO tested and works in both parity regimes. iserial removed and 
  !          substituted with Print_Screen, i.e, record results only when Nsh>0
  ! ver#125: Implemented and tested reflection symmetry as optional. When broken the
  !          computer time per iteration is almost 5 times bigger:
  !             parity good:   Time per iteration:  4.5276 seconds
  !             parity broken: Time per iteration: 19.4646 seconds
  !             the difference in total energiy is  0.001  keV 
  ! ver#124: Rewrite to prepare for braking the reflection symmetry
  !          Preprocessing directives included: #ifndef hide_qrpa, hide_tho, hide_dme 
  !          For preprocessing one needs -Dhide_qrpa -Dhide_tho -Dhide_dme
  !          If no preprocessing or -Uhide_qrpa -Uhide_tho -Uhide_dme then
  !          all modules are included
  ! ver#123: Playground for QRPA calculations [16/12/2010]:
  !          The changes are: 
  !          - all variables in Module HFBTHO are Public 
  !          - include_qrpa=0 is added to Module HFBTHO with asssociated
  !            declarations eventualy used by qrpa
  !          - Subroutine ByNucleus moved to PTHO_PROGRAM where is its place
  !            and it should be done long ago. Call Do_QRPA() is used only there.
  !          - So if the program is compiled with -Uhide_qrpa one can use Do_QRPA()
  !            to do qrpa calculations.  
  ! ver#123: Fixed crash after iterations limit. 
  !          Tested against anl version hfbtho101spt15.f90 - itterations go differently 
  !          but the final results are identical. 
  ! ver#122: (MK) Added CExPar for coulomb exchange. Parameter read from UNEDF module
  ! ver#121: (MK) Added possibility to use zero particle number for droplet calculations
  ! ver#120: (MK) Added external field, and all channels to direct Hartree. e^2 for Coulomb
  !          now read from the UNEDF module. Direct Hartree now always calculated based on
  !          module function regardless of the value of DMEorder parameter
  ! ver#117: Direct Hartree added when DME_order>-1
  ! ver#115: (MK) added use_cm_cor variable to hb0 calculation and (nabla rho)^2 terms to
  !          the calculate_U_parameters function calls
  ! ver#114: Name list, new tho.dat file, proton/neutron fields, confirms all results of 
  !          recent published version after ANL optimization ptho101b_last_tested.f90,
  !          public/Public variables            
  ! ver#113: Cleaning
  ! ver#112: No parameter functions
  ! ver#111: Main program detached from the file as PTHO_MAIN_PROGRAM.f90 which will not be
  !          versioned. ptho becomes jus a HFBTHO module.  Pairing constants V0(2),V1(2)
  !          replaced by CpV0(0:1), CpV1(0:1) coming as public from defined in UNEDF module
  !          Removed dalf and ippforce form the pairing. For compatibility, ippforce
  !          stays in the input file by now but the kind of pairing is given by CpV1 only
  !          Dropped corrections 'ecmcpavpj', 'erotcorrection' which should be added later.
  !          For compatibility inputfile stays the same. Added IDEUB.
  !          THO part in 'densit' (not densitpj), 'gamdel' commented HO/THO for speed
  ! ver#109: All public variables, expectpj works with a jump: not clear how UNEDF can work 
  !          with complex numbers, just skip this part by now bu write results data
  !          New thodefh(iw1)
  ! ver#108: Removed all programs not used in ver#107
  !          expect contains a key DO_FITT: 
  !            =0 calculare energy, delta, def  & rms only 
  !            =1 the same+all integrals for the regression optimization
  !          V0,V1 pairing constant separated for neutrons and protons: v0(2),v1(2)
  !          HFBTHO collected in MODULE HFBTHO
  !          KOP3 removed
  ! ver#107: Towards UNEDF: complete rewrite based on Marcus to include N2LO
  !          LN for ZR110 at prolate solution with SLY4, mixed pairing and tensor terms:
  !          -agreement with previouse Skyrme implemetation to the last significant digit
  !          -agreement with previouse LO+LDA implemetation to the last significant digit
  !          -agreement with previouse LO+CB  implemetation to the last significant digit
  ! ver#106: Towards UNEDF: the standard functional rewritten in terms of UNEDF U-amplitudes
  !          The assumption U=U(tau_0,Delta rho_0,rho_0,rho_1) becomes possible after 
  !          adding Nabla rho_ij terms  (STANDARD FUNCTIONAL ONLY)
  ! ver#105: Towards UNEDF: the standard functional rewritten in terns of UNEDF U-amplitudes
  !          The assumption U=U(rho_0,rho_1) becomes possible after adding Delta rho_ij terms  
  ! ver#104: Broyden improved with linear search at negative curvature
  !          Implemented Agumented Lagrangian method for constraint calculations
  !          Manual blocking included and tested, key: manualBlocking
  ! ver#103: From this version on-no more support for VAP (VAP completely removed)
  !          The whole program in terms of C-parameters (including tenzor terms)
  ! ver#102: The whole program in terms of C-parameters (without tenzor terms)
  ! ver#101: Optimization in terms of nuclear matter: 'FITS' regime
  ! ver#100: Toward isovector pairing following Sagawa and Yamagami
  ! ver# 99: Subroutine HFBiterations. The isotopic line in tho.dat removed.
  !          Subroutines byNucleus, byConstraint, FitPairing, HFBTHO_HFODD isolated 
  !          at the end and could be ported if necessary. skyrme='FITS' assumes the skyrme
  !          parameters as explicitely given. -N00 supresses completely the output and only
  !          hodef.dat and thodef.dat are charged (if iserial=0 even these files are supressed)
  !          HFBTHO_HFODD updated (think further about a constraint in Q2 terms)   
  ! ver# 98: INOUT modified and added interface subroutine HFBTHO_HFODD
  !          Pairing fitted with MIX/(LN-NOLN) for SLY4,SKP,SKM* forces 
  ! ver# 97: Corrected blocking candidates criteria
  ! ver# 96: Extended so term W0,W1 and SKLY4T forces
  ! ver# 95: Pairing regularization removed, linear HFBDIAG mixing for Lambda 
  !          when blocking, DSYEV replaced by the faster DSYEVD, hfbdiag caculates
  !          canonical basis only at the last hfbdiag iteration, expect optimized
  !          new subroutines HFBiterations, FitPairing, byConstraint
  !          Work around a bug in LAPAK related to DSYEVD           
  ! ver# 94: Misprints 
  ! ver# 93: Removed hh and de matrices and related manipulations
  !          Broyden_min now escape maximum and inflex points
  ! ver# 92: Bug in blocking while no pairing cleaned
  ! ver# 91: BLAST  & LAPACK diagonalization 
  ! ver# 90: If applying Broyden method to matrix elements then
  !          at 20 shells the total number of matrix elements
  !          is 2x2x65307=261228 or about 2.1 Mb and if one keeps
  !          8 iterations it will be  about 17 Mb-not too much
  !          This is the only way one can mix Lipkin-Nogami
  !
  !          If one uses the potentials at 30x30 grid points 
  !          the numbers are 8x2x900=14400 or 115 Kb and if one keeps
  !          8 iterations it will be about 1 Mb-much better
  !          but Lipkin-Nogami is out of this scheme (?!)
  !          If one uses densities at 30x30 grid points they 
  !          are 9x2x900-almost like the potential case.
  !          (sent to George)  
  ! ver# 89: Reduced printout (no anymore lprinter)
  !          LN del+ala2 printed during the iterations
  !          Strength in the initial constraint calculations reduced to 0.3, requested
  !          deformation+/-0.3, untill si<1.1
  !          If too slow convergence (1000 iterations) and Lambda>0 interrupt iterations
  !          Odd nucleus right away from the even-even (even) solution
  !          When even solution missing/corrupt (even at inin<0) calculate it first 
  !          and then odd one
  ! ver# 88: Synchronization for the parallel run
  !          FileLabel subroutine added.
  !          Modified inin control 
  !            inin<0: Always start from a file if it exists, not corrupted and correct
  !                    otherwise inin=iabs(inin) and start from scratch
  !            inin>0: Always start from scratch
  !          SCRATCH calculations start with initial 20 constrained iteration if constrain
  !          is not requested (icstr=0). When icstr#0 standard constraint calculations. 
  !          BY CHAIN calculations temporary removed due to blocking complications           
  ! ver# 87: Approximate Blocking keeping time-reversal symmetry to PNP PAV
  ! ver# 86: Approximate Blocking keeping time-reversal symmetry and tested agings HFODD
  ! ver# 85: LN in canonical basis. Benchmark to HFODD 
  ! ver# 84: Testing HFODD LN again HFB-HO
  ! ver# 83: Cleaning, SKM* mixed volume LN pairing fitted
  ! ver# 82: As 81 but prepared for jaguar
  ! ver# 81: Pairing regularization/renormalization. PAV done with unprojected v_k
  !          V0(Nsh=20,pwi=50) fitted for SLY4,SKP,Renormalized,Regularized,Mixed,Volume
  !          Removed delta and gamdel0 completely 
  ! ver# 80: Accuracy for large number of shells increased by the number of gauss points
  !          Gaussian points now calculated
  !          Initial guest now from deformed Wood-Saxon
  !          Initial run now starts with requested shell number n00
  ! ver# 79: Cranking rotational correction implemented:
  !          Printed to screen, thoout.dat, hodef.dat and thodef.dat 
  !          but not added to the energy
  ! ver# 78: Full CM correction implemented in HFB  & HFB(PAV) 
  !          Printed to screen, thoout.dat, hodef.dat and thodef.dat 
  !          but not added to the energy
  !          NB : ilpjnp(2) removed
  ! ver# 77: Automated Blocking: 
  !          First is calculated N,Z without blocking, remembered *.hel *.tel files and
  !           determined blocking candidates according to pwablo criteria. 
  !           if N(Z) is odd we have neutron(proton) blocking candidates.
  !           if both, N and Z, are odd we have both, neutron and proton, block.candidates
  !          Then we block state after state among the blocking candidates and calculate 
  !           starting from the recorded unblocked (N,Z) solution 
  !           if only N ( or Z) is odd then all neutron (or proton) blocking candidates
  !           are calculated
  !           if both N and Z are odd all pairs of proton and neutron blocking candidates 
  !           are calculated                                               (PAV  &LN unclear)
  ! ver# 76: Manual Blocking for a particular state in a particular block 
  !          overlap criteria used to avoid the level crossing problem    (PAV  &LN unclear)
  ! ver# 75: Manual Blocking for the minimal qpe within a given block. 
  ! ver# 74: Bulgac procedure .. not done
  ! ver# 73: Cleanup, introducing the cpc notations, beyond unit circle removed (MARK 1)
  ! ver# 72: PNP: still valid version for integration over the unit circle
  ! ver# 71: PNP: towards beyond unit circle integration 
  ! ver# 70: PNP: detached neutron from proton projection
  ! ver# 69: PNP: quadrupole constraint fixed to converge
  ! ver# 68: equivalent to var#67
  ! ver# 67: PNP: corrections to the tensor term and initial dumping factor
  ! ver# 66: PNP: V0 fitted to PLN energies at Sn126 for SKP mixed and volume at Nsh=20,HO 
  ! ver# 65: deformed HO basis implemented and tested
  ! ver# 64: fixed byChain to go not 2 beyond the forced break
  ! ver# 63: iasswrong(3) fixed for correct multiprocessor run
  ! ver# 61/62: pthotop for Cheetah added at the end
  ! ver# 60: Thodef.dat header line fixed (added U:). Fixed P/N in ByChain calculations
  !          The 'Stop' is removed from  mishmatch conditions with alternative to use old one.
  !          Consistent pairing for SLY4 and SKP forces
  ! ver# 59: LST modified to accept negative aa-values. SLY4 And SKP with consistent 
  !          (high densities regime) pairing for all cases. Old asymptotic prescription 
  !          is used in the case of Mishmatch asymptotic. Temporary,still new SLY4 
  !          pairing constants are commented.   ! ver# 58: Proton line in byChain 
  !          calculations goes vertically. In case of wrong 
  !          asymptotic parameter 'kindhfb' is recorded as 'kindhfb+100' in thodef.dat 
  !          file where the results for this nucleus are substituted with HO results
  !          Only Nsh=20 pairing constants are already fitted to the higher density  
  !          asymptotic prescription which is already enforced. (Temporary, still old 
  !          SLY4 pairing constants are in the code).  
  ! ver# 57: Partially refitted pairing constants according to the new asymptotic 
  !          prescription. NB! old constants are still for the SLY4 force and not all of 
  !          the cases with SKP are fitted. Old ass. regime is temporary enforced.
  ! ver# 56: Code optimization and checks
  ! ver# 55: Tensor term J.J implemented and tested 
  ! ver# 54: Back to *.hel *.tel files; Including new hodef.dat like thodef.dat file
  ! ver# 53: LST is choosing the higher density in the asymptotic region
  !=======================================================================================
End Module HFBTHO_VERSION
!===================================================================================================================================
!#END HFBTHO_VERSION MODULE
!===================================================================================================================================
!
!===================================================================================================================================
!#START UNEDF MODULE
!===================================================================================================================================
Module UNEDF 
  !=================================================================== 
  ! M.Kortelainen & M.Stoitsov, 2009-2011
  ! UNEDF interface for Skyrme, DME(LO,NLO,N2LO) and other DFT solvers 
  !=================================================================== 
  Implicit None
  !
  Character(16), Private :: Version='17'
  ! 
  ! Version History
  !-------------------------------------------------------------------
  ! ver#17:(Mario)   use_TMR_pairing=0/1 standard/TMR pairing added 
  !                  to Namelist. Using: 
  !                  CpV0(0)=G,    CpV0(1)=a
  !                  CpV1(0)=vfacn,CpV1(1)=vfacp
  ! ver#16:(Mario)   #ifndef hide_dme preprocessing directive included    
  ! ver#15:(Markus)  Added parameter CExPar, used in Coul. excange term.
  !                  Also, all the channels included in direct Hartree
  ! ver#14:(Markus)  Added function Vexternal for the external field,
  !                  and use_j2terms to switch off tensor terms. 
  !                  Direct Hartree set to zero.
  ! Ver#13:(Mario)   Added ac2,ac3,acoord
  ! ver#12:(Mario)   hartree term temprorary dropped. rDr NNN terms taken
  !                  with a factor of 1/2 
  ! ver#11:(Mario)   Gaussian approximation to the Hartree term added,
  ! [3/10/2010]      hatree_rc removed. NB! Function HartreeV is an
  !                  elemental function with possible array arguments
  ! ver#10: (Markus) Added e2charg (e^2 for Coulomb) to the public variables
  ! ver#9: (Mario)   Hartree 'CHrho' calculated in INM with rc='hatree_rc'
  ! [2/2/2010]       is subtracted from Crho(0)at DMEorder >= 0. 
  !                  CHrho added to the public list, 'hatree_rc' added
  !                  to interaction parameters and the namelist.   
  !                  In the case DMEorder=-1 (standard Skyrme) 
  !                  both, 'CHrho' and 'hatree_rc', do not play.
  !                  New function HartreeV(u) defines Hatree energy as
  !                  E(Hartree)=(1/2)*Int[rho_0(r)*V(|r-r'|)*rho_0(r')]
  !                  HartreeV(u) is zero for u=<'hatree_rc'
  ! ver#8: (Markus)  Hartree DME terms dropped out.
  ! ver#7: (Markus)  Added switch to turn off the 3N terms.
  !        (Mario)   Added Abs to density and gradient dependent LDA       
  !                  Public :: DMEorder,DMElda,use_DME3N_terms   
  ! ver#6: (Mario)   Skyrme transformation added.
  ! ver#5: (Mario)   Print_Namelist=T/F added to the namelist 
  ! ver#4: (Markus)  Added natural units to the module. Used only for printing.
  ! ver#3: (Mario)   Uamplitudes(0:3,0:7) in normal order
  ! t for Uamplitudes(t,*)
  ! 0 -> 0,0
  ! 1 -> 1,1
  ! 2 -> 0,1
  ! 3 -> 1,0
  ! n for Uamplitudes(*,n)
  ! 0 -> U
  ! 1 -> dU/dRHO_0
  ! 2 -> dU/dRHO_1
  ! 3 -> d2U/(dRHO_0*dRHO_0)
  ! 4 -> d2U/(dRHO_1*dRHO_1) 
  ! 5 -> d2U/(dRHO_0*dRHO_1) 
  ! 6 -> dU/d(TAU_0)   
  ! 7 -> dU/d(Delta RHO_0)   
  ! TESTED MATTHEMETICA<=>BIRUC & SCOTT; MATTHEMETICA<=>Module UNEDF (energy amplitudes only)
  ! ver#2: (Mario) Pairing included
  !  - set_functional_parameters(fname,lpr)
  !  - pairing incorporated into CpV0(0:1),CpV1(0:1)
  !    as public variables also serving two public amplitudes
  !     Urhorhopr(0:1,0)=CpV0(0:1)+CpV1(0:1)*rho(0)
  !     Urhorhopr(0:1,1)=CpV1(0:1)
  !     so, they can be used with appropriate values by the DME solver
  !  -need improvement later, 
  !      currently HFBTHO uses CpV0(0:1), CpV0(0:1)  as before
  !      just substituting V0,V1 in pn-representation 
  !      CpV0*(1-CpV1/0.16*rho_0)and this defines
  !      the default values in the module CpV0=V0,CpV1=1/2)
  !  -NAMELIST and input/output modified. RESERVED NAMES ARE:
  !      -namelist forbiden:
  !          'UNRDF'  - best UNEDF
  !          'SKYRME' - best SKYRME
  !      -namelist inforced but not for C-parameters (use_INM=F) 
  !       or NM-parameters (use_INM=T) defined by the solver 
  !          'FITS'   
  !      -namelist inforced (one can overwrite all):
  !          'ANY OTHER NAME' 
  !       i.e., the solver defines C-/NM- only using 'FITS'       
  ! ver#1: (Mario) Complete rewrite consistent with HFBTHO
  !  -CB-LDA added
  !  -INM added 
  !  -HFBTHO BENCHMARK: LN, ZR(110) prolate solution with SLY4, 
  !   mixed pairing and tensor terms. Agreement with previouse 
  !   implemetation to the last significant digit in the cases:
  !      - Standard Skyrme  
  !      - LO+LDA
  !      - LO+CB-LDA
  !      - (NrNr=0,rDj=0), (rDr=0,jDr=0), 0.5(NrNr=-rDr,jDr=-rDj)  
  !   -use_j2terms removed, i.e., in the SKYRME case CJ=0 removes all
  !    tensor terms, while in DME tensor terms are always present     
  ! ver#0: (Marcus) Basic coding from scratch
  !   -DME(u) consistent with Mathematica numbers 
  !   -including small 'u' approximation  
  !-------------------------------------------------------------------
  !
  ! === PUBLIC VARIABLES ===
  !
  ! Use pointers to prevent conflicts with UNEDF public variabes
  ! Example: Use UNEDF, pr=>my_pr, ipr=>my_ipr, Crho=>my_Crho ...
  !
  Integer, Parameter, Public :: ipr=Kind(1)                                   ! to set the precision of the DFT solver
  Integer, Parameter, Public :: pr=Kind(1.000D0)                              ! to set the precision of the DFT solver
  Logical, Public :: use_charge_density, use_cm_cor           
  Real(pr), Public, Dimension(0:3,0:7) :: Urhorho,Urhotau,UrhoDrho,Unablarho  ! ph DME amplitudes   
  Real(pr), Public, Dimension(0:3,0:7) :: UJnablarho,UrhonablaJ,UJJ
  Real(pr), Public, Dimension(0:3,0:7) :: Urhorhopr                           ! pp amplitudes   
  Real(pr), Public, Dimension(0:1) :: UEnonstdr,UFnonstdr,URnonstdr           ! Other amplitudes  
  Real(pr), Public :: hbzero,sigma,e2charg,CExPar                             ! hbr^2/2m, DD sigma, e^2 charge, coul. exchange parameter
  Real(pr), Public, Dimension(0:1) :: Crho,Cdrho,Ctau,CrDr,CrdJ,CJ,CpV0,CpV1  ! basic coupling constants
  Real(pr), Public :: E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM,P_NM,KA_NM
  Real(pr), Public :: CHrho                                                   ! Crho(0) from the Hartree term in NM
  Logical, Public :: use_DME3N_terms,use_j2terms
  Integer(ipr), Public :: DMEorder,DMElda,use_TMR_pairing
  Real(pr), Public :: mpi,gA,fpi,c1,c3,c4,cd,ce,LambdaX
  Logical, Public  :: use_INM,use_Namelist,Print_Namelist  
  !  
  ! === PRIVATE VARIABLES ===
  !
  Real(pr), Private, Dimension(0:1) :: nuCrho,nuCdrho,nuCtau,nuCrDr           ! basic coupling constants in natural units
  Real(pr), Private, Dimension(0:1) :: nuCrdJ,nuCJ,nuCpV0,nuCpV1              !
  Real(pr), Private :: t0,t1,t2,t3,x0,x1,x2,x3,b4,b4p,te,to
  Real(pr), Private :: nuLambda,nufpi                                         ! parameters associated to natural units
  Real(pr), Private, Dimension(0:1) :: Cnrho,CJdr                             ! hidden and always zero
  Integer(ipr), Private :: i_cut                                              ! dmeorder: -1=Standard Skyrme, 0=LO, 1=NLO, 2=N2LO  
  Real(pr), Private :: Pi,eps                                                 ! dmelda: 0=Kf-LDA, 1=CB-LDA
  Real(pr), Private :: kfconst,CK                                             ! (3Pi^2/2)^(1/3)
  Real(pr), Parameter, Private :: mevfm=197.30_pr; 
  Real(pr), Private :: rho(0:1),tau(0:1),nrho2(0:1),lrho(0:1)
  Real(pr), Private :: mpi2,fpi2,fpi4,gA2,gA4,gA6,CHartree
  Real(pr), Private :: arhorho,brhorho,arhodrho,brhodrho,arhotau,brhotau,ajj,bjj,adrdr,bdrdr
  Real(pr), Private :: darhorho,dbrhorho,darhodrho,dbrhodrho,darhotau,dbrhotau,dajj,dbjj,dadrdr,dbdrdr
  Real(pr), Private :: ddarhodrho,ddbrhodrho,ddarhotau,ddbrhotau,ddarhorho,ddbrhorho
  Real(pr), Private :: hrho0rho0,hrho1rho1,hdr0dr0,hdr1dr1,hrho0Drho0,hrho1Drho0, &
       hrho1Drho1,hrho0tau0,hrho1tau0,hrho1tau1,hJ0dr0,hrho0DJ0,hJ1dr1,hrho1DJ1, &
       hJ0dr1,hrho1DJ0,hJ1dr0,hJ0J0,hJ0J1,hJ1J1     
  Real(pr), Private :: dhrho0rho0,dhrho1rho1,dhdr0dr0,dhdr1dr1,dhrho0Drho0, &
       dhrho1Drho0,dhrho1Drho1,dhrho0tau0,dhrho1tau0,dhrho1tau1,dhJ0dr0,dhrho0DJ0, &
       dhJ1dr1,dhrho1DJ1,dhJ0dr1,dhrho1DJ0,dhJ1dr0,dhJ0J0,dhJ0J1,dhJ1J1     
  Real(pr), Private :: ddhrho0rho0,ddhrho1rho1,ddhrho0Drho0,ddhrho1Drho0, &
       ddhrho1Drho1,ddhrho0tau0,ddhrho1tau0,ddhrho1tau1
  Real(pr), Private, Dimension(3,3,33) :: ctr0r0,ctr1r1,ctdr0dr0,ctdr1dr1, & ! coefficients for 3N part   
       ctr0Dr0,ctr1Dr0,ctr1Dr1,ctr0t0,ctr1t0,ctr1t1,ctJ0dr0,ctr0dJ0,ctJ1dr1, &
       ctr1dJ1,ctJ0dr1,ctr1dJ0,ctJ1dr0,ctJ0J0,ctJ0J1,ctJ1J1
  Real(pr), Private :: u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12
  Real(pr), Private :: ual,lual,atu,asqu,asqu4 
  Real(pr), Private :: ac2,ac3,acoord
  Parameter(acoord=0.50_pr,ac2=4.0_pr*(acoord**2-acoord+0.50_pr),ac3=2.0_pr*(acoord**2-acoord+0.50_pr))
  Character (30) :: FunctionalName 
  ! 
  Real(pr), Private :: A1_1,A1_2,A1_3,A1_4,A1_5,b1_1,b1_2,b1_3,b1_4,b1_5
  Real(pr), Private :: A3_1,A3_2,A3_3,A3_4,A3_5,b3_1,b3_2,b3_3,b3_4,b3_5    
  Real(pr), Private :: h0mpi6,h0mpi6c1,h0mpi6c3,h0mpi6NM,h0mpi6c1NM,h0mpi6c3NM
  ! 
  Namelist /UNEDF_NAMELIST/ FunctionalName,DMEorder,DMElda,use_INM,hbzero,use_TMR_pairing, &
       Crho,Cdrho,Ctau,CrDr,CrdJ,CJ,sigma,CpV0,CpV1,e2charg, &
       E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM, &
       mpi,gA,fpi,c1,c3,c4,cd,ce,LambdaX, &
       use_cm_cor,use_charge_density,use_DME3N_terms,use_j2terms,CExPar, &
       Print_Namelist     
Contains
  !
  !=======================================================================================================
  Subroutine calculate_U_parameters(rho0_in,rho1_in,tau0_in,tau1_in,laprho0,laprho1,nablarho0s,nablarho1s)
    !-----------------------------------------------------------------------------------------------------
    Implicit None
    Real(pr), Intent(in) :: rho0_in,rho1_in,tau0_in,tau1_in
    Real(pr), Intent(in), Optional :: &
         nablarho0s,nablarho1s,laprho0,laprho1
    Integer(ipr) :: t,i,j,k,l
    Real(pr) :: u,du,ddu,dtu,dlu
    Real(pr) :: ph,aux,daux,ddaux
    Real(pr) :: y,dy,ddy,marc,dmarc,ddmarc,mlog,dmlog,ddmlog
    Real(pr) :: ucut,ucut3n
    !
    ucut=0.1_pr; ucut3n=0.6_pr
    !    
    rho(0)=rho0_in; rho(1)=rho1_in; 
    tau(0)=tau0_in; tau(1)=tau1_in; 
    !
    lrho=0.0_pr; nrho2=0.0_pr; 
    If (Present(laprho0)) lrho(0)=laprho0 
    If (Present(laprho1)) lrho(1)=laprho1
    If (Present(nablarho0s)) nrho2(0)=nablarho0s 
    If (Present(nablarho1s)) nrho2(1)=nablarho1s
    !
    arhorho=0.0_pr; darhorho=0.0_pr; ddarhorho=0.0_pr
    brhorho=0.0_pr; dbrhorho=0.0_pr; ddbrhorho=0.0_pr
    arhodrho=0.0_pr; darhodrho=0.0_pr; ddarhodrho=0.0_pr
    brhodrho=0.0_pr; dbrhodrho=0.0_pr; ddbrhodrho=0.0_pr
    arhotau=0.0_pr; darhotau=0.0_pr; ddarhotau=0.0_pr
    brhotau=0.0_pr; dbrhotau=0.0_pr; ddbrhotau=0.0_pr  
    adrdr=0.0_pr; dadrdr=0.0_pr
    bdrdr=0.0_pr; dbdrdr=0.0_pr
    ajj=0.0_pr; dajj=0.0_pr
    bjj=0.0_pr; dbjj=0.0_pr    
    !
    hrho0rho0=0.0_pr; hrho1rho1=0.0_pr; hdr0dr0=0.0_pr; hdr1dr1=0.0_pr
    hrho0Drho0=0.0_pr; hrho1Drho0=0.0_pr; hrho1Drho1=0.0_pr 
    hrho0tau0=0.0_pr; hrho1tau0=0.0_pr; hrho1tau1=0.0_pr 
    hJ0dr0=0.0_pr; hrho0DJ0=0.0_pr; hJ1dr1=0.0_pr; hrho1DJ1=0.0_pr
    hJ0dr1=0.0_pr; hrho1DJ0=0.0_pr; hJ1dr0=0.0_pr 
    hJ0J0=0.0_pr; hJ0J1=0.0_pr; hJ1J1=0.0_pr
    dhrho0rho0=0.0_pr; dhrho1rho1=0.0_pr; dhdr0dr0=0.0_pr; dhdr1dr1=0.0_pr
    dhrho0Drho0=0.0_pr; dhrho1Drho0=0.0_pr; dhrho1Drho1=0.0_pr 
    dhrho0tau0=0.0_pr; dhrho1tau0=0.0_pr; dhrho1tau1=0.0_pr 
    dhJ0dr0=0.0_pr; dhrho0DJ0=0.0_pr; dhJ1dr1=0.0_pr; dhrho1DJ1=0.0_pr
    dhJ0dr1=0.0_pr; dhrho1DJ0=0.0_pr; dhJ1dr0=0.0_pr 
    dhJ0J0=0.0_pr; dhJ0J1=0.0_pr; dhJ1J1=0.0_pr
    ddhrho0rho0=0.0_pr; ddhrho1rho1=0.0_pr 
    ddhrho0Drho0=0.0_pr; ddhrho1Drho0=0.0_pr; ddhrho1Drho1=0.0_pr
    ddhrho0tau0=0.0_pr; ddhrho1tau0=0.0_pr; ddhrho1tau1=0.0_pr
    !
#ifndef hide_dme    
    ! U and partial derivatives with respect of rho_0 in Thomas Fermi approximation
    !
    u=0.0_pr; du=0.0_pr; ddu=0.0_pr; dtu=0.0_pr; dlu=0.0_pr
    If (dmeorder.Ge.0) Then
       If (dmelda.Eq.0) Then
          ! density dependent LDA
          u=(kfconst/mpi)*rho(0)**(1.0_pr/3.0_pr)
          du=(1.0_pr/3.0_pr)*u/(rho(0)+eps)                 ! u'(RHO_0)   
          ddu=-(2.0_pr/9.0_pr)*u/(rho(0)**2+eps)            ! u''(RHO_0)  
          dtu=0.0_pr                                        ! u'(TAU_0)   
          dlu=0.0_pr                                        ! u'(DeltaRHO_0) 
       Else 
          ! density and gradient dependent LDA       
          u=Sqrt(Abs((5.0_pr/3.0_pr)*(tau(0)-0.250_pr*lrho(0))/(rho(0)+eps)))/mpi 
          du=-0.50_pr*u/(rho(0)+eps)                        ! u'(RHO_0)   
          ddu=0.750_pr*u/(rho(0)**2+eps)                    ! u''(RHO_0)
          dtu=0.50_pr*u/(Abs(tau(0)-0.250_pr*lrho(0))+eps)  ! u'(TAU_0)   
          dlu=-0.250_pr*dtu                                 ! u'(DeltaRHO_0) 
       Endif
    Endif
    !
    ! Partial optimiztion
    u2=u*u; u3=u2*u; u4=u3*u; u5=u4*u; u6=u5*u; u7=u6*u; u8=u7*u; u9=u8*u; u10=u9*u; u11=u10*u; u12=u11*u;
    ual=1.0_pr+4.0_pr*u2; lual=Log(ual); atu=Atan(2.0_pr*u); asqu=Sqrt(1.0_pr+u2); asqu4=Sqrt(4.0_pr+u2)
    !
    ! A and B functions and their partial derivatives with respect of u
    !
    !  LO, 2N terms
    !
    If (dmeorder.Ge.0) Then
       If (u.Gt.ucut) Then          
          arhorho=(gA2*(4.0_pr*u2*(21.0_pr-498.0_pr*u2-64.0_pr*u4+16.0_pr*u6)+48.0_pr*u3 &
               *(35.0_pr+4.0_pr*u2)*atu+3.0_pr*(-7.0_pr+16.0_pr*u2*(-8.0_pr+9.0_pr*u2))*lual)) &
               /(1024.0_pr*fpi2*u8)
          darhorho=(gA2*(4.0_pr*u2*(-21.0_pr+279.0_pr*u2+16.0_pr*u4)-6.0_pr*u3*(175.0_pr &
               +12.0_pr*u2)*atu+3.0_pr*(7.0_pr+96.0_pr*u2-72.0_pr*u4)*lual))/(128.0_pr*fpi2*u9)
          ddarhorho=(-3.0_pr*gA2*(4.0_pr*u2*(-63.0_pr+504.0_pr*u2+16.0_pr*u4)-12.0_pr*u3*(175.0_pr &
               +8.0_pr*u2)*atu+(63.0_pr+672.0_pr*u2-360.0_pr*u4)*lual))/(128.0_pr*fpi2*u10)
          brhorho=2.0_pr*arhorho   
          dbrhorho=2.0_pr*darhorho   
          ddbrhorho=2.0_pr*ddarhorho
          !
          arhodrho=(35.0_pr*gA2*(-4.0_pr*u2*(-3.0_pr+72.0_pr*u2+4.0_pr*u4-60.0_pr*u*atu) &
               +3.0_pr*(-1.0_pr-18.0_pr*u2+24.0_pr*u4)*lual))/(12288.0_pr*fpi2*mpi2*u10)
          darhodrho=(35.0_pr*gA2*(4.0_pr*u2*(-15.0_pr+234.0_pr*u2+8.0_pr*u4 &
               -210.0_pr*u*atu)+3.0_pr*(5.0_pr+72.0_pr*u2-72.0_pr*u4)*lual)) &
               /(6144.0_pr*fpi2*mpi2*u11)
          ddarhodrho=(35.0_pr*gA2*(4.0_pr*u2*(165.0_pr-1746.0_pr*u2-40.0_pr*u4+1680.0_pr*u &
               *atu)+3.0_pr*(-55.0_pr+72.0_pr*u2*(-9.0_pr+7.0_pr*u2))*lual))/(6144.0_pr*fpi2*mpi2*u12)
          !
          brhodrho=2.0_pr*arhodrho
          dbrhodrho=2.0_pr*darhodrho
          ddbrhodrho=2.0_pr*ddarhodrho
          ! 
          arhotau=(35.0_pr*gA2*(4.0_pr*u2*(-3.0_pr+72.0_pr*u2+4.0_pr*u4)-240.0_pr*u3 &
               *atu+(3.0_pr+54.0_pr*u2-72.0_pr*u4)*lual))/(3072.0_pr*fpi2*mpi2*u10)
          darhotau=(35.0_pr*gA2*(4.0_pr*u2*(15.0_pr-234.0_pr*u2-8.0_pr*u4+210.0_pr*u &
               *atu)+3.0_pr*(-5.0_pr+72.0_pr*u2*(-1.0_pr+u2))*lual))/(1536.0_pr*fpi2*mpi2*u11)
          ddarhotau=(35.0_pr*gA2*(4.0_pr*u2*(-165.0_pr+1746.0_pr*u2+40.0_pr*u4-1680.0_pr*u*atu) &
               +3.0_pr*(55.0_pr+72.0_pr*u2*(9.0_pr-7.0_pr*u2))*lual))/(1536.0_pr*fpi2*mpi2*u12)
          brhotau=2.0_pr*arhotau
          dbrhotau=2.0_pr*darhotau
          ddbrhotau=2.0_pr*ddarhotau
          !
          ajj=(3.0_pr*gA2*(-4.0_pr*u2+8.0_pr*u4+lual))/(512.0_pr*fpi2*mpi2*u6)
          dajj=(-3.0_pr*gA2*(4.0_pr*u2*(-3.0_pr-6.0_pr*u2+8.0_pr*u4)+3.0_pr*(ual)*lual)) &
               /(256.0_pr*fpi2*mpi2*u7*(ual))
          bjj=2.0_pr*ajj
          dbjj=2.0_pr*dajj            
       Else
          Arhorho=(gA2*u4*(121.0_pr-448.0_pr*u2))/(1155.0_pr*fpi2)
          dArhorho=(4.0_pr*gA2*u3*(121.0_pr-672.0_pr*u2))/(1155.0_pr*fpi2)
          ddArhorho=(4.0_pr*gA2*u2*(1573.0_pr+560.0_pr*u2*(-26.0_pr+153.0_pr*u2)))/(5005.0_pr*fpi2)
          Brhorho=2.0_pr*Arhorho
          dBrhorho=2.0_pr*dArhorho
          ddBrhorho=2.0_pr*ddArhorho
          !     
          Arhotau=(gA2*(1287.0_pr-4004.0_pr*u2+11232.0_pr*u4-31680.0_pr*u6))/(10296.0_pr*fpi2*mpi2)
          dArhotau=(gA2*u*(-1001.0_pr+432.0_pr*u2*(13.0_pr-55.0_pr*u2)))/(1287.0_pr*fpi2*mpi2)
          ddArhotau=(gA2*(-1001.0_pr+16.0_pr*u2*(1053.0_pr-7425.0_pr*u2+40040.0_pr*u4)))/(1287.0_pr*fpi2*mpi2)
          Brhotau=2.0_pr*Arhotau
          dBrhotau=2.0_pr*dArhotau
          ddBrhotau=2.0_pr*ddArhotau
          !
          ArhoDrho=-Arhotau/4.0_pr
          dArhoDrho=-dArhotau/4.0_pr
          ddArhoDrho=-ddArhotau/4.0_pr
          BrhoDrho=2.0_pr*ArhoDrho
          dBrhoDrho=2.0_pr*dArhoDrho
          ddBrhoDrho=2.0_pr*ddArhoDrho
          !
          AJJ=(gA2*(5.0_pr-15.0_pr*u2+48.0_pr*u4-160.0_pr*u6))/(40.0_pr*fpi2*mpi2)
          dAJJ=(-3.0_pr*gA2*u*(5.0_pr-32.0_pr*u2+160.0_pr*u4))/(20.0_pr*fpi2*mpi2)
          Bjj=2.0_pr*Ajj; 
          dBjj=2.0_pr*dAjj   
       End If
    End If
    !
    ! NLO, 2N terms
    !
    If (dmeorder.Ge.1) Then
       If (u.Gt.ucut) Then
          arhorho=arhorho+(mpi2*(4.0_pr*u2*(11025.0_pr*(-1.0_pr + 10.0_pr*gA2 + 127.0_pr*gA4) &
               - 525.0_pr*(-839.0_pr + 3014.0_pr*gA2 + 50489.0_pr*gA4)*u2 - 4200.0_pr*(-17.0_pr &
               - 6.0_pr*gA2 + 495.0_pr*gA4)*u4 - 140.0_pr*(-133.0_pr - 718.0_pr*gA2 &
               + 107.0_pr*gA4)*u6 + 1536.0_pr*(1.0_pr + 10.0_pr*gA2 + 13.0_pr*gA4)*u8) &
               - 525.0_pr*Log(1.0_pr + 2.0_pr*u*(u + asqu))*(4.0_pr*u*asqu*(21.0_pr*(-1.0_pr &
               + 10.0_pr*gA2 + 127.0_pr*gA4) - 2.0_pr*(-167.0_pr + 518.0_pr*gA2 + 9305.0_pr*gA4)*u2 &
               + 8.0_pr*(7.0_pr + 10.0_pr*gA2 - 153.0_pr*gA4)*u4 + 16.0_pr*(1.0_pr + 6.0_pr*gA2 &
               + gA4)*u6) + 3.0_pr*(7.0_pr - 70.0_pr*gA2 - 889.0_pr*gA4 - 64.0_pr*(-1.0_pr + 6.0_pr*gA2 &
               + 83.0_pr*gA4)*u2 + 48.0_pr*(-1.0_pr + 2.0_pr*gA2 + 47.0_pr*gA4)*u4)*Log(1.0_pr &
               + 2.0_pr*u*(u + asqu)))))/(1.72032e7_pr*fpi4*Pi**2*u8)

          darhorho=darhorho+(mpi2*(4.0_pr*u2*(-3675.0_pr*(-1.0_pr+10.0_pr*gA2+127.0_pr*gA4) &
               +175.0_pr*(-503.0_pr+1766.0_pr*gA2+29897.0_pr*gA4)*u2+2800.0_pr*(-3.0_pr-2.0_pr*gA2 &
               +81.0_pr*gA4)*u4-700.0_pr*(1.0_pr+6.0_pr*gA2+gA4)*u6+128.0_pr*(1.0_pr+10.0_pr*gA2 &
               +13.0_pr*gA4)*u8)*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu))) &
               -700.0_pr*u*(1.0_pr+u2)*(21.0_pr-2*u2*(95.0_pr+8.0_pr*u2)+gA2*(-210.0_pr+556.0_pr*u2 &
               -32.0_pr*u4)+gA4*(-2667.0_pr+10306.0_pr*u2+304.0_pr*u4))*(1.0_pr+2.0_pr*u*(u+Sqrt(1.0_pr &
               +u2)))*Log(1.0_pr+2.0_pr*u*(u+asqu))+525.0_pr*(7.0_pr-70.0_pr*gA2 &
               -889.0_pr*gA4-48.0_pr*(-1.0_pr+6.0_pr*gA2+83.0_pr*gA4)*u2+24.0_pr*(-1.0_pr+2.0_pr*gA2 &
               +47.0_pr*gA4)*u4)*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu)))*Log(1.0_pr &
               +2.0_pr*u*(u+asqu))**2))/(716800.0_pr*fpi4*Pi**2*u9*asqu*(1.0_pr &
               +2.0_pr*u*(u+asqu)))

          ddarhorho=ddarhorho+(mpi2*(4.0_pr*u2*(33075.0_pr*(-1.0_pr+10.0_pr*gA2+127.0_pr*gA4) &
               -525.0_pr*(-965.0_pr+3314.0_pr*gA2+56699.0_pr*gA4)*u2-2800.0_pr*(-11.0_pr-10.0_pr*gA2 &
               +281.0_pr*gA4)*u4+700.0_pr*(1.0_pr+6.0_pr*gA2+gA4)*u6+128.0_pr*(1.0_pr+10.0_pr*gA2 &
               +13.0_pr*gA4)*u8)*(asqu+4.0_pr*(u+u3)*(1.0_pr+2.0_pr*u*(u &
               +asqu)))+525.0_pr*Log(1.0_pr+2.0_pr*u*(u+asqu))*(4.0_pr*u*(1.0_pr &
               +u2)*(63.0_pr-346.0_pr*u2-16.0_pr*u4+gA2*(-630.0_pr+964.0_pr*u2-32.0_pr*u4) &
               +gA4*(-8001.0_pr+18406.0_pr*u2+304.0_pr*u4))*(1.0_pr+4.0_pr*u*(asqu &
               +2.0_pr*u*(1.0_pr+u*(u+asqu))))-3.0_pr*(-21.0_pr*(-1.0_pr+10.0_pr*gA2 &
               +127.0_pr*gA4)-112.0_pr*(-1.0_pr+6.0_pr*gA2+83.0_pr*gA4)*u2+40.0_pr*(-1.0_pr &
               +2.0_pr*gA2+47.0_pr*gA4)*u4)*(asqu+4.0_pr*(u+u3)*(1.0_pr+2.0_pr*u*(u &
               +asqu)))*Log(1.0_pr+2.0_pr*u*(u+asqu)))))/(716800.0_pr*fpi4*Pi**2*u10*asqu &
               *(1.0_pr+2.0_pr*u*(u+asqu))**2)

          brhorho=brhorho+(mpi2*(4.0_pr*u2*(11025.0_pr*(-1.0_pr + 10.0_pr*gA2 + 43.0_pr*gA4) &
               - 525.0_pr*(-839.0_pr + 3014.0_pr*gA2 + 12269.0_pr*gA4)*u2 - 4200.0_pr*(-17.0_pr &
               - 6.0_pr*gA2 + 123.0_pr*gA4)*u4 - 140.0_pr*(-133.0_pr - 718.0_pr*gA2 &
               + 1223.0_pr*gA4)*u6 - 1536.0_pr*(-1.0_pr - 10.0_pr*gA2 + 23.0_pr*gA4)*u8) &
               + 525.0_pr*Log(1.0_pr + 2.0_pr*u*(u + asqu))*(4.0_pr*u*asqu*(21.0_pr &
               - 2.0_pr*u2*(167.0_pr + 28.0_pr*u2 + 8.0_pr*u4) - 2.0_pr*gA2*(105.0_pr &
               - 518.0_pr*u2 + 40.0_pr*u4 + 48.0_pr*u6) + gA4*(-903.0_pr + 4378.0_pr*u2 &
               + 360.0_pr*u4 + 176.0_pr*u6)) - 3.0_pr*(7.0_pr - 70.0_pr*gA2 - 301.0_pr*gA4 &
               - 64.0_pr*(-1.0_pr + 6.0_pr*gA2 + 23.0_pr*gA4)*u2 + 48.0_pr*(-1.0_pr + 2.0_pr*gA2 &
               + 11.0_pr*gA4)*u4)*Log(1.0_pr + 2.0_pr*u*(u + asqu)))))/(8.6016e6_pr*fpi4*Pi**2*u8)

          dbrhorho=dbrhorho+(mpi2*(-4.0_pr*u2*(3675.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4) &
               -175.0_pr*(-503.0_pr+1766.0_pr*gA2+7229.0_pr*gA4)*u2-2800.0_pr*(-3.0_pr-2.0_pr*gA2 &
               +21.0_pr*gA4)*u4-700.0_pr*(-1.0_pr-6.0_pr*gA2+11.0_pr*gA4)*u6+128.0_pr*(-1.0_pr &
               -10.0_pr*gA2+23.0_pr*gA4)*u8)*(asqu+2.0_pr*u*(1.0_pr+u*(u+Sqrt(1.0_pr &
               +u2))))-700.0_pr*u*(1.0_pr+u2)*(21.0_pr-2.0_pr*u2*(95.0_pr+8.0_pr*u2)+gA2*(-210.0_pr &
               +556.0_pr*u2-32.0_pr*u4)+gA4*(-903.0_pr+2410.0_pr*u2+112.0_pr*u4))*(1.0_pr+2.0_pr*u*(u &
               +asqu))*Log(1.0_pr+2.0_pr*u*(u+asqu))+525.0_pr*(7.0_pr-70.0_pr*gA2 &
               -301.0_pr*gA4-48.0_pr*(-1.0_pr+6.0_pr*gA2+23.0_pr*gA4)*u2+24.0_pr*(-1.0_pr+2.0_pr*gA2 &
               +11.0_pr*gA4)*u4)*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu)))*Log(1.0_pr &
               +2.0_pr*u*(u+asqu))**2))/(358400.0_pr*fpi4*Pi**2*u9*asqu*(1.0_pr &
               +2.0_pr*u*(u+asqu)))

          ddbrhorho=ddbrhorho+(mpi2*(-4.0_pr*u2*(-33075.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4) &
               +525.0_pr*(-965.0_pr+3314.0_pr*gA2+13655.0_pr*gA4)*u2+2800.0_pr*(-11.0_pr-10.0_pr*gA2 &
               +77.0_pr*gA4)*u4+700.0_pr*(-1.0_pr-6.0_pr*gA2+11.0_pr*gA4)*u6+128.0_pr*(-1.0_pr &
               -10.0_pr*gA2+23.0_pr*gA4)*u8)*(asqu+4.0_pr*(u+u3)*(1.0_pr+2.0_pr*u*(u &
               +asqu)))+525.0_pr*Log(1.0_pr+2.0_pr*u*(u+asqu))*(4.0_pr*u*(1.0_pr &
               +u2)*(-63.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4)+2.0_pr*(-173.0_pr+482.0_pr*gA2 &
               +2159.0_pr*gA4)*u2+16.0_pr*(-1.0_pr-2.0_pr*gA2+7.0_pr*gA4)*u4)*(1.0_pr &
               +4.0_pr*u*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu))))-3.0_pr*( &
               -21.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4)-112.0_pr*(-1.0_pr+6.0_pr*gA2 &
               +23.0_pr*gA4)*u2+40.0_pr*(-1.0_pr+2.0_pr*gA2+11.0_pr*gA4)*u4)*(asqu &
               +4.0_pr*(u+u3)*(1.0_pr+2.0_pr*u*(u+asqu)))*Log(1.0_pr+2.0_pr*u*(u &
               +asqu)))))/(358400.0_pr*fpi4*Pi**2*u10*asqu*(1.0_pr &
               +2.0_pr*u*(u+asqu))**2)

          arhodrho=arhodrho -(-44100.0_pr*(-1.0_pr+10.0_pr*gA2+127.0_pr*gA4)*u2+4200.0_pr*(-421.0_pr &
               +1516.0_pr*gA2+25369.0_pr*gA4)*u4+420.0_pr*(-283.0_pr-30.0_pr*gA2+8641.0_pr*gA4)*u6 &
               +1680.0_pr*(-19.0_pr-106.0_pr*gA2+5.0_pr*gA4)*u8+464.0_pr*(1.0_pr+10.0_pr*gA2 &
               +13.0_pr*gA4)*u10-420.0_pr*u*asqu*(105.0_pr-1685.0_pr*u2-134.0_pr*u4 &
               -8.0_pr*u6+16.0_pr*u8+10.0_pr*gA2*(-105.0_pr+527.0_pr*u2-22.0_pr*u4-8.0_pr*u6 &
               +16.0_pr*u8)+gA4*(-13335.0_pr+94295.0_pr*u2+2738.0_pr*u4-104.0_pr*u6+208.0_pr*u8)) &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))+11025.0_pr*(1.0_pr-10.0_pr*gA2-127.0_pr*gA4 &
               -9.0_pr*(-1.0_pr+6.0_pr*gA2+83.0_pr*gA4)*u2+8.0_pr*(-1.0_pr+2.0_pr*gA2 &
               +47.0_pr*gA4)*u4)*Log(1.0_pr+2.0_pr*u*(u+asqu))**2-128.0_pr*u3*(33600.0_pr &
               *(-1.0_pr+4.0_pr*gA2+16.0_pr*gA4)*u+560.0_pr*(-29.0_pr-10.0_pr*gA2+212.0_pr*gA4)*u3 &
               +420.0_pr*(-3.0_pr-16.0_pr*gA2+27.0_pr*gA4)*u5+33.0_pr*(-1.0_pr-10.0_pr*gA2 &
               +23.0_pr*gA4)*u7+16800.0_pr*(-1.0_pr+4.0_pr*gA2+16.0_pr*gA4)*asqu4*Log(2.0_pr) &
               +1680.0_pr*(-4.0_pr-5.0_pr*gA2+22.0_pr*gA4)*u2*asqu4*Log(2.0_pr)+210.0_pr &
               *(-1.0_pr-10.0_pr*gA2+23.0_pr*gA4)*u4*asqu4*Log(2.0_pr)+105.0_pr*(1.0_pr &
               +10.0_pr*gA2-23.0_pr*gA4)*u6*asqu4*Log(2.0_pr)+105.0_pr*asqu4 &
               *(-((-10.0_pr+u2)*(4.0_pr+u2)**2)-10.0_pr*gA2*(64.0_pr-8.0_pr*u2-2.0_pr*u4+u6) &
               +gA4*(-2560.0_pr-352.0_pr*u2-46.0_pr*u4+23.0_pr*u6))*Log(2.0_pr+u*(u+Sqrt(4.0_pr &
               +u2)))))/(4.128768e7_pr*fpi4*Pi**2*u10) 

          darhodrho=darhodrho+(-420.0_pr*u*(1.0_pr+u2)*asqu4*(15.0_pr-158.0_pr*u2 &
               -8.0_pr*u4-2.0_pr*gA2*(75.0_pr-238.0_pr*u2+8.0_pr*u4)+gA4*(-1905.0_pr+8690.0_pr*u2 &
               +152.0_pr*u4))*Log(1.0_pr+2.0_pr*u*(u+asqu))+315.0_pr*asqu &
               *asqu4*(5.0_pr-50.0_pr*gA2-635.0_pr*gA4-36.0_pr*(-1.0_pr+6.0_pr*gA2 &
               +83.0_pr*gA4)*u2+24.0_pr*(-1.0_pr+2.0_pr*gA2+47.0_pr*gA4)*u4)*Log(1.0_pr+2.0_pr*u*(u &
               +asqu))**2+4.0_pr*u2*asqu*(1575.0_pr*asqu4+gA4 & 
               *(-200025.0_pr*asqu4+u*(u*(845775.0_pr*asqu4+4.0_pr*u*(u*(-47705.0_pr &
               *asqu4+3.0_pr*u*(u*asqu4*(-1253.0_pr+236.0_pr*u2) &
               -7840.0_pr*Log(2.0_pr)))-268800.0_pr*Log(2.0_pr)))-3440640.0_pr*Log(2.0_pr))) &
               +u*(64575.0_pr*u*asqu4+33740.0_pr*u3*asqu4+1092.0_pr*u5*asqu4 &
               -48.0_pr*u7*asqu4+13440.0_pr*(4.0_pr+u2)**2*Log(2.0_pr))-2.0_pr*gA2*(7875.0_pr &
               *asqu4+u*(138915.0_pr*u*asqu4-8540.0_pr*u3*asqu4 &
               -3276.0_pr*u5*asqu4+240.0_pr*u7*asqu4-13440.0_pr*(-8.0_pr+u2) &
               *(4.0_pr+u2)*Log(2.0_pr)))+13440.0_pr*u*(-(4.0_pr+u2)**2+gA2*(64.0_pr+8.0_pr*u2 &
               -2.0_pr*u4)+gA4*(256.0_pr+80.0_pr*u2+7.0_pr*u4))*Log(2.0_pr+u*(u+Sqrt(4.0_pr &
               +u2)))))/(589824.0_pr*fpi4*Pi**2*u11*asqu*asqu4)

          ddarhodrho=ddarhodrho+(420.0_pr*u*asqu*(4.0_pr+u2)**1.5_pr*(-165.0_pr*(-1.0_pr &
               +10.0_pr*gA2+127.0_pr*gA4)+2.0_pr*(-593.0_pr+1730.0_pr*gA2+32183.0_pr*gA4)*u2 &
               +40.0_pr*(-1.0_pr-2.0_pr*gA2+19.0_pr*gA4)*u4)*Log(1.0_pr+2.0_pr*u*(u+asqu)) &
               -315.0_pr*(4.0_pr+u2)**1.5_pr*(-55.0_pr*(-1.0_pr+10.0_pr*gA2+127.0_pr*gA4)-324.0_pr &
               *(-1.0_pr+6.0_pr*gA2+83.0_pr*gA4)*u2+168.0_pr*(-1.0_pr+2.0_pr*gA2+47.0_pr*gA4)*u4) &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))**2-4.0_pr*u2*(-((4.0_pr+u2)*(-17325.0_pr*asqu4 &
               +u*(u*(-526365.0_pr*asqu4+4.0_pr*u*(u*(-48475.0_pr*asqu4+3.0_pr*u &
               *(u*asqu4*(-273.0_pr+4.0_pr*u2)-5600.0_pr*Log(2.0_pr)))-174720.0_pr*Log(2.0_pr))) &
               -1720320.0_pr*Log(2.0_pr))))+gA4*(-8801100.0_pr*asqu4+21899745.0_pr*u2*Sqrt(4.0_pr &
               +u2)+1906205.0_pr*u4*asqu4-1290772.0_pr*u6*asqu4-33780.0_pr*u8 &
               *asqu4+2832.0_pr*u10*asqu4-13440.0_pr*u*(8192.0_pr+4224.0_pr*u2 &
               +672.0_pr*u4+35.0_pr*u6)*Log(2.0_pr))-2.0_pr*gA2*(4.0_pr+u2)*(86625.0_pr*asqu4 &
               +u*(3440640.0_pr*Log(2.0_pr)+u*(1137465.0_pr*asqu4+4.0_pr*u*(u*(-16975.0_pr*Sqrt(4.0_pr &
               +u2)+3.0_pr*u*(u*asqu4*(-819.0_pr+20.0_pr*u2)-5600.0_pr*Log(2.0_pr))) &
               +107520.0_pr*Log(2.0_pr)))))+13440.0_pr*u*(512.0_pr*(-1.0_pr+4.0_pr*gA2+16.0_pr*gA4) &
               +48.0_pr*(-7.0_pr+16.0_pr*gA2+88.0_pr*gA4)*u2+24.0_pr*(-3.0_pr+gA2+28.0_pr*gA4)*u4 &
               +5.0_pr*(-1.0_pr-2.0_pr*gA2+7.0_pr*gA4)*u6)*Log(2.0_pr+u*(u+asqu4)))) &
               /(589824.0_pr*fpi4*Pi**2*u12*(4.0_pr+u2)**1.5_pr)

          brhodrho=brhodrho -(-44100.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4)*u2+4200.0_pr*(-421.0_pr &
               +1516.0_pr*gA2+6169.0_pr*gA4)*u4+420.0_pr*(-283.0_pr-30.0_pr*gA2+2029.0_pr*gA4)*u6 &
               +1680.0_pr*(-19.0_pr-106.0_pr*gA2+185.0_pr*gA4)*u8-464.0_pr*(-1.0_pr-10.0_pr*gA2 &
               +23.0_pr*gA4)*u10+420.0_pr*u*asqu*(-105.0_pr+1685.0_pr*u2+134.0_pr*u4 &
               +8.0_pr*u6-16.0_pr*u8+10.0_pr*gA2*(105.0_pr-527.0_pr*u2+22.0_pr*u4+8.0_pr*u6 &
               -16.0_pr*u8)+gA4*(4515.0_pr-22235.0_pr*u2-842.0_pr*u4-184.0_pr*u6+368.0_pr*u8)) &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))+11025.0_pr*(1.0_pr+9.0_pr*u2-8.0_pr*u4 &
               +2.0_pr*gA2*(-5.0_pr-27.0_pr*u2+8.0_pr*u4)+gA4*(-43.0_pr-207.0_pr*u2+88.0_pr*u4)) &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))**2+64.0_pr*u3*(33600.0_pr*(-1.0_pr+4.0_pr*gA2 &
               +16.0_pr*gA4)*u+560.0_pr*(-29.0_pr-10.0_pr*gA2+212.0_pr*gA4)*u3+420.0_pr*(-3.0_pr &
               -16.0_pr*gA2+27.0_pr*gA4)*u5+33.0_pr*(-1.0_pr-10.0_pr*gA2+23.0_pr*gA4)*u7 &
               +16800.0_pr*(-1.0_pr+4.0_pr*gA2+16.0_pr*gA4)*asqu4*Log(2.0_pr)+1680.0_pr*(-4.0_pr &
               -5.0_pr*gA2+22.0_pr*gA4)*u2*asqu4*Log(2.0_pr)+210.0_pr*(-1.0_pr-10.0_pr*gA2 &
               +23.0_pr*gA4)*u4*asqu4*Log(2.0_pr)+105.0_pr*(1.0_pr+10.0_pr*gA2 &
               -23.0_pr*gA4)*u6*asqu4*Log(2.0_pr)+105.0_pr*asqu4*(-((-10.0_pr &
               +u2)*(4.0_pr+u2)**2)-10.0_pr*gA2*(64.0_pr-8*u2-2*u4+u6)+gA4*(-2560.0_pr &
               -352.0_pr*u2-46.0_pr*u4+23.0_pr*u6))*Log(2.0_pr+u*(u+asqu4)))) &
               /(2.064384e7_pr*fpi4*Pi**2*u10)

          dbrhodrho=dbrhodrho+(-420.0_pr*u*(1.0_pr+u2)*asqu4*(15.0_pr-158.0_pr*u2 &
               -8.0_pr*u4-2.0_pr*gA2*(75.0_pr-238.0_pr*u2+8.0_pr*u4)+gA4*(-645.0_pr+2042.0_pr*u2 &
               +56.0_pr*u4))*Log(1.0_pr+2.0_pr*u*(u+asqu))+315.0_pr*asqu &
               *asqu4*(5.0_pr-50.0_pr*gA2-215.0_pr*gA4-36.0_pr*(-1.0_pr+6.0_pr*gA2 & 
               +23.0_pr*gA4)*u2+24.0_pr*(-1.0_pr+2.0_pr*gA2+11.0_pr*gA4)*u4)*Log(1.0_pr+2.0_pr*u*(u &
               +asqu))**2-4.0_pr*u2*asqu*(-1575.0_pr*asqu4 &
               +gA4*(67725.0_pr*asqu4+u*(u*(-1482075.0_pr*asqu4+4.0_pr*u*(u &
               *(-35035.0_pr*asqu4+3.0_pr*u*(u*asqu4*(-847.0_pr+184.0_pr*u2) &
               -3920.0_pr*Log(2.0_pr)))-134400.0_pr*Log(2.0_pr)))-1720320.0_pr*Log(2.0_pr))) &
               +u*(96705.0_pr*u*asqu4+20020.0_pr*u3*asqu4+924.0_pr*u5*asqu4 &
               -96.0_pr*u7*asqu4+6720.0_pr*(4.0_pr+u2)**2*Log(2.0_pr))+2.0_pr*gA2*(7875.0_pr &
               *asqu4+u*(-183645.0_pr*u*asqu4+4900.0_pr*u3*asqu4 &
               +2772.0_pr*u5*asqu4-480.0_pr*u7*asqu4+6720.0_pr*(-8.0_pr+u2) &
               *(4.0_pr+u2)*Log(2.0_pr)))+6720.0_pr*u*(-(4.0_pr+u2)**2+gA2*(64.0_pr+8.0_pr*u2 &
               -2.0_pr*u4)+gA4*(256.0_pr+80.0_pr*u2+7.0_pr*u4))*Log(2.0_pr+u*(u+asqu4)))) &
               /(294912.0_pr*fpi4*Pi**2*u11*asqu*asqu4)

          ddbrhodrho=ddbrhodrho+(420.0_pr*u*asqu*(4.0_pr+u2)**1.5_pr*(-165.0_pr*(-1.0_pr &
               +10.0_pr*gA2+43.0_pr*gA4)+2.0_pr*(-593.0_pr+1730.0_pr*gA2+7571.0_pr*gA4)*u2 &
               +40.0_pr*(-1.0_pr-2.0_pr*gA2+7.0_pr*gA4)*u4)*Log(1.0_pr+2.0_pr*u*(u+asqu)) &
               -315.0_pr*(4.0_pr+u2)**1.5_pr*(-55.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4)-324.0_pr*(-1.0_pr &
               +6.0_pr*gA2+23.0_pr*gA4)*u2+168.0_pr*(-1.0_pr+2.0_pr*gA2+11.0_pr*gA4)*u4)*Log(1.0_pr &
               +2.0_pr*u*(u+asqu))**2+4.0_pr*u2*(-((4.0_pr+u2)*(17325.0_pr*asqu4 &
               +u*(u*(-763875.0_pr*asqu4+4.0_pr*u*(u*(-28805.0_pr*asqu4+3.0_pr*u*(u &
               *asqu4*(-231.0_pr+8*u2)-2800.0_pr*Log(2.0_pr)))-87360.0_pr*Log(2.0_pr))) &
               -860160.0_pr*Log(2.0_pr))))+ gA4*(2979900.0_pr*asqu4-45909045.0_pr*u2*Sqrt(4.0_pr &
               +u2)-14728385.0_pr*u4*asqu4-928508.0_pr*u6*asqu4-21660.0_pr*u8 &
               *asqu4+2208.0_pr*u10*asqu4-6720.0_pr*u*(8192.0_pr+4224.0_pr*u2 &
               +672.0_pr*u4+35.0_pr*u6)*Log(2.0_pr))-2.0_pr*gA2*(4.0_pr+u2)*(-86625.0_pr*asqu4 &
               +u*(1720320.0_pr*Log(2.0_pr)+u*(1443015.0_pr*asqu4+4.0_pr*u*(u*(-9905.0_pr*Sqrt(4.0_pr &
               +u2)+3.0_pr*u*(u*asqu4*(-693.0_pr+40.0_pr*u2)-2800.0_pr*Log(2.0_pr))) &
               +53760.0_pr*Log(2.0_pr)))))+6720.0_pr*u*(512.0_pr*(-1.0_pr+4.0_pr*gA2+16.0_pr*gA4) &
               +48.0_pr*(-7.0_pr+16.0_pr*gA2+88.0_pr*gA4)*u2+24.0_pr*(-3.0_pr+gA2+28.0_pr*gA4)*u4 &
               +5.0_pr*(-1.0_pr-2.0_pr*gA2+7.0_pr*gA4)*u6)*Log(2.0_pr+u*(u+asqu4)))) &
               /(294912.0_pr*fpi4*Pi**2*u12*(4.0_pr+u2)**1.5_pr) 

          adrdr=adrdr +0.0_pr
          dadrdr=dadrdr +0.0_pr
          bdrdr=bdrdr +0.0_pr
          dbdrdr=dbdrdr +0.0_pr

          arhotau=arhotau+(-44100.0_pr*(-1.0_pr+10.0_pr*gA2+127.0_pr*gA4)*u2+4200.0_pr*(-421.0_pr &
               +1516.0_pr*gA2+25369.0_pr*gA4)*u4+420.0_pr*(-283.0_pr-30.0_pr*gA2+8641.0_pr*gA4)*u6 &
               +1680.0_pr*(-19.0_pr-106.0_pr*gA2+5.0_pr*gA4)*u8+464.0_pr*(1.0_pr+10.0_pr*gA2 &
               +13.0_pr*gA4)*u10-420.0_pr*u*asqu*(105.0_pr-1685.0_pr*u2-134.0_pr*u4 &
               -8.0_pr*u6+16.0_pr*u8+10.0_pr*gA2*(-105.0_pr+527.0_pr*u2-22.0_pr*u4-8.0_pr*u6 &
               +16.0_pr*u8)+gA4*(-13335.0_pr+94295.0_pr*u2+2738.0_pr*u4-104.0_pr*u6+208.0_pr*u8)) &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))+11025.0_pr*(1.0_pr-10.0_pr*gA2-127.0_pr*gA4 &
               -9.0_pr*(-1.0_pr+6.0_pr*gA2+83.0_pr*gA4)*u2+8.0_pr*(-1.0_pr+2.0_pr*gA2 &
               +47.0_pr*gA4)*u4)*Log(1.0_pr+2.0_pr*u*(u+asqu))**2)/(1.032192e7_pr*fpi4*Pi**2*u10)

          darhotau=darhotau+(-4.0_pr*u2*(-525.0_pr*(-1.0_pr+10.0_pr*gA2+127.0_pr*gA4)+35.0_pr &
               *(-409.0_pr+1450.0_pr*gA2+24439.0_pr*gA4)*u2+140.0_pr*(-5.0_pr-2.0_pr*gA2 &
               +143.0_pr*gA4)*u4-84.0_pr*(1.0_pr+6.0_pr*gA2+gA4)*u6+16.0_pr*(1.0_pr+10.0_pr*gA2 &
               +13.0_pr*gA4)*u8)*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu)))+35.0_pr &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))*(4.0_pr*u*(1.0_pr+u2)*(15.0_pr-158.0_pr*u2 &
               -8.0_pr*u4-2.0_pr*gA2*(75.0_pr-238.0_pr*u2+8.0_pr*u4)+gA4*(-1905.0_pr+8690.0_pr*u2 &
               +152.0_pr*u4))*(1.0_pr+2.0_pr*u*(u+asqu))-3.0_pr*(5.0_pr-50.0_pr*gA2 &
               -635.0_pr*gA4-36.0_pr*(-1.0_pr+6.0_pr*gA2+83.0_pr*gA4)*u2+24.0_pr*(-1.0_pr+2.0_pr*gA2 &
               +47.0_pr*gA4)*u4)*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu)))*Log(1.0_pr &
               +2.0_pr*u*(u+asqu))))/(49152.0_pr*fpi4*Pi**2*u11*asqu*(1.0_pr &
               +2.0_pr*u*(u+asqu)))

          ddarhotau=ddarhotau+(4.0_pr*u2*(-5775.0_pr*(-1.0_pr+10.0_pr*gA2+127.0_pr*gA4)+35.0_pr &
               *(-3179.0_pr+11102.0_pr*gA2+188453.0_pr*gA4)*u2+140.0_pr*(-29.0_pr-18.0_pr*gA2 &
               +791.0_pr*gA4)*u4-252.0_pr*(1.0_pr+6.0_pr*gA2+gA4)*u6+16.0_pr*(1.0_pr+10.0_pr*gA2 &
               +13.0_pr*gA4)*u8)*(asqu+4.0_pr*(u+u3)*(1.0_pr+2.0_pr*u*(u+Sqrt(1.0_pr &
               +u2))))-35.0_pr*Log(1.0_pr+2.0_pr*u*(u+asqu))*(4.0_pr*u*(1.0_pr+u2)*(-165.0_pr &
               *(-1.0_pr+10.0_pr*gA2+127.0_pr*gA4)+2.0_pr*(-593.0_pr+1730.0_pr*gA2+32183.0_pr*gA4)*u2 &
               +40.0_pr*(-1.0_pr-2.0_pr*gA2+19.0_pr*gA4)*u4)*(1.0_pr+4.0_pr*u*(asqu &
               +2.0_pr*u*(1.0_pr+u*(u+asqu))))-3.0_pr*(-55.0_pr*(-1.0_pr+10.0_pr*gA2 &
               +127.0_pr*gA4)-324.0_pr*(-1.0_pr+6.0_pr*gA2+83.0_pr*gA4)*u2+168.0_pr*(-1.0_pr &
               +2.0_pr*gA2+47.0_pr*gA4)*u4)*(asqu+4.0_pr*(u+u3)*(1.0_pr+2.0_pr*u*(u &
               +asqu)))*Log(1.0_pr+2.0_pr*u*(u+asqu))))/(49152.0_pr*fpi4*Pi**2*u12 &
               *asqu*(1.0_pr+2.0_pr*u*(u+asqu))**2)

          brhotau=brhotau+(-44100.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4)*u2+4200.0_pr*(-421.0_pr &
               +1516.0_pr*gA2+6169.0_pr*gA4)*u4+420.0_pr*(-283.0_pr-30.0_pr*gA2+2029.0_pr*gA4)*u6 &
               +1680.0_pr*(-19.0_pr-106.0_pr*gA2+185.0_pr*gA4)*u8-464.0_pr*(-1.0_pr-10.0_pr*gA2 &
               +23.0_pr*gA4)*u10+105.0_pr*Log(1.0_pr+2.0_pr*u*(u+asqu))*(4.0_pr*u*Sqrt(1.0_pr &
               +u2)*(-105.0_pr+1685.0_pr*u2+134.0_pr*u4+8.0_pr*u6-16.0_pr*u8+10.0_pr*gA2*(105.0_pr &
               -527.0_pr*u2+22.0_pr*u4+8.0_pr*u6-16.0_pr*u8)+gA4*(4515.0_pr-22235.0_pr*u2 &
               -842.0_pr*u4-184.0_pr*u6+368.0_pr*u8))+105.0_pr*(1.0_pr+9.0_pr*u2-8.0_pr*u4 &
               +2.0_pr*gA2*(-5.0_pr-27.0_pr*u2+8.0_pr*u4)+gA4*(-43.0_pr-207.0_pr*u2+88.0_pr*u4)) &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))))/(5.16096e6*fpi4*Pi**2*u10)

          dbrhotau=dbrhotau+(4.0_pr*u2*(525.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4)-35.0_pr*(-409.0_pr &
               +1450.0_pr*gA2+5923.0_pr*gA4)*u2+140.0_pr*(5.0_pr+2.0_pr*gA2-35.0_pr*gA4)*u4 &
               +84.0_pr*(1.0_pr+6.0_pr*gA2-11.0_pr*gA4)*u6+16.0_pr*(-1.0_pr-10.0_pr*gA2 &
               +23.0_pr*gA4)*u8)*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu)))+35.0_pr &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))*(4.0_pr*u*(1.0_pr+u2)*(15.0_pr-158.0_pr*u2 &
               -8.0_pr*u4-2.0_pr*gA2*(75.0_pr-238.0_pr*u2+8.0_pr*u4)+gA4*(-645.0_pr+2042.0_pr*u2 &
               +56.0_pr*u4))*(1.0_pr+2.0_pr*u*(u+asqu))-3.0_pr*(5.0_pr-50.0_pr*gA2 &
               -215.0_pr*gA4-36.0_pr*(-1.0_pr+6.0_pr*gA2+23.0_pr*gA4)*u2+24.0_pr*(-1.0_pr+2.0_pr*gA2 &
               +11.0_pr*gA4)*u4)*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu)))*Log(1.0_pr &
               +2.0_pr*u*(u+asqu))))/(24576.0_pr*fpi4*Pi**2*u11*asqu*(1.0_pr &
               +2.0_pr*u*(u+asqu)))

          ddbrhotau=ddbrhotau+(-4.0_pr*u2*(5775.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4)-35.0_pr &
               *(-3179.0_pr+11102.0_pr*gA2+45545.0_pr*gA4)*u2+140.0_pr*(29.0_pr+18.0_pr*gA2 &
               -203.0_pr*gA4)*u4-252.0_pr*(-1.0_pr-6.0_pr*gA2+11.0_pr*gA4)*u6+16.0_pr*(-1.0_pr &
               -10.0_pr*gA2+23.0_pr*gA4)*u8)*(asqu+4.0_pr*(u+u3)*(1.0_pr+2.0_pr*u*(u &
               +asqu)))-35.0_pr*Log(1.0_pr+2.0_pr*u*(u+asqu))*(4.0_pr*u*(1.0_pr &
               +u2)*(-165.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4)+2.0_pr*(-593.0_pr+1730.0_pr*gA2 &
               +7571.0_pr*gA4)*u2+40.0_pr*(-1.0_pr-2.0_pr*gA2+7.0_pr*gA4)*u4)*(1.0_pr+4.0_pr*u &
               *(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu))))-3.0_pr*(-55.0_pr*(-1.0_pr &
               +10.0_pr*gA2+43.0_pr*gA4)-324.0_pr*(-1.0_pr+6.0_pr*gA2+23.0_pr*gA4)*u2 &
               +168.0_pr*(-1.0_pr+2.0_pr*gA2+11.0_pr*gA4)*u4)*(asqu+4.0_pr*(u+u3)*(1.0_pr &
               +2.0_pr*u*(u+asqu)))*Log(1.0_pr+2.0_pr*u*(u+asqu)))) &
               /(24576.0_pr*fpi4*Pi**2*u12*asqu*(1.0_pr+2.0_pr*u*(u+asqu))**2)

          ajj=ajj+(4.0_pr*u2*(-9.0_pr+24.0_pr*u2+u4+gA4*(18.0_pr-363.0_pr*u2-113.0_pr*u4) &
               +2.0_pr*gA2*(9.0_pr+66.0_pr*u2+5.0_pr*u4))-12.0_pr*u*asqu*(-3.0_pr+u2 &
               -2.0_pr*u4+2.0_pr*gA2*(3.0_pr+5.0_pr*u2-10.0_pr*u4)+gA4*(6.0_pr-59.0_pr*u2 &
               +10.0_pr*u4))*Log(1.0_pr+2.0_pr*u*(u+asqu))+9.0_pr*(-1.0_pr+2.0_pr*gA2 &
               +2.0_pr*gA4+3.0_pr*(-1.0_pr-2.0_pr*gA2+gA4)*u2)*Log(1.0_pr+2.0_pr*u*(u &
               +asqu))**2)/(73728.0_pr*fpi4*Pi**2*u6)

          dajj=dajj+(-4.0_pr*u2*(-9.0_pr+9.0_pr*u2-2.0_pr*u4+2.0_pr*gA2*(9.0_pr+27.0_pr*u2 &
               -10.0_pr*u4)+2.0_pr*gA4*(9.0_pr-90.0_pr*u2+5.0_pr*u4))*(asqu+2.0_pr*u &
               *(1.0_pr+u*(u+asqu)))-9.0_pr*Log(1.0_pr+2.0_pr*u*(u+asqu)) &
               *(4.0_pr*u*(1.0_pr+u2-2.0_pr*gA2*(1.0_pr+u2)+gA4*(-2.0_pr+8.0_pr*u2+6.0_pr*u4)) &
               *(1.0_pr+2.0_pr*u*(u+asqu))+(-1.0_pr-2.0_pr*u2+2.0_pr*gA2*(1.0_pr+gA2 &
               +(-2.0_pr+gA2)*u2))*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu))) &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))))/(12288.0_pr*fpi4*Pi**2*u7*asqu &
               *(1.0_pr+2.0_pr*u*(u+asqu)))

          bjj=bjj+(4.0_pr*u2*(-9.0_pr+24.0_pr*u2+u4+gA4*(99.0_pr-228.0_pr*u2-23.0_pr*u4) &
               +2.0_pr*gA2*(9.0_pr+66.0_pr*u2+5.0_pr*u4))-12.0_pr*u*asqu*(-3.0_pr+u2 &
               -2.0_pr*u4+2.0_pr*gA2*(3.0_pr+5.0_pr*u2-10.0_pr*u4)+gA4*(33.0_pr-23.0_pr*u2 &
               +46.0_pr*u4))*Log(1.0_pr+2.0_pr*u*(u+asqu))+9.0_pr*(-1.0_pr+2.0_pr*gA2 &
               +11.0_pr*gA4+3.0_pr*(-1.0_pr-2.0_pr*gA2+7.0_pr*gA4)*u2)*Log(1.0_pr+2.0_pr*u*(u &
               +asqu))**2)/(36864.0_pr*fpi4*Pi**2*u6)

          dbjj=dbjj+(-4.0_pr*u2*(-9.0_pr+9.0_pr*u2-2.0_pr*u4+2.0_pr*gA2*(9.0_pr+27.0_pr*u2 &
               -10.0_pr*u4)+gA4*(99.0_pr-99.0_pr*u2+46.0_pr*u4))*(asqu+2.0_pr*u*(1.0_pr &
               +u*(u+asqu)))+36.0_pr*u*(-1.0_pr-u2+2.0_pr*gA2*(1.0_pr+u2)+gA4*(11.0_pr &
               +7.0_pr*u2))*(1.0_pr+2.0_pr*u*(u+asqu))*Log(1.0_pr+2.0_pr*u*(u+Sqrt(1.0_pr &
               +u2)))-9.0_pr*(-1.0_pr+2.0_pr*gA2+11.0_pr*gA4+2.0_pr*(-1.0_pr-2.0_pr*gA2 &
               +7.0_pr*gA4)*u2)*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu)))*Log(1.0_pr &
               +2.0_pr*u*(u+asqu))**2)/(6144.0_pr*fpi4*Pi**2*u7*asqu*(1.0_pr &
               +2.0_pr*u*(u+asqu)))
       Else    !! u smaller than ucut
          arhorho=arhorho +((-1.0_pr - 4.0_pr*gA2 + 8.0_pr*gA**4)*mpi2)/(384.0_pr*fpi4*Pi**2) &
               + (11.0_pr*(1.0_pr + 14.0_pr*gA2 + 17.0_pr*gA**4)*mpi2*u4)/(201600.0_pr*fpi4*Pi**2) &
               - ((1.0_pr + 18.0_pr*gA2 + 13.0_pr*gA**4)*mpi2*u6)/(69300.0_pr*fpi4*Pi**2)
          darhorho=darhorho+(mpi2*u3*(121.0_pr*(1.0_pr+14.0_pr*gA2+17.0_pr*gA4)-48.0_pr*(1.0_pr &
               +18.0_pr*gA2+13.0_pr*gA4)*u2))/(554400.0_pr*fpi4*Pi**2)
          ddarhorho=ddarhorho+(mpi2*u2*(1573.0_pr*(1.0_pr+14.0_pr*gA2+17.0_pr*gA4)-1040.0_pr*(1.0_pr &
               +18.0_pr*gA2+13.0_pr*gA4)*u2+680.0_pr*(1.0_pr+22.0_pr*gA2+gA4)*u4)) &
               /(2.4024e6_pr*fpi4*Pi**2)

          brhorho=brhorho +((-1.0_pr - 4.0_pr*gA2 + 8.0_pr*gA**4)*mpi2)/(192.0_pr*fpi4*Pi**2) &
               - (11.0_pr*(-1.0_pr - 14.0_pr*gA2 + 43.0_pr*gA**4)*mpi2*u4)/(100800.0_pr*fpi4*Pi**2) &
               + ((-1.0_pr - 18.0_pr*gA2 + 71.0_pr*gA**4)*mpi2*u6)/(34650.0_pr*fpi4*Pi**2)
          dbrhorho=dbrhorho+(mpi2*u3*(121.0_pr*(1.0_pr+14.0_pr*gA2-43.0_pr*gA4)+48.0_pr*(-1.0_pr &
               -18.0_pr*gA2+71.0_pr*gA4)*u2))/(277200.0_pr*fpi4*Pi**2)
          ddbrhorho=ddbrhorho+(mpi2*u2*(-1573.0_pr*(-1.0_pr-14.0_pr*gA2+43.0_pr*gA4)+1040.0_pr &
               *(-1.0_pr-18.0_pr*gA2+71.0_pr*gA4)*u2-680.0_pr*(-1.0_pr-22.0_pr*gA2 &
               +107.0_pr*gA4)*u4))/(1.2012e6_pr*fpi4*Pi**2)

          arhodrho=arhodrho+(120120.0_pr*(-2.0_pr-17.0_pr*gA2+88.0_pr*gA4)+8008.0_pr*(1.0_pr &
               +14.0_pr*gA2+167.0_pr*gA4)*u2-39.0_pr*(109.0_pr+1962.0_pr*gA2+4357.0_pr*gA4)*u4 &
               +540.0_pr*(3.0_pr+66.0_pr*gA2+31.0_pr*gA4)*u6)/(5.5351296e8_pr*fpi4*Pi**2)
          darhodrho=darhodrho+(4004.0_pr*(1.0_pr+14.0_pr*gA2+167.0_pr*gA4)*u-39.0_pr*(109.0_pr &
               +1962.0_pr*gA2+4357.0_pr*gA4)*u3+810.0_pr*(3.0_pr+66.0_pr*gA2+31.0_pr*gA4)*u5) &
               /(1.3837824e8_pr*fpi4*Pi**2)
          ddarhodrho=ddarhodrho+(12012.0_pr*(1.0_pr+14.0_pr*gA2+167.0_pr*gA4)-351.0_pr*(109.0_pr &
               +1962.0_pr*gA2+4357.0_pr*gA4)*u2+12150.0_pr*(3.0_pr+66.0_pr*gA2+31.0_pr*gA4)*u4 &
               +35.0_pr*(-811.0_pr-21086.0_pr*gA2+12637.0_pr*gA4)*u6)/(4.1513472e8_pr*fpi4*Pi**2)

          brhodrho=brhodrho -(480480.0_pr*(-2.0_pr-17.0_pr*gA2+34.0_pr*gA4)+76076.0_pr*(-1.0_pr &
               -14.0_pr*gA2+43.0_pr*gA4)*u2-12597.0_pr*(-1.0_pr-18.0_pr*gA2+71.0_pr*gA4)*u4 &
               +3660.0_pr*(-1.0_pr-22.0_pr*gA2+107.0_pr*gA4)*u6)/(5.5351296e8_pr*fpi4*Pi**2) 
          dbrhodrho=dbrhodrho+(-38038.0_pr*(-1.0_pr-14.0_pr*gA2+43.0_pr*gA4)*u+12597.0_pr*(-1.0_pr &
               -18.0_pr*gA2+71.0_pr*gA4)*u3-5490.0_pr*(-1.0_pr-22.0_pr*gA2+107.0_pr*gA4)*u5) &
               /(1.3837824e8_pr*fpi4*Pi**2)
          ddbrhodrho=ddbrhodrho+(-114114.0_pr*(-1.0_pr-14.0_pr*gA2+43.0_pr*gA4)+113373.0_pr*(-1.0_pr &
               -18.0_pr*gA2+71.0_pr*gA4)*u2-82350.0_pr*(-1.0_pr-22.0_pr*gA2+107.0_pr*gA4)*u4 &
               +58975.0_pr*(-1.0_pr-26.0_pr*gA2+151.0_pr*gA4)*u6)/(4.1513472e8_pr*fpi4*Pi**2)

          adrdr=adrdr +0.0_pr
          dadrdr=dadrdr +0.0_pr
          bdrdr=bdrdr +0.0_pr
          dbdrdr=dbdrdr +0.0_pr

          arhotau=arhotau -(30030.0_pr*(2.0_pr+17.0_pr*gA2+20.0_pr*gA4)+7007.0_pr*(1.0_pr+14.0_pr*gA2 &
               +17.0_pr*gA4)*u2-1404.0_pr*(1.0_pr+18.0_pr*gA2+13.0_pr*gA4)*u4+440.0_pr*(1.0_pr &
               +22.0_pr*gA2+gA4)*u6)/(3.459456e7_pr*fpi4*Pi**2)
          darhotau=darhotau+(-7007.0_pr*(1.0_pr+14.0_pr*gA2+17.0_pr*gA4)*u+2808.0_pr*(1.0_pr &
               +18.0_pr*gA2+13.0_pr*gA4)*u3-1320.0_pr*(1.0_pr+22.0_pr*gA2+gA4)*u5) &
               /(1.729728e7_pr*fpi4*Pi**2)
          ddarhotau=ddarhotau -(21021.0_pr*(1.0_pr+14.0_pr*gA2+17.0_pr*gA4)-25272.0_pr*(1.0_pr &
               +18.0_pr*gA2+13.0_pr*gA4)*u2+19800.0_pr*(1.0_pr+22.0_pr*gA2+gA4)*u4+14560.0_pr &
               *(-1.0_pr-26.0_pr*gA2+19.0_pr*gA4)*u6)/(5.189184e7_pr*fpi4*Pi**2) 

          brhotau=brhotau+(30030.0_pr*(-2.0_pr-17.0_pr*gA2+34.0_pr*gA4)+7007.0_pr*(-1.0_pr &
               -14.0_pr*gA2+43.0_pr*gA4)*u2-1404.0_pr*(-1.0_pr-18.0_pr*gA2+71.0_pr*gA4)*u4 &
               +440.0_pr*(-1.0_pr-22.0_pr*gA2+107.0_pr*gA4)*u6)/(1.729728e7_pr*fpi4*Pi**2)
          dbrhotau=dbrhotau+(7007.0_pr*(-1.0_pr-14.0_pr*gA2+43.0_pr*gA4)*u-2808.0_pr*(-1.0_pr &
               -18.0_pr*gA2+71.0_pr*gA4)*u3+1320.0_pr*(-1.0_pr-22.0_pr*gA2+107.0_pr*gA4)*u5) &
               /(8.64864e6_pr*fpi4*Pi**2)
          ddbrhotau=ddbrhotau+(21021.0_pr*(-1.0_pr-14.0_pr*gA2+43.0_pr*gA4)-25272.0_pr*(-1.0_pr &
               -18.0_pr*gA2+71.0_pr*gA4)*u2+19800.0_pr*(-1.0_pr-22.0_pr*gA2+107.0_pr*gA4)*u4 &
               -14560.0_pr*(-1.0_pr-26.0_pr*gA2+151.0_pr*gA4)*u6)/(2.594592e7_pr*fpi4*Pi**2)

          ajj=ajj+(-1050.0_pr*(-2.0_pr-17.0_pr*gA2+7.0_pr*gA4)+315.0_pr*(1.0_pr+14.0_pr*gA2 &
               -28.0_pr*gA4)*u2+81.0_pr*(-1.0_pr-18.0_pr*gA2+57.0_pr*gA4)*u4+16.0_pr*(2.0_pr &
               +44.0_pr*gA2-187.0_pr*gA4)*u6)/(2.4192e6_pr*fpi4*Pi**2)
          dajj=dajj+(105.0_pr*(1.0_pr+14.0_pr*gA2-28.0_pr*gA4)*u+54.0_pr*(-1.0_pr-18.0_pr*gA2 &
               +57.0_pr*gA4)*u3+16.0_pr*(2.0_pr+44.0_pr*gA2-187.0_pr*gA4)*u5)/(403200.0_pr*fpi4*Pi**2)

          bjj=bjj+(-1050.0_pr*(-2.0_pr-17.0_pr*gA2+34.0_pr*gA4)-315.0_pr*(-1.0_pr-14.0_pr*gA2 &
               +43.0_pr*gA4)*u2+81.0_pr*(-1.0_pr-18.0_pr*gA2+71.0_pr*gA4)*u4+32.0_pr*(1.0_pr &
               +22.0_pr*gA2-107.0_pr*gA4)*u6)/(1.2096e6_pr*fpi4*Pi**2)
          dbjj=dbjj+(105.0_pr*(1.0_pr+14.0_pr*gA2-43.0_pr*gA4)*u+54.0_pr*(-1.0_pr-18.0_pr*gA2 &
               +71.0_pr*gA4)*u3+32.0_pr*(1.0_pr+22.0_pr*gA2-107.0_pr*gA4)*u5)/(201600.0_pr*fpi4*Pi**2)

       End If     !! if (u.gt.ucut...
    End If    !! if (dmeorder.ge.1....
    !
    ! N2LO, 2N terms
    !
    If (dmeorder.Ge.2) Then
       If (u.Gt.ucut) Then
          arhorho=arhorho +(gA2*mpi**3*(u2*(210.0_pr*(-129.0_pr*c3 + 56.0_pr*c4) + 6.0_pr*(38607.0_pr*c3 &
               - 27325.0_pr*c4)*u2 + 7.0_pr*(3165.0_pr*c3 - 2014.0_pr*c4)*u4 + 2714.0_pr*(3.0_pr*c3 &
               + c4)*u6 + 570.0_pr*(3.0_pr*c3 + c4)*u8 - 66.0_pr*c1*(1050.0_pr - 6366.0_pr*u2 + 299.0_pr*u4 &
               + 237.0_pr*u6)) + 3.0_pr*u3*(-242550.0_pr*c1 + 6237.0_pr*(2.0_pr*c1 - c3)*u4 &
               + 572.0_pr*(3.0_pr*c1 - 3.0_pr*c3 - c4)*u6 + 190.0_pr*(3.0_pr*c3 + c4)*u8 - 1155.0_pr*(3.0_pr*c3 &
               - 2.0_pr*c4)*(35.0_pr + 2.0_pr*u2))*Atan(u) - 6.0_pr*(-35.0_pr*(330.0_pr*c1 + 129.0_pr*c3 &
               - 56.0_pr*c4) - 352.0_pr*(162.0_pr*c1 + 69.0_pr*c3 - 40.0_pr*c4)*u2 + 891.0_pr*(14.0_pr*c1 &
               + 11.0_pr*c3 - 8.0_pr*c4)*u4)*Log(1.0_pr + u2)))/(110880.0_pr*fpi4*Pi*u8)

          darhorho=darhorho+(gA2*mpi**3*(1680.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)*u2 &
               -3.0_pr*(620862.0_pr*c1+351039.0_pr*c3-250970.0_pr*c4)*u4+(71742.0_pr*c1-56163.0_pr*c3 &
               +30482.0_pr*c4)*u6+18.0_pr*(286.0_pr*c1-127.0_pr*(3.0_pr*c3+c4))*u8+1710.0_pr*(3.0_pr*c3 &
               +c4)*u10+3.0_pr*u3*(1212750.0_pr*c1-6237.0_pr*(2.0_pr*c1-c3)*u4+572.0_pr*(3.0_pr*c1 &
               -3.0_pr*c3-c4)*u6+570.0_pr*(3.0_pr*c3+c4)*u8+1155.0_pr*(3.0_pr*c3-2.0_pr*c4)*(175.0_pr &
               +6.0_pr*u2))*Atan(u)+24.0_pr*(-70.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)-528.0_pr*(162.0_pr*c1 &
               +69.0_pr*c3-40.0_pr*c4)*u2+891.0_pr*(14.0_pr*c1+11.0_pr*c3-8.0_pr*c4)*u4)*Log(1.0_pr+u2))) &
               /(110880.0_pr*fpi4*Pi*u9)

          ddarhorho=ddarhorho+(gA2*mpi**3*(-2520.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)*u2 &
               +21.0_pr*(79002.0_pr*c1+45789.0_pr*c3-32950.0_pr*c4)*u4-9.0_pr*(4774.0_pr*c1-3657.0_pr*c3 &
               +1630.0_pr*c4)*u6-190.0_pr*(3.0_pr*c3+c4)*u8+570.0_pr*(3.0_pr*c3+c4)*u10 &
               +3.0_pr*(u3*(-1212750.0_pr*c1-2079.0_pr*(-2.0_pr*c1+c3)*u4+190.0_pr*(3.0_pr*c3+c4)*u8 &
               -1155.0_pr*(3.0_pr*c3-2.0_pr*c4)*(175.0_pr+4.0_pr*u2))*Atan(u)+4.0_pr*(210.0_pr*(330.0_pr*c1 &
               +129.0_pr*c3-56.0_pr*c4)+1232.0_pr*(162.0_pr*c1+69.0_pr*c3-40.0_pr*c4)*u2-1485.0_pr*(14.0_pr*c1 &
               +11.0_pr*c3-8.0_pr*c4)*u4)*Log(1.0_pr+u2))))/(18480.0_pr*fpi4*Pi*u10)

          brhorho=brhorho +(c4*gA2*mpi**3*(u2*(5880.0_pr - 81975.0_pr*u2 - 7049.0_pr*u4 + 1357.0_pr*u6 &
               + 285.0_pr*u8) + 3.0_pr*u3*(40425.0_pr + 2310.0_pr*u2 - 286.0_pr*u6 + 95.0_pr*u8)*Atan(u) &
               + 24.0_pr*(-245.0_pr - 1760.0_pr*u2 + 891.0_pr*u4)*Log(1.0_pr + u2)))/(27720.0_pr*fpi4*Pi*u8)

          dbrhorho=dbrhorho+(c4*gA2*mpi**3*(u2*(-47040.0_pr+376455.0_pr*u2+15241.0_pr*u4 &
               -1143.0_pr*u6+855.0_pr*u8)+3.0_pr*u3*(-202125.0_pr-6930.0_pr*u2-286.0_pr*u6 &
               +285.0_pr*u8)*Atan(u)-96.0_pr*(-490.0_pr-2640.0_pr*u2+891.0_pr*u4)*Log(1.0_pr+u2))) &
               /(27720.0_pr*fpi4*Pi*u9)

          ddbrhorho=ddbrhorho+(c4*gA2*mpi**3*(u2*(14112.0_pr-69195.0_pr*u2-1467.0_pr*u4-19.0_pr*u6 &
               +57.0_pr*u8)+3.0_pr*u3*(40425.0_pr+924.0_pr*u2+19.0_pr*u8)*Atan(u)+48.0_pr*(-294.0_pr &
               -1232.0_pr*u2+297.0_pr*u4)*Log(1.0_pr+u2)))/(924.0_pr*fpi4*Pi*u10)

          arhodrho=arhodrho+(gA2*mpi*(10395.0_pr*u3*(4.0_pr+u2)*(8.0_pr*c1*(40.0_pr-10.0_pr*u2+u4) &
               +c3*(160.0_pr+8.0_pr*u2+4.0_pr*u4-3.0_pr*u6))*Atan(0.5_pr*u)-2.0_pr*(480.0_pr*(330.0_pr*c1 &
               +129.0_pr*c3-56.0_pr*c4)*u2+48.0_pr*(49038.0_pr*c1+23511.0_pr*c3+7870.0_pr*c4)*u4 &
               +32.0_pr*(-7953.0_pr*c1+10545.0_pr*c3+479.0_pr*c4)*u6-12.0_pr*(5874.0_pr*c1-2769.0_pr*c3 &
               +232.0_pr*c4)*u8+27.0_pr*(147.0_pr*c3+16.0_pr*c4)*u10-48.0_pr*u3*(-5775.0_pr*(6.0_pr*c1 &
               +3.0_pr*c3-2.0_pr*c4)+297.0_pr*(2.0_pr*c1-c3)*u4+88.0_pr*(-3.0_pr*c1+3.0_pr*c3+c4)*u6 &
               +90.0_pr*(3.0_pr*c3+c4)*u8)*Atan(u)+48.0_pr*(-1290.0_pr*c3+560.0_pr*c4+33.0_pr*(-100.0_pr*c1 &
               -3.0_pr*(162.0_pr*c1+69.0_pr*c3-40.0_pr*c4)*u2+9.0_pr*(14.0_pr*c1+11.0_pr*c3 &
               -8.0_pr*c4)*u4))*Log(1.0_pr+u2))))/(1.216512e6_pr*fpi4*Pi*u10)

          darhodrho=darhodrho+(gA2*mpi*(-10395.0_pr*u3*(8960.0_pr*c1+72.0_pr*(-2.0_pr*c1+c3)*u4 &
               +8.0_pr*(c1-c3)*u6+3.0_pr*c3*u8+320.0_pr*c3*(14.0_pr+3.0_pr*u2))*Atan(0.5_pr*u) &
               +2.0_pr*(4800.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)*u2+48.0_pr*(354354.0_pr*c1+169473.0_pr*c3 &
               +51970.0_pr*c4)*u4+8.0_pr*(-226050.0_pr*c1+177369.0_pr*c3+7676.0_pr*c4)*u6 &
               -12.0_pr*(5874.0_pr*c1-8979.0_pr*c3+472.0_pr*c4)*u8+135.0_pr*(-135.0_pr*c3+32.0_pr*c4)*u10 &
               +48.0_pr*u3*(40425.0_pr*(6.0_pr*c1+3.0_pr*c3-2.0_pr*c4)+891.0_pr*(-2.0_pr*c1+c3)*u4 &
               +88.0_pr*(3.0_pr*c1-3.0_pr*c3-c4)*u6+90.0_pr*(3.0_pr*c3+c4)*u8)*Atan(u) &
               +96.0_pr*(-50.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)-396.0_pr*(162.0_pr*c1+69.0_pr*c3 &
               -40.0_pr*c4)*u2+891.0_pr*(14.0_pr*c1+11.0_pr*c3-8.0_pr*c4)*u4)*Log(1.0_pr+u2)))) &
               /(1.216512e6_pr*fpi4*Pi*u11)

          ddarhodrho=ddarhodrho+(gA2*mpi*(-2400.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)*u2 &
               -24.0_pr*(273126.0_pr*c1+129687.0_pr*c3+34280.0_pr*c4)*u4-4.0_pr*(229530.0_pr*c1+264813.0_pr*c3 &
               +56912.0_pr*c4)*u6+16.0_pr*(8052.0_pr*c1-6513.0_pr*c3-188.0_pr*c4)*u8+6.0_pr*(267.0_pr*(c1 &
               -c3)+16.0_pr*c4)*u10+3.0_pr*(4.0_pr+u2)*(315.0_pr*u3*(2240.0_pr*(2.0_pr*c1+c3) &
               +360.0_pr*c3*u2+18.0_pr*(-2.0_pr*c1+c3)*u4+(c1-c3)*u6)*Atan(0.5_pr*u) &
               -4.0_pr*(2.0_pr*u3*(7350.0_pr*(6.0_pr*c1+3.0_pr*c3-2.0_pr*c4)+81.0_pr*(-2.0_pr*c1+c3)*u4 &
               -4.0_pr*(-3.0_pr*c1+3.0_pr*c3+c4)*u6)*Atan(u)+(-50.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4) &
               -324.0_pr*(162.0_pr*c1+69.0_pr*c3-40.0_pr*c4)*u2+567.0_pr*(14.0_pr*c1+11.0_pr*c3 &
               -8.0_pr*c4)*u4)*Log(1.0_pr+u2)))))/(6912.0_pr*fpi4*Pi*u12*(4.0_pr+u2))

          brhodrho=brhodrho+(gA2*mpi*(-99.0_pr*u3*(33600.0_pr*(2.0_pr*c1+c3)*u+560.0_pr*(-10.0_pr*c1 &
               +13.0_pr*c3)*u3+840.0_pr*(-2.0_pr*c1+c3)*u5+54.0_pr*c3*u7-105.0_pr*(4.0_pr+u2) &
               *(8.0_pr*c1*(40.0_pr-10.0_pr*u2+u4)+c3*(160.0_pr+8.0_pr*u2+4.0_pr*u4-3.0_pr*u6)) &
               *Atan(0.5_pr*u))+64.0_pr*c4*(u2*(1680.0_pr-23610.0_pr*u2-958.0_pr*u4+174.0_pr*u6 &
               -27.0_pr*u8)+6.0_pr*u3*(5775.0_pr+44.0_pr*u6+45.0_pr*u8)*Atan(u)+24.0_pr*(-70.0_pr &
               +99.0_pr*u2*(-5.0_pr+3.0_pr*u2))*Log(1.0_pr+u2))))/(1.216512e6_pr*fpi4*Pi*u10)

          dbrhodrho=dbrhodrho+(gA2*mpi*(2.0_pr*u2*(11642400.0_pr*(2.0_pr*c1+c3)*u2+64.0_pr*c4*( &
               -8400.0_pr+77955.0_pr*u2+1919.0_pr*u4-177.0_pr*u6+135.0_pr*u8)-3465.0_pr*u4 &
               *(8.0_pr*c1*(70+3*u2)+c3*(-440.0_pr+9.0_pr*u2*(-4.0_pr+u2))))-3.0_pr*(3465.0_pr*u3 &
               *(8960.0_pr*c1+72.0_pr*(-2.0_pr*c1+c3)*u4+8.0_pr*(c1-c3)*u6+3.0_pr*c3*u8+320.0_pr*c3 &
               *(14.0_pr+3.0_pr*u2))*Atan(0.5_pr*u)-128.0_pr*c4*(u3*(-40425.0_pr-44.0_pr*u6+45.0_pr*u8) &
               *Atan(u)+8.0_pr*(350.0_pr+99.0_pr*u2*(20.0_pr-9.0_pr*u2))*Log(1.0_pr+u2))))) &
               /(1.216512e6_pr*fpi4*Pi*u11)

          ddbrhodrho=ddbrhodrho+(gA2*mpi*(2.0_pr*u2*(32.0_pr*c4*(4.0_pr+u2)*(1050.0_pr-6690.0_pr*u2 &
               -106.0_pr*u4+3.0_pr*u6)+315.0_pr*u2*(-6720.0_pr*(2.0_pr*c1+c3)-40.0_pr*(56.0_pr*c1 &
               +55.0_pr*c3)*u2+2.0_pr*(110.0_pr*c1-89.0_pr*c3)*u4+3.0_pr*(c1-c3)*u6))+3.0_pr*(4.0_pr &
               +u2)*(315.0_pr*u3*(2240.0_pr*(2.0_pr*c1+c3)+360.0_pr*c3*u2+18.0_pr*(-2.0_pr*c1+c3)*u4+(c1 &
               -c3)*u6)*Atan(0.5_pr*u)+64.0_pr*c4*(u3*(3675.0_pr+u6)*Atan(u)+(-350.0_pr+81.0_pr*u2 &
               *(-20.0_pr+7.0_pr*u2))*Log(1.0_pr+u2)))))/(6912.0_pr*fpi4*Pi*u12*(4.0_pr+u2))

          adrdr=adrdr+0.0_pr
          dadrdr=dadrdr+0.0_pr
          bdrdr=bdrdr+0.0_pr
          dbdrdr=dbdrdr+0.0_pr

          arhotau=arhotau+(gA2*mpi*(60.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)*u2-6.0_pr*(20262.0_pr*c1 &
               +11139.0_pr*c3-7870.0_pr*c4)*u4+(2838.0_pr*c1-2865.0_pr*c3+1916.0_pr*c4)*u6 &
               +12.0_pr*(132.0_pr*c1-29.0_pr*(3.0_pr*c3+c4))*u8+54.0_pr*(3.0_pr*c3+c4)*u10 &
               -6.0_pr*u3*(-5775.0_pr*(6.0_pr*c1+3.0_pr*c3-2.0_pr*c4)+297.0_pr*(2.0_pr*c1-c3)*u4 &
               +88.0_pr*(-3.0_pr*c1+3.0_pr*c3+c4)*u6+90.0_pr*(3.0_pr*c3+c4)*u8)*Atan(u) &
               +6.0_pr*(-1290.0_pr*c3+560.0_pr*c4+33.0_pr*(-100.0_pr*c1-3.0_pr*(162.0_pr*c1+69.0_pr*c3 &
               -40.0_pr*c4)*u2+9.0_pr*(14.0_pr*c1+11.0_pr*c3-8.0_pr*c4)*u4))*Log(1.0_pr+u2))) &
               /(19008.0_pr*fpi4*Pi*u10)

          darhotau=darhotau -(gA2*mpi*(300.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)*u2 &
               -3.0_pr*(130746.0_pr*c1+73077.0_pr*c3-51970.0_pr*c4)*u4+(8250.0_pr*c1-6603.0_pr*c3 &
               +3838.0_pr*c4)*u6+6.0_pr*(132.0_pr*c1-59.0_pr*(3.0_pr*c3+c4))*u8+270.0_pr*(3.0_pr*c3 &
               +c4)*u10+3.0_pr*u3*(40425.0_pr*(6.0_pr*c1+3.0_pr*c3-2.0_pr*c4)+891.0_pr*(-2.0_pr*c1 &
               +c3)*u4+88.0_pr*(3.0_pr*c1-3.0_pr*c3-c4)*u6+90.0_pr*(3.0_pr*c3+c4)*u8)*Atan(u) &
               +6.0_pr*(-50.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)-396.0_pr*(162.0_pr*c1+69.0_pr*c3 &
               -40.0_pr*c4)*u2+891.0_pr*(14.0_pr*c1+11.0_pr*c3-8.0_pr*c4)*u4)*Log(1.0_pr+u2))) &
               /(9504.0_pr*fpi4*Pi*u11)

          ddarhotau=ddarhotau+(gA2*mpi*(150.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)*u2 &
               -3.0_pr*(43962.0_pr*c1+24969.0_pr*c3-17840.0_pr*c4)*u4+(2154.0_pr*c1-1695.0_pr*c3 &
               +848.0_pr*c4)*u6-24.0_pr*(-3.0_pr*c1+3.0_pr*c3+c4)*u8+6.0_pr*u3*(7350.0_pr*(6.0_pr*c1 &
               +3.0_pr*c3-2.0_pr*c4)+81.0_pr*(-2.0_pr*c1+c3)*u4-4.0_pr*(-3.0_pr*c1+3.0_pr*c3+c4)*u6) &
               *Atan(u)+3.0_pr*(-50.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)-324.0_pr*(162.0_pr*c1+69.0_pr*c3 &
               -40.0_pr*c4)*u2+567.0_pr*(14.0_pr*c1+11.0_pr*c3-8.0_pr*c4)*u4)*Log(1.0_pr+u2))) &
               /(432.0_pr*fpi4*Pi*u12)

          brhotau=brhotau+(c4*gA2*mpi*(u2*(-1680.0_pr+23610.0_pr*u2+958.0_pr*u4-174.0_pr*u6 &
               +27.0_pr*u8)-6.0_pr*u3*(5775.0_pr+44.0_pr*u6+45.0_pr*u8)*Atan(u)+24.0_pr*(70.0_pr &
               +99.0_pr*u2*(5.0_pr-3.0_pr*u2))*Log(1.0_pr+u2)))/(4752.0_pr*fpi4*Pi*u10)

          dbrhotau=dbrhotau+(c4*gA2*mpi*(u2*(8400.0_pr-77955.0_pr*u2-1919.0_pr*u4+177.0_pr*u6 &
               -135.0_pr*u8)+3.0_pr*u3*(40425.0_pr+44.0_pr*u6-45.0_pr*u8)*Atan(u)+24.0_pr*(-350.0_pr &
               +99.0_pr*u2*(-20.0_pr+9.0_pr*u2))*Log(1.0_pr+u2)))/(2376.0_pr*fpi4*Pi*u11)

          ddbrhotau=ddbrhotau+(c4*gA2*mpi*(-1050.0_pr*u2+6690.0_pr*u4+106.0_pr*u6-3.0_pr*u8 &
               -3.0_pr*u3*(3675.0_pr+u6)*Atan(u)+3.0_pr*(350.0_pr+81.0_pr*u2*(20.0_pr-7.0_pr*u2)) &
               *Log(1.0_pr+u2)))/(27.0_pr*fpi4*Pi*u12)

          ajj=ajj+(gA2*mpi*(u2*(-84.0_pr*(c1+2.0_pr*c1*u2)+6.0_pr*c3*(-11.0_pr+16.0_pr*u2+u4) &
               -c4*(12.0_pr+11.0_pr*u2*(2.0_pr+u2)))+2.0_pr*u3*(210.0_pr*c1-105.0_pr*c3 &
               +28.0_pr*(-3.0_pr*c1+3.0_pr*c3+c4)*u2+12.0_pr*(9.0_pr*c3+c4)*u4)*Atan(u)+(84.0_pr*c1 &
               +66.0_pr*c3+12.0_pr*c4+7.0_pr*(-30.0_pr*c1+21.0_pr*c3+4.0_pr*c4)*u2)*Log(1.0_pr+u2))) &
               /(4480.0_pr*fpi4*Pi*u6)

          dajj=dajj+(gA2*mpi*(u3*(315.0_pr*(-2.0_pr*c1+c3)+28.0_pr*(3.0_pr*c1-3.0_pr*c3-c4)*u2 &
               +12.0_pr*(9.0_pr*c3+c4)*u4)*Atan(u)+2.0_pr*(9.0_pr*(14.0_pr*c1+11.0_pr*c3+2.0_pr*c4)*u2 &
               +(42.0_pr*c1-60.0_pr*c3+19.0_pr*c4)*u4+6.0_pr*(9.0_pr*c3+c4)*u6+(-9.0_pr*(14.0_pr*c1 &
               +11.0_pr*c3+2.0_pr*c4)+7.0_pr*(30.0_pr*c1-21.0_pr*c3-4.0_pr*c4)*u2)*Log(1.0_pr+u2)))) &
               /(2240.0_pr*fpi4*Pi*u7)

          bjj=bjj+(c4*gA2*mpi*(-(u2*(12.0_pr+11.0_pr*u2*(2.0_pr+u2)))+8.0_pr*u5*(7.0_pr &
               +3.0_pr*u2)*Atan(u)+4.0_pr*(3.0_pr+7.0_pr*u2)*Log(1.0_pr+u2)))/(2240.0_pr*fpi4*Pi*u6)

          dbjj=dbjj+(c4*gA2*mpi*(u2*(18.0_pr+19.0_pr*u2+6.0_pr*u4)+2.0_pr*u5*(-7.0_pr &
               +3.0_pr*u2)*Atan(u)-2.0_pr*(9.0_pr+14.0_pr*u2)*Log(1.0_pr+u2)))/(560.0_pr*fpi4*Pi*u7)

       Else  !! if (u.gt.ucut...
          arhorho=arhorho +(3.0_pr*(2.0_pr*c1 - c3)*gA2*mpi**3)/(64.0_pr*fpi4*Pi) + (11.0_pr*(42.0_pr*c1 &
               + 129.0_pr*c3 + 40.0_pr*c4)*gA2*mpi**3*u4)/(100800.0_pr*fpi4*Pi) - ((162.0_pr*c1 &
               + 213.0_pr*c3 + 56.0_pr*c4)*gA2*mpi**3*u6)/(69300.0_pr*fpi4*Pi)
          darhorho=darhorho+(gA2*mpi**3*u3*(121.0_pr*(42.0_pr*c1+129.0_pr*c3+40.0_pr*c4) &
               -24.0_pr*(162.0_pr*c1+213.0_pr*c3+56.0_pr*c4)*u2))/(277200.0_pr*fpi4*Pi)
          ddarhorho=ddarhorho+(gA2*mpi**3*u2*(1573.0_pr*(42.0_pr*c1+129.0_pr*c3+40.0_pr*c4) &
               -520.0_pr*(162.0_pr*c1+213.0_pr*c3+56.0_pr*c4)*u2+765.0_pr*(110.0_pr*c1+107.0_pr*c3 &
               +24.0_pr*c4)*u4))/(1.2012e6_pr*fpi4*Pi)

          brhorho=brhorho +(11.0_pr*c4*gA2*mpi**3*u4)/(1260.0_pr*fpi4*Pi) &
               - (4.0_pr*c4*gA2*mpi**3*u6)/(2475.0_pr*fpi4*Pi)
          dbrhorho=dbrhorho+(c4*gA2*mpi**3*u3*(605.0_pr-168.0_pr*u2))/(17325.0_pr*fpi4*Pi)
          ddbrhorho=ddbrhorho+(c4*gA2*mpi**3*u2*(1573.0_pr-728.0_pr*u2+459.0_pr*u4)) &
               /(15015.0_pr*fpi4*Pi)

          arhodrho=arhodrho+(gA2*mpi*(1081080.0_pr*(10.0_pr*c1-11.0_pr*c3+4.0_pr*c4)+16016.0_pr &
               *(42.0_pr*c1+129.0_pr*c3+140.0_pr*c4)*u2-117.0_pr*(5886.0_pr*c1+7739.0_pr*c3+2688.0_pr*c4)*u4 &
               +135.0_pr*(2970.0_pr*c1+2889.0_pr*c3+704.0_pr*c4)*u6))/(5.5351296e8_pr*fpi4*Pi)
          darhodrho=darhodrho+(gA2*mpi*u*(16016.0_pr*(42.0_pr*c1+129.0_pr*c3+140.0_pr*c4)-234.0_pr &
               *(5886.0_pr*c1+7739.0_pr*c3+2688.0_pr*c4)*u2+405.0_pr*(2970.0_pr*c1+2889.0_pr*c3 &
               +704.0_pr*c4)*u4))/(2.7675648e8_pr*fpi4*Pi)
          ddarhodrho=ddarhodrho+(gA2*mpi*(192192.0_pr*(42.0_pr*c1+129.0_pr*c3+140.0_pr*c4)-8424.0_pr &
               *(5886.0_pr*c1+7739.0_pr*c3+2688.0_pr*c4)*u2+24300.0_pr*(2970.0_pr*c1+2889.0_pr*c3 &
               +704.0_pr*c4)*u4-175.0_pr*(442806.0_pr*c1+367383.0_pr*c3+73216.0_pr*c4)*u6)) &
               /(3.32107776e9_pr*fpi4*Pi)

          brhodrho=brhodrho+(gA2*mpi*(2162160.0_pr*(10.0_pr*c1-11.0_pr*c3+4*c4)-40040.0_pr*(42.0_pr*c1 &
               +129.0_pr*c3-112.0_pr*c4)*u2+819.0_pr*(270.0_pr*c1+355.0_pr*c3-768.0_pr*c4)*u4 &
               -45.0_pr*(770.0_pr*c1+749.0_pr*c3-4224.0_pr*c4)*u6))/(5.5351296e8_pr*fpi4*Pi)
          dbrhodrho=dbrhodrho+(gA2*mpi*u*(-40040.0_pr*(42.0_pr*c1+129.0_pr*c3-112.0_pr*c4)+1638.0_pr &
               *(270.0_pr*c1+355.0_pr*c3-768.0_pr*c4)*u2-135.0_pr*(770.0_pr*c1+749.0_pr*c3 &
               -4224.0_pr*c4)*u4))/(2.7675648e8_pr*fpi4*Pi)
          ddbrhodrho=ddbrhodrho+(gA2*mpi*(-480480.0_pr*(42.0_pr*c1+129.0_pr*c3-112.0_pr*c4)+58968.0_pr &
               *(270.0_pr*c1+355.0_pr*c3-768.0_pr*c4)*u2-8100.0_pr*(770.0_pr*c1+749.0_pr*c3 &
               -4224.0_pr*c4)*u4+175.0_pr*(11466.0_pr*c1+9513.0_pr*c3-146432.0_pr*c4)*u6)) &
               /(3.32107776e9_pr*fpi4*Pi)

          adrdr=adrdr+0.0_pr
          dadrdr=dadrdr+0.0_pr
          bdrdr=bdrdr+0.0_pr
          dbdrdr=dbdrdr+0.0_pr

          arhotau=arhotau+(gA2*mpi*(135135.0_pr*(10.0_pr*c1-11.0_pr*c3-4.0_pr*c4)-7007.0_pr*(42.0_pr*c1 &
               +129.0_pr*c3+40.0_pr*c4)*u2+702.0_pr*(162.0_pr*c1+213.0_pr*c3+56.0_pr*c4)*u4-495.0_pr &
               *(110.0_pr*c1+107.0_pr*c3+24.0_pr*c4)*u6))/(1.729728e7_pr*fpi4*Pi)
          darhotau=darhotau+(gA2*mpi*u*(-7007.0_pr*(42.0_pr*c1+129.0_pr*c3+40.0_pr*c4)+1404.0_pr &
               *(162.0_pr*c1+213.0_pr*c3+56.0_pr*c4)*u2-1485.0_pr*(110.0_pr*c1+107.0_pr*c3+24.0_pr*c4)*u4)) &
               /(8.64864e6_pr*fpi4*Pi)
          ddarhotau=ddarhotau+(gA2*mpi*(-21021.0_pr*(42.0_pr*c1+129.0_pr*c3+40.0_pr*c4)+12636.0_pr &
               *(162.0_pr*c1+213.0_pr*c3+56.0_pr*c4)*u2-22275.0_pr*(110.0_pr*c1+107.0_pr*c3+24.0_pr*c4)*u4 &
               +4550.0_pr*(546.0_pr*c1+453.0_pr*c3+88.0_pr*c4)*u6))/(2.594592e7_pr*fpi4*Pi)

          brhotau=brhotau -(c4*gA2*mpi*(135135.0_pr+70070.0_pr*u2-9828.0_pr*u4+2970.0_pr*u6)) &
               /(2.16216e6_pr*fpi4*Pi)
          dbrhotau=dbrhotau -(c4*gA2*mpi*u*(35035.0_pr-9828.0_pr*u2+4455.0_pr*u4))/(540540.0_pr*fpi4*Pi)
          ddbrhotau=ddbrhotau+(c4*gA2*mpi*(-105105.0_pr+88452.0_pr*u2-66825.0_pr*u4+50050.0_pr*u6)) &
               /(1.62162e6_pr*fpi4*Pi)

          ajj=ajj+(gA2*mpi*(-1050.0_pr*(10.0_pr*c1-11.0_pr*c3-2.0_pr*c4)+70.0_pr*(42.0_pr*c1 &
               +129.0_pr*c3+10.0_pr*c4)*u2-3.0_pr*(486.0_pr*c1+639.0_pr*c3+28.0_pr*c4)*u4+8.0_pr &
               *(110.0_pr*c1+107.0_pr*c3+3.0_pr*c4)*u6))/(268800.0_pr*fpi4*Pi)
          dajj=dajj+(gA2*mpi*u*(35.0_pr*(42.0_pr*c1+129.0_pr*c3+10.0_pr*c4)-3.0_pr*(486.0_pr*c1 &
               +639.0_pr*c3+28.0_pr*c4)*u2+12.0_pr*(110.0_pr*c1+107.0_pr*c3+3.0_pr*c4)*u4)) &
               /(67200.0_pr*fpi4*Pi)

          bjj=bjj+(c4*gA2*mpi*(525.0_pr+175.0_pr*u2-21.0_pr*u4+6.0_pr*u6))/(33600.0_pr*fpi4*Pi)
          dbjj=dbjj+(c4*gA2*mpi*u*(175.0_pr-42.0_pr*u2+18.0_pr*u4))/(16800.0_pr*fpi4*Pi)
       End If  !! if (u.gt.ucut...
       !
       ! N2LO, 3N terms
       If (use_DME3N_terms) Then
          If (u.Gt.ucut3n) Then
             y=1.0_pr/u
             dy=-y*y
             ddy=2.0_pr*y*y*y
             marc=atu
             dmarc=2.0_pr/(1.0_pr+ 4.0_pr*u2)
             ddmarc=-16.0_pr*u/((1.0_pr+ 4.0_pr*u2)**2)
             mlog=lual
             dmlog=8.0_pr*u/(1.0_pr+ 4.0_pr*u2)
             ddmlog=(8.0_pr-32.0_pr*u2)/((1.0_pr+ 4.0_pr*u2)**2)
             Do l=1,33
                Do j=1,3
                   Do k=1,3
                      aux=((y**(l-1))*(atu**(j-1))*(lual**(k-1)))*mevfm
                      daux=(Real(-1+j,pr)*marc**(-2+j)*mlog**(-1+k)*y**(-1+l)*dmarc &     
                           +Real(-1+k,pr)*marc**(-1+j)*mlog**(-2+k)*y**(-1+l)*dmlog &
                           +Real(-1+l,pr)*marc**(-1+j)*mlog**(-1+k)*y**(-2+l)*dy)*mevfm
                      ddaux=((2.0_pr*Real(-1+l,pr)*y**(-2+l)*(Real(-1+j,pr) &
                           *marc**(-2+j)*mlog**(-1+k)*dmarc+Real(-1+k,pr)*marc**(-1+j)*mlog**(-2+k) &
                           *dmlog)*dy+y**(-1+l)*(2.0_pr*Real((-1+j)*(-1+k),pr)*marc**(-2+j)*mlog**(-2+k) &
                           *dmarc*dmlog+mlog**(-1+k)*(Real((-2+j)*(-1+j),pr)*marc**(-3+j)*dmarc**2 &
                           +Real(-1+j,pr)*marc**(-2+j)*ddmarc)+marc**(-1+j)*(Real((-2+k)*(-1+k),pr) &
                           *mlog**(-3+k)*dmlog**2+Real(-1+k,pr)*mlog**(-2+k)*ddmlog))+marc**(-1+j) &
                           *mlog**(-1+k)*(Real((-2+l)*(-1+l),pr)*y**(-3+l)*dy**2+Real(-1+l,pr) &
                           *y**(-2+l)*ddy)))*mevfm
                      !
                      hrho0rho0 =hrho0rho0 +ctr0r0(k,j,l)*aux
                      hrho1rho1 =hrho1rho1 +ctr1r1(k,j,l)*aux
                      hdr0dr0   =hdr0dr0+ctdr0dr0(k,j,l)*aux
                      hdr1dr1   =hdr1dr1+ctdr1dr1(k,j,l)*aux
                      hrho0Drho0=hrho0Drho0+ctr0Dr0(k,j,l)*aux
                      hrho1Drho0=hrho1Drho0+ctr1Dr0(k,j,l)*aux
                      hrho1Drho1=hrho1Drho1+ctr1Dr1(k,j,l)*aux
                      hrho0tau0 =hrho0tau0+ctr0t0(k,j,l)*aux
                      hrho1tau0 =hrho1tau0+ctr1t0(k,j,l)*aux
                      hrho1tau1 =hrho1tau1+ctr1t1(k,j,l)*aux
                      hJ0dr0    =hJ0dr0+ctJ0dr0(k,j,l)*aux
                      hrho0DJ0  =hrho0DJ0+ctr0dJ0(k,j,l)*aux
                      hJ1dr1    =hJ1dr1+ctJ1dr1(k,j,l)*aux
                      hrho1DJ1  =hrho1DJ1+ctr1dJ1(k,j,l)*aux
                      hJ0dr1    =hJ0dr1+ctJ0dr1(k,j,l)*aux
                      hrho1DJ0  =hrho1DJ0+ctr1dJ0(k,j,l)*aux
                      hJ1dr0    =hJ1dr0+ctJ1dr0(k,j,l)*aux
                      hJ0J0     =hJ0J0+ctJ0J0(k,j,l)*aux
                      hJ0J1     =hJ0J1+ctJ0J1(k,j,l)*aux        
                      hJ1J1     =hJ1J1+ctJ1J1(k,j,l)*aux
                      !
                      dhrho0rho0 =dhrho0rho0 +ctr0r0(k,j,l)*daux
                      dhrho1rho1 =dhrho1rho1 +ctr1r1(k,j,l)*daux
                      dhdr0dr0   =dhdr0dr0+ctdr0dr0(k,j,l)*daux
                      dhdr1dr1   =dhdr1dr1+ctdr1dr1(k,j,l)*daux
                      dhrho0Drho0=dhrho0Drho0+ctr0Dr0(k,j,l)*daux
                      dhrho1Drho0=dhrho1Drho0+ctr1Dr0(k,j,l)*daux
                      dhrho1Drho1=dhrho1Drho1+ctr1Dr1(k,j,l)*daux
                      dhrho0tau0 =dhrho0tau0+ctr0t0(k,j,l)*daux
                      dhrho1tau0 =dhrho1tau0+ctr1t0(k,j,l)*daux
                      dhrho1tau1 =dhrho1tau1+ctr1t1(k,j,l)*daux
                      dhJ0dr0    =dhJ0dr0+ctJ0dr0(k,j,l)*daux
                      dhrho0DJ0  =dhrho0DJ0+ctr0dJ0(k,j,l)*daux
                      dhJ1dr1    =dhJ1dr1+ctJ1dr1(k,j,l)*daux
                      dhrho1DJ1  =dhrho1DJ1+ctr1dJ1(k,j,l)*daux
                      dhJ0dr1    =dhJ0dr1+ctJ0dr1(k,j,l)*daux
                      dhrho1DJ0  =dhrho1DJ0+ctr1dJ0(k,j,l)*daux
                      dhJ1dr0    =dhJ1dr0+ctJ1dr0(k,j,l)*daux
                      dhJ0J0     =dhJ0J0+ctJ0J0(k,j,l)*daux
                      dhJ0J1     =dhJ0J1+ctJ0J1(k,j,l)*daux        
                      dhJ1J1     =dhJ1J1+ctJ1J1(k,j,l)*daux
                      !
                      ddhrho0rho0 =ddhrho0rho0 +ctr0r0(k,j,l)*ddaux
                      ddhrho1rho1 =ddhrho1rho1 +ctr1r1(k,j,l)*ddaux
                      ddhrho0Drho0=ddhrho0Drho0+ctr0Dr0(k,j,l)*ddaux
                      ddhrho1Drho0=ddhrho1Drho0+ctr1Dr0(k,j,l)*ddaux
                      ddhrho1Drho1=ddhrho1Drho1+ctr1Dr1(k,j,l)*ddaux
                      ddhrho0tau0 =ddhrho0tau0+ctr0t0(k,j,l)*ddaux
                      ddhrho1tau0 =ddhrho1tau0+ctr1t0(k,j,l)*ddaux
                      ddhrho1tau1 =ddhrho1tau1+ctr1t1(k,j,l)*ddaux
                   End Do
                End Do
             End Do  !! end of l,j,k loops
          Else  !! if (u.gt.ucut3n..
             hrho0rho0=mevfm*(2.0555896e-8_pr*ce*(-15.058202_pr + u)*(13.851257_pr + u)*(218.02762_pr &
                  - 15.661674_pr*u + u2)*(200.58157_pr + 13.247785_pr*u + u2) + gA*(-0.042893802_pr*cd*(-1.5948724_pr &
                  + u)*(-0.01657426_pr + u)*(1.8884833_pr - 2.1431293_pr*u + u2)*(0.029641853_pr - 0.28437579_pr*u &
                  + u2) + gA*LambdaX*(6.0619551_pr*c3*(-1.1935357_pr + u)*(-0.30214838_pr + u)*(0.9281019_pr &
                  - 1.5964577_pr*u + u2)*(0.15916821_pr - 0.63445995_pr*u + u2) + 3.9928539_pr*c1*(-1.2431741_pr &
                  + u)*(-0.19155888_pr + u)*(1.1456958_pr - 1.8138841_pr*u + u2)*(0.19109424_pr - 0.36771486_pr*u &
                  + u2) - 0.06402411_pr*c4*(-1.3982905_pr + u)*(0.1480126_pr + u)*(1.6779335_pr - 2.2114154_pr*u &
                  + u2)*(0.041603529_pr - 0.35726143_pr*u + u2))))/(fpi4*LambdaX)

             hrho1rho1=mevfm*(2.0555896e-8_pr*ce*(-15.058202_pr + u)*(13.851257_pr + u)*(218.02762_pr &
                  - 15.661674_pr*u + u2)*(200.58157_pr + 13.247785_pr*u + u2) + gA*(-0.042893802_pr*cd*(-1.5948724_pr &
                  + u)*(-0.01657426_pr + u)*(1.8884833_pr - 2.1431293_pr*u + u2)*(0.029641853_pr - 0.28437579_pr*u &
                  + u2) + gA*LambdaX*(6.0619551_pr*c3*(-1.1935357_pr + u)*(-0.30214838_pr + u)*(0.9281019_pr &
                  - 1.5964577_pr*u + u2)*(0.15916821_pr - 0.63445995_pr*u + u2) + 3.9928539_pr*c1*(-1.2431741_pr &
                  + u)*(-0.19155888_pr + u)*(1.1456958_pr - 1.8138841_pr*u + u2)*(0.19109424_pr - 0.36771486_pr*u &
                  + u2) - 0.06402411_pr*c4*(-1.3982905_pr + u)*(0.1480126_pr + u)*(1.6779335_pr - 2.2114154_pr*u &
                  + u2)*(0.041603529_pr - 0.35726143_pr*u + u2))))/(fpi4*LambdaX)

             hdr0dr0=0.0_pr
             hdr1dr1=0.0_pr

             hrho0Drho0=mevfm*(gA*(0.12666864_pr*cd*(-1.2705742_pr + u)*(0.21652154_pr + u)*(1.4119605_pr &
                  - 2.1693922_pr*u + u2)*(0.89941341_pr - 1.1919922_pr*u + u2) + gA*LambdaX*(29.501949_pr*c1*( &
                  -1.0215543_pr + u)*(-0.15813736_pr + u)*(0.84479502_pr - 1.6998479_pr*u + u2)*(0.34169853_pr &
                  - 0.81578924_pr*u + u2) - 7.1108952_pr*c3*(1.1461026_pr - 2.0823268_pr*u + u2)*(0.62212787_pr &
                  - 1.1642433_pr*u + u2)*(0.074190685_pr - 0.4633705_pr*u + u2) + 0.18892679_pr*c4*(1.7714304_pr &
                  - 2.5914415_pr*u + u2)*(1.2438813_pr - 1.6683493_pr*u + u2)*(0.004127797_pr - 0.089415137_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             hrho1Drho0=mevfm*(gA*(-0.042222876_pr*cd*(-1.2705742_pr + u)*(0.21652155_pr + u)*(1.4119605_pr &
                  - 2.1693922_pr*u + u2)*(0.89941343_pr - 1.1919922_pr*u + u2) + gA*LambdaX*(-0.064338673_pr*c1*( &
                  -1.9302829_pr + u)*(0.24441857_pr + u)*(2.1104361_pr - 2.7720023_pr*u + u2)*(1.3417224_pr &
                  - 1.5814849_pr*u + u2) + 0.049739131_pr*c3*(1.8967607_pr - 2.6742154_pr*u + u2)*(1.316844_pr &
                  - 1.6740809_pr*u + u2)*(0.011061538_pr - 0.13317506_pr*u + u2) - 0.073720916_pr*c4*(1.6783675_pr &
                  - 2.5203342_pr*u + u2)*(1.1903356_pr - 1.6163511_pr*u + u2)*(0.0074741816_pr - 0.10943007_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             hrho1Drho1=mevfm*(gA*(-0.084445761_pr*cd*(-1.2705742_pr + u)*(0.21652154_pr + u)*(1.4119605_pr &
                  - 2.1693922_pr*u + u2)*(0.89941341_pr - 1.1919922_pr*u + u2) + gA*LambdaX*(-7.0472043_pr*c1*( &
                  -1.0402373_pr + u)*(-0.15610243_pr + u)*(0.87447681_pr - 1.7182557_pr*u + u2)*(0.33211898_pr &
                  - 0.74186157_pr*u + u2) + 8.0775475_pr*c3*(0.95424896_pr - 1.9113384_pr*u + u2)*(0.50394896_pr &
                  - 1.1611245_pr*u + u2)*(0.11279874_pr - 0.61804716_pr*u + u2) - 0.1259512_pr*c4*(1.7714304_pr &
                  - 2.5914415_pr*u + u2)*(1.2438813_pr - 1.6683493_pr*u + u2)*(0.004127797_pr - 0.089415137_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             hrho0tau0=mevfm*(ga*(-0.25333726_pr*cd*(-1.2705742_pr+u)*(0.21652155_pr+u)*(1.4119605_pr &
                  +(-2.1693922_pr+u)*u)*(0.89941343_pr+(-1.1919922_pr+u)*u)+ga*LambdaX*(-4.0039062_pr*c1 &
                  *(-1.2352256_pr+u)*(0.064279963_pr+u)*(1.2516998_pr+(-1.9866561_pr+u)*u)*(0.62481483_pr &
                  +(-0.6277641_pr+u)*u)-0.43884277_pr*c4*(1.682325_pr+(-2.523347_pr+u)*u)*(1.1930406_pr &
                  +(-1.6187124_pr+u)*u)*(0.0072756559_pr+(-0.10828832_pr+u)*u)-1.3046875_pr*c3*(-1.3097655_pr &
                  +u)*(-0.12575226_pr+u)*(1.337577_pr+(-1.9544714_pr+u)*u)*(0.15931999_pr+u*(0.30346213_pr+u))))) &
                  /(fpi4*LambdaX*mpi2)

             hrho1tau0=mevfm*(ga*(0.084445755_pr*cd*(-1.2705742_pr+u)*(0.21652155_pr+u)*(1.4119605_pr &
                  +(-2.1693922_pr+u)*u)*(0.89941343_pr+(-1.1919922_pr+u)*u)+ga*LambdaX*(0.12695312_pr*c1 &
                  *(-1.9629389_pr+u)*(0.24499103_pr+u)*(2.1016515_pr+(-2.7692262_pr+u)*u)*(1.3407016_pr &
                  +(-1.5830183_pr+u)*u)-0.10644531_pr*c3*(1.8516394_pr+(-2.6412277_pr+u)*u)*(1.2914728_pr &
                  +(-1.6504171_pr+u)*u)*(0.012431285_pr+(-0.14092406_pr+u)*u)+0.14648438_pr*c4*(1.6816422_pr &
                  +(-2.5228331_pr+u)*u)*(1.1925378_pr+(-1.6182745_pr+u)*u)*(0.007307093_pr+(-0.10847576_pr &
                  +u)*u))))/(fpi4*LambdaX*mpi2)

             hrho1tau1=mevfm*(ga*(0.16889151_pr*cd*(-1.2705742_pr+u)*(0.21652155_pr+u)*(1.4119605_pr &
                  +(-2.1693922_pr+u)*u)*(0.89941343_pr+(-1.1919922_pr+u)*u)+ga*LambdaX*(1.0390625_pr*c1 &
                  *(-1.2278564_pr+u)*(0.017559313_pr+u)*(1.2287066_pr+(-1.9487414_pr+u)*u)*(0.66412923_pr &
                  +(-0.34133743_pr+u)*u)-0.19921875_pr*c3*(1.86082_pr+(-2.5829231_pr+u)*u)*(1.9047944_pr &
                  +(-1.441557_pr+u)*u)*(0.025980569_pr+(-0.2323826_pr+u)*u)+0.29296875_pr*c4*(1.6816422_pr &
                  +(-2.5228331_pr+u)*u)*(1.1925378_pr+(-1.6182745_pr+u)*u)*(0.007307093_pr+(-0.10847576_pr &
                  +u)*u))))/(fpi4*LambdaX*mpi2)

             hJ0dr0=0.0_pr

             hrho0DJ0=mevfm*(gA2*(0.87488139_pr*c1*(-1.5075509_pr+u)*(0.26375366_pr+u)*(1.6807042_pr &
                  +(-2.4330383_pr+u)*u)*(0.96407512_pr+(-1.3041604_pr+u)*u)+0.10945568_pr*c4*(1.6868327_pr &
                  +(-2.5235789_pr+u)*u)*(1.16087_pr+(-1.5798009_pr+u)*u)*(0.0058674843_pr+(-0.099090543_pr+u)*u) &
                  -0.65769729_pr*c3*(1.7722509_pr+(-2.5823827_pr+u)*u)*(1.1921659_pr+(-1.56672_pr+u)*u) &
                  *(0.0031379608_pr+(-0.079284979_pr+u)*u)))/(fpi4*mpi2)

             hJ1dr1=0.0_pr

             hrho1DJ1=mevfm*(gA2*(-0.43565926_pr*c1*(-1.2316322_pr+u)*(0.20486088_pr+u)*(1.3126963_pr &
                  +(-2.0776323_pr+u)*u)*(0.80436474_pr+(-1.0724984_pr+u)*u)-0.07297045_pr*c4*(1.6868327_pr &
                  +(-2.5235789_pr+u)*u)*(1.1608701_pr+(-1.5798009_pr+u)*u)*(0.0058674831_pr+(-0.099090535_pr &
                  +u)*u)+0.1942317_pr*c3*(1.7488022_pr+(-2.5635214_pr+u)*u)*(1.1781067_pr+(-1.5454099_pr+u)*u) &
                  *(0.0037519521_pr+(-0.082717466_pr+u)*u)))/(fpi4*mpi2)

             hJ0dr1=0.0_pr

             hrho1DJ0=mevfm*(gA2*(-0.49880616_pr*c1*(-1.2658085_pr+u)*(0.20574562_pr+u)*(1.4174942_pr &
                  +(-2.1861549_pr+u)*u)*(0.9423041_pr+(-1.257487_pr+u)*u)-0.036485225_pr*c4*(1.6868327_pr &
                  +(-2.5235789_pr+u)*u)*(1.1608701_pr+(-1.5798009_pr+u)*u)*(0.0058674831_pr+(-0.099090536_pr &
                  +u)*u)+0.10945568_pr*c3*(1.6868327_pr+(-2.5235789_pr+u)*u)*(1.1608701_pr+(-1.5798009_pr+u)*u) &
                  *(0.0058674829_pr+(-0.099090534_pr+u)*u)))/(fpi4*mpi2)

             hJ1dr0=0.0_pr

             hJ0J0=mevfm*(ga*(-0.1536163_pr*cd*(-1.3967122_pr+u)*(0.25293929_pr+u)*(1.6675202_pr+(-2.3021441_pr &
                  +u)*u)*(0.98775885_pr+(-1.1088032_pr+u)*u)+ga*LambdaX*(0.37890625_pr*c1*(-1.9915809_pr+u) &
                  *(-1.0930869_pr+u)*(-0.66138682_pr+u)*(0.62037829_pr+u)*(0.88993989_pr+(-1.6302515_pr+u)*u) &
                  -0.82470703_pr*c3*(-0.1151617_pr+u)*(0.19519195_pr+u)*(1.3633551_pr+(-2.2799728_pr+u)*u) &
                  *(0.92255094_pr+(-1.4901995_pr+u)*u)-0.31591797_pr*c4*(1.6836216_pr+(-2.5209696_pr+u)*u) &
                  *(1.1565142_pr+(-1.5751323_pr+u)*u)*(0.0053910501_pr+(-0.095764354_pr+u)*u)))) &
                  /(fpi4*LambdaX*mpi2)

             hJ0J1=mevfm*(ga*(0.10241087_pr*cd*(-1.3967122_pr+u)*(0.25293929_pr+u)*(1.6675202_pr+(-2.3021441_pr &
                  +u)*u)*(0.98775885_pr+(-1.1088032_pr+u)*u)+ga*LambdaX*(0.19000244_pr*c3*(-0.078125417_pr+u) &
                  *(0.032069683_pr+u)*(1.5685634_pr+(-2.4365485_pr+u)*u)*(1.0625904_pr+(-1.537248_pr+u)*u) &
                  -0.69067383_pr*c1*(-1.1098034_pr+u)*(0.18815441_pr+u)*(1.1047942_pr+(-1.9744345_pr+u)*u) &
                  *(0.7614037_pr+(-1.2677199_pr+u)*u)+0.21069336_pr*c4*(1.6833299_pr+(-2.5207351_pr+u)*u) &
                  *(1.1564102_pr+(-1.5750045_pr+u)*u)*(0.0053966658_pr+(-0.095801508_pr+u)*u)))) &
                  /(fpi4*LambdaX*mpi2)

             hJ1J1=mevfm*(ga*(0.051205435_pr*cd*(-1.3967122_pr+u)*(0.25293929_pr+u)*(1.6675202_pr+(-2.3021441_pr &
                  +u)*u)*(0.98775885_pr+(-1.1088032_pr+u)*u)+ga*LambdaX*(0.19219971_pr*c3*(-0.12213139_pr+u) &
                  *(0.2951145_pr+u)*(2.0521898_pr+(-2.8194623_pr+u)*u)*(1.2084447_pr+(-1.7226538_pr+u)*u) &
                  -0.54418945_pr*c1*(-1.7111589_pr+u)*(0.30206023_pr+u)*(1.421294_pr+(-2.2931824_pr+u)*u) &
                  *(0.85528206_pr+(-1.3180868_pr+u)*u)+0.10534668_pr*c4*(1.6833299_pr+(-2.5207351_pr+u)*u) &
                  *(1.1564102_pr+(-1.5750045_pr+u)*u)*(0.0053966658_pr+(-0.095801508_pr+u)*u)))) &
                  /(fpi4*LambdaX*mpi2)

             dhrho0rho0=mevfm*(-3.0696805e-6_pr*ce*(-0.63337387_pr + u)*(-0.6_pr + u)*(0.37421487_pr &
                  - 1.2218226_pr*u + u2)*(0.33405954_pr - 1.1537321_pr*u + u2) + gA*(0.27912917_pr*cd*(&
                  1.7611632_pr - 2.5821217_pr*u + u2)*(1.2266141_pr - 1.6507188_pr*u + u2)*(0.029143369_pr &
                  - 0.28917996_pr*u + u2) + gA*LambdaX*(-1491.7558_pr*c1*(0.67170782_pr - 1.62489_pr*u &
                  + u2)*(0.4525586_pr - 1.2543926_pr*u + u2)*(0.17617549_pr - 0.80139711_pr*u + u2) &
                  - 536.43232_pr*c3*(0.76316673_pr - 1.7218477_pr*u + u2)*(0.43361641_pr - 1.1835313_pr*u &
                  + u2)*(0.16367811_pr - 0.78054496_pr*u + u2) + 0.30876998_pr*c4*(-1.5784873_pr + u)*( &
                  -1.0937068_pr + u)*(1.340452_pr - 2.04474_pr*u + u2)*(0.019873806_pr - 0.25332432_pr*u &
                  + u2))))/(fpi4*LambdaX)

             dhrho1rho1=mevfm*(-3.0696805e-6_pr*ce*(-0.63337387_pr + u)*(-0.6_pr + u)*(0.37421487_pr &
                  - 1.2218226_pr*u + u2)*(0.33405954_pr - 1.1537321_pr*u + u2) + gA*(0.27912917_pr*cd*(&
                  1.7611632_pr - 2.5821217_pr*u + u2)*(1.2266141_pr - 1.6507188_pr*u + u2)*(0.029143369_pr &
                  - 0.28917996_pr*u + u2) + gA*LambdaX*(-1491.7558_pr*c1*(0.67170782_pr - 1.62489_pr*u &
                  + u2)*(0.4525586_pr - 1.2543926_pr*u + u2)*(0.17617549_pr - 0.80139711_pr*u + u2) &
                  - 536.43232_pr*c3*(0.76316673_pr - 1.7218477_pr*u + u2)*(0.43361641_pr - 1.1835313_pr*u &
                  + u2)*(0.16367811_pr - 0.78054496_pr*u + u2) + 0.30876998_pr*c4*(-1.5784873_pr + u)*( &
                  -1.0937068_pr + u)*(1.340452_pr - 2.04474_pr*u + u2)*(0.019873806_pr - 0.25332432_pr*u &
                  + u2))))/(fpi4*LambdaX)

             dhdr0dr0=0.0_pr
             dhdr1dr1=0.0_pr

             dhrho0Drho0=mevfm*(gA*(-0.52109887_pr*cd*(-1.3082043_pr + u)*(-0.071686644_pr + u)*(1.4677685_pr &
                  - 2.2766738_pr*u + u2)*(0.93721725_pr - 1.4019146_pr*u + u2) + gA*LambdaX*(-1.0505013_pr*c4*( &
                  -0.6756734_pr + u)*(-0.093758568_pr + u)*(1.59727_pr - 2.4899591_pr*u + u2)*(1.0589336_pr &
                  - 1.6598738_pr*u + u2) + 410.5668_pr*c3*(-0.82531999_pr + u)*(-0.33506221_pr + u)*(0.60822284_pr &
                  - 1.5136238_pr*u + u2)*(0.31121003_pr - 0.98997717_pr*u + u2) - 1232.118_pr*c1*(-0.8439295_pr &
                  + u)*(-0.39116927_pr + u)*(0.58509589_pr - 1.4743171_pr*u + u2)*(0.28266076_pr - 0.97625946_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             dhrho1Drho0=mevfm*(gA*(0.17369909_pr*cd*(-1.308206_pr + u)*(-0.071686588_pr + u)*(1.4677702_pr &
                  - 2.2766756_pr*u + u2)*(0.9372174_pr - 1.4019153_pr*u + u2) + gA*LambdaX*(0.0083633968_pr*c3*( &
                  -1.1273967_pr + u)*(-0.75552081_pr + u)*(-0.07564609_pr + u)*(36.068809_pr + u)*(1.033597_pr &
                  - 1.7618699_pr*u + u2) + 0.36997868_pr*c4*(-0.67567294_pr + u)*(-0.091577817_pr + u)*(1.4806242_pr &
                  - 2.3819186_pr*u + u2)*(1.075619_pr - 1.644617_pr*u + u2) - 0.68655766_pr*c1*(-0.063343379_pr &
                  + u)*(0.97369203_pr + u)*(1.3870134_pr - 2.3075684_pr*u + u2)*(1.0431426_pr - 1.6344504_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             dhrho1Drho1=mevfm*(gA*(0.34739925_pr*cd*(-1.3082043_pr + u)*(-0.071686644_pr + u)*(1.4677685_pr &
                  - 2.2766738_pr*u + u2)*(0.93721725_pr - 1.4019146_pr*u + u2) + gA*LambdaX*(0.7003342_pr*c4*( &
                  -0.6756734_pr + u)*(-0.093758568_pr + u)*(1.59727_pr - 2.4899591_pr*u + u2)*(1.0589336_pr &
                  - 1.6598738_pr*u + u2) - 352.88094_pr*c3*(-0.80527266_pr + u)*(-0.38439748_pr + u)*(0.56781414_pr &
                  - 1.468113_pr*u + u2)*(0.30457752_pr - 1.0189296_pr*u + u2) + 388.57352_pr*c1*(-0.83954784_pr &
                  + u)*(-0.4308896_pr + u)*(0.57229209_pr - 1.4557589_pr*u + u2)*(0.25241536_pr - 0.93218873_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             dhrho0tau0=mevfm*(ga*(1.0421949_pr*cd*(-1.3082058_pr+u)*(-0.071686595_pr+u)*(1.46777_pr &
                  +(-2.2766754_pr+u)*u)*(0.93721738_pr+(-1.4019152_pr+u)*u)+ga*LambdaX*(2.1171875_pr*c4 &
                  *(-0.6756739_pr+u)*(-0.090082162_pr+u)*(1.5215283_pr+(-2.4175155_pr+u)*u)*(1.0854982_pr &
                  +(-1.6625771_pr+u)*u)+224.0_pr*c1*(-0.94326159_pr+u)*(-0.29152041_pr+u)*(0.71677264_pr &
                  +(-1.5877363_pr+u)*u)*(0.25328282_pr+(-0.77748171_pr+u)*u)+48.0_pr*c3*(0.88506991_pr &
                  +(-1.8677526_pr+u)*u)*(0.70765726_pr+(-1.4865692_pr+u)*u)*(0.079022134_pr+(-0.37067813_pr &
                  +u)*u))))/(fpi4*LambdaX*mpi2)

             dhrho1tau0=mevfm*(ga*(-0.34739831_pr*cd*(-1.3082058_pr+u)*(-0.071686595_pr+u)*(1.46777_pr &
                  +(-2.2766754_pr+u)*u)*(0.93721738_pr+(-1.4019152_pr+u)*u)+ga*LambdaX*(0.3046875_pr*c3 &
                  *(-0.75536572_pr+u)*(-0.092756738_pr+u)*(2.1347117_pr+(-2.8953677_pr+u)*u)*(1.2457416_pr &
                  +(-1.8629201_pr+u)*u)-0.6953125_pr*c4*(-0.67567303_pr+u)*(-0.089866675_pr+u)*(1.5384947_pr &
                  +(-2.4324091_pr+u)*u)*(1.0870855_pr+(-1.6682029_pr+u)*u)+1.65625_pr*c1*(-0.058760606_pr+u) &
                  *(0.81179488_pr+u)*(1.3491153_pr+(-2.2757768_pr+u)*u)*(1.0158947_pr+(-1.6126348_pr+u)*u)))) &
                  /(fpi4*LambdaX*mpi2)

             dhrho1tau1=mevfm*(ga*(-0.69479662_pr*cd*(-1.3082058_pr+u)*(-0.071686595_pr+u)*(1.46777_pr &
                  +(-2.2766754_pr+u)*u)*(0.93721738_pr+(-1.4019152_pr+u)*u)+ga*LambdaX*(-1.390625_pr*c4 &
                  *(-0.67567303_pr+u)*(-0.089866675_pr+u)*(1.5384947_pr+(-2.4324091_pr+u)*u)*(1.0870855_pr &
                  +(-1.6682029_pr+u)*u)-12.0_pr*c3*(-0.2307474_pr+u)*(0.76898229_pr+(-1.7297015_pr+u)*u) &
                  *(0.4881705_pr+(-1.0916344_pr+u)*u)-192.0_pr*c1*(-0.89170462_pr+u)*(-0.40511316_pr+u) &
                  *(0.6397114_pr+(-1.5139956_pr+u)*u)*(0.21220628_pr+(-0.78918661_pr+u)*u)))) &
                  /(fpi4*LambdaX*mpi2)

             dhJ0dr0=0.0_pr

             dhrho0DJ0=mevfm*(gA2*(-0.59709561_pr*c4*(-0.69515033_pr+u)*(-0.089020409_pr+u)*(1.4150964_pr &
                  +(-2.3274531_pr+u)*u)*(1.0076893_pr+(-1.5882571_pr+u)*u)+4.125_pr*c3*(-0.74030576_pr+u) &
                  *(-0.084829197_pr+u)*(1.3684104_pr+(-2.2855723_pr+u)*u)*(0.95092906_pr+(-1.51202_pr+u)*u) &
                  -20.0_pr*c1*(-1.0733108_pr+u)*(-0.097980712_pr+u)*(0.99914461_pr+(-1.8457061_pr+u)*u) &
                  *(0.66726597_pr+(-1.0330025_pr+u)*u)))/(fpi4*mpi2)

             dhJ1dr1=0.0_pr

             dhrho1DJ1=mevfm*(gA2*(0.39806361_pr*c4*(-0.69515033_pr+u)*(-0.089020402_pr+u)*(1.4150967_pr &
                  +(-2.3274533_pr+u)*u)*(1.0076893_pr+(-1.5882572_pr+u)*u)- 1.125_pr*c3*(-0.74849274_pr+u) &
                  *(-0.080617806_pr+u)*(1.3642178_pr+(-2.2788157_pr+u)*u)*(0.96130355_pr+(-1.5129071_pr+u)*u) &
                  -5.0_pr*c1*(1.0009125_pr+(-1.9641805_pr+u)*u)*(0.6986326_pr+(-1.3789259_pr+u)*u)*(0.041493177_pr &
                  +u*(0.16810641_pr+u))))/(fpi4*mpi2)

             dhJ0dr1=0.0_pr

             dhrho1DJ0=mevfm*(gA2*(0.19903181_pr*c4*(-0.69515033_pr+u)*(-0.089020402_pr+u)*(1.4150967_pr &
                  +(-2.3274533_pr+u)*u)*(1.0076893_pr+(-1.5882572_pr+u)*u)-0.59709555_pr*c3*(-0.69515033_pr+u) &
                  *(-0.089020406_pr+u)*(1.4150964_pr+(-2.3274531_pr+u)*u)*(1.0076893_pr+(-1.5882571_pr+u)*u) &
                  +1.4395216_pr*c1*(-1.6910704_pr+u)*(-0.073466262_pr+u)*(1.5355865_pr+(-2.3975331_pr+u)*u) &
                  *(0.99872048_pr+(-1.5169799_pr+u)*u)))/(fpi4*mpi2)

             dhJ1dr0=0.0_pr

             dhJ0J0=mevfm*(ga*(0.45817098_pr*cd*(-1.6176653_pr+u)*(-0.066866022_pr+u)*(1.7318817_pr &
                  +(-2.5016073_pr+u)*u)*(1.0139312_pr+(-1.425551_pr+u)*u)+ga*LambdaX*(1.7695313_pr*c4 &
                  *(-0.69501265_pr+u)*(-0.088068636_pr+u)*(1.3968124_pr+(-2.3120281_pr+u)*u)*(0.99464383_pr &
                  +(-1.5741511_pr+u)*u)+28.0_pr*c3*(-0.52322774_pr+u)*(-0.14878162_pr+u)*(0.83449184_pr &
                  +(-1.7952477_pr+u)*u)*(0.5636391_pr+(-1.2398858_pr+u)*u)-48.0_pr*c1*(-0.86033423_pr+u) &
                  *(-0.098046912_pr+u)*(0.71197564_pr+(-1.6379394_pr+u)*u)*(0.50452901_pr+(-1.1286795_pr &
                  +u)*u))))/(fpi4*LambdaX*mpi2)

             dhJ0J1=mevfm*(ga*(-0.30544732_pr*cd*(-1.6176653_pr+u)*(-0.066866022_pr+u)*(1.7318817_pr &
                  +(-2.5016073_pr+u)*u)*(1.0139312_pr+(-1.425551_pr+u)*u)+ga*LambdaX*(-1.1796875_pr*c4 &
                  *(-0.69501214_pr+u)*(-0.088062144_pr+u)*(1.3949349_pr+(-2.3100216_pr+u)*u)*(0.99576969_pr &
                  +(-1.5747849_pr+u)*u)-1.5332031_pr*c3*(-0.63533554_pr+u)*(-0.083219793_pr+u)*(1.2047515_pr &
                  +(-2.1485114_pr+u)*u)*(0.85574473_pr+(-1.4759269_pr+u)*u)+4.515625_pr*c1*(-1.0974449_pr+u) &
                  *(-0.076646725_pr+u)*(1.0942329_pr+(-1.9814273_pr+u)*u)*(0.7814348_pr+(-1.3562458_pr &
                  +u)*u))))/(fpi4*LambdaX*mpi2)

             dhJ1J1=mevfm*(ga*(-0.15272366_pr*cd*(-1.6176653_pr+u)*(-0.066866022_pr+u)*(1.7318817_pr &
                  +(-2.5016073_pr+u)*u)*(1.0139312_pr+(-1.425551_pr+u)*u)+ga*LambdaX*(-0.58984375_pr*c4 &
                  *(-0.69501214_pr+u)*(-0.088062144_pr+u)*(1.3949349_pr+(-2.3100216_pr+u)*u)*(0.99576969_pr &
                  +(-1.5747849_pr+u)*u)-4.0_pr*c3*(-0.59708253_pr+u)*(-0.031303824_pr+u)*(1.029716_pr &
                  +(-1.979428_pr+u)*u)*(0.81581869_pr+(-1.3671856_pr+u)*u)-16.0_pr*c1*(0.93008679_pr &
                  +(-1.8984052_pr+u)*u)*(0.65743897_pr+(-1.3677632_pr+u)*u)*(0.048712336_pr+u*(0.22866833_pr &
                  +u)))))/(fpi4*LambdaX*mpi2)

             ddhrho0rho0=mevfm*(3.5082063e-6_pr*ce*(-4.2415409_pr + u)*(-0.63600549_pr + u)*(0.36430595_pr &
                  - 1.2064733_pr*u + u2)*(0.36384534_pr - 1.2034803_pr*u + u2) + gA*(-0.10526469_pr*cd*(-15.717486_pr &
                  + u)*(-1.1063208_pr + u)*(-0.79219406_pr + u)*(-0.15821661_pr + u)*(0.99686003_pr - 1.7359229_pr*u &
                  + u2) + gA*LambdaX*(-23.140197_pr*c4*(-0.54937402_pr + u)*(-0.25617028_pr + u)*(0.87042847_pr &
                  - 1.8342813_pr*u + u2)*(0.55052287_pr - 1.1703549_pr*u + u2) + 34417.155_pr*c3*(-0.75220388_pr &
                  + u)*(-0.45847411_pr + u)*(0.49790159_pr - 1.3890715_pr*u + u2)*(0.31420079_pr - 1.0830145_pr*u &
                  + u2) + 150844.82_pr*c1*(0.5097462_pr - 1.4229463_pr*u + u2)*(0.38532363_pr - 1.2176864_pr*u &
                  + u2)*(0.27148791_pr - 1.0367952_pr*u + u2))))/(fpi4*LambdaX)

             ddhrho1rho1=mevfm*(3.5082063e-6_pr*ce*(-4.2415409_pr + u)*(-0.63600549_pr + u)*(0.36430595_pr &
                  - 1.2064733_pr*u + u2)*(0.36384534_pr - 1.2034803_pr*u + u2) + gA*(-0.10526469_pr*cd*(-15.717486_pr &
                  + u)*(-1.1063208_pr + u)*(-0.79219406_pr + u)*(-0.15821661_pr + u)*(0.99686003_pr - 1.7359229_pr*u &
                  + u2) + gA*LambdaX*(-23.140197_pr*c4*(-0.54937402_pr + u)*(-0.25617028_pr + u)*(0.87042847_pr &
                  - 1.8342813_pr*u + u2)*(0.55052287_pr - 1.1703549_pr*u + u2) + 34417.155_pr*c3*(-0.75220388_pr &
                  + u)*(-0.45847411_pr + u)*(0.49790159_pr - 1.3890715_pr*u + u2)*(0.31420079_pr - 1.0830145_pr*u &
                  + u2) + 150844.82_pr*c1*(0.5097462_pr - 1.4229463_pr*u + u2)*(0.38532363_pr - 1.2176864_pr*u &
                  + u2)*(0.27148791_pr - 1.0367952_pr*u + u2))))/(fpi4*LambdaX)

             ddhrho0Drho0=mevfm*(gA*(-3.1120954_pr*cd*(-0.30605132_pr + u)*(1.2680842_pr + u)*(1.0861033_pr &
                  - 2.0506502_pr*u + u2)*(0.79950278_pr - 1.5067227_pr*u + u2) + gA*LambdaX*(-24514.431_pr*c3*( &
                  -0.76031522_pr + u)*(-0.46096893_pr + u)*(0.49710047_pr - 1.3844361_pr*u + u2)*(0.29795986_pr &
                  - 1.055072_pr*u + u2) + 22229.743_pr*c1*(-0.74615955_pr + u)*(-0.61188783_pr + u)*(0.46154959_pr &
                  - 1.3315221_pr*u + u2)*(0.25673397_pr - 0.99266282_pr*u + u2) + 110.78144_pr*c4*(-0.8714267_pr &
                  + u)*(-0.31340402_pr + u)*(0.6979966_pr - 1.6133393_pr*u + u2)*(0.37679963_pr - 0.97187797_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             ddhrho1Drho0=mevfm*(gA*(1.0374387_pr*cd*(-0.30605117_pr + u)*(1.2679976_pr + u)*(1.0860961_pr &
                  - 2.0506436_pr*u + u2)*(0.79949758_pr - 1.5067185_pr*u + u2) + gA*LambdaX*(7.4928249_pr*c1*( &
                  -0.37847001_pr + u)*(0.20959238_pr + u)*(1.3329648_pr - 2.2776684_pr*u + u2)*(0.96918292_pr &
                  - 1.6804594_pr*u + u2) + 0.50449085_pr*c4*(-0.25756953_pr + u)*(4.642167_pr + u)*(1.0022598_pr &
                  - 1.987984_pr*u + u2)*(0.89301348_pr - 1.6569412_pr*u + u2) - 8.2338865_pr*c3*(-0.24176677_pr &
                  + u)*(-0.021931235_pr + u)*(0.86956532_pr - 1.8470984_pr*u + u2)*(0.70044277_pr - 1.4791249_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             ddhrho1Drho1=mevfm*(gA*(2.0747303_pr*cd*(-0.30605132_pr + u)*(1.2680842_pr + u)*(1.0861033_pr &
                  - 2.0506502_pr*u + u2)*(0.79950278_pr - 1.5067227_pr*u + u2) + gA*LambdaX*(15284.793_pr*c3*( &
                  -0.74604143_pr + u)*(-0.48901435_pr + u)*(0.48038865_pr - 1.3647474_pr*u + u2)*(0.30523896_pr &
                  - 1.0777475_pr*u + u2) - 73.854292_pr*c4*(-0.8714267_pr + u)*(-0.31340402_pr + u)*(0.6979966_pr &
                  - 1.6133393_pr*u + u2)*(0.37679963_pr - 0.97187797_pr*u + u2) - 14523.513_pr*c1*(0.51615934_pr &
                  - 1.4339319_pr*u + u2)*(0.42639291_pr - 1.2743751_pr*u + u2)*(0.23481515_pr - 0.95144731_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             ddhrho0tau0=mevfm*(ga*(6.2246138_pr*cd*(-0.30605117_pr+u)*(1.2680013_pr+u)*(1.0860964_pr &
                  +(-2.0506438_pr+u)*u)*(0.79949779_pr+(-1.5067187_pr+u)*u)+ga*LambdaX*(7.25_pr*c4*(-0.25490416_pr &
                  +u)*(1.9314619_pr+u)*(0.97416501_pr+(-1.9598149_pr+u)*u)*(0.85402263_pr+(-1.6270877_pr+u)*u) &
                  -4096.0_pr*c1*(-0.82877841_pr+u)*(-0.61205527_pr+u)*(0.55294169_pr+(-1.4322212_pr+u)*u) &
                  *(0.21133499_pr+(-0.85194512_pr+u)*u)-1024.0_pr*c3*(-1.5110662_pr+u)*(-0.81113397_pr+u) &
                  *(0.48220899_pr+(-1.3083223_pr+u)*u)*(0.13559828_pr+(-0.71947752_pr+u)*u))))/(fpi4*LambdaX*mpi2)

             ddhrho1tau0=mevfm*(ga*(-2.0748713_pr*cd*(-0.30605117_pr+u)*(1.2680013_pr+u)*(1.0860964_pr &
                  +(-2.0506438_pr+u)*u)*(0.79949779_pr+(-1.5067187_pr+u)*u)+ga*LambdaX*(-2.1875_pr*c4 &
                  *(-0.25526539_pr+u)*(2.1267165_pr+u)*(0.98082879_pr+(-1.9667092_pr+u)*u)*(0.86027632_pr &
                  +(-1.6333133_pr+u)*u)+4.0_pr*c3*(-0.27578203_pr+u)*(0.73913647_pr+u)*(0.98670704_pr &
                  +(-1.9658775_pr+u)*u)*(0.85430872_pr+(-1.6131019_pr+u)*u)-29.0_pr*c1*(-0.37708308_pr+u) &
                  *(0.063860656_pr+u)*(1.1337991_pr+(-2.0962155_pr+u)*u)*(0.84554078_pr+(-1.5461655_pr+u)*u)))) &
                  /(fpi4*LambdaX*mpi2)

             ddhrho1tau1=mevfm*(ga*(-4.1497425_pr*cd*(-0.30605117_pr+u)*(1.2680013_pr+u)*(1.0860964_pr &
                  +(-2.0506438_pr+u)*u)*(0.79949779_pr+(-1.5067187_pr+u)*u)+ga*LambdaX*(-4.375_pr*c4 &
                  *(-0.25526539_pr+u)*(2.1267165_pr+u)*(0.98082879_pr+(-1.9667092_pr+u)*u)*(0.86027632_pr &
                  +(-1.6333133_pr+u)*u)-9728.0_pr*c1*(-0.72811509_pr+u)*(-0.40683586_pr+u)*(0.51101775_pr &
                  +(-1.4100026_pr+u)*u)*(0.32463579_pr+(-1.0747833_pr+u)*u)+1536.0_pr*c3*(0.61725539_pr &
                  +(-1.5571261_pr+u)*u)*(0.39877632_pr+(-1.1610973_pr+u)*u)*(0.15498418_pr+(-0.77760994_pr &
                  +u)*u))))/(fpi4*LambdaX*mpi2)

          End If !! if (u.gt.ucut3n..
       End If !! if(use_dme3N_terms
       !
    End If !! if (dmeorder.ge.2....
#endif
    !
    Urhorho=0.0_pr   ; Urhotau=0.0_pr
    UrhoDrho=0.0_pr  ; Unablarho=0.0_pr
    UJnablarho=0.0_pr; UrhonablaJ=0.0_pr
    Urhorhopr=0.0_pr ; UJJ=0.0_pr
    !
    ! Notations for Uamplitudes(0:3,0:7)
    ! t for Uamplitudes(t,*)
    ! 0 -> 0,0
    ! 1 -> 1,1
    ! 2 -> 0,1
    ! 3 -> 1,0
    ! n for Uamplitudes(*,n)
    ! 0 -> U
    ! 1 -> dU/dRHO_0
    ! 2 -> dU/dRHO_1
    ! 3 -> d2U/(dRHO_0*dRHO_0)
    ! 4 -> d2U/(dRHO_1*dRHO_1) 
    ! 5 -> d2U/(dRHO_0*dRHO_1) 
    ! 6 -> dU/d(TAU_0)   
    ! 7 -> dU/d(Delta RHO_0)   
    !
    !! 2N terms
    Do t=0,1
       ph=1.0_pr
       If(t.Eq.1) ph=-1.0_pr
       Urhorho(t,0)=Crho(t)+Cdrho(t)*rho(0)**sigma &
            +0.50_pr*(arhorho+ph*brhorho)*mevfm
       Urhotau(t,0)=Ctau(t)+0.50_pr*(arhotau+ph*brhotau)*mevfm
       UrhoDrho(t,0)=Crdr(t)+ac2*0.50_pr*(arhoDrho+ph*brhoDrho)*mevfm
       UJJ(t,0)=CJ(t)+0.50_pr*(ajj+ph*bjj)*mevfm
       Unablarho(t,0)=Cnrho(t)+0.50_pr*(adrdr+ph*bdrdr)*mevfm
       UrhonablaJ(t,0)=Crdj(t)
       UJnablarho(t,0)=Cjdr(t)
       !!       Urhorhopr(t,0)=CpV0(t)+CpV1(t)*rho(0)

       Urhorho(t,1)=sigma*Cdrho(t)*(rho(0)**sigma)/(rho(0)+eps) &
            +0.50_pr*(darhorho+ph*dbrhorho)*du*mevfm
       Urhotau(t,1)=0.50_pr*(darhotau+ph*dbrhotau)*du*mevfm
       UrhoDrho(t,1)=ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*du*mevfm
       UJJ(t,1)=0.50_pr*(dajj+ph*dbjj)*du*mevfm
       Unablarho(t,1)=0.50_pr*(dadrdr+ph*dbdrdr)*du*mevfm

       Urhorho(t,6)=0.50_pr*(darhorho+ph*dbrhorho)*dtu*mevfm
       Urhotau(t,6)=0.50_pr*(darhotau+ph*dbrhotau)*dtu*mevfm
       UrhoDrho(t,6)=ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*dtu*mevfm
       UJJ(t,6)=0.50_pr*(dajj+ph*dbjj)*dtu*mevfm
       Unablarho(t,6)=0.50_pr*(dadrdr+ph*dbdrdr)*dtu*mevfm

       Urhorho(t,7)=0.50_pr*(darhorho+ph*dbrhorho)*dlu*mevfm
       Urhotau(t,7)=0.50_pr*(darhotau+ph*dbrhotau)*dlu*mevfm
       UrhoDrho(t,7)=ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*dlu*mevfm
       UJJ(t,7)=0.50_pr*(dajj+ph*dbjj)*dlu*mevfm
       Unablarho(t,7)=0.50_pr*(dadrdr+ph*dbdrdr)*dlu*mevfm


       Urhorho(t,3)=sigma*(sigma-1.0_pr)*Cdrho(t)*(rho(0)**sigma)/(rho(0)**2+eps) &
            +0.50_pr*(darhorho+ph*dbrhorho)*ddu*mevfm &
            +0.50_pr*(ddarhorho+ph*ddbrhorho)*du*du*mevfm
       Urhotau(t,3)=0.50_pr*(darhotau+ph*dbrhotau)*ddu*mevfm &
            +0.50_pr*(ddarhotau+ph*ddbrhotau)*du*du*mevfm
       UrhoDrho(t,3)=ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*ddu*mevfm &
            +ac2*0.50_pr*(ddarhoDrho+ph*ddbrhoDrho)*du*du*mevfm

       !!       Urhorhopr(t,1)=CpV1(t)
    End Do
    Urhorhopr(0,0) = CpV0(0)*(1.0_pr-CpV1(0)*rho(0)/0.16_pr)          &
         +CpV0(1)*(1.0_pr-CpV1(1)*rho(0)/0.16_pr)
    Urhorhopr(1,0) = CpV0(0)*(1.0_pr-CpV1(0)*rho(0)/0.16_pr)          &
         +CpV0(1)*(1.0_pr-CpV1(1)*rho(0)/0.16_pr)
    Urhorhopr(2,0) = (CpV0(0)*(1.0_pr-CpV1(0)*rho(0)/0.16_pr)         &
         -CpV0(1)*(1.0_pr-CpV1(1)*rho(0)/0.16_pr))*2.0_pr
    Urhorhopr(0,1) = (-CpV0(0)*CpV1(0)-CpV0(1)*CpV1(1))/0.16_pr
    Urhorhopr(1,1) = (-CpV0(0)*CpV1(0)-CpV0(1)*CpV1(1))/0.16_pr
    Urhorhopr(2,1) = 2.0_pr*(-CpV0(0)*CpV1(0)+CpV0(1)*CpV1(1))/0.16_pr
    Urhorhopr=Urhorhopr/16.0_pr
#ifndef hide_dme    
    !! 3N terms
    If (use_DME3N_terms) Then 
       Urhorho(0,0)=Urhorho(0,0)+hrho0rho0*rho(0)
       Urhorho(1,0)=Urhorho(1,0)+hrho1rho1*rho(0)
       Unablarho(0,0)=Unablarho(0,0)+hdr0dr0*rho(0)
       Unablarho(1,0)=Unablarho(1,0)+hdr1dr1*rho(0)
       UrhoDrho(0,0)=UrhoDrho(0,0)+hrho0Drho0*rho(0)*ac3
       UrhoDrho(1,0)=UrhoDrho(1,0)+hrho1Drho1*rho(0)*ac3
       UrhoDrho(3,0)=UrhoDrho(3,0)+hrho1Drho0*rho(1)*ac3
       Urhotau(0,0)=Urhotau(0,0)+hrho0tau0*rho(0)
       Urhotau(1,0)=Urhotau(1,0)+hrho1tau1*rho(0)
       Urhotau(3,0)=Urhotau(3,0)+hrho1tau0*rho(1)
       UJnablarho(0,0)=UJnablarho(0,0)+hJ0dr0*rho(0)
       UJnablarho(1,0)=UJnablarho(1,0)+hJ1dr1*rho(0)
       UJnablarho(2,0)=UJnablarho(2,0)+hJ0dr1*rho(1)
       UJnablarho(3,0)=UJnablarho(3,0)+hJ1dr0*rho(1)
       UrhonablaJ(0,0)=UrhonablaJ(0,0)+hrho0dJ0*rho(0)
       UrhonablaJ(1,0)=UrhonablaJ(1,0)+hrho1dJ1*rho(0)
       UrhonablaJ(3,0)=UrhonablaJ(3,0)+hrho1dJ0*rho(1)
       UJJ(0,0)=UJJ(0,0)+hJ0J0*rho(0)
       UJJ(1,0)=UJJ(1,0)+hJ1J1*rho(0)
       UJJ(2,0)=UJJ(2,0)+hJ0J1*rho(1)

       Urhorho(0,1)=Urhorho(0,1)+dhrho0rho0*du*rho(0) +hrho0rho0
       Urhorho(1,1)=Urhorho(1,1)+dhrho1rho1*du*rho(0) +hrho1rho1
       Unablarho(0,1)=Unablarho(0,1)+dhdr0dr0*du*rho(0) +hdr0dr0
       Unablarho(1,1)=Unablarho(1,1)+dhdr1dr1*du*rho(0) +hdr1dr1
       UrhoDrho(0,1)=UrhoDrho(0,1)+(dhrho0Drho0*du*rho(0) +hrho0Drho0)*ac3
       UrhoDrho(1,1)=UrhoDrho(1,1)+(dhrho1Drho1*du*rho(0) +hrho1Drho1)*ac3
       UrhoDrho(3,1)=UrhoDrho(3,1)+dhrho1Drho0*du*rho(1)*ac3
       Urhotau(0,1)=Urhotau(0,1)+dhrho0tau0*du*rho(0) +hrho0tau0
       Urhotau(1,1)=Urhotau(1,1)+dhrho1tau1*du*rho(0) +hrho1tau1
       Urhotau(3,1)=Urhotau(3,1)+dhrho1tau0*du*rho(1)
       UJnablarho(0,1)=UJnablarho(0,1)+dhJ0dr0*du*rho(0) +hJ0dr0
       UJnablarho(1,1)=UJnablarho(1,1)+dhJ1dr1*du*rho(0) +hJ1dr1
       UJnablarho(2,1)=UJnablarho(2,1)+dhJ0dr1*du*rho(1)
       UJnablarho(3,1)=UJnablarho(3,1)+dhJ1dr0*du*rho(1)
       UrhonablaJ(0,1)=UrhonablaJ(0,1)+dhrho0dJ0*du*rho(0) +hrho0dJ0
       UrhonablaJ(1,1)=UrhonablaJ(1,1)+dhrho1dJ1*du*rho(0) +hrho1dJ1
       UrhonablaJ(3,1)=UrhonablaJ(3,1)+dhrho1dJ0*du*rho(1)
       UJJ(0,1)=UJJ(0,1)+dhJ0J0*du*rho(0) +hJ0J0
       UJJ(1,1)=UJJ(1,1)+dhJ1J1*du*rho(0) +hJ1J1
       UJJ(2,1)=UJJ(2,1)+dhJ0J1*du*rho(1)

       Urhorho(0,6)=Urhorho(0,6)+dhrho0rho0*dtu*rho(0)
       Urhorho(1,6)=Urhorho(1,6)+dhrho1rho1*dtu*rho(0)
       Unablarho(0,6)=Unablarho(0,6)+dhdr0dr0*dtu*rho(0)
       Unablarho(1,6)=Unablarho(1,6)+dhdr1dr1*dtu*rho(0)
       UrhoDrho(0,6)=UrhoDrho(0,6)+dhrho0Drho0*dtu*rho(0)*ac3
       UrhoDrho(1,6)=UrhoDrho(1,6)+dhrho1Drho1*dtu*rho(0)*ac3
       UrhoDrho(3,6)=UrhoDrho(3,6)+dhrho1Drho0*dtu*rho(1)*ac3
       Urhotau(0,6)=Urhotau(0,6)+dhrho0tau0*dtu*rho(0)
       Urhotau(1,6)=Urhotau(1,6)+dhrho1tau1*dtu*rho(0)
       Urhotau(3,6)=Urhotau(3,6)+dhrho1tau0*dtu*rho(1)
       UJnablarho(0,6)=UJnablarho(0,6)+dhJ0dr0*dtu*rho(0)
       UJnablarho(1,6)=UJnablarho(1,6)+dhJ1dr1*dtu*rho(0)
       UJnablarho(2,6)=UJnablarho(2,6)+dhJ0dr1*dtu*rho(1)
       UJnablarho(3,6)=UJnablarho(3,6)+dhJ1dr0*dtu*rho(1)
       UrhonablaJ(0,6)=UrhonablaJ(0,6)+dhrho0dJ0*dtu*rho(0)
       UrhonablaJ(1,6)=UrhonablaJ(1,6)+dhrho1dJ1*dtu*rho(0)
       UrhonablaJ(3,6)=UrhonablaJ(3,6)+dhrho1dJ0*dtu*rho(1)
       UJJ(0,6)=UJJ(0,6)+dhJ0J0*dtu*rho(0)
       UJJ(1,6)=UJJ(1,6)+dhJ1J1*dtu*rho(0)
       UJJ(2,6)=UJJ(2,6)+dhJ0J1*dtu*rho(1)

       Urhorho(0,7)=Urhorho(0,7)+dhrho0rho0*dlu*rho(0)
       Urhorho(1,7)=Urhorho(1,7)+dhrho1rho1*dlu*rho(0)
       Unablarho(0,7)=Unablarho(0,7)+dhdr0dr0*dlu*rho(0)
       Unablarho(1,7)=Unablarho(1,7)+dhdr1dr1*dlu*rho(0)
       UrhoDrho(0,7)=UrhoDrho(0,7)+dhrho0Drho0*dlu*rho(0)*ac3
       UrhoDrho(1,7)=UrhoDrho(1,7)+dhrho1Drho1*dlu*rho(0)*ac3
       UrhoDrho(3,7)=UrhoDrho(3,7)+dhrho1Drho0*dlu*rho(1)*ac3
       Urhotau(0,7)=Urhotau(0,7)+dhrho0tau0*dlu*rho(0)
       Urhotau(1,7)=Urhotau(1,7)+dhrho1tau1*dlu*rho(0)
       Urhotau(3,7)=Urhotau(3,7)+dhrho1tau0*dlu*rho(1)
       UJnablarho(0,7)=UJnablarho(0,7)+dhJ0dr0*dlu*rho(0)
       UJnablarho(1,7)=UJnablarho(1,7)+dhJ1dr1*dlu*rho(0)
       UJnablarho(2,7)=UJnablarho(2,7)+dhJ0dr1*dlu*rho(1)
       UJnablarho(3,7)=UJnablarho(3,7)+dhJ1dr0*dlu*rho(1)
       UrhonablaJ(0,7)=UrhonablaJ(0,7)+dhrho0dJ0*dlu*rho(0)
       UrhonablaJ(1,7)=UrhonablaJ(1,7)+dhrho1dJ1*dlu*rho(0)
       UrhonablaJ(3,7)=UrhonablaJ(3,7)+dhrho1dJ0*dlu*rho(1)
       UJJ(0,7)=UJJ(0,7)+dhJ0J0*dlu*rho(0)
       UJJ(1,7)=UJJ(1,7)+dhJ1J1*dlu*rho(0)
       UJJ(2,7)=UJJ(2,7)+dhJ0J1*dlu*rho(1)


       Urhorho(0,3)=Urhorho(0,3)+2.0_pr*dhrho0rho0*du &
            +ddhrho0rho0*du*du*rho(0) +dhrho0rho0*ddu*rho(0)
       Urhorho(1,3)=Urhorho(1,3)+2.0_pr*dhrho1rho1*du &
            +ddhrho1rho1*du*du*rho(0) +dhrho1rho1*ddu*rho(0)
       UrhoDrho(0,3)=UrhoDrho(0,3)+(2.0_pr*dhrho0Drho0*du &
            +ddhrho0Drho0*du*du*rho(0) +dhrho0Drho0*ddu*rho(0))*ac3
       UrhoDrho(1,3)=UrhoDrho(1,3)+(2.0_pr*dhrho1Drho1*du &
            +ddhrho1Drho1*du*du*rho(0)+dhrho1Drho1*ddu*rho(0))*ac3
       Urhotau(0,3)=Urhotau(0,3)+2.0_pr*dhrho0tau0*du &
            +ddhrho0tau0*du*du*rho(0) +dhrho0tau0*ddu*rho(0)
       Urhotau(1,3)=Urhotau(1,3)+2.0_pr*dhrho1tau1*du &
            +ddhrho1tau1*du*du*rho(0) +dhrho1tau1*ddu*rho(0)
       UrhoDrho(3,3)=UrhoDrho(3,3) +(ddhrho1Drho0*du*du*rho(1) &
            +dhrho1Drho0*ddu*rho(1))*ac3
       Urhotau(3,3)=Urhotau(3,3) +ddhrho1tau0*du*du*rho(1) &
            +dhrho1tau0*ddu*rho(1) 


       UrhoDrho(3,2)=UrhoDrho(3,2)+hrho1Drho0*ac3
       Urhotau(3,2)=Urhotau(3,2)+hrho1tau0
       UJnablarho(3,2)=UJnablarho(3,2)+hJ1dr0
       UJnablarho(2,2)=UJnablarho(2,2)+hJ0dr1
       UrhonablaJ(3,2)=UrhonablaJ(3,2)+hrho1dJ0
       UJJ(2,2)=UJJ(2,2)+hJ0J1

       UrhoDrho(3,5)=UrhoDrho(3,5) +dhrho1Drho0*du*ac3
       Urhotau(3,5)=Urhotau(3,5) +dhrho1tau0*du

    End If   !! if(use3Nterm..

    UEnonstdr=0.0_pr; UFnonstdr=0.0_pr; URnonstdr=0.0_pr    
#endif   
    !
    If (.Not.use_j2terms) Then
       UJJ=0.0_pr
    End If
    !
  End Subroutine calculate_U_parameters
  !==================================================================================
  !
#ifndef hide_dme    
  !==================================================================================
  Subroutine load_tables
    !==================================================================================
    Implicit None
    !
    !  Urho0rho0
    !
    ctr0r0=0.0_pr
    ctr0r0(1,1,1) = (-0.964246762695464_pr*c3*gA2)/fpi4 - &
         (0.0263462308309038_pr*c4*gA2)/fpi4 - (0.1875_pr*ce)/(fpi4*LambdaX) + &
         (0.046875_pr*cd*ga)/(fpi4*LambdaX)
    ctr0r0(1,1,3) = (-6.84195361981036_pr*c1*gA2)/fpi4 + &
         (5.01410456224811_pr*c3*gA2)/fpi4 + (0.0732776804846939_pr*c4*gA2)/fpi4 &
         - (0.06328125_pr*cd*ga)/(fpi4*LambdaX)
    ctr0r0(2,1,3) = (0.147698198400372_pr*c3*gA2)/fpi4 + &
         (0.0703526785714286_pr*c4*gA2)/fpi4 - &
         (0.02109375_pr*cd*ga)/(fpi4*LambdaX)
    ctr0r0(1,2,4) = (-0.1875_pr*c3*gA2)/fpi4 - (0.375_pr*c4*gA2)/fpi4 + &
         (0.140625_pr*cd*ga)/(fpi4*LambdaX)
    ctr0r0(1,1,5) = (-23.1434396482781_pr*c1*gA2)/fpi4 + &
         (12.2465870768616_pr*c3*gA2)/fpi4 - (0.18795405731824_pr*c4*gA2)/fpi4 + &
         (0.017578125_pr*cd*ga)/(fpi4*LambdaX)
    ctr0r0(2,1,5) = (35.0618960855328_pr*c1*gA2)/fpi4 - &
         (16.8781682906031_pr*c3*gA2)/fpi4 + (0.0543163105867347_pr*c4*gA2)/fpi4 &
         - (0.06328125_pr*cd*ga)/(fpi4*LambdaX)
    ctr0r0(1,2,6) = (-325.455281192602_pr*c1*gA2)/fpi4 + &
         (83.7334611805583_pr*c3*gA2)/fpi4 + (1.3574300625_pr*c4*gA2)/fpi4
    ctr0r0(2,2,6) = (32.9911875098339_pr*c3*gA2)/fpi4 + &
         (0.16875_pr*c4*gA2)/fpi4
    ctr0r0(1,1,7) = (-5761.43588948783_pr*c1*gA2)/fpi4 - &
         (869.283572934387_pr*c3*gA2)/fpi4 - (13.4068326170281_pr*c4*gA2)/fpi4
    ctr0r0(2,1,7) = (2192.10779645869_pr*c1*gA2)/fpi4 + &
         (263.842952909344_pr*c3*gA2)/fpi4 + (3.99388560267857_pr*c4*gA2)/fpi4 - &
         (0.00439453125_pr*cd*ga)/(fpi4*LambdaX)
    ctr0r0(3,1,7) = (-198.153843790565_pr*c1*gA2)/fpi4 - &
         (58.742085950407_pr*c3*gA2)/fpi4 - (0.4904296875_pr*c4*gA2)/fpi4
    ctr0r0(1,3,7) = (-526.213527565944_pr*c3*gA2)/fpi4 - &
         (2.84068526785714_pr*c4*gA2)/fpi4
    ctr0r0(1,2,8) = (11194.2459766741_pr*c1*gA2)/fpi4 - &
         (901.889210563564_pr*c3*gA2)/fpi4 + (1.90516445471939_pr*c4*gA2)/fpi4
    ctr0r0(2,2,8) = (-751.043571428571_pr*c1*gA2)/fpi4 + &
         (1336.68785503856_pr*c3*gA2)/fpi4 + (4.15185267857143_pr*c4*gA2)/fpi4
    ctr0r0(1,1,9) = (-19577.5085802329_pr*c1*gA2)/fpi4 - &
         (20061.3662796935_pr*c3*gA2)/fpi4 - (402.203416845285_pr*c4*gA2)/fpi4
    ctr0r0(2,1,9) = (6927.60941249884_pr*c1*gA2)/fpi4 + &
         (10204.912096572_pr*c3*gA2)/fpi4 + (157.45919188058_pr*c4*gA2)/fpi4
    ctr0r0(3,1,9) = (-1688.44574624066_pr*c1*gA2)/fpi4 - &
         (1848.73711933536_pr*c3*gA2)/fpi4 - (17.251846749442_pr*c4*gA2)/fpi4
    ctr0r0(1,3,9) = (-14639.9109308036_pr*c1*gA2)/fpi4 - &
         (2073.4495546875_pr*c3*gA2)/fpi4 + (12.406921875_pr*c4*gA2)/fpi4
    ctr0r0(1,2,10) = (-6051.26397950415_pr*c1*gA2)/fpi4 + &
         (25318.5711035127_pr*c3*gA2)/fpi4 + (1111.77481741071_pr*c4*gA2)/fpi4
    ctr0r0(2,2,10) = (20263.462438058_pr*c1*gA2)/fpi4 - &
         (645.573171294818_pr*c3*gA2)/fpi4 - (220.10902734375_pr*c4*gA2)/fpi4
    ctr0r0(1,1,11) = (-117612.589511906_pr*c1*gA2)/fpi4 - &
         (22681.6420225835_pr*c3*gA2)/fpi4 + (1494.47871529855_pr*c4*gA2)/fpi4
    ctr0r0(2,1,11) = (78131.8227221231_pr*c1*gA2)/fpi4 + &
         (11306.3195229362_pr*c3*gA2)/fpi4 - (1044.42429449456_pr*c4*gA2)/fpi4
    ctr0r0(3,1,11) = (-19244.1407277278_pr*c1*gA2)/fpi4 - &
         (4031.58809797695_pr*c3*gA2)/fpi4 + (143.952209472656_pr*c4*gA2)/fpi4
    ctr0r0(1,3,11) = (-12782.365790625_pr*c1*gA2)/fpi4 - &
         (19898.5801884375_pr*c3*gA2)/fpi4 - (763.061715_pr*c4*gA2)/fpi4
    ctr0r0(1,2,12) = (112532.892207631_pr*c1*gA2)/fpi4 + &
         (2673.93101274094_pr*c3*gA2)/fpi4 - (3075.39951940848_pr*c4*gA2)/fpi4
    ctr0r0(2,2,12) = (-4174.98147321429_pr*c1*gA2)/fpi4 + &
         (16959.475747232_pr*c3*gA2)/fpi4 + (1249.48156138393_pr*c4*gA2)/fpi4
    ctr0r0(1,1,13) = (-62220.7516704385_pr*c1*gA2)/fpi4 - &
         (37645.4468013683_pr*c3*gA2)/fpi4 - (1458.88122707611_pr*c4*gA2)/fpi4
    ctr0r0(2,1,13) = (49392.5657264478_pr*c1*gA2)/fpi4 + &
         (35938.6848928339_pr*c3*gA2)/fpi4 + (1614.83538151428_pr*c4*gA2)/fpi4
    ctr0r0(3,1,13) = (-21590.8338635122_pr*c1*gA2)/fpi4 - &
         (12971.4689072793_pr*c3*gA2)/fpi4 - (402.76586151123_pr*c4*gA2)/fpi4
    ctr0r0(1,3,13) = (-66293.2611893973_pr*c1*gA2)/fpi4 - &
         (8811.70601861049_pr*c3*gA2)/fpi4 + (1385.87636601563_pr*c4*gA2)/fpi4
    ctr0r0(1,2,14) = (-8634.9144059954_pr*c1*gA2)/fpi4 + &
         (13584.0523217795_pr*c3*gA2)/fpi4 + (1941.14579905622_pr*c4*gA2)/fpi4
    ctr0r0(2,2,14) = (46092.3086405378_pr*c1*gA2)/fpi4 + &
         (6421.96926354563_pr*c3*gA2)/fpi4 - (1169.36919344657_pr*c4*gA2)/fpi4
    ctr0r0(1,1,15) = (-54823.1734859566_pr*c1*gA2)/fpi4 - &
         (21478.1701524643_pr*c3*gA2)/fpi4 + (30.6758968031093_pr*c4*gA2)/fpi4
    ctr0r0(2,1,15) = (79796.8451847046_pr*c1*gA2)/fpi4 + &
         (29550.3190922257_pr*c3*gA2)/fpi4 - (265.149560200718_pr*c4*gA2)/fpi4
    ctr0r0(3,1,15) = (-34549.4910141965_pr*c1*gA2)/fpi4 - &
         (12487.0361763301_pr*c3*gA2)/fpi4 + (129.866342965262_pr*c4*gA2)/fpi4
    ctr0r0(1,3,15) = (-6841.585265625_pr*c1*gA2)/fpi4 - &
         (8394.50421527972_pr*c3*gA2)/fpi4 - (637.487162214007_pr*c4*gA2)/fpi4
    ctr0r0(1,2,16) = (-966.144667675781_pr*c1*gA2)/fpi4 - &
         (1635.80437776282_pr*c3*gA2)/fpi4 - (20.2837476928711_pr*c4*gA2)/fpi4
    ctr0r0(2,2,16) = (10019.4643315377_pr*c1*gA2)/fpi4 + &
         (5806.54125301064_pr*c3*gA2)/fpi4 + (162.342586296387_pr*c4*gA2)/fpi4
    ctr0r0(1,1,17) = (-27328.2660786696_pr*c1*gA2)/fpi4 - &
         (11536.8965413_pr*c3*gA2)/fpi4 - (0.153939957765852_pr*c4*gA2)/fpi4
    ctr0r0(2,1,17) = (51084.8335772702_pr*c1*gA2)/fpi4 + &
         (21317.5006264206_pr*c3*gA2)/fpi4 - (5.04213459735325_pr*c4*gA2)/fpi4
    ctr0r0(3,1,17) = (-23176.714608858_pr*c1*gA2)/fpi4 - &
         (9591.77377013768_pr*c3*gA2)/fpi4 - (7.01374572372437_pr*c4*gA2)/fpi4
    ctr0r0(1,3,17) = (-839.926762520926_pr*c1*gA2)/fpi4
    ctr0r0(1,2,18) = (-71.6191333138602_pr*c1*gA2)/fpi4 - &
         (29.5562121168625_pr*c3*gA2)/fpi4
    ctr0r0(2,2,18) = (1117.2724960968_pr*c1*gA2)/fpi4 + &
         (423.729200499137_pr*c3*gA2)/fpi4 + (5.07093692321777_pr*c4*gA2)/fpi4
    ctr0r0(1,1,19) = (-9253.22855094582_pr*c1*gA2)/fpi4 - &
         (4005.49423137719_pr*c3*gA2)/fpi4
    ctr0r0(2,1,19) = (21870.1711895189_pr*c1*gA2)/fpi4 + &
         (9310.62118848288_pr*c3*gA2)/fpi4 + (0.0769699788829259_pr*c4*gA2)/fpi4
    ctr0r0(3,1,19) = (-11369.3247545611_pr*c1*gA2)/fpi4 - &
         (4767.97117814042_pr*c3*gA2)/fpi4 - (0.656709900856018_pr*c4*gA2)/fpi4
    ctr0r0(2,2,20) = (17.9047833284651_pr*c1*gA2)/fpi4 + &
         (7.38905302921563_pr*c3*gA2)/fpi4
    ctr0r0(1,1,21) = (-2139.34177995088_pr*c1*gA2)/fpi4 - &
         (951.965547431256_pr*c3*gA2)/fpi4
    ctr0r0(2,1,21) = (6585.92097742191_pr*c1*gA2)/fpi4 + &
         (2872.41393631472_pr*c3*gA2)/fpi4
    ctr0r0(3,1,21) = (-4206.97567815619_pr*c1*gA2)/fpi4 - &
         (1804.69962149452_pr*c3*gA2)/fpi4 - &
         (0.00962124736036573_pr*c4*gA2)/fpi4
    ctr0r0(1,1,23) = (-320.600347402137_pr*c1*gA2)/fpi4 - &
         (147.828203709359_pr*c3*gA2)/fpi4
    ctr0r0(2,1,23) = (1373.03667928303_pr*c1*gA2)/fpi4 + &
         (615.472966679238_pr*c3*gA2)/fpi4
    ctr0r0(3,1,23) = (-1139.80667854591_pr*c1*gA2)/fpi4 - &
         (500.614173373695_pr*c3*gA2)/fpi4
    ctr0r0(1,1,25) = (-28.0192319249075_pr*c1*gA2)/fpi4 - &
         (13.7137479594838_pr*c3*gA2)/fpi4
    ctr0r0(2,1,25) = (187.596940698279_pr*c1*gA2)/fpi4 + &
         (87.2234077809484_pr*c3*gA2)/fpi4
    ctr0r0(3,1,25) = (-216.193884091206_pr*c1*gA2)/fpi4 - &
         (97.5972978305762_pr*c3*gA2)/fpi4
    ctr0r0(1,1,27) = (-1.08369739154545_pr*c1*gA2)/fpi4 - &
         (0.619027285041159_pr*c3*gA2)/fpi4
    ctr0r0(2,1,27) = (15.0933133539992_pr*c1*gA2)/fpi4 + &
         (7.47177985304339_pr*c3*gA2)/fpi4
    ctr0r0(3,1,27) = (-27.1326378098226_pr*c1*gA2)/fpi4 - &
         (12.7192853287926_pr*c3*gA2)/fpi4
    ctr0r0(1,1,29) = (-0.00618211760948692_pr*c3*gA2)/fpi4
    ctr0r0(2,1,29) = (0.541848695772724_pr*c1*gA2)/fpi4 + &
         (0.315695760130066_pr*c3*gA2)/fpi4
    ctr0r0(3,1,29) = (-2.02212634319308_pr*c1*gA2)/fpi4 - &
         (1.01238124519548_pr*c3*gA2)/fpi4
    ctr0r0(2,1,31) = (0.00309105880474346_pr*c3*gA2)/fpi4
    ctr0r0(3,1,31) = (-0.0677310869715905_pr*c1*gA2)/fpi4 - &
         (0.0402347347174442_pr*c3*gA2)/fpi4
    ctr0r0(3,1,33) = (-0.000386382350592933_pr*c3*gA2)/fpi4
    !
    !  Urho1rho1
    !
    ctr1r1=0.0_pr
    ctr1r1(1,1,1) = (-0.964246762695464_pr*c3*gA2)/fpi4 - &
         (0.0263462308309038_pr*c4*gA2)/fpi4 - (0.1875_pr*ce)/(fpi4*LambdaX) + &
         (0.046875_pr*cd*ga)/(fpi4*LambdaX)
    ctr1r1(1,1,3) = (-6.84195361981036_pr*c1*gA2)/fpi4 + &
         (5.01410456224811_pr*c3*gA2)/fpi4 + (0.0732776804846939_pr*c4*gA2)/fpi4 &
         - (0.06328125_pr*cd*ga)/(fpi4*LambdaX)
    ctr1r1(2,1,3) = (0.147698198400372_pr*c3*gA2)/fpi4 + &
         (0.0703526785714286_pr*c4*gA2)/fpi4 - &
         (0.02109375_pr*cd*ga)/(fpi4*LambdaX)
    ctr1r1(1,2,4) = (-0.1875_pr*c3*gA2)/fpi4 - (0.375_pr*c4*gA2)/fpi4 + &
         (0.140625_pr*cd*ga)/(fpi4*LambdaX)
    ctr1r1(1,1,5) = (-23.1434396482781_pr*c1*gA2)/fpi4 + &
         (12.2465870768616_pr*c3*gA2)/fpi4 - (0.18795405731824_pr*c4*gA2)/fpi4 + &
         (0.017578125_pr*cd*ga)/(fpi4*LambdaX)
    ctr1r1(2,1,5) = (35.0618960855328_pr*c1*gA2)/fpi4 - &
         (16.8781682906031_pr*c3*gA2)/fpi4 + (0.0543163105867347_pr*c4*gA2)/fpi4 &
         - (0.06328125_pr*cd*ga)/(fpi4*LambdaX)
    ctr1r1(1,2,6) = (-325.455281192602_pr*c1*gA2)/fpi4 + &
         (83.7334611805583_pr*c3*gA2)/fpi4 + (1.3574300625_pr*c4*gA2)/fpi4
    ctr1r1(2,2,6) = (32.9911875098339_pr*c3*gA2)/fpi4 + &
         (0.16875_pr*c4*gA2)/fpi4
    ctr1r1(1,1,7) = (-5761.43588948783_pr*c1*gA2)/fpi4 - &
         (869.283572934387_pr*c3*gA2)/fpi4 - (13.4068326170281_pr*c4*gA2)/fpi4
    ctr1r1(2,1,7) = (2192.10779645869_pr*c1*gA2)/fpi4 + &
         (263.842952909344_pr*c3*gA2)/fpi4 + (3.99388560267857_pr*c4*gA2)/fpi4 - &
         (0.00439453125_pr*cd*ga)/(fpi4*LambdaX)
    ctr1r1(3,1,7) = (-198.153843790565_pr*c1*gA2)/fpi4 - &
         (58.742085950407_pr*c3*gA2)/fpi4 - (0.4904296875_pr*c4*gA2)/fpi4
    ctr1r1(1,3,7) = (-526.213527565944_pr*c3*gA2)/fpi4 - &
         (2.84068526785714_pr*c4*gA2)/fpi4
    ctr1r1(1,2,8) = (11194.2459766741_pr*c1*gA2)/fpi4 - &
         (901.889210563564_pr*c3*gA2)/fpi4 + (1.90516445471939_pr*c4*gA2)/fpi4
    ctr1r1(2,2,8) = (-751.043571428571_pr*c1*gA2)/fpi4 + &
         (1336.68785503856_pr*c3*gA2)/fpi4 + (4.15185267857143_pr*c4*gA2)/fpi4
    ctr1r1(1,1,9) = (-19577.5085802329_pr*c1*gA2)/fpi4 - &
         (20061.3662796935_pr*c3*gA2)/fpi4 - (402.203416845285_pr*c4*gA2)/fpi4
    ctr1r1(2,1,9) = (6927.60941249884_pr*c1*gA2)/fpi4 + &
         (10204.912096572_pr*c3*gA2)/fpi4 + (157.45919188058_pr*c4*gA2)/fpi4
    ctr1r1(3,1,9) = (-1688.44574624066_pr*c1*gA2)/fpi4 - &
         (1848.73711933536_pr*c3*gA2)/fpi4 - (17.251846749442_pr*c4*gA2)/fpi4
    ctr1r1(1,3,9) = (-14639.9109308036_pr*c1*gA2)/fpi4 - &
         (2073.4495546875_pr*c3*gA2)/fpi4 + (12.406921875_pr*c4*gA2)/fpi4
    ctr1r1(1,2,10) = (-6051.26397950415_pr*c1*gA2)/fpi4 + &
         (25318.5711035127_pr*c3*gA2)/fpi4 + (1111.77481741071_pr*c4*gA2)/fpi4
    ctr1r1(2,2,10) = (20263.462438058_pr*c1*gA2)/fpi4 - &
         (645.573171294818_pr*c3*gA2)/fpi4 - (220.10902734375_pr*c4*gA2)/fpi4
    ctr1r1(1,1,11) = (-117612.589511906_pr*c1*gA2)/fpi4 - &
         (22681.6420225835_pr*c3*gA2)/fpi4 + (1494.47871529855_pr*c4*gA2)/fpi4
    ctr1r1(2,1,11) = (78131.8227221231_pr*c1*gA2)/fpi4 + &
         (11306.3195229362_pr*c3*gA2)/fpi4 - (1044.42429449456_pr*c4*gA2)/fpi4
    ctr1r1(3,1,11) = (-19244.1407277278_pr*c1*gA2)/fpi4 - &
         (4031.58809797695_pr*c3*gA2)/fpi4 + (143.952209472656_pr*c4*gA2)/fpi4
    ctr1r1(1,3,11) = (-12782.365790625_pr*c1*gA2)/fpi4 - &
         (19898.5801884375_pr*c3*gA2)/fpi4 - (763.061715_pr*c4*gA2)/fpi4
    ctr1r1(1,2,12) = (112532.892207631_pr*c1*gA2)/fpi4 + &
         (2673.93101274094_pr*c3*gA2)/fpi4 - (3075.39951940848_pr*c4*gA2)/fpi4
    ctr1r1(2,2,12) = (-4174.98147321429_pr*c1*gA2)/fpi4 + &
         (16959.475747232_pr*c3*gA2)/fpi4 + (1249.48156138393_pr*c4*gA2)/fpi4
    ctr1r1(1,1,13) = (-62220.7516704385_pr*c1*gA2)/fpi4 - &
         (37645.4468013683_pr*c3*gA2)/fpi4 - (1458.88122707611_pr*c4*gA2)/fpi4
    ctr1r1(2,1,13) = (49392.5657264478_pr*c1*gA2)/fpi4 + &
         (35938.6848928339_pr*c3*gA2)/fpi4 + (1614.83538151428_pr*c4*gA2)/fpi4
    ctr1r1(3,1,13) = (-21590.8338635122_pr*c1*gA2)/fpi4 - &
         (12971.4689072793_pr*c3*gA2)/fpi4 - (402.76586151123_pr*c4*gA2)/fpi4
    ctr1r1(1,3,13) = (-66293.2611893973_pr*c1*gA2)/fpi4 - &
         (8811.70601861049_pr*c3*gA2)/fpi4 + (1385.87636601563_pr*c4*gA2)/fpi4
    ctr1r1(1,2,14) = (-8634.9144059954_pr*c1*gA2)/fpi4 + &
         (13584.0523217795_pr*c3*gA2)/fpi4 + (1941.14579905622_pr*c4*gA2)/fpi4
    ctr1r1(2,2,14) = (46092.3086405378_pr*c1*gA2)/fpi4 + &
         (6421.96926354563_pr*c3*gA2)/fpi4 - (1169.36919344657_pr*c4*gA2)/fpi4
    ctr1r1(1,1,15) = (-54823.1734859566_pr*c1*gA2)/fpi4 - &
         (21478.1701524643_pr*c3*gA2)/fpi4 + (30.6758968031093_pr*c4*gA2)/fpi4
    ctr1r1(2,1,15) = (79796.8451847046_pr*c1*gA2)/fpi4 + &
         (29550.3190922257_pr*c3*gA2)/fpi4 - (265.149560200718_pr*c4*gA2)/fpi4
    ctr1r1(3,1,15) = (-34549.4910141965_pr*c1*gA2)/fpi4 - &
         (12487.0361763301_pr*c3*gA2)/fpi4 + (129.866342965262_pr*c4*gA2)/fpi4
    ctr1r1(1,3,15) = (-6841.585265625_pr*c1*gA2)/fpi4 - &
         (8394.50421527972_pr*c3*gA2)/fpi4 - (637.487162214007_pr*c4*gA2)/fpi4
    ctr1r1(1,2,16) = (-966.144667675781_pr*c1*gA2)/fpi4 - &
         (1635.80437776282_pr*c3*gA2)/fpi4 - (20.2837476928711_pr*c4*gA2)/fpi4
    ctr1r1(2,2,16) = (10019.4643315377_pr*c1*gA2)/fpi4 + &
         (5806.54125301064_pr*c3*gA2)/fpi4 + (162.342586296387_pr*c4*gA2)/fpi4
    ctr1r1(1,1,17) = (-27328.2660786696_pr*c1*gA2)/fpi4 - &
         (11536.8965413_pr*c3*gA2)/fpi4 - (0.153939957765852_pr*c4*gA2)/fpi4
    ctr1r1(2,1,17) = (51084.8335772702_pr*c1*gA2)/fpi4 + &
         (21317.5006264206_pr*c3*gA2)/fpi4 - (5.04213459735325_pr*c4*gA2)/fpi4
    ctr1r1(3,1,17) = (-23176.714608858_pr*c1*gA2)/fpi4 - &
         (9591.77377013768_pr*c3*gA2)/fpi4 - (7.01374572372437_pr*c4*gA2)/fpi4
    ctr1r1(1,3,17) = (-839.926762520926_pr*c1*gA2)/fpi4
    ctr1r1(1,2,18) = (-71.6191333138602_pr*c1*gA2)/fpi4 - &
         (29.5562121168625_pr*c3*gA2)/fpi4
    ctr1r1(2,2,18) = (1117.2724960968_pr*c1*gA2)/fpi4 + &
         (423.729200499137_pr*c3*gA2)/fpi4 + (5.07093692321777_pr*c4*gA2)/fpi4
    ctr1r1(1,1,19) = (-9253.22855094582_pr*c1*gA2)/fpi4 - &
         (4005.49423137719_pr*c3*gA2)/fpi4
    ctr1r1(2,1,19) = (21870.1711895189_pr*c1*gA2)/fpi4 + &
         (9310.62118848288_pr*c3*gA2)/fpi4 + (0.0769699788829259_pr*c4*gA2)/fpi4
    ctr1r1(3,1,19) = (-11369.3247545611_pr*c1*gA2)/fpi4 - &
         (4767.97117814042_pr*c3*gA2)/fpi4 - (0.656709900856018_pr*c4*gA2)/fpi4
    ctr1r1(2,2,20) = (17.9047833284651_pr*c1*gA2)/fpi4 + &
         (7.38905302921563_pr*c3*gA2)/fpi4
    ctr1r1(1,1,21) = (-2139.34177995088_pr*c1*gA2)/fpi4 - &
         (951.965547431256_pr*c3*gA2)/fpi4
    ctr1r1(2,1,21) = (6585.92097742191_pr*c1*gA2)/fpi4 + &
         (2872.41393631472_pr*c3*gA2)/fpi4
    ctr1r1(3,1,21) = (-4206.97567815619_pr*c1*gA2)/fpi4 - &
         (1804.69962149452_pr*c3*gA2)/fpi4 - &
         (0.00962124736036573_pr*c4*gA2)/fpi4
    ctr1r1(1,1,23) = (-320.600347402137_pr*c1*gA2)/fpi4 - &
         (147.828203709359_pr*c3*gA2)/fpi4
    ctr1r1(2,1,23) = (1373.03667928303_pr*c1*gA2)/fpi4 + &
         (615.472966679238_pr*c3*gA2)/fpi4
    ctr1r1(3,1,23) = (-1139.80667854591_pr*c1*gA2)/fpi4 - &
         (500.614173373695_pr*c3*gA2)/fpi4
    ctr1r1(1,1,25) = (-28.0192319249075_pr*c1*gA2)/fpi4 - &
         (13.7137479594838_pr*c3*gA2)/fpi4
    ctr1r1(2,1,25) = (187.596940698279_pr*c1*gA2)/fpi4 + &
         (87.2234077809484_pr*c3*gA2)/fpi4
    ctr1r1(3,1,25) = (-216.193884091206_pr*c1*gA2)/fpi4 - &
         (97.5972978305762_pr*c3*gA2)/fpi4
    ctr1r1(1,1,27) = (-1.08369739154545_pr*c1*gA2)/fpi4 - &
         (0.619027285041159_pr*c3*gA2)/fpi4
    ctr1r1(2,1,27) = (15.0933133539992_pr*c1*gA2)/fpi4 + &
         (7.47177985304339_pr*c3*gA2)/fpi4
    ctr1r1(3,1,27) = (-27.1326378098226_pr*c1*gA2)/fpi4 - &
         (12.7192853287926_pr*c3*gA2)/fpi4
    ctr1r1(1,1,29) = (-0.00618211760948692_pr*c3*gA2)/fpi4
    ctr1r1(2,1,29) = (0.541848695772724_pr*c1*gA2)/fpi4 + &
         (0.315695760130066_pr*c3*gA2)/fpi4
    ctr1r1(3,1,29) = (-2.02212634319308_pr*c1*gA2)/fpi4 - &
         (1.01238124519548_pr*c3*gA2)/fpi4
    ctr1r1(2,1,31) = (0.00309105880474346_pr*c3*gA2)/fpi4
    ctr1r1(3,1,31) = (-0.0677310869715905_pr*c1*gA2)/fpi4 - &
         (0.0402347347174442_pr*c3*gA2)/fpi4
    ctr1r1(3,1,33) = (-0.000386382350592933_pr*c3*gA2)/fpi4
    !
    !  UDrho0Drho0
    !
    ctdr0dr0=0.0_pr
    !
    !  UDrho1Drho1
    !
    ctdr1dr1=0.0_pr
    !
    !  Urho0DDrho0
    !
    ctr0Dr0=0.0_pr
    ctr0Dr0(1,1,3) = (-0.682809650873706_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0239679050048591_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,5) = (-5.25041224998774_pr*c1*gA2)/(fpi4*mpi2) + &
         (3.77245681259736_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.166505020727041_pr*c4*gA2)/(fpi4*mpi2) + &
         (0.03515625_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr0Dr0(2,1,5) = (0.12308183200031_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0586272321428571_pr*c4*gA2)/(fpi4*mpi2) - &
         (0.017578125_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr0Dr0(1,1,7) = (-26.0336741183036_pr*c1*gA2)/(fpi4*mpi2) + &
         (12.9213033479161_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0724766691246811_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,7) = (29.218246737944_pr*c1*gA2)/(fpi4*mpi2) - &
         (13.9991785771689_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.0572336176658163_pr*c4*gA2)/(fpi4*mpi2) - &
         (0.0087890625_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr0Dr0(1,2,8) = (-252.855350558036_pr*c1*gA2)/(fpi4*mpi2) + &
         (56.1626339963032_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.49978546875_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,2,8) = (27.4926562581949_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.140625_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,9) = (-4737.41799564353_pr*c1*gA2)/(fpi4*mpi2) - &
         (732.936270957556_pr*c3*gA2)/(fpi4*mpi2) - &
         (10.8826855699936_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,9) = (1812.39349332203_pr*c1*gA2)/(fpi4*mpi2) + &
         (230.621448631588_pr*c3*gA2)/(fpi4*mpi2) + &
         (3.48044224330357_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,9) = (-165.128203158804_pr*c1*gA2)/(fpi4*mpi2) - &
         (48.9517382920058_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.40869140625_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,3,9) = (-412.495287388393_pr*c3*gA2)/(fpi4*mpi2) - &
         (1.89848772321429_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,2,10) = (8987.88523861607_pr*c1*gA2)/(fpi4*mpi2) - &
         (631.430929909134_pr*c3*gA2)/(fpi4*mpi2) + &
         (1.32764820631378_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,2,10) = (-625.869642857143_pr*c1*gA2)/(fpi4*mpi2) + &
         (1045.17490521998_pr*c3*gA2)/(fpi4*mpi2) + &
         (3.10831473214286_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,11) = (-17022.0529432446_pr*c1*gA2)/(fpi4*mpi2) - &
         (15425.6442847704_pr*c3*gA2)/(fpi4*mpi2) - &
         (325.368422118792_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,11) = (6144.54246637738_pr*c1*gA2)/(fpi4*mpi2) + &
         (7927.36933047731_pr*c3*gA2)/(fpi4*mpi2) + &
         (127.556172677176_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,11) = (-1407.03812186722_pr*c1*gA2)/(fpi4*mpi2) - &
         (1455.27334812581_pr*c3*gA2)/(fpi4*mpi2) - &
         (13.9398324148996_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,3,11) = (-11514.015703125_pr*c1*gA2)/(fpi4*mpi2) - &
         (1727.87462890625_pr*c3*gA2)/(fpi4*mpi2) + &
         (10.3391015625_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,2,12) = (-694.446624441964_pr*c1*gA2)/(fpi4*mpi2) + &
         (18247.3109917788_pr*c3*gA2)/(fpi4*mpi2) + &
         (899.399912946429_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,2,12) = (15092.8155175781_pr*c1*gA2)/(fpi4*mpi2) - &
         (116.351373250285_pr*c3*gA2)/(fpi4*mpi2) - &
         (178.470087890625_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,13) = (-84816.595261174_pr*c1*gA2)/(fpi4*mpi2) - &
         (21005.6538813155_pr*c3*gA2)/(fpi4*mpi2) + &
         (1207.97363186907_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,13) = (56932.8731897877_pr*c1*gA2)/(fpi4*mpi2) + &
         (11075.3299701944_pr*c3*gA2)/(fpi4*mpi2) - &
         (844.741128183943_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,13) = (-14559.3117129491_pr*c1*gA2)/(fpi4*mpi2) - &
         (3504.49103169193_pr*c3*gA2)/(fpi4*mpi2) + &
         (116.596984863281_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,3,13) = (-12860.985515625_pr*c1*gA2)/(fpi4*mpi2) - &
         (14398.998496875_pr*c3*gA2)/(fpi4*mpi2) - &
         (616.763278125_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,2,14) = (71603.420010498_pr*c1*gA2)/(fpi4*mpi2) + &
         (8297.07769839308_pr*c3*gA2)/(fpi4*mpi2) - &
         (2483.1196515904_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,2,14) = (2438.09865025112_pr*c1*gA2)/(fpi4*mpi2) + &
         (10960.7457584461_pr*c3*gA2)/(fpi4*mpi2) + &
         (1009.19850167411_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,15) = (-55860.3376653923_pr*c1*gA2)/(fpi4*mpi2) - &
         (25835.4946732797_pr*c3*gA2)/(fpi4*mpi2) - &
         (1174.87146717911_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,15) = (47592.1515045526_pr*c1*gA2)/(fpi4*mpi2) + &
         (24786.4294824231_pr*c3*gA2)/(fpi4*mpi2) + &
         (1301.50458021698_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,15) = (-19018.8059498995_pr*c1*gA2)/(fpi4*mpi2) - &
         (9429.01449770894_pr*c3*gA2)/(fpi4*mpi2) - &
         (324.792194366455_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,3,15) = (-44074.7843378906_pr*c1*gA2)/(fpi4*mpi2) - &
         (9754.29660498047_pr*c3*gA2)/(fpi4*mpi2) + &
         (1118.35204980469_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,2,16) = (1466.81527313232_pr*c1*gA2)/(fpi4*mpi2) + &
         (4924.95415195852_pr*c3*gA2)/(fpi4*mpi2) + &
         (1563.11953388149_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,2,16) = (29864.464335022_pr*c1*gA2)/(fpi4*mpi2) + &
         (7829.64002841679_pr*c3*gA2)/(fpi4*mpi2) - &
         (942.380516836984_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,17) = (-42578.4126891706_pr*c1*gA2)/(fpi4*mpi2) - &
         (17202.0887780327_pr*c3*gA2)/(fpi4*mpi2) + &
         (24.7012874970572_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,17) = (61948.0344894364_pr*c1*gA2)/(fpi4*mpi2) + &
         (24461.7278379897_pr*c3*gA2)/(fpi4*mpi2) - &
         (213.463917761666_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,17) = (-26551.291288457_pr*c1*gA2)/(fpi4*mpi2) - &
         (10248.0190226489_pr*c3*gA2)/(fpi4*mpi2) + &
         (104.612325701032_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,3,17) = (-7839.31645019531_pr*c1*gA2)/(fpi4*mpi2) - &
         (4409.58359893799_pr*c3*gA2)/(fpi4*mpi2) - &
         (513.293134852818_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,2,18) = (-2406.80499142456_pr*c1*gA2)/(fpi4*mpi2) - &
         (982.897425834421_pr*c3*gA2)/(fpi4*mpi2) - &
         (16.3321086730957_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,2,18) = (8676.01512762451_pr*c1*gA2)/(fpi4*mpi2) + &
         (3669.79377386557_pr*c3*gA2)/(fpi4*mpi2) + &
         (130.679305718994_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,19) = (-21929.7613258823_pr*c1*gA2)/(fpi4*mpi2) - &
         (9003.06173155003_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.123949706677028_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,19) = (41084.26511068_pr*c1*gA2)/(fpi4*mpi2) + &
         (16712.9487662382_pr*c3*gA2)/(fpi4*mpi2) - &
         (4.06063970530374_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,19) = (-18477.8744897258_pr*c1*gA2)/(fpi4*mpi2) - &
         (7453.07953912234_pr*c3*gA2)/(fpi4*mpi2) - &
         (5.64185216903687_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,2,20) = (601.70124785614_pr*c1*gA2)/(fpi4*mpi2) + &
         (245.724356458605_pr*c3*gA2)/(fpi4*mpi2) + &
         (4.08302716827393_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,21) = (-7429.54966559941_pr*c1*gA2)/(fpi4*mpi2) - &
         (3087.22338867281_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,21) = (17522.7958314409_pr*c1*gA2)/(fpi4*mpi2) + &
         (7223.97599005467_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0619748533385141_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,21) = (-9051.54263718816_pr*c1*gA2)/(fpi4*mpi2) - &
         (3700.96455501676_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.528670542240143_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,23) = (-1722.74624112667_pr*c1*gA2)/(fpi4*mpi2) - &
         (722.002478964633_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,23) = (5292.28471362693_pr*c1*gA2)/(fpi4*mpi2) + &
         (2204.46546742604_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,23) = (-3370.28377443923_pr*c1*gA2)/(fpi4*mpi2) - &
         (1394.125888934_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.00774685666731426_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,25) = (-258.704980303594_pr*c1*gA2)/(fpi4*mpi2) - &
         (108.950617217712_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,25) = (1106.14939761973_pr*c1*gA2)/(fpi4*mpi2) + &
         (464.078201799952_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,25) = (-916.533683372315_pr*c1*gA2)/(fpi4*mpi2) - &
         (382.506711741246_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,27) = (-22.6474359092818_pr*c1*gA2)/(fpi4*mpi2) - &
         (9.54811915326822_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,27) = (151.415132381586_pr*c1*gA2)/(fpi4*mpi2) + &
         (63.7775488277401_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,27) = (-174.235171472089_pr*c1*gA2)/(fpi4*mpi2) - &
         (73.1584858358234_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,29) = (-0.877190519238951_pr*c1*gA2)/(fpi4*mpi2) - &
         (0.368818401576571_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,29) = (12.2009084738798_pr*c1*gA2)/(fpi4*mpi2) + &
         (5.14287797821068_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,29) = (-21.9040194562317_pr*c1*gA2)/(fpi4*mpi2) - &
         (9.22717823122213_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,31) = (0.438595259619476_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.184409200788286_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,31) = (-1.63476237413985_pr*c1*gA2)/(fpi4*mpi2) - &
         (0.688962047473406_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,33) = (-0.0548244074524344_pr*c1*gA2)/(fpi4*mpi2) - &
         (0.0230511500985357_pr*c3*gA2)/(fpi4*mpi2)
    !
    !  Urho1DDrho0
    !
    ctr1Dr0=0.0_pr
    ctr1Dr0(1,1,3) = (0.0185237085762877_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.00798930166828636_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,1,5) = (0.0956155117984694_pr*c1*gA2)/(fpi4*mpi2) - &
         (0.0234414046556122_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0555016735756803_pr*c4*gA2)/(fpi4*mpi2) - &
         (0.01171875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1Dr0(2,1,5) = (-0.00389508928571429_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.0195424107142857_pr*c4*gA2)/(fpi4*mpi2) + &
         (0.005859375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1Dr0(1,1,7) = (-0.382716029575893_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.159635060586735_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.024158889708227_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,1,7) = (-0.0256136798469388_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.00456883769132653_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0190778725552721_pr*c4*gA2)/(fpi4*mpi2) + &
         (0.0029296875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1Dr0(1,2,8) = (0.706878683035714_pr*c1*gA2)/(fpi4*mpi2) - &
         (0.423425189732143_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.16659515625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,2,8) = (-0.0234375_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.046875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,1,9) = (11.4067307059152_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.609889416055485_pr*c3*gA2)/(fpi4*mpi2) + &
         (3.62756185666454_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,1,9) = (-3.54638253348214_pr*c1*gA2)/(fpi4*mpi2) - &
         (0.0995982142857143_pr*c3*gA2)/(fpi4*mpi2) - &
         (1.16014741443452_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(3,1,9) = (0.25048828125_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.04833984375_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.13623046875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,3,9) = (0.949243861607143_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.632829241071429_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,2,10) = (-22.8885294642857_pr*c1*gA2)/(fpi4*mpi2) + &
         (18.9939639150191_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.442549402104592_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,2,10) = (2.53108258928571_pr*c1*gA2)/(fpi4*mpi2) - &
         (4.85298549107143_pr*c3*gA2)/(fpi4*mpi2) - &
         (1.03610491071429_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,1,11) = (220.437076035854_pr*c1*gA2)/(fpi4*mpi2) + &
         (298.906979832589_pr*c3*gA2)/(fpi4*mpi2) + &
         (108.456140706264_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,1,11) = (-68.8102811802455_pr*c1*gA2)/(fpi4*mpi2) - &
         (121.491437360491_pr*c3*gA2)/(fpi4*mpi2) - &
         (42.5187242257254_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(3,1,11) = (6.13759068080357_pr*c1*gA2)/(fpi4*mpi2) + &
         (12.8441554478237_pr*c3*gA2)/(fpi4*mpi2) + &
         (4.64661080496652_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,3,11) = (11.318203125_pr*c1*gA2)/(fpi4*mpi2) - &
         (23.38048828125_pr*c3*gA2)/(fpi4*mpi2) - &
         (3.4463671875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,2,12) = (-732.624496372768_pr*c1*gA2)/(fpi4*mpi2) - &
         (842.4480234375_pr*c3*gA2)/(fpi4*mpi2) - &
         (299.799970982143_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,2,12) = (117.953701171875_pr*c1*gA2)/(fpi4*mpi2) + &
         (175.760551757812_pr*c3*gA2)/(fpi4*mpi2) + &
         (59.490029296875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,1,13) = (-1470.77295530744_pr*c1*gA2)/(fpi4*mpi2) - &
         (1177.66729674159_pr*c3*gA2)/(fpi4*mpi2) - &
         (402.65787728969_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,1,13) = (906.809157364328_pr*c1*gA2)/(fpi4*mpi2) + &
         (809.938923974958_pr*c3*gA2)/(fpi4*mpi2) + &
         (281.580376061314_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(3,1,13) = (-110.580871582031_pr*c1*gA2)/(fpi4*mpi2) - &
         (112.623596191406_pr*c3*gA2)/(fpi4*mpi2) - &
         (38.8656616210938_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,3,13) = (585.845296875_pr*c1*gA2)/(fpi4*mpi2) + &
         (588.08008125_pr*c3*gA2)/(fpi4*mpi2) + &
         (205.587759375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,2,14) = (3251.26034692383_pr*c1*gA2)/(fpi4*mpi2) + &
         (2508.40495063477_pr*c3*gA2)/(fpi4*mpi2) + &
         (827.706550530134_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,2,14) = (-1184.58378557478_pr*c1*gA2)/(fpi4*mpi2) - &
         (995.956718052455_pr*c3*gA2)/(fpi4*mpi2) - &
         (336.399500558036_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,1,15) = (1877.99099063873_pr*c1*gA2)/(fpi4*mpi2) + &
         (1378.17489546073_pr*c3*gA2)/(fpi4*mpi2) + &
         (391.623822393036_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,1,15) = (-1930.62984776132_pr*c1*gA2)/(fpi4*mpi2) - &
         (1432.9768855957_pr*c3*gA2)/(fpi4*mpi2) - &
         (433.834860072327_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(3,1,15) = (440.209643554687_pr*c1*gA2)/(fpi4*mpi2) + &
         (341.656391143799_pr*c3*gA2)/(fpi4*mpi2) + &
         (108.264064788818_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,3,15) = (-1549.60388085937_pr*c1*gA2)/(fpi4*mpi2) - &
         (1159.75179052734_pr*c3*gA2)/(fpi4*mpi2) - &
         (372.784016601562_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,2,16) = (-2582.26988983154_pr*c1*gA2)/(fpi4*mpi2) - &
         (1875.22398799962_pr*c3*gA2)/(fpi4*mpi2) - &
         (521.039844627162_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,2,16) = (1459.93622772217_pr*c1*gA2)/(fpi4*mpi2) + &
         (1058.4680922808_pr*c3*gA2)/(fpi4*mpi2) + &
         (314.126838945661_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,1,17) = (-56.9118636474609_pr*c1*gA2)/(fpi4*mpi2) - &
         (35.9791338062831_pr*c3*gA2)/(fpi4*mpi2) - &
         (8.23376249901908_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,1,17) = (396.575549217224_pr*c1*gA2)/(fpi4*mpi2) + &
         (276.498193133763_pr*c3*gA2)/(fpi4*mpi2) + &
         (71.1546392538888_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(3,1,17) = (-182.383999263218_pr*c1*gA2)/(fpi4*mpi2) - &
         (125.990329873221_pr*c3*gA2)/(fpi4*mpi2) - &
         (34.8707752336775_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,3,17) = (871.035161132812_pr*c1*gA2)/(fpi4*mpi2) + &
         (629.940514133998_pr*c3*gA2)/(fpi4*mpi2) + &
         (171.097711617606_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,2,18) = (38.1077882995605_pr*c1*gA2)/(fpi4*mpi2) + &
         (24.4976645050049_pr*c3*gA2)/(fpi4*mpi2) + &
         (5.44403622436523_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,2,18) = (-244.521582824707_pr*c1*gA2)/(fpi4*mpi2) - &
         (173.383349386597_pr*c3*gA2)/(fpi4*mpi2) - &
         (43.559768572998_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,1,19) = (0.340248109817505_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.240601839556013_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0413165688923427_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,1,19) = (9.06178956413269_pr*c1*gA2)/(fpi4*mpi2) + &
         (5.50013281108311_pr*c3*gA2)/(fpi4*mpi2) + &
         (1.35354656843458_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(3,1,19) = (10.8524751663208_pr*c1*gA2)/(fpi4*mpi2) + &
         (8.40921586990356_pr*c3*gA2)/(fpi4*mpi2) + &
         (1.88061738967896_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,2,20) = (-9.52694707489014_pr*c1*gA2)/(fpi4*mpi2) - &
         (6.12441612625122_pr*c3*gA2)/(fpi4*mpi2) - &
         (1.36100905609131_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,1,21) = (-0.170124054908752_pr*c1*gA2)/(fpi4*mpi2) - &
         (0.120300919778006_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.0206582844461714_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(3,1,21) = (1.29154408693314_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.873662660121918_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.176223514080048_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(3,1,23) = (0.0212655068635941_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.0150376149722508_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.00258228555577142_pr*c4*gA2)/(fpi4*mpi2)
    !
    !  Urho1DDrho1
    !
    ctr1Dr1=0.0_pr
    ctr1Dr1(1,1,3) = (0.225058111718187_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.0159786033365727_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,5) = (1.76926855346105_pr*c1*gA2)/(fpi4*mpi2) - &
         (1.25129754669487_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.111003347151361_pr*c4*gA2)/(fpi4*mpi2) - &
         (0.0234375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1Dr1(2,1,5) = (-0.0449223666191509_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.0390848214285714_pr*c4*gA2)/(fpi4*mpi2) + &
         (0.01171875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1Dr1(1,1,7) = (8.62173784319196_pr*c1*gA2)/(fpi4*mpi2) - &
         (4.30128936533428_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.0483177794164541_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,7) = (-9.7579980091616_pr*c1*gA2)/(fpi4*mpi2) + &
         (4.67405890125441_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0381557451105442_pr*c4*gA2)/(fpi4*mpi2) + &
         (0.005859375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1Dr1(1,2,8) = (84.4669955357143_pr*c1*gA2)/(fpi4*mpi2) - &
         (18.660361268857_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.3331903125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,2,8) = (-9.18765625273164_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.09375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,9) = (1573.18444149334_pr*c1*gA2)/(fpi4*mpi2) + &
         (246.450668847167_pr*c3*gA2)/(fpi4*mpi2) + &
         (7.25512371332908_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,9) = (-602.738093849159_pr*c1*gA2)/(fpi4*mpi2) - &
         (77.5124918057675_pr*c3*gA2)/(fpi4*mpi2) - &
         (2.32029482886905_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,9) = (54.959238292518_pr*c1*gA2)/(fpi4*mpi2) + &
         (16.3655859410853_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.2724609375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,3,9) = (137.182014508929_pr*c3*gA2)/(fpi4*mpi2) + &
         (1.26565848214286_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,2,10) = (-2950.89558816964_pr*c1*gA2)/(fpi4*mpi2) + &
         (192.771436581372_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.885098804209184_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,2,10) = (202.294921875_pr*c1*gA2)/(fpi4*mpi2) - &
         (345.017723242968_pr*c3*gA2)/(fpi4*mpi2) - &
         (2.07220982142857_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,11) = (6073.45819058092_pr*c1*gA2)/(fpi4*mpi2) + &
         (4996.23922290246_pr*c3*gA2)/(fpi4*mpi2) + &
         (216.912281412528_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,11) = (-2212.39051736854_pr*c1*gA2)/(fpi4*mpi2) - &
         (2578.84096274467_pr*c3*gA2)/(fpi4*mpi2) - &
         (85.0374484514509_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,11) = (484.396391719877_pr*c1*gA2)/(fpi4*mpi2) + &
         (478.331459408286_pr*c3*gA2)/(fpi4*mpi2) + &
         (9.29322160993304_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,3,11) = (3776.4796875_pr*c1*gA2)/(fpi4*mpi2) + &
         (595.892330729167_pr*c3*gA2)/(fpi4*mpi2) - &
         (6.892734375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,2,12) = (-1047.22431947545_pr*c1*gA2)/(fpi4*mpi2) - &
         (5677.84872382209_pr*c3*gA2)/(fpi4*mpi2) - &
         (599.599941964286_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,2,12) = (-4730.6293359375_pr*c1*gA2)/(fpi4*mpi2) - &
         (56.7874345025092_pr*c3*gA2)/(fpi4*mpi2) + &
         (118.98005859375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,13) = (26158.9112133749_pr*c1*gA2)/(fpi4*mpi2) + &
         (7519.86242027358_pr*c3*gA2)/(fpi4*mpi2) - &
         (805.315754579381_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,13) = (-17547.6432543683_pr*c1*gA2)/(fpi4*mpi2) - &
         (4066.12243819318_pr*c3*gA2)/(fpi4*mpi2) + &
         (563.160752122628_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,13) = (4641.24227101558_pr*c1*gA2)/(fpi4*mpi2) + &
         (1225.12479783611_pr*c3*gA2)/(fpi4*mpi2) - &
         (77.7313232421875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,3,13) = (5176.02796875_pr*c1*gA2)/(fpi4*mpi2) + &
         (4535.422284375_pr*c3*gA2)/(fpi4*mpi2) + &
         (411.17551875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,2,14) = (-19673.8226904297_pr*c1*gA2)/(fpi4*mpi2) - &
         (3790.5753035822_pr*c3*gA2)/(fpi4*mpi2) + &
         (1655.41310106027_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,2,14) = (-2415.02747628348_pr*c1*gA2)/(fpi4*mpi2) - &
         (3234.1744717979_pr*c3*gA2)/(fpi4*mpi2) - &
         (672.799001116071_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,15) = (20527.2675777519_pr*c1*gA2)/(fpi4*mpi2) + &
         (8144.9046506711_pr*c3*gA2)/(fpi4*mpi2) + &
         (783.247644786072_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,15) = (-17967.065175939_pr*c1*gA2)/(fpi4*mpi2) - &
         (7747.02659701377_pr*c3*gA2)/(fpi4*mpi2) - &
         (867.669720144653_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,15) = (6843.01542324125_pr*c1*gA2)/(fpi4*mpi2) + &
         (3014.94458858161_pr*c3*gA2)/(fpi4*mpi2) + &
         (216.528129577637_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,3,15) = (12800.9049609375_pr*c1*gA2)/(fpi4*mpi2) + &
         (3707.04542871094_pr*c3*gA2)/(fpi4*mpi2) - &
         (745.568033203125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,2,16) = (-3092.52618530273_pr*c1*gA2)/(fpi4*mpi2) - &
         (1011.19649671026_pr*c3*gA2)/(fpi4*mpi2) - &
         (1042.07968925432_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,2,16) = (-8378.82125244141_pr*c1*gA2)/(fpi4*mpi2) - &
         (2987.23517760232_pr*c3*gA2)/(fpi4*mpi2) + &
         (628.253677891323_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,17) = (14135.5592743402_pr*c1*gA2)/(fpi4*mpi2) + &
         (5746.10194483582_pr*c3*gA2)/(fpi4*mpi2) - &
         (16.4675249980382_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,17) = (-20256.3592693034_pr*c1*gA2)/(fpi4*mpi2) - &
         (8245.07712593439_pr*c3*gA2)/(fpi4*mpi2) + &
         (142.309278507778_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,17) = (8660.61555253689_pr*c1*gA2)/(fpi4*mpi2) + &
         (3459.42204756108_pr*c3*gA2)/(fpi4*mpi2) - &
         (69.7415504673549_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,3,17) = (3484.14064453125_pr*c1*gA2)/(fpi4*mpi2) + &
         (1259.881028268_pr*c3*gA2)/(fpi4*mpi2) + &
         (342.195423235212_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,2,18) = (840.376118774414_pr*c1*gA2)/(fpi4*mpi2) + &
         (319.466587109805_pr*c3*gA2)/(fpi4*mpi2) + &
         (10.8880724487305_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,2,18) = (-3131.19715759277_pr*c1*gA2)/(fpi4*mpi2) - &
         (1166.81519764533_pr*c3*gA2)/(fpi4*mpi2) - &
         (87.1195371459961_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,19) = (7310.26069007057_pr*c1*gA2)/(fpi4*mpi2) + &
         (3000.94037657016_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0826331377846854_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,19) = (-13685.5267014613_pr*c1*gA2)/(fpi4*mpi2) - &
         (5572.85595346119_pr*c3*gA2)/(fpi4*mpi2) + &
         (2.70709313686916_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,19) = (6169.2187169198_pr*c1*gA2)/(fpi4*mpi2) + &
         (2481.7781746784_pr*c3*gA2)/(fpi4*mpi2) + &
         (3.76123477935791_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,2,20) = (-210.094029693604_pr*c1*gA2)/(fpi4*mpi2) - &
         (79.8666467774513_pr*c3*gA2)/(fpi4*mpi2) - &
         (2.72201811218262_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,21) = (2476.5165551998_pr*c1*gA2)/(fpi4*mpi2) + &
         (1029.07446289094_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,21) = (-5841.10206786854_pr*c1*gA2)/(fpi4*mpi2) - &
         (2407.9518963783_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.0413165688923427_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,21) = (3018.45160491616_pr*c1*gA2)/(fpi4*mpi2) + &
         (1233.36858750781_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.352447028160095_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,23) = (574.248747042225_pr*c1*gA2)/(fpi4*mpi2) + &
         (240.667492988211_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,23) = (-1764.09490454231_pr*c1*gA2)/(fpi4*mpi2) - &
         (734.821822475346_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,23) = (1123.44919031994_pr*c1*gA2)/(fpi4*mpi2) + &
         (464.703617106343_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.00516457111154284_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,25) = (86.2349934345314_pr*c1*gA2)/(fpi4*mpi2) + &
         (36.3168724059041_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,25) = (-368.716465873243_pr*c1*gA2)/(fpi4*mpi2) - &
         (154.692733933317_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,25) = (305.511227790772_pr*c1*gA2)/(fpi4*mpi2) + &
         (127.502237247082_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,27) = (7.54914530309393_pr*c1*gA2)/(fpi4*mpi2) + &
         (3.18270638442274_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,27) = (-50.4717107938621_pr*c1*gA2)/(fpi4*mpi2) - &
         (21.25918294258_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,27) = (58.0783904906964_pr*c1*gA2)/(fpi4*mpi2) + &
         (24.3861619452745_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,29) = (0.292396839746317_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.12293946719219_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,29) = (-4.06696949129328_pr*c1*gA2)/(fpi4*mpi2) - &
         (1.71429265940356_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,29) = (7.30133981874389_pr*c1*gA2)/(fpi4*mpi2) + &
         (3.07572607707404_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,31) = (-0.146198419873159_pr*c1*gA2)/(fpi4*mpi2) - &
         (0.0614697335960952_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,31) = (0.54492079137995_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.229654015824469_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,33) = (0.0182748024841448_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.00768371669951191_pr*c3*gA2)/(fpi4*mpi2)
    !
    !  Urho0tau0
    !
    ctr0t0=0.0_pr
    ctr0t0(1,1,3)=(1.36561930174741_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0479358100097182_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,5)=(10.5008244999755_pr*c1*gA2)/(fpi4*mpi2)- &
         (7.54491362519471_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.333010041454082_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.0703125_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr0t0(2,1,5)=(-0.24616366400062_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.117254464285714_pr*c4*gA2)/(fpi4*mpi2)+ &
         (0.03515625_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr0t0(1,1,7)=(52.0673482366071_pr*c1*gA2)/(fpi4*mpi2)- &
         (25.8426066958322_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.144953338249362_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,7)=(-58.436493475888_pr*c1*gA2)/(fpi4*mpi2)+ &
         (27.9983571543377_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.114467235331633_pr*c4*gA2)/(fpi4*mpi2)+ &
         (0.017578125_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr0t0(1,2,8)=(505.710701116071_pr*c1*gA2)/(fpi4*mpi2)- &
         (112.325267992606_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.9995709375_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,2,8)=(-54.9853125163899_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.28125_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,9)=(9474.83599128706_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1465.87254191511_pr*c3*gA2)/(fpi4*mpi2)+ &
         (21.7653711399872_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,9)=(-3624.78698664406_pr*c1*gA2)/(fpi4*mpi2)- &
         (461.242897263176_pr*c3*gA2)/(fpi4*mpi2)- &
         (6.96088448660714_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,9)=(330.256406317608_pr*c1*gA2)/(fpi4*mpi2)+ &
         (97.9034765840117_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.8173828125_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,3,9)=(824.990574776786_pr*c3*gA2)/(fpi4*mpi2)+ &
         (3.79697544642857_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,2,10)=(-17975.7704772321_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1262.86185981827_pr*c3*gA2)/(fpi4*mpi2)- &
         (2.65529641262755_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,2,10)=(1251.73928571429_pr*c1*gA2)/(fpi4*mpi2)- &
         (2090.34981043995_pr*c3*gA2)/(fpi4*mpi2)- &
         (6.21662946428571_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,11)=(34044.1058864892_pr*c1*gA2)/(fpi4*mpi2)+ &
         (30851.2885695409_pr*c3*gA2)/(fpi4*mpi2)+ &
         (650.736844237584_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,11)=(-12289.0849327548_pr*c1*gA2)/(fpi4*mpi2)- &
         (15854.7386609546_pr*c3*gA2)/(fpi4*mpi2)- &
         (255.112345354353_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,11)=(2814.07624373444_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2910.54669625161_pr*c3*gA2)/(fpi4*mpi2)+ &
         (27.8796648297991_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,3,11)=(23028.03140625_pr*c1*gA2)/(fpi4*mpi2)+ &
         (3455.7492578125_pr*c3*gA2)/(fpi4*mpi2)- &
         (20.678203125_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,2,12)=(1388.89324888393_pr*c1*gA2)/(fpi4*mpi2)- &
         (36494.6219835576_pr*c3*gA2)/(fpi4*mpi2)- &
         (1798.79982589286_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,2,12)=(-30185.6310351562_pr*c1*gA2)/(fpi4*mpi2)+ &
         (232.70274650057_pr*c3*gA2)/(fpi4*mpi2)+ &
         (356.94017578125_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,13)=(169633.190522348_pr*c1*gA2)/(fpi4*mpi2)+ &
         (42011.3077626309_pr*c3*gA2)/(fpi4*mpi2)- &
         (2415.94726373814_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,13)=(-113865.746379575_pr*c1*gA2)/(fpi4*mpi2)- &
         (22150.6599403888_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1689.48225636789_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,13)=(29118.6234258982_pr*c1*gA2)/(fpi4*mpi2)+ &
         (7008.98206338386_pr*c3*gA2)/(fpi4*mpi2)- &
         (233.193969726563_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,3,13)=(25721.97103125_pr*c1*gA2)/(fpi4*mpi2)+ &
         (28797.99699375_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1233.52655625_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,2,14)=(-143206.840020996_pr*c1*gA2)/(fpi4*mpi2)- &
         (16594.1553967862_pr*c3*gA2)/(fpi4*mpi2)+ &
         (4966.2393031808_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,2,14)=(-4876.19730050223_pr*c1*gA2)/(fpi4*mpi2)- &
         (21921.4915168923_pr*c3*gA2)/(fpi4*mpi2)- &
         (2018.39700334821_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,15)=(111720.675330785_pr*c1*gA2)/(fpi4*mpi2)+ &
         (51670.9893465594_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2349.74293435822_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,15)=(-95184.3030091052_pr*c1*gA2)/(fpi4*mpi2)- &
         (49572.8589648463_pr*c3*gA2)/(fpi4*mpi2)- &
         (2603.00916043396_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,15)=(38037.6118997991_pr*c1*gA2)/(fpi4*mpi2)+ &
         (18858.0289954179_pr*c3*gA2)/(fpi4*mpi2)+ &
         (649.58438873291_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,3,15)=(88149.5686757812_pr*c1*gA2)/(fpi4*mpi2)+ &
         (19508.5932099609_pr*c3*gA2)/(fpi4*mpi2)- &
         (2236.70409960937_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,2,16)=(-2933.63054626465_pr*c1*gA2)/(fpi4*mpi2)- &
         (9849.90830391705_pr*c3*gA2)/(fpi4*mpi2)- &
         (3126.23906776297_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,2,16)=(-59728.9286700439_pr*c1*gA2)/(fpi4*mpi2)- &
         (15659.2800568336_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1884.76103367397_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,17)=(85156.8253783411_pr*c1*gA2)/(fpi4*mpi2)+ &
         (34404.1775560654_pr*c3*gA2)/(fpi4*mpi2)- &
         (49.4025749941145_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,17)=(-123896.068978873_pr*c1*gA2)/(fpi4*mpi2)- &
         (48923.4556759794_pr*c3*gA2)/(fpi4*mpi2)+ &
         (426.927835523333_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,17)=(53102.5825769139_pr*c1*gA2)/(fpi4*mpi2)+ &
         (20496.0380452978_pr*c3*gA2)/(fpi4*mpi2)- &
         (209.224651402065_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,3,17)=(15678.6329003906_pr*c1*gA2)/(fpi4*mpi2)+ &
         (8819.16719787598_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1026.58626970564_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,2,18)=(4813.60998284912_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1965.79485166884_pr*c3*gA2)/(fpi4*mpi2)+ &
         (32.6642173461914_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,2,18)=(-17352.030255249_pr*c1*gA2)/(fpi4*mpi2)- &
         (7339.58754773114_pr*c3*gA2)/(fpi4*mpi2)- &
         (261.358611437988_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,19)=(43859.5226517645_pr*c1*gA2)/(fpi4*mpi2)+ &
         (18006.1234631001_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.247899413354056_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,19)=(-82168.5302213599_pr*c1*gA2)/(fpi4*mpi2)- &
         (33425.8975324765_pr*c3*gA2)/(fpi4*mpi2)+ &
         (8.12127941060747_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,19)=(36955.7489794515_pr*c1*gA2)/(fpi4*mpi2)+ &
         (14906.1590782447_pr*c3*gA2)/(fpi4*mpi2)+ &
         (11.2837043380737_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,2,20)=(-1203.40249571228_pr*c1*gA2)/(fpi4*mpi2)- &
         (491.448712917211_pr*c3*gA2)/(fpi4*mpi2)- &
         (8.16605433654785_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,21)=(14859.0993311988_pr*c1*gA2)/(fpi4*mpi2)+ &
         (6174.44677734561_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,21)=(-35045.5916628818_pr*c1*gA2)/(fpi4*mpi2)- &
         (14447.9519801093_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.123949706677028_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,21)=(18103.0852743763_pr*c1*gA2)/(fpi4*mpi2)+ &
         (7401.92911003352_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1.05734108448029_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,23)=(3445.49248225335_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1444.00495792927_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,23)=(-10584.5694272539_pr*c1*gA2)/(fpi4*mpi2)- &
         (4408.93093485208_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,23)=(6740.56754887846_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2788.251777868_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0154937133346285_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,25)=(517.409960607188_pr*c1*gA2)/(fpi4*mpi2)+ &
         (217.901234435425_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,25)=(-2212.29879523946_pr*c1*gA2)/(fpi4*mpi2)- &
         (928.156403599904_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,25)=(1833.06736674463_pr*c1*gA2)/(fpi4*mpi2)+ &
         (765.013423482493_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,27)=(45.2948718185636_pr*c1*gA2)/(fpi4*mpi2)+ &
         (19.0962383065364_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,27)=(-302.830264763173_pr*c1*gA2)/(fpi4*mpi2)- &
         (127.55509765548_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,27)=(348.470342944178_pr*c1*gA2)/(fpi4*mpi2)+ &
         (146.316971671647_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,29)=(1.7543810384779_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.737636803153143_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,29)=(-24.4018169477597_pr*c1*gA2)/(fpi4*mpi2)- &
         (10.2857559564214_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,29)=(43.8080389124633_pr*c1*gA2)/(fpi4*mpi2)+ &
         (18.4543564624443_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,31)=(-0.877190519238951_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.368818401576571_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,31)=(3.2695247482797_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1.37792409494681_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,33)=(0.109648814904869_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0461023001970714_pr*c3*gA2)/(fpi4*mpi2)
    !
    !  Urho1tau0
    !
    ctr1t0=0.0_pr
    ctr1t0(1,1,3)=(-0.0370474171525753_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0159786033365727_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,1,5)=(-0.191231023596939_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0468828093112245_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.111003347151361_pr*c4*gA2)/(fpi4*mpi2)+ &
         (0.0234375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1t0(2,1,5)=(0.00779017857142857_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0390848214285714_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.01171875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1t0(1,1,7)=(0.765432059151786_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.319270121173469_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0483177794164541_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,1,7)=(0.0512273596938776_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.00913767538265306_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0381557451105442_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.005859375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1t0(1,2,8)=(-1.41375736607143_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.846850379464286_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.3331903125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,2,8)=(0.046875_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.09375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,1,9)=(-22.8134614118304_pr*c1*gA2)/(fpi4*mpi2)- &
         (1.21977883211097_pr*c3*gA2)/(fpi4*mpi2)- &
         (7.25512371332908_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,1,9)=(7.09276506696429_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.199196428571429_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2.32029482886905_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(3,1,9)=(-0.5009765625_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0966796875_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.2724609375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,3,9)=(-1.89848772321429_pr*c3*gA2)/(fpi4*mpi2)- &
         (1.26565848214286_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,2,10)=(45.7770589285714_pr*c1*gA2)/(fpi4*mpi2)- &
         (37.9879278300383_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.885098804209184_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,2,10)=(-5.06216517857143_pr*c1*gA2)/(fpi4*mpi2)+ &
         (9.70597098214286_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2.07220982142857_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,1,11)=(-440.874152071708_pr*c1*gA2)/(fpi4*mpi2)- &
         (597.813959665179_pr*c3*gA2)/(fpi4*mpi2)- &
         (216.912281412528_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,1,11)=(137.620562360491_pr*c1*gA2)/(fpi4*mpi2)+ &
         (242.982874720982_pr*c3*gA2)/(fpi4*mpi2)+ &
         (85.0374484514509_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(3,1,11)=(-12.2751813616071_pr*c1*gA2)/(fpi4*mpi2)- &
         (25.6883108956473_pr*c3*gA2)/(fpi4*mpi2)- &
         (9.29322160993304_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,3,11)=(-22.63640625_pr*c1*gA2)/(fpi4*mpi2)+ &
         (46.7609765625_pr*c3*gA2)/(fpi4*mpi2)+ &
         (6.892734375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,2,12)=(1465.24899274554_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1684.896046875_pr*c3*gA2)/(fpi4*mpi2)+ &
         (599.599941964286_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,2,12)=(-235.90740234375_pr*c1*gA2)/(fpi4*mpi2)- &
         (351.521103515625_pr*c3*gA2)/(fpi4*mpi2)- &
         (118.98005859375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,1,13)=(2941.54591061489_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2355.33459348319_pr*c3*gA2)/(fpi4*mpi2)+ &
         (805.315754579381_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,1,13)=(-1813.61831472866_pr*c1*gA2)/(fpi4*mpi2)- &
         (1619.87784794992_pr*c3*gA2)/(fpi4*mpi2)- &
         (563.160752122628_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(3,1,13)=(221.161743164063_pr*c1*gA2)/(fpi4*mpi2)+ &
         (225.247192382813_pr*c3*gA2)/(fpi4*mpi2)+ &
         (77.7313232421875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,3,13)=(-1171.69059375_pr*c1*gA2)/(fpi4*mpi2)- &
         (1176.1601625_pr*c3*gA2)/(fpi4*mpi2)- &
         (411.17551875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,2,14)=(-6502.52069384766_pr*c1*gA2)/(fpi4*mpi2)- &
         (5016.80990126953_pr*c3*gA2)/(fpi4*mpi2)- &
         (1655.41310106027_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,2,14)=(2369.16757114955_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1991.91343610491_pr*c3*gA2)/(fpi4*mpi2)+ &
         (672.799001116071_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,1,15)=(-3755.98198127747_pr*c1*gA2)/(fpi4*mpi2)- &
         (2756.34979092146_pr*c3*gA2)/(fpi4*mpi2)- &
         (783.247644786072_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,1,15)=(3861.25969552264_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2865.95377119141_pr*c3*gA2)/(fpi4*mpi2)+ &
         (867.669720144653_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(3,1,15)=(-880.419287109375_pr*c1*gA2)/(fpi4*mpi2)- &
         (683.312782287598_pr*c3*gA2)/(fpi4*mpi2)- &
         (216.528129577637_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,3,15)=(3099.20776171875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2319.50358105469_pr*c3*gA2)/(fpi4*mpi2)+ &
         (745.568033203125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,2,16)=(5164.53977966309_pr*c1*gA2)/(fpi4*mpi2)+ &
         (3750.44797599923_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1042.07968925432_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,2,16)=(-2919.87245544434_pr*c1*gA2)/(fpi4*mpi2)- &
         (2116.93618456159_pr*c3*gA2)/(fpi4*mpi2)- &
         (628.253677891323_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,1,17)=(113.823727294922_pr*c1*gA2)/(fpi4*mpi2)+ &
         (71.9582676125663_pr*c3*gA2)/(fpi4*mpi2)+ &
         (16.4675249980382_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,1,17)=(-793.151098434448_pr*c1*gA2)/(fpi4*mpi2)- &
         (552.996386267526_pr*c3*gA2)/(fpi4*mpi2)- &
         (142.309278507778_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(3,1,17)=(364.767998526437_pr*c1*gA2)/(fpi4*mpi2)+ &
         (251.980659746443_pr*c3*gA2)/(fpi4*mpi2)+ &
         (69.7415504673549_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,3,17)=(-1742.07032226563_pr*c1*gA2)/(fpi4*mpi2)- &
         (1259.881028268_pr*c3*gA2)/(fpi4*mpi2)- &
         (342.195423235212_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,2,18)=(-76.2155765991211_pr*c1*gA2)/(fpi4*mpi2)- &
         (48.9953290100098_pr*c3*gA2)/(fpi4*mpi2)- &
         (10.8880724487305_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,2,18)=(489.043165649414_pr*c1*gA2)/(fpi4*mpi2)+ &
         (346.766698773193_pr*c3*gA2)/(fpi4*mpi2)+ &
         (87.1195371459961_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,1,19)=(-0.68049621963501_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.481203679112026_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0826331377846854_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,1,19)=(-18.1235791282654_pr*c1*gA2)/(fpi4*mpi2)- &
         (11.0002656221662_pr*c3*gA2)/(fpi4*mpi2)- &
         (2.70709313686916_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(3,1,19)=(-21.7049503326416_pr*c1*gA2)/(fpi4*mpi2)- &
         (16.8184317398071_pr*c3*gA2)/(fpi4*mpi2)- &
         (3.76123477935791_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,2,20)=(19.0538941497803_pr*c1*gA2)/(fpi4*mpi2)+ &
         (12.2488322525024_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2.72201811218262_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,1,21)=(0.340248109817505_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.240601839556013_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0413165688923427_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(3,1,21)=(-2.58308817386627_pr*c1*gA2)/(fpi4*mpi2)- &
         (1.74732532024384_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.352447028160095_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(3,1,23)=(-0.0425310137271881_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0300752299445016_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00516457111154284_pr*c4*gA2)/(fpi4*mpi2)
    !
    !  Urho1tau1
    !
    ctr1t1=0.0_pr
    ctr1t1(1,1,3)=(-0.450116223436374_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0319572066731454_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,5)=(-3.5385371069221_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2.50259509338973_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.222006694302721_pr*c4*gA2)/(fpi4*mpi2)+ &
         (0.046875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1t1(2,1,5)=(0.0898447332383019_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0781696428571429_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.0234375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1t1(1,1,7)=(-17.2434756863839_pr*c1*gA2)/(fpi4*mpi2)+ &
         (8.60257873066855_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0966355588329082_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,7)=(19.5159960183232_pr*c1*gA2)/(fpi4*mpi2)- &
         (9.34811780250883_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0763114902210884_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.01171875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1t1(1,2,8)=(-168.933991071429_pr*c1*gA2)/(fpi4*mpi2)+ &
         (37.320722537714_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.666380625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,2,8)=(18.3753125054633_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.1875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,9)=(-3146.36888298668_pr*c1*gA2)/(fpi4*mpi2)- &
         (492.901337694334_pr*c3*gA2)/(fpi4*mpi2)- &
         (14.5102474266582_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,9)=(1205.47618769832_pr*c1*gA2)/(fpi4*mpi2)+ &
         (155.024983611535_pr*c3*gA2)/(fpi4*mpi2)+ &
         (4.6405896577381_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,9)=(-109.918476585036_pr*c1*gA2)/(fpi4*mpi2)- &
         (32.7311718821706_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.544921875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,3,9)=(-274.364029017857_pr*c3*gA2)/(fpi4*mpi2)- &
         (2.53131696428571_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,2,10)=(5901.79117633929_pr*c1*gA2)/(fpi4*mpi2)- &
         (385.542873162743_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1.77019760841837_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,2,10)=(-404.58984375_pr*c1*gA2)/(fpi4*mpi2)+ &
         (690.035446485936_pr*c3*gA2)/(fpi4*mpi2)+ &
         (4.14441964285714_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,11)=(-12146.9163811618_pr*c1*gA2)/(fpi4*mpi2)- &
         (9992.47844580492_pr*c3*gA2)/(fpi4*mpi2)- &
         (433.824562825056_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,11)=(4424.78103473708_pr*c1*gA2)/(fpi4*mpi2)+ &
         (5157.68192548934_pr*c3*gA2)/(fpi4*mpi2)+ &
         (170.074896902902_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,11)=(-968.792783439754_pr*c1*gA2)/(fpi4*mpi2)- &
         (956.662918816572_pr*c3*gA2)/(fpi4*mpi2)- &
         (18.5864432198661_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,3,11)=(-7552.959375_pr*c1*gA2)/(fpi4*mpi2)- &
         (1191.78466145833_pr*c3*gA2)/(fpi4*mpi2)+ &
         (13.78546875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,2,12)=(2094.44863895089_pr*c1*gA2)/(fpi4*mpi2)+ &
         (11355.6974476442_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1199.19988392857_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,2,12)=(9461.258671875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (113.574869005018_pr*c3*gA2)/(fpi4*mpi2)- &
         (237.9601171875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,13)=(-52317.8224267498_pr*c1*gA2)/(fpi4*mpi2)- &
         (15039.7248405472_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1610.63150915876_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,13)=(35095.2865087366_pr*c1*gA2)/(fpi4*mpi2)+ &
         (8132.24487638635_pr*c3*gA2)/(fpi4*mpi2)- &
         (1126.32150424526_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,13)=(-9282.48454203116_pr*c1*gA2)/(fpi4*mpi2)- &
         (2450.24959567222_pr*c3*gA2)/(fpi4*mpi2)+ &
         (155.462646484375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,3,13)=(-10352.0559375_pr*c1*gA2)/(fpi4*mpi2)- &
         (9070.84456875_pr*c3*gA2)/(fpi4*mpi2)- &
         (822.3510375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,2,14)=(39347.6453808594_pr*c1*gA2)/(fpi4*mpi2)+ &
         (7581.1506071644_pr*c3*gA2)/(fpi4*mpi2)- &
         (3310.82620212054_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,2,14)=(4830.05495256696_pr*c1*gA2)/(fpi4*mpi2)+ &
         (6468.34894359579_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1345.59800223214_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,15)=(-41054.5351555039_pr*c1*gA2)/(fpi4*mpi2)- &
         (16289.8093013422_pr*c3*gA2)/(fpi4*mpi2)- &
         (1566.49528957214_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,15)=(35934.130351878_pr*c1*gA2)/(fpi4*mpi2)+ &
         (15494.0531940275_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1735.33944028931_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,15)=(-13686.0308464825_pr*c1*gA2)/(fpi4*mpi2)- &
         (6029.88917716322_pr*c3*gA2)/(fpi4*mpi2)- &
         (433.056259155273_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,3,15)=(-25601.809921875_pr*c1*gA2)/(fpi4*mpi2)- &
         (7414.09085742188_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1491.13606640625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,2,16)=(6185.05237060547_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2022.39299342052_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2084.15937850865_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,2,16)=(16757.6425048828_pr*c1*gA2)/(fpi4*mpi2)+ &
         (5974.47035520464_pr*c3*gA2)/(fpi4*mpi2)- &
         (1256.50735578265_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,17)=(-28271.1185486804_pr*c1*gA2)/(fpi4*mpi2)- &
         (11492.2038896716_pr*c3*gA2)/(fpi4*mpi2)+ &
         (32.9350499960763_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,17)=(40512.7185386068_pr*c1*gA2)/(fpi4*mpi2)+ &
         (16490.1542518688_pr*c3*gA2)/(fpi4*mpi2)- &
         (284.618557015555_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,17)=(-17321.2311050738_pr*c1*gA2)/(fpi4*mpi2)- &
         (6918.84409512217_pr*c3*gA2)/(fpi4*mpi2)+ &
         (139.48310093471_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,3,17)=(-6968.2812890625_pr*c1*gA2)/(fpi4*mpi2)- &
         (2519.76205653599_pr*c3*gA2)/(fpi4*mpi2)- &
         (684.390846470424_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,2,18)=(-1680.75223754883_pr*c1*gA2)/(fpi4*mpi2)- &
         (638.933174219611_pr*c3*gA2)/(fpi4*mpi2)- &
         (21.7761448974609_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,2,18)=(6262.39431518555_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2333.63039529067_pr*c3*gA2)/(fpi4*mpi2)+ &
         (174.239074291992_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,19)=(-14620.5213801411_pr*c1*gA2)/(fpi4*mpi2)- &
         (6001.88075314032_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.165266275569371_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,19)=(27371.0534029225_pr*c1*gA2)/(fpi4*mpi2)+ &
         (11145.7119069224_pr*c3*gA2)/(fpi4*mpi2)- &
         (5.41418627373832_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,19)=(-12338.4374338396_pr*c1*gA2)/(fpi4*mpi2)- &
         (4963.55634935681_pr*c3*gA2)/(fpi4*mpi2)- &
         (7.52246955871582_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,2,20)=(420.188059387207_pr*c1*gA2)/(fpi4*mpi2)+ &
         (159.733293554903_pr*c3*gA2)/(fpi4*mpi2)+ &
         (5.44403622436523_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,21)=(-4953.03311039961_pr*c1*gA2)/(fpi4*mpi2)- &
         (2058.14892578187_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,21)=(11682.2041357371_pr*c1*gA2)/(fpi4*mpi2)+ &
         (4815.90379275659_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0826331377846854_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,21)=(-6036.90320983232_pr*c1*gA2)/(fpi4*mpi2)- &
         (2466.73717501561_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.70489405632019_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,23)=(-1148.49749408445_pr*c1*gA2)/(fpi4*mpi2)- &
         (481.334985976422_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,23)=(3528.18980908462_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1469.64364495069_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,23)=(-2246.89838063988_pr*c1*gA2)/(fpi4*mpi2)- &
         (929.407234212685_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0103291422230857_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,25)=(-172.469986869063_pr*c1*gA2)/(fpi4*mpi2)- &
         (72.6337448118083_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,25)=(737.432931746486_pr*c1*gA2)/(fpi4*mpi2)+ &
         (309.385467866635_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,25)=(-611.022455581544_pr*c1*gA2)/(fpi4*mpi2)- &
         (255.004474494164_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,27)=(-15.0982906061879_pr*c1*gA2)/(fpi4*mpi2)- &
         (6.36541276884548_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,27)=(100.943421587724_pr*c1*gA2)/(fpi4*mpi2)+ &
         (42.51836588516_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,27)=(-116.156780981393_pr*c1*gA2)/(fpi4*mpi2)- &
         (48.7723238905489_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,29)=(-0.584793679492634_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.245878934384381_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,29)=(8.13393898258656_pr*c1*gA2)/(fpi4*mpi2)+ &
         (3.42858531880712_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,29)=(-14.6026796374878_pr*c1*gA2)/(fpi4*mpi2)- &
         (6.15145215414809_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,31)=(0.292396839746317_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.12293946719219_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,31)=(-1.0898415827599_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.459308031648938_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,33)=(-0.0365496049682896_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0153674333990238_pr*c3*gA2)/(fpi4*mpi2)
    !
    !  UJ0Drho0
    !
    ctJ0dr0=0.0_pr
    !
    !  Urho0DJ0
    !
    ctr0DJ0=0.0_pr
    ctr0DJ0(1,1,3)=(-0.244419766459837_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0234375000151846_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,5)=(-1.03158597015212_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.780470479470306_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0351562492028061_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,7)=(-21.0022652575419_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2.27019651550668_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00439427564971301_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,7)=(6.64452901902785_pr*c1*gA2)/(fpi4*mpi2)- &
         (1.38867403439195_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0175780851403061_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,2,8)=(-0.000291666666666667_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0000478125_pr*c3*gA2)/(fpi4*mpi2)- &
         (4.01785714285714e-7_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,9)=(-307.281664803321_pr*c1*gA2)/(fpi4*mpi2)- &
         (17.4248922644585_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0131892671795281_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,9)=(167.072976616587_pr*c1*gA2)/(fpi4*mpi2)+ &
         (4.91318494493787_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00878976004464286_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,9)=(-20.6458593811462_pr*c1*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,2,10)=(0.07489375_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.000520145089285714_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0000100207270408163_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,2,10)=(-0.009375_pr*c1*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,11)=(-1417.24568001336_pr*c1*gA2)/(fpi4*mpi2)- &
         (148.086566637738_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00426528076171875_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,11)=(1006.30895614752_pr*c1*gA2)/(fpi4*mpi2)+ &
         (93.3959309491168_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0102003609793527_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,11)=(-165.352324283285_pr*c1*gA2)/(fpi4*mpi2)- &
         (12.9178818130112_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00332205636160714_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,3,11)=(-0.0875_pr*c1*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,2,12)=(-0.663125_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.09263671875_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00322935267857143_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,2,12)=(0.22640625_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.01634765625_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00052734375_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,13)=(-3324.31750521588_pr*c1*gA2)/(fpi4*mpi2)- &
         (384.460747217616_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00614578334263393_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,13)=(2946.4643595763_pr*c1*gA2)/(fpi4*mpi2)+ &
         (332.165470161666_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00216459803989955_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,13)=(-602.212375765244_pr*c1*gA2)/(fpi4*mpi2)- &
         (64.5897786600606_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.001190185546875_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,3,13)=(0.37125_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0631546875_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0026578125_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,2,14)=(1.72600911458333_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.135301904296875_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0143205845424107_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,2,14)=(-0.62890625_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0793212890625_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00530482700892857_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,15)=(-4612.5195280616_pr*c1*gA2)/(fpi4*mpi2)- &
         (513.154765943116_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00963807367597307_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,15)=(5043.02746283907_pr*c1*gA2)/(fpi4*mpi2)+ &
         (573.956470083127_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00939304275512695_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,15)=(-1259.59323184652_pr*c1*gA2)/(fpi4*mpi2)- &
         (143.759310262187_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00231170654296875_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,3,15)=(-0.71465625_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0487265625_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0069609375_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,2,16)=(-0.13190234375_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.084516064453125_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0133438685825893_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,2,16)=(0.345556640625_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0366448974609375_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0070806884765625_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,17)=(-4081.8063003269_pr*c1*gA2)/(fpi4*mpi2)- &
         (405.584681811802_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.000288566207885742_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,17)=(5511.17933250622_pr*c1*gA2)/(fpi4*mpi2)+ &
         (582.019749761264_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00207645034790039_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,17)=(-1676.42040558742_pr*c1*gA2)/(fpi4*mpi2)- &
         (184.111758668806_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.000897445678710938_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,3,17)=(-0.031904296875_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00455775669642857_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,2,18)=(-0.0106455078125_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.00226593017578125_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00019940185546875_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,2,18)=(0.03829833984375_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0119082458496094_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00132149047851563_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,19)=(-2389.93945013977_pr*c1*gA2)/(fpi4*mpi2)- &
         (200.758956740934_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2.18095779418945e-6_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,19)=(4026.89040316843_pr*c1*gA2)/(fpi4*mpi2)+ &
         (373.930398299718_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0000424012184143066_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,19)=(-1499.34103770493_pr*c1*gA2)/(fpi4*mpi2)- &
         (149.680535490279_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0000683426856994629_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,2,20)=(0.002661376953125_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.000566482543945312_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0000498504638671875_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,21)=(-937.343180466548_pr*c1*gA2)/(fpi4*mpi2)- &
         (62.9319292613857_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,21)=(2012.17746809803_pr*c1*gA2)/(fpi4*mpi2)+ &
         (156.764429796522_pr*c3*gA2)/(fpi4*mpi2)- &
         (1.09047889709473e-6_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,21)=(-928.899230533649_pr*c1*gA2)/(fpi4*mpi2)- &
         (80.7120141836349_pr*c3*gA2)/(fpi4*mpi2)+ &
         (7.43508338928223e-6_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,23)=(-243.60536380909_pr*c1*gA2)/(fpi4*mpi2)- &
         (12.1417483739015_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,23)=(690.008776700697_pr*c1*gA2)/(fpi4*mpi2)+ &
         (42.8121424298444_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,23)=(-403.794812461118_pr*c1*gA2)/(fpi4*mpi2)- &
         (29.2881400946876_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1.36309862136841e-7_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,25)=(-40.2494438142034_pr*c1*gA2)/(fpi4*mpi2)- &
         (1.31603730589479_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,25)=(159.713027711956_pr*c1*gA2)/(fpi4*mpi2)+ &
         (7.34601767825556_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,25)=(-122.828777785615_pr*c1*gA2)/(fpi4*mpi2)- &
         (7.0783524478032_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,27)=(-3.82761492020625_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0613407218849266_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,27)=(23.8460141906355_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.719359374832321_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,27)=(-25.6066641016362_pr*c1*gA2)/(fpi4*mpi2)- &
         (1.09298032666628_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,29)=(-0.159483955008594_pr*c1*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,29)=(2.07329141511172_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0306703609424633_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,29)=(-3.48578429802332_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0975875120896559_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,31)=(0.0797419775042969_pr*c1*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,31)=(-0.279096921265039_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.00383379511780791_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,33)=(-0.00996774718803712_pr*c1*gA2)/(fpi4*mpi2)
    !
    !  UJ1Drho1
    !
    ctJ1dr1=0.0_pr
    !
    !  Urho1DJ1
    !
    ctr1dJ1=0.0_pr
    ctr1dJ1(1,1,3)=(0.0736607554815508_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0156250000101231_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,5)=(0.308705740050707_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.248438076755833_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0234374994685374_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,7)=(6.96559883584731_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.755267413285654_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00292951709980867_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,7)=(-2.17968675634262_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.45703198308388_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0117187234268707_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,2,8)=(0.0000972222222222222_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0000158035714285714_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2.67857142857143e-7_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,9)=(102.418432538607_pr*c1*gA2)/(fpi4*mpi2)+ &
         (5.80390099909298_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00879284478635204_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,9)=(-55.6646250180288_pr*c1*gA2)/(fpi4*mpi2)- &
         (1.63479839496441_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0058598400297619_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,9)=(6.87316406454873_pr*c1*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,2,10)=(-0.0249645833333333_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.000170041454081633_pr*c3*gA2)/(fpi4*mpi2)+ &
         (6.68048469387755e-6_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,2,10)=(0.003125_pr*c1*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,11)=(472.415226671119_pr*c1*gA2)/(fpi4*mpi2)+ &
         (49.3607671189922_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0028435205078125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,11)=(-335.431924184589_pr*c1*gA2)/(fpi4*mpi2)- &
         (31.1285768627125_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00680024065290179_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,11)=(55.1130468965117_pr*c1*gA2)/(fpi4*mpi2)+ &
         (4.30485325221653_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00221470424107143_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,3,11)=(0.0291666666666667_pr*c1*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,2,12)=(0.221041666666667_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0298024553571429_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00215290178571429_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,2,12)=(-0.07546875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0052734375_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0003515625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,13)=(1108.10583507196_pr*c1*gA2)/(fpi4*mpi2)+ &
         (128.155631000319_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00409718889508929_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,13)=(-982.154786525432_pr*c1*gA2)/(fpi4*mpi2)- &
         (110.722544919902_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00144306535993304_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,13)=(200.736909272009_pr*c1*gA2)/(fpi4*mpi2)+ &
         (21.5295294915046_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00079345703125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,3,13)=(-0.12375_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.020165625_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.001771875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,2,14)=(-0.575336371527778_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0403271065848214_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00954705636160714_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,2,14)=(0.209635416666667_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0246721540178571_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00353655133928571_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,15)=(1537.50650935387_pr*c1*gA2)/(fpi4*mpi2)+ &
         (171.04837595648_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00642538245064872_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,15)=(-1681.00915427969_pr*c1*gA2)/(fpi4*mpi2)- &
         (191.31569234679_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00626202850341797_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,15)=(419.864410615508_pr*c1*gA2)/(fpi4*mpi2)+ &
         (47.9189995185479_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0015411376953125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,3,15)=(0.23821875_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.013921875_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.004640625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,2,16)=(0.0439674479166667_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0237240652901786_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00889591238839286_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,2,16)=(-0.115185546875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.009854736328125_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.004720458984375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,17)=(1360.60210010897_pr*c1*gA2)/(fpi4*mpi2)+ &
         (135.194990126003_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.000192377471923828_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,17)=(-1837.05977750207_pr*c1*gA2)/(fpi4*mpi2)- &
         (194.007275403871_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00138430023193359_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,17)=(558.806801862472_pr*c1*gA2)/(fpi4*mpi2)+ &
         (61.370885371495_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.000598297119140625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,3,17)=(0.00911551339285714_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00303850446428571_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,2,18)=(0.00354850260416667_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0006888427734375_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0001329345703125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,2,18)=(-0.01276611328125_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.00352891845703125_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00088099365234375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,19)=(796.646483379922_pr*c1*gA2)/(fpi4*mpi2)+ &
         (66.919651519992_pr*c3*gA2)/(fpi4*mpi2)- &
         (1.45397186279297e-6_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,19)=(-1342.29680105614_pr*c1*gA2)/(fpi4*mpi2)- &
         (124.643480233646_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0000282674789428711_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,19)=(499.780345901644_pr*c1*gA2)/(fpi4*mpi2)+ &
         (49.8934890491979_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0000455617904663086_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,2,20)=(-0.000887125651041667_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.000172210693359375_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.000033233642578125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,21)=(312.447726822183_pr*c1*gA2)/(fpi4*mpi2)+ &
         (20.9773097537952_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,21)=(-670.725822699343_pr*c1*gA2)/(fpi4*mpi2)- &
         (52.2548095686812_pr*c3*gA2)/(fpi4*mpi2)+ &
         (7.26985931396484e-7_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,21)=(309.63307684455_pr*c1*gA2)/(fpi4*mpi2)+ &
         (26.9040022495172_pr*c3*gA2)/(fpi4*mpi2)- &
         (4.95672225952148e-6_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,23)=(81.2017879363635_pr*c1*gA2)/(fpi4*mpi2)+ &
         (4.04724945796715_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,23)=(-230.002925566899_pr*c1*gA2)/(fpi4*mpi2)- &
         (14.2707141432815_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,23)=(134.598270820373_pr*c1*gA2)/(fpi4*mpi2)+ &
         (9.76271331945924_pr*c3*gA2)/(fpi4*mpi2)- &
         (9.08732414245605e-8_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,25)=(13.4164812714011_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.43867910196493_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,25)=(-53.2376759039853_pr*c1*gA2)/(fpi4*mpi2)- &
         (2.44867255941852_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,25)=(40.9429259285383_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2.3594508159344_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,27)=(1.27587164006875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0204469072949755_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,27)=(-7.94867139687851_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.23978645827744_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,27)=(8.53555470054539_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.364326775555427_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,29)=(0.053161318336198_pr*c1*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,29)=(-0.691097138370574_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0102234536474878_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,29)=(1.16192809934111_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.032529170696552_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,31)=(-0.026580659168099_pr*c1*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,31)=(0.0930323070883464_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.00127793170593597_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,33)=(0.00332258239601237_pr*c1*gA2)/(fpi4*mpi2)
    !
    !  UJ0Drho1
    !
    ctJ0dr1=0.0_pr
    !
    !  Urho1DJ0
    !
    ctr1dJ0=0.0_pr
    ctr1dJ0(1,1,3)=(0.0234375000151846_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00781250000506155_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,1,5)=(0.10546875_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0351562492028061_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0117187497342687_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,1,7)=(0.10546875_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.00439427564971301_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00146475854990434_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,1,7)=(-0.10546875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0175780851403061_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00585936171343537_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,2,8)=(-4.01785714285714e-7_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1.33928571428571e-7_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,1,9)=(0.0263671875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0131892671795281_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00439642239317602_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,1,9)=(-0.0791015625_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.00878976004464286_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00292992001488095_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(3,1,9)=(0.0263671875_pr*c1*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,2,10)=(-0.0000100207270408163_pr*c3*gA2)/(fpi4*mpi2)+ &
         (3.34024234693878e-6_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,1,11)=(0.00426528076171875_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00142176025390625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,1,11)=(-0.01318359375_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0102003609793527_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00340012032645089_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(3,1,11)=(0.01318359375_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.00332205636160714_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00110735212053571_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,2,12)=(-0.00322935267857143_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00107645089285714_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,2,12)=(0.00052734375_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00017578125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,1,13)=(-0.00614578334263393_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00204859444754464_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,1,13)=(0.00216459803989955_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.000721532679966518_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(3,1,13)=(0.00164794921875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.001190185546875_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.000396728515625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,3,13)=(0.0026578125_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0008859375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,2,14)=(0.0143205845424107_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00477352818080357_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,2,14)=(-0.00530482700892857_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00176827566964286_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,1,15)=(0.00963807367597307_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00321269122532436_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,1,15)=(-0.00939304275512695_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00313101425170898_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(3,1,15)=(0.00231170654296875_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00077056884765625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,3,15)=(-0.0069609375_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0023203125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,2,16)=(-0.0133438685825893_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00444795619419643_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,2,16)=(0.0070806884765625_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0023602294921875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,1,17)=(-0.000288566207885742_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0000961887359619141_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,1,17)=(0.00207645034790039_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.000692150115966797_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(3,1,17)=(-0.000897445678710938_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.000299148559570312_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,3,17)=(0.00455775669642857_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00151925223214286_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,2,18)=(0.00019940185546875_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00006646728515625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,2,18)=(-0.00132149047851563_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.000440496826171875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,1,19)=(2.18095779418945e-6_pr*c3*gA2)/(fpi4*mpi2)- &
         (7.26985931396484e-7_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,1,19)=(0.0000424012184143066_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0000141337394714355_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(3,1,19)=(0.0000683426856994629_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0000227808952331543_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,2,20)=(-0.0000498504638671875_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0000166168212890625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,1,21)=(-1.09047889709473e-6_pr*c3*gA2)/(fpi4*mpi2)+ &
         (3.63492965698242e-7_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(3,1,21)=(7.43508338928223e-6_pr*c3*gA2)/(fpi4*mpi2)- &
         (2.47836112976074e-6_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(3,1,23)=(1.36309862136841e-7_pr*c3*gA2)/(fpi4*mpi2)- &
         (4.54366207122803e-8_pr*c4*gA2)/(fpi4*mpi2)
    !
    !  UJ1Drho0
    !
    ctJ1dr0=0.0_pr
    !
    !  UJ0J0
    !
    ctJ0J0=0.0_pr
    ctJ0J0(1,1,3)=(-0.0388151948706056_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0689519291560374_pr*c4*gA2)/(fpi4*mpi2)+ &
         (0.03515625_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ0J0(1,1,5)=(-0.148030393439646_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.125249890279611_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.124165184550383_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.017578125_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ0J0(1,2,6)=(0.266832652280755_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0685045280612245_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,7)=(-50.9400757481906_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1.49639517552504_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00598016955317283_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,7)=(9.4436896594477_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.387304024062572_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00211338089923469_pr*c4*gA2)/(fpi4*mpi2)+ &
         (0.00439453125_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ0J0(1,2,8)=(60.4392047247024_pr*c1*gA2)/(fpi4*mpi2)+ &
         (21.2791752723318_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0653900223214286_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,9)=(-6519.80035123157_pr*c1*gA2)/(fpi4*mpi2)+ &
         (316.045390495029_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2.16832575424506_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,9)=(2094.21472214372_pr*c1*gA2)/(fpi4*mpi2)- &
         (79.8751918140941_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.347175990513393_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,9)=(-174.923789099377_pr*c1*gA2)/(fpi4*mpi2)
    ctJ0J0(1,3,9)=(-136.232770647321_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.632829241071429_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,2,10)=(22280.2381832589_pr*c1*gA2)/(fpi4*mpi2)- &
         (3725.75432559614_pr*c3*gA2)/(fpi4*mpi2)- &
         (40.7846229073661_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,2,10)=(-3482.03035714286_pr*c1*gA2)/(fpi4*mpi2)+ &
         (857.700259443584_pr*c3*gA2)/(fpi4*mpi2)+ &
         (7.79614955357143_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,11)=(52284.4157057857_pr*c1*gA2)/(fpi4*mpi2)- &
         (22649.0853151035_pr*c3*gA2)/(fpi4*mpi2)- &
         (389.622252285779_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,11)=(-29244.8802258893_pr*c1*gA2)/(fpi4*mpi2)+ &
         (11240.6895935096_pr*c3*gA2)/(fpi4*mpi2)+ &
         (170.39136171439_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,11)=(3224.34088479694_pr*c1*gA2)/(fpi4*mpi2)- &
         (1377.64430659016_pr*c3*gA2)/(fpi4*mpi2)- &
         (18.1685616629464_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,3,11)=(-20240.849921875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (3674.53847265625_pr*c3*gA2)/(fpi4*mpi2)+ &
         (49.9097109375_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,2,12)=(-142538.810602121_pr*c1*gA2)/(fpi4*mpi2)+ &
         (51959.7177006278_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1056.48083035714_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,2,12)=(47556.1071386719_pr*c1*gA2)/(fpi4*mpi2)- &
         (13973.2343116668_pr*c3*gA2)/(fpi4*mpi2)- &
         (251.396513671875_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,13)=(-141614.511333235_pr*c1*gA2)/(fpi4*mpi2)+ &
         (35973.0682685409_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1263.90306398228_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,13)=(120423.751289578_pr*c1*gA2)/(fpi4*mpi2)- &
         (33366.8468838635_pr*c3*gA2)/(fpi4*mpi2)- &
         (938.242939340534_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,13)=(-23741.2058036637_pr*c1*gA2)/(fpi4*mpi2)+ &
         (6070.38539035745_pr*c3*gA2)/(fpi4*mpi2)+ &
         (145.681594848633_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,3,13)=(76403.207390625_pr*c1*gA2)/(fpi4*mpi2)- &
         (27942.3675271205_pr*c3*gA2)/(fpi4*mpi2)- &
         (671.186349441964_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,2,14)=(204907.907447184_pr*c1*gA2)/(fpi4*mpi2)- &
         (56211.484659242_pr*c3*gA2)/(fpi4*mpi2)- &
         (2422.99460657436_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,2,14)=(-95398.9531424386_pr*c1*gA2)/(fpi4*mpi2)+ &
         (30129.2960691746_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1018.34175868443_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,15)=(9978.05620114632_pr*c1*gA2)/(fpi4*mpi2)- &
         (6474.12892378917_pr*c3*gA2)/(fpi4*mpi2)- &
         (968.323984416389_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,15)=(-36487.3307862002_pr*c1*gA2)/(fpi4*mpi2)+ &
         (13937.926974747_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1128.44835742705_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,15)=(13562.8790229223_pr*c1*gA2)/(fpi4*mpi2)- &
         (5370.52492057074_pr*c3*gA2)/(fpi4*mpi2)- &
         (292.012623596191_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,3,15)=(-76061.0427714844_pr*c1*gA2)/(fpi4*mpi2)+ &
         (21292.048190918_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1052.55900878906_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,2,16)=(-26755.8098873901_pr*c1*gA2)/(fpi4*mpi2)+ &
         (11051.0158225074_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1284.3833242048_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,2,16)=(33104.1130819702_pr*c1*gA2)/(fpi4*mpi2)- &
         (10790.492615363_pr*c3*gA2)/(fpi4*mpi2)- &
         (813.688233476911_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,17)=(-13525.0298722901_pr*c1*gA2)/(fpi4*mpi2)+ &
         (3920.12706464951_pr*c3*gA2)/(fpi4*mpi2)+ &
         (21.2402742906025_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,17)=(18744.3022469218_pr*c1*gA2)/(fpi4*mpi2)- &
         (5272.85370859068_pr*c3*gA2)/(fpi4*mpi2)- &
         (174.676497182301_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,17)=(-7131.99279047117_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1930.52177440425_pr*c3*gA2)/(fpi4*mpi2)+ &
         (88.4860045787266_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,3,17)=(7839.31645019531_pr*c1*gA2)/(fpi4*mpi2)- &
         (2869.73837755476_pr*c3*gA2)/(fpi4*mpi2)- &
         (419.970217895508_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,2,18)=(355.059586151123_pr*c1*gA2)/(fpi4*mpi2)+ &
         (183.100210884983_pr*c3*gA2)/(fpi4*mpi2)- &
         (13.999078478132_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,2,18)=(-1327.89377142334_pr*c1*gA2)/(fpi4*mpi2)+ &
         (15.4343164854196_pr*c3*gA2)/(fpi4*mpi2)+ &
         (105.873926083374_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,19)=(-8867.21170169755_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2936.70582294006_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.109368735306604_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,19)=(13819.8747331675_pr*c1*gA2)/(fpi4*mpi2)- &
         (4438.03614592496_pr*c3*gA2)/(fpi4*mpi2)- &
         (3.51122917092868_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,19)=(-4725.00077792637_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1532.04177968309_pr*c3*gA2)/(fpi4*mpi2)- &
         (4.33029295563698_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,2,20)=(-88.7648965377808_pr*c1*gA2)/(fpi4*mpi2)- &
         (45.7750527212457_pr*c3*gA2)/(fpi4*mpi2)+ &
         (3.49976961953299_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,21)=(-3950.27983470918_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1347.42215988918_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,21)=(7835.22431087985_pr*c1*gA2)/(fpi4*mpi2)- &
         (2624.64064001171_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0546843676533018_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,21)=(-3337.91716486204_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1111.14574092555_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.449709850430489_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,23)=(-1142.83677022112_pr*c1*gA2)/(fpi4*mpi2)+ &
         (401.209483659988_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,23)=(3005.12399409001_pr*c1*gA2)/(fpi4*mpi2)- &
         (1034.46527332476_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,23)=(-1635.837501786_pr*c1*gA2)/(fpi4*mpi2)+ &
         (553.405150406837_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00683554595666272_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,25)=(-207.173615193637_pr*c1*gA2)/(fpi4*mpi2)+ &
         (74.5943777265557_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,25)=(765.620232650329_pr*c1*gA2)/(fpi4*mpi2)- &
         (270.43438082633_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,25)=(-549.802354657625_pr*c1*gA2)/(fpi4*mpi2)+ &
         (190.708530684211_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,27)=(-21.371516474291_pr*c1*gA2)/(fpi4*mpi2)+ &
         (7.86480912750357_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,27)=(124.32036907328_pr*c1*gA2)/(fpi4*mpi2)- &
         (44.92276431339_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,27)=(-125.00166164342_pr*c1*gA2)/(fpi4*mpi2)+ &
         (44.3795879210135_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,29)=(-0.956932496743248_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.358850516087132_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,29)=(11.6426907338887_pr*c1*gA2)/(fpi4*mpi2)- &
         (4.29125507983892_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,29)=(-18.3709744429036_pr*c1*gA2)/(fpi4*mpi2)+ &
         (6.65825509945955_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,31)=(0.478466248371624_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.179425258043566_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,31)=(-1.574952903829_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.581263199490756_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,33)=(-0.059808281046453_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0224281572554458_pr*c3*gA2)/(fpi4*mpi2)
    !
    !  UJ0J1
    !
    ctJ0J1=0.0_pr
    ctJ0J1(1,1,3)=(0.0285278380152656_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0459679527706916_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.0234375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ0J1(1,1,5)=(-0.00310602359693878_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0516361625744048_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0827767897002551_pr*c4*gA2)/(fpi4*mpi2)+ &
         (0.01171875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ0J1(1,2,6)=(0.0475446853741497_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0456696853741497_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,1,7)=(0.147463309151786_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0358081265611713_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00398677970211522_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,1,7)=(0.00200860969387755_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0111190874787415_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0014089205994898_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.0029296875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ0J1(1,2,8)=(-0.363757366071429_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.177937700892857_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0435933482142857_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,1,9)=(11.9185698381696_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.18312424087213_pr*c3*gA2)/(fpi4*mpi2)- &
         (1.44555050283004_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,1,9)=(-2.81250837053571_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0791064453125_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.231450660342262_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(3,1,9)=(0.17578125_pr*c1*gA2)/(fpi4*mpi2)
    ctJ0J1(1,3,9)=(-0.210943080357143_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.421886160714286_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,2,10)=(-90.1323160714286_pr*c1*gA2)/(fpi4*mpi2)+ &
         (9.21932595663265_pr*c3*gA2)/(fpi4*mpi2)+ &
         (27.1897486049107_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,2,10)=(12.6565848214286_pr*c1*gA2)/(fpi4*mpi2)- &
         (1.095703125_pr*c3*gA2)/(fpi4*mpi2)- &
         (5.19743303571429_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,1,11)=(-798.881085665458_pr*c1*gA2)/(fpi4*mpi2)+ &
         (160.074257447634_pr*c3*gA2)/(fpi4*mpi2)+ &
         (259.748168190519_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,1,11)=(328.414995954241_pr*c1*gA2)/(fpi4*mpi2)- &
         (61.5935256713867_pr*c3*gA2)/(fpi4*mpi2)- &
         (113.594241142927_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(3,1,11)=(-30.7629743303571_pr*c1*gA2)/(fpi4*mpi2)+ &
         (5.57625906808036_pr*c3*gA2)/(fpi4*mpi2)+ &
         (12.1123744419643_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,3,11)=(123.05109375_pr*c1*gA2)/(fpi4*mpi2)- &
         (15.2580234375_pr*c3*gA2)/(fpi4*mpi2)- &
         (33.273140625_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,2,12)=(2557.41305524554_pr*c1*gA2)/(fpi4*mpi2)- &
         (494.193056919643_pr*c3*gA2)/(fpi4*mpi2)- &
         (704.320553571429_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,2,12)=(-600.61833984375_pr*c1*gA2)/(fpi4*mpi2)+ &
         (105.762021484375_pr*c3*gA2)/(fpi4*mpi2)+ &
         (167.59767578125_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,1,13)=(4226.57441403285_pr*c1*gA2)/(fpi4*mpi2)- &
         (759.654990867397_pr*c3*gA2)/(fpi4*mpi2)- &
         (842.602042654855_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,1,13)=(-2859.96228445522_pr*c1*gA2)/(fpi4*mpi2)+ &
         (519.521237971814_pr*c3*gA2)/(fpi4*mpi2)+ &
         (625.495292893689_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(3,1,13)=(423.723815917969_pr*c1*gA2)/(fpi4*mpi2)- &
         (74.1513366699219_pr*c3*gA2)/(fpi4*mpi2)- &
         (97.1210632324219_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,3,13)=(-1778.06559375_pr*c1*gA2)/(fpi4*mpi2)+ &
         (349.290122544643_pr*c3*gA2)/(fpi4*mpi2)+ &
         (447.457566294643_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,2,14)=(-8387.96795947266_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1557.82783909738_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1615.32973771624_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,2,14)=(3204.65585239955_pr*c1*gA2)/(fpi4*mpi2)- &
         (610.780761369978_pr*c3*gA2)/(fpi4*mpi2)- &
         (678.894505789621_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,1,15)=(-3814.31004524231_pr*c1*gA2)/(fpi4*mpi2)+ &
         (782.002794581876_pr*c3*gA2)/(fpi4*mpi2)+ &
         (645.54932294426_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,1,15)=(4206.02934884295_pr*c1*gA2)/(fpi4*mpi2)- &
         (834.729188615308_pr*c3*gA2)/(fpi4*mpi2)- &
         (752.298904951368_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(3,1,15)=(-1006.82687988281_pr*c1*gA2)/(fpi4*mpi2)+ &
         (200.402084350586_pr*c3*gA2)/(fpi4*mpi2)+ &
         (194.675082397461_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,3,15)=(3781.37963671875_pr*c1*gA2)/(fpi4*mpi2)- &
         (709.003098632812_pr*c3*gA2)/(fpi4*mpi2)- &
         (701.706005859375_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,2,16)=(5207.17552185059_pr*c1*gA2)/(fpi4*mpi2)- &
         (1065.41645752128_pr*c3*gA2)/(fpi4*mpi2)- &
         (856.255549469866_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,2,16)=(-3152.00038513184_pr*c1*gA2)/(fpi4*mpi2)+ &
         (618.68313205392_pr*c3*gA2)/(fpi4*mpi2)+ &
         (542.458822317941_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,1,17)=(114.489910766602_pr*c1*gA2)/(fpi4*mpi2)- &
         (21.7326762047359_pr*c3*gA2)/(fpi4*mpi2)- &
         (14.1601828604017_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,1,17)=(-785.971121017456_pr*c1*gA2)/(fpi4*mpi2)+ &
         (158.614487102509_pr*c3*gA2)/(fpi4*mpi2)+ &
         (116.450998121534_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(3,1,17)=(379.629753897531_pr*c1*gA2)/(fpi4*mpi2)- &
         (73.6791959381104_pr*c3*gA2)/(fpi4*mpi2)- &
         (58.9906697191511_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,3,17)=(-1742.07032226563_pr*c1*gA2)/(fpi4*mpi2)+ &
         (357.745064784459_pr*c3*gA2)/(fpi4*mpi2)+ &
         (279.980145263672_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,2,18)=(-76.2155765991211_pr*c1*gA2)/(fpi4*mpi2)+ &
         (14.7764228733608_pr*c3*gA2)/(fpi4*mpi2)+ &
         (9.33271898542132_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,2,18)=(478.384230102539_pr*c1*gA2)/(fpi4*mpi2)- &
         (98.7791618408203_pr*c3*gA2)/(fpi4*mpi2)- &
         (70.582617388916_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,1,19)=(-0.68049621963501_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.150680578790392_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0729124902044024_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,1,19)=(-18.4566708641052_pr*c1*gA2)/(fpi4*mpi2)+ &
         (3.32744608689717_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2.34081944728579_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(3,1,19)=(-19.8544406890869_pr*c1*gA2)/(fpi4*mpi2)+ &
         (4.63924562215805_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2.88686197042465_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,2,20)=(19.0538941497803_pr*c1*gA2)/(fpi4*mpi2)- &
         (3.69410571834019_pr*c3*gA2)/(fpi4*mpi2)- &
         (2.33317974635533_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,1,21)=(0.340248109817505_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0753402893951961_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0364562451022012_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(3,1,21)=(-2.54145170688629_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.526430741071701_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.299806566953659_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(3,1,23)=(-0.0425310137271881_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.00941753617439951_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00455703063777515_pr*c4*gA2)/(fpi4*mpi2)
    !
    !  UJ1J1
    !
    ctJ1J1=0.0_pr
    ctJ1J1(1,1,3)=(0.0496322121163276_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0229839763853458_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.01171875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ1J1(1,1,5)=(0.196140339479882_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0113169982018643_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0413883948501276_pr*c4*gA2)/(fpi4*mpi2)+ &
         (0.005859375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ1J1(1,2,6)=(-0.0889442174269183_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0228348426870748_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,7)=(16.7237752493969_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.34790459899251_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00199338985105761_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,7)=(-3.22524030314923_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.137722860274259_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.000704460299744898_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.00146484375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ1J1(1,2,8)=(-19.6214015749008_pr*c1*gA2)/(fpi4*mpi2)- &
         (7.57700061161061_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0217966741071429_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,9)=(2190.64598296261_pr*c1*gA2)/(fpi4*mpi2)- &
         (106.868359765483_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.722775251415019_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,9)=(-703.063761547905_pr*c1*gA2)/(fpi4*mpi2)+ &
         (27.158281478954_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.115725330171131_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,9)=(58.6594921997924_pr*c1*gA2)/(fpi4*mpi2)
    ctJ1J1(1,3,9)=(46.67658203125_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.210943080357143_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,2,10)=(-7494.70074858631_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1278.61760582159_pr*c3*gA2)/(fpi4*mpi2)+ &
         (13.5948743024554_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,2,10)=(1169.53616071429_pr*c1*gA2)/(fpi4*mpi2)- &
         (294.126983802623_pr*c3*gA2)/(fpi4*mpi2)- &
         (2.59871651785714_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,11)=(-17607.1420353921_pr*c1*gA2)/(fpi4*mpi2)+ &
         (7994.24713374194_pr*c3*gA2)/(fpi4*mpi2)+ &
         (129.87408409526_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,11)=(9843.68403362975_pr*c1*gA2)/(fpi4*mpi2)- &
         (3932.01024951878_pr*c3*gA2)/(fpi4*mpi2)- &
         (56.7971205714634_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,11)=(-1084.01759961981_pr*c1*gA2)/(fpi4*mpi2)+ &
         (478.820795649102_pr*c3*gA2)/(fpi4*mpi2)+ &
         (6.05618722098214_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,3,11)=(6819.79372395833_pr*c1*gA2)/(fpi4*mpi2)- &
         (1268.16076692708_pr*c3*gA2)/(fpi4*mpi2)- &
         (16.6365703125_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,2,12)=(48059.0188986235_pr*c1*gA2)/(fpi4*mpi2)- &
         (18566.944349986_pr*c3*gA2)/(fpi4*mpi2)- &
         (352.160276785714_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,2,12)=(-16034.3911816406_pr*c1*gA2)/(fpi4*mpi2)+ &
         (4929.07689946187_pr*c3*gA2)/(fpi4*mpi2)+ &
         (83.798837890625_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,13)=(47847.3513627872_pr*c1*gA2)/(fpi4*mpi2)- &
         (13686.6719432792_pr*c3*gA2)/(fpi4*mpi2)- &
         (421.301021327427_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,13)=(-40664.4224147226_pr*c1*gA2)/(fpi4*mpi2)+ &
         (12306.5684431232_pr*c3*gA2)/(fpi4*mpi2)+ &
         (312.747646446845_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,13)=(8015.01712823947_pr*c1*gA2)/(fpi4*mpi2)- &
         (2193.04572012566_pr*c3*gA2)/(fpi4*mpi2)- &
         (48.5605316162109_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,3,13)=(-25770.923296875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (10166.4482434152_pr*c3*gA2)/(fpi4*mpi2)+ &
         (223.728783147321_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,2,14)=(-69245.3594485406_pr*c1*gA2)/(fpi4*mpi2)+ &
         (22270.458788223_pr*c3*gA2)/(fpi4*mpi2)+ &
         (807.664868858119_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,2,14)=(32217.3951881045_pr*c1*gA2)/(fpi4*mpi2)- &
         (11458.4663920128_pr*c3*gA2)/(fpi4*mpi2)- &
         (339.44725289481_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,15)=(-3355.18276569786_pr*c1*gA2)/(fpi4*mpi2)+ &
         (4003.15120252837_pr*c3*gA2)/(fpi4*mpi2)+ &
         (322.77466147213_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,15)=(12334.8284220602_pr*c1*gA2)/(fpi4*mpi2)- &
         (6594.07536966716_pr*c3*gA2)/(fpi4*mpi2)- &
         (376.149452475684_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,15)=(-4584.16347069416_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2259.89314979311_pr*c3*gA2)/(fpi4*mpi2)+ &
         (97.3375411987305_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,3,15)=(25694.7668613281_pr*c1*gA2)/(fpi4*mpi2)- &
         (8712.71905517578_pr*c3*gA2)/(fpi4*mpi2)- &
         (350.853002929688_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,2,16)=(8939.92116689046_pr*c1*gA2)/(fpi4*mpi2)- &
         (6189.35971202371_pr*c3*gA2)/(fpi4*mpi2)- &
         (428.127774734933_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,2,16)=(-11150.7683255005_pr*c1*gA2)/(fpi4*mpi2)+ &
         (5032.6588526575_pr*c3*gA2)/(fpi4*mpi2)+ &
         (271.22941115897_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,17)=(4508.67638249921_pr*c1*gA2)/(fpi4*mpi2)- &
         (1354.76069989184_pr*c3*gA2)/(fpi4*mpi2)- &
         (7.08009143020085_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,17)=(-6244.51076026543_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2125.28532690204_pr*c3*gA2)/(fpi4*mpi2)+ &
         (58.225499060767_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,17)=(2384.7618078426_pr*c1*gA2)/(fpi4*mpi2)- &
         (812.91389298321_pr*c3*gA2)/(fpi4*mpi2)- &
         (29.4953348595755_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,3,17)=(-2613.10548339844_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1796.50318320138_pr*c3*gA2)/(fpi4*mpi2)+ &
         (139.990072631836_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,2,18)=(-118.353195383708_pr*c1*gA2)/(fpi4*mpi2)- &
         (28.3697180204175_pr*c3*gA2)/(fpi4*mpi2)+ &
         (4.66635949271066_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,2,18)=(437.301789367676_pr*c1*gA2)/(fpi4*mpi2)- &
         (234.978396185244_pr*c3*gA2)/(fpi4*mpi2)- &
         (35.291308694458_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,19)=(2955.73723389918_pr*c1*gA2)/(fpi4*mpi2)- &
         (978.581137073306_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0364562451022012_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,19)=(-4606.79145692376_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1486.71857443532_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1.17040972364289_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,19)=(1575.92551413057_pr*c1*gA2)/(fpi4*mpi2)- &
         (499.689660100295_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1.44343098521233_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,2,20)=(29.5882988459269_pr*c1*gA2)/(fpi4*mpi2)+ &
         (7.09242950510437_pr*c3*gA2)/(fpi4*mpi2)- &
         (1.16658987317766_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,21)=(1316.75994490306_pr*c1*gA2)/(fpi4*mpi2)- &
         (449.140719963059_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,21)=(-2611.74143695995_pr*c1*gA2)/(fpi4*mpi2)+ &
         (874.719811383881_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0182281225511006_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,21)=(1112.6598731875_pr*c1*gA2)/(fpi4*mpi2)- &
         (369.221981860559_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.14990328347683_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,23)=(380.945590073708_pr*c1*gA2)/(fpi4*mpi2)- &
         (133.736494553329_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,23)=(-1001.70799803_pr*c1*gA2)/(fpi4*mpi2)+ &
         (344.82175777492_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,23)=(545.279167261999_pr*c1*gA2)/(fpi4*mpi2)- &
         (184.448333224776_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00227851531888757_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,25)=(69.0578717312122_pr*c1*gA2)/(fpi4*mpi2)- &
         (24.8647925755186_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,25)=(-255.206744216776_pr*c1*gA2)/(fpi4*mpi2)+ &
         (90.1447936087768_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,25)=(183.267451552542_pr*c1*gA2)/(fpi4*mpi2)- &
         (63.5695102280704_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,27)=(7.12383882476366_pr*c1*gA2)/(fpi4*mpi2)- &
         (2.62160304250119_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,27)=(-41.4401230244268_pr*c1*gA2)/(fpi4*mpi2)+ &
         (14.97425477113_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,27)=(41.6672205478068_pr*c1*gA2)/(fpi4*mpi2)- &
         (14.7931959736712_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,29)=(0.318977498914416_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.119616838695711_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,29)=(-3.88089691129625_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1.43041835994631_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,29)=(6.12365814763455_pr*c1*gA2)/(fpi4*mpi2)- &
         (2.21941836648652_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,31)=(-0.159488749457208_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0598084193478553_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,31)=(0.524984301276333_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.193754399830252_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,33)=(0.019936093682151_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.00747605241848192_pr*c3*gA2)/(fpi4*mpi2)
  End Subroutine load_tables
  !==================================================================================
#endif
  !
  !==================================================================================
  Subroutine read_UNEDF_NAMELIST(fname)
    !--------------------------------------------------------------------------------
    ! RESERVED NAMES ARE:
    !  -namelist forbiden:
    !          'UNRDF'  - best UNEDF
    !          'SKYRME' - best SKYRME
    !  -namelist inforced but not for C-parameters (use_INM=F) 
    !   or NM-parameters (use_INM=T) defined by the solver 
    !          'FITS'   
    !  -namelist inforced (one can overwrite all):
    !          'ANY OTHER NAME' 
    ! i.e., the DME solver defines C-/NM- only using 'FITS'
    !--------------------------------------------------------------------------------
    Implicit None
    Character (30), Intent(inout) :: fname 
    Character (30) :: inforcedname     
    Logical        :: regularization
    Integer(ipr)   :: ios,LOUT=6,lnamelist=16 
    !
    ! parameters
    use_DME3N_terms = .True.
    use_TMR_pairing=0
    eps=Spacing(1.0_pr)
    Pi=4.0_pr*Atan(1.0_pr)  
    kfconst=(1.50_pr*Pi**2)**(1.0_pr/3.0_pr)    ! (3Pi^2/2)^(1/3)
    CK=3.0_pr/5.0_pr*kfconst**2
    !
    Open(lnamelist,file='UNEDF_NAMELIST.DAT',DELIM='APOSTROPHE') ! 'QUOTE' 
    Do  
       !==================================  
       ! SOME DEFAULT VALUSE FOR ALL CASES  
       !==================================  
       Print_Namelist=.True.
       FunctionalName="Bla-Bla"
       use_charge_density=.False.
       regularization=.False.
       use_cm_cor=.True. 
       use_j2terms = .True.
       DMEorder=2 
       DMElda=0
       Crho(0)=-727.0933239596374733_pr; Crho(1)=474.8709969984467989_pr 
       CDrho(0)=612.1037411660222460_pr; CDrho(1)=-705.7204872069220301_pr
       Ctau(0)=33.8846741217252401_pr;   Ctau(1)=32.4047409594248919_pr
       CrDr(0)=-76.9962031249999939_pr;  CrDr(1)=15.6571351249999999_pr
       CrdJ(0)=-92.2500000000000000_pr;  CrdJ(1)=-30.7500000000000000_pr
       CJ(0)=17.2096115000000012_pr;     CJ(1)=64.5758124999999978_pr
       Cnrho=0.0_pr;                     CJdr=0.0_pr
       CpV0=-258.20_pr;                  CpV1=0.50_pr
       sigma=0.3062227576210547_pr;     
       hbzero=20.7355300000000007_pr;    
       e2charg = 1.43997840_pr ; CExPar = 1.0_pr
       mpi=138.03_pr/197.3_pr; fpi=92.4_pr/197.3_pr; gA=1.29_pr
       c1=-0.81_pr/1000.0_pr*197.3_pr
       c3=-3.4_pr/1000.0_pr*197.3_pr
       c4= 3.4_pr/1000.0_pr*197.3_pr
       cd=-2062.0_pr/1000.0_pr
       ce=-625.0_pr/1000.0_pr
       LambdaX=700.0_pr/197.3_pr
       E_NM=-15.972149141444596410_pr;  RHO_NM=0.159538756711733398_pr 
       K_NM=229.900964482603626493_pr;  SMASS_NM=1.439546988976078357_pr 
       ASS_NM=32.004302815052007247_pr; LASS_NM=45.961751480461579433_pr 
       VMASS_NM=1.249838547196253424_pr
       nuLambda=700.0_pr ; nufpi = 93.0_pr
       !======================
       ! SELECT THE FUNCTIONAL
       !======================
       FUNCTIONAL: Select Case (Trim(fname))
       Case ('UNEDF')
          inforcedname='UNEDF'
          use_INM=.True. 
          use_Namelist=.False.
          ! kind of the functional
          use_charge_density=.False.
          regularization=.False.
          use_cm_cor=.True. 
          DMEorder=2
          DMElda=0
          use_TMR_pairing=0
          ! Best DME parameters up to date
          Crho(0)=-727.0933239596374733_pr;  Crho(1)=474.8709969984467989_pr
          CDrho(0)=612.1037411660222460_pr;  CDrho(1)=-705.7204872069220301_pr
          Ctau(0)=33.8846741217252401_pr;    Ctau(1)=32.4047409594248919_pr
          CrDr(0)=-76.9962031249999939_pr;   CrDr(1)=15.6571351249999999_pr
          CrdJ(0)=-92.2500000000000000_pr;   CrdJ(1)=-30.7500000000000000_pr
          CJ(0)=17.2096115000000012_pr;      CJ(1)=64.5758124999999978_pr
          Cnrho=0.0_pr;                      CJdr=0.0_pr      
          CpV0=-258.20_pr;                   CpV1=0.50_pr
          sigma=0.3062227576210547_pr;     
          hbzero=20.7355300000000007_pr;    
          e2charg = 1.43997840_pr
          ! Associated NN and NNN parameters
          mpi=138.03_pr/197.3_pr; fpi=92.4_pr/197.3_pr; gA=1.29_pr
          c1=-0.81_pr/1000.0_pr*197.3_pr
          c3=-3.4_pr/1000.0_pr*197.3_pr
          c4= 3.4_pr/1000.0_pr*197.3_pr
          cd=-2062.0_pr/1000.0_pr
          ce=-625.0_pr/1000.0_pr
          LambdaX=700.0_pr/197.3_pr
          ! Associated INM parameters
          E_NM=-15.972149141444596410_pr;  RHO_NM=0.159538756711733398_pr  
          K_NM=229.900964482603626493_pr;  SMASS_NM=1.439546988976078357_pr       
          ASS_NM=32.004302815052007247_pr; LASS_NM=45.961751480461579433_pr 
          VMASS_NM=1.249838547196253424_pr 
       Case ("SKYRME")
          inforcedname='SKYRME'
          use_Namelist=.False.
          ! kind of the functional
          use_charge_density=.False.
          regularization=.False.
          use_cm_cor=.True. 
          DMEorder=-1
          DMElda=0
          use_TMR_pairing=0		  
          ! Best SKYRME parameters up to date
          hbzero=20.735530_pr;               sigma=0.1666666666666667_pr    
          Crho(0)=-933.3423749999999472_pr;  Crho(1)=830.0524855000001025_pr
          CDrho(0)=861.0625000000000000_pr;  CDrho(1)=-1064.2732499999999618_pr
          Ctau(0)=57.1286874999999981_pr;    Ctau(1)=24.6567364999999974_pr
          CrDr(0)=-76.9962031249999939_pr;   CrDr(1)=15.6571351249999999_pr 
          CrdJ(0)=-92.2500000000000000_pr;   CrdJ(1)=-30.7500000000000000_pr
          CJ(0)=17.2096115000000012_pr;      CJ(1)=64.5758124999999978_pr
          Cnrho=0.0_pr;                      CJdr=0.0_pr       
          CpV0=-258.20_pr;                   CpV1=0.50_pr
          sigma=0.1666666666666667_pr;     
          hbzero=20.7355300000000007_pr;    
          e2charg = 1.43997840_pr
          ! Associated INM parameters
          E_NM=-15.972149141444596410_pr;   RHO_NM=0.159538756711733398_pr 
          K_NM=229.900964482603626493_pr;   SMASS_NM=1.439546988976078357_pr 
          ASS_NM=32.004302815052007247_pr;  LASS_NM=45.961751480461579433_pr 
          VMASS_NM=1.249838547196253424_pr 
       Case default
          inforcedname=fname
          use_Namelist=.True.
       End Select FUNCTIONAL
       !======================
       ! READ NAMELISTS     
       !======================
       If(.Not.use_Namelist) Exit
       Read(UNIT=lnamelist,NML=UNEDF_NAMELIST,iostat=ios)
       If(ios.Ne.0) Then
          ! WRONG ENTRY WITHIN UNEDF_NAMELIST.DAT FILE
          Write(LOUT,'(1X,/,A)') 'ATTENTION: WRONG INPUT!'      
          Write(LOUT,*) 'THE INPUT DATA WITH LABEL FUNCTIONALNAME=''',Trim(INFORCEDNAME),''''
          Write(LOUT,*) 'INSIDE THE UNEDF_NAMELIST.DAT FILE IS WRONG.'      
          Write(LOUT,*) 'PLESE CORECT AND TRY AGAIN!'      
          Stop 'PROGRAM STOP IN read_UNEDF_NAMELIST'  
       Endif
       If(Trim(FunctionalName).Eq.Trim(inforcedname)) Exit
    Enddo
    Close(lnamelist)
    !
    !================================
    ! SEE WHAT THE NAMELISTS MODIFIED    
    !================================
    INFORCED_FUNCTIONAL: Select Case (Trim(inforcedname))
    Case ("UNEDF")
       FunctionalName='UNEDF'
    Case ("SKYRME")
       FunctionalName='SKYRME'
    Case default
       ! MISSING ENTRY WITHIN UNEDF_NAMELIST.DAT FILE
       If(Trim(FunctionalName).Ne.Trim(inforcedname)) Then
          Write(LOUT,'(1X,/,A)') 'ATTENTION: MISSING INPUT!'      
          Write(LOUT,*) 'THE INPUT DATA WITH LABEL FUNCTIONALNAME=''',Trim(INFORCEDNAME),''''
          Write(LOUT,*) 'IS MISSING INSIDE THE UNEDF_NAMELIST.DAT FILE.'      
          Write(LOUT,*) 'PLESE CORECT AND TRY AGAIN!'      
          Stop 'PROGRAM STOP IN SET_FUNCTIONAL_PARAMETERS'  
       Endif
    End Select INFORCED_FUNCTIONAL
    !  
  End Subroutine read_UNEDF_NAMELIST
  !==================================================================================
  !          
  !==================================================================================
  Subroutine set_functional_parameters(fname,lpr)
    !--------------------------------------------------------------------------------
    ! set functional parameters
    !--------------------------------------------------------------------------------
    Implicit None
    Logical, Intent(in) :: lpr
    Character (30), Intent(inout) :: fname     
    Logical :: regularization
    Integer(ipr), Parameter :: lin=15,lout=6
    !
    ! parameters
    FunctionalName=fname    
    eps=Spacing(1.0_pr)
    Pi=4.0_pr*Atan(1.0_pr)  
    kfconst=(1.50_pr*Pi**2)**(1.0_pr/3.0_pr)    ! (3Pi^2/2)^(1/3)
    CK=3.0_pr/5.0_pr*kfconst**2
    nuLambda=700.0_pr ; nufpi = 93.0_pr    
    !
    Call Make_Parameter_Free_Useful_Combintions()
    !
    ! exact Hartree CHrho from INM
    CHrho=0.0_pr; !!!!If (dmeorder.eq.3) Call CHrho_from_NM()
    !
    If(use_INM) Then 
     Call calculate_C_form_NM(E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM)
    Else
     Crho(0)=Crho(0)+CHrho                       !*(0.0_pr) !mario
    End If
    Call calculate_NM_properties()   
    !
    Crho(0)=Crho(0)-CHrho                        !*(0.00_pr) !mario
    !
    Call calculate_natural_units()
    !
    ! Print output
    If(lpr) Then
     Call print_functional_parameters(lout)
    Endif 
    !
  End Subroutine set_functional_parameters       
  !==================================================================================
  !          
  !==================================================================================
  Subroutine print_functional_parameters(lout)         
    !--------------------------------------------------------------------------------
    Implicit None
    Integer(ipr), Intent(in) :: lout 
    !
    Write(lout,'(a)')   '  ---------------------------------------'
    Write(lout,'(a,a)') '           UNEDF Module Version:', &
                                         Trim(Version)    
    Write(lout,'(a)')   '         M.Kortelainen & M.Stoitsov ' 
    Write(lout,'(a)')   '  ---------------------------------------'
    
    
    Write(lout,'(a)')    
    Write(lout,'(100(2x,a,a,f15.8))') Trim(FunctionalName),' functional'
    Write(lout,'(100(2x,a,f15.8))') '----------------------------------'
    Write(lout,'("  Crho(0)= ",g26.18,"; Crho(1)= ",g26.18)') Crho
    Write(lout,'("  CDrho(0)=",g26.18,"; CDrho(1)=",g26.18)') CDrho
    Write(lout,'("  Ctau(0)= ",g26.18,"; Ctau(1)= ",g26.18)') Ctau
    Write(lout,'("  CrDr(0)= ",g26.18,"; CrDr(1)= ",g26.18)') Crdr
    Write(lout,'("  CrdJ(0)= ",g26.18,"; CrdJ(1)= ",g26.18)') CrdJ
    Write(lout,'("  CJ(0)=   ",g26.18,"; CJ(1)=   ",g26.18)') CJ
    Write(lout,'("  CpV0(0)= ",g26.18,"; CpV0(1)= ",g26.18)') CpV0
    Write(lout,'("  CpV1(0)= ",g26.18,"; CpV1(1)= ",g26.18)') CpV1
    Write(lout,'("  sigma=   ",g26.18,"; hbzero=  ",g26.18)') sigma,hbzero 
    Write(lout,'("  e^2 chrg=",g26.18,"; CExPar=  ",g26.18)') e2charg,CExPar
    Write(lout,'("  c.m. correction: ",L1,", chr. density in direct Coul: ",L1)') use_cm_cor,use_charge_density
    Write(lout,'("  use tensor terms: ",L1)') use_j2terms
    Write(lout,'("  use TMR pairing:  ",I1)') use_TMR_pairing
        
    Write(lout,'(100(2x,a,f15.8))') 
    Write(lout,'(100(2x,a,f15.8))') 'Coupling constants in natural units'
    Write(lout,'(100(2x,a,f15.8))') '-----------------------------------'
    Write(lout,'("  Crho_nu(0)= ",g26.18,"; Crho_nu(1)= ",g26.18)') nuCrho
    Write(lout,'("  CDrho_nu(0)=",g26.18,"; CDrho_nu(1)=",g26.18)') nuCDrho
    Write(lout,'("  Ctau_nu(0)= ",g26.18,"; Ctau_nu(1)= ",g26.18)') nuCtau
    Write(lout,'("  CrDr_nu(0)= ",g26.18,"; CrDr_nu(1)= ",g26.18)') nuCrdr
    Write(lout,'("  CrdJ_nu(0)= ",g26.18,"; CrdJ_nu(1)= ",g26.18)') nuCrdJ
    Write(lout,'("  CJ_nu(0)=   ",g26.18,"; CJ_nu(1)=   ",g26.18)') nuCJ
    Write(lout,'("  CpV0_nu(0)= ",g26.18,"; CpV0_nu(1)= ",g26.18)') nuCpV0
    Write(lout,'("  CpV1_nu(0)= ",g26.18,"; CpV1_nu(1)= ",g26.18)') nuCpV1
    Write(lout,'("  fpi_nu=     ",g26.18,"; Lambda_nu=  ",g26.18)') nufpi,nuLambda
    
    If (dmeorder.Ge.0) Then
      Write(lout,'(100(2x,a,f15.8))') 
      Write(lout,'(100(2x,a,f15.8))') 'DME parameters'
      Write(lout,'(100(2x,a,f15.8))') '----------------------------------'
      Write(lout,'("       gA=",f12.6," mpi [1/fm]=",f12.6," fpi [1/fm]=",f12.6)') gA,mpi,fpi
      Write(lout,'("  c1 [fm]=",f12.6,"    c3 [fm]=",f12.6,"    c4 [fm]=",f12.6)') c1,c3,c4
      Write(lout,'("       cd=",f12.6,"         ce=",f12.6," LamX[1/fm]=",f12.6)') cd,ce,LambdaX
      Write(lout,'("  ->CHrho=",f12.6)') CHrho
      If (dmeorder.Ge.2) Write(lout,'("  use 3N terms: ",L1)') use_DME3N_terms
    End If
    
    Write(lout,'(100(2x,a,f15.8))') 
    Write(lout,'(100(2x,a,f15.8))') 'Nuclear matter properties'
    Write(lout,'(100(2x,a,f15.8))') '----------------------------------'
    Write(lout,'(100(2x,a9,f25.16))') 'E_NM=',E_NM,'K_NM=',K_NM
    Write(lout,'(100(2x,a9,f25.16))') 'P_NM=',P_NM,'RHO_NM=',RHO_NM
    Write(lout,'(100(2x,a9,f25.16))') 'ASS_NM=',ASS_NM,'LASS_NM=',LASS_NM
    Write(lout,'(100(2x,a9,f25.16))') 'SMASS_NM=',SMASS_NM,'VMASS_NM=',VMASS_NM
    
    Call t_from_C()
    Write(lout,'(100(2x,a,f15.8))') 
    Write(lout,'(100(2x,a,f15.8))') 'Associated (t,x)-coupling constants' 
    Write(lout,'(100(2x,a,f15.8))') '-----------------------------------'
    Write(lout,'("  t0=    ",g26.18,"; x0=     ",g26.18)') t0,x0
    Write(lout,'("  t1=    ",g26.18,"; x1=     ",g26.18)') t1,x1
    Write(lout,'("  t2=    ",g26.18,"; x2=     ",g26.18)') t2,x2
    Write(lout,'("  t3=    ",g26.18,"; x3=     ",g26.18)') t3,x3
    Write(lout,'("  b4=    ",g26.18,"; b4p=    ",g26.18)') b4,b4p
    Write(lout,'("  te=    ",g26.18,"; to=     ",g26.18)') te,to
    Write(lout,'("  sigma= ",g26.18,"; hbzero= ",g26.18)') sigma,hbzero
    
    If(Print_Namelist) Then
     Write(lout,'(100(2x,a,f15.8))') 
     SELECTED_FUNCTIONAL: Select Case (Trim(FunctionalName))
      Case ("UNEDF","SKYRME")
       Write(lout,'(100(2x,a,f15.8))') 'NAMELIST CONTENT (cannot be modified for functional names UNEDF,SKYRME)'
       Write(lout,'(100(2x,a,f15.8))') '-----------------------------------------------------------------------'
      Case ("FITS")
       Write(lout,'(100(2x,a,f15.8))') 'NAMELIST CONTENT (Advane usage: modify all but not C-, NM-, and more...)'
       Write(lout,'(100(2x,a,f15.8))') '-----------------------------------------------------------------------'
      Case default
       Write(lout,'(100(2x,a,f15.8))') 'NAMELIST CONTENT (copy/past to UNEDF_NAMELIST.DAT and modify)'
       Write(lout,'(100(2x,a,f15.8))') '-------------------------------------------------------------'
     End Select SELECTED_FUNCTIONAL
     Write(lout,'(100(a,f15.8))')    ' !NB: FUNCTIONALNAME should be always in quotations'
     Write (lout,UNEDF_NAMELIST)     
    End If
  End Subroutine print_functional_parameters
  !==================================================================================
  !
  !==================================================================================
  Subroutine calculate_natural_units
    ! Calculates coupling constants in natural units
    !-----------------------------------------------------------------------------------
    Implicit None
    nuCrho = Crho*(nufpi**2)/(mevfm**3)
    nuCdrho = Cdrho*(nufpi**2)*((nuLambda*nufpi*nufpi)**sigma)/(mevfm**(3.0_pr*(1.0_pr+sigma)))
    nuCtau = Ctau*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCrDr = CrDr*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCrdJ = CrdJ*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCJ = CJ*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCpV0 = CpV0*(nufpi**2)/(mevfm**3)
    nuCpV1 = CpV1*(nufpi**4)*nuLambda/(mevfm**6)
  End Subroutine calculate_natural_units
  !==================================================================================
  ! 
  !==================================================================================
  Subroutine calculate_C_form_NM(E,K,SMASS,RHO,ASS,LASS,VMASS,sigma_NM)
    ! Calculates volume C-constants (and sigma) form NM properties
    ! Interface usage: 
    !  hbzero,CK,kfconst,mpi,sigma
    !  aRhoRho,bRhoRho...
    !  hRho0Rho0,dhRho0Rho0...
    !  Crho(0),Crho(1),Cdrho(0),Cdrho(1),Ctau(0),Ctau(0)
    !  subroutine calculate_U_parameters
    !
    !  input: E,K,SMASS,RHO,ASS,LASS,VMASS,sigma_NM(optional)
    ! output: Crho(0),Crho(1),Cdrho(0),Cdrho(1),Ctau(0),Ctau(0),sigma(optional)
    ! 
    ! Options:
    !  When sigma_NM exists then 'sigma'=sigma_NM
    !  When sigma_NM does not exist then 'sigma' is defined from NM
    !------------------------------------------------------------------
    Implicit None
    Real(pr), Intent(in) :: E,K,SMASS,RHO,ASS,LASS,VMASS
    Real(pr), Intent(in), Optional :: sigma_NM
    Real(pr) :: aRho0Rho0,daRho0Rho0,ddaRho0Rho0,aRho1Rho1,daRho1Rho1,ddaRho1Rho1
    Real(pr) :: aRho0Tau0,daRho0Tau0,ddaRho0Tau0,aRho1Tau1,daRho1Tau1,ddaRho1Tau1 
    Real(pr) :: u,tauc,rho2
    Real(pr),Parameter :: c13=1.0_pr/3.0_pr,c23=2.0_pr/3.0_pr
    !
    tauc=CK*RHO**c23; u=(kfconst/mpi)*RHO**c13; rho2=rho**2
    Call calculate_U_parameters(RHO,RHO,tauc*RHO,tauc*RHO,0.0_pr,0.0_pr)
    aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm
    aRho1Rho1=0.50_pr*(aRhoRho-bRhoRho)*mevfm
    aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm
    aRho1Tau1=0.50_pr*(aRhoTau-bRhoTau)*mevfm
    daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
    daRho1Rho1=0.50_pr*(daRhoRho-dbRhoRho)*mevfm
    daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
    daRho1Tau1=0.50_pr*(daRhoTau-dbRhoTau)*mevfm
    ddaRho0Rho0=0.50_pr*(ddaRhoRho+ddbRhoRho)*mevfm
    ddaRho1Rho1=0.50_pr*(ddaRhoRho-ddbRhoRho)*mevfm
    ddaRho0Tau0=0.50_pr*(ddaRhoTau+ddbRhoTau)*mevfm
    ddaRho1Tau1=0.50_pr*(ddaRhoTau-ddbRhoTau)*mevfm
    !
    ! set/calculate sigma
    If (Present(sigma_NM)) Then 
     sigma=sigma_NM 
    Else 
     sigma=((1.0_pr/3.0_pr)*(-K+tauc*hbzero*(-3.0_pr+4.0_pr*SMASS)-9.0_pr*E+9.0_pr*RHO2*hRho0Rho0 & 
             +21.0_pr*tauc*RHO2*hRho0Tau0+u*RHO*(daRho0Rho0+5.0_pr*tauc*daRho0Tau0 & 
             +7.0_pr*RHO*dhRho0Rho0+11.0_pr*tauc*RHO*dhRho0Tau0+u*ddaRho0Rho0 & 
             +u*tauc*ddaRho0Tau0+u*RHO*ddhRho0Rho0+u*tauc*RHO*ddhRho0Tau0))) & 
             /(tauc*hbzero*(-3.0_pr+2.0_pr*SMASS)+3.0_pr*E+3.0_pr*RHO2*hRho0Rho0 & 
             +3.0_pr*tauc*RHO2*hRho0Tau0+u*RHO*(daRho0Rho0+tauc*daRho0Tau0 & 
             + RHO*dhRho0Rho0+tauc*RHO*dhRho0Tau0)) 
    End If
    !  
    Crho(0)=(c13*(tauc*hbzero*(-3.0_pr+(2.0_pr-3.0_pr*sigma)*SMASS) &
        +3.0_pr*(1.0_pr+sigma)*E-3.0_pr*sigma*RHO*aRho0Rho0 &
        +3.0_pr*(1.0_pr-sigma)*RHO2*hRho0Rho0+3.0_pr*tauc*RHO2*hRho0Tau0 &
        +u*RHO*(daRho0Rho0+tauc*daRho0Tau0+RHO*dhRho0Rho0 &
        +tauc*RHO*dhRho0Tau0)))/(sigma*RHO)
    Cdrho(0)=(c13*RHO**(-1.0_pr-sigma)*(tauc*hbzero*(3.0_pr-2.0_pr*SMASS)&
        -3.0_pr*E-3.0_pr*RHO**2*hRho0Rho0-3.0_pr*tauc*RHO2*hRho0Tau0&
        -u*RHO*(daRho0Rho0+tauc*daRho0Tau0+RHO*dhRho0Rho0 &
        +tauc*RHO*dhRho0Tau0)))/sigma
    Ctau(0)=(hbzero*(SMASS-1.0_pr)-RHO*(aRho0Tau0+RHO*hRho0Tau0))/RHO
    !
    Crho(1)=(27.0_pr*ASS*(1.0_pr+sigma)-9.0_pr*LASS &
        +5.0_pr*tauc*hbzero*(5.0_pr-6.0_pr*VMASS+3.0_pr*sigma*(-4.0_pr+3.0_pr*VMASS)) &
        +20.0_pr*tauc*(2.0_pr-3.0_pr*sigma)*RHO*aRho0Tau0 &
        +RHO*(-27.0_pr*sigma*aRho1Rho1+5.0_pr*tauc*(11.0_pr-12.0_pr*sigma)*RHO*hRho0Tau0 &
        -27.0_pr*(-1.0_pr+sigma)*RHO*hRho1Rho1+9.0_pr*tauc*(5.0_pr-3.0_pr*sigma)*RHO*hRho1Tau0 &
        +45.0_pr*tauc*RHO*hRho1Tau1+40.0_pr*tauc*Ctau(0)-60.0_pr*tauc*sigma*Ctau(0) &
        +5.0_pr*u*tauc*daRho0Tau0+9.0_pr*u*daRho1Rho1+15.0_pr*u*tauc*daRho1Tau1 &
        +5.0_pr*u*tauc*RHO*dhRho0Tau0+9.0_pr*u*RHO*dhRho1Rho1+9.0_pr*u*tauc*RHO*dhRho1Tau0 &
        +15.0_pr*u*tauc*RHO*dhRho1Tau1))/(27.0_pr*sigma*RHO)
     Cdrho(1)=-(RHO**(-1.0_pr-sigma)*(27.0_pr*ASS-9.0_pr*LASS &
        +5.0_pr*tauc*hbzero*(5.0_pr-6.0_pr*VMASS)+40.0_pr*tauc*RHO*aRho0Tau0 &
        +55.0_pr*tauc*RHO2*hRho0Tau0+27.0_pr*RHO**2*hRho1Rho1+45.0_pr*tauc*RHO2*hRho1Tau0 &
        +45.0_pr*tauc*RHO2*hRho1Tau1+40.0_pr*tauc*RHO*Ctau(0) +5.0_pr*u*tauc*RHO*daRho0Tau0 &
        +9.0_pr*u*RHO*daRho1Rho1+15.0_pr*u*tauc*RHO*daRho1Tau1 &
        +5.0_pr*u*tauc*RHO2*dhRho0Tau0+9.0_pr*u*RHO2*dhRho1Rho1 &
        +9.0_pr*u*tauc*RHO2*dhRho1Tau0 +15.0_pr*u*tauc*RHO2*dhRho1Tau1))/(27.0_pr*sigma)
     Ctau(1)=(hbzero-hbzero*VMASS+RHO*(aRho0Tau0-aRho1Tau1+RHO*hRho0Tau0-RHO*hRho1Tau1+Ctau(0)))/RHO
     !
  End Subroutine calculate_C_form_NM
  !==================================================================================
  !
  !==================================================================================
  Subroutine calculate_NM_properties()
    ! Calculates INM properties
    ! Interface usage: 
    !  hbzero,CK,kfconst,mpi,sigma
    !  aRhoRho,bRhoRho...
    !  hRho0Rho0,dhRho0Rho0...
    !  Crho(0),Crho(1),Cdrho(0),Cdrho(1),Ctau(0),Ctau(0)
    !  E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM,sigma,P_NM,KA_NM
    !  function find_NM_RHOC()
    ! input:  Crho(0),Crho(1),Cdrho(0),Cdrho(1),Ctau(0),Ctau(0),sigma
    ! output: E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM,sigma,P_NM,KA_NM
    ! Using:
    !  RHO_NM=find_NM_RHOC()
    !------------------------------------------------------------------ 
    Implicit None
    Real(pr) :: aRho0Rho0,daRho0Rho0,ddaRho0Rho0,aRho1Rho1,daRho1Rho1,ddaRho1Rho1
    Real(pr) :: aRho0Tau0,daRho0Tau0,ddaRho0Tau0,aRho1Tau1,daRho1Tau1,ddaRho1Tau1 
    Real(pr) :: u,tauc,rho_NM2
    Real(pr), Parameter :: c13=1.0_pr/3.0_pr,c23=2.0_pr/3.0_pr
    !
    RHO_NM=find_NM_RHOC()
    !
    aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm
    aRho1Rho1=0.50_pr*(aRhoRho-bRhoRho)*mevfm
    aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm
    aRho1Tau1=0.50_pr*(aRhoTau-bRhoTau)*mevfm
    daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
    daRho1Rho1=0.50_pr*(daRhoRho-dbRhoRho)*mevfm
    daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
    daRho1Tau1=0.50_pr*(daRhoTau-dbRhoTau)*mevfm
    ddaRho0Rho0=0.50_pr*(ddaRhoRho+ddbRhoRho)*mevfm
    ddaRho1Rho1=0.50_pr*(ddaRhoRho-ddbRhoRho)*mevfm
    ddaRho0Tau0=0.50_pr*(ddaRhoTau+ddbRhoTau)*mevfm
    ddaRho1Tau1=0.50_pr*(ddaRhoTau-ddbRhoTau)*mevfm
    tauc=CK*RHO_NM**c23; u=(kfconst/mpi)*RHO_NM**c13; rho_NM2=rho_NM**2
    !
    ! Symmetric Nuclear Matter 
    E_NM=tauc*hbzero+RHO_NM*(aRho0Rho0+RHO_NM*hRho0Rho0+Crho(0)+RHO_NM**sigma*Cdrho(0)) &
      +tauc*RHO_NM*(aRho0Tau0+RHO_NM*hRho0Tau0+Ctau(0))
    P_NM=c13*RHO_NM**2*((2.0_pr*tauc*hbzero)/RHO_NM+3.0_pr*aRho0Rho0+5.0_pr*tauc*aRho0Tau0 &
      +6.0_pr*RHO_NM*hRho0Rho0+8.0_pr*tauc*RHO_NM*hRho0Tau0+3.0_pr*Crho(0) &
      +3.0_pr*(1+sigma)*RHO_NM**sigma*Cdrho(0)+5.0_pr*tauc*Ctau(0)+u*daRho0Rho0 &
      +u*tauc*daRho0Tau0+u*RHO_NM*dhRho0Rho0+u*tauc*RHO_NM*dhRho0Tau0)
    SMASS_NM=1.0_pr+(RHO_NM*(aRho0Tau0+RHO_NM*hRho0Tau0+Ctau(0)))/hbzero
    K_NM=9.0_pr*sigma*(1+sigma)*RHO_NM**(1+sigma)*Cdrho(0) &
      +(-2.0_pr*tauc*hbzero+10.0_pr*tauc*RHO_NM*aRho0Tau0+18.0_pr*RHO_NM2*hRho0Rho0 &
      +40.0_pr*tauc*RHO_NM**2*hRho0Tau0+4.0_pr*u*RHO_NM*daRho0Rho0 & 
      +RHO_NM*(10.0_pr*tauc*Ctau(0)+u*(8.0_pr*tauc*daRho0Tau0+u*ddaRho0Rho0 &
      +(10.0_pr*RHO_NM*dhRho0Rho0+14.0_pr*tauc*RHO_NM*dhRho0Tau0 &
      +(u*tauc*ddaRho0Tau0+u*RHO_NM*ddhRho0Rho0+u*tauc*RHO_NM*ddhRho0Tau0)))))
    ! 
    ! Asymmetric Nuclear Matter 
    ASS_NM=RHO_NM2*hRho1Rho1+RHO_NM*(aRho1Rho1+Crho(1)+RHO_NM**sigma*Cdrho(1)) &
       +(tauc*(5.0_pr*hbzero+RHO_NM*(5.0_pr*aRho0Tau0+15.0_pr*aRho1Tau1+5.0_pr*RHO_NM*hRho0Tau0 &
       +9.0_pr*RHO_NM*hRho1Tau0+5.0_pr*(3.0_pr*RHO_NM*hRho1Tau1+Ctau(0)+3.0_pr*Ctau(1)))))/9.0_pr
    VMASS_NM=(hbzero+RHO_NM*(aRho0Tau0-aRho1Tau1+RHO_NM*hRho0Tau0-RHO_NM*hRho1Tau1+Ctau(0)-Ctau(1)))/hbzero
    LASS_NM=6.0_pr*RHO_NM2*hRho1Rho1+3.0_pr*RHO_NM*(aRho1Rho1+Crho(1)+(1.0_pr+sigma)*RHO_NM**sigma*Cdrho(1)) &
       +u*RHO_NM*daRho1Rho1 +u*RHO_NM2*dhRho1Rho1 &
       +(tauc*(10.0_pr*hbzero+8.0_pr*RHO_NM2*(5.0_pr*hRho0Tau0+9.0_pr*hRho1Tau0+15.0_pr*hRho1Tau1) &
       +25.0_pr*RHO_NM*(aRho0Tau0+3.0_pr*aRho1Tau1+Ctau(0)+3*Ctau(1)) &
       +5.0_pr*u*RHO_NM*(daRho0Tau0+3.0_pr*daRho1Tau1) &
       +u*RHO_NM2*(5.0_pr*dhRho0Tau0+9.0_pr*dhRho1Tau0+15.0_pr*dhRho1Tau1)))/9.0_pr
    KA_NM=18.0_pr*RHO_NM2*hRho1Rho1+9.0_pr*sigma*(1.0_pr+sigma)*RHO_NM**(1.0_pr+sigma)*Cdrho(1) &
       +4.0_pr*u*RHO_NM*daRho1Rho1 +10.0_pr*u*RHO_NM2*dhRho1Rho1 &
       + u**2*RHO_NM*ddaRho1Rho1+u**2*RHO_NM2*ddhRho1Rho1 &
       +(tauc*(-10.0_pr*hbzero+40.0_pr*RHO_NM2*(5.0_pr*hRho0Tau0+9.0_pr*hRho1Tau0+15.0_pr*hRho1Tau1) &
       +50.0_pr*RHO_NM*(aRho0Tau0+3.0_pr*aRho1Tau1+Ctau(0)+3*Ctau(1)) &
       +40.0_pr*u*RHO_NM*(daRho0Tau0+3.0_pr*daRho1Tau1) &
       +14.0_pr*u*RHO_NM2*(5.0_pr*dhRho0Tau0+9.0_pr*dhRho1Tau0 &
       +15.0_pr*dhRho1Tau1)+5.0_pr*u**2*RHO_NM*(ddaRho0Tau0 &
       +3.0_pr*ddaRho1Tau1)+u**2*RHO_NM2*(5.0_pr*ddhRho0Tau0+9*ddhRho1Tau0+15*ddhRho1Tau1)))/9.
     !       
  End Subroutine calculate_NM_properties
  !==================================================================================  
  !
  !=======================================================================
  Function find_NM_RHOC()
    ! Search for the INM saturation density RHO_NM using the Secant Method
    !====================================================================
    Implicit None
    !Integer(pr) intent(out) :: ierr
    Integer(pr) :: iter
    Real(pr) :: aRho0Rho0,daRho0Rho0,ddaRho0Rho0,aRho1Rho1,daRho1Rho1,ddaRho1Rho1
    Real(pr) :: aRho0Tau0,daRho0Tau0,ddaRho0Tau0,aRho1Tau1,daRho1Tau1,ddaRho1Tau1 
    Real(pr) :: find_NM_RHOC,kfconstmpi,u,tauc
    Real(pr) :: rhom0,rhom,rhom2,w,w0,step,energy
    Real(pr),Parameter :: c13=1.0_pr/3.0_pr,c23=2.0_pr/3.0_pr
    !
    kfconstmpi=kfconst/mpi; step=-0.0010_pr; iter=0
    ! initial point
    rhom=0.170_pr; tauc=CK*rhom**c23; u=kfconstmpi*rhom**c13; rhom2=rhom**2
    Call calculate_U_parameters(rhom,rhom,tauc*rhom,tauc*rhom,0.0_pr,0.0_pr)
    aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm; daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
    aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm; daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
    w0=c13*rhom2*((2.0_pr*tauc*hbzero)/rhom+3.0_pr*aRho0Rho0+5.0_pr*tauc*aRho0Tau0 &
      +6.0_pr*rhom*hRho0Rho0+8.0_pr*tauc*rhom*hRho0Tau0+3.0_pr*Crho(0) &
      +3.0_pr*(1.0_pr+sigma)*rhom**sigma*Cdrho(0)+5.0_pr*tauc*Ctau(0)+u*daRho0Rho0 &
      +u*tauc*daRho0Tau0+u*rhom*dhRho0Rho0+u*tauc*rhom*dhRho0Tau0)
    rhom0=rhom; rhom=rhom+step
    !	
    ! secant method
    Do While(Abs(step).Ge.eps*100.0_pr)
       iter=iter+1
       tauc=CK*rhom**c23; u=kfconstmpi*rhom**c13; rhom2=rhom**2
       Call calculate_U_parameters(rhom,rhom,tauc*rhom,tauc*rhom,0.0_pr,0.0_pr)
       aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm; daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
       aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm; daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
       w=c13*rhom2*((2.0_pr*tauc*hbzero)/rhom+3.0_pr*aRho0Rho0+5.0_pr*tauc*aRho0Tau0 &
         +6.0_pr*rhom*hRho0Rho0+8.0_pr*tauc*rhom*hRho0Tau0+3.0_pr*Crho(0) &
         +3.0_pr*(1.0_pr+sigma)*rhom**sigma*Cdrho(0)+5.0_pr*tauc*Ctau(0)+u*daRho0Rho0 &
         +u*tauc*daRho0Tau0+u*rhom*dhRho0Rho0+u*tauc*rhom*dhRho0Tau0)
       step=-w*(rhom-rhom0)/(w-w0)
       rhom0=rhom; w0=w; rhom=rhom+step       
       If(iter.Gt.100) Stop 'STOP(In find_NM_RHOC)'
       !energy=tauc*hbzero+rhom*(aRho0Rho0+rhom*hRho0Rho0+Crho(0)+rhom**sigma*Cdrho(0)) &
       ! +tauc*rhom*(aRho0Tau0+rhom*hRho0Tau0+Ctau(0))        
       !write(6,'(a,15(1pg12.4))') ' rhom0,rhom,step,e,w=',rhom0,rhom,step,energy,w
    Enddo
    find_NM_RHOC=rhom
  End Function find_NM_RHOC
  !==================================================================================
  !
  !==================================================================================
  Subroutine C_from_t()
    !--------------------------------------------------------------------------------
    ! C- from (t,x)-
    !--------------------------------------------------------------------------------
    Implicit None
     Crho(0) =  3.0_pr/8.0_pr*t0                            
     Cdrho(0) =  (1.0_pr/16.0_pr)*t3
     Crho(1) = -(1.0_pr/4.0_pr)*t0*(0.50_pr+x0)   
     Cdrho(1) = -(1.0_pr/24.0_pr)*t3*(0.50_pr+x3)   
     Ctau(0)  =  (3.0_pr/16.0_pr)*t1+(1.0_pr/4.0_pr)*t2*(5.0_pr/4.0_pr+x2)   
     Ctau(1)  = -(1.0_pr/8.0_pr)*t1*(0.5+x1)+(1.0_pr/8.0_pr)*t2*(0.50_pr+x2)   
     CrDr(0) =  (1.0_pr/16.0_pr)*t2*(5.0_pr/4.0_pr+x2)-(9.0_pr/64.0_pr)*t1
     CrDr(1) =  (3.0_pr/32.0_pr)*t1*(0.5+x1)+(1.0_pr/32.0_pr)*t2*(0.50_pr+x2)   
     CJ(0)    = -(1.0_pr/16.0_pr)*(t1*(2.0_pr*x1-1.0_pr)+t2*(2.0_pr*x2+1)-5*te-15*to)
     CJ(1)    = -(1.0_pr/16.0_pr)*(t2 -t1 + 5.0_pr*te -5.0_pr*to )
     CrdJ(0)   = -b4-(0.50_pr)*b4p              
     CrdJ(1)   = -0.50_pr*b4p                       
  End Subroutine C_from_t
  !==================================================================================
  !
  !==================================================================================
  Subroutine t_from_C()
    !--------------------------------------------------------------------------------
    ! (t,x)- from C-
    !--------------------------------------------------------------------------------
    Implicit None
    t0     =  (8.0_pr/3)*Crho(0)
    t1     =  4.0_pr/3.0_pr*(Ctau(0)-4.0_pr*CrDr(0))         
    t2     =  4.0_pr/3.0_pr*(3.0_pr*Ctau(0)-6.0_pr*Ctau(1)+4.0_pr*CrDr(0)-8.0_pr*CrDr(1))  
    t3     =  16.0_pr*Cdrho(0)
    x0     = -0.50_pr*(3.0_pr*Crho(1)/Crho(0)+1.0_pr)   
    x1     =  2.0_pr*(-Ctau(0)-3.0_pr*Ctau(1)+4.0_pr*CrDr(0)+12.0_pr*CrDr(1))/t1/3.0_pr
    x2     = -2.0_pr*(3.0_pr*Ctau(0)-15.0_pr*Ctau(1)+4.0_pr*CrDr(0)-20.0_pr*CrDr(1))/t2/3.0_pr
    x3     = -0.50_pr*(3.0_pr*Cdrho(1)/Cdrho(0)+1.0_pr) 
    b4     =  CrdJ(1)-CrdJ(0)
    b4p    = -2.0_pr*CrdJ(1)
    te     = (4.0_pr/15.0_pr)*(3.0_pr*CJ(0)-9.0_pr*CJ(1)-4.0_pr*CrDr(0)+12.0_pr*CrDr(1)-2.0_pr*Ctau(0)+6.0_pr*Ctau(1))
    to     = (4.0_pr/15.0_pr)*(3.0_pr*CJ(0)+3.0_pr*CJ(1)+4.0_pr*CrDr(0)+4.0_pr*CrDr(1))
  End Subroutine t_from_C
  !==================================================================================
  !   
  !==================================================================================
  Subroutine CHrho_from_NM()
    !--------------------------------------------------------------------------------
    ! CHrho from NM, E_NM(Hartree)=CHrho*RHO_NM
    !--------------------------------------------------------------------------------
    Implicit None
	Real(pr) :: z3=1.50_pr
        !
        !!CHrho= & 
        !!+h0mpi6c3NM*(A3_1/b3_1**z3+A3_2/b3_2**z3+A3_3/b3_3**z3+A3_4/b3_4**z3+A3_5/b3_5**z3) &
        !!+h0mpi6c1NM*(A1_1/b1_1**z3+A1_2/b1_2**z3+A1_3/b1_3**z3+A1_4/b1_4**z3+A1_5/b1_5**z3)
        CHrho = 0.0_pr
        !
  End Subroutine CHrho_from_NM
  !==================================================================================
  !
  !==================================================================================
  Elemental Function HartreeV00(u)
    !--------------------------------------------------------------------------------
    ! HartreeV(u), E(Hartree)=(1/2)*Int[rho_0(r)*V(|r-r'|)*rho_0(r')]
    !--------------------------------------------------------------------------------
    Implicit None
	Real(pr), Intent(in) :: u 
	Real(pr)             :: x2,HartreeV00 
        !
        !!x2=(u*mpi)**2
	!
        !!HartreeV=h0mpi6c1*(Exp(-x2*b1_1)*A1_1+Exp(-x2*b1_2)*A1_2+Exp(-x2*b1_3)*A1_3+Exp(-x2*b1_4)*A1_4+Exp(-x2*b1_5)*A1_5)+&
        !!h0mpi6c3*(Exp(-x2*b3_1)*A3_1+Exp(-x2*b3_2)*A3_2+Exp(-x2*b3_3)*A3_3+Exp(-x2*b3_4)*A3_4+Exp(-x2*b3_5)*A3_5)
	!
        HartreeV00=0.0_pr
        !
  End Function HartreeV00
  !
  Elemental Function HartreeV01(u)
    !--------------------------------------------------------------------------------
    ! HartreeV(u), E(Hartree)=(1/2)*Int[rho_0(r)*V(|r-r'|)*rho_1(r')]
    !--------------------------------------------------------------------------------
    Implicit None
	Real(pr), Intent(in) :: u 
	Real(pr)             :: x2,HartreeV01 
        !
        HartreeV01=0.0_pr
        !
  End Function HartreeV01
  !
  Elemental Function HartreeV11(u)
    !--------------------------------------------------------------------------------
    ! HartreeV(u), E(Hartree)=(1/2)*Int[rho_1(r)*V(|r-r'|)*rho_1(r')]
    !--------------------------------------------------------------------------------
    Implicit None
	Real(pr), Intent(in) :: u 
	Real(pr)             :: x2,HartreeV11 
        !
        HartreeV11=0.0_pr
        !
  End Function HartreeV11
  !==================================================================================
  !
  !==================================================================================
  Elemental Function ThetaFunction2(u)
    !--------------------------------------------------------------------------------
    ! ThetaFunction2(u)=0 or 1  when x2<2  or x2>2
    !--------------------------------------------------------------------------------
    Implicit None
	Real(pr), Intent(in) :: u 
	Real(pr)             :: x2,ThetaFunction2 
    !
	x2=(u*mpi)
	!
	ThetaFunction2=0.0_pr
	If(x2.Gt.2.0_pr) ThetaFunction2=1.0_pr
	!
  End Function ThetaFunction2
  !================================================================================== 
  !
  !==================================================================================
  Subroutine Make_Parameter_Free_Useful_Combintions()
    !--------------------------------------------------------------------------------
    ! Make Useful combintions
    !--------------------------------------------------------------------------------
    Implicit None
    !
    If (dmeorder.Ge.0) Then
       mpi2=mpi**2
       gA2=gA**2; gA4=gA2**2; gA6=gA2**3; 
       fpi2=fpi**2; fpi4=fpi2**2;        
       CHartree =mevfm*(3.0_pr*gA2)/(32.0_pr*fpi4*Pi**2) 
	   h0mpi6=197.30_pr*(mpi**6)*(3.0_pr*gA*gA)/(32.0_pr*fpi**4*Pi**2)
	   h0mpi6c1=h0mpi6*c1;         h0mpi6c3=h0mpi6*c3	 
       !
	   h0mpi6NM=197.30_pr*(3.0_pr*(mpi**3)*gA2)/(64.0_pr*fpi**4*Sqrt(Pi))
	   h0mpi6c1NM=h0mpi6NM*c1;     h0mpi6c3NM=h0mpi6NM*c3	 
	   !
	   A3_1=42.7132145164590_pr;   A3_2=0.670441422115440_pr; A3_3=0.0525713896514650_pr; 
	   A3_4=0.0012545731701320_pr; A3_5=5.81008627207380_pr
	   b3_1=3.0809379008590_pr;    b3_2=0.905186811964580_pr; b3_3=0.474514509597610_pr; 
	   b3_4=0.228138177966090_pr;  b3_5=1.66931540698090_pr;
	   !
	   A1_1=2.5000830618386_pr;    A1_2=0.619542286897850_pr; A1_3=0.169682589033730_pr;    
	   A1_4=0.0276112113725470_pr; A1_5=0.00108164458809540_pr
	   b1_1=1.75854210706510_pr;   b1_2=0.88882524524657_pr;  b1_3=0.46377235143756_pr;   
	   b1_4=0.247665887704790_pr;  b1_5=0.132222413002680_pr
	   !
#ifndef hide_dme    
	   Call load_tables()
#endif
	   !
    Endif
    !
  End Subroutine Make_Parameter_Free_Useful_Combintions
  !==================================================================================
  !  
  !==================================================================================
  Elemental Function Vexternal(t,x,y,z)
    !
    Implicit None
    Integer(ipr), Intent(in) :: t  !! isospin index: 0=isoscalar, 1=isovector
    Real(pr), Intent(in) :: x,y,z  !! position in cartesian basis
    Real(pr) :: Vexternal
    !
    Vexternal = 0.0_pr
    !
  End Function Vexternal
  !
End Module UNEDF
!===================================================================================================================================
!#END UNEDF MODULE
!===================================================================================================================================
!
!===================================================================================================================================
!#START HFBTHO MODULE
!===================================================================================================================================
Module HFBTHO
  Use  HFBTHO_VERSION
  Use  UNEDF  
  !
  ! Input for HFBiterations
  Integer(ipr)  :: n00_INI,iLST_INI,inin_INI,icou_INI
  Integer(ipr)  :: npr_INI(3),kindhfb_INI
  Integer(ipr)  :: keyblo1_INI,keyblo2_INI,IDEBUG_INI
  Real(pr)      :: b0_INI,q_INI,cdef_INI,cqad_INI
  Character(30) :: skyrme_INI
  Real(pr)      :: epsi_INI
  Logical       :: Add_Pairing_INI,Print_HFBTHO_Namelist_INI,DO_FITT_INI                 
  ! Output for regression optimization
  Real(pr)  :: efit_0
  Real(pr), Dimension(0:1) :: efit_rhorho,efit_rhorhoD,efit_rhotau,efit_rhoDrho
  Real(pr), Dimension(0:1) :: efit_rhonablaJ,efit_JJ,efitV0,dfitV0,efV_0
  ! serial output (1:on/0:off)
  Integer(ipr)      :: IDEBUG
  Logical           :: DO_FITT
  ! For loop over used particle types. For normal nuclei min=1, max=2. For n droplets min=max=1.
  Integer(ipr) :: itmin,itmax,irestart
  ! Global numbers
  Real(pr), Parameter :: zero=0.00_pr,half=0.50_pr,one=1.00_pr,two=2.0_pr,three=3.0_pr &
       ,four=4.0_pr,five=5.0_pr,six=6.0_pr,seven=7.0_pr,eight=8.0_pr,nine=9.0_pr,ten=10.0_pr
  ! Whole global numbers pp# 
  Real(pr), Parameter :: pp12=12.0_pr,pp16=16.0_pr,pp15=15.0_pr,pp20=20.0_pr &
       ,pp24=24.0_pr,pp27=27.0_pr,pp32=32.0_pr,pp64=64.0_pr,pp40=40.0_pr
  ! Fractional global numbers p#
  Real(pr), Parameter :: p12=one/two,p13=one/three,p14=0.250_pr,p23=two/three &
       ,p43=four/three,p32=1.50_pr,p34=three/four,p53=five/three,p18=one/eight,p38=three/eight &
       ,p59=five/nine,p52=2.50_pr,p54=five/four,p74=seven/four
  ! Frequent Constants
  Real(pr) :: PI,ffdef3,ffdef4,ffdef5,ffdef6,ffdef7  
  ! Single constants
  Real(pr) :: bet,beta0,q,bp,bpp,bz,hom,hb0,b0,etot,coex,t0s,t0a,drs,dra,ts  &
       ,ta,t3alp,t3al0,t3alm,alp,alm,wla0,wla1,cex,cdef,cqad,ty20,vin,rin  &
       ,ain,qin,pwi,si,siold,epsi,xmix,xmix0,xmax,alst,clst,sklst,alphi,amas  &
       ,skass,varmas,v0ws,akv,hqc,amu,r0,r00,r02,r04,decay,rmm3,amm3,bmm3,cmm3  &
       ,chargee2,EBASECUT
  Real(pr) :: rho_c       
  Integer(ipr) :: lfile,lin,lout,lwin,lwou,lplo,lwel,lres,icstr,icou,ncut  &
       ,iLST1,iLST,maxi,iiter,inin,nzm,nrm,icacou,iqrpa,icacoupj,icahartree &
	   ,nlm,nb,nt,n00,itass,kindhfb,iappend,iError_in_HO,iError_in_THO,ierest
  Integer(ipr) :: n00max=50,esu   
  ! Results 
  Integer(ipr), Parameter :: ieresu=50,ieresl=20,ieresj=50,ieresbl=6 
  Integer(ipr), Parameter :: ieres=ieresu+ieresl+ieresj+ieresbl 
  Real(pr) :: eres(ieres)
  Character(13) :: ereslbl(2)
  Character(2) :: nucname
  Real(pr)    :: eresu(ieresu),eresl(ieresl),eresbl(ieresbl),eresj(ieresj)
  Character(13) :: hlabels(ieres+3)  
  ! Common small arrays
  Real(pr) :: alast(2),ala(2),ala1(2),tz(2),ass(2),drhoi(2),del(2),vso(2),r0v(2)  &
       ,av(2),rso(2),aso(2),Sumnz(2),Dispersion(2),v2min(2),v2minv(2),rms(3),ept(3),q2(3)  &
       ,Dnfactor(3),varmasnz(2),pjmassnz(2)
  Integer(ipr) :: npr(3),inz(2),ldel(2),nk(2),itbl(2),kbl(2),tpar(2),ipbl(2),nbl(2),ibbl(2)  &
       ,klmax(2),inner(2),iasswrong(3),lcc !remove
  ! Lipkin-Nogami
  Real(pr) :: ala2(2),etr(3),ssln(3,2),Geff(2)
  ! Blocking
  Real(pr)     :: pwiblo=2.0_pr, eqpmin(2)=0.0_pr 
  Integer(ipr) :: bloall; Parameter(bloall=200)
  Integer(ipr), Dimension(0:bloall,2) :: bloblo,blo123=0,blok1k2=0
  Real(pr),     Dimension(0:bloall,2) :: bloqpdif
  Integer(ipr) :: iparenti(2),keyblo(3),nkblo_INI(2,5),nkblo(2,5)=0
  Integer(ipr) :: blocross(2),blomax(2),blo123d(2),blok1k2d(2),blocanon(2)
  ! manualBlocking
  Integer(ipr):: manualBlocking=0 
  ! Logical and character variables
  Character(1) :: tq,tp(2),tl(0:20),tis(2)
  Character(30) :: skyrme
  Character(8) :: tit(2)
  Character(7) :: protn(2)
  Data  protn/'neutron','proton '/
  ! Allocatable arrays
  Character(13), Allocatable  :: tb(:)
  Character(25), Allocatable  :: txb(:)
  Real(pr), Allocatable, Target  :: rk(:,:),ak(:,:),hh0(:,:),de0(:,:)  &
       ,ddc(:,:,:),ddc1(:,:,:),qh(:,:),qh1(:,:),ql(:,:,:),ql1(:,:,:)  &
       ,ek(:,:),dk(:,:),vk(:,:),vk1(:,:),uk(:,:),hfb1(:,:),vkmax(:,:)
  Real(pr), Allocatable  :: fdsx(:),fdsy(:),fdsy1(:),fdsy2(:),fdsy3(:),fspb0(:)  &
       ,fspc0(:),fspd0(:),fspb1(:),fspc1(:),fspd1(:),fspb2(:),fspc2(:),fspd2(:),fspb3(:)  &
       ,fspc3(:),fspd3(:),fak(:),fi(:),sq(:),sqi(:),wf(:),wfi(:),rkass(:,:)
  Integer(ipr), Allocatable  :: id(:),ia(:),ikb(:),ipb(:),nz(:),nr(:),nl(:),ns(:),npar(:)  &
       ,ka(:,:),kd(:,:),numax(:,:),iv(:), lcanon(:,:)
  Real(pr), Allocatable  :: AN(:),ANK(:),PFIU(:),PFID(:)
  Real(pr), Allocatable  :: FIU(:),FID(:),FIUR(:),FIDR(:)
  Real(pr), Allocatable  :: FIUD2N(:),FIDD2N(:),FIUZ(:),FIDZ(:)
  ! optimization arrays
  Real(pr), Allocatable  :: QHLA_opt(:,:),FI1R_opt(:,:),FI1Z_opt(:,:),FI2D_opt(:,:),y_opt(:)   
  ! Arrays depending on mesh points
  Integer(ipr)  :: ngh,ngl,nleg,nghl,nbx,ntx,nzx,nrx,nlx,ndx,ndx2,ndxs,nqx
  Integer(ipr)  :: nhfbqx,nb2x,nhfbx,nkx,nzrlx,iqqmax
  Real(pr), Allocatable :: xh(:),wh(:),xl(:),sxl(:),wl(:),xleg(:),wleg(:),vc(:,:)
  Real(pr), Allocatable :: vhbn(:),vn(:),vrn(:),vzn(:),vdn(:),vsn(:),dvn(:)  
  Real(pr), Allocatable :: vhbp(:),vp(:),vrp(:),vzp(:),vdp(:),vsp(:),dvp(:)    
  Real(pr), Allocatable :: vSZFIn(:),vSFIZn(:),vSRFIn(:),vSFIRn(:)
  Real(pr), Allocatable :: vSZFIp(:),vSFIZp(:),vSRFIp(:),vSFIRp(:)
  Real(pr), Allocatable, Target :: aka(:,:),ro(:,:),tau(:,:),dro(:,:),dj(:,:) &
       ,SZFI(:,:),SFIZ(:,:),SRFI(:,:),SFIR(:,:),NABLAR(:,:),NABLAZ(:,:)
  Real(pr), Allocatable  :: fl(:),fli(:),fh(:),fd(:),fp1(:),fp2(:),fp3(:),fp4(:),fp5(:),fp6(:)  &
       ,fs1(:),fs2(:),fs3(:),fs4(:),fs5(:),fs6(:),wdcor(:),wdcori(:),cou(:)
  Real(pr), Allocatable  :: vDHartree(:,:),vhart00(:,:),vhart01(:,:),vhart11(:,:)
  ! PAV Projection
  Integer(ipr)  :: keypj,ilpj,ilpj2,ilnqx,ilnghl
  Real(pr)  :: rehfbcan,ehfb,retotpj,depnp,iproj,npr1pj,npr2pj
  Complex(pr)  :: onei=(0.0_pr,1.0_pr)
  Complex(pr), Allocatable, Target  :: phypj(:),sinphy(:),exp1iphy(:)  &
       ,exp1iphym(:),exp2iphy(:),exp2iphym(:),coupj(:,:),ropj(:,:,:)  &
       ,taupj(:,:,:),dropj(:,:,:),djpj(:,:,:),akapj(:,:,:),pjk(:,:)  &
       ,SZFIpj(:,:,:),SFIZpj(:,:,:),SRFIpj(:,:,:),SFIRpj(:,:,:),epj(:,:)  &
       ,ddepj(:,:,:),cpj(:,:,:),ypj(:,:,:),rpj(:,:,:)
  Real(pr)  :: polem(2),polep(2)
  ! CMC
  Integer(ipr)  :: ICMinput
  Real(pr)     :: ECMHFB(3),ECMPAV(3)
  ! CRC
  Integer(ipr)  :: ICRinput 
  Real(pr) :: DEROT(3),SQUJ(3),CRAN(3),ERIGHFB(3)
  ! hfbdiagonal
  Real(pr), Allocatable :: erhfb(:),drhfb(:),erhfb1(:),drhfb1(:)
  Real(pr), Allocatable :: hfb(:,:),zhfb(:),evvk(:),hfbcan(:,:),evvkcan(:)
  ! Broyden
  Character(1)   ::  bbroyden
  Integer(ipr)  :: nbroyden=7
  Real(pr)     :: alphamix=0.70_pr
  Integer(ipr)  :: nhhdim,nhhdim2,nhhdim3,nhhdim4,ialwork,ilwork
  Real(pr),    Allocatable, Target  :: brout(:),brin(:)
  Real(pr),     Allocatable :: alwork(:)
  Integer(ipr), Allocatable :: lwork(:)
  ! new keys
  Logical :: Parity,Parity_INI  
  Logical :: Print_Screen=.True.
  Logical :: Add_Pairing,Print_HFBTHO_Namelist
  Integer(ipr)  :: MAX_ITER_INI,keypj_INI,iproj_INI,npr1pj_INI,npr2pj_INI
  ! Eqp U,V   
  Integer(ipr) :: nuv,nqp
  Real(pr), Allocatable, Target :: RVqpN(:),RVqpP(:),RUqpN(:),RUqpP(:),REqpN(:),REqpP(:)
  Integer(ipr), Allocatable, Target :: KpwiN(:),KpwiP(:),KqpN(:),KqpP(:)
  ! error indicator
  Integer(ipr)  :: ierror_flag=0
  Character(60) :: ierror_info(0:10)
  ! namelist
  Namelist /HFBTHO_NAMELIST/ MAX_ITER_INI,epsi_INI,Add_Pairing_INI &
       ,icou_INI,iLST_INI,keypj_INI,iproj_INI,npr1pj_INI,npr2pj_INI & 
       ,DO_FITT_INI,IDEBUG_INI,Parity_INI,Print_HFBTHO_Namelist_INI

  !
  ! mpi setup
  integer(ipr),Public :: iam_mpi  !! mpi id number of this process
  Integer, Public     :: num_processors_mpi,num_nodes_mpi,icom_err_mpi  !! for mpi
End Module  HFBTHO
!===================================================================================================================================
!#END HFBTHO MODULE
!===================================================================================================================================
!
!===================================================================================================================================
!#START pairing_HFBTHO MODULE
!===================================================================================================================================
Module pairing_HFBTHO 
  Use HFBTHO 
  Implicit None 
  Character(6) :: pairing_Version='9' 
  ! Version History 
  !=================================================================================================================================  
  !Integer(ipr)             :: nlimtmr
  Integer(ipr)              :: NNN0,mv,ivpair=0 
  Real(pr)                  :: gl0,gal,gl(2),vfac(2),e_pair(2),d_pair(2)
  Real(pr),     Allocatable :: wnn(:,:),bin(:,:),rk_pair(:,:),ak_pair(:,:)
  Integer(ipr), Allocatable :: jsort(:)
  !=================================================================================================================================
  !
Contains 
  !
  !=================================================================================================================================
  Subroutine vpair
    !-----------------------------------------------------------------------
    ! calculates matrix elements V_nz and V_nr' for TMR-pairing in the axially deformed oscillator basis
    !-----------------------------------------------------------------------
    Implicit None
    Real(pr), Allocatable :: v12z(:,:),v12p(:,:),pnosc(:),wn(:)
	Real(pr)     :: eps,s,smax
	Integer(ipr) :: i,i1,i2,j,j0,j1,ib,nd,n1,n2,nn,nnz,nnr
    Data eps/1.0d-8/
    !----------------------------------------------------------------------
	! TMR pairing input use_TMR_pairing=1
    !----------------------------------------------------------------------
	vfac(1)=CpV1(0);     vfac(2)=CpV1(1) !  1.100d0     
    gl0 = CpV0(0)                        ! -728.0d0
    gal = CpV0(1)                        !  0.6442049360d0
	gl = vfac * gl0                      !  nlimtmr=Infty (was 8)
    Write(*,'(a,2g16.4)') '  TMR pairing: Tian,Ma,Ring, PRB 676, 44 (2009)'
    Write(*,'(a,2i6)'   ) '  TMR pairing        nlimtmr =:   No Limit'
    Write(*,'(a,2g16.8)') '  TMR pairing strength     G =: ',gl0
    Write(*,'(a,2g16.8)') '  TMR pairing width    a,a^2 =: ',gal,gal**2
    Write(*,'(a,2g16.8)') '  TMR pairing reduction vfac =: ',vfac
    Write(*,'(a,2g16.8)') '  TMR pairing strength    gl =: ',gl
    !----------------------------------------------------------------------
    ! determines mv value (and quantum numbers for output)
    !----------------------------------------------------------------------
    i = 0
    Do ib = 1,nb
       nd = id(ib)
       Do n2   = 1,nd
          Do n1 = n2,nd
             i  = i + 1        
          Enddo   ! n1
       Enddo   ! n2
    Enddo   ! ib
    mv = i
    If(Allocated(v12z)) Deallocate(v12z,v12p,pnosc,wn,rk_pair,ak_pair,jsort)
    Allocate(v12z(1:MV,0:N00),v12p(1:MV,0:N00),pnosc(0:N00),wn(MV),rk_pair(MV,2),ak_pair(MV,2),jsort(MV))
    !----------------------------------------------------------------------
	! calculate binomial coefficients
    !----------------------------------------------------------------------
	Call  binom
    !----------------------------------------------------------------------
    ! calculate V12Z and V12P
    !----------------------------------------------------------------------
    Call v12nz(v12z,pnosc)
    Call v12np(v12p,pnosc)
    !----------------------------------------------------------------------
    ! calculate the size NNN0 of the s.p. matrix elements wnn, alocate it 
    !----------------------------------------------------------------------
    nn = 0
    Do nnz = 0,n00
       Do nnr = 0,n00
          ! If (2*nnz+2*nnr.Gt.nlimtmr) Cycle
          nn = nn + 1
          smax = zero
          Do i = 1,mv
             smax = Max(smax,Abs(v12z(i,nnz)*v12p(i,nnr)))
          Enddo   
          If (smax.Lt.eps) nn = nn - 1
       Enddo   
    Enddo   
    NNN0 = nn; If(Allocated(wnn)) Deallocate(wnn); Allocate(wnn(MV,NNN0))
    !----------------------------------------------------------------------
    ! calculate the single particle matrix elements wnn
    !----------------------------------------------------------------------
    nn = 0
    Do nnz = 0,n00
       Do nnr = 0,n00
          !If (2*nnz+2*nnr.Gt.nlimtmr) Cycle
          nn = nn + 1
          smax = zero
          Do i = 1,mv
             wn(i) = v12z(i,nnz)*v12p(i,nnr)
             smax = Max(smax,Abs(wn(i)))
          Enddo  
          If (smax.Lt.eps) Then
             nn = nn - 1
          Else
             Do i = 1,mv
                wnn(i,nn) = wn(i)
             Enddo   
          Endif  ! smax < eps
       Enddo   
    Enddo   
    If(nn.Ne.NNN0) Stop 'Something wrong in vpair!'
    !----------------------------------------------------------------------
    ! Reordering array jpint
    !----------------------------------------------------------------------
	i = 0; j0=0; 
    Do ib = 1,nb
	 nd = id(ib); i1=0
	 Do n2 = 1,nd
	  i1=i1+n2; i2=n2-2
	  Do n1 = n2,nd
	   i2=i2+1
	   If(n1.Eq.n2) Then
	    j1=i1
	   Else
	    j1=j1+i2
	   Endif
	   j=j0+j1
	   i = i + 1
	   jsort(i)=j
	  Enddo   ! n1
     Enddo   ! n2
	 j0=j
    Enddo   ! ib
	!
    If(IDEBUG.Gt.10) Write(6,'(10(a,1x,i6,2x))')  ' VPAIR: mv =', mv,'NNN0 = ',NNN0
	!
	Deallocate(v12z,v12p,pnosc,wn,bin)
	!
  End  Subroutine vpair
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine v12nz(v12z,pnosc)
    !-----------------------------------------------------------------------
    ! calculates matrix elements V_nz for TMR-pairing
    !-----------------------------------------------------------------------
    Implicit None
	Real(pr)     :: eps
    Real(pr)     :: v12z(1:MV,0:N00),pnosc(0:N00)	
    Integer(ipr) :: i,il,ib,nd,i0,n2,nz2,n1,nz1,nn,n,nh,nnh,nnzm,nx 
    Data eps/1.0d-6/
    !
    Call get_CPU_time('v12nz',0)	
    !----------------------------------------------------------------------
    !    calculate V12NZ 
    !----------------------------------------------------------------------
    Call pnoscz(n00,gal,bz,pnosc)
	v12z=zero
    il = 0
    Do ib = 1,nb
       nd = id(ib)
       i0= ia(ib)
       Do n2   = 1,nd
          nz2 = nz(i0+n2)
          Do n1 = n2,nd
             nz1 = nz(i0+n1)
             il  = il + 1        
             If (Mod(nz1+nz2,2).Ne.0) Cycle
             Do nn = 0,nz1+nz2,2
                n  = nz1 + nz2 - nn
                nh  = n/2
                nnh = nn/2
                If (nh .Gt.n00) Stop 'in VPAIR: n too large'
                If (nnh.Gt.n00) Stop 'in VPAIR: nn too large'
                v12z(il,nnh) = pnosc(nh)*talmos1(nz1,nz2,nn,n)
             Enddo   ! nn 	
          Enddo   ! n1
       Enddo   ! n2
    Enddo   ! ib
    !----------------------------------------------------------------------
    ! calculate the maximal nnz:  nnzm 
    !----------------------------------------------------------------------
    nnzm =0
    Do nn = 0,n00
       nx = 0
       Do i = 1,mv
          If (Abs(v12z(i,nn)).Gt.eps) nx = nn
       Enddo   ! i
       nnzm = Max(nnzm,nx)
    Enddo   ! nn
    If(nnzm.Gt.n00) Stop 'in V12NZ: nnzm too large'
    !
    If(IDEBUG.Gt.10) Write(6,*) 'nnzm =',nnzm
    !
    Call get_CPU_time('v12nz',1)	
	!
  End Subroutine v12nz
  !=================================================================================================================================
  !  
  !=================================================================================================================================
  Subroutine v12np(v12p,pnosc)  
    !-----------------------------------------------------------------------
    !     matrix elements V_np (perpendicular direction) for TMR-pairing
    !-----------------------------------------------------------------------
    Implicit None
	Real(pr)     :: eps
    Real(pr)     :: v12p(1:MV,0:N00),pnosc(0:N00)	
	Integer(ipr) :: i,il,ib,nd,i0,n,nn,n1,nr1,ml1,n2,nr2,ml2,nnrm,nx  
    Data eps/1.0d-10/
    !
    Call get_CPU_time('v12np',0)	
	!
    !----------------------------------------------------------------------
    !    calculate V12NP
    !----------------------------------------------------------------------
    Call pnoscp(n00,gal,bp,pnosc)
	v12p=zero
    il = 0
    Do ib = 1,nb
       nd = id(ib); i0= ia(ib)
       Do n2 = 1,nd
          nr2 = nr(i0+n2); ml2 = nl(i0+n2)
          Do n1 = n2,nd
             nr1 = nr(i0+n1); ml1 = nl(i0+n1)
             il  = il + 1        
             If (ml1.Ne.ml2) Cycle
             Do nn = 0,nr1+nr2+ml1
                n = nr1+nr2+ml1-nn
                If (nn.Gt.n00) Stop 'in VPAIR: nn too large'
                If (n .Gt.n00) Stop 'in VPAIR: n  too large'
                v12p(il,nn) = pnosc(n)*talmos20(ml1,nr1,nr2,0,nn,n) 
             Enddo   ! nn
          Enddo   ! n1
       Enddo   ! n2
    Enddo   ! ib
    !----------------------------------------------------------------------
    !---- calculate the maximal nnr:  nnrm 
    !----------------------------------------------------------------------
    nnrm = 0
    Do nn = 0,n00
       nx = 0
       Do i = 1,mv
          If (Abs(v12p(i,nn)).Gt.eps) nx = nn
       Enddo   ! i
       nnrm = Max(nnrm,nx)
    Enddo   ! nn
    If(nnrm.Gt.n00) Stop 'in V12NP: nnrm too large'
    !
    If(IDEBUG.Gt.10) Write(6,*) 'nnrm =',nnrm
    !
    Call get_CPU_time('v12np',1)	
	!
  End Subroutine v12np
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine pnoscz(nm,a,b,pnosc)	    
    Implicit None
	Integer(IPR) :: nm,nh,n
    Real(pr)     :: a,b,a2,b2,alpha,br,s0,s1
    Real(pr)     :: pnosc(0:N00)
    !
    br = b*sq(2); a2=a*a; b2=b*b
    s0 = (two*pi)**(-0.25d0)*Sqrt(b/(a2+b2))
    s1 = (a2-b2)/(a2+b2)
    Do nh = 0,nm
       n = nh + nh
       pnosc(nh) = s0 * s1**nh * wf(n) / 2**nh * fi(nh)
    Enddo    ! n
    If(IDEBUG.Gt.10) Write(6,*) 'PNZ',nm,pnosc(0:nm)
  End Subroutine pnoscz
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine pnoscp(nm,a,b,pnosc)
    Implicit None
	Integer(IPR) :: i,n,nm
    Real(pr)     :: a,b,br,alpha,s,x,y,al2
    Real(pr)     :: pnosc(0:N00)			
    !
    If(nm.Gt.n00) Stop ' in PNOSCP: nm too large'
    !
    br = b*sq(2); alpha = a/br; al2=alpha*alpha
    Do n = 0,nm
       s = zero
       Do i = 1,NGL
          x = sxl(i); y = xl(i)/(four*al2)
          s = s + osc2(n,0,x)*Exp(-y)*wl(i)   
       Enddo  
       pnosc(n) = s/(four*br*al2*Sqrt(two*pi))
    Enddo    ! n
    If(IDEBUG.Gt.10)  Write(6,*) 'PNP',nm,pnosc(0:nm)	
  End Subroutine pnoscp
  !=================================================================================================================================
  ! 
  !=================================================================================================================================
  Function talmos1(n1,n2,n3,n4)
    !-----------------------------------------------------------------------
    !	1-dimensional Talmi-Moshinsky bracket: <n1 , n2 | n3 , n4 >
    !	quantum number n_i start from zero: n = 0,1,2,....	
    !-----------------------------------------------------------------------
    Implicit None
	Real(pr)     :: talmos1,f,s
	Integer(ipr) :: n1,n2,n3,n4,m4,m3,n14
	!
    talmos1 = zero
    If (n1+n2.Ne.n3+n4) Return
    f = wf(n1)*wfi(n3)*wf(n2)*wfi(n4)/Sqrt(two**(n3+n4))
    s = zero; n14=n1-n4
    Do m4 = 0,n4
       m3 = n14 + m4
       If(m3.Ge.0.And.n3.Ge.m3) s = s + iv(m4)*bin(n3,m3)*bin(n4,m4)
    Enddo
    talmos1 = f*s
  End Function talmos1
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Function talmos2(im1,in1,in2,im3,in3,in4)
    !-----------------------------------------------------------------------
    !	2d-moshinsky bracket:	
    !	<n1 m1, n2 -m1 | n3 m3, n4 -m3>
    !	radial quantum number start from zero: n=0,1,2,....	
    !	orbital angular momentum m1=-m2>0,m3=-m4>0
    !-----------------------------------------------------------------------
    !     iv(n)  =  (-1)**n
    !     fak(n) =  n!
    !     fi(n)  =  1/n!
    !     wf(n)  =  sqrt(n!)
    !     wfi(n) =  1/sqrt(n!)
    !-----------------------------------------------------------------------
    Implicit None
    !
	Real(pr)     :: talmos2,s34,sn3,sn4,sm3,sm4,prout
	Integer(ipr) :: n1,n2,n3,n4,nn1,nn2,nn3,nn4,m1,m2,m3,m4,im1,in1,in2,im3,in3,in4,im2,im4,nn12,nn34
	Integer(ipr) :: i3,i4,j3,j4,k3,k4,l3,l4,it3,it4	
    !
    im2 = -im1; im4 = -im3; talmos2 = zero
    nn12 = 2*in1+iabs(im1)+2*in2+iabs(im2)
    nn34 = 2*in3+iabs(im3)+2*in4+iabs(im4)
    If (im1.Lt.0)     Stop 'in TALMOS2: m1 < 0'
    If (im3.Lt.0)     Stop 'in TALMOS2: m3 < 0'
    If (nn12.Ne.nn34) Return
    s34 = zero
    If (im1.Lt.im2) Then
       n1 = in2; n2 = in1; m1 = im2; m2 = im1
    Else
       n1 = in1; n2 = in2; m1 = im1; m2 = im2
    Endif
    n3 = in3; n4 = in4
    m3 = im3; m4 = im4
    nn1 = 2*n1 + Abs(m1)
    nn2 = 2*n2 + Abs(m2)
    nn3 = 2*n3 + Abs(m3)
    nn4 = 2*n4 + Abs(m4)
    If (m3.Gt.m4) Then
       If (n1+n2.Ne.n3+n4+m2-m4) Return
    Else
       If (n1+n2.Ne.n3+n4+m2-m3) Return
    Endif
    prout = iv(n3+n4-n1-n2)/Sqrt(two**(nn3+nn4))* wf(n1)*wf(n1+Abs(m1))* &
	         wf(n2)*wf(n2+Abs(m2))*wfi(n3)*wfi(n3+Abs(m3))*wfi(n4)*wfi(n4+Abs(m4))
    sn3 = zero; sn4 = zero; sm3 = zero; sm4 = zero
    Do i3 = 0,n3
       Do j3 = 0,n3
          Do k3 = 0,n3
             l3 = n3 - i3 - j3 - k3
             Do i4 = 0,n4
                Do j4 = 0,n4
                   Do k4 = 0,n4
                      l4 = n4 - i4 - j4 - k4
                      If (l3.Lt.0.Or.l4.Lt.0) Cycle ! k2
                      Do it3 = 0,Abs(m3)
                         Do it4 = 0,Abs(m4)
                            If (m3.Gt.m4) Then
                               If (i3+i4+j3+j4+it3.Ne.n2) Cycle
                               If (j3+j4.Ne.m2+k3+k4-it3+it4) Cycle
                               If (l3+l4.Ne.n3+n4-n2-k3-k4+it3) Cycle
                               sn3 = fak(n3)*fi(l3)*fi(i3)*fi(j3)*fi(k3)
                               sn4 = iv(j4+k4)*fak(n4)*fi(l4)*fi(i4)*fi(j4)*fi(k4)
                               sm3 = fak(Abs(m3))*fi(it3)*fi(Abs(m3)-it3)
                               sm4 = iv(it4)*fak(Abs(m4))*fi(it4)*fi(Abs(m4)-it4)
                               s34 = s34 + sn3*sn4*sm3*sm4
                            Else
                               If (i3+i4+j3+j4+it4.Ne.n2) Cycle
                               If (j3+j4.Ne.m2+k3+k4+it3-it4) Cycle
                               If (l3+l4.Ne.n3+n4-n2-k3-k4+it4) Cycle
                               sn3 = fak(n3)*fi(l3)*fi(i3)*fi(j3)*fi(k3)
                               sn4 = iv(j4+k4)*fak(n4)*fi(l4)*fi(i4)*fi(j4)*fi(k4)
                               sm3 = fak(Abs(m3))*fi(it3)*fi(Abs(m3)-it3)
                               sm4 = iv(it4)*fak(Abs(m4))*fi(it4)*fi(Abs(m4)-it4)
                               s34 = s34 + sn3*sn4*sm3*sm4
                            Endif ! m3,m4
                         Enddo   ! it4
                        Enddo   ! it3
                      Enddo   ! k4
                   Enddo   ! j2
                Enddo   ! i2
             Enddo   ! k1
          Enddo   ! j1
       Enddo   ! i1
       talmos2 = s34*prout
     End Function talmos2
     !=================================================================================================================================
     !
     !=================================================================================================================================
     Function talmos20(im1,in1,in2,im3,in3,in4)
      !-----------------------------------------------------------------------
      !	2d-moshinsky bracket:	(CASE im3=0)
      !	<n1 m1, n2 -m1 | n3 m3, n4 -m3>
      !	radial quantum number start from zero: n=0,1,2,....	
      !	orbital angular momentum m1=-m2>0,m3=-m4>0
      !-----------------------------------------------------------------------
      !     iv(n)  =  (-1)**n
      !     fak(n) =  n!
      !     fi(n)  =  1/n!
      !     wf(n)  =  sqrt(n!)
      !     wfi(n) =  1/sqrt(n!)
      !-----------------------------------------------------------------------
      Implicit None
      !
	  Real(pr)     :: talmos20,s34,sn3,sn4,sm3,sm4,prout,wfi34
	  Integer(ipr) :: n1,n2,n3,n4,nn1,nn2,nn3,nn4,m1,m2,m3,m4,im1,in1,in2,im3,in3,in4,im2,im4,nn12,nn34,m22
	  Integer(ipr) :: i3,i4,j3,j4,k3,k4,l3,l4,it3,it4	
      !
      If(im1.Lt.0)     Stop 'in TALMOS20: m1 < 0'
      im2 = -im1; 
	  talmos20 = zero
      nn12 = in1+in1+Abs(im1)+in2+in2+Abs(im2)
      nn34 = in3+in3+in4+in4
      If(nn12.Ne.nn34) Return
      If (im1.Lt.im2) Then
       n1 = in2; n2 = in1; m1 = im2; m2 = im1
      Else
       n1 = in1; n2 = in2; m1 = im1; m2 = im2
      Endif
      n3 = in3; n4 = in4
	  If (n1+n2.Ne.n3+n4+m2) Return
	  m1=Abs(m1); m22=m2; m2=Abs(m2)
      nn1 = n1 + n1 + m1
      nn2 = n2 + n2 + m2
      nn3 = n3 + n3 
      nn4 = n4 + n4 
	  wfi34=(wfi(n3)*wfi(n4))**2
      prout = iv(n3+n4-n1-n2)*wf(n1)*wf(n1+m1)*wf(n2)*wf(n2+m2)*wfi34/Sqrt(two**(nn3+nn4))
      s34 = zero
      Do i3 = 0,n3
       Do j3 = 0,n3
		Do k3 = 0,n3
		 l3 = n3-i3-j3-k3
		 If (l3.Lt.0) Cycle 
		 sn3 = fak(n3)*fi(l3)*fi(i3)*fi(j3)*fi(k3)			 
		 Do i4 = 0,n4
		  j4=n2-i3-i4-j3
		  k4=j3+j4-m22-k3
		  l4 = n4-i4-j4-k4
		  If (l4.Lt.0) Cycle 
		  sn4 = iv(j4+k4)*fak(n4)*fi(l4)*fi(i4)*fi(j4)*fi(k4)
		  s34 = s34 + sn3*sn4
		 Enddo ! i4
		Enddo ! k3
	   Enddo ! j3
	  Enddo ! i3
	  talmos20 = s34*prout
     End Function talmos20
     !=================================================================================================================================
     !
     !=================================================================================================================================
     Subroutine binom
      !---------------------------------------------------------------------
      ! THE ARRAY OF BINOMIAL COEFFICIENTS BIN(I,J)= = I!/J!/(I-J)! 
      !---------------------------------------------------------------------
	  Implicit None
      Integer(ipr) :: i,k,igfv
      Parameter(igfv=170)               !maximal number for GFV
	  !
      If(Allocated(bin)) Deallocate(bin)
      Allocate(bin(0:igfv,0:igfv))
	  bin = zero
      Do i = 0,IGFV
         bin(i,0)=one; bin(i,i)=one
         Do k = 1,i/2
            bin(i,k)   = dnint(bin(i,k-1)/dfloat(k)*dfloat(i-k+1))
            bin(i,i-k) = bin(i,k)
         Enddo
      Enddo
     End Subroutine binom
     !=================================================================================================================================
     !
     !=================================================================================================================================
      Function osc2(n,m,r)
      !---------------------------------------------------------------------
      ! calculates the radial wavefunctions for the zylindrical oscillator
      ! osc2(n,m,r) = N_(n,m)*sqrt(2) *  r^|m| * L_n^|m|(r*r) * exp(-r*r/2)
      ! N_(n,m) = sqrt( n! / (n+|m|)! )
      !---------------------------------------------------------------------
      Implicit None
      Real(pr)     :: r,x,w0,s0,s1,s2,osc2
	  Integer(ipr) :: n,m,i
      !
      If (m.Lt.0) Stop 'in OSC2:  m < 0'
      !
      x = r*r;        w0 = sq(2)*Exp(-half*x) * r**m
      s0 = w0*wfi(m); s1 = zero
      Do i = 1,n
         s2 = s1; s1 = s0
         s0 = ((i+i+m-1-x)*s1-sq(i-1)*sq(i-1+m)*s2)*sqi(i)*sqi(i+m)
      Enddo   ! i
      osc2 = s0 
      End Function osc2
     !=================================================================================================================================
     !	 
     !=================================================================================================================================
      Subroutine delta
      !---------------------------------------------------------------------
	  ! pairing in oscillator basis
      ! ivpair = 0   initializing
      !        > 0   during the iteration
      !---------------------------------------------------------------------
      Implicit None
      Integer(ipr)          :: i,j,it,nn,ib,nd,m,nq,il,n1,n2,i1,i2,ne,nnr,n12,n21
	  Real(pr)              :: sk,sr,g
      Real(pr), Allocatable :: Pknn(:,:),Prnn(:,:)	
      !---------------------------------------------------------------------
	  ! Calculate only once wnn
      !---------------------------------------------------------------------
	  If(ivpair.Eq.0) Then
	    Call vpair
		ivpair=1
	  Endif		  
      !
      Call get_CPU_time('delta',0)	
	  !
      If(Allocated(Pknn)) Deallocate(Pknn,Prnn); Allocate(Pknn(NNN0,2),Prnn(NNN0,2))
	  !
      !---------------------------------------------------------------------
	  ! HFBTHO densities in HO space	
      !---------------------------------------------------------------------
	  i=0; 
	  Do ib=1,nb
	   nd=id(ib); m=ib+(it-1)*nbx	
	   Do n2=1,nd
        Do n1=n2,nd
	     i=i+1
	     n12=n1+(n2-1)*nd  ! n21=n2+(n1-1)*nd
		 ak_pair(i,1)=ak(n12,ib)
		 ak_pair(i,2)=ak(n12,ib+nbx)
		 rk_pair(i,1)=rk(n12,ib)
		 rk_pair(i,2)=rk(n12,ib+nbx)
	    End Do !n1
       End Do !n2
	  End Do !ib
	  !
      Do it = 1,2
       g = 0.50d0*gl(it)
       !---------------------------------------------------------------------
	   !calculation of knn,rnn
       !---------------------------------------------------------------------
       Do nn = 1,NNN0
         sk  = zero; sr  = zero 
         Do i = 1,mv
            sk = sk + wnn(i,nn) * ak_pair(i,it)
            sr = sr + wnn(i,nn) * rk_pair(i,it)
         Enddo   ! i
         Pknn(nn,it) = sk
         Prnn(nn,it) = sr
       Enddo   ! nn
       !---------------------------------------------------------------------
       !calculation of the pairing energy and avarage delta
       !---------------------------------------------------------------------
       sk = zero; sr = zero 
       Do nn = 1,NNN0
         sk = sk + Pknn(nn,it)**2
         sr = sr + Pknn(nn,it)*Prnn(nn,it)
       Enddo   ! i
	   d_pair(it) = -g*sr
       e_pair(it) =  g*sk
       !---------------------------------------------------------------------
	   !calculation of pairing matrix delta
       !---------------------------------------------------------------------
	   i=0
       Do ib = 1,nb
	    nd = id(ib)
        Do n2 = 1,nd
         Do n1 = n2,nd
          i = i + 1
          sk = zero
          Do nn = 1,NNN0
             sk = sk + wnn(i,nn)*Pknn(nn,it)
          Enddo   ! nn
          sk =  g*sk
		  If(it.Eq.1) Then
		   i1=jsort(i)+nhhdim2
		   !write(*,*) brout(i1),sk,brout(i1)-sk; pause
		   brout(i1)=sk 
		  Else
		   i2=jsort(i)+nhhdim3
		   brout(i2)=sk 
		  Endif
         Enddo   ! n1
        Enddo   ! n2
       Enddo   ! ib
	  Enddo   ! it
	  !
	  Deallocate(Pknn,Prnn)
	  !
      Call get_CPU_time('delta',1)	
      !
      End Subroutine delta
     !=================================================================================================================================
     !
End Module pairing_HFBTHO 
!===================================================================================================================================
!#END pairing_HFBTHO MODULE
!===================================================================================================================================
!
#ifndef hide_tho
!===================================================================================================================================
!#START THO_MODULE MODULE
!===================================================================================================================================
Module THO_MODULE 
  Use HFBTHO, Only: pr,ipr,nghl,nzrlx,fh,fl,fli,fd,fp1,fp2,fp3,fp4,fp5,fp6,fs1,fs2,fs3,fs4,fs5,fs6, &
       wdcor,wdcori,one,bp,bz,ngl,xl,ngh,xh,ilst,iasswrong,half,two,four,six,three,five,pi,wh,wl,lout, &
       lfile,ilst1,zero,fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,fspd0,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2, &
       fspb3,fspc3,fspd3,bmm3,iqqmax,rmm3,amm3,cmm3,bpp,iv,fak,fi,sq,nb,id,ia,nz,nr,nl,qh,ql,ns,nbx,rk, &
       qh,sq,nzm,sqi,nlm,ql,wfi,nrm,b0,itmin,itmax,itass,decay,ass,rms,npr,hb0,dnfactor, &
	   ierror_flag,ierror_info  
  Implicit None 
  !
Contains 
  !  
  !=======================================================================
  Subroutine f01234(lpr)
    !---------------------------------------------------------------------
    ! calculates lst-function 'f' its derivatives 'f1,f2,f3,f4',
    ! the jacobian 'fd=(f^2 f1/r^2)^(1/2)' and its first and second
    ! derivatives 'fd1' and 'fd2' at the point 'r'.
    ! The lst-function is defined as in 'thofun'
    !
    ! All integrations are performed in the dimensionless variables
    ! 'u=xh(ih)' and 'v=xl(il)' being the dimensionless gausss mesh
    ! points in cylindrical coordinates. Thus 'r' is found as an invers
    ! function of the input variable 'qq=Sqrt[xh(ih)^2+xl(il)]', i.e.,
    ! solving f(r)=qq. The original 'z' and 'rho' variables in 'fm' are
    ! defined as 'fh(ihli)=zz*bz' and 'fl(ihli)=rr*bp'.
    !
    ! All required quantities for evaluating w.f., its first derivative
    ! and the associated laplasian are calculated
    ! -------------------------------------------------
    ! definitions used:
    ! ' xh(ih),xl(il)'  gauss mesh points
    ! ' bp, bz       '  oscillator param. in fm
    ! ' iLST1       '  0->ho, #0->tho
    !-----
    ! in this way:
    ! (1) \delta \rho devided by j^2 is:
    ! dro=(cz*fs1+cr*fs2+fs3)*qhab*qlab
    !+ fs1*qha1b1*qlab
    !+ fs2*qhab*qla1b1/(4.*v**2)
    !+ fs4*(qha1b+qhab1)*(qla1b+qlab1)/(2.*v)
    !+ fs5*(qha1b+qhab1)*qlab
    !+ fs6*qhab*(qla1b+qlab1)/(2.*v)
    ! with
    ! cr=1./4.-(nr_a+nr_b+m+1)/(2.*v)+(m/(2.*v))**2
    ! cz=u*u-(nz_a+nz_b+1)
    ! (2) the first (r,z) derivatives are:
    ! fi'_z=fp1*(qh*ql)+fp2*(qh1*ql)+fp3*(qh*ql1)/(2.*v)->HO->qh1*ql/bpz
    ! fi'_r=fp4*(qh*ql)+fp5*(qh1*ql)+fp6*(qh*ql1)/(2.*v)->HO->qh*ql1/Sqrt(v)/bp
    !
    !---------------------------------------------------------------------
    Implicit None
    Logical :: lpr
    Integer(ipr), Save :: i,il,ih,ihli,iw,key0=0,key1=1
    Real(pr), Save :: bri,bri2,bzi,bzi2,wv1,u,qq,f,f1,f2,f3,rhoi,fd1,fd2,  &
         fd12,fdd,r,r2,r3,r4,r5,r6,rr,rr2,rr4,zz,zz2,drr1,drr2,drz1,  &
         drz2,drr12,drz12,g,g1,g2,g3,gg,gg1,gg2,g1g1,uz1,ur1,vz1,vr1,  &
         uz2,ur2,vz2,vr2,ur12,vr12,uz12,vz12
    !
    bri=one/bp; bri2=bri*bri; bzi=one/bz; bzi2=bzi*bzi
    Do il=1,ngl
       wv1=xl(il)
       Do ih=1,ngh
          ihli=ih+(il-1)*ngh
          u=xh(ih); qq=Sqrt(u*u+wv1)
          If(iLST1.Eq.0) Then
             ! ho-case
             r=qq; f=r; f1=one; f2=zero; f3=zero
          Else
             ! tho-case: initial run
             If(ih*il.Eq.1.And.iLST.Lt.0)  &
                  Call thofun(key0,g,f,f1,f2,f3,g1,.True.)  !.false.) !lpr)
             If(iasswrong(3).Ne.0) Then
                ! reinforce ho results
                r=qq; f=r; f1=one; f2=zero; f3=zero
             Else
                ! tho-case: f(r)=qq,f'(r),f''(r),f'''(r), r=Invers_f(r)
                Call thofun(key1,qq,f,f1,f2,f3,r,.False.)
             End If
          End If
          ! Jacobian calculations
          r2=r*r; r3=r2*r; r4=r2*r2 ; r5=r3*r2 ; r6=r3*r3
          ! fdd=(f(r)^2 f'(r)/r^2)^(1/2),fd1=fdd'/fdd, fd2=fdd''/fdd
          fd1=f1/f-one/r+half*f2/f1
          fdd=f*Sqrt(f1)/r
          fd2=two*(f2/f-fd1/r)-(f2/f1)**2/four+half*f3/f1
          fd12=fd1**2
          ! g=(f/r)-derivatives
          g=f/r; g1=-(f-f1*r)/r2
          g2=(two*(f-f1*r)+f2*r2)/r3
          g3=(six*(f1*r-f)-three*f2*r2+f3*r3)/r4
          ! g4=(24.0_pr*(f-f1*r+half*f2*r2)-four*f3*r3+f4*r4)/r5
          gg=g*g; gg1=g*g1; gg2=g*g2; g1g1=g1*g1
          ! (rr,zz)-definitions
          rr=Sqrt(wv1)/g; rhoi=bri/rr; zz=u/g
          rr2=rr*rr; rr4=rr2*rr2 ; zz2=zz*zz
          ! (r,z)-derivatives
          drr1=bri*rr/r;           drz1=bzi*zz/r
          drr2=(bri2-drr1**2)/r; drz2=(bzi2-drz1**2)/r
          drr12=drr1**2;            drz12=drz1**2
          ! (u,v)-derivatives
          uz1=bzi*(g+g1*zz2/r);    ur1=bri*g1*rr*zz/r
          vz1=bzi*two*gg1*rr2*zz/r
          vr1=bri*two*rr*(gg+gg1*rr2/r)
          uz2=bzi2*zz*(three*g1*r2-g1*zz2+g2*r*zz2)/r3
          ur2=bri2*(g1*r2-g1*rr2+g2*r*rr2)*zz/r3
          vz2=bzi2*two*rr2*(gg1*r2-gg1*zz2+g1g1*r*zz2+gg2*r*zz2)/r3
          vr2=gg*r3+five*gg1*r2*rr2-gg1*rr4+g1g1*r*rr4+gg2*r*rr4
          vr2=vr2*bri2*two/r3; ur12=ur1**2; vr12=vr1**2; uz12=uz1**2
          vz12=vz1**2
          ! storage
          fh(ihli)=zz*bz; fl(ihli)=rr*bp; fli(ihli)=one/fl(ihli); fd(ihli)=fdd*fdd
          ! for first derivatives
          fp1(ihli)=fd1*drz1; fp2(ihli)=uz1; fp3(ihli)=vz1
          fp4(ihli)=fd1*drr1; fp5(ihli)=ur1; fp6(ihli)=vr1
          ! for the laplasian
          fs1(ihli)=two*(ur12+uz12); fs2(ihli)=two*(vr12+vz12)
          fs3(ihli)=two*(fd1*(drr1*rhoi+drr2+drz2)+  &
               (fd12+fd2)*(drr12+drz12))
          fs4(ihli)=two*(ur1*vr1+uz1*vz1)
          fs5(ihli)=four*fd1*(drr1*ur1+drz1*uz1)+ &
               ur1*rhoi+ur2+uz2
          fs6(ihli)=four*fd1*(drr1*vr1+drz1*vz1)+ &
               vr1*rhoi+vr2+vz2-(vr12+vz12)/wv1
       End Do   !ihs
    End Do   !il
    !
    ! Associated (z,r)-weights
    Do il=1,ngl
       Do ih=1,ngh
          i=ih+(il-1)*ngh
          wdcor(i)=pi*wh(ih)*wl(il)*bz*bp*bp/fd(i)
          wdcori(i)=one/wdcor(i)
       End Do
    End Do
    !
    If(lpr) Then
       Do iw=lout,lfile
          Write(iw,*)
          If(iLST1.Eq.0) Then
             Write(iw,*) ' ### HO case: wdcor charged'
          Else
             Write(iw,*) ' ### THO case: wdcor charged'
          End If
          Write(iw,*)
       End Do
    End If
    Return
  End Subroutine f01234
  !=======================================================================
  !
  !=======================================================================
  Subroutine thofun(key,r,f,f1,f2,f3,fj,lpr)
    !---------------------------------------------------------------------
    ! Calculates LST-function 'f' its derivatives 'f1,f2,f3'
    ! at the point 'r' (all dimensionless).
    !---------------------------------------------------------------------
    Implicit None
    Logical :: lpr
    Integer(ipr) :: key,msw
    Integer(ipr) :: it,iter,ir,iqq,irmax,irmsit,immho,imm1,  &
         imm2,immm,immmax,imm3
    Real(pr) ::  r,f,f1,f2,f3,fj
    Real(pr), Allocatable :: dsx(:),dsy(:),dsyT(:),dsyi(:),dsyii(:),dsy1(:),  &
         dsy1i(:),spb0(:),spc0(:),spd0(:),spbi(:),spci(:),spdi(:)
    Real(pr) :: h,hhb,pihhb,c00,snorm,snorm1,assm,asm1,asm2,asm3,  &
         rmsit,rmmho,z1,z10,aaa,bex,rend,fj1,fj2,fj3,  &
         s,s1,sN,sP,sT,qq,qqup,qqdn,zqq,zqqi,df,zfj1,zfj1i,fjb,aa,bb,yyy,  &
         sqsq,rmmm,rmmmb0,z1mmm,rmm1,rmm1b0,z1mm1,rmm2,rmm2b0,z1mm2,  &
         rmmmax,z1mmmax,rmmmaxb0,rmmx,z1mmx,z1mmxx,alaex,aldsy1,decay2
    Real(pr) :: denm1(2),denm2(2),rdenm(2)
    Real(pr) :: epsf=1.0d-14,epsdsy=1.0d-16
    Complex(pr) :: yyy1,bbb1,aac !for the e3rd order equation
    !
    !===========================
    ! KEY=0 INITIAL CALCULATIONS
    !===========================
    qq=r
    If(key.Eq.0.And.iLST.Le.0) Then
       Write(*,*)
       Write(*,*) ' LST transformation...'
       !
       ! steps
       h=0.01_pr; hhb=b0*h; pihhb=4.0_pr*pi*hhb
       !
       ! correct density asymptotic
       !
       !================================
       ! Test neutron/proton asymptotic
       !================================
       Do it=itmin,itmax 
          ! neutron/proton density decay constant
          itass=it; decay=ass(itass); decay2=decay**2; 
          rmsit=rms(itass)+one; irmsit=rmsit/hhb
          bb=Dble((itass-1)*npr(2))*Sqrt(1.440_pr)/hb0
          ! correct density Rend
          Rend=40.0_pr; irmax=Rend/hhb
          ! DeAllocate/Allocate
          If(Allocated(dsx)) Deallocate(dsx,dsy,dsyT,dsyi,dsyii,dsy1,dsy1i  &
               ,spb0,spc0,spd0,spbi,spci,spdi)
          Allocate(dsx(irmax),dsy(irmax),dsyT(irmax),dsyi(irmax),dsyii(irmax),  &
               dsy1(irmax),dsy1i(irmax),spb0(irmax),spc0(irmax),spd0(irmax),  &
               spbi(irmax),spci(irmax),spdi(irmax))
          ! HFB+HO_{L=0} density 'dsy' and its normalization
          ! integral 'dsyi' at points 'dsx' with step 'hhb=h*b0'
          ! up to the point where 'dsy*dsx*dsx < epsdsy'
          msw=0; snorm=zero     
          Do While(Abs(snorm-Dble(npr(itass))).Gt.0.01.And.msw.Lt.25)
             msw=msw+6 ! increase for good norm of HFB+HO_{L=0}
             itass=-itass; s1=zero
             Do  ir=1,irmax
                rmmho=hhb*Dble(ir); immho=ir
                Call densitr(itass,rmmho,sN,sP,msw)
                If(itass.Eq.1) Then
                   s=sN; sT=sP
                Else
                   s=sP; sT=sN
                End If
                s=s*Dnfactor(itass)
                z1=s*rmmho**2; s1=s1+z1
                dsyT(ir)=sT; 
                dsy(ir)=s; dsyi(ir)=hhb*s1 !p-ho density and its integral
                immho=ir                   !up to the point immho
                If(z1.Lt.epsdsy) Exit
             End Do
             snorm=pihhb*s1  ! HFB+HO_{L=0} norm
          End Do
          Call deri(hhb,immho,dsy,spb0)
          ! MIN: Find 'rmm1', the first minimum of Ln(HFB+HO_{L=0})'
          z10=1.0d10
          Do ir=irmsit,immho-5
             denm1(it)=dsy(ir); denm2(it)=dsyT(ir)
             z1=spb0(ir)/dsy(ir)
             If(z1.Le.z10) Then
                imm1=ir; rmm1=hhb*Dble(ir);  z1mm1=z1; Else; Exit
             End If
             z10=z1
          End Do
          rdenm(itass)=rmm1/b0
          ! no minimum of Ln(HFB+HO_{L=0})'
          If(rmm1.Ge.hhb*Dble(immho-5)) Then
             Write(*,*)
             Write(*,*) '#####################################'
             Write(*,*) 'Please increase Nsh NB!!!(NO THO RUN)'
             Write(*,*) '#####################################'
             Write(*,*)
             iasswrong(itass)=-1
	         ierror_flag=ierror_flag+1
	         ierror_info(ierror_flag)='STOP: in thofun: Please increase Nsh NB!!!(NO THO RUN)'
			 Return  
             If(lpr) Then
                !Open(1110,file='dat0.dat')
                !Write(1110,*) ' r rhoh Log(rhoh)'
                Do ir=5,immho-5
                   !Write(1110,'(14(4x,e13.6))') hhb*Dble(ir),dsy(ir),spb0(ir)/dsy(ir)
                End Do
                !Close(1110)
             End If
             Return
          End If
       End Do
       !
       ! Asymptotics
       If( (denm1(1)-denm2(1))*(denm2(2)-denm1(2)).Le.zero ) Then 
          itass=1; If(ass(1).Gt.ass(2)) itass=2                 ! mismatch: use old asymptotic (lower decay) 
       Else        
          itass=2; If(denm1(1).Gt.denm2(1)) itass=1             ! use new asymptotic (higher density)
       End If
       !
       iasswrong(3)=0
       If(iasswrong(itass).Ne.0) iasswrong(3)=iasswrong(itass)  ! wrong assymptotic=> reinforce HO results 
       !
       Write(*,*) '                          min.point         neutron density       proton density'
       Write(*,*) '  1. Neutron min.point ',rdenm(1),denm1(1),denm2(1)
       Write(*,*) '  2. Protons min.point ',rdenm(2),denm2(2),denm1(2)
       Write(*,*) '     Neutron/Proton decay',ass(1),'/',ass(2)
       Write(*,*) '     Chosen Case=',itass
       !
       !===================
       ! Actual asymptotic
       !===================
       ! neutron/proton density decay constant
       decay=ass(itass); decay2=decay**2; 
       rmsit=rms(itass)+one; irmsit=rmsit/hhb
       bb=Dble((itass-1)*npr(2))*Sqrt(1.440_pr)/hb0
       ! correct density Rend
       Rend=40.0_pr; irmax=Rend/hhb
       ! DeAllocate/Allocate
       If(Allocated(dsx)) Deallocate(dsx,dsy,dsyT,dsyi,dsyii,dsy1,dsy1i  &
            ,spb0,spc0,spd0,spbi,spci,spdi)
       Allocate(dsx(irmax),dsy(irmax),dsyi(irmax),dsyii(irmax),  &
            dsy1(irmax),dsy1i(irmax),spb0(irmax),spc0(irmax),spd0(irmax),  &
            spbi(irmax),spci(irmax),spdi(irmax))
       ! HFB+HO_{L=0} density 'dsy' and its normalization
       ! integral 'dsyi' at points 'dsx' with step 'hhb=h*b0'
       ! up to the point where 'dsy*dsx*dsx < epsdsy'
       msw=0; snorm=zero     
       Do While(Abs(snorm-Dble(npr(itass))).Gt.0.01.And.msw.Lt.25)
          msw=msw+6 ! increase for good norm of HFB+HO_{L=0}
          itass=-itass; s1=zero
          Do  ir=1,irmax
             rmmho=hhb*Dble(ir); immho=ir
             Call densitr(itass,rmmho,sN,sP,msw)
             If(itass.Eq.1) Then
                s=sN; sT=sP
             Else
                s=sP; sT=sN
             End If
             s=s*Dnfactor(itass)
             z1=s*rmmho**2; s1=s1+z1
             dsy(ir)=s; dsyi(ir)=hhb*s1 !p-ho density and its integral
             immho=ir                   !up to the point immho
             If(z1.Lt.epsdsy) Exit
          End Do
          snorm=pihhb*s1  ! HFB+HO_{L=0} norm
       End Do
       Call deri(hhb,immho,dsy,spb0)
       ! MIN: Find 'rmm1', the first minimun of Ln(HFB+HO_{L=0})'
       z10=1.0d10
       Do ir=irmsit,immho-5
          z1=spb0(ir)/dsy(ir)
          If(z1.Le.z10) Then
             imm1=ir; rmm1=hhb*Dble(ir);  z1mm1=z1; Else; Exit
          End If
          z10=z1
       End Do
       rmm1b0=rmm1/b0
       ! MAX: Find 'rmmmax', the first maximum of Ln(HFB+HO_{L=0})'
       z10=z1mm1
       Do ir=imm1,immho-5
          z1=spb0(ir)/dsy(ir)
          If(z1.Ge.z10) Then
             immmax=ir; rmmmax=hhb*Dble(ir);  z1mmmax=z1; Else; Exit
          End If
          z10=z1
       End Do
       rmmmaxb0=rmmmax/b0
       !
       ! END: Find 'rmm2', the last point of Ln(HFB+HO_{L=0}) at the level of the first minimum
       z10=z1mm1
       Do ir=immmax,immho-5
          z1=spb0(ir)/dsy(ir)
          If(z1.Ge.z10) Then
             imm2=ir; rmm2=hhb*Dble(ir);  z1mm2=z1
          End If
       End Do
       rmm2b0=rmm2/b0
       !
       ! MID: Find 'rmmm', the MinMaX mid point
       z10=half*(z1mmmax+z1mm1)
       Do ir=imm1,immho
          z1=spb0(ir)/dsy(ir)
          If(z1.Ge.z10) Exit
          immm=ir; rmmm=hhb*Dble(ir);  z1mmm=z1
       End Do
       rmmmb0=rmmm/b0
       ! -------------------------------------------------------------
       ! Important points required:
       ! Minimum point  'rmm1'   and its log.density 'z1mm1'
       ! First maximum  'rmmmax' and its log.density 'z1mmmax'
       ! Last acceptable point 'rmm2'   and its log.density 'z1mm2'
       ! Mid point      'rmmm'   and its log.density 'z1mmm'
       ! -------------------------------------------------------------
       !
       ! fit 'aa' from the mid match point 'rmmm','z1mmm'
       If(z1mmm.Ge.-decay) z1mmm=(-decay+z1mm1)/two   !just in case
       ! the 3rd order equation
       sqsq=rmmm*(two*bb+decay2*rmmm)
       bbb1=one+rmmm*z1mmm
       yyy1=(2.0_pr*bbb1**6+18.0_pr*bbb1**3*sqsq+27.0_pr*sqsq**2+ &
            3.0_pr*Sqrt(3.0_pr)*sqsq**1.5_pr*Sqrt(4.0_pr*bbb1**3+ &
            27.0_pr*sqsq))**(1.0_pr/3.0_pr)
       !NB it was DReal but the lincker complains
       aa=Real((2.*bbb1**2*yyy1+2.**(2.0_pr/3.0_pr)*yyy1**2-&
            2.*2.**(1.0_pr/3.0_pr)*(-bbb1**4-&
            6.0_pr*bbb1*sqsq))/(6.*yyy1))
       aa=(-4.0_pr*bb*rmmm-decay2*rmmm**2+aa)*0.250_pr
       !write(*,*)  ' aa=',aa !If(aa.Le.zero) aa=zero  ! in this case take l=0
       !
       ! matching logder at Rmin='rmm1' and 'rmmx'
       ! log density for rmm1<r<rmmx is z1mm1+aaa*(r-rmm1)**2/r**assm
       ! rmmx=rmm2     !rmm2 is last acceptable point
       rmmx=rmmmax      !rmmmax is the firts maximum point
       sqsq=Sqrt(decay2+(4.0_pr*(aa+bb*rmmx))/rmmx**2)
       z1mmx=-one/rmmx-sqsq-(two*bb+decay2*rmmx)/  &
            (4.0_pr*aa+rmmx*(4.0_pr*bb+decay2*rmmx))
       z1mmxx=(two*(8.0_pr*aa**2+16.0_pr*aa*bb*rmmx+ &
            12.0_pr*bb**2*rmmx**2+two*aa*decay2*rmmx**2+ &
            6.0_pr*bb*decay2*rmmx**3+decay2**2*rmmx**4+ &
            rmmx*(two*aa+bb*rmmx)*sqsq*(4.0_pr*aa+ &
            rmmx*(4.0_pr*bb+decay2*rmmx))))/(rmmx**2*(4.0_pr*aa+ &
            rmmx*(4.0_pr*bb+decay2*rmmx))**2)
       assm=rmmx*(two/(rmmx-rmm1)-z1mmxx/(z1mmx-z1mm1))
       asm1=one-assm; asm2=two-assm; asm3=three-assm
       aaa=-rmmx*rmmx**assm*z1mmxx/  &
            ((rmm1-rmmx)*(assm*(rmm1-rmmx)+two*rmmx))
       alaex=Log(dsy(imm1))-&
            (z1mm1*rmm1+two*aaa*rmm1**asm3/(asm3*asm2*asm1))
       ! correct density 'dsy1'
       Do ir=1,irmax
          qq=hhb*Dble(ir); dsx(ir)=qq
          If(ir.Le.imm1) Then
             ! inner region (r < rmm1)
             dsy1(ir)=dsy(ir)
          Else
             sqsq=Sqrt(decay2+4.0_pr*(aa/qq**2+bb/qq))
             yyy=-qq*sqsq-Log(qq*qq*sqsq)-&
                  (two*bb/decay)*Log((two*bb+decay2*qq+decay*qq*sqsq)/decay)
             ! take complex in the case of negative aa
             aac=aa
             yyy=yyy+two*Sqrt(aac)*  &    
                  Log((two*aac+bb*qq+Sqrt(aac)*qq*sqsq)/(two*aac**1.50_pr*qq))
             If(qq.Le.rmmx) Then
                ! region (rmm1 < r < rmmx)
                aldsy1=alaex+z1mm1*qq+aaa*qq**asm1*  &
                     (asm3*asm2*rmm1**2-two*asm3*asm1*rmm1*qq+asm2*asm1*qq*qq)/  &
                     (asm3*asm2*asm1)
                dsy1(ir)=Exp(aldsy1)
                bex=aldsy1-yyy
             Else
                ! region (r > rmmx)
                dsy1(ir)=Exp(bex+yyy)
             End If
          End If
       End Do
       ! correct density norm
       snorm1=pihhb*Sum(dsy1*dsx*dsx)
       ! normalized correct density 'dsy1'
       dsy1=snorm*dsy1/snorm1
       ! zero constant
       c00=(dsy1(1)/dsy(1))**(1.0_pr/3.0_pr)
       ! splining correct density and its integral
       s1=zero
       Do ir=1,irmax
          s1=s1+ dsy1(ir)*dsx(ir)**2
          dsy1i(ir)=hhb*s1
       End Do
       !
       ! correct density dsy1 and its integral dsy1i known up to irmax
       Call csplin(irmax,dsx,dsy1 ,spb0,spc0,spd0)
       Call csplin(irmax,dsx,dsy1i,spbi,spci,spdi)
       !
       ! print 'dat1.dat' with HFB+HO and 'correct' densities
       ! and their Log derivatives at lpr=.true.
       If(lpr) Then
          !Open(1110,file='dat0.dat')
          !Write(1110,*) ' r rhoh rhoc Log(rhoh)'' Log(rhoc)'' '
          Do ir=5,immho-5
             ! ho density derivative
             s=(45.0_pr*( dsy(ir+1)-dsy(ir-1))-9.0_pr*  &
                  (dsy(ir+2)-dsy(ir-2))+dsy(ir+3)-dsy(ir-3))/(60.0_pr*hhb)/dsy(ir)
             !correct density derivative
             s1=(45.0_pr*(dsy1(ir+1)-dsy1(ir-1))-9.0_pr*(dsy1(ir+2)-&
                  dsy1(ir-2))+dsy1(ir+3)-dsy1(ir-3))/(60.0_pr*hhb)/dsy1(ir)
             !Write(1110,'(14(4x,e13.6))') dsx(ir),dsy(ir),dsy1(ir),s,s1
          End Do
          !Close(1110)
          !Open(1111,file='dat1.dat')
          !Write(1111,*) ' Dimensionless_qq  Invers_f Invers_f1 Invers_f2 Invers_f3 '
       End If
       !
       !=======================================
       ! Calculations at given dimensionless 'qq'
       !=======================================
       ! f(R->0)=c00*R, therefore Invers_f(R)=R/c00
       fj=h/c00; bmm3=zero; z1=zero; ir=0
       Do iqq=1,immho
          qq=Dble(iqq)*h
          ! HFB+HO density and integral at 'b0*qq'
          zqq=dsy(iqq); zqqi=dsyi(iqq)
          ! Iterations to find 'fj=Invers_f(qq)'
          ! NB! f[Invers_f(qq)]=qq
          iter=0; df=0.00010_pr
          Do While(Abs(df).Ge.epsf.And.iter.Le.500)
             iter=iter+1; fjb=fj*b0
             Call cseval(irmax,fjb,dsx,dsy1 ,spb0,spc0,spd0,zfj1 )
             Call cseval(irmax,fjb,dsx,dsy1i,spbi,spci,spdi,zfj1i)
             qqup=(Log(zfj1i/zqqi))*zfj1i; qqdn=zfj1*b0*fjb**2
             ! Secant  & Newton
             If(zfj1i.Le.zqqi.And.df.Le.0.0_pr) df=-half*df
             If(zfj1i.Gt.zqqi.And.df.Gt.0.0_pr) df=-half*df
             If(Abs(qqdn).Gt.Abs(qqup).And.iter.Le.20) df=-qqup/qqdn
             fj=fj+df
          End Do
          fj1=(zqq*qq*qq)/(zfj1*fj*fj)
          fj2=(fj1-z1)/h; z1=fj1
          If(qq.Gt.rmm2b0) Then
             If(fj1.Ge.bmm3) Then
                bmm3=fj1; Else; Exit
             End If
          End If
          dsy(iqq)=fj
          dsyi(iqq)=fj1
          dsyii(iqq)=fj2
          iqqmax=iqq
       End Do
       imm3=iqqmax-50
       rmm3=Dble(imm3)*h
       amm3=dsy(imm3)
       bmm3=dsyi(imm3)
       cmm3=dsyii(imm3)
       ! second and p13 derivatives up to 'iqqmax''
       Call Deri(h,iqqmax,dsyi,spb0)
       Call Deri(h,iqqmax,spb0,spc0)
       !
       If(Allocated(fdsx)) Deallocate(fdsx,fdsy,fdsy1,fdsy2,fdsy3,  &
            fspb0,fspc0,fspd0,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,  &
            fspb3,fspc3,fspd3)
       Allocate(fdsx(iqqmax),fdsy(iqqmax),fdsy1(iqqmax),fdsy2(iqqmax),  &
            fdsy3(iqqmax),fspb0(iqqmax),fspc0(iqqmax),fspd0(iqqmax),  &
            fspb1(iqqmax),fspc1(iqqmax),fspd1(iqqmax),fspb2(iqqmax),  &
            fspc2(iqqmax),fspd2(iqqmax),fspb3(iqqmax),fspc3(iqqmax),  &
            fspd3(iqqmax))
       Do iqq=1,iqqmax
          fdsx(iqq)=Dble(iqq)*h
          fdsy(iqq)=dsy(iqq)
          fdsy1(iqq)=dsyi(iqq)
          fdsy2(iqq)=spb0(iqq)
          fdsy3(iqq)=spc0(iqq)
          !
          ! print 'dat1.dat' with fj=Inverse_f(qq) and its derivatives
          ! fj1..3 at no smoothing when lpr=.true.
          If(lpr) Then
             !Write(1111,'(14(4x,e13.6))') fdsx(iqq),fdsy(iqq),  &
             !     fdsy1(iqq),fdsy2(iqq),fdsy3(iqq)
          End If
       End Do
       !
       If(Allocated(dsx)) Deallocate(dsx,dsy,dsyi,dsyii,dsy1,dsy1i,  &
            spb0,spc0,spd0,spbi,spci,spdi)
       !
       Call csplin(iqqmax,fdsx,fdsy ,fspb0,fspc0,fspd0)
       Call csplin(iqqmax,fdsx,fdsy1,fspb1,fspc1,fspd1)
       Call csplin(iqqmax,fdsx,fdsy2,fspb2,fspc2,fspd2)
       Call csplin(iqqmax,fdsx,fdsy3,fspb3,fspc3,fspd3)
       !
       Do ir=lout,lfile
          Write(ir,*)
          Write(ir,*) ' Legandre points= ',msw
          Write(ir,*) ' b0, decay=       ',b0,decay
          Write(ir,*) ' h, hhb=          ',h,hhb
          Write(ir,*) ' rms, rmsit=      ',rms(itass),rmsit
          Write(ir,*) ' Rend,   irmax=   ',Rend,irmax
          Write(ir,*) ' HORend, immho=   ',rmmho,immho
          Write(ir,*) ' snorm,snorm1=    ',snorm,snorm1
          Write(ir,*) ' snorm/snorm1=    ',snorm/snorm1
          Write(ir,*) ' min:rmm1,/b0=    ',rmm1,rmm1b0
          Write(ir,*) ' max:rmmmax,/b0=  ',rmmmax,rmmmaxb0
          Write(ir,*) ' last:rmm2,/b0=   ',rmm2,rmm2b0
          Write(ir,*) ' num:rmm3*b0,rmm3=',rmm3*b0,rmm3
          Write(ir,*) ' rmmho, rmmho/b0= ',rmmho,rmmho/b0
          Write(ir,*) ' alaex,bex=       ',alaex,bex
          Write(ir,*) ' aa,bb=           ',aa,bb
          Write(ir,*) ' L_eff=           ',(Sqrt(one+4.0_pr*aac)-one)/two
          Write(ir,*) ' amm3, bmm3=      ',amm3,bmm3
          Write(ir,*) ' cmm3, one/c00=   ',cmm3,one/c00
          Write(ir,*) ' z1mm1,z1mmm=     ',z1mm1,z1mmm
          Write(ir,*)
       End Do
       ! print 'dat2..4.dat' after when lpr=.true.
       ! dat2.dat: 'r=b0*qq',correct density
       ! dat3.dat   qq, Invers_f,Invers_f1...3
       ! dat4.dat   qq, f,f1,f2,f3; 'qq' are the Gauss points
       If(lpr) Then
          !Close(1111)
          !Open(1112,file='dat2.dat')
          !Open(1113,file='dat3.dat')
          !Open(1114,file='dat4.dat')
          !Write(1112,*) ' r  den_correct'
          !Write(1113,*) ' qq Invers_f Invers_f1 Invers_f2 Invers_f3'
          !Write(1114,*) ' qq f f1 f2 f3'
          Do iqq=1,ngh
             Do ir=1,ngl
                qq=Sqrt(xh(iqq)**2+xl(ir))
                Call densitr(itass,b0*qq,sN,sP,msw)
                If(itass.Eq.1) Then
                   s=sN; sT=sP
                Else
                   s=sP; sT=sN
                End If
                s=s*Dnfactor(itass)
                If(qq.Le.rmm3) Then
                   Call cseval(iqqmax,qq,fdsx,fdsy ,fspb0,fspc0,fspd0,fj)
                   Call cseval(iqqmax,qq,fdsx,fdsy1,fspb1,fspc1,fspd1,fj1)
                   Call cseval(iqqmax,qq,fdsx,fdsy2,fspb2,fspc2,fspd2,fj2)
                   Call cseval(iqqmax,qq,fdsx,fdsy3,fspb3,fspc3,fspd3,fj3)
                Else
                   fj=amm3+bmm3*(qq-rmm3)+cmm3*(qq-rmm3)**2/two
                   fj1=bmm3+cmm3*(qq-rmm3)
                   fj2=cmm3; fj3=zero
                End If
                f=qq; f1=one/fj1; f2=-fj2*f1**3
                f3=three*fj2**2*f1**5-fj3*f1**4
                s=s*qq*qq/(fj*fj*fj1)
                !Write(1112,'(14(4x,e13.6))') b0*fj,s
                !Write(1113,'(14(4x,e13.6))') qq,fj,fj1,fj2,fj3
                !Write(1114,'(14(4x,e13.6))') qq,f,f1,f2,f3
             End Do
          End Do
          !Close(1112); Close(1113); Close(1114)
       End If
       !=========================
    Else  !KEY=1 CALCULATIONS
       !=========================
       !
       ! Calculations of Invers_f(qq),_f'(qq),_f''(qq),_f'''(qq)
       If(qq.Le.rmm3) Then
          Call cseval(iqqmax,qq,fdsx,fdsy ,fspb0,fspc0,fspd0,fj)
          Call cseval(iqqmax,qq,fdsx,fdsy1,fspb1,fspc1,fspd1,fj1)
          Call cseval(iqqmax,qq,fdsx,fdsy2,fspb2,fspc2,fspd2,fj2)
          Call cseval(iqqmax,qq,fdsx,fdsy3,fspb3,fspc3,fspd3,fj3)
       Else
          fj=amm3+bmm3*(qq-rmm3)+cmm3*(qq-rmm3)**2/two
          fj1=bmm3+cmm3*(qq-rmm3); fj2=cmm3; fj3=zero
       End If
       ! Calculations of f(fj),f'(fj),f''(fj),f'''(fj)
       f=qq; f1=one/fj1; f2=-fj2*f1**3
       f3=three*fj2**2*f1**5-fj3*f1**4
    End If
    Return
  End Subroutine thofun
  !=======================================================================
  !
  !=======================================================================
  Subroutine printLST(ffn,ffp)
    !---------------------------------------------------------------------
    ! prints ffn, ffp
    !---------------------------------------------------------------------
    Implicit None
    Integer(ipr), Save :: ihli,Ifle
    Real(pr), Save :: sn,sp,rrf,rfn,rfp
    Real(pr) :: ffn(nghl),ffp(nghl)
    !
    Ifle=76+iLST1
    If(Ifle.Eq.76) Open(Ifle,file='hodennp.dat',status='unknown')
    If(Ifle.Eq.77) Open(Ifle,file='thodenp.dat',status='unknown')
    sn=zero; sp=zero
    Write(Ifle,*) 'r=Sqrt[fh**2+fl**2)    den_N    den_P    fh    fl '
    Do ihli=1,nghl
       rrf=Sqrt(fh(ihli)**2+fl(ihli)**2)
       rfn=ffn(ihli);            rfp=ffp(ihli)
       sn=sn+rfn*wdcor(ihli); sp=sp+rfp*wdcor(ihli)
       Write(Ifle,'(12(1x,e16.8))') rrf,rfn,rfp,fh(ihli),fl(ihli),wdcor(ihli)
    End Do
    Write(Ifle,*) ' neutron and proton norms:',sn,sp
    Close(Ifle)
    Return
  End Subroutine printLST
  !=======================================================================
  !
  !=======================================================================
  Subroutine densitr(it,xr,yr,yrP,msw)
    !---------------------------------------------------------------------
    ! Calculates Legandre decomposition of neutron(proton) 'it=1(2)'
    ! HFB+HO_{L=0}(r) density 'yr' at point 'xr' (in fm)
    !---------------------------------------------------------------------
    Implicit None
    Integer(ipr) :: it,msw
    Integer(ipr), Save :: iw,ik,il,i0,i02,jk,nsa,nsb,nrb,ny,nyy,ib,nd,  &
         n1,n2,n1n2nd,ibit,ibitnb,nzb,mlb,ngh1,ngl1
    ! msw=20 test for protons at the crasy case of U 212 120 92
    ! msw=3  s,a1=107.665062 80.82396  s/s1=1.33209
    ! msw=6  s,s1=88.973061 80.99903  s/s1=1.09844
    ! msw=12 s,s1=92.025592 80.99595  s/s1=1.13617
    ! msw=18 s,s1=92.000017 80.99595  s/s1=1.13585
    ! msw=24 s,s1=92.000010 80.99595  s/s1=1.13585
    ! Taken up to msw=24
    Real(pr), Allocatable, Save :: xmw(:),yi(:,:)
    Real(pr) :: phy(msw,nzrlx),anl(msw,msw),yl(msw,msw)
    Real(pr), Save :: sl,w,hw,ct2,s,frit,fritP,wdcorin,bzi,bri,ct,st,z,t
    Real(pr) :: xr,yr,yrp
    !
    ngh1=ngh+1; ngl1=ngl+1
    If(it.Lt.0) Then
       it=-it
       If(Allocated(xmw)) Deallocate(xmw,yi)
       Allocate(xmw(msw),yi(msw,msw))
       wdcorin=one/Sqrt(pi*bz*bpp); bzi=one/bz; bri=one/bp
       ! 'msw' mesh-points in 'angle' space
       xmw(1)=zero; hw=half*pi/Dble(msw-1)
       Do  il=2,msw
          xmw(il)=hw*Dble(il-1)
       End Do
       ! coefficients for the L-decomposition
       sl=4.0_pr
       Do  il=1,msw
          sl=0.250_pr*sl
          Do ny=1,il
             anl(ny,il)=iv(il-ny)*fak(2*(il+ny-2))*  &
                  fi(il-ny)*fi(il+ny-2)*fi(2*ny-2)*sl
          End Do
       End Do
       Do iw=1,msw
          w=xmw(iw); ct2=Cos(w)**2
          Do il=1,msw
             yi(iw,il)=zero; s=zero
             Do nyy=1,il
                ny=il+1-nyy; s=s*ct2+anl(ny,il)
             End Do
             yl(iw,il)=s*sq(4*il-3)
          End Do
          yi(iw,iw)=one
       End Do
       Call lingd(msw,msw,msw,msw,yl,yi,s,il)
    End If
    ! 'xr/yr' calculations
    Do iw=1,msw
       w=xmw(iw); ct=Cos(w); st=Sin(w)
       z=ct*xr*bzi; t=(st*xr*bri)**2
       Call gaupolr(z,t)
       ik=0
       Do ib=1,nb
          nd=id(ib); i0=ia(ib)
          Do n2=1,nd
             ik=ik+1; i02=i0+n2
             nzb=nz(i02);  nrb=nr(i02); mlb=nl(i02)
             phy(iw,ik)=qh(nzb,ngh1)*ql(nrb,mlb,ngl1)*wdcorin
          End Do
       End Do
    End Do
    ! 'yr' over the blocks
    yr=zero; yrP=zero; ik=0
    Do ib=1,nb
       nd=id(ib); i0=ia(ib)
       Do n2=1,nd
          jk=ik; ik=ik+1; i02=i0+n2; nsb=ns(i02)
          Do n1=n2,nd
             jk=jk+1; i02=i0+n1; nsa=ns(i02)
             If(nsa.Eq.nsb) Then
                ibit=ib; ibitnb=ib+nbx; n1n2nd=n1+(n2-1)*nd
                frit=rk(n1n2nd,ibit); fritP=rk(n1n2nd,ibitnb)
                If(n1.Ne.n2) Then
                   frit=two*frit; fritP=two*fritP
                End If
                s=zero
                Do iw=1,msw
                   s=s+yi(1,iw)*phy(iw,ik)*phy(iw,jk)
                End Do
                yr=yr+frit*s; yrP=yrP+fritP*s
             End If
          End Do !n2
       End Do !n1
    End Do !ib
    !
    Return
  End Subroutine densitr
  !=======================================================================
  !
  !=======================================================================
  Subroutine gaupolr(z,x)
    !---------------------------------------------------------------------
    ! see 'gaupol'
    !---------------------------------------------------------------------
    Implicit None
    Real(pr) :: z,x
    Real(pr) :: w0,w00,w4pii,dsq,d1,d2
    Integer(ipr) :: N,L,NGH1,NGL1
    !
    NGH1=NGH+1; NGL1=NGL+1
    W4PII=PI**(-0.250_pr); W0=W4PII*Exp(-HALF*Z*Z)
    ! W0 =W0*SQRT(Z) NOT MULTIPLIED BY WDCOR
    QH(0,NGH1)=W0; QH(1,NGH1)=SQ(2)*W0*Z
    Do N=2,NZM
       QH(N,NGH1)=SQI(N)*(SQ(2)*Z*QH(N-1,NGH1)-SQ(N-1)*QH(N-2,NGH1))
    End Do
    W00=SQ(2)*Exp(-HALF*X)
    Do L=0,NLM
       If(L.Eq.0) Then
          W0=W00*Sqrt(HALF)
       Else
          W0=W00*Sqrt(HALF*X**L)
       End If
       QL(0,L,NGL1)=WFI(L)*W0; QL(1,L,NGL1)=(Dble(L+1)-X)*WFI(L+1)*W0
       Do N=2,NRM
          DSQ=SQ(N)*SQ(N+L); D1=Dble(2*N+L-1)-X
          D2=SQ(N-1)*SQ(N-1+L)
          QL(N,L,NGL1)=(D1*QL(N-1,L,NGL1)-D2*QL(N-2,L,NGL1))/DSQ
       End Do
    End Do
    Return
  End Subroutine gaupolr
  !=======================================================================
  !
  !=======================================================================
  Subroutine lingd(ma,mx,n,m,a,x,d,Ifl)
    !---------------------------------------------------------------------
    ! Solves a system of linear equations A*X=B
    ! At the beginning the matrix B is stored in X
    ! During the calculation it will be overwritten
    ! D is the determinant of A
    !---------------------------------------------------------------------
    Implicit None
    Integer(ipr) :: ma,mx,n,m,Ifl
    Integer(ipr), Save :: i,j,k,l,k1,n1
    Real(pr) ::  a(ma,m),x(mx,m),d
    Real(pr), Save :: tollim,one,zero,p,q,tol,cp,cq
    Data tollim/1.d-10/,one/1._pr/,zero/0._pr/
    Ifl=1; p=zero
    Do i=1,n
       q=zero
       Do j=1,n
          q=q+Abs(a(i,j))
       End Do
       If(q.Gt.p) p=q
    End Do
    tol=tollim*p; d=one
    Do k=1,n
       p=zero
       Do j=k,n
          q=Abs(a(j,k))
          If(q.Lt.p) Cycle
          p=q; i=j
       End Do
       If (p.Lt.tol) Then
          Write (6,200) ('-',j=1,80),tol,i,k,a(i,k),('-',j=1,80)
200       Format (/1x,80a1/' *****  ERROR IN LINGD , TOLERANZ=',e10.4,  &
               ' VALUE OF A(',i3,',',i3,') IS ',e10.4/1x,80a1)
          Ifl=-1
          Return
       End If
       cp=one/a(i,k)
       If(i.Ne.k) Then
          d=-d
          Do l=1,m
             cq=x(i,l); x(i,l)=x(k,l); x(k,l)=cq
          End Do
          Do l=k,n
             cq=a(i,l); a(i,l)=a(k,l); a(k,l)=cq
          End Do
       End If
       d=d*a(k,k)
       If(k.Eq.n) Exit
       k1=k+1
       Do i=k1,n
          cq=a(i,k)*cp
          Do l=1,m
             x(i,l)=x(i,l)-cq*x(k,l)
          End Do
          Do l=k1,n
             a(i,l)=a(i,l)-cq*a(k,l)
          End Do
       End Do
    End Do
    Do l=1,m
       x(n,l)=x(n,l)*cp
    End Do
    If(n.Eq.1) Return
    n1=n-1
    Do k=1,n1
       cp=one/a(n-k,n-k)
       Do l=1,m
          cq=x(n-k,l)
          Do i=1,k
             cq=cq-a(n-k,n+1-i)*x(n+1-i,l)
          End Do
          x(n-k,l)=cq*cp
       End Do
    End Do
    Return
  End Subroutine lingd
  !=======================================================================
  !
  !=======================================================================
  Subroutine csplin (n, x, y, b, c, d)
    !---------------------------------------------------------------------
    ! file: csplin.for  (from slac)
    ! the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
    ! for a cubic interpolating spline
    ! s(x)=y(i)+b(i)*(x-x(i))+c(i)*(x-x(i))**2+d(i)*(x-x(i))**3
    ! for  x(i) .le. x .le. x(i+1)
    ! input..
    ! n=the number of data points or knots (n.ge.2)
    ! x=the abscissas of the knots in strictly increasing order
    ! y=the ordinates of the knots
    ! output..
    ! b, c, d=arrays of spline coefficients as defined above.
    ! using  p  to denote dIfferentiation,
    ! y(i)=s(x(i))
    ! b(i)=sp(x(i))
    ! c(i)=spp(x(i))/2
    ! d(i)=sppp(x(i))/6  (derivative from the right)
    ! the accompanying function subprogram  cseval  can be used
    ! to evaluate the spline, its derivative or even its 2nd derivative.
    !=====================================================================
    Implicit None
    Integer(ipr), Save :: nm1,i,ib
    Integer(ipr) :: n
    Real(pr)   :: x(n), y(n), b(n), c(n), d(n)
    Real(pr), Save    :: t,zero=0.0_pr,two=2.0_pr,tr=3.0_pr
    ! check input for consistency
    If(n.Lt.2) Stop '-n < 2 in csplin call--'
    nm1=n-1
    Do i=1, nm1
       If(x(i).Ge.x(i+1)) Stop 'x not strictly ascending in csplin call'
    End Do
    If (n.Ne.2) Then
       ! set up tridiagonal system
       ! b=diagonal, d=offdiagonal, c=right hand side.
       d(1)=x(2)-x(1); c(2)=(y(2)-y(1))/d(1)
       Do i=2, nm1
          d(i)=x(i+1)-x(i); b(i)=two*(d(i-1)+d(i))
          c(i+1)=(y(i+1)-y(i))/d(i); c(i)=c(i+1)-c(i)
       End Do
       ! end conditions.  p13 derivatives at  x(1)  and  x(n)
       ! obtained from divided dIfferences
       b(1)=-d(1); b(n)=-d(n-1); c(1)=zero; c(n)=zero
       If (n.Ne.3) Then
          c(1)=c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
          c(n)=c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
          c(1)=c(1)*d(1)**2/(x(4)-x(1))
          c(n)=-c(n)*d(n-1)**2/(x(n)-x(n-3))
          ! forward elimination
       Else
          Do i=2, n
             t=d(i-1)/b(i-1); b(i)=b(i)-t*d(i-1); c(i)=c(i)-t*c(i-1)
          End Do
       End If
       ! back substitution
       c(n)=c(n)/b(n)
       Do ib=1, nm1
          i=n-ib
          c(i)=(c(i)-d(i)*c(i+1))/b(i)
       End Do
       ! compute polynomial coefficients
       b(n)=(y(n)-y(nm1))/d(nm1)+d(nm1)*(c(nm1)+two*c(n))
       Do i=1, nm1
          b(i)=(y(i+1)-y(i))/d(i)-d(i)*(c(i+1)+two*c(i))
          d(i)=(c(i+1)-c(i))/d(i); c(i)=tr*c(i)
       End Do
       c(n)=tr*c(n); d(n)=d(n-1)
       Return
    Else
       b(1)=(y(2)-y(1))/(x(2)-x(1)); c(1)=zero; d(1)=zero
       Return
    End If
  End Subroutine csplin
  !=======================================================================
  !
  !=======================================================================
  Subroutine cseval(n,u,x,y,b,c,d,splf0)
    !---------------------------------------------------------------------
    ! this subroutine is a copy of 'cseva' but only for the function
    !---------------------------------------------------------------------
    Implicit None
    Integer(ipr) :: n
    Integer(ipr), Save :: i=1,j,k
    Real(pr) :: x(n),y(n),b(n),c(n),d(n),u,splf0
    Real(pr), Save :: dx
    If(i.Ge.n)      i=1
    If(u.Lt.x(i))   Go To 10
    If(u.Le.x(i+1)) Go To 30
    ! binary search
10  i=1
    j=n+1
20  k=(i+j)/2
    If(u.Lt.x(k)) j=k
    If(u.Ge.x(k)) i=k
    If(j.Gt.i+1) Go To 20
    ! evaluate splf0
30  dx=u-x(i)
    splf0=y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
    Return
  End Subroutine cseval
  !=======================================================================
  !
  !=======================================================================
  Subroutine deri(h,n,f1,dunl)
    !---------------------------------------------------------------------
    ! first derivative of 'f1' if the the step is 'h'
    !---------------------------------------------------------------------
    Implicit None
    Integer(ipr) :: n
    Integer(ipr), Save :: k
    Real(pr) :: h,f1(n),dunl(n)
    Real(pr), Save :: t60,t12
    Real(pr), Save :: t8=8.0_pr,t45=45.0_pr,t9=9.0_pr
    t60=1.0_pr/(h*60.0_pr); t12=1.0_pr/(h*12.0_pr)
    !
    dunl(1)=(t8*f1(2)-f1(3)+f1(1))*t12
    dunl(2)=(t45*(f1(3)-f1(1))-t9*f1(4)+f1(5)-f1(1))*t60
    dunl(3)=(t45*(f1(4)-f1(2))-t9*(f1(5)-f1(1))+f1(6))*t60
    dunl(n)=(-t8*f1(n-1)+f1(n)+f1(n-2))*t12
    dunl(n-1)=(t45*(f1(n)-f1(n-2))+t9*f1(n-3)-f1(n)-f1(n-4))*t60
    dunl(n-2)=(t45*(f1(n-1)-f1(n-3))-t9*(f1(n)-f1(n-4))-f1(n-5))*t60
    Do k=4,n-3
       dunl(k)=(t45*(f1(k+1)-f1(k-1))-&
            t9*(f1(k+2)-f1(k-2))+f1(k+3)-f1(k-3))*t60
    End Do
    Return
  End Subroutine deri
  !=======================================================================
  !
End Module THO_MODULE
!===================================================================================================================================
!#END THO_MODULE MODULE
!===================================================================================================================================
#endif
!
!
#ifndef hide_mainprog
!===================================================================================================================================
!#START MAIN PROGRAM
!===================================================================================================================================
!=========================================================================
! PARALLEL RUN ex script uses these two lines. Do not modify them!
Program hfbthoprog; Use HFBTHO; Call Main_Program; End Program hfbthoprog
  !=========================================================================
  !
  Subroutine Main_Program
    Use HFBTHO
#ifndef hide_qrpa
    Use qrpa_HFBTHO
#endif	
    Implicit None
#ifndef hide_mpi_qrpa
    Include 'mpif.h'
#endif
    Logical      :: lpr
    Integer(ipr) :: iblocase(2),nkblocase(2,5)
    Integer(ipr) :: i,it,iodd,ios,envstat,envlen
    integer(ipr) :: qrpa_points_thiscore,qrpa_start_thiscore,idivm
    integer(ipr) :: arrayTaskID,arrayTaskCount,arrayTaskMin,arrayTaskMax
    CHARACTER    :: envRead*20
    real(pr),allocatable :: qrpa_omega_points_re(:),qrpa_omega_points_im(:)
    !
    iam_mpi = 0; num_nodes_mpi = 1 !! set default value when no MPI
    ! MPI init
#ifndef hide_qrpa
#ifndef hide_mpi_qrpa       
#ifndef hide_array_qrpa
#error "mpi and array jobs can not be compiled to the same executable"
#endif	
    Call mpi_init(icom_err_mpi)
    Call mpi_comm_rank(mpi_comm_world,iam_mpi,icom_err_mpi)
    Call mpi_comm_size(mpi_comm_world,num_nodes_mpi,icom_err_mpi)
#endif	

    !! for array task, use the same variables as with mpi
#ifndef hide_array_qrpa
    CALL get_environment_variable("SLURM_ARRAY_TASK_ID", envRead,envlen,envstat,.true.);
    if(envstat.ne.0) then; write(*,*) "SLURM_ARRAY_TASK_ID env status ",envstat; stop; end if;    read(envRead(1:envlen),*) arrayTaskID;
    CALL get_environment_variable("SLURM_ARRAY_TASK_COUNT", envRead,envlen,envstat,.true.);
    if(envstat.ne.0) then; write(*,*) "SLURM_ARRAY_TASK_COUNT env status ",envstat; stop; end if; read(envRead(1:envlen),*) arrayTaskCount;
    CALL get_environment_variable("SLURM_ARRAY_TASK_MAX", envRead,envlen,envstat,.true.);
    if(envstat.ne.0) then; write(*,*) "SLURM_ARRAY_TASK_MAX env status ",envstat; stop; end if;   read(envRead(1:envlen),*) arrayTaskMax;
    CALL get_environment_variable("SLURM_ARRAY_TASK_MIN", envRead,envlen,envstat,.true.);
    if(envstat.ne.0) then; write(*,*) "SLURM_ARRAY_TASK_MIN env status ",envstat; stop; end if;   read(envRead(1:envlen),*) arrayTaskMin;
    if(arrayTaskCount.lt.1) stop "arrayTaskCount.lt.1"
    if(arrayTaskMax-arrayTaskMin +1 .ne. arrayTaskCount) then
      write(*,*) "stopping due to condition: arrayTaskMax-arrayTaskMin +1 .ne. arrayTaskCount" 
      write(*,*) "arrayTaskMax: ",arrayTaskMax," arrayTaskMin: ",arrayTaskMin," arrayTaskCount: ",arrayTaskCount
    STOP
    end if 
    Iam_mpi = arrayTaskID-arrayTaskMin
    num_nodes_mpi = arrayTaskCount
#endif	
#endif	

    Open(901,file='tho.dat',status='old')         ! open the input file 
    Do
       !-------------------------------------------------------------------
       ! Read input and namelist data for the requested nucleus 
       !-------------------------------------------------------------------
       Read(901,*,iostat=ios) n00_INI,b0_INI,q_INI,npr_INI(1),npr_INI(2)   &
            ,skyrme_INI,kindhfb_INI,inin_INI,cdef_INI,cqad_INI &
			,nkblocase(1,1:5),nkblocase(2,1:5)
       If(ios.Ne.0)     Stop 'HFBTHO TERMINATED: ERROR in tho.dat file!'
       If(n00_INI.Eq.0) then
         ! finalize mpi          
#ifndef hide_mpi_qrpa
         Call mpi_finalize(icom_err_mpi)
#endif
         !
         Stop 'HFBTHO TERMINATED: No more nuclei!'
       end if
       Call read_UNEDF_NAMELIST(skyrme_INI)  
	   !memo: Namelist /UNEDF_NAMELIST/ FunctionalName,DMEorder,DMElda,use_INM,use_TMR_pairing,hbzero, &
	   !      Crho,Cdrho,Ctau,CrDr,CrdJ,CJ,sigma,CpV0,CpV1,e2charg, &
	   !      E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM, &
	   !      mpi,gA,fpi,c1,c3,c4,cd,ce,LambdaX, &
	   !      use_cm_cor,use_charge_density,use_DME3N_terms,use_j2terms,CExPar, &
	   !      Print_Namelist     
       Call read_HFBTHO_NAMELIST	   
	   !memo: Namelist /HFBTHO_NAMELIST/ MAX_ITER_INI,epsi_INI,Add_Pairing_INI, &
	   !      icou_INI,iLST_INI,keypj_INI,iproj_INI,npr1pj_INI,npr2pj_INI, & 
	   !      DO_FITT_INI,IDEBUG_INI,Parity_INI,Print_HFBTHO_Namelist_INI
	   !	   
       iblocase=0; blomax=0; !blomax will be cahrged from the previouse solution
	   ! come back here for the next blocked level if requested
11     Continue
	   !=========================================================================== 
	   ! GROUND STATE BLOCKING WALKER: blocking candidates are predefined 
	   ! by the parent nucleus and we block them one by one 
	   !---------------------------------------------------------------------------
       Do it=1,2  
         If(nkblocase(it,1).Ne.0.And.nkblocase(it,2).Eq.0) Then
		  If(it.Eq.1) Then
		   iblocase(1)=iblocase(1)+1
		   If(iblocase(1).Gt.blomax(1)) iblocase(1)=1
		  Else
		   If(iblocase(1).Le.1) iblocase(2)=iblocase(2)+1
		  Endif
		  nkblo_INI(it,1)=Sign(iblocase(it),nkblocase(it,1))
		  nkblo_INI(it,2)=0
          Else
		  ! case of external blocking
		  nkblo_INI(it,:)=nkblocase(it,:) 		  
          Endif				
        Enddo
	   !---------------------------------------------------------------------------
	   ! MANUAL BLOCKING: manualBlocking=1 in the module (ocasionaly used)
	   ! One types which level to be blocked referencing the parent nucleus 
	   !---------------------------------------------------------------------------
        If(manualBlocking.Ne.0) Then
        Write(*,'(a,5(1pg12.4))') 'Please print the number of the neutron level to block, num='
        Read(*,*) nkblo_INI(1,1)
        Write(*,'(a,5(1pg12.4))') 'Neutron blocked level num=',nkblo_INI(1,1)
		nkblo_INI(1,2)=0
        Write(*,'(a,5(1pg12.4))') 'Please print the number of the proton level to block, num='
        Read(*,*) nkblo_INI(2,1)
        Write(*,'(a,5(1pg12.4))') 'Proton blocked level num=',nkblo_INI(2,1)
		nkblo_INI(2,2)=0
       End If
       !
       If(Trim(skyrme_INI).Eq.'FITS') Then
 	    !--------------------------------------------------------------------
	    ! Calculations for 'FITS' functional (Modifies some values if needed)
 	    !--------------------------------------------------------------------
        !
        ! JASON 15 September
        DMEORDER=-1; DMELDA=0; use_cm_cor=.False.; use_TMR_pairing=0
		!
        HBZERO=20.735530D0; E2CHARG=1.439978400000000; 
        CRHO(0)=-731.222785829509800000d0;  CDRHO(0)= 855.690051584978500000d0;   CTAU(0)=-0.543988860905982100d0;
        CRHO(1)= 263.710305524676100000d0;  CDRHO(1)=-176.864195604041100000d0;   CTAU(1)=-33.361881866521340000d0;
		CRDR(0)=-43.29008975531540d0;       CJ(0)=0.000000000000000000d0;         CRDJ(0)=-75.260870048289410000d0
		CRDR(1)=-164.1379857135440d0;       CJ(1)=0.000000000000000000d0;         CRDJ(1)=-22.652819964871300000d0
        CPV0(0)=-186.192246596249d0;        CPV1(0)=0.5000000000000000d0;         SIGMA=0.298783982778235700d0;         
		CPV0(1)=-206.746416898386d0;        CPV1(1)=0.500000000000000d0;          CEXPAR=0.639129523762364d0
		!
		RHO_NM=0.157327162963946800d0;
        E_NM=-15.800000484870580000d0;  K_NM=225.943394883389600000d0;   SMASS_NM=0.995872580822852000d0;   
        ASS_NM=28.348338556986590000d0; LASS_NM=40.001962979089330000d0; VMASS_NM=1.248999953269958000d0; 
        !		
       Endif
	   !
       !===========================================================================
       ! run the solver in all cases EVEN/ODDS, FITS/NO-FITS	
       !===========================================================================
	   ! below are listed all variables used by the solver. They are NEVER CHANGED
	   ! if the sover runs on a slave processor they should be trensfered to it 
       ! - from tho.dat
       !      n00_INI=20;        npr1_INI=70; npr2_INI=50;    kindhfb_INI=-1;  inin_INI=-1; 
       !      b0_INI=2.234776;   q_INI=0.0;   cdef_INI=0.0;   cqad_INI=0.5;    skyrme_INI='SLY4'; 
       ! - from read_HFBTHO_NAMELIST
       !      ILST_INI=0;              keypj_INI=1;         iproj_INI=0;    npr1pj_INI=0; 
       !      icou_INI=2;              IDEBUG_INI=0;        npr2pj_INI=0;   
       !      Parity_INI=.False.;      epsi_INI=0.00001_pr; MAX_ITER_INI=101         
       !      Add_Pairing_INI=.False.; DO_FITT_INI=.False.; Print_PTHO_Namelist_INI=.True.;
       !- from read_UNEDF_NAMELIST
       !      DMEORDER=-1; DMELDA=0; use_TMR_pairing=0
       !      HBZERO=20.73553000000000;    E2CHARG=1.439978400000000; CRHO(0)=-933.3423749999999;  CRHO(1)=830.0524855000001;
       !      CDRHO(0)=861.0625000000000;  CDRHO(1)=-1064.2732500000; CTAU(0)=57.12868750000000;   CTAU(1)=24.65673650000000;
       !      CRDR(0)=-76.99620312499999;  CRDR(1)=15.65713512500000; CRDJ(0)=-92.25000000000000;  CRDJ(1)=-30.7500000000000;
       !      CJ(0)=17.20961150000000;     CJ(1)=64.57581250000000;   CPV0(0)=-258.2000000000000;  CPV0(1)=-258.2000000000000;
       !      CPV1(0)=0.5000000000000000;  CPV1(1)=0.500000000000000; SIGMA=0.1666666666666667;    CEXPAR=1.000000000000000;
       !      E_NM=-15.97214914144462;     K_NM=229.9009644826037;    SMASS_NM =1.439546988976078; RHO_NM =0.1595387567117334;
       !      ASS_NM =32.00430281505202;   LASS_NM=45.96175148046161; VMASS_NM =1.249838547196253;                
       !      MPI=0.6995945261023822;      GA=1.290000000000000;      FPI=0.4683223517486062;      C1=-0.1598130000000000; 
       !      C3 =-0.6708200000000;        C4 =0.6708200000000000;    CD =-2.062000000000000;      CE=-0.6250000000000;       
       !      LAMBDAX =3.547896604156107;        
       !      USE_INM=.false.;             USE_CM_COR =.true.;        USE_DME3N_TERMS=.true.;      USE_J2TERMS =.true.; 
       !      USE_CHARGE_DENSITY=.false.;  PRINT_NAMELIST=.true.;
       !   
       Call HFBTHO_SOLVER
       !
       If(ierror_flag.Ne.0) Then
		 Write(*,*)
		 Write(*,'(a)') ' ERRORS IN HFBTHO_SOLVER'	   	 
	     Do i=1,ierror_flag
		 Write(*,'(a,i2,2x,a)') ' error_flag=',i,ierror_info(i)
		 Enddo
		 Write(*,*)
	   Else
		 Write(*,*)
		 Write(*,'(a)') ' HFBTHO_SOLVER ended without errors'	   	 
		 Write(*,*)
	   Endif	   
	   !
	   ! at this point all HFB results are available
	   ! if the sover runs on a slave processor all results should be trensfered to the master processor 
       !===========================================================================
       !
	   ! go for new blocked level
	   inin_INI=-Abs(inin) ! Always start from old *.hel files when blocking
	   If(iblocase(1).Ne.0.And.iblocase(1).Lt.blomax(1)) go to 11 
	   If(iblocase(2).Ne.0.And.iblocase(2).Lt.blomax(2)) go to 11    
	   !
	   If(lout.Lt.lfile) Close(lfile)                 ! close the output
       !
#ifndef hide_qrpa
       !=============================================================================================
       ! FAM QRPA
       !=============================================================================================
       ! READ qrpa.inp FILE
       Open(112,file='qrpa.inp')
       Read(112,*) 
       Read(112,*) qrpa_calc_mode,qrpa_file_io
       Read(112,*) QRPA_Operator_type, T_qrpa_responce,L_qrpa_responce,K_qrpa_responce       
       Read(112,*) qrpa_eps;         
       Read(112,*) max_iter_qrpa;    
       Read(112,*) qrpa_nbroyden
       Read(112,*) qrpa_alphamix;    
       Read(112,*) qrpa_Romega,qrpa_Iomega,qrpa_Romega_step,qrpa_Iomega_step  !! line params
       Read(112,*) qrpa_circRe,qrpa_circIm,qrpa_circRad;                      !! circle params
       Read(112,*) qrpa_HcircRe,qrpa_HcircRad;                                !! half-circle params
       Read(112,*) max_qrpa_points
       Close(112)   
       !
       ! allocate table for all points to be calculated
       if(allocated(qrpa_omega_points_re)) deallocate(qrpa_omega_points_re,qrpa_omega_points_im)
       allocate(qrpa_omega_points_re(max_qrpa_points),qrpa_omega_points_im(max_qrpa_points))
       !
       ! default settings without mpi
       qrpa_start_thiscore = 1                   !! start by default from first item on the table
       qrpa_points_thiscore = max_qrpa_points    !! by default this core calculates all points
       !
       ! set all omega points in tables. all mpi threads will have the same table
       select case(qrpa_calc_mode)
       case(1) !! line
         Do i=1,max_qrpa_points 
           qrpa_omega_points_re(i) = qrpa_Romega + real(i-1,kind=pr)*qrpa_Romega_step
           qrpa_omega_points_im(i) = qrpa_Iomega + real(i-1,kind=pr)*qrpa_Iomega_step
         end do

       case(2)  !! circle
         Do i=1,max_qrpa_points 
           qrpa_omega_points_re(i) = qrpa_circRe + qrpa_circRad * cos( 2.0_pr*pi*real(i-1,kind=pr)/real(max_qrpa_points,kind=pr) )
           qrpa_omega_points_im(i) = qrpa_circIm + qrpa_circRad * sin( 2.0_pr*pi*real(i-1,kind=pr)/real(max_qrpa_points,kind=pr) )
         end do

       case(3)  !! half-circle
         Do i=1,max_qrpa_points 
           qrpa_omega_points_re(i) = qrpa_HcircRe + qrpa_HcircRad * sin( pi*real(i-1,kind=pr)/real(max_qrpa_points-1,kind=pr) )
           qrpa_omega_points_im(i) = -qrpa_HcircRad * cos( pi*real(i-1,kind=pr)/real(max_qrpa_points-1,kind=pr) )
         end do
       case default
         STOP "unknown qrpa calculation mode"
       end select
       !
       ! check file io mode. 0: noting is read or written, 1: write amplitudes, 2: read operators
       if((qrpa_file_io.gt.2) .or. (qrpa_file_io.lt.0)) stop "wrong qrpa_file_io"

       ! check that more points than nodes
       if(max_qrpa_points.lt.num_nodes_mpi) then
         WRITE(*,*) "less qrpa points than nodes!"
#ifndef hide_mpi_qrpa
         Call mpi_finalize(icom_err_mpi)
#endif
         STOP
       end if
       !
       ! with mpi or array job, distribute points to various cores. modify starting position and number of points to be calculated in each task
       ! writing of amplitudes is not supported with MPI
       ! NB: improve MPI later! Add dynamic distribution of omega points to MPI tasks in order to reduce reserved CPU hours
       idivm = mod(max_qrpa_points,num_nodes_mpi)
#ifndef hide_mpi_qrpa
       if(qrpa_file_io.eq.1) then
         WRITE(*,*) "Writing of amplitudes in MPI mode not supported"
         Call mpi_finalize(icom_err_mpi)
         STOP
       end if       
       if(iam_mpi.ge.idivm) then  !! last qrpa points. this is also always the case if number of points is exact multiple of number of mpi nodes
          qrpa_points_thiscore = max_qrpa_points/num_nodes_mpi
          qrpa_start_thiscore = (qrpa_points_thiscore+1)*idivm +qrpa_points_thiscore*(iam_mpi-idivm) +1 

       else                       !! first qrpa points
          qrpa_points_thiscore = max_qrpa_points/num_nodes_mpi +1  !! +1 extra point due to modulus
          qrpa_start_thiscore = qrpa_points_thiscore*iam_mpi + 1   !! +1 since table starts from index 1
       end if
#endif

#ifndef hide_array_qrpa
       if(qrpa_file_io.eq.1) then
         WRITE(*,*) "Writing of amplitudes in array mode not supported"
         STOP
       end if
       if(iam_mpi.ge.idivm) then  !! last qrpa points. this is also always the case if number of points is exact multiple of number of mpi nodes
          qrpa_points_thiscore = max_qrpa_points/num_nodes_mpi
          qrpa_start_thiscore = (qrpa_points_thiscore+1)*idivm +qrpa_points_thiscore*(iam_mpi-idivm) +1 

       else                       !! first qrpa points
          qrpa_points_thiscore = max_qrpa_points/num_nodes_mpi +1  !! +1 extra point due to modulus
          qrpa_start_thiscore = qrpa_points_thiscore*iam_mpi + 1   !! +1 since table starts from index 1
       end if
#endif

       Do i=qrpa_start_thiscore,qrpa_start_thiscore+qrpa_points_thiscore-1
        qrpa_omega=Cmplx(qrpa_omega_points_re(i),qrpa_omega_points_im(i),kind=pr) ! set omega as complex number 
        Call qrpa_HFBTHO_solver                                                   ! solve qrpa for this omega
       End Do
       deallocate(qrpa_omega_points_re,qrpa_omega_points_im)
#endif
    Enddo
    !
    Close(901)  !  close the input tho.dat file  
    !
    !help     
    !=============================================================================================
    ! Output variables nucname,ereslbl(1:2),eres(1:ierest) contain:
    !  LBL,BLKN,BLKZ,Jsi,JININ,A,N',Z,Efn,Efp,JEtot,Jbett,Jbetn,Jbetp,JQt,JQn,JQp             
    !  JpEn,JpEp,JpDn,JpDp,JAsn,JAsp,Jrt,Jrn,Jrp,Jrc,Jht,Jhn,Jhp,Jqht,Jqhn,Jqhp,
    !  JKINt,JKINn,JKINp,JSO,JCDIR,JCEX,JDisn,JDisp,JV2Mn,JV2Mp,JILST,JKIND,JL,
    !  JECMPAV1,JECMPAV2,JECMPAV3,JA,JN,JZ,ITER,UEtot,Ubett,Ubetn,Ubetp,UQt,UQn,
    !  UQp,Uln,Ulp,UpEn,UpEp,UpDn,UpDp,UAsn,UAsp,Urt,Urn,Urp,Urc,Uht,Uhn,Uhp,
    !  Uqht,Uqhn,Uqhp,UKINT,UKINN,UKINP,USO,UCDIR,UCEX,UDisn,UDisp,UV2Mn,UV2Mp,
    !  UECMT,UECMN,UECMP,UROTT,UROTN,UROTP,USQUJT,USQUJN,USQUJP,UCRANT,UCRANN,
    !  UCRANP,UERIGT,UERIGN,UERIGP,EHFBLN,EHFB,LNbet,LNben,LNbep,LNQt,LNQn,LNQp,
    !  LNpEn,LNpEp,LNpDn,LNpDp,LNrt,LNrn,LNrp,LNrC,LNam2n,LNam2p,LNe2n,LNe2p,
    !  BlEqpN,BlDEqpN,BlOvrN,BlEqpZ,BlDEqpZ,BlOvrZ
    !=============================================================================================
    !
    ! Example:
    !write(*,*) nucname,ereslbl(1:2),eres(1:ierest)
    !
  End Subroutine Main_Program
!===================================================================================================================================
!#END MAIN PROGRAM
!===================================================================================================================================
#endif
!
!===================================================================================================================================
!#START HFBTHO_SOLVER
!===================================================================================================================================
!===================================================================================================================================
Subroutine HFBTHO_SOLVER  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Universal HFBTHO_SOLVER: 
  !      
  !      Axially deformed configurational constrained and/or unconstrained Hartree-Fock-Bogoliubov
  !      calculations with SkyrmePLUS functionals, and delta pairing using Harmonic-Oscillator (HO), 
  !      and/or Transformed HO (THO) bassis with or without reflection symmetry imposed, with or 
  !      without Lipkin-Nogami procedure. The solver can handle all Skyrme-like functionals, 
  !      DME-functionals, Fayans-functionals, calculating infinite nuclear matter properties, 
  !      finite nuclei (even-even, odd-even, odd-odd), and neutron drops. Isoscalar and isovector 
  !      monopole FAM QRPA calculations for spherical and deformed nuclei.
  !
  !      THE SOLVER ASSUMES THAT ALL *_INI VARIABLES ARE ALREADY AVAILABLE!
  !
  !      The *_INI variables are one-to-one related to the numbers that may come from the input 
  !      file tho.dat and from both name-lists from UNEDF_NAMELIST.DAT
  !       
  !      THE SOVER NEVER CHANGE THEM!
  !
  !      THEREFORE, THE SOLVER CAN BE USED WITHIN TWO DIFFERENT MODELS:
  !
  !      I) MAIN PROGRAM READING ALL REQUIRED *_INI VARIABLES AND THEN CALLING THE SOLVER, i.e.,
  !              Open(901,file='tho.dat',status='old')          
  !              Read(901,*,iostat=ios) n00_INI,b0_INI,q_INI,npr_INI(1),npr_INI(2),skyrme_INI,&
  !                                     kindhfb_INI,inin_INI,cdef_INI,cqad_INI, &
  !                                    ,nkblo_INI(1,1:5),nkblo_INI(2,1:5)
  ! 	         Call read_UNEDF_NAMELIST(skyrme_INI)  
  !	             Call read_HFBTHO_NAMELIST	   
  ! 	         Call HFBTHO_SOLVER
  !              if(ierror_flag.ne.0) write(*,*) ierror_info(1:ierror_flag)
  !
  !      II) SUBROUTINE ACCEPTING *_INI VARIABLES FORM A BUFFER AND THEN CALLING THE SOLVER
  !
  !      In a case of a parralel run the master processor should read tho.dat and BOTH NAMELISTS,
  !      and then should pass all *_INI variables trough a buffer to the slave processors,
  !      or one can pass manually these #_INI parameters trough a buffer.
  !
  !      The complete list of tthe *_INI variables with some exampe values (buffer content)
  !      are the following
  !
  !      ======== from tho.dat 
  !      n00_INI=20;      npr1_INI=70; npr2_INI=50;  kindhfb_INI=-1; inin_INI=-1; 
  !      b0_INI=2.234776; q_INI=0.0;   cdef_INI=0.0; cqad_INI=0.5;   skyrme_INI='SLY4'; nkblo_INI=0
  !
  !      ======== from read_HFBTHO_NAMELIST 
  !      ILST_INI=0;              keypj_INI=1;         iproj_INI=0;    npr1pj_INI=0; 
  !      icou_INI=2;              IDEBUG_INI=0;        npr2pj_INI=0;   
  !      Parity_INI=.False.;      epsi_INI=0.00001_pr; MAX_ITER_INI=101         
  !      Add_Pairing_INI=.False.; DO_FITT_INI=.False.; Print_PTHO_Namelist_INI=.True.;
  !
  !      ======== from read_UNEDF_NAMELIST 
  !      DMEORDER=-1; DMELDA=0; use_TMR_pairing=0
  !      HBZERO=20.73553000000000;    E2CHARG=1.439978400000000;  CRHO(0)=-933.3423749999999;  CRHO(1)=830.0524855000001;
  !      CDRHO(0)=861.0625000000000;  CDRHO(1)=-1064.2732500000;  CTAU(0)=57.12868750000000;   CTAU(1)=24.65673650000000;
  !      CRDR(0)=-76.99620312499999;  CRDR(1)=15.65713512500000;  CRDJ(0)=-92.25000000000000;  CRDJ(1)=-30.7500000000000;
  !      CJ(0)=17.20961150000000;     CJ(1)=64.57581250000000;    CPV0(0)=-258.2000000000000;  CPV0(1)=-258.2000000000000;
  !      CPV1(0)=0.5000000000000000;  CPV1(1)=0.500000000000000;  SIGMA=0.1666666666666667;    CEXPAR=1.000000000000000;
  !      E_NM=-15.97214914144462;     K_NM=229.9009644826037;     SMASS_NM =1.439546988976078; RHO_NM =0.1595387567117334;
  !      ASS_NM =32.00430281505202;   LASS_NM=45.96175148046161;  VMASS_NM =1.249838547196253;                
  !      MPI=0.6995945261023822;      GA=1.290000000000000;       FPI=0.4683223517486062;      C1=-0.1598130000000000; 
  !      C3 =-0.6708200000000;        C4 =0.6708200000000000;     CD =-2.062000000000000;      CE=-0.6250000000000;       
  !      LAMBDAX =3.547896604156107;  USE_INM=.false.;            USE_CM_COR =.true.;          USE_DME3N_TERMS=.true.;      
  !      USE_J2TERMS =.true.;         USE_CHARGE_DENSITY=.false.; PRINT_NAMELIST=.true.;
  !
  ! memo:
  !  -  inin_INI forces scratch unconstraint (inin=1,2,3) or constraint (inin=100,200,300) calculations 
  !     with the value of 'cdef' used. Unconstraint case start with small number constrained iterations.
  !  -  inin_INI forces unconstraint (inin=-1,-2,-3) or constraint (inin=-100,-200,-300) calculations from 
  !     previouse solotion if it exists. If not, the solver sets inin=Abs(inin) and using 'cdef' starts from scratch.
  !  -  The same holds for odd nuclei. If even-even solution for the odd nucleus does not exists it is calculated first.
  !  -  Print_Screen=T/F for n00_INI=+/-  
  !
  !     OUPUT IS GENERATED TO thoout.dat file ONLY IF n00_INI>0. At n00_INI<0 number of shells abs(n00_INI) is used but
  !     all output is supressed.
  !
  !     AT THE END OF THE SOLUTION, THE SOLVER CONIANS ALL RESULTS IN THE ARRAYS nucname,ereslbl(1:2),eres(1:ierest)
  ! memo:     
  !     Output variables nucname,ereslbl(1:2),eres(1:ierest) contain:
  !     LBL,BLKN,BLKZ,Jsi,JININ,A,N',Z,Efn,Efp,JEtot,Jbett,Jbetn,Jbetp,JQt,JQn,JQp             
  !     JpEn,JpEp,JpDn,JpDp,JAsn,JAsp,Jrt,Jrn,Jrp,Jrc,Jht,Jhn,Jhp,Jqht,Jqhn,Jqhp,
  !     JKINt,JKINn,JKINp,JSO,JCDIR,JCEX,JDisn,JDisp,JV2Mn,JV2Mp,JILST,JKIND,JL,
  !     JECMPAV1,JECMPAV2,JECMPAV3,JA,JN,JZ,ITER,UEtot,Ubett,Ubetn,Ubetp,UQt,UQn,
  !     UQp,Uln,Ulp,UpEn,UpEp,UpDn,UpDp,UAsn,UAsp,Urt,Urn,Urp,Urc,Uht,Uhn,Uhp,
  !     Uqht,Uqhn,Uqhp,UKINT,UKINN,UKINP,USO,UCDIR,UCEX,UDisn,UDisp,UV2Mn,UV2Mp,
  !     UECMT,UECMN,UECMP,UROTT,UROTN,UROTP,USQUJT,USQUJN,USQUJP,UCRANT,UCRANN,
  !     UCRANP,UERIGT,UERIGN,UERIGP,EHFBLN,EHFB,LNbet,LNben,LNbep,LNQt,LNQn,LNQp,
  !     LNpEn,LNpEp,LNpDn,LNpDp,LNrt,LNrn,LNrp,LNrC,LNam2n,LNam2p,LNe2n,LNe2p,
  !     BlEqpN,BlDEqpN,BlOvrN,BlEqpZ,BlDEqpZ,BlOvrZ
  !     Example: write(*,*) nucname,ereslbl(1:2),eres(1:ierest)
  !
  ! memo:
  !     Input reads from the file         'tho.dat'
  !     Final solution stored to files '  '*.hel' and/or '*.tel'
  !     Output is written to files        'thoout.dat', 'thodef.dat' and 'hodef.dat'
  !     Output files *.dat may exist as   'thoout','thores','hodenp',
  !                                       'thodenp','thodef','thoene',
  !                                       'thoprc','dat0.1.2.3.4'
  !     External accuracy pr/ipr always come from UNEDF module
  !     INTEGER, PARAMETER :: ipr=KIND(1), pr=KIND(1.0_pr)
  !
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Input file 'tho.dat' has the following srtucture:
  !.................................................................................................................................
  !  20    -2.234776    0.    20   20   'UNEDF1'              -1   -1   0.0    0.5    5 -1 5 1 2    7 -1 3 0 3
  !  20    -2.234776    0.    70   50   'UNEDF0'              -1   -1   0.0    0.5    5 -1 5 1 2    7 -1 3 0 3
  !   0    -2.234776    0.    82   50   '-------------'       -1   -1   0.0    0.5    0 -1 5 1 2    0 -1 3 0 3
  !.................................................................................................................................
  ! (a)      (b)       (c)    (d)  (e)   (f)                  (g)  (h) (i)    (j)     (k)           (l)  
  ! n00       b0      beta0    N    Z    force                LN    C  Cbeta  cquad   blNeutrons    blProtons
  !---------------------------------------------------------------------------------------------------------------------------------
  ! memo: 
  !  (a) n00    Number of oscillator shells
  !              n00>0 prints to thoout.dat & screen
  !              n00<0 no print at all 
  !              n00=0 program stops (NB!)
  !  (b) b0     Oscillator Basis parameter b0>0 (If b0<0 it takes a default value)
  !  (c) beta0  Value of Basis deformation parameter
  !  (d) AN     Number of neutrons N 
  !  (e) AZ     Number of protons Z
  !  (f)'FTST'  Fayance forces label
  !  (g) kind   Kind of calculations 1: noLN, -1:LN
  !  (h) inin   Unconstraint Iterations from peviouse solution
  !              -1: (from spherical *.hel or *.tel file)
  !              -2: (from prolate   *.hel or *.tel file)
  !              -3: (from oblate    *.hel or *.tel file)
  !             Unconstraint Iterations fom scratch with
  !             preliminary constraint at deformation Cbeta, (i):
  !               1: Spherical scratch 
  !               2: Prolate scratch
  !               3: Oblate scratch
  !             Constrained calculation (icstr) at Cbeta, see (i)
  !               100, 200, 300 fom scratch
  !              -100,-200,-300 fom previouse solution
  !  (i) Cbeta Requested deformation used in (h)
  !  (j) Cqad  Spring constraint constant  used in (h)
  !  (k) blNeutrons: a group responsible for blocking a particular neutron level
  !      The group consists of 5 numbers, e.g., for 7-[ 3, 0, 3]: 7 -1 3 0 3
  !       k1  2 \times \Omega 
  !          =0: the whole group (k) is disregarded (n0 blocking)
  !          >0: blocking in N+1 nucleus
  !          <0: blocking in N-1 nucleus
  !       k2  parity (+1 or -1); NB! when k2=0, the ground state walker is applied
  !       k3,k4,k5  Nilson quantum numbers
  !  (l) blProtons: exactly the same as (k) but for protons 
  !=================================================================================================================================
  Use HFBTHO 
  Implicit None 
  Integer(ipr) :: iw,ib,j,i,it 
  Real(pr)     :: epsi0
  !-------------------------------------------------------------
  ! inatializing all according to *_INI values
  !-------------------------------------------------------------
  Call iniialize_HFBTHO_SOLVER
  If(ierror_flag.Ne.0) Return
  If(lout.Lt.lfile) Open(lfile,file='thoout.dat',status='unknown')  
  Call Constraint_or_not(inin_INI,inin,icstr)  
  If(ierror_flag.Ne.0) Return
  !-------------------------------------------------------------------------
  ! Loop recalculating eventually the even-even solution for an odd nucleus
  !-------------------------------------------------------------------------
  Do
     irestart=0
     n00=Abs(n00_INI);  b0=b0_INI;           q=q_INI;           iLST=iLST_INI;
     maxi=MAX_ITER_INI; npr(1)=npr_INI(1);   npr(2)=npr_INI(2); npr(3)=npr(1)+npr(2);
     skyrme=skyrme_INI; kindhfb=kindhfb_INI; cdef=cdef_INI;     cqad=cqad_INI;     
     keypj=keypj_INI;   iproj=iproj_INI;     npr1pj=npr1pj_INI; npr2pj=npr2pj_INI; 
     nkblo=nkblo_INI	
     ! Blocking
     Do it=1,2
	  If(nkblo(it,1).Ne.0) Then
	   If(nkblo(it,1).Gt.0) Then
         ! particle state
	     npr(it)=npr(it)+1
	     iparenti(it)=-1 
	   Else	
         ! hole state
	     npr(it)=npr(it)-1
		 iparenti(it)=+1 
	   End If
	   nkblo(it,1)=Abs(nkblo(it,1)) 
	   If(nkblo(it,2).Eq.0) Then
         ! ground state walker
		 keyblo(it)=nkblo(it,1)
	   End If
	  End If
     End Do
     !-------------------------------------------------------------
     ! HFB+HO calculations
     !-------------------------------------------------------------
     If(ILST.Le.0) Then
        icacou=0; icahartree=0
        Call preparer(.True.)
        If(ierror_flag.Ne.0) Return
        Call inout(1)
		!
        If(ierror_flag.Ne.0) Return
        !-------------------------------------------------------------
        ! Preliminary constrained calculations
        !-------------------------------------------------------------
        If(inin.Gt.0.And.icstr.Eq.0) Then
           epsi0=epsi            ! remember accuracu
           icstr=1               ! constraint true 
           epsi=1.0_pr           ! small accuracu
           Do iw=lout,lfile
              If(Parity) Then
                 Write(iw,'(/,a,i3,a,i2,a,/)') '  ### INITIAL STAGE(constrained calculations, reflection symmetry used)' 
              Else
                 Write(iw,'(/,a,i3,a,i2,a,/)') '  ### INITIAL STAGE(constrained calculations, no reflection symmetry used)' 
              Endif
           Enddo
           Call iter(.True.)     ! small constraint iterations            
           If(ierror_flag.Ne.0) Return
           icstr=0;              ! requested unconstraint calculation
           epsi=epsi0            ! return to requested requested accuracy
        Endif
        !-------------------------------------------------------------
        ! REGULAR HFB+HO ITERATIONS
        !-------------------------------------------------------------
        Do iw=lout,lfile
           If(Parity) Then
              Write(iw,'(/,a,i3,a,i2,a,/)')    '  ### REGULAR STAGE (reflection symmetry imposed)'
           Else
              Write(iw,'(/,a,i3,a,i2,a,/)')    '  ### REGULAR STAGE (no reflection symmetry imposed)'           
           Endif
        Enddo
        Call iter(.True.)  
        If(ierror_flag.Ne.0) Return
        Call resu(1)        
        If(ierror_flag.Ne.0) Return
     End If
#ifndef hide_tho
     !-------------------------------------------------------------
     ! HFB+THO calculations from HFB+HO
     !-------------------------------------------------------------
     If(ILST.Lt.0) Then
        ILST1=1; icacou=0; icahartree=0
        Call coordinateLST(.False.) ! THO basis
        If(ierror_flag.Ne.0) Return
        Call densit                 ! THO densities
        If(ierror_flag.Ne.0) Return
        Call field                  ! Nuclear fields
        If(ierror_flag.Ne.0) Return		
        Call iter(.True.)           ! HFB+THO iterations
        If(ierror_flag.Ne.0) Return
        Call resu(1)                ! print/record results
        If(ierror_flag.Ne.0) Return
     End If
     !-------------------------------------------------------------
     ! HFB+THO calculations from *.tel
     !-------------------------------------------------------------
     If(ILST.Gt.0) Then
        If(inin.Gt.0) Then
           ierror_flag=ierror_flag+1
           ierror_info(ierror_flag)=' Stop: Forbidden iLST>0, inin>0 (try inin<0 if the file *.tel exists)' 
           Return 
        Endif
        icacou=0; icahartree=0
        Call preparer(.True.)
        If(ierror_flag.Ne.0) Return
        Call inout(1)               ! reading HFB matrices
        If(ierror_flag.Ne.0) Return
        Call iter(.True.)           ! HFB+THO iterations
        If(ierror_flag.Ne.0) Return
        Call resu(1)                ! print/record results
        If(ierror_flag.Ne.0) Return
     End If
#endif
     !-------------------------------------------------------------
     ! go for the requested blocking state in a case of odd nuclei
     ! if restarted due to corrupted/missing previouse solution
     !-------------------------------------------------------------
     inin=-Abs(inin)                
     If(irestart.Eq.0) Exit
  End Do
  ! 
End Subroutine HFBTHO_SOLVER
!=======================================================================
!
!=======================================================================
Subroutine heading
  !---------------------------------------------------------------------
  ! print heading to screen 'lout' and to tape thoout.dat 'lfile'
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: iw,idt(8)
  Character(len=12) rcl(3)
  Character(len=50) today
  !
  Call Date_and_time(rcl(1),rcl(2),rcl(3),idt)
  Write(today,'(a,i2,a,i2,a,i4,a,i2,a,i2,a)')'(',idt(2),'/',idt(3),'/',idt(1),', ',idt(5),':',idt(6),')'
  Do iw=lout,lfile
     Write(iw,'(a)')
     Write(iw,'(a)')      '  ======================================='
     Write(iw,'(a,i2,a)') '           FORTRAN 95 CODE (KIND=',pr,') '
     Write(iw,'(a,a)')    '               Version: ',Version
     Write(iw,'(a)')      '  ======================================='
     Write(iw,'(a)')      '       AXIALLY DEFORMED CONFIGURATIONAL  '
     Write(iw,'(a)')      '     HARTREE-FOCK-BOGOLIUBOV CALCULATIONS'
     Write(iw,'(a)')      '                     WITH                '
     Write(iw,'(a)')      '            UNEDF AND DELTA PAIRING      '
     Write(iw,'(a)')      '                     USING               '
     Write(iw,'(a)')      '             HARMONIC-OSCILLATOR         '
     Write(iw,'(a)')      '                   AND/OR                '
     Write(iw,'(a)')      '       TRANSFORMED HARMONIC-OSCILLATOR   '
     Write(iw,'(a)')      '                    BASIS                '
     Write(iw,'(a)')      '                     ---                 '
     Write(iw,'(a)')      '   Stoitsov,Dobaczewski,Nazarewicz,Ring  ' 
     Write(iw,'(a)')      '  ======================================='
     Write(iw,'(a,a,a,i4,a,i3,a,i3,a)')'    Nucleus: ',nucname,' (A=',npr(1)+npr(2),', N=',npr(1),', Z=',npr(2),')'
     If(Parity) Then
        Write(iw,'(a)')     '       Reflection Symmetry Imposed       '
     Else
        Write(iw,'(a)')     '      No Reflection Symmetry Imposed     '
     Endif
     Write(iw,'(a,a)')    '            ',today
     Write(iw,'(a)')      '  ======================================='
     Write(iw,'(a)')
#ifdef hide_qrpa    
     Write(iw,'(a)')      '  Preprocessing: hide_qrpa'
#endif     
#ifdef hide_tho    
     Write(iw,'(a)')      '  Preprocessing: hide_tho'
#endif     
#ifdef hide_dme    
     Write(iw,'(a)')      '  Preprocessing: hide_dme'
#endif     
     Write(iw,'(a)')
  End Do
End Subroutine heading
!=======================================================================
!
!=======================================================================
Subroutine thodefh(iw1)
  !---------------------------------------------------------------------
  ! print labels to hodef.dat or/and thodef.dat files
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None 
  Integer(ipr) :: iw1
  hlabels(1)='LBL';       hlabels(11)='JEtot';    hlabels(21)='JpDp'; 
  hlabels(2)='BLKN';      hlabels(12)='Jbett';    hlabels(22)='JAsn'; 
  hlabels(3)='BLKZ';      hlabels(13)='Jbetn';    hlabels(23)='JAsp';  
  hlabels(4)='Jsi';       hlabels(14)='Jbetp';    hlabels(24)='Jrt'; 
  hlabels(5)='JININ';     hlabels(15)='JQt';      hlabels(25)='Jrn';  
  hlabels(6)='A';         hlabels(16)='JQn';      hlabels(26)='Jrp';    
  hlabels(7)='N';         hlabels(17)='JQp';      hlabels(27)='Jrc';   
  hlabels(8)='Z';         hlabels(18)='JpEn';     hlabels(28)='Jht';    
  hlabels(9)='Efn';       hlabels(19)='JpEp';     hlabels(29)='Jhn';    
  hlabels(10)='Efp';      hlabels(20)='JpDn';     hlabels(30)='Jhp';    
  !
  hlabels(31)='Jqht';    hlabels(41)='JDisp';    hlabels(51)='JN';          
  hlabels(32)='Jqhn';    hlabels(42)='JV2Mn';    hlabels(52)='JZ';            
  hlabels(33)='Jqhp';    hlabels(43)='JV2Mp';    hlabels(53)='ITER';          
  hlabels(34)='JKINt';   hlabels(44)='JILST';    hlabels(54)='UEtot';          
  hlabels(35)='JKINn';   hlabels(45)='JKIND';    hlabels(55)='Ubett';   
  hlabels(36)='JKINp';   hlabels(46)='JL';       hlabels(56)='Ubetn';          
  hlabels(37)='JSO';     hlabels(47)='JECMPAV1'; hlabels(57)='Ubetp';           
  hlabels(38)='JCDIR';   hlabels(48)='JECMPAV2'; hlabels(58)='UQt';            
  hlabels(39)='JCEX';    hlabels(49)='JECMPAV3'; hlabels(59)='UQn';     
  hlabels(40)='JDisn';   hlabels(50)='JA';       hlabels(60)='UQp';            
  !
  hlabels(61)='Uln';     hlabels(71)='Urp';      hlabels(81)='UKINP';             
  hlabels(62)='Ulp';     hlabels(72)='Urc';      hlabels(82)='USO';       
  hlabels(63)='UpEn';    hlabels(73)='Uht';      hlabels(83)='UCDIR'; 
  hlabels(64)='UpEp';    hlabels(74)='Uhn';      hlabels(84)='UCEX';    
  hlabels(65)='UpDn';    hlabels(75)='Uhp';      hlabels(85)='UDisn';      
  hlabels(66)='UpDp';    hlabels(76)='Uqht';     hlabels(86)='UDisp';   
  hlabels(67)='UAsn';    hlabels(77)='Uqhn';     hlabels(87)='UV2Mn'; 
  hlabels(68)='UAsp';    hlabels(78)='Uqhp';     hlabels(88)='UV2Mp';    
  hlabels(69)='Urt';     hlabels(79)='UKINT';    hlabels(89)='UECMT';    
  hlabels(70)='Urn';     hlabels(80)='UKINN';    hlabels(90)='UECMN'; 
  !
  hlabels(91)='UECMP';   hlabels(101)='UERIGT';  hlabels(111)='LNQp';     
  hlabels(92)='UROTT';   hlabels(102)='UERIGN';  hlabels(112)='LNpEn';         
  hlabels(93)='UROTN';   hlabels(103)='UERIGP';  hlabels(113)='LNpEp';         
  hlabels(94)='UROTP';   hlabels(104)='EHFBLN';  hlabels(114)='LNpDn';        
  hlabels(95)='USQUJT';  hlabels(105)='EHFB';    hlabels(115)='LNpDp';   
  hlabels(96)='USQUJN';  hlabels(106)='LNbet';   hlabels(116)='LNrt';          
  hlabels(97)='USQUJP';  hlabels(107)='LNben';   hlabels(117)='LNrn';          
  hlabels(98)='UCRANT';  hlabels(108)='LNbep';   hlabels(118)='LNrp';          
  hlabels(99)='UCRANN';  hlabels(109)='LNQt';    hlabels(119)='LNrC';    
  hlabels(100)='UCRANP'; hlabels(110)='LNQn';    hlabels(120)='LNam2n';        
  !
  hlabels(121)='LNam2p';          
  hlabels(122)='LNe2n';          
  hlabels(123)='LNe2p'; 
  hlabels(124)='BlEqpN';        
  hlabels(125)='BlDEqpN';        
  hlabels(126)='BlOvrN';       
  hlabels(127)='BlEqpZ'; 
  hlabels(128)='BlDEqpZ';        
  hlabels(129)='BlOvrZ'; 
  ! 
  Write(iw1,'(660(a,1x))') hlabels
  !
  ! HELP
  !Do i=1,129
  ! write(iw1,'(1x,i3,a,a)',advance='NO') i,':',trim(hlabels(i))
  !Enddo
  ! 1:LBL  2:BLKN  3:BLKZ  4:Jsi  5:JININ  6:A  7:N  8:Z  9:Efn 10:Efp 
  ! 11:JEtot 12:Jbett 13:Jbetn 14:Jbetp 15:JQt 16:JQn 17:JQp 18:JpEn 19:JpEp 20:JpDn 
  ! 21:JpDp 22:JAsn 23:JAsp 24:Jrt 25:Jrn 26:Jrp 27:Jrc 28:Jht 29:Jhn 30:Jhp 
  ! 31:Jqht 32:Jqhn 33:Jqhp 34:JKINt 35:JKINn 36:JKINp 37:JSO 38:JCDIR 39:JCEX 40:JDisn 
  ! 41:JDisp 42:JV2Mn 43:JV2Mp 44:JILST 45:JKIND 46:JL 47:JECMPAV1 48:JECMPAV2 49:JECMPAV3 50:JA 
  ! 51:JN 52:JZ 53:ITER 54:UEtot 55:Ubett 56:Ubetn 57:Ubetp 58:UQt 59:UQn 60:UQp 
  ! 61:Uln 62:Ulp 63:UpEn 64:UpEp 65:UpDn 66:UpDp 67:UAsn 68:UAsp 69:Urt 70:Urn 
  ! 71:Urp 72:Urc 73:Uht 74:Uhn 75:Uhp 76:Uqht 77:Uqhn 78:Uqhp 79:UKINT 80:UKINN 
  ! 81:UKINP 82:USO 83:UCDIR 84:UCEX 85:UDisn 86:UDisp 87:UV2Mn 88:UV2Mp 89:UECMT 90:UECMN 
  ! 91:UECMP 92:UROTT 93:UROTN 94:UROTP 95:USQUJT 96:USQUJN 97:USQUJP 98:UCRANT 99:UCRANN 100:UCRANP
  ! 101:UERIGT 102:UERIGN 103:UERIGP 104:EHFBLN 105:EHFB 106:LNbet 107:LNben 108:LNbep 109:LNQt 110:LNQn
  ! 111:LNQp 112:LNpEn 113:LNpEp 114:LNpDn 115:LNpDp 116:LNrt 117:LNrn 118:LNrp 119:LNrC 120:LNam2n
  ! 121:LNam2p 122:LNe2n 123:LNe2p 124:BlEqpN 125:BlDEqpN 126:BlOvrN 127:BlEqpZ 128:BlDEqpZ 129:BlOvrZ
End Subroutine thodefh
!=======================================================================
!
!=======================================================================
Subroutine thoalloc
  !---------------------------------------------------------------------
  ! Allocates arrays at given number of oscillator shells 'n00'
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer :: ier,ib,ND
  !    
  ! number of int.points      
  If(Parity) Then
!     ngh=30; ngl=30; nleg=30                                    !Yesp
     !value from CPC(2013)
     ngh=40; ngl=40; nleg=80                                    !Nop
  Else  
!     ngh=60; ngl=30; nleg=30                                    !Nop
     !value from CPC(2013)
     ngh=80; ngl=40; nleg=80                                    !Nop
  Endif
  !
  !nbx=2*n00+1                   ! maximal number of k-blocks
  !ntx=(n00+1)*(n00+2)*(n00+3)/6 ! max.num. p/n levels
  !nzx=n00                       ! maximal nz-quantum number
  !nrx=n00/2+1                   ! maximal nr-quantum number
  !nlx=n00                       ! maximal ml-quantum number
  !ndx=(n00+2)*(n00+2)/4         ! maximal dim. of one k-block
  !nhhdim=number of nonzero HH matrix elements 
  !
  nzrlx=(nzx+1)*(nrx+1)*(nlx+1)   ! phy(:,:,nzrlx)
  nghl=ngh*ngl                    ! nghl=ngh*ngl 
  nqx=ndx*ndx; nb2x=nbx+nbx; ndx2=ndx+ndx
  ilnqx=ilpj*nqx; ilnghl=ilpj*nghl
  nhfbx=ndx+ndx; nhfbqx=nhfbx*nhfbx; nkx=ntx; ndxs=ndx*(ndx+1)/2 
  !-----------------------------------------
  !Arrays depending on gauss points
  !-----------------------------------------
  If(Allocated(xh)) Deallocate(xh,wh,xl,sxl,wl,xleg,wleg,vc &
       ,vhbn,vn,vrn,vzn,vdn,vsn,dvn,vhbp,vp,vrp,vzp,vdp,vsp,dvp  &
       ,vSZFIn,vSFIZn,vSRFIn,vSFIRn,vSZFIp,vSFIZp,vSRFIp,vSFIRp &
       ,fl,fli,fh,fd,fp1,fp2,fp3,fp4,fp5,fp6  &
       ,fs1,fs2,fs3,fs4,fs5,fs6,wdcor,wdcori,cou,vDHartree,vhart00,vhart01,vhart11)
  Allocate(xh(ngh),wh(ngh),xl(ngl),sxl(ngl),wl(ngl),xleg(nleg),wleg(nleg),vc(nghl,nghl))
  Allocate(vhbn(nghl),vn(nghl),vrn(nghl),vzn(nghl),vdn(nghl),vsn(nghl),dvn(nghl)  &
       ,vhbp(nghl),vp(nghl),vrp(nghl),vzp(nghl),vdp(nghl),vsp(nghl),dvp(nghl)  &         
       ,vSZFIn(nghl),vSFIZn(nghl),vSRFIn(nghl),vSFIRn(nghl)  &
       ,vSZFIp(nghl),vSFIZp(nghl),vSRFIp(nghl),vSFIRp(nghl))
  Allocate(fl(nghl),fli(nghl),fh(nghl),fd(nghl),fp1(nghl),fp2(nghl),fp3(nghl)  &
       ,fp4(nghl),fp5(nghl),fp6(nghl),fs1(nghl),fs2(nghl),fs3(nghl),fs4(nghl)  &
       ,fs5(nghl),fs6(nghl),wdcor(nghl),wdcori(nghl),cou(nghl),vDHartree(nghl,2) &
       ,vhart00(nghl,nghl),vhart01(nghl,nghl),vhart11(nghl,nghl))
  If(Allocated(aka)) Deallocate(aka,ro,tau,dro,dj,NABLAR,NABLAZ,SZFI,SFIZ,SRFI,SFIR)
  Allocate(aka(nghl,2),ro(nghl,2),tau(nghl,2),dro(nghl,2),dj(nghl,2)  & 
       ,SZFI(nghl,2),SFIZ(nghl,2),SRFI(nghl,2),SFIR(nghl,2)  &
       ,NABLAR(nghl,2),NABLAZ(nghl,2))
  !-----------------------------------------
  ! Arrays depending on configurations
  !-----------------------------------------
  If(Allocated(rk)) Deallocate(rk,ak,qh,qh1,ql,ql1,nz,nr,nl,ns,npar,id  &
       ,ia,ikb,ipb,ka,kd,tb,txb,numax,ek,dk,vk,vk1,uk,vkmax,ddc,ddc1,hfb1,lcanon)
  Allocate(rk(nqx,nb2x),ak(nqx,nb2x),qh(0:nzx,1:ngh+1)  &
       ,qh1(0:nzx,1:ngh+1),ql(0:nrx,0:nlx,1:ngl+1),ql1(0:nrx,0:nlx,1:ngl+1)  &
       ,nz(ntx),nr(ntx),nl(ntx),ns(ntx),npar(ntx),id(nbx),ia(nbx),ikb(nbx),lcanon(0:nbx,2)  &
       ,ipb(nbx),ka(nbx,2),kd(nbx,2),tb(ntx),txb(nbx),numax(0:nkx,2)  &
       ,ek(nkx,2),dk(nkx,2),vk(nkx,2),vk1(nkx,2),uk(nkx,2),vkmax(nkx,2)  &
       ,ddc(ndx,nkx,2),ddc1(ndx,nkx,2),hfb1(nhfbx,2))
  !-----------------------------------------
  ! HFB Arrays
  !-----------------------------------------
  If(Allocated(erhfb)) Deallocate(erhfb,drhfb,erhfb1,drhfb1)
  Allocate(erhfb(nkx),drhfb(nkx),erhfb1(nkx),drhfb1(nkx))
  If(Allocated(hfb)) Deallocate(hfb,zhfb,evvk,hfbcan,evvkcan)
  Allocate(hfb(ndx2,ndx2),zhfb(ndx2),evvk(ndx2),hfbcan(ndx,ndx),evvkcan(ndx))
  If(Allocated(AN)) Deallocate(AN,ANk,PFIU,PFID,FIU,FID,FIUR,FIDR,FIUD2N,FIDD2N,FIUZ,FIDZ)
  Allocate(AN(nqx),ANk(nqx),PFIU(ndx),PFID(ndx),FIU(ndx),FID(ndx)  &
       ,FIUR(ndx),FIDR(ndx),FIUD2N(ndx),FIDD2N(ndx),FIUZ(ndx),FIDZ(ndx))
  !-----------------------------------------
  ! Optimal LAPACK storage
  !-----------------------------------------
  ialwork=1; ilwork=1; 
  If(Allocated(alwork)) Deallocate(alwork,lwork)
  Allocate(alwork(ialwork),lwork(ilwork))
  ier=0; Call DSYEVD('V','L',ndx2,hfb,ndx2,evvk,ALWORK,-1,LWORK,-1,ier)
  If(ier.Ne.0) Then
     ierror_flag=ierror_flag+1
     ierror_info(ierror_flag)='STOP: FATAL ERROR CONDITION IN DSYEVD'
     Return
  Endif
  ialwork=Int(alwork(1)); ilwork=lwork(1)
  If(Allocated(alwork)) Deallocate(alwork,lwork)
  Allocate(alwork(ialwork),lwork(ilwork))
  !-----------------------------------------
  ! Eqp, U,V 
  !-----------------------------------------
  If(Allocated(RVqpN)) Deallocate(RVqpN,RVqpP,RUqpN,RUqpP,REqpN,REqpP)
  Allocate(RVqpN(nuv),RVqpP(nuv),RUqpN(nuv),RUqpP(nuv),REqpN(nqp),REqpP(nqp))
  If(Allocated(KpwiP)) Deallocate(KpwiP,KpwiN,KqpN,KqpP)
  Allocate(KpwiN(nqp),KpwiP(nqp),KqpN(nqp),KqpP(nqp))
  !-----------------------------------------
  ! PNP ARRAYS: CONF. AND GAUGE ANGLE
  !-----------------------------------------
  If(Allocated(exp1iphy))Deallocate(ropj,taupj,dropj,djpj,akapj,coupj,pjk  &
       ,SZFIpj,SFIZpj,SRFIpj,SFIRpj,epj,cpj,ypj,rpj,ddepj,phypj,sinphy  &
       ,exp1iphy,exp2iphy,exp1iphym,exp2iphym)
  Allocate(ropj(nghl,ilpj,2),taupj(nghl,ilpj,2),dropj(nghl,ilpj,2)  &
       ,djpj(nghl,ilpj,2),akapj(nghl,ilpj,2),coupj(nghl,ilpj),pjk(ilpj,2)  &
       ,SZFIpj(nghl,ilpj,2),SFIZpj(nghl,ilpj,2),SRFIpj(nghl,ilpj,2)  &
       ,SFIRpj(nghl,ilpj,2),epj(ilpj,2),cpj(nkx,ilpj,2),ypj(nkx,ilpj,2)  &
       ,rpj(nkx,ilpj,2),ddepj(nqx,ilpj,nb2x),phypj(ilpj),sinphy(ilpj),  &
       exp1iphy(ilpj),exp2iphy(ilpj),exp1iphym(ilpj),exp2iphym(ilpj))
  !-----------------------------------------
  ! FIELDS INITIALIZATION (NB! optimize)
  !-----------------------------------------
  ro=zero;     tau=zero;    dro=zero;    dj=zero;  aka=zero; rk=zero
  vn=zero;     vsn=zero;    vhbn=zero;   vrn=zero; vzn=zero; vdn=zero; 
  vp=zero;     vsp=zero;    vhbp=zero;   vrp=zero; vzp=zero; vdp=zero; 
  vSFIZn=zero; vSZFIn=zero; vSFIRn=zero; vSRFIn=zero;  vDHartree=zero;
  vSFIZp=zero; vSZFIp=zero; vSFIRp=zero; vSRFIp=zero; 
End Subroutine thoalloc
!=======================================================================
!
!=======================================================================
Subroutine preparer(lpr)
  !---------------------------------------------------------------------
  ! setup routine
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Logical      :: lpr
  Integer(ipr) :: iw
  !
  If(n00.Eq.0) Then
     ierror_flag=ierror_flag+1
     ierror_info(ierror_flag)=' STOP: No more nuclei pass to the solver'
     Return
  Endif
  !-----------------------------------------
  ! select the symbol of the nucleus
  !-----------------------------------------
  Call nucleus(1,npr(2),nucname)   
  If(ierror_flag.Ne.0) Return  
  !-----------------------------------------
  ! print headings to screen/'thoout.dat'
  !-----------------------------------------
  If(lpr) Then
     Call heading             
     Do iw=lout,lfile
        Call print_functional_parameters(iw)
        If(ierror_flag.Ne.0) Return		
        If(Print_HFBTHO_Namelist) Then
           Write(iw,'(100(2x,a,f15.8))')
           Write(iw,'(100(2x,a,f15.8))') 'NAMELIST CONTENT (copy/past to UNEDF_NAMELIST.DAT and modify)'
           Write(iw,'(100(2x,a,f15.8))') '-------------------------------------------------------------'      
           Write(iw,HFBTHO_NAMELIST) 
        End If
     End Do
  End If
  !-----------------------------------------
  ! pairing parameters (NB! modify later)
  !-----------------------------------------
  !rho_c=0.160_pr; pwi=2000.0_pr !60.0_pr; 
  rho_c=0.160_pr; pwi=60.0_pr; 
  !-----------------------------------------
  ! particle number as real variable
  !-----------------------------------------
  tz(1)=Dble(npr(1)); tz(2)=Dble(npr(2)); amas=tz(1)+tz(2)
  !-----------------------------------------
  ! default combinations
  !-----------------------------------------
  r00=r0*amas**p13; r02=r00**2; r04=r02**2
  chargee2 = e2charg
  coex=-chargee2*(3.0_pr/pi)**p13; cex=-0.750_pr*coex
  hom=41.0_pr*amas**(-p13)*r0
  !-----------------------------------------
  ! hbzero from forces [hqc**2/(two*amu)] 
  !-----------------------------------------
  hb0=hbzero; If (use_cm_cor) hb0=hb0*(one-one/amas) 
  !-----------------------------------------
  ! basis parameter q
  !-----------------------------------------
  beta0=q; q=Exp((3.0_pr*Sqrt(5.0_pr/(16.0_pr*pi)))*beta0)
  !-----------------------------------------
  ! basis parameters b0,bp,bz
  !-----------------------------------------
  If(b0.Le.zero) b0=Sqrt(two*hbzero/hom)
  If(b0.Le.zero) b0=Sqrt(hqc**2/(hom*amu))
  bp=b0*q**(-one/6.0_pr); bz=b0*q**(one/3.0_pr); bpp=bp*bp
  !-----------------------------------------
  ! constraint in terms of beta
  !-----------------------------------------
  ty20=Sqrt(5.0_pr/pi)*hom/b0**2/two     
  !-----------------------------------------
  ! projection: number of grid points 
  !-----------------------------------------
  keypj=Max(1,keypj); ilpj=keypj; ilpj2=ilpj**2; 
  !-----------------------------------------
  ! projecting on different nucleus
  !-----------------------------------------
  If(iproj.Eq.0) Then
     npr1pj=npr(1); npr2pj=npr(2)
  Else
     npr1pj=npr(1)+npr1pj; npr2pj=npr(2)+npr2pj
  End If
  !-----------------------------------------
  ! blocking window
  !-----------------------------------------
  pwiblo=Min(Max(25.0_pr/Sqrt(Float(npr(1)+npr(2))),2.0_pr),8.0_pr)
  !---------------------------------------------------------
  ! statistics to screen('lout')/file('lfile')
  !---------------------------------------------------------
  If(lpr) Then
     Do iw=lout,lfile
        Write(iw,*)
        Write(iw,'(a,i5)')          '  output file device for lfile: ',lfile
        Write(iw,'(a,2i5)')         '  number of oscillator shells : ',n00
        Write(iw,'(a,2f10.3)')      '  ho  length b0 (fm)          : ',b0
        Write(iw,'(a,f8.3,a,f8.3)') '  basis deformation           :  beta0=',beta0,' q=',q
        If(iLST.Eq.0)  Then         ! HFB+HO case only
           iLST1=0
           Write(iw,'(a,2f10.3)')   '  tho basis is                :  off'
        Else                        ! HFB+THO case
           Write(iw,'(a,2f10.3)')   '  tho basis is                :  on'
           If(iLST.Gt.0) Then       ! HFB+THO  only
              iLST1=1
              If(inin.Gt.0) Then
                 ierror_flag=ierror_flag+1
                 ierror_info(ierror_flag)=' Stop: Forbidden iLST>0, inin>0 combination. Try iLST>0, inin<0 If the file *.hel exists'
                 Return
              End If
              Write(iw,'(a,2f10.3)')'  THO parameters from tholst.wel'
           Else                     ! HFB+THO after HFB+HO
              iLST1=0
              Write(iw,'(a,2f10.3)')'  HFB+THO after a HFB+HO run    '
           End If
        End If
        Write(iw,'(a,2i5)')         '  maximal number of iterations: ',maxi
        Write(iw,'(a,2f10.3)')      '  initial mixing parameter    : ',xmix
        If(inin.Eq.1)  Then
           Write(iw,'(a,2i5)')      '  initial w.f.                :  from  &
                & spherical scratch'
        End If
        If(inin.Eq.2)  Then
           Write(iw,'(a,2i5)')      '  initial w.f.                :  from prolate scratch'
        End If
        If(inin.Eq.3)  Then
           Write(iw,'(a,2i5)')      '  initial w.f.                :  from oblate scratch'
        End If
        If(inin.Lt.0) Then
           Write(iw,'(a,2i5)')      '  initial wave functions from :  tape'
        End If
        Write(iw,'(a,2x,a2,i4)')    '  nucleus:                    : ',nucname,npr(1)+npr(2)
        Write(iw,'(a,3x,a)')        '  skyrme forces               : ',skyrme
        If(kindhfb.Lt.0) Then
           Write(iw,'(a,2i5)')      '  Lipkin-Nogami procedure is  :  ON'
        Else
           Write(iw,'(a,2i5)')      '  Lipkin-Nogami procedure is  :  OFF'
        End If
        If(ilpj-1.Eq.0) Then
           Write(iw,'(a,2i5)')      '  PAV procedure is OFF, lpj   : ',keypj
        Else
           Write(iw,'(a,2i5)')      '  PAV(ON) #gauge-points, lpj  : ',keypj
        End If
        If(icstr.Eq.0) Then
           Write(iw,'(a,2i5)')      '  constraint calculation is   :  OFF'
        Else
           Write(iw,'(a,2i5)')      '  constraint calculation is   :   ON'
           Write(iw,'(a,2f10.3)')   '  requested def value   cdef  :  ',cdef
           Write(iw,'(a,2f10.3)')   '  quadratic.constr.const.cqad :  ',cqad
        End If
        If(keyblo(1).Ne.0) Then
           Write(iw,'(a,2f10.3)')   '         Neutron blocking is  :   ON'
        End If
        If(keyblo(2).Ne.0) Then
           Write(iw,'(a,2f10.3)')   '          Proton blocking is  :   ON'
        End If
        Write(iw,'(a,I3)')          '                       ILST   :  ',iLST
        Write(iw,'(a,I3)')          '                       ININ   :  ',inin
        If(nbroyden.Eq.0) Then
           Write(iw,'(a,I3)')       '                Linear mixing :  ',nbroyden
        Else
           Write(iw,'(a,I3)')       '               Broyden mixing :  ',nbroyden
        Endif
        Write(iw,'(a,I3)')   
        !
        Write(iw,'(a,a,i4,i6,i4)') '  nucleus: ',nucname,npr(1)+npr(2),npr(1),npr(2)
        Write(iw,'(a,5(1x,e14.8))')'  b0,1/b0,bp,bz,q=',b0,one/b0,bp,bz,q
        Write(iw,'(a,5(1x,e14.8))')'  hb0=h**2/(2m)  & cmc  & e**2=',hbzero,hb0,chargee2
        Write(iw,'(a,5(1x,e14.8))')'  hom=r0*41.0_pr*A^{-1/3}=',hom
        If (icou.Eq.0) Write(iw,'(a,5(1x,e14.8))')    '  without coulomb forces'
        If (icou.Eq.1) Write(iw,'(a,5(1x,e14.8))')    '  with direct coulomb force only'
        If (icou.Eq.2) Write(iw,'(a,5(1x,e14.8))')    '  with coulomb force with exchange'
        Write(iw,'(a,5(1x,e14.8))')
     End Do
  Endif
  !-----------------------------------------
  ! BASIS, GAUSS POINTS, HOWF
  !-----------------------------------------
  Call gfv                      ! factorials    
  If(ierror_flag.Ne.0) Return     
  Call base0(lpr)               ! basis space (calculate configurational space)
  If(ierror_flag.Ne.0) Return  
  Call thoalloc                 ! global allocation
  If(ierror_flag.Ne.0) Return  
  Call gausspoints              ! GAUSS mesh points
  If(ierror_flag.Ne.0) Return  
  Call base(lpr)                ! oscillator configurations (set up quantum numbers)
  If(ierror_flag.Ne.0) Return  
  Call gaupol(lpr)              ! basis wf at gauss mesh points   
  If(ierror_flag.Ne.0) Return  
  !
End Subroutine preparer
!====================================================================
!
!====================================================================
Subroutine coordinateLST(lpr)
  !------------------------------------------------------------------
  ! HO/THO
  !------------------------------------------------------------------
  Use HFBTHO
#ifndef hide_tho
  Use THO_MODULE, Only: f01234    
#endif
  Implicit None
  Logical :: lpr   
  Integer(ipr) :: i,il,ih
  If(iLST1.Eq.0) Then
     ! HO-basis
     Do il=1,ngl
        Do ih=1,ngh
           i=ih+(il-1)*ngh
           fh(i)=bz*xh(ih)
           fl(i)=bp*Sqrt(xl(il))
           wdcor(i)=pi*wh(ih)*wl(il)*bz*bp*bp
           wdcori(i)=one/wdcor(i)
        End Do
     End Do
#ifndef hide_tho 
  Else
     ! THO basis
     Call f01234(.False.)
     If(ierror_flag.Ne.0) Return	 
#endif               
  End If
  !
  Call optHFBTHO                ! optimal HO/THO combinations
  If(ierror_flag.Ne.0) Return  
  !  
End Subroutine coordinateLST
!====================================================================
!
!====================================================================
Subroutine iter(lpr) 
  !------------------------------------------------------------------
  ! Iterations through successive diagonalisation
  !------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Logical :: lpr
  Real(pr)       :: assprn,delln(2)
  Real(pr), Save :: time
  Integer(ipr)   :: iw,it,ite
  Real(pr)       :: time1,time2,time3
  !---------------------------------------------------
  ! print to screen('lout')/thoout.dat('lfile')
  !---------------------------------------------------
  Do iw=lout,lfile
     If(iLST.Eq.0) Then
        Write(iw,'(a,f7.3,4(a,i3),a)')  & 
             '  |HFB+HO> iterations(b0=',b0,', Nsh=',n00,  & 
             ', inin=',inin,', N=',npr(1),', Z=',npr(2),')...'
     Else
        If(iLST1.Eq.0.Or.iasswrong(3).Ne.0) Then
           If(iasswrong(3).Ne.0) Then
              Write(iw,'(a,f7.3,a,i3,a)')  & 
                   '  |HFB+THO substituted by HFB+HO> iterations (b0=',  & 
                   b0,', Nsh=',n00,')...'
           Else
              Write(iw,'(a,f7.3,a)')'  towards |hfb+tho> iterations...'
              Write(iw,'(a,f7.3,a)')
              Write(iw,'(a,f7.3,a,i3,a)')  & 
                   '  |Preliminary HFB+HO> iterations (b0=',b0,', Nsh=',n00,')...'
           End If
        Else
           If(itass.Eq.1) Then
              Write(iw,'(2(a,f7.3),a,i3,a)')  & 
                   '  |HFB+THO> iterations(b0=',b0,', neutron density decay=',  & 
                   decay,', Nsh=',n00,')...'
           Else
              Write(iw,'(2(a,f7.3),a,i3,a)')  & 
                   '  |HFB+THO> iterations(b0=',b0,', proton density decay=',  & 
                   decay,', Nsh=',n00,')...'
           End If
        End If
     End If
     Write(iw,1)
     Write(iw,'(20(a))')'  i','          si ','    mix ','  beta', '  &
          &  Etot ','    A ','      rn','      rp ','        En', '   &
          &    Dn','      Ep','      Dp','        Ln  ','   Lp ', '   &
          &    time' !Idro '
     Write(iw,1)
1    Format(2x,130('-'))       
  End Do
  !---------------------------------------------------------------------
  ! main hfb iteration loop
  !---------------------------------------------------------------------
  iError_in_HO=0; iError_in_THO=0; time=0.0_pr; time3=0.0_pr 
  Do ite=1,maxi   
     Call Cpu_time(time1)
     !
     iiter=ite
     !
     If (lpr.Or.iiter.Eq.1) Then
        assprn=ass(1); If(assprn.Gt.ass(2)) assprn=-ass(2) ! protons come with '-'                  
        delLN=del; If(kindhfb.Lt.0) delLN=del+ala2         ! LN case
        ! during iterations print
        Do iw=lout,lfile
           If(Max(Abs(drhoi(1)),Abs(drhoi(2))).Gt.1.0e-10_pr) Then
              Write(*,*) '  WARNING! Int(Dro)=',Max(Abs(drhoi(1)),Abs(drhoi(2)))
           Endif
           Write(iw,2) iiter,bbroyden,si,xmix,bet,etot,varmas,rms(1),rms(2),ept(1),delLN(1), &
                ept(2),delLN(2),alast(1),alast(2),time 
        End Do
     End If
     !-------------------------------------------------
     ! HFBDIAG
     !-------------------------------------------------        
     Do it=itmin,itmax
        Call hfbdiag(it,0)   ! hfb diagonalization with minimal canonical
        If(ierror_flag.Ne.0) Return
     End Do
     !-------------------------------------------------
     ! EXPECT, DENSIT, COULOMB, FILED, GAMDEL
     !-------------------------------------------------        
     Call expect(.False.)    ! expectation values
     If(ierror_flag.Ne.0) Return	 
     Call field              ! new fields
     If(ierror_flag.Ne.0) Return	 
     Call gamdel             ! hf-matrix
     If(ierror_flag.Ne.0) Return	 
     !-------------------------------------------------
     ! Dumping control (old linear mixing)
     !-------------------------------------------------
     xmix0=0.10 !original 0.1
     If(si.Lt.siold) Then
        xmix=Min(xmax,xmix * 1.130_pr);  !old value 1.13
     Else
        xmix=xmix0
     End If
     siold=si
     !-------------------------------------------------
     ! time per iteration
     !-------------------------------------------------
     Call Cpu_time(time2)
     time=time2-time1; time3=time3+time
     !-------------------------------------------------
     ! Solution is OK within the iteration limit
     !-------------------------------------------------
     If(iiter.Ge.2.And.si.Lt.epsi) Then
        If(iLST1.Eq.0) Then
           iError_in_HO=0 
        Else 
           iError_in_THO=0
        End If
        ! iteration interrupted print
        If(.Not.lpr) Then
           delLN=del; If(kindhfb.Lt.0) delLN=del+ala2 
           Do iw=lout,lfile       
              Write(iw,3) iiter,bbroyden,si,xmix,bet,etot,varmas,rms(1),rms(2),ept(1),delLN(1), &
                   ept(2),delLN(2),alast(1),alast(2),time !Max(Abs(drhoi(1)),Abs(drhoi(2)))
              Write(iw,'(a,f8.3,a)') '  Total CPU time=',time3/60.0_pr,' minutes'
           Enddo
        End If
        ! converged print
        Do iw=lout,lfile       
           Write(iw,4) iiter,si,iError_in_HO,iError_in_THO
           Write(iw,'(a,f8.3,a)') '  Total CPU time=',time3/60.0_pr,' minutes'
        Enddo
        iiter=iiter+1
        Return       
     End If
     !-------------------------------------------------
     ! Slow convergence and lambda >0 (stop iterations)
     !-------------------------------------------------
     If(iiter.Ge.1000.And.(alast(1).Gt.zero.Or.alast(2).Gt.zero)) Exit
     !
  End Do    ! ite
  iiter=iiter+1
  !-------------------------------------------------
  ! Solution interrupted due to iterations limit
  !-------------------------------------------------
  If(iLST1.Eq.0) Then 
     iError_in_HO=-1
  Else
     iError_in_THO=-1
  End If
  delLN=del; If(kindhfb.Lt.0) delLN=del+ala2 
  ! iterations limit print
  Do iw=lout,lfile       
     Write(iw,2) iiter,bbroyden,si,xmix,bet,etot,varmas,rms(1),rms(2),ept(1),delLN(1), &
          ept(2),delLN(2),alast(1),alast(2),Max(Abs(drhoi(1)),Abs(drhoi(2)))
     Write(iw,5) iiter,si,iError_in_HO,iError_in_THO
     Write(iw,'(a,f8.3,a)') '  Total CPU time=',time3/60.0_pr,' minutes'
  Enddo
  !-------------------------------------------------
2 Format(i4,a,1x,f12.8,f5.2,f7.3,f11.3,1x,f6.1,2(f8.3),' | ',4(f8.3),' | ',20(f8.3))    
3 Format(2x,130('-'),/,'  *   iteration interrupted after',i4,' steps   si=',f17.10,' ho=',i3,' tho=',i3,/,2x,130('-'))    
4 Format(2x,130('-'),/,'  *   iteration converged   after',i4,' steps   si=',f17.10,' ho=',i3,' tho=',i3,/,2x,130('-'))
5 Format(2x,130('-'),/,'  *   iterations limit interrupt after',i4,' steps   si=',f17.10,' ho=',i3,' tho=',i3,/,2x,130('-'))    
  !-------------------------------------------------
End Subroutine iter
!====================================================================
! 
!====================================================================
Subroutine hfbdiag(it,icanon)
  !------------------------------------------------------------------
  ! Skyrme-HFB diagonalization in axial HO/THO basis
  !------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Logical :: lpr_pwi,norm_to_improve    
  Real(pr) :: al,al2,emin,hla,dla,pn,eqpe,ela,enb,enb1,ekb,  &
       s1,s2,s3,alnorm,sitest
  Integer(ipr) :: iw,it,i0,icanon,ibiblo,ier,i,j,k,k0,kl,lc,ib,nd,  &
       nhfb,n1,n2,kaib,m,ndk,nd1,nd2,kdib,k1,k2,id1,id2,  &
       n12,n21,ntz,nhhph,nhhpp,ibro,ibroib,i_uv,i_eqp
  Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
  Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
  !
  Call get_CPU_time('hfbdiag',0)
  !	 
  If(it.Eq.1) Then
     EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN
  Else
     EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP
  Endif
  KpwiPo=0; KqpPo=0  
  !
  nhhph=(it-1)*nhhdim; nhhpp=(it+1)*nhhdim 
  !
  !------------------------------------------------------------------
  ! Loop the internal normalization
  !------------------------------------------------------------------
  sitest=Max(Min(0.10_pr,si*0.010_pr),0.000010_pr)
  norm_to_improve=.True.; inner(it)=-1; sumnz(it)=one
  Do While(norm_to_improve)
     !
     inner(it)=inner(it)+1
     !
     If(Abs(sumnz(it)).Lt.sitest.Or.inner(it).Eq.20) norm_to_improve=.False.
     !
     sumnz(it)=zero; v2min(it)=one; Dispersion(it)=zero
     !
     kl=0; emin=1000.0_pr; al=ala(it)
     !
     ! blocking
     If(iparenti(it).Eq.0) blomax(it)=0
     blo123d(it)=0; blok1k2d(it)=0; blocanon(it)=0; 
     ibiblo=bloblo(keyblo(it),it)
     !------------------------------------------------------------------
     ! Runs over blocks
     !------------------------------------------------------------------
     i_uv=0; i_eqp=0
     lc=0; lcanon(0,it)=0; klmax=0; ibro=0
     Do ib=1,nb
        nd=id(ib); nhfb=nd+nd; i0=ia(ib); m=ib+(it-1)*nbx; ibroib=ibro
        !------------------------------------------------------------------
        !  hfb-matrix
        !------------------------------------------------------------------
        Do n1=1,nd
           nd1=n1+nd
           Do n2=1,n1
              nd2=n2+nd; ibro=ibro+1
              hla=brin(nhhph+ibro);     dla=brin(nhhpp+ibro)     
              hfb(n1,n2)=hla;           hfb(nd2,n1)=dla 
              hfb(nd1,n2)=dla;          hfb(nd1,nd2)=-hla
           End Do
           hfb(n1,n1)=hfb(n1,n1)-al; hfb(nd1,nd1)=hfb(nd1,nd1)+al
        End Do
        ier=0; 
        Call DSYEVD('V','L',nhfb,hfb,ndx2,evvk,ALWORK,ialwork,LWORK,ilwork,ier)  
        !Call dsyev('V','L',nhfb,hfb,ndx2,evvk,ALWORK,ialwork,ier)        
        !------------------------------------------------------------------
        ! NB! Diagonalization bug in LAPACK 
        !------------------------------------------------------------------
        If(ier.Gt.0) Then
           Do iw=lout,lfile
              Write(iw,*) 'FATAL ERROR CONDITION IN HFBDIAG DSYEVD, ier=',ier,'(RECOVERED)' 
           End Do
           ibro=ibroib
           Do n1=1,nd
              nd1=n1+nd
              Do n2=1,n1
                 nd2=n2+nd; ibro=ibro+1
                 hla=brin(nhhph+ibro);     dla=brin(nhhpp+ibro)     
                 hfb(n1,n2)=hla;           hfb(nd2,n1)=dla 
                 hfb(nd1,n2)=dla;          hfb(nd1,nd2)=-hla
              End Do
              hfb(n1,n1)=hfb(n1,n1)-al; hfb(nd1,nd1)=hfb(nd1,nd1)+al
           End Do
           Call sdiag(ndx2,nhfb,hfb,evvk,hfb,zhfb,+1)
        End If
        !------------------------------------------------------------------
        ! Blocking
        !------------------------------------------------------------------
        ! external blocking
        If(iiter.Eq.1.And.inner(it).Eq.0) Then
           If(iparenti(it).Ne.0.And.keyblo(it).Eq.0) Then
              ! eventually charging 
              !   keyblo(it)=1
              !   bloblo(keyblo(it),it)=ib
              !   blo123(keyblo(it),it)=requested level (k0)
              Call requested_blocked_level(ib,it)
              If(ierror_flag.Ne.0) Return	 		  
  	      ibiblo=bloblo(keyblo(it),it)		  
           Endif
        Endif
        ! general blocking
        k0=0
        If(ibiblo.Eq.ib) Then
           If(iiter.Eq.1.And.inner(it).Eq.0) Then 
              ! blocked level as in the even-even nucleus           
              k0=blo123(keyblo(it),it); ndk=k0+nd
              Do n2=1,nd
                 nd2=n2+nd
                 hfb1(n2,it)=hfb(n2,ndk)    !U
                 hfb1(nd2,it)=hfb(nd2,ndk)  !V
              Enddo
              ! number of states in the block to be tested
              blocross(it)=Min(blomax(it)+10,nd)
           Endif
           ! overlap between new and old blocked levels
           s3=zero
           Do n1=1,blocross(it)
              ndk=n1+nd; s1=zero  
              Do n2=1,nd 
                 nd2=n2+nd
                 s1=s1+Abs(hfb1(nd2,it)*hfb(nd2,ndk)) !VV
                 s1=s1+Abs(hfb1(n2,it)*hfb(n2,ndk))   !UU 
              Enddo
              If(s1.Gt.s3) Then
                 s3=s1; k0=n1
              Endif
           Enddo
           blo123d(it)=k0 
           If(.Not.norm_to_improve) Then
              ! find maximal HO component
              ndk=k0+nd
              s1=zero
              Do n1=1,nd
                 nd1=n1+nd
                 hfb1(n1,it)=hfb(n1,ndk); hfb1(nd1,it)=hfb(nd1,ndk)
                 s2=Max(s1,Abs(hfb(n1,ndk)),Abs(hfb(nd1,ndk)))
                 If(s2.Gt.s1) Then
                    s1=s2; i=n1+i0  ! labels in k[k1,k2] numbering
                 End If
              End Do
              ! print blocked state
              Do iw=lout,lfile
                 Write(iw,'(4x,a,2(a,i3),2x,3(a,1x,f12.8,1x),(i3,a,i3,1x),a)')  &
                      protn(it),' Blocking: block=',ib,  &
                      ' state=',k0,  &
                      ' Eqp=',evvk(k0+nd),  &
                      ' Dqpe=',evvk(k0+nd)-eqpmin(it),  &
                      ' Ovlp=',s3  &                                  
                      , keyblo(it),'/',blomax(it)  &
                      , tb(i)
              Enddo
              ! ieresbl=6, 'BLKN','BLKZ'
              ereslbl(it)=tb(i)
              If(it.Eq.1) Then
                 ! 'BlEqpN','BlDEqpN','BlOvrN'
                 eresbl(1)=evvk(k0+nd); eresbl(2)=evvk(k0+nd)-eqpmin(it); eresbl(3)=s1
              Else
                 ! 'BlEqpZ','BlDEqpZ','BlOvrZ'
                 eresbl(4)=evvk(k0+nd); eresbl(5)=evvk(k0+nd)-eqpmin(it); eresbl(6)=s1
              Endif
           End If
        End If
        !------------------------------------------------------------------
        ! Run over all qp states k in the block
        !------------------------------------------------------------------
        kaib=kl
        Do k=1,nd
           ndk=k+nd
           ! referent spectra
           pn=zero
           Do i=1,nd
              hla=hfb(i+nd,ndk)**2; pn=pn+hla
              If(i.Eq.1) Then                         
                 !write(*,*) hfb(i+nd,ndk),hfb(i,ndk),evvk(nd+k)         ! Vak, Uak, Sign( V*U)  Ek 
                 !write(*,*) hfb(i,nd-k+1),hfb(i+nd,nd-k+1),evvk(nd-k+1) ! Vak, Uak, Sign(-V*U) -Ek
              Endif
           Enddo
           ! Blocking 
           If(k.Eq.k0) Then 
              n1=k0+nd
              Do i=1,nd 
                 hla=hfb(i+nd,n1)**2; dla=hfb(i,n1)**2; pn=pn-half*(hla-dla)
              Enddo
           Endif
           eqpe=evvk(nd+k); ela=eqpe*(one-two*pn) 
           enb=ela+al;      ekb=Sqrt(Abs(eqpe*eqpe-ela**2))
           !------------------------------------------------------------------
           ! cut-off condition: energy pwi + Fermi cut-off function
           !------------------------------------------------------------------
           lpr_pwi=enb.Le.pwi.Or.Abs(one/(one+Exp(100.0_pr*(enb-pwi)))).Gt.1.0d-6              
#ifndef hide_qrpa
           lpr_pwi=enb.Le.pwi !cristina sharp cut off for qrpa
#endif
           If(use_TMR_pairing.Ne.0) lpr_pwi=.True. !no pairing window with TMR pairing
           !------------------------------------------------------------------
           ! Remember the whole qp solution
           !------------------------------------------------------------------
           If(.Not.norm_to_improve) Then
              i_eqp=i_eqp+1
              EqpPo(i_eqp)=evvk(nd+k)                       ! Eqp_k
              If(lpr_pwi) KqpPo(kl+1)=i_eqp                 ! below pwi otherwise zero
              If(lpr_pwi) KpwiPo(kl+1)=i_uv                 ! below pwi otherwise zero
              Do n2=1,nd
                 nd2=n2+nd; i_uv=i_uv+1
                 UqpPo(i_uv)=hfb(n2,ndk)                    ! U_ak
                 VqpPo(i_uv)=hfb(nd2,ndk)                   ! V_ak
              End Do
           Endif
           !------------------------------------------------------------------
           ! Pairing window
           !------------------------------------------------------------------
           If(lpr_pwi) Then  
              kl=kl+1                                      !number of active states
              If(k0.Eq.k) blok1k2d(it)=kl                  !blocking: dynamic #: k[k1,k2] numbering
              If((eqpe.Le.emin).And.(pn.Gt.0.0001)) Then   !to avoid unocc at magic numbers
                 emin=eqpe; alnorm=pn                      !min qpe and its occupation
              End If
              erhfb(kl)=enb; drhfb(kl)=ekb; uk(kl,it)=pn   !ref.s.p. energies, deltas, occupancies
              sumnz(it)=sumnz(it)+two*pn                   !internal normalization
           End If
        End Do
        If(norm_to_improve) Cycle
        !------------------------------------------------------------------
        !  Density matrices
        !------------------------------------------------------------------
        kdib=kl-kaib; ka(ib,it)=kaib; kd(ib,it)=kdib
        k1=kaib+1; k2=kaib+kdib
        eqpe=0.
        Do n2=1,nd
           Do n1=n2,nd
              s1=zero; s2=zero
              If(k1.Le.k2) Then
                 Do k=k1,k2
                    nd1=KpwiPo(k)+n1; nd2=KpwiPo(k)+n2
                    s1=s1+VqpPo(nd1)*VqpPo(nd2) 
                    s2=s2+UqpPo(nd1)*VqpPo(nd2) +VqpPo(nd1)*UqpPo(nd2)
                 End Do
			     !if(it.eq.1) write(200,*) n1,n2,s1-s2
                 s1=two*s1; s2=half*s2               ! two:due to m-projection, half:due to symmetrization
                 ! blocking
                 If(ibiblo.Eq.ib) Then
                    i=blok1k2d(it); id1=KpwiPo(i)+n1; id2=KpwiPo(i)+n2					
                    s1=s1-VqpPo(id1)*VqpPo(id2)+UqpPo(id1)*UqpPo(id2) 
                    s2=s2-half*(UqpPo(id1)*VqpPo(id2)+VqpPo(id1)*UqpPo(id2))
                 Endif
              End If
              n12=n1+(n2-1)*nd; n21=n2+(n1-1)*nd
              rk(n12,m)=s1; rk(n21,m)=s1              !  V V'
              ak(n12,m)=-s2; ak(n21,m)=-s2            !- U V', ak=half*(pairing density)             
              hfbcan(n1,n2)=s1; hfb(n1,n2)=s1
           End Do !n1
        End Do !n2
        !------------------------------------------------------------------
        ! Canonical basis 
        !------------------------------------------------------------------
        If(k1.Le.k2) Then
           Call Canonical(it,icanon,k2,k1,nd,i0,lc,ib,ibiblo,m,ibroib)
           If(ierror_flag.Ne.0) Return	 
        Endif
        lcanon(ib,it)=lc
     End Do !ib	 
     If(kl.Eq.0) Then
        ierror_flag=ierror_flag+1
        ierror_info(ierror_flag)=' STOP: kl=zero, no states below pwi!!!'       
        Return
     Endif
     If(iparenti(it).Ne.0.And.ibiblo.Eq.0) Then
        ierror_flag=ierror_flag+1
        ierror_info(ierror_flag)='STOP: No blocking candidate found!!!'	 
        Return
     Endif
     eqpmin(it)=emin; klmax(it)=kl; sumnz(it)=sumnz(it)-tz(it) 
     !------------------------------------------------------------------
     ! Lambda search
     !------------------------------------------------------------------
     Call ALambda(al,it,kl)
     If(ierror_flag.Ne.0) Return	 	 
     If(keyblo(it).Eq.0) Then
        ala(it)=al 
     Else
        ala(it)=ala(it)+0.50_pr*(al-ala(it))   
     End If     
     ! NB! 'alast' instead of 'al' at small pairing
     alast(it)=al
     If(Abs(ept(it)).Lt.0.0001) Then
        ntz=tz(it)+0.1; ntz=ntz/2
        Do k=1,kl
           drhfb(k)=erhfb(k)
        End Do
        Call ord(kl,drhfb)
        alast(it)=drhfb(ntz)  !last bound s.p. energy
     End If
     !------------------------------------------------------------------
     ! THO asymptotic decay
     !------------------------------------------------------------------
     ! density asymptotic decay \rho(r)->Exp(-ass(it)*r)
     ! ass(it)=2*Sqrt((E_min-\lambda)/((A-1)/A)*hbar**2/(2*m)))
     al2=zero
     If(kindhfb.Lt.0) Then
        al2=al+two*ala2(it)*(one-two*alnorm) ! al=al+two*ala2(it)
     Endif
     al2=(emin-al2)/hb0
     ! wrong asymptotic
     iasswrong(it)=0; If(al2.Le.zero) iasswrong(it)=1; ass(it)=two*Sqrt(Abs(al2))
     !
  End Do ! While(norm_to_improve)
  !
  If(it.Eq.1) Then
     Call get_CPU_time('hfbdiag(n)',1)
  Else
     Call get_CPU_time('hfbdiag(p)',1)
  Endif	 
  !
End Subroutine hfbdiag
!=======================================================================
!
!=======================================================================
Subroutine ALambda(al,it,kl)
  !---------------------------------------------------------------------
  ! Adjusting Fermi energy
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Real(pr), Save :: fm7=1.0d-7,fm10=1.0d-10
  Real(pr)  :: al,vh ,xinf,xsup,esup,ez,dez,dvh,y,a,b,einf,absez,sn
  Integer(ipr) :: it,i,k,kl,icze,lit,ntz
  !-------------------------------------------------
  ! Chemical potential without pairing
  !-------------------------------------------------
  If(CpV0(it-1).Eq.zero) Then
     ntz=tz(it)+0.1; ntz=ntz/2
     Do k=1,kl
        drhfb(k)=erhfb(k)
     End Do
     Call ord(kl,drhfb)
     If (ntz.Lt.kl) Then
        al=half*(drhfb(ntz)+drhfb(ntz+1))
     Else
        al=drhfb(ntz)+0.001
     End If
     Return
  End If
  !-------------------------------------------------
  ! Chemical potential with pairing
  !-------------------------------------------------
  xinf=-100.; xsup=80.; esup=one; icze=0
  Do lit=1,500
     sn=zero; dez=zero
     Do i=1,kl
        vh=zero; dvh=zero 
        y=erhfb(i)-al; a=y*y+drhfb(i)**2; b=Sqrt(a)
        If(b.Gt.zero)  vh=half*(one-y/b)
        If(b.Lt.fm7.And.icze.Eq.1) vh=-einf/(esup-einf) !no pairing
        If(vh.Lt.zero) vh=zero
        If(vh.Gt.one)  vh=one
        If(b.Gt.zero) dvh=drhfb(i)**2/(a*b)         ! D[ez,al](i)
        ! blocking
        If(i.Eq.blok1k2d(it)) Then
           vh=half; dvh=zero
        End If
        sn=sn+two*vh; dez=dez+dvh                   ! D[ez,al]
     End Do
     ez=sn-tz(it); absez=Abs(ez)
     !-------------------------------------------------
     ! Correcting bounds
     !-------------------------------------------------
     If(ez.Lt.zero) Then
        xinf=Max(xinf,al); einf=ez
     Else
        xsup=Min(xsup,al); esup=ez
     End If
     If(lit.Eq.1) Then
        If(absez.Le.0.10_pr) Then
           al=al-ez
        Else
           al=al-0.10_pr*Sign(one,ez)
        End If
     Else
        al=al-ez/(dez+1.d-20)                      ! newton method
     End If
     If(xsup-xinf.Lt.fm7) icze=1                      ! low/upp close
     If(al.Lt.xinf.Or.al.Gt.xsup) al=half*(xinf+xsup) ! mean upp/low
     If(absez.Le.fm10) Return
  End Do
  !-------------------------------------------------
  ! Low accuracy warning
  !-------------------------------------------------
  Write(lout,'(a,2(e12.5,2x),a,2(2x,f8.4),a,i2)') ' Low accuracy=',sn,ez,' for N,Z=',tz,' it=',it
End Subroutine Alambda
!=========================================================================================
!
!=========================================================================================
Subroutine Canonical(it,icanon,k2,k1,nd,i0,lc,ib,ibiblo,m,ibroib)
  !---------------------------------------------------------------------------------------
  ! Canonical diagonalization
  !---------------------------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: it,i0,icanon,ibiblo,i,iw,k,kk,lc,ib,nd,n1,n2,m,nd1,k1,k2,n12,ier
  Integer(ipr) :: nhhph,nhhpp,ibro,ibroib
  Real(pr) :: s1,s2,vx,h1,d1,h2,d2,ddn1,ddn2
  Real(pr), Allocatable :: hh(:,:),de(:,:)
  Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
  Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
  !
  If(it.Eq.1) Then
     EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN
  Else
     EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP
  Endif
  !
  If(Abs(ept(it)).Lt.0.0001) Then
     !------------------------------------------------------
     ! No pairing => just taking the HF states
     !------------------------------------------------------
     Do k=1,nd 
        kk=k1+k-1; lc=lc+1                              ! total number of the canonical states
        ddc(1:nd,lc,it)=zero; vk(lc,it)=zero            ! zeros: nd could be larger then k2-k1+1 
        If(kk.Gt.k2) Cycle
        vx=zero
        Do i=1,nd
           h1=VqpPo(KpwiPo(kk)+i)**2; vx=vx+h1
        Enddo
        If (vx.Le.zero) vx=zero                         ! roundoff errors
        If (vx.Ge.one ) vx=one
        Do i=1,nd
           If(vx.Ge.half) Then
              ddc(i,lc,it)=VqpPo(KpwiPo(kk)+i)          ! (ph) s.p. orbitals in conf.space
           Else
              ddc(i,lc,it)=UqpPo(KpwiPo(kk)+i)          ! (ph) s.p. orbitals in conf.space
           Endif
        Enddo
        Dispersion(it)=Dispersion(it)+four*vx*(one-vx)  ! internal P/N Dispersion
        If(Abs(vx-half).Le.v2min(it)) Then
           v2min(it)=Abs(vx-half); v2minv(it)=vx        ! divergent condition             
           lcc=lc
        Endif
        vk(lc,it)=vx                                    ! (ph) s.p. occupations v^2
        !------------------------------------------------------
        ! RESU only
        !------------------------------------------------------
        If(icanon.Ne.0) Then
           ek(lc,it)=EqpPo(KqpPo(kk))*(one-two*vx)+ala(it)      ! (ph) s.p. energies
           dk(lc,it)=zero                                       ! (ph) s.p. deltas
        Endif
     Enddo !k
  Else
     !------------------------------------------------------
     ! Pairing => calculate canonical basis
     !------------------------------------------------------
     ier=0; Call DSYEVD('V','L',nd,hfbcan,ndx,evvkcan,ALWORK,ialwork,LWORK,ilwork,ier)
     ! bug in LAPAK
     If(ier.Gt.0) Then 
        Do iw=lout,lfile
           Write(iw,*) 'FATAL ERROR CONDITION IN CANONICAL DSYEVD, ier=',ier,'(RECOVERED)'
        End Do
        Do n2=1,nd
           Do n1=n2,nd
              vx=hfb(n1,n2)
              hfbcan(n2,n1)=vx; hfbcan(n1,n2)=vx
           End Do
        End Do
        Call sdiag(ndx,nd,hfbcan,evvkcan,hfbcan,zhfb,+1)
     End If
     !------------------------------------------------------
     ! Eigenvalues and wavefunctions
     !------------------------------------------------------
     Do k=1,nd
        lc=lc+1                                      ! total number of the canonical states
        Do i=1,nd
           ddc(i,lc,it)=hfbcan(i,k)                    ! (ph) canon orbitals in conf.space
        End Do
        vx=evvkcan(k)*half
        If (vx.Le.zero) vx=zero                        ! roundoff errors
        If (vx.Ge.one ) vx=one
        ! blocking 
        If(ibiblo.Eq.ib.And.vx.Gt.0.49.And.vx.Le.0.51) blocanon(it)=lc
        Dispersion(it)=Dispersion(it)+four*vx*(one-vx) ! internal P/N Dispersion
        If(Abs(vx-half).Le.v2min(it)) Then
           v2min(it)=Abs(vx-half); v2minv(it)=vx       ! divergent condition             
           lcc=lc
        Endif
        vk(lc,it)=vx                                   ! (ph) canon occupations v^2
        !------------------------------------------------------
        ! RESU only
        !------------------------------------------------------
        If(icanon.Ne.0) Then
           ! canonical energies and deltas (no physical meaning in PNP)
           nhhph=(it-1)*nhhdim; nhhpp=(it+1)*nhhdim
           Allocate(hh(nd,nd),de(nd,nd))           
           ibro=ibroib  
           Do n1=1,nd
              Do n2=1,n1
                 ibro=ibro+1
                 vx=brin(nhhph+ibro); hh(n2,n1)=vx; hh(n1,n2)=vx   
                 vx=brin(nhhpp+ibro); de(n2,n1)=vx; de(n1,n2)=vx   
              End Do
           End Do
           h1=zero; d1=zero
           Do n2=1,nd
              h2=zero; d2=zero
              Do n1=1,nd
                 ddn1=hfbcan(n1,k)
                 h2=h2+ddn1*hh(n1,n2)
                 d2=d2+ddn1*de(n1,n2)
              End Do
              ddn2=hfbcan(n2,k)
              h1=h1+h2*ddn2
              d1=d1+d2*ddn2
           End Do
           ek(lc,it)=h1                              ! (ph) canon s.p. energies
           dk(lc,it)=d1                              ! (ph) canon s.p. deltas
           Deallocate(hh,de)
        End If
        !
     End Do !k
  Endif
  !------------------------------------------------------
  ! RESU only
  !------------------------------------------------------
  If(icanon.Ne.0) Then
     !------------------------------------------------------
     ! Find maximal HO components of all qp states
     !------------------------------------------------------
     Do k=k1,k2
        s1=zero
        Do n1=1,nd
           nd1=nd+n1
           s2=Max(s1,Abs(VqpPo(KpwiPo(k)+n1)),Abs(UqpPo(KpwiPo(k)+n1)))
           If(s2.Gt.s1) Then
              s1=s2
              vkmax(k,it)=s1                       ! maximal overlap
              numax(k,it)=n1+i0                    ! its number in k[k1,k2] numbering
           End If
        End Do
     End Do
  End If
  !------------------------------------------------------
  ! Searching for possible blocking candidates
  !------------------------------------------------------
  If(icanon.Ne.0) Then
     ! Search within |(1-2*N)*Eqpe| lover than 'pwiblo' 
     ! The levels number n1 is 1,2,3,... for the given block ([123] numbering)
     If(iparenti(it).Eq.0) Then
        n1=0
        Do k=k1,k2
           n1=n1+1
           If(Abs(EqpPo(KqpPo(k))-eqpmin(it)).Le.pwiblo) Then
              blomax(it)=blomax(it)+1                         ! blocked state #, maximel # of block candidates
              If(blomax(it).Gt.bloall) Then
                 ierror_flag=ierror_flag+1
                 ierror_info(ierror_flag)='Too many blocking candidates! Increase bloall and run again'
                 Return
              Endif
              bloblo(blomax(it),it)=ib                        ! block where to block
              blo123(blomax(it),it)=n1                        ! state # [123] numbering
              blok1k2(blomax(it),it)=k                        ! state # k[k1,k2] numbering
              bloqpdif(blomax(it),it)=Abs(EqpPo(KqpPo(k))-eqpmin(it))
           End If
        End Do
     End If
  End If
  !
End  Subroutine canonical
!=======================================================================
!
!=======================================================================
Subroutine resu(irecord)
  !---------------------------------------------------------------------
  ! prints results: single particle energies, densities, fields
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: it,iw,ib,im,m,nd,k,k0,k1,k2,j,n,imax,nhfb,irecord
  Real(pr)     :: sum,eqpe,pn,ela,enb,ek0,vk0,ekk,delb,ovmax,s,uuvv,  &
       dk0,skk,summ(4),vvs,vvc,enjacek
  Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
  Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
  !
  !--------------------------------------------
  ! last HFB run for full canon.calculations
  !--------------------------------------------
  Do it=itmin,itmax
     Call hfbdiag(it,1)       ! hfb with maximal canonical
     If(ierror_flag.Ne.0) Return	 	 
  End Do
  !PAV
  Call expect(.False.)        ! expectation values
  If(ierror_flag.Ne.0) Return	   
  Call field                  ! new fields
  If(ierror_flag.Ne.0) Return	   
  Call gamdel                 ! hf-matrix
  If(ierror_flag.Ne.0) Return	   
  ! inout(2): HFB matrices if nucleus is even-even
  If(npr(1).Eq.2*(npr(1)/2).And.npr(2).Eq.2*(npr(2)/2)) Then
     Call inout(2)               
     If(ierror_flag.Ne.0) Return	   
  Endif
  !--------------------------------------------
  ! Printing densities and fields
  !--------------------------------------------
  ! Call printLST(ro(1,1),ro(1,2))      !need fix
  ! call printLST(tau(1,1),tau(1,2))
  ! call printLST(dro(1,1),dro(1,2))
  ! call printLST(dj (1,1),dj (1,2))
  ! printing of fields
  ! call printLST(vhb(1,1),vhb(1,2))
  ! call printLST(v (1,1),v (1,2))
  ! call printLST(vs(1,1),vs(1,2))
  !--------------------------------------------
  ! Printing quasiparticle states
  !--------------------------------------------
  Do it=itmin,itmax
     If(it.Eq.1) Then
        EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN
     Else
        EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP
     Endif
     !
     If(Print_Screen) Then
        iw=lfile       
        Write(iw,200) tit(it)
        Write(iw,*) ' eqp(k) -> q.p. energy '
        Write(iw,*) ' e(k)   -> referent s.p. energy '
        Write(iw,*) ' p(k)   -> occ.probability '
        Write(iw,*) ' del(k) -> referent s.p. gap '
        Write(iw,*) ' fermi energy alast=',alast(it)
        Write(iw,'(a,a)')  &
             '  #k  block#    eqp(k)     e(k)       (1-2N)E      decay        p(k)',  &
             '        del(k)    overl      labels'
200     Format(//,' #quasiparticle energies ',a,/,1x,32('-'))
     Endif
     sum=zero
     Do ib=1,nb
        nd=id(ib); im=ia(ib); m=ib+(it-1)*nbx; nhfb=nd+nd
        k1=ka(ib,it)+1
        k2=ka(ib,it)+kd(ib,it)
        If(k1.Le.k2) Then
           Do k=k1,k2                                 ! print active states only
              pn=uk(k,it)                             ! qp probabilities
              j=k
              If(pn.Gt.-1.d-14) Then                   ! print If signIficant pn
                 ! main oscillator component
                 ovmax=vkmax(k,it)                    ! maximal overlap
                 imax=numax(k,it)                     ! its number
                 ! printing
                 eqpe=EqpPo(KqpPo(k))                 ! qp energies
                 skk=two*Sqrt(Abs(eqpe-ala(it))/hb0)  ! qp decay
                 ela=eqpe*(one-two*pn)
                 enb=ela+ala(it)                      ! ref. s.p. energies
                 delb=Sqrt(Abs(eqpe**2-ela**2))       ! ref. s.p. delta
                 sum=sum+two*pn                       ! particle number
                 If(Print_Screen) Then
                    iw=lfile       
                    Write(iw,201) k,ib,eqpe,enb,(one-two*pn)*eqpe,skk,pn,delb,ovmax,tb(imax)
201                 Format(i4,2x,i3,1x,f12.6,f12.6,f12.6,f12.6,2x,f12.8,  &
                         2(2x,f7.4),' ',a13)
                 Endif
              End If
           End Do
        End If
     End Do !ib
     !--------------------------------------------
     ! Printing canonical single particle states
     !--------------------------------------------
     If(Print_Screen) Then        
        iw=lfile                 
        Write(iw,'(a,i4,a,i4)')  &
             '#all active are ',j,' q.p. states out of ',nt
        Write(iw,'(a,f6.1)') '#since the cut off is pwi=',pwi
        Write(iw,'(3a,f6.1)')'#check: number of ',tit(it),'=',sum
        Write(iw,100) tit(it)
        Write(iw,*) ' labels -> {2*omega}{parity}[nn=nz+2*nr+nl,nz,nl]'
        Write(iw,*) ' cqpe   -> canonical q.p. energies'
        Write(iw,*) ' ce     -> canonical s.p. energies'
        Write(iw,*) ' fermi energy=',alast(it)
        Write(iw,*) ' average cdelt=',del(it)
        Write(iw,'(a,a)')'  k0      ceqp        ce         v*v',  &
             '       u*v        cdel     overl      labels'
100     Format(//,' #canonical s.p. energies ',a,/,1x,33('-'),//)
     Endif
     k0=0
     summ=zero; enjacek=zero
     Do ib=1,nb
        nd=id(ib); im=ia(ib)
        k1=ka(ib,it)+1; k2=ka(ib,it)+kd(ib,it)
        If(k1.Le.k2) Then
           Do k=1,nd
              k0=k0+1
              ! for Lipkin Nogami
              vvs=two*Sqrt(vk(k0,it))*Sqrt(one-vk(k0,it))    !2vu
              vvc=two*vk(k0,it)-one                          !2v^2-1
              summ(1)=summ(1)+vvs**2
              summ(2)=summ(2)+vvs**2*vvc
              summ(3)=summ(3)+vvs**4
              summ(4)=summ(4)+(vvs*vvc)**2
              ! search for main oscillator component
              ovmax=zero
              Do n=1,nd
                 s=Abs(ddc(n,k0,it))                         !canon orbitals in conf.space
                 If (s.Ge.ovmax) Then
                    ovmax=s; imax=n
                 End If
              End Do
              ! printing
              ek0=ek(k0,it)                                  !canon s.p. energies
              enjacek=enjacek+ek0*vk(k0,it)
              If(ek0.Lt.pwi) Then                            !print up to 'pwi'
                 vk0=vk(k0,it)                               !canon occupations v^2
                 If(vk0.Gt.-1.d-4) Then                       !print If signIficant v^2
                    dk0=-dk(k0,it)                           !canon s.p. deltas
                    ekk=Sqrt((ek0-ala(it))**2+dk(k0,it)**2)  !resulting cqpe
                    uuvv=Sqrt(Abs(vk0*(one-vk0)))            !resulting u*v
                    If(Print_Screen) Then
                       iw=lfile       
                       Write(iw,101) k0,ekk+ala(it),ek0,vk0,uuvv,dk0,ovmax,tb(im+imax)
101                    Format(i4,2f12.6,2(1x,f12.8),2(2x,f7.4),' ',a13)
                    Endif
                 End If
              End If
           End Do !k0
        End If
     End Do !ib
     !--------------------------------------------
     ! Lipkin-Nogami
     !--------------------------------------------
     ssln(1,it)=summ(1)
     ssln(2,it)=summ(2)
     ssln(3,it)=summ(4)*summ(1)-summ(2)**2+summ(1)**3/4.0_pr-half*summ(3)*summ(1)
     If(Print_Screen) Then
        iw=lfile       
        Write(iw,*) ' Sum canonical e_v*V^2_k=',two*enjacek  
     Endif
  End Do !it
  !--------------------------------------------
  ! To thoout.dat, thodef.dat and hodef.dat 
  !--------------------------------------------
  If(irecord.Ne.0) Then
     iappend=1
     Call expect(.True.)    !print  & record HFB+PAV results
     If(ierror_flag.Ne.0) Return	   	 
     iappend=0
  Else
     Call expect(.True.)    !print HFB+PAV results
     If(ierror_flag.Ne.0) Return	   	 
  End If
  !
End Subroutine resu
!=======================================================================
!
!=======================================================================
Subroutine read_HFBTHO_NAMELIST  
  Use HFBTHO
  Implicit None
  Integer(ipr) :: ios,lnamelist=16 
  !------------------------------------
  ! Namelist (default values)
  !------------------------------------   
  iLST_INI=0                      ! 0:HO, -1:HO->THO, 1:THO
  MAX_ITER_INI=100                    ! max number of iterations
  keypj_INI=1                     ! number of gaouge points 
  iproj_INI=0                     ! projecting on different nucleus 
  npr1pj_INI=0                    !  its neutron number
  npr2pj_INI=0                    !  its proton number 
  !epsi_INI=0.00000000001_pr      ! stop criteria: original
  !epsi_INI=0.0000001_pr          ! tests
  epsi_INI=0.00001_pr             ! relaxing tests
  !epsi_INI=0.0001_pr             ! enough mario
  !epsi_INI=0.01_pr               ! fit
  Add_Pairing_INI=.False.         ! add pairing starting from file
  icou_INI=2                      ! coul: no-(0), dir.only-(1), plus exchange-(2)
  DO_FITT_INI=.False.             ! calculates quantities for reg.optimization
  IDEBUG_INI=0                    ! debug
  Parity_INI=.True.               ! reflection symmetry  
  Print_HFBTHO_Namelist_INI=.True.  ! Print Namelist 
  !------------------------------------
  ! Namelist (handling)
  !------------------------------------   
  Open(lnamelist,file='UNEDF_NAMELIST.DAT',DELIM='APOSTROPHE') ! 'QUOTE'   
  Read(UNIT=lnamelist,NML=HFBTHO_NAMELIST,iostat=ios)
  If (ios.Ne.0) Then
     ierror_flag=ierror_flag+1
     ierror_info(ierror_flag)='Error in HFBTHO_NAMELIST read'    
     Return
  Endif
  Close(lnamelist)     
  !
  Call iniialize_HFBTHO_SOLVER
  !
End Subroutine read_HFBTHO_NAMELIST
!=======================================================================
!
!=======================================================================
Subroutine iniialize_HFBTHO_SOLVER
  !---------------------------------------------------------------------
  ! default parameters
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  !------------------------------------
  ! tapes
  !------------------------------------   
  lwin=21; lwou=22;  lwel=12; lres=17; lin=3; lout=6   
  !------------------------------------
  ! From Namelist or default values
  !------------------------------------   
  epsi=epsi_INI                                   ! stop criteria
  Add_Pairing=Add_Pairing_INI                     ! add pairing starting from file
  icou=icou_INI                                   ! coul: no-(0), dir.only-(1), plus exchange-(2)
  DO_FITT=DO_FITT_INI                             ! calculates quantities for reg.optimization
  IDEBUG=IDEBUG_INI                               ! debug
  Parity=Parity_INI                               ! reflection symmetry
  Print_HFBTHO_Namelist=Print_HFBTHO_Namelist_INI ! Print Namelist 
  !------------------------------------
  ! output control
  !------------------------------------   
  If(n00_INI.Gt.0) Then    
     Print_Screen=.True. 
     lfile=lout+1                                 ! output to screen & thoout.dat
  Else
     Print_Screen=.False.
     lfile=lout-1                                 ! no output to screen & thoout.dat
  Endif
  !------------------------------------
  ! Pi
  !------------------------------------
  PI=four*Atan(one)
  !------------------------------------
  ! blocking
  !------------------------------------
  bloblo=0; blo123=0; blok1k2=0;  keyblo=0 
  blomax=0; nkblo=0;  iparenti=0; irestart=0
  blocanon=zero;      eqpmin=zero   
  !------------------------------------
  ! buffers
  !------------------------------------
  eres=zero;  eresu=zero;  eresl=zero;  
  eresj=zero; eresbl=zero; ereslbl=' 00[00,00,00]'
  !------------------------------------
  ! def parameters
  !------------------------------------
  ffdef3=Sqrt(five/(four*pi))/two
  ffdef4=Sqrt(117.0_pr)/(four*pi)
  ffdef5=Sqrt(nine/(four*pi))/eight
  ffdef6=Sqrt(five*pi)/three
  ffdef7=Sqrt(pi)/four
  !------------------------------------
  ! former linear mixing
  !------------------------------------
  xmix0=0.30_pr            ! lowest mixing parameter  (redefined later)
  xmix=0.3_pr              ! initial mixing parameter (changes every iteration)
  xmax=1.0_pr              ! mario
  !------------------------------------
  ! misck (redefined later)
  !------------------------------------
  rehfbcan=0.0_pr; depnp=0.0_pr; ala2=0.00_pr             
  ept=-2.0_pr; del=1.0_pr; ala=-7.0_pr 
  ala1(1)=-14.6851; ala1(2)=-3.7522; si=1.0_pr
  iqrpa=0; icacou=0;  icacoupj=0; icahartree=0; iasswrong=0
  iError_in_HO=0;  iError_in_THO=0
  ECMHFB=0.0_pr; ECMPAV=0.0_pr
  !------------------------------------
  ! Saxon-Woods: von koepf und ring, z.phys. (1991)
  !------------------------------------
  v0ws=-71.28;   akv=0.4616; r0v=1.2334; av=0.6150
  vso=11.1175; rso=1.1443; aso=0.6476
  !------------------------------------
  ! fixed text
  !------------------------------------
  tp(1)='+'; tp(2)='-'; tis(1)='n'; tis(2)='p'; 
  tit(1)='neutrons'; tit(2)='protons '
  tl(0)='s'; tl(1)='p'; tl(2)='d'; tl(3)='f'; tl(4)='g'
  tl(5)='h'; tl(6)='i'; tl(7)='j'; tl(8)='k'; tl(9)='l'
  tl(10)='m'; tl(11)='n'; tl(12)='o'; tl(13)='p'; tl(14)='q'
  tl(15)='r'; tl(16)='s'; tl(17)='t'; tl(18)='u'; tl(19)='v'; tl(20)='w'
  !------------------------------------
  ! fixed parity sign
  !------------------------------------
  tpar(1)=+1; tpar(2)=-1;
  !------------------------------------
  ! physical constants
  !------------------------------------
  amu=938.90590_pr; r0=1.20_pr; 
  alphi=137.036020_pr; hqc=197.328910_pr ! hqc=197.3282849_pr; 
  !------------------------------------
  ! e2 for protons (set now in elsewhere)
  !------------------------------------
  !chargee2=hqc/alphi
  !chargee2=1.43997840_pr 
  !-----------------------------------
  ! set the loops over particle types
  !-----------------------------------
  itmin=1 ; itmax = 2;
  If ((npr_INI(1).Lt.1).And.(npr_INI(2).Lt.1)) Then
     ierror_flag=ierror_flag+1
     ierror_info(ierror_flag)='STOP: Nonsense in particle numbers N<0 and Z<0.'    
     Return  
  Endif
  If (npr_INI(1).Eq.0) itmin = 2
  If (npr_INI(2).Eq.0) itmax = 1
  !-----------------------------------
  ! error flag and info
  !-----------------------------------
  ierror_flag=0
  ierror_info(ierror_flag)='No errors in the solver!'  
  !  
  Call set_functional_parameters(skyrme_INI,.False.)
  !
End Subroutine iniialize_HFBTHO_SOLVER
!=======================================================================
!
!=======================================================================
Subroutine base0(lpr)
  !---------------------------------------------------------------------
  ! selects HO basis configurations in cylindrical coordinates
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Logical :: lpr
  Integer(ipr) :: iw,k,nre,nze,la,le,ip,ir,iz,il,is,Iall,ilauf,jlauf,ib,nd
  Integer(ipr) :: NOSCIL
  Real(pr), Allocatable :: e(:)
  Real(pr) :: hbz,hbp,ee
  !
  If(n00.Gt.n00max) Then
     ierror_flag=ierror_flag+1
     ierror_info(ierror_flag)='STOP: too large n00 versus n00max'
     Return  
  Endif
  !-----------------------------------------------
  ! MAXIMUM NUMBER OF THE HO SHELLS (n00,NOSCIL) 
  ! (7,120),(8,165),(9,220),(10,286),(11,364)
  ! (12,455),(14,680),(16,969),(18,1330),(20,1771)
  !-----------------------------------------------
  NOSCIL=(n00+1)*(n00+2)*(n00+3)/6
  !-----------------------------------------------
  ! count all states for n00max
  !-----------------------------------------------
  nze=n00max; nre=n00max/2; 
  Iall=0; 
  Do k=1,n00max+1
     la=k-1; le=min0(n00max,k)
     Do ip=1,2
        Do ir=0,nre
           Do iz=0,nze
              Do il=la,le
                 Do is=+1,-1,-2
                    If (iz+2*ir+il.Gt.n00max)    Cycle
                    If (il+(is+1)/2.Ne.k)        Cycle
                    If (Mod(iz+il,2).Ne.ip-1)    Cycle
                    Iall=Iall+1
                 End Do
              End Do
           End Do
        End Do
     End Do
  End Do
  !-----------------------------------------------
  ! charge all energies for n00max
  !-----------------------------------------------
  Allocate(e(Iall))
  hbz=two*hbzero/bz**2; hbp=two*hbzero/bp**2; 
  Iall=0; 
  Do k=1,n00max+1
     la=k-1; le=min0(n00max,k)
     Do ip=1,2
        Do ir=0,nre
           Do iz=0,nze
              Do il=la,le
                 Do is=+1,-1,-2
                    If (iz+2*ir+il.Gt.n00max)    Cycle
                    If (il+(is+1)/2.Ne.k)        Cycle
                    If (Mod(iz+il,2).Ne.ip-1)    Cycle
                    Iall=Iall+1
                    e(Iall)=hbz*(float(iz)+half)+hbp*(two*float(ir)+float(il)+one)
                 End Do
              End Do
           End Do
        End Do
     End Do
  End Do
  !-----------------------------------------------
  ! sort energies  and derive base cut-off energy
  !-----------------------------------------------
  Call ord(Iall,e); 
  EBASECUT=E(NOSCIL)+1.0D-5
  Deallocate(e) 
  !-----------------------------------------------
  ! calculate the actual states
  !-----------------------------------------------
  nze=n00max; nre=n00max/2; ib=0; ilauf=0; nzx=0; nrx=0; nlx=0; nqp=0; nuv=0
  ! loop over k-quantum number
  Do k=1,n00max+1
     la=k-1; le=min0(n00max,k)
     ! loop over parity
     If(.Not.Parity) jlauf=ilauf !Nop 
     Do ip=1,2
        If(Parity) jlauf=ilauf !Yesp 
        Do ir=0,nre
           Do iz=0,nze
              Do il=la,le
                 Do is=+1,-1,-2
                    If (iz+2*ir+il.Gt.n00max)    Cycle
                    If (il+(is+1)/2.Ne.k)        Cycle
                    If (Mod(iz+il,2).Ne.ip-1)    Cycle
                    ee=hbz*(float(iz)+half)+hbp*(two*float(ir)+float(il)+one)
                    If(ee.Le.EBASECUT) Then
                       ilauf=ilauf+1
                       nzx=Max(nzx,iz); nrx=Max(nrx,ir); nlx=Max(nlx,il)
                    Endif
                 End Do
              End Do
           End Do
        End Do
        If(Parity) Then                !Yesp
           If (ilauf.Gt.jlauf) Then    
              ib=ib+1
              nd=ilauf-jlauf
              ndx=Max(ndx,nd)
              nqp=nqp+nd; nuv=nuv+nd*nd       
           End If
        Endif
     End Do
     If(.Not.Parity) Then              !Nop 
        If(ilauf.Gt.jlauf) Then  
           ib=ib+1
           nd=ilauf-jlauf
           ndx=Max(ndx,nd)
           nqp=nqp+nd; nuv=nuv+nd*nd       
        End If
     Endif
  End Do
  nbx=ib; ntx=ilauf
  !-----------------------------------------------
  ! print statistics 
  !-----------------------------------------------
  If(lpr) Then 
     Do iw=lout,lfile
        Write(iw,'(a,2(i6,2x),a)') '  NUV, NQP:                      ',nuv,nqp
        Write(iw,'(a,2(i6,2x),a)') '  Comparison with bookkeeping spherical basis:'
        Write(iw,'(a,2(i6,2x),a)') '  n00:                           ',n00,n00,  &
             'maximal number of shells'
        Write(iw,'(a,2(i6,2x),a)') '  nbx, 2*n00+1:                  ',nbx,2*n00+1,  &
             'maximal number of k-blocks'
        Write(iw,'(a,2(i6,2x),a)') '  ntx, (n00+1)*(n00+2)*(n00+3)/6 ',ntx,(n00+1)*(n00+2)*(n00+3)/6,  &
             'max.num. p/n levels'
        Write(iw,'(a,2(i6,2x),a)') '  nzx, n00:                      ',nzx,n00,  &
             'maximal nz-quantum number'
        Write(iw,'(a,2(i6,2x),a)') '  nrx, n00/2  :                  ',nrx,n00/2,  &
             'maximal nr-quantum number'
        Write(iw,'(a,2(i6,2x),a)') '  nlx, n00:                      ',nlx,n00,  &
             'maximal ml-quantum number'
        Write(iw,'(a,2(i6,2x),a)') '  ndx, (n00+2)*(n00+2)/4:        ',ndx,(n00+2)*(n00+2)/4,  &
             'maximal dim. of one k-block'
        Write(iw,*)
     Enddo
  Endif
  !
End Subroutine base0
!=======================================================================
!
!=======================================================================
Subroutine base(lpr)
  !---------------------------------------------------------------------
  ! set HO basis configurations in cylindrical coordinates
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Logical :: lpr
  Integer(ipr) :: nze,nre,ib,ilauf,jlauf,nom,nnm,  &
       k,la,le,ip,ir,iz,il,is,nn,ND,IBX,N1,N2,iw
  Real(pr) :: hbz,hbp,ee
  !
  hbz=two*hbzero/bz**2; hbp=two*hbzero/bp**2; 
  !
  nze=n00max; nre=n00max/2; ib=0; ilauf=0; nzm=0; nrm=0; nlm=0; nom=0; nnm=0
  !-----------------------------------------------
  ! loop over k-quantum number
  !-----------------------------------------------
  Do k=1,n00max+1
     la=k-1; le=min0(n00max,k)
     ! loop over parity
     If(.Not.Parity) jlauf=ilauf !Nop 
     Do ip=1,2
        If(Parity) jlauf=ilauf   !Yesp 
        Do ir=0,nre
           Do iz=0,nze
              Do il=la,le
                 Do is=+1,-1,-2
                    If (iz+2*ir+il.Gt.n00max) Cycle
                    If (il+(is+1)/2.Ne.k)     Cycle
                    If (Mod(iz+il,2).Ne.ip-1) Cycle
                    ee=hbz*(float(iz)+half)+hbp*(two*float(ir)+float(il)+one)
                    If(ee.Le.EBASECUT) Then
                       ilauf=ilauf+1
                       If (ilauf.Gt.ntx) Then
                          ierror_flag=ierror_flag+1
                          ierror_info(ierror_flag)='STOP: in base: ntx too small'
                          Return  
                       Endif
                       nz(ilauf)=iz; nr(ilauf)=ir; nl(ilauf)=il; ns(ilauf)=is; npar(ilauf)=ip					   
                       nn =iz+2*ir+il
                       Write(tb(ilauf),100) 2*k-1,tp(ip),nn,iz,il
100                    Format(i2,a1,'[',i2,',',i2,',',i2,']')
                       Do iw=lout,lfile
			              If(lpr.And.IDEBUG.Gt.10) Write(iw,'(i4,a,i2,a,i2,a,i2,a,i2,a,i2,a,2x,a,1x,a,f14.8)')  &
                               ilauf,'   nn=',nn,'   nz=',iz,'   nr=',ir,  &
                               '   ml=',il,'  ms=',is,' /2',tb(ilauf),'e=',ee
                       Enddo
                       nzm=Max(nzm,iz); nrm=Max(nrm,ir); nlm=Max(nlm,il)
                       nom=Max(nom,2*k-1); nnm=Max(nnm,iz+2*ir+il)
                    Endif
                 End Do
              End Do
           End Do
        End Do
        !-----------------------------------------------
        ! Block memory
        !-----------------------------------------------
        If(Parity) Then                !Yesp        
           If (ilauf.Gt.jlauf) Then   
              ib=ib+1
              ia(ib)=jlauf; id(ib)=ilauf-jlauf
              ikb(ib)=k; ipb(ib)=ip
              Write(txb(ib),'(i3,a,i2,a,a1)') ib,'. block:  k=',k+k-1,'/2',tp(ip)              
              !ir=(ib+1)/2
              !write(*,*)  ib,2*k-1,'2*Omega=',2*ir - 1  
              Do iw=lout,lfile
                 If(lpr.And.IDEBUG.Gt.10) Write(iw,'(/,a,i3,a,a1)')'  For the above block:  k=',k+k-1,'/2',tp(ip)
              End Do
           End If
           If(id(ib).Eq.0) Then
              ierror_flag=ierror_flag+1
              ierror_info(ierror_flag)='STOP: in base Block Memory(1)'
              Return  
           Endif
        Endif
     End Do
     !-----------------------------------------------
     ! Block memory
     !-----------------------------------------------
     If(.Not.Parity) Then               !Nop
        If (ilauf.Gt.jlauf) Then   
           ib=ib+1
           ia(ib)=jlauf; id(ib)=ilauf-jlauf
           ikb(ib)=k; ipb(ib)=ip           
!           Write(txb(ib),'(i3,a,i2,a,a1)') ib,'. block:  k=',k+k-1,'/2',tp(ip)
           Write(txb(ib),'(i3,a,i2,a)') ib,'. block:  k=',k+k-1,'/2'
           Do iw=lout,lfile
              If(lpr.And.IDEBUG.Gt.10) Write(iw,'(/,a,i3,a,a1)')'  For the above block:  k=',k+k-1,'/2',tp(ip)
           End Do
        End If
        If(id(ib).Eq.0) Then
           ierror_flag=ierror_flag+1
           ierror_info(ierror_flag)='STOP: in base Block Memory(2)'
           Return  
        Endif
     Endif
  End Do
  nb=ib;  nt=ilauf
  !-----------------------------------------------
  ! broyden/linear mixing (storage)
  !-----------------------------------------------
  nhhdim=0
  Do ib=1,NB
     ND=ID(ib)
     Do N1=1,ND
        Do N2=1,N1
           nhhdim=nhhdim+1
        Enddo
     Enddo
  Enddo
  nhhdim2=2*nhhdim; nhhdim3=3*nhhdim; nhhdim4=4*nhhdim
  If(Allocated(brin)) Deallocate(brin,brout)
  Allocate(brin(nhhdim4),brout(nhhdim4))    
  !-----------------------------------------------
  ! Print statistics
  !-----------------------------------------------
  If(lpr) Then
     Do iw=lout,lfile
        Write(iw,'(a,i4)')   '  actual basis used:'
        Write(iw,'(a,i4)')   '  number of blocks: nb=',nb
        Write(iw,'(a,i4)')   '  number of levels: nt=',nt
        Write(iw,'(a,i4)')   '  maximal 2*omega : nom=',nom
        Write(iw,'(a,i4)')   '  maximal nz:       nzm=',nzm
        Write(iw,'(a,i4)')   '  maximal nr:       nrm=',nrm
        Write(iw,'(a,i4)')   '  maximal ml:       nlm=',nlm
        Write(iw,'(a,i4)')   '  maximal nn=nz+2*nr+nl=',nnm
        Write(iw,'(a,i4)')   '  2 x bigest block dim.=',ndx2
        Write(iw,'(a,i8)')   '  Nonzero elemets of HH=',nhhdim
        Write(iw,'(a,i8)')   '  Number Broyden elemens=',nhhdim4
        Write(iw,'(a,i4)') 
     Enddo
  Endif
  If(nzm.Ge.n00max.Or.(nom-1)/2.Eq.n00max) Then
     Write(*,*) 'nzm=',nzm,'  (nom-1)/2=',(nom-1)/2,'  n00max=',n00max
     ierror_flag=ierror_flag+1
     ierror_info(ierror_flag)='STOP: Please increase n00max to have correct basis'
  Endif
End Subroutine base
!=======================================================================
!
!=======================================================================
Subroutine gaupol(lpr)
  !---------------------------------------------------------------------
  ! HO wave functions in cylindrical coordinates
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Logical :: lpr
  Real(pr) :: w0,z,x,s,s0,s1,w00,w4pii,dsq,d1,d2,d3,d4,hs0,hs1
  Integer(ipr) :: ih,il,iw,ix,n,l,n1,n2
  !-----------------------------------------------
  ! Hermit
  !-----------------------------------------------
  w4pii=pi**(-0.250_pr)
  Do ih=1,ngh
     z=xh(ih); w0=w4pii*Exp(-half*z*z)
     ! functions qh, qh1; norm: \sum_{ih} qh(n1,ih)*qh(n2,ih)=\delta_{n1,n2}
     w0 =w0*Sqrt(wh(ih))
     qh(0,ih)=w0;      qh(1,ih)=sq(2)*w0*z
     qh1(0,ih)=-w0*z;   qh1(1,ih)=sq(2)*w0*(one-z*z)
     Do n=2,nzm
        qh(n,ih)=sqi(n)*(sq(2)*z*qh(n-1,ih)-sq(n-1)*qh(n-2,ih))
        qh1(n,ih)=sq(n+n)*qh(n-1,ih)-z*qh(n,ih)
     End Do
  End Do
  !-----------------------------------------------
  ! Laguerre
  !-----------------------------------------------
  Do il=1,ngl
     x=xl(il); w00=sq(2)*Exp(-half*x)
     Do l=0,nlm
        ! functions ql, ql1; norm: \sum_{il} ql(n1,l,il)*ql(n2,l,il)=\delta_{n1,n2}
        w0=w00*Sqrt(half*wl(il)*x**l)
        ql(0,l,il)=wfi(l)*w0;         ql(1,l,il)=(l+1-x)*wfi(l+1)*w0
        ql1(0,l,il)=(l-x)*wfi(l)*w0;  ql1(1,l,il)=(Dble(l*l+l)-x*Dble(l+l+3)+x*x)*wfi(l+1)*w0
        Do n=2,nrm
           dsq=sq(n)*sq(n+l); d1=Dble(n+n+l-1)-x
           d2=sq(n-1)*sq(n-1+l); d3=n+n+l-x; d4=two*dsq
           ql(n,l,il)=(d1*ql(n-1,l,il)-d2*ql(n-2,l,il))/dsq
           ql1(n,l,il)=d3*ql(n,l,il)-d4*ql(n-1,l,il)
        End Do
     End Do
  End Do
  !-----------------------------------------------
  ! test accuracy for hermit orthonormalization
  !-----------------------------------------------
  hs0=0.0_pr; hs1=2.0_pr
  Do n1=0,nzm
     Do n2=0,n1
        If (Mod(n1-n2,2).Eq.0) Then
           s=zero
           Do ih=1,ngh
              s=s+qh(n1,ih)*qh(n2,ih)
           End Do
           If(n1.Ne.n2) Then
              hs0=Max(s,hs0)
           Else
              hs1=Min(s,hs1)
           End If
        End If
     End Do
  End Do
  !-----------------------------------------------
  ! test accuracy for laguerre orthonormalization
  !-----------------------------------------------
  s0=0.0_pr; s1=2.0_pr
  Do l=0,nlm
     Do n1=0,nrm
        Do n2=0,n1
           s=zero
           Do il=1,ngl
              s=s+ql(n1,l,il)*ql(n2,l,il)
           End Do
           If(n1.Ne.n2) Then
              s0=Max(s,s0)
           Else
              s1=Min(s,s1)
           End If
        End Do
     End Do
  End Do
  !-----------------------------------------------
  ! print accuracy
  !-----------------------------------------------
  If(lpr) Then
     Do iw=lout,lfile
        Write(iw,'(a,i3)')  &
             '  number of gauss-hermite  meshpoints ngh=',ngh
        Write(iw,'(2(a,f12.8))')  &
             '  from xh(1)=',xh(1),' to xh(ngh)=',xh(ngh)
        Write(iw,'(a,i3)')  &
             '  number of gauss-laguerre meshpoints ngl=',ngl
        Write(iw,'(2(a,f12.8))')  &
             '  from xl(1)=',xl(1),' to xl(ngl)=',xl(ngl)
        Write(iw,*)  &
             ' max.dev.in:    orthogonality           normalization'
        Write(iw,*) ' hermit   ',hs0,Abs(1.0_pr-hs1)
        Write(iw,*) ' laguerre ',s0,Abs(1.0_pr-s1)
     End Do
  Endif
  !-----------------------------------------------
  ! debug
  !-----------------------------------------------
  If (lpr.And.IDEBUG.Gt.20) Then
     ix=3
     Do iw=lout,lfile
        Write(iw,*) ' nz    qh(nz,ih=1,...)'
        Do n=0,nzm
           Write(iw,'(i4,3f15.8)') n,(qh(n,ih),ih=1,ix)
           Write(iw,'(i4,3f15.8)') n,(qh1(n,ih),ih=1,ix)
           Write(iw,*) ' '
        End Do
        Do l=0,nlm
           Write(iw,*) ' nr ml    ql(nr,l,il=1,...)'
           Do n=0,nrm
              Write(iw,'(i4,i3,3f15.8)') n,l,(ql(n,l,il),il=1,ix)
              Write(iw,'(i4,i3,3f15.8)') n,l,(ql1(n,l,il),il=1,ix)
              Write(iw,*) ' '
           End Do
        End Do
     End Do
     !-----------------------------------------------
     ! test for hermit integration
     !-----------------------------------------------
     Do n1=0,nzm
        Do n2=0,n1
           If (Mod(n1-n2,2).Eq.0) Then
              s=zero
              Do ih=1,ngh
                 s=s+qh(n1,ih)*qh(n2,ih)
              End Do
              Do iw=lout,lfile
                 Write(iw,100) n1,n2,s
              End Do
100           Format(' gauss hermit: n1=',i3,'  n2=',i3,f20.8)
           End If
        End Do
     End Do
     !-----------------------------------------------
     ! test for orthonormalization
     !-----------------------------------------------
     Do l=0,nlm
        Do n1=0,nrm
           Do n2=0,n1
              s=zero
              Do il=1,ngl
                 s=s+ql(n1,l,il)*ql(n2,l,il)
              End Do
              Do iw=lout,lfile
                 Write(iw,101) l,n1,n2,s
101              Format(' gauss laguerre: l='  &
                      ,i2,' n1=',i3,'  n2=',i3,f20.8)
              End Do
           End Do
        End Do
     End Do
  End If
  !
  Call coordinateLST(.False.)  ! coordinate LST  
  !
End Subroutine gaupol
!=======================================================================
!
!=======================================================================
Subroutine FileLabels(NPRI,ININL,FILELABEL)
  !---------------------------------------------------------------------
  ! file labels, e.g., filelabel='s070_040' 
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: it,ininabs,ininl,nprt(2),NPRI(2)
  Character(1)  :: sinin
  Character(3)  :: snpr(2)
  Character(8)  :: filelabel
  !
  ininabs=iabs(ininl)
  If(ininabs.Eq. 3) sinin='o'
  If(ininabs.Eq. 2) sinin='p'
  If(ininabs.Eq. 1) sinin='s'
  If(ininabs.Eq. 4) sinin='t'
  !
  nprt=npri  
  Do it=itmin,itmax  
     If(npri(it).Ne.2*(npri(it)/2)) nprt(it)=nprt(it)+iparenti(it)  !iparent=-/+ menas particles/holes
     If(nprt(it).Lt.10  ) Then
        Write(snpr(it),'(a2,i1)') '00',nprt(it)
     Else
        If(nprt(it).Lt.100 ) Then
           Write(snpr(it),'(a1,i2)') '0',nprt(it)
        Else
           Write(snpr(it),'(i3)') nprt(it)
        End If
     End If
  Enddo
  !
  Write(filelabel,'(a1,a3,a1,a3)')  sinin,snpr(1),'_',snpr(2)
  !
End Subroutine FileLabels
!=======================================================================
!
!=======================================================================
Subroutine inout(is)
  !---------------------------------------------------------------------
  ! is=1: reads matrix elements from tape and exit
  ! is=2: writes matrix elements to tape and exit
  ! NB! if the welfile is missing or corrupt call start
  !     to restart calculations from scratch   
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: is,iw,n1,n2,nd,ib,bloall1
  Character(2)  :: nucname1
  Character(8)  :: filelabel
  Character(12) :: welfile
  Real(pr)    :: tz1(2),b01,beta1,v0r(2),v1r(2),pwir
  Integer(ipr) :: npr1,npr11,ngh1,ngl1,n001,nt1
  Integer(ipr) :: ntx1,nb1,nhhdim1,NLANSA0,NLANSA1,NZA2NRA,NZA1,NLA1
  Integer(ipr) :: ID1(nbx)
  !==== HFBODD interface
  Integer(ipr) :: i,nza,nra,nla,nsa,ibasis
  Integer(ipr) :: ibro
  !====
  ! label organization
  Call FileLabels(NPR,ININ,FILELABEL)
  If(ierror_flag.Ne.0) Return	   	   
  If(iLST1.Le.0) Write(welfile,'(a8,a4)')  FILELABEL,'.hel'
  If(iLST1.Gt.0) Write(welfile,'(a8,a4)')  FILELABEL,'.tel'
  !
  If (is.Eq.1) Then
     !---------------------------------------------------------------------
     ! read matrix elements from 'welfile' file or start from scratch 
     !---------------------------------------------------------------------
     If(inin.Gt.0) Then
        Call Start
        Return
     End If
     !read matrix elements
     Open(lwin,file=welfile,status='unknown',form='unformatted', ERR=100)
     !hel
     Read(lwin,Err=100,End=100) nucname1,npr11,npr1,ngh1,ngl1,n001,nb1,nt1
     If(Abs(n001).Ne.Abs(n00).And.nb1.Ne.nb) go to 100
     Read(lwin,Err=100,End=100) b01,beta1,si,etot,rms,bet,xmix,v0r,v1r,pwir  &
          ,del,ept,ala,ala2,alast,tz1,varmas,varmasNZ,pjmassNZ,ass,skass
     brin=zero; si=one; bbroyden='L'
     Read(lwin,Err=100,End=100) ntx1,nb1,nhhdim1
     Read(lwin,Err=100,End=100) id1
     Read(lwin,Err=100,End=100) brin
     !
     ! Add small pairing de=de+0.1 in the no-LN 
     ! case to prevent pairing collaps  
     If(kindhfb.Eq.1.And.Add_Pairing) Then
        ibro=0
        Do ib=1,NB
           ND=ID1(ib)
           I=ibro
           Do N1=1,ND
              Do N2=1,N1
                 I=I+1
                 brin(i+nhhdim2)=brin(i+nhhdim2)+0.10_pr 
                 brin(i+nhhdim3)=brin(i+nhhdim3)+0.10_pr 
              End Do !N2
           End Do !N1
           ibro=i
        End Do !IB
     Endif
     Do ib=1,NB
        ND=ID1(ib)
        Do N1=1,ND
           Read(lwin,Err=100,End=100) NLANSA0,NLANSA1,NZA2NRA,NZA1,NLA1 
        End Do
     End Do
     ! blocking
     Read(lwin,Err=100,End=100)  bloall1
     Read(lwin,Err=100,End=100)  bloblo,blo123,blok1k2,blomax,bloqpdif
     If(bloall1.Ne.bloall) go to 100
     !tel
     If(iLST.Gt.0) Then
        Read(lwin,Err=100,End=100) decay,rmm3,cmm3,amm3,bmm3,itass,iqqmax
        If(Allocated(fdsx)) Deallocate(fdsx,fdsy,fdsy1,fdsy2,fdsy3,  &
             fspb0,fspc0,fspd0,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,  &
             fspb3,fspc3,fspd3)
        Allocate(fdsx(iqqmax),fdsy(iqqmax),fdsy1(iqqmax),  &
             fdsy2(iqqmax),fdsy3(iqqmax),fspb0(iqqmax),fspc0(iqqmax),  &
             fspd0(iqqmax),fspb1(iqqmax),fspc1(iqqmax),fspd1(iqqmax),  &
             fspb2(iqqmax),fspc2(iqqmax),fspd2(iqqmax),fspb3(iqqmax),  &
             fspc3(iqqmax),fspd3(iqqmax))
        Read(lwin,Err=100,End=100) fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,fspd0  &
             ,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
     End If
     Do iw=lout,lfile
        Write(iw,*)
        Write(iw,*) ' Reading from wel_file: ',welfile
        Write(iw,*)
     Enddo
     Close(lwin)
     !---------------------------------------------------------------------
     !for HFODD interface 
     !ib=abs(inin); npr_temp(1)=npr11; npr_temp(2)=npr1; 
     !Call HFBTHO_HFODD(1616,npr_temp,ib)
     !---------------------------------------------------------------------
     Return
     !
100  Continue
     !---------------------------------------------------------------------
     ! missing or corrupt 'welfile' file 
     !---------------------------------------------------------------------
     Close(lwin)
     Do iw=lout,lfile          
        Write(iw,'(1x,a,a,a)')
        Write(iw,'(1x,a,a,a)')   ' The file ',welfile,' is corrupted, missing, or wrong!' 
        Write(iw,'(1x,a,a,a)')   ' STARTING FROM SCRATCH WITH ININ=IABS(ININ)!' 
        Write(iw,'(1x,a,a,a)')
     Enddo
     Call start
     Return
  End If
  !---------------------------------------------------------------------
  ! write matrix elements to 'welfile' file 
  !---------------------------------------------------------------------
  If (is.Eq.2.And.iasswrong(3).Eq.0) Then
  ! don't write with mpi qrpa
#ifdef hide_mpi_qrpa
     Open(lwou,file=welfile,status='unknown',form='unformatted')
     !hel
     npr11=npr(1); npr1=npr(2)
     Write(lwou) nucname,npr11,npr1,ngh,ngl,n00,nb,nt 
     Write(lwou) b0,beta0,si,etot,rms,bet,xmix,CpV0,CpV1,pwi,del,ept  &
          ,ala,ala2,alast,tz,varmas,varmasNZ,pjmassNZ,ass,skass
     Write(lwou) ntx,nb,nhhdim
     Write(lwou) id
     Write(lwou) brin
     ibasis=0
     Do ib=1,NB
        ND=ID(ib)
        Do N1=1,ND
           ibasis=ibasis+1
           NLA=NL(ibasis); NRA=NR(ibasis); NZA=NZ(ibasis); NSA=NS(ibasis); NLANSA1=(-1)**(NZA+NLA)
           Write(lwou) 2*NLA+NSA,NLANSA1,NZA+2*NRA+NLA,NZA,NLA             
        End Do
     End Do
     !---------------------------------------------------------------------
     ! blocking: sort blocking candidates first
     !---------------------------------------------------------------------
     Do ib=1,2
        Call blosort(ib,blomax(ib))
     Enddo
     Write(lwou) bloall
     Write(lwou) bloblo,blo123,blok1k2,blomax,bloqpdif
     !tel
     If(Allocated(fdsx)) Then
        Write(lwou) decay,rmm3,cmm3,amm3,bmm3,itass,iqqmax
        Write(lwou) fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,fspd0  &
             ,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
     End If
     Close(lwou)
     Do iw=lout,lfile
        Write(iw,'(a,a,a)')
        Write(iw,'(a,a,a)') '  Writing to wel_file: ',welfile
        Write(iw,'(a,a,a)') ' __________________________________  '
        Write(iw,'(a,a,a)') '  The tape ',welfile,' recorded:     '
        Write(iw,'(a,a,a)') '  nucname,npr,ngh,ngl,n00,nb,nt      '
        Write(iw,'(a,a,a)') '  b0,beta0,si,etot,rms,bet,xmix      '
        Write(iw,'(a,a,a)') '  pairing:     CpV0,CpV1,pwi         '
        Write(iw,'(a,a,a)') '  delta:       del,ept               '
        Write(iw,'(a,a,a)') '  lambda:      ala,ala2,alast,tz     '
        Write(iw,'(a,a,a)') '  asymptotic:  varmas,ass,skass      '
        Write(iw,'(a,a,a)') '  ntx,nb,nhhdim,id,N_rz,n_r,n_z      '
        Write(iw,'(a,a,a)') '  Omega2,Sigma2,Parirty,Lambda       '
        Write(iw,'(a,a,a)') '  matrices(inbro):    hh,de          '
        Write(iw,'(a,a,a)') '  *all blocking candidates           '
        If(Allocated(fdsx)) Write(iw,'(a,a,a)') '  *all THO arrays                    '
        Write(iw,'(a,a,a)') ' __________________________________  '
        Write(iw,'(a,a,a)')
     Enddo
#else 
   ! with qrpa no hel file: do nothing here
#endif
  Endif
  !
  !ib=abs(inin); npr_temp(1)=npr11; npr_temp(2)=npr1; 
  !Call HFBTHO_HFODD(1717,npr_temp,ib)
  !
End Subroutine inout
!=======================================================================
!
!=======================================================================
Subroutine start
  !---------------------------------------------------------------------
  ! initializes scratch Saxon-Woods potentials
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: iw,i,ih,il,ihl,it,ita
  Real(pr) :: zb(ngh),rrb(ngl),rb(ngl),rav,rao,vpws,vls,betas,gamma  &
       ,fac,facb,zz,rr,r,ctet,cphi,p2,p20,p22,s,u,w,f,rc,c,beta00
  !--------------------------------------------------------------------
  ! initializing
  !--------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! reinatializing all again since scratch calculation
  !----------------------------------------------------------------------------
  Call iniialize_HFBTHO_SOLVER
  If(ierror_flag.Ne.0) Return	   	   
  !mario
  Do it=itmin,itmax  
     If(npr(it).Ne.2*(npr(it)/2)) Then
        irestart=irestart+1; npr(it)=npr_INI(it)     
     Endif
  Enddo
  If(irestart.Ne.0) Then
     ! odd nucleus requested but no even-even solution
     ! recalculate the even-even nucleus from scratch
     Do iw=lout,lfile         
        Write(iw,'(1x,a,2i4)') 
        Write(iw,'(1x,a,2i4)') ' Initialization for the even-even core (N,Z)=: ',npr(1:2) 
     Enddo
  Else
     ! scratch for the even-even nucleus requested
     Do iw=lout,lfile
        Write(iw,'(1x,a,2i4)')
        Write(iw,'(a,a,3i4)')    '  Scratch initialization for the nucleus: ',nucname,npr(1:2) 
        Write(iw,'(1x,a,2i4)')
     Enddo
  Endif
  n00=Abs(n00_INI);       b0=b0_INI;           q=q_INI; iLST=iLST_INI
  maxi=MAX_ITER_INI;      inin=inin_INI;      
  skyrme=skyrme_INI;      kindhfb=kindhfb_INI
  cdef=cdef_INI;          cqad=cqad_INI;       iproj=iproj_INI 
  npr1pj=npr1pj_INI;      npr2pj=npr2pj_INI;   
  Call Constraint_or_not(inin_INI,inin,icstr) 
  If(ierror_flag.Ne.0) Return	   	             
  Call preparer(.False.)
  If(ierror_flag.Ne.0) Return	   	             
  inin=Abs(inin)          ! positive even if inin_INI is not
  bet=cdef
  !-----------------------------------
  ! Saxon-Woods potentials
  !-----------------------------------
  Do iw=lout,lfile
     Write(iw,'(/,a)') '  Initial potentials of Saxon-Woods shape '
  End Do
  beta00=bet     ! wf to requested deformation
  Do iw=lout,lfile          
     Write(iw,'(a,f10.4)')  '  v0ws   =',v0ws
     Write(iw,'(a,f10.4)')  '  kappa  =',akv
     Write(iw,'(a,2f10.4)') '  vs0    =',vso
     Write(iw,'(a,2f10.4)') '  r0     =',r0v
     Write(iw,'(a,2f10.4)') '  a      =',av
     Write(iw,'(a,2f10.4)') '  r0-so  =',rso
     Write(iw,'(a,2f10.4)') '  a-so   =',aso
     Write(iw,'(a,f10.4)')  '  beta00 =',beta00
  Enddo
  !-----------------------------------
  ! Densities
  !-----------------------------------
  Do it=itmin,itmax
     ita=3-it; rav=r0v(it)*amas**p13; rao=rso(it)*amas**p13
     vpws=v0ws*(one-akv*(npr(it)-npr(ita))/amas)
     vls=half*(hqc/amu)**2*vpws*vso(it)
     betas=beta00 * Sqrt(5.0_pr/(16.0_pr*pi))
     gamma=zero
     fac=one+betas*Cos(gamma*pi/180.0_pr)
     fac=(one+betas*Cos((gamma+120.0_pr)*pi/180.0_pr))*fac
     fac=(one+betas*Cos((gamma-120.0_pr)*pi/180.0_pr))*fac
     fac=fac**(-p13)
     ! z,r-coordinates in fm
     zb=xh*bz; rrb=xl*bp*bp; rb=Sqrt(rrb)
     Do ih=1,ngh
        zz=zb(ih)**2
        Do il=1,ngl
           rr=rrb(il)+zz; r=Sqrt(rr)
           ! woods saxon
           ctet=zz/rr; cphi=zero
           p20=3.0_pr*ctet-one; p22=Sqrt(3.0_pr)*cphi
           p2=p20*Cos(gamma*pi/180)+p22*Sin(gamma*pi/180.0_pr)
           facb=fac*(one+betas*p2)
           u=vpws/( one+Exp( (r-rav*facb) / av(it) ))
           w=-vls /( one+Exp( (r-rao*facb) / aso(it)))
           ihl=ih+(il-1)*ngh
           If(it.Eq.1) Then
              vhbn(ihl)=hb0; vn(ihl)=u; vsn(ihl)=w;              
              vrn(ihl)=zero; vzn(ihl)=zero; vdn(ihl)=zero; 
              vSFIZn(ihl)=zero; vSZFIn(ihl)=zero; 
              vSFIRn(ihl)=zero; vSRFIn(ihl)=zero; 
           Else 
              vhbp(ihl)=hb0; vp(ihl)=u; vsp(ihl)=w;              
              vrp(ihl)=zero; vzp(ihl)=zero; vdp(ihl)=zero; 
              vSFIZp(ihl)=zero; vSZFIp(ihl)=zero; 
              vSFIRp(ihl)=zero; vSRFIp(ihl)=zero; 
           End If
           ro(ihl,it)=u
           aka(ihl,it)=1.d-3*Exp((r-rav*facb)/2._pr)
        End Do
     End Do
     s=npr(it)/Sum(ro(:,it))
     Do il=1,ngl
        Do ih=1,ngh
           ihl=ih+(il-1)*ngh
           f=s/(pi*wh(ih)*wl(il)* bz*bp*bp); ro(ihl,it)=f*ro(ihl,it)
        End Do
     End Do
     !-----------------------------------
     ! pairing
     !-----------------------------------
     Do il=1,nghl
        f=(ro(il,1)+ro(il,2))/rho_c
        If(it.Eq.1) Then
           dvn(il)=CpV0(it-1)*(one-CpV1(it-1)*f)*aka(il,it)
        Else
           dvp(il)=CpV0(it-1)*(one-CpV1(it-1)*f)*aka(il,it)
        End If
     End Do
  End Do
  !-----------------------------------
  ! coulomb
  !-----------------------------------
  If (icou.Eq.0) Then
     cou=zero
  Else
     rc=r0v(2)*amas**p13
     Do il=1,ngl
        Do ih=1,ngh
           r=Sqrt(zb(ih)**2+rrb(il))
           If (r.Lt.rc) Then
              c=half*(3/rc-r*r/(rc**3))
           Else
              c=one/r
           End If
           cou(ih+(il-1)*ngh)=c*npr(2)/alphi
        End Do
     End Do
  End If
  !-----------------------------------
  ! initial ph+pp matrix elements
  !-----------------------------------
  ak=0.1_pr; rk=0.1_pr ! initial density matrix elements (improve later)
  brin=zero            ! initial matrix elements to zero
  iiter=0              ! iteration number iiter to zero 
  Call gamdel    
  !
End Subroutine start
!=======================================================================
!
!=======================================================================
Subroutine printRHO
  !---------------------------------------------------------------------
  ! prints rho, aka
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr), Save :: ihli,Ifle
  Ifle=76+iLST1
  If(Ifle.Eq.76) Open(Ifle,file='ho_den.dat',status='unknown')
  If(Ifle.Eq.77) Open(Ifle,file='tho_den.dat',status='unknown')
  Write(Ifle,*) 'r  denN  denP  akaN  akaP '
  Do ihli=1,nghl
     Write(Ifle,'(12(1x,e16.8))') Sqrt(fh(ihli)**2+fl(ihli)**2)  &
          ,ro(ihli,1),ro(ihli,2),aka(ihli,1),aka(ihli,2)
  End Do
  Close(Ifle)
End Subroutine printRHO
!=======================================================================
!
!=======================================================================
Subroutine gfv
  !---------------------------------------------------------------------
  ! Calculates sign, Sqrt, factorials, etc. of integers and half int.
  ! iv(n)=(-1)**n, sq(n)=Sqrt(n), sqi(n)=1/Sqrt(n)
  ! fak(n)=n!; wf(n)=Sqrt(n!); wfi(n)=1/Sqrt(n!)
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: i,igfv
  Parameter(igfv=170)               !maximal number for GFV
  If(Allocated(iv)) Deallocate(iv,fak,fi,sq,sqi,wf,wfi)
  Allocate(iv(-igfv:igfv),fak(0:igfv),fi(0:igfv),sq(0:igfv),sqi(0:igfv))
  Allocate(wf(0:igfv),wfi(0:igfv))
  iv(0)=1; sq(0)=zero; sqi(0)=1.0d30
  fak(0)=one; fi(0)=one; wf(0)=one; wfi(0)=one
  Do i=1,igfv
     iv(i)=-iv(i-1)
	 iv(-i) = iv(i)	 
     sq(i)=Sqrt(Dble(i)); sqi(i)=one/sq(i)
     fak(i)=Dble(i)*fak(i-1); fi(i)=one/fak(i)
     wf(i)=sq(i)*wf(i-1); wfi(i)=one/wf(i)
  End Do
End Subroutine gfv
!=======================================================================
!
!=======================================================================
Subroutine sdiag(nmax,n,a,d,x,e,is)
  !---------------------------------------------------------------------
  ! A   matrix to be diagonalized
  ! D   eigenvalues,  X   eigenvectors, E   auxiliary field
  ! IS=1  eigenvalues are ordered (major component of X is positive)
  ! 0  eigenvalues are not ordered
  !---------------------------------------------------------------------
  Use HFBTHO, Only: pr,ipr
  Implicit None
  Integer(ipr), Save :: i,j,j1,k,l,im
  Integer(ipr)       :: n,nmax,is
  Real(pr), Save :: f,g,h,hi,s,p,b,r,pra,c
  Real(pr) :: a(nmax,nmax),x(nmax,nmax),e(n),d(n)
  Real(pr), Save :: tol=1.0D-32,eps=9.0D-12,one=1.0_pr,zero=0.0_pr
  !
  If (n.Le.1) Then
     d(1)=a(1,1); x(1,1)=one
     Return
  End If
  Do i=1,n
     Do j=1,i
        x(i,j)=a(i,j)
     End Do
  End Do
  ! householder-reduktion
  i=n
15 Continue
  If (i.Ge.2) Then
     l=i-2
     f=x(i,i-1); g=f; h=zero
     If (l.Gt.0) Then
        Do k=1,l
           h=h+x(i,k)*x(i,k)
        End Do
     Endif
     s=h+f*f
     If (s.Lt.tol) Then
        h=zero
        Go To 100
     Endif
     If (h.Gt.zero) Then
        l=l+1; g=Sqrt(s)
        If (f.Ge.zero) g=-g
        h=s-f*g; hi=one/h; x(i,i-1)=f-g; f=zero
        If (l.Gt.0) Then
           Do j=1,l
              x(j,i)=x(i,j)*hi
              s=zero
              Do k=1,j
                 s=s+x(j,k)*x(i,k)
              End Do
              j1=j+1
              If (l.Ge.j1) Then
                 Do k=j1,l
                    s=s+x(k,j)*x(i,k)
                 End Do
              Endif
              e(j)=s*hi; f=f+s*x(j,i)
           End Do
        Endif
        f=f*hi*0.50_pr
        If (l.Gt.0) Then
           Do j=1,l
              s=x(i,j); e(j)=e(j)-f*s; p=e(j)
              Do  k=1,j
                 x(j,k)=x(j,k)-s*e(k)-x(i,k)*p
              End Do
           End Do
        Endif
     Endif
100  Continue
     d(i)=h; e(i-1)=g; i=i-1
     Go To 15
     ! Bereitstellen der Transformationmatrix
  Endif
  d(1)=zero; e(n)=zero; b=zero; f=zero
  Do i=1,n
     l=i-1
     If (d(i).Eq.0.) Go To 221
     If (l.Gt.0) Then
        Do J=1,L
           s=zero
           Do k=1,l
              s=s+x(i,k)*x(k,j)
           End Do
           Do k=1,l
              x(k,j)=x(k,j)-s*x(k,i)
           End Do
        End Do
     Endif
221  Continue
     d(i)=x(i,i)
     x(i,i)=one
     If (l.Gt.0) Then
        Do j=1,l
           x(i,j)=zero; x(j,i)=zero
        End Do
     Endif
  End Do
  ! Diagonalisieren der Tri-Diagonal-Matrix
  Do l=1,n
     h=eps*(Abs(d(l))+ Abs(e(l)))
     If (h.Gt.b) b=h
     ! Test fuer Splitting
     Do  j=l,n
        If (Abs(e(j)).Le.b) Exit
     End Do
     ! test fuer konvergenz
     If (j.Eq.l) Go To 300
340  p=(d(l+1)-d(l))/(2.0_pr*e(l))
     r=Sqrt(p*p+one); pra=p+r
     If (p.Lt.zero) pra=p-r
     h=d(l)-e(l)/pra
     Do i=l,n
        d(i)=d(i)-h
     End Do
     f=f+h
     ! QR-transformation
     p=d(j); c=one; s=zero; i=j
360  i=i-1
     If (i.Lt.l) Go To 362
     g=c*e(i); h=c*p
     If ( Abs(p).Ge.Abs(e(i))) Then
        c=e(i)/p
        r=Sqrt(c*c+one); e(i+1)=s*p*r; s=c/r; c=one/r
        Go To 365
     Endif
     c=p/e(i)
     r=Sqrt(c*c+one); e(i+1)=s*e(i)*r; s=one/r; c=c/r
365  p=c*d(i)-s*g
     d(i+1)=h+s*(c*g+s*d(i))
     Do k=1,n
        h=x(k,i+1); x(k,i+1)=x(k,i)*s+h*c
        x(k,i)=x(k,i)*c-h*s
     End Do
     Go To 360
362  e(l)=s*p
     d(l)=c*p
     If ( Abs(e(l)).Gt.b) Go To 340
     ! konvergenz
300  d(l)=d(l)+f
  End Do
  If (is.Eq.0) Return
  ! ordnen der eigenwerte
  Do i=1,n
     k=i; p=d(i); j1=i+1
     If (j1.Le.n) Then
        Do j=j1,n
           If (d(j).Ge.p) Cycle
           k=j; p=d(j)
        End Do
        If (k.Eq.i) Cycle
        d(k)=d(i); d(i)=p
        Do j=1,n
           p=x(j,i); x(j,i)=x(j,k)
           x(j,k)=p
        End Do
     Endif
  End Do
  ! signum
  Do  k=1,n
     s=zero
     Do i=1,n
        h=Abs(x(i,k))
        If (h.Gt.s) Then 
           s=h; im=i
        End If
     End Do
     If (x(im,k).Lt.zero) Then
        Do i=1,n
           x(i,k)=-x(i,k)
        End Do
     End If
  End Do
End Subroutine sdiag
!=======================================================================
!
!=======================================================================
Subroutine nucleus(is,npr2,te)
  !---------------------------------------------------------------------
  ! is=1 determines the symbol for a given proton number npr2
  ! 2 determines the proton number for a given symbol te
  !---------------------------------------------------------------------
  Use HFBTHO, Only: pr,ipr,ierror_flag,ierror_info
  Implicit None
  Integer(ipr) :: is,npr2,np
  Integer(ipr) :: maxz
  Parameter (maxz=133)
  Character(2) te
  Character(2*maxz+2) t
  T(  1: 40)=' n HHELIBE B C N O FNENAMGALSI P SCLAR K'
  T( 41: 80)='CASCTI VCRMNFECONICUZNGAGEASSEBRKRRBSR Y'
  T( 81:120)='ZRNBMOTCRORHPDAGCDINSNSBTE IXECSBALACEPR'
  T(121:160)='NDPMSMEUGDTBDYHOERTMYBLUHFTA WREOSIRPTAU'
  T(161:200)='HGTLPBBIPOATRNFRRAACTHPA UNPPUAMCMBKCFES'
  T(201:220)='FMMDNOLR040506070809'
  T(221:265)='101112131415161718192021222324252627282930313233'
  If (is.Eq.1) Then
     If (npr2.Lt.0.Or.npr2.Gt.maxz) Then
        ierror_flag=ierror_flag+1
        ierror_info(ierror_flag)='STOP: in nucleus npr2 is wrong:'
        Return  
     Endif
     te=t(2*npr2+1:2*npr2+2) 
     Return
  Else
     Do np=0,maxz
        If (te.Eq.t(2*np+1:2*np+2)) Then
           npr2=np 
           Return
        End If
     End Do
  End If
  ierror_flag=ierror_flag+1
  ierror_info(ierror_flag)='STOP: in nucleus the nucleus is unknown!'
End Subroutine nucleus
!=======================================================================
!
!=======================================================================
Subroutine stab(npr2,npr3)
  !---------------------------------------------------------------------
  ! given 'Z' returns mass number 'A' on the stability line
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: npr2,npr3
  Real(pr), Save :: sn,sz,dsn,c,c5
  c=0.0060_pr; c5=5.0_pr*c/3.0_pr; sz=npr2; sn=npr2
  Do While(Abs(dsn).Lt.1.0d-5)
     dsn=sz-sn+c*(sn+sz)**(5.0_pr/3.0_pr)
     dsn=dsn/(-1.0_pr+c5*(sn+sz)**(2.0_pr/3.0_pr))
     sn=sn-dsn
  End Do
  npr3=sn
End Subroutine stab
!=======================================================================
!
!=======================================================================
Subroutine ord(n,e)
  !---------------------------------------------------------------------
  ! orders a set of numbers according to their size
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: n,i,k,j
  Real(pr), Save :: p
  Real(pr) :: e(n)
  Do i=1,n
     k=i; p=e(i)
     If (i.Lt.n) Then
        Do j=i+1,n
           If (e(j).Lt.p) Then
              k=j; p=e(j)
           End If
        End Do
        If (k.Ne.i) Then
           e(k)=e(i); e(i)=p
        End If
     End If
  End Do
End Subroutine ord
!=======================================================================
!
!=======================================================================
Subroutine blosort(it,n)
  !---------------------------------------------------------------------
  ! sorting blocking candidates
  ! Integer(ipr) :: iblocking,bloall; Parameter(bloall=200)
  ! Integer(ipr), Dimension(0:bloall,2) :: bloblo,blo123=0,blok1k2=0
  ! Real(pr),    Dimension(0:bloall,2) :: bloqpdif
  ! Integer(ipr), Dimension(3) :: keyblo
  ! Integer(ipr), Dimension(2) :: blocross,blomax,blo123d,blok1k2d,blocanon
  ! Write(lwou) bloblo,blo123,blok1k2,blomax,bloqpdif
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: it,ip,n,i,k,j
  Real(pr) :: p
  Do i=1,n
     k=i; p=bloqpdif(i,it)
     If (i.Lt.n) Then
        Do j=i+1,n
           If (bloqpdif(j,it).Lt.p) Then
              k=j; p=bloqpdif(j,it)
           End If
        End Do
        If (k.Ne.i) Then
           bloqpdif(k,it)=bloqpdif(i,it); bloqpdif(i,it)=p
           ip=bloblo(k,it);  bloblo(k,it)=bloblo(i,it);  bloblo(i,it)=ip   
           ip=blo123(k,it);  blo123(k,it)=blo123(i,it);  blo123(i,it)=ip   
           ip=blok1k2(k,it); blok1k2(k,it)=blok1k2(i,it); blok1k2(i,it)=ip   
        End If
     End If
  End Do
End Subroutine blosort
!=======================================================================
!
!=======================================================================
Subroutine tracesln
  !---------------------------------------------------------------------
  !        CALCULATING THE LIPKIN-NOGAMI SUMS IN CANONICAL BASIS
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: iw,it,ib,k1,k2,kkk,k
  Real(pr)    :: AAV,SNtor,SDtor
  Real(pr)    :: S_U1V1,S_U1V3,S_U2V2,S_U3V1,S_U4V4
  Real(pr)    :: U_ACTU,U_ACTU2,U_ACTU3,U_ACTU4
  Real(pr)    :: V_ACTU,V_ACTU2,V_ACTU3,V_ACTU4
  !
  etr=zero
  Do it=itmin,itmax
     S_U1V1=ZERO; S_U1V3=ZERO; S_U2V2=ZERO; S_U3V1=ZERO; S_U4V4=ZERO
     Do ib=1,nb
        k1=ka(ib,it)+1; k2=ka(ib,it)+kd(ib,it)
        If(k1.Le.k2) Then
           kkk=lcanon(ib-1,it)
           Do k=1,id(ib)
              kkk=kkk+1; aav=vk(kkk,it)               ! v^2
              U_ACTU=Sqrt(AAV);       U_ACTU2=U_ACTU*U_ACTU 
              U_ACTU3=U_ACTU2*U_ACTU; U_ACTU4=U_ACTU2*U_ACTU2
              V_ACTU=Sqrt(ONE-AAV);   V_ACTU2=V_ACTU*V_ACTU 
              V_ACTU3=V_ACTU2*V_ACTU; V_ACTU4=V_ACTU2*V_ACTU2
              S_U1V1=S_U1V1+U_ACTU  * V_ACTU
              S_U1V3=S_U1V3+U_ACTU  * V_ACTU3
              S_U2V2=S_U2V2+U_ACTU2 * V_ACTU2     !Tr r (1-r)
              S_U3V1=S_U3V1+U_ACTU3 * V_ACTU        
              S_U4V4=S_U4V4+U_ACTU4 * V_ACTU4     !Tr (1-r)^2 r^2
           End Do
        End If
     End Do !ib
     SNtor=8.0_pr*(S_U3V1*S_U1V3-S_U4V4)
     SDtor=32.0_pr*(S_U2V2*S_U2V2-S_U4V4)
     Geff(it)=del(it)**2/ept(it)
     ala2(it)=-Geff(it)*(SNtor/SDtor)
     If(ala2(it).Ge.10.) ala2(it)=4.                   ! ala2 goes to hell
     etr(it)=-four*ala2(it)*S_U2V2                  ! to total energy
  End Do !it
  etr(3)=etr(1)+etr(2)                                   !to total energy
  Do iw=lout,lfile
     Write(iw,'(26x,a,2(1x,f7.3),a,3(1x,f9.3),a,2(1x,f7.3))')  &
          '  #LN: ala2(n,p)=',ala2,' #eln(n,p,t)=',etr,' #del+ala2=',del+ala2
  Enddo
End Subroutine tracesln
!=======================================================================
!
!=======================================================================
Subroutine tracesln_qp
  !---------------------------------------------------------------------
  !        CALCULATING THE LIPKIN-NOGAMI TRACES IN QP SPACE
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: iw,nd,ib,i1,i2,n2,ibitnb,it,i1n2nd,i2n2nd,i1i2nd
  Real(pr) :: frit,frit2,ftit
  Real(pr) :: etr2(2),trk(2),trk1(2),SNtor(2),SDtor(2),Sum(2)
  !
  ! initialization
  etr=zero; etr2=zero; trk=zero; trk1=zero
  ! loop over the blocks
  Do ib=1,nb
     nd=id(ib)
     ! Traces for neutrons and protons
     Do i2=1,nd           ! index alpha
        Do i1=i2,nd         ! index beta.ge.alpha
           sum=zero
           Do n2=1,nd
              i1n2nd=Max(i1,n2)+(Min(i1,n2)-1)*nd
              i2n2nd=Max(i2,n2)+(Min(i2,n2)-1)*nd
              Do it=itmin,itmax
                 ibitnb=ib+(it-1)*nbx
                 Sum(it)=Sum(it)+rk(i1n2nd,ibitnb)*rk(i2n2nd,ibitnb)*p14
              End Do !it
           End Do !n2
           i1i2nd=i1+(i2-1)*nd
           Do it=itmin,itmax
              ibitnb=ib+(it-1)*nbx
              frit=rk(i1i2nd,ibitnb)*half
              ftit=ak(i1i2nd,ibitnb)
              frit2=Sum(it)
              If(i1.Eq.i2) Then                                        
                 etr(it)=etr(it)+frit-frit**2                  ! Tr r (1-r)
                 etr2(it)=etr2(it)+(one-two*frit+frit2)*frit2      ! Tr (1-r)^2 r^2
                 trk(it)=trk(it)+frit*ftit                       ! Tr r k
                 trk1(it)=trk1(it)+ftit -frit*ftit                 ! Tr k (1-r)
              Else
                 etr(it)=etr(it) -two*frit**2                     ! Tr r (1-r) 
                 etr2(it)=etr2(it)+two*(-two*frit+frit2)*frit2     ! Tr (1-r)^2 r^2
                 trk(it)=trk(it)+two*frit*ftit                   ! Tr r k
                 trk1(it)=trk1(it)-two*ftit*frit                   ! Tr k (1-r)
              End If
           End Do !it
        End Do !i1
     End Do !i2
  End Do !ib
  ! total traces
  Do it=itmin,itmax
     SNtor(it)=8.0_pr*(trk1(it)*trk(it)-etr2(it))
     SDtor(it)=32.0_pr*(etr(it)**2     -etr2(it))
     Geff(it)=del(it)**2/ept(it)
     ala2(it)=-( SNtor(it)/SDtor(it) )*Geff(it)     
     If(ala2(it).Ge.10.) ala2(it)=4.  ! in case ala2 goes to hell
     etr(it)=-four*ala2(it)*etr(it) ! to total energy
  End Do
  etr(3)=etr(1)+etr(2)         !to total energy  
  Do iw=lout,lfile
     Write(iw,'(26x,a,2(1x,f7.3),a,3(1x,f9.3),a,2(1x,f7.3))')  &
          '  #LN: ala2(n,p)=',ala2,' #eln(n,p,t)=',etr,' #del+ala2=',del+ala2
  Enddo
End Subroutine tracesln_qp  
!=======================================================================
!
!=======================================================================
Subroutine densitln
  !---------------------------------------------------------------------
  ! calculates the densities in r-space at gauss-meshpoints
  ! corrected due to Lipkin-Nogami
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: ih,il,ib,nd,i0,i01,i02,n1,n2,nza,nzb,nra,nrb,nla,  &
       nlb,nsa,nsb,it,ml,ihli,k,k0(2),k00(2),k1,k2
  Real(pr) :: fr(2),vvs,vvc,ssln1,ssln2,ssln3,vks
  Real(pr) :: qla,qlb,qlab,qha,qhb,qhlab,qhab,sro
  !
  k0=0; k00=0
  ro=zero
  ! loop over the blocks
  Do ib=1,nb
     k00=k0
     nd=id(ib); i0=ia(ib)
     Do n2=1,nd
        i02=i0+n2; nzb=nz(i02); nrb=nr(i02); 
        nlb=nl(i02); nsb=ns(i02)
        Do n1=1,n2
           i01=i0+n1; nza=nz(i01); nra=nr(i01)
           nla=nl(i01); nsa=ns(i01)
           k0=k00
           Do it=itmin,itmax
              k1=ka(ib,it)+1
              k2=ka(ib,it)+kd(ib,it)
              fr(it)=zero
              If(k1.Le.k2) Then
                 Do k=1,nd
                    k0(it)=k0(it)+1
                    ssln1=ssln(1,it)
                    ssln2=ssln(2,it)
                    ssln3=ssln(3,it)
                    vks=vk(k0(it),it)
                    vvc=vks
                    vvs=Abs(one-vks)
                    If(vvs.Ge.1.0d-40) Then
                       vvs=two*Sqrt(vks*vvs)   !2vu
                       vvc=vks+vvs**2*p14*ssln1*((two*vks-one)*ssln1-ssln2)/ssln3
                    End If
                    fr(it)=fr(it)+two*ddc(n2,k0(it),it)*ddc(n1,k0(it),it)*vvc
                 End Do
                 If (n1.Ne.n2) Then
                    fr(it)=two*fr(it)
                 End If
              End If
           End Do
           !---diagonal in spin
           If (nsa.Eq.nsb) Then
              ml=nla
              Do il=1,ngl
                 qla=ql (nra,ml,il);    qlb=ql (nrb,ml,il)
                 qlab=qla*qlb
                 Do ih=1,ngh
                    ihli=ih+(il-1)*ngh
                    qha=qh (nza,ih); qhb=qh (nzb,ih)
                    qhab=qha*qhb
                    qhlab=qhab*qlab; sro=qhlab
                    ro(ihli,:)=ro(ihli,:)+fr(:)*sro
                 End Do   !ih
              End Do  !il
           End If
        End Do !n2
     End Do !n1
  End Do !ib
  ! set the THO weights
  Do ihli=1,nghl
     ro(ihli,:)=ro(ihli,:)*wdcori(ihli)
  End Do
End Subroutine densitln
!=======================================================================
!
!=======================================================================
Subroutine coulom1
  !---------------------------------------------------------------------
  ! Coulom-field (direct part) Vautherin prescription
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr), Save :: i,k
  Real(pr), Save :: zd2,rhl,y1,y2,xx1,xx2,s1,s2,e1,e2,vik,f,r,r1,r4,  &
       rr2,z,z1,zd1
  Real(pr) x1,x2
  Real(pr), Save :: bb1=0.443251414630_pr,g1=0.249983683100_pr,  &
       bb2=0.062606012200_pr,g2=0.092001800370_pr,  &
       bb3=0.047573835460_pr,g3=0.040696975260_pr,  &
       bb4=0.017365064510_pr,g4=0.005264496390_pr
  Real(pr), Save :: small=1.0d-14
  !
  If(icacou.Eq.0) Then
     icacou=1  
     f=half*chargee2/pi
     Do i=1,nghl
        r=fl(i); z=fh(i)
        r4=four*r
        Do k=1,i
           r1=fl(k); z1=fh(k)
           rhl=r4*r1
           rr2=(r+r1)**2
           zd1=(z-z1)**2;  zd2=(z+z1)**2
           y1=zd1+rr2;     y2=zd2+rr2
           xx1=rhl/y1;     xx2=rhl/y2
           x1=one-xx1;     x2=one-xx2
           s1=Sqrt(y1);    s2=Sqrt(y2)
           e1=one+x1*(bb1+x1*(bb2+x1*(bb3+x1*bb4)))
           If(x1.Gt.small) e1=e1-x1*Log(x1)*(g1+x1*(g2+x1*(g3+x1*g4)))
           e2=one+x2*(bb1+x2*(bb2+x2*(bb3+x2*bb4)))
           If(x2.Gt.small) e2=e2-x2*Log(x2)*(g1+x2*(g2+x2*(g3+x2*g4)))
           vik=f*(e1*s1+e2*s2)
           vc(i,k)=vik*wdcor(k)  !wdcor=pi*wh*wl*bz*bp*bp/fd
           vc(k,i)=vik*wdcor(i)
        End Do  !k
     End Do  !i
  End If
  ! calculation of the coulomb field
  !cou=zero
  !Do i=1,nghl
  !   cou(:)=cou(:)+vc(:,i)*dro(i,2)
  !End Do
  Call dgemm('n','n',nghl,1,nghl,1.0_pr,vc,nghl,ro(:,2),nghl,0.0_pr,cou,nghl)
End Subroutine coulom1
!=======================================================================
!
!=======================================================================
Subroutine coulom
  !---------------------------------------------------------------------
  ! Coulomb field (direct part), Gogny prescription
  ! Ref.: Phys. Rev. C 27, 2317 (1983)
  ! taken from HFBTHO v200d 
  !---------------------------------------------------------------------
  !  Use HFBTHO_utilities
!  USE UNEDF, Only: ipr, pr
  USE HFBTHO !, Only: vc
!  USE HFBTHO_VERSION, Only: ipr, pr
!  Use HFBTHO Only vc
  Use bessik
  Implicit None
  Integer(ipr), Save :: i,j,k
  Real(pr), Save :: zd2,y1,y2,xx1,s1,vik,f,r,r1,fac1,fac2,rr2,z,z1,zd1,t,  &
                    bb,r2,r12,rrr,rz1,rz2,rrz1,rrz2,xx,rk1,rip1,rkp1,alpha,&
                    beta,xxx
  !
  If (IDEBUG.Eq.1) Call get_CPU_time('coulom',0)
  !
  If(icacou.Eq.0) Then
     !
     icacou=1
     !
     ! For parity-breaking shapes, the Coulomb potential was incorrectly
     ! calculated by assuming the two intervals [0,+\infty[ and ]-infty,0]
     ! were equivalent (see also below). This bug was corrected in version
     ! 139a
     If(Parity) Then
        fac1 = one;  fac2 = one
     Else
        fac1 = zero; fac2 = two
     End If
     ! Notes:
     !   - Missing factor 2 compared to Eq. (58) CPC paper because the density
     !     ro(:,it) already contains it (see routine DENSIT) due to T-invariance
     !   - Missing factor 1/2 when applying Gauss-Legendre quadrature (from [0,1]
     !     to the proper [-1,1] interval because it will be put back in subroutine
     !     expect() and is cancelled by a factor 2 in the HF field
     !   - For conserved parity, Gauss-Hermite points are all positive, the full
     !     integral over z' is split in z'<0 and z'>0, values of z and z1 below
     !     refer to the absolute values of z' (=-z' if z'<0)
     !
     bb=50.0_pr          ! Length scale L
     beta=2.00_pr
     alpha=one/beta
     f=chargee2/Sqrt(pi) ! e^2/Sqrt(pi)
     !
!$OMP PARALLEL DO        &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(nghl,fl,fh,nleg,xleg,bb,fac1,fac2,wleg,wdcor,vc,f,alpha,beta) &
!$OMP& PRIVATE(i,r,z,k,r1,z1,rrr,rr2,zd1,zd2,rz1,rz2,rrz1,rrz2, &
!$OMP&         xx1,j,xx,y1,s1,t,y2,vik,xxx)
     Do i=1,nghl
        r = fl(i); z = fh(i)
        Do k=1,i
           !
           r1 = fl(k); z1 = fh(k)
           rrr = two*r*r1; rr2 = (r - r1)**2
           ! z>0 part
           zd1 = (z - z1)**2
           rz1 = rr2 + zd1
           ! z<0 part
           zd2 = (z + z1)**2
           rz2 = rr2 + zd2
           ! Gauss-Legendre integration over u from 0 to D
           xx1=zero
           Do j=1,nleg
              xx=(one-xleg(j)**beta)**alpha ! change of variable to 0 <= u <= 1
              xxx=(one-xleg(j)**beta)**(alpha+one)
              y1=(xleg(j)/(bb*xx))**2 ! u^2
              s1=y1*rrr               ! 2 u^2 r r'
              y2=besei0(s1)           ! I0( 2 u^2 r r' ) * exp(-2 u^2 r r')
              xx1=xx1+fac2*wleg(j)*y2*(Exp(-rz1*y1) + fac1*Exp(-rz2*y1)) / xxx
           End Do
           vik=f*xx1/bb
           !
           vc(i,k)=vik*wdcor(k)  !wdcor=pi*wh*wl*bz*bp*bp
           vc(k,i)=vik*wdcor(i)  !wdcor=pi*wh*wl*bz*bp*bp
           !
        End Do  !k
     End Do  !i
!$OMP End Parallel Do
     !
  End If

  ! Calculation of the Coulomb field
  cou=zero
  Call dgemm('n','n',nghl,1,nghl,1.0_pr,vc,nghl,ro(:,2),nghl,0.0_pr,cou,nghl)
  !
  If (IDEBUG.Eq.1) Call get_CPU_time('coulom',1)
  !
End Subroutine coulom
!=======================================================================
!

!=======================================================================
!Subroutine coulom
!  !---------------------------------------------------------------------
!  ! Coulom-field (direct part), Gogny prescription
!  !---------------------------------------------------------------------
!  Use HFBTHO
!  Implicit None
!  Integer(ipr), Save :: i,j,k
!  Real(pr), Save :: zd2,y1,y2,xx1,s1,vik,f,r,r1,  &
!       rr2,z,z1,zd1,t,bb,r2,r12,rrr,rz1,rz2,rrz1,rrz2,xx
!  Real(pr), Save ::  &
!       bb1=3.51562290_pr,g1=0.398942280_pr,g5=0.009162810_pr,  &
!       bb2=3.08994240_pr,g2=0.013285920_pr,g6=0.020577060_pr,  &
!       bb3=1.20674920_pr,g3=0.002253190_pr,g7=0.026355370_pr,  &
!       bb4=0.26597320_pr,g4=0.001575650_pr,g8=0.016476330_pr,  &
!       bb5=0.03607680_pr,g9=0.003923770_pr,  &
!       bb6=0.00458130_pr,bbxx=3.750_pr
!  !
!  Call get_CPU_time('coulom',0)
!  !
!  If(icacou.Eq.0) Then
!     icacou=1
!     bb=Max(bz,bp)**4;  f=chargee2/Sqrt(pi) ! e^2/Sqrt(pi)
!     bb=50.0d0
!     Do i=1,nghl
!        r=fl(i); z=fh(i); r2=r*r
!        Do k=1,i
!           r1=fl(k); z1=fh(k); r12=r1*r1
!           !
!           rrr=two*r*r1;    rr2=(r-r1)**2
!           zd1=(z-z1)**2;   zd2=(z+z1)**2
!           rz1=r2+r12+zd1;  rz2=r2+r12+zd2
!           rrz1=rr2+zd1;    rrz2=rr2+zd2
!           !
!           xx1=zero
!           Do j=1,nleg
!              xx=Sqrt(one-xleg(j)**2)
!              y1=(xleg(j)/(bb*xx))**2
!              s1=y1*rrr
!              If(s1.Le.bbxx) Then
!                 t=(s1/bbxx)**2
!                 y2=one+t*(bb1+t*(bb2+t*(bb3+t*(bb4+t*(bb5+t*bb6)))))
!                 IF(Parity) THEN
!                    y2=y2*(Exp(-rz1*y1)+Exp(-rz2*y1))
!                 ELSE
!                    y2=y2* Exp(-rz1*y1)  !+Exp(-rz2*y1))
!                 END IF
!              Else
!                 t=(bbxx/s1)
!                 y2=g1+t*(g2+t*(g3+t*(-g4+t*(g5+t*(-g6+t*(g7+t*(-g8+t*g9)))))))
!                 IF(Parity) THEN
!                    y2=y2/Sqrt(s1)*(Exp(-rrz1*y1)+Exp(-rrz2*y1))
!                 ELSE
!                    y2=y2/Sqrt(s1)* Exp(-rrz1*y1) !+Exp(-rrz2*y1))
!                 END IF
!              End If
!              xx1=xx1+wleg(j)*y2/(bb*xx**3)
!           End Do
!           vik=f*xx1
!           vc(i,k)=vik*wdcor(k)  !wdcor=pi*wh*wl*bz*bp*bp/fd
!           vc(k,i)=vik*wdcor(i)
!        End Do  !k
!     End Do  !i
!  End If
!  ! calculation of the coulomb field
!  cou=zero
!  Do i=1,nghl
!     cou(:)=cou(:)+vc(:,i)*ro(i,2)
!  End Do
!  ! Call dgemm('n','n',nghl,1,nghl,1.0_pr,vc,nghl,ro(:,2),nghl,0.0_pr,cou,nghl)
!  !
!  Call get_CPU_time('coulom',1)
!  !
!End Subroutine coulom
!=======================================================================
!
!=======================================================================
Subroutine HartreeDir
  !---------------------------------------------------------------------
  ! Hartree-field (direct part)
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: i,j,k
  Real(pr) :: xx1,vik00,vik01,vik11
  Real(pr) :: r,rr,rrr,r1,r2,rr1,rr2
  Real(pr) :: z,z1,zdm,zdp,rzm,rzp
  Real(pr), Allocatable :: u(:)
  If(icahartree.Eq.0) Then
     icahartree=1
     !
     If(Allocated(u)) Deallocate(u); Allocate(u(nleg))
     u=Cos(HALF*Pi*xleg)
     Do i=1,nghl
        r=fl(i); z=fh(i); rr=r*r
        Do k=1,i
           r1=fl(k); z1=fh(k); rr1=r1*r1;     rr2=two*r*r1; rrr=rr+rr1;
           zdm=(z-z1)**2;      zdp=(z+z1)**2; rzm=rrr+zdm;  rzp=rrr+zdp
           vik00=0.250_pr*Sum(wleg*( &
                + HartreeV00(Sqrt(rzp-rr2*u)) &
                + HartreeV00(Sqrt(rzm-rr2*u)) &
                + HartreeV00(Sqrt(rzp+rr2*u)) &
                + HartreeV00(Sqrt(rzm+rr2*u))))
           vik01=0.250_pr*Sum(wleg*( &
                + HartreeV01(Sqrt(rzp-rr2*u)) &
                + HartreeV01(Sqrt(rzm-rr2*u)) &
                + HartreeV01(Sqrt(rzp+rr2*u)) &
                + HartreeV01(Sqrt(rzm+rr2*u))))
           vik11=0.250_pr*Sum(wleg*( &
                + HartreeV11(Sqrt(rzp-rr2*u)) &
                + HartreeV11(Sqrt(rzm-rr2*u)) &
                + HartreeV11(Sqrt(rzp+rr2*u)) &
                + HartreeV11(Sqrt(rzm+rr2*u))))

           vhart00(i,k)=vik00*wdcor(k)         ! wdcor=pi*wh*wl*bz*bp*bp/fd
           vhart00(k,i)=vik00*wdcor(i)
           vhart01(i,k)=vik01*wdcor(k)
           vhart01(k,i)=vik01*wdcor(i)
           vhart11(i,k)=vik11*wdcor(k)
           vhart11(k,i)=vik11*wdcor(i)
        End Do  !k
     End Do  !i
     Deallocate(u)
  End If
  ! calculation of the Hartree field
  vDHartree=0.0_pr
  Do i=1,nghl
     vDHartree(:,1)=vDHartree(:,1)+vhart00(:,i)*(ro(i,1)+ro(i,2))+vhart01(:,i)*(ro(i,1)-ro(i,2))
     vDHartree(:,2)=vDHartree(:,2)+vhart11(:,i)*(ro(i,1)-ro(i,2))+vhart01(:,i)*(ro(i,1)+ro(i,2))
  End Do
  !write(*,*) ' A=',sum((ro(:,1)+ro(:,2))*wdcor(:)),hatree_rc
  !write(*,*) ' EHartree=',sum(half*vDHartree(:)*(ro(:,1)+ro(:,2))*wdcor(:))
End Subroutine HartreeDir
!=======================================================================
!
!=======================================================================
Subroutine optHFBTHO
  !---------------------------------------------------------------------
  ! oprimization arrays
  ! NB FI2D_opt(JA,ihil) == Laplassisan(r,z) HOwf
  !    FID2D-xlamy2*FID  == Laplassisan(r,z,phy) FID
  !    FIU2D-xlapy2*FIU  == Laplassisan(r,z,phy) FIU
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: i,ih,il,ib,ibx,nd,nza,nra,nla,nsa
  Integer(ipr) :: ihil,laplus,im,JA,N1,N2,ndnd,n12,n21
  Real(pr)    :: qla,v2,v4,yi,y,y2,qha,qhla,xmi,u,u2,un,up,xxx
  Real(pr)    :: sml2,cnzaa,cnraa,a,b
  Real(pr)    :: FITW1,FITW2,FITW3,FITW4
  Real(pr)    :: fi1r,fi1z,fi2d,QHL1A,QH1LA,vh,vdh,vsh,hbh
  Real(pr)    :: SRFIh,SFIRh,SFIZh,SZFIh,SNABLARh,SNABLAZh
  Real(pr)    :: xlam,xlam2,xlamy,xlamy2,xlap,xlap2,xlapy,xlapy2,XLAMPY
  Real(pr)    :: bpi,bpi2,bzi,bzi2,xh2
  !
  bpi=one/bp; bpi2=bpi*bpi; bzi=one/bz; bzi2=bzi*bzi
  !
  !-----------------------------------------
  ! Allocate the optimization arrays
  !-----------------------------------------
  If(Allocated(QHLA_opt)) Deallocate(QHLA_opt,FI1R_opt,FI1Z_opt,FI2D_opt,y_opt)
  Allocate(QHLA_opt(ntx,nghl),FI1R_opt(ntx,nghl),FI1Z_opt(ntx,nghl),FI2D_opt(ntx,nghl),y_opt(nghl))  
  !----------------------------------------------
  ! START BLOCKS
  !----------------------------------------------
  Do ib=1,NB
     ND=ID(ib); IM=ia(ib)
     If(Parity) Then
        LAPLUS=(ib+1)/2 !Yesp
     Else
        LAPLUS=ib       !Nop
     Endif
     XLAP=LAPLUS; XLAM=XLAP-ONE; xlap2=xlap*xlap; xlam2=xlam*xlam
     !----------------------------------------------
     ! SUM OVER GAUSS INTEGRATION POINTS
     !----------------------------------------------    
     Do IL=1,ngl
        v2=half/xl(il); v4=v2*v2 
        Do IH=1,ngh
           ihil=ih+(il-1)*ngh; xh2=xh(ih)**2
           If(iLST1.Eq.0) Then
              ! HO-basis
              yi=Sqrt(xl(il))*bp; y=one/yi; y2=y*y
              xlamy=xlam*y; xlamy2=xlam2*y2; xlapy=xlap*y; xlapy2=xlap2*y2; XLAMPY=XLAMY+XLAPY
#ifndef hide_tho 
           Else
              ! THO-basis
              y=fli(ihil); y2=y*y; xlamy=xlam*y; u=xh(ih); u2=u*u; 
              xlamy2=xlam2*y2; xlapy=xlap*y; xlapy2=xlap2*y2; XLAMPY=XLAMY+XLAPY
#endif    
           End If
           y_opt(ihil)=y
           !----------------------------------------------
           ! SCAN OVER BASIS STATES
           !----------------------------------------------
           Do N1=1,ND
              JA=N1+IM; NLA=NL(JA); NRA=NR(JA); NZA=NZ(JA); NSA=NS(JA)
              SML2=NLA*NLA; CNZAA=NZA+NZA+1; CNRAA=NRA+NRA+NLA+1
              QHA=QH(NZA,IH); QLA=QL(NRA,NLA,IL); QHLA=QHA*QLA
              QHL1A=QHA*QL1(NRA,NLA,IL)*V2; QH1LA=QH1(NZA,IH)*QLA   
              If(iLST1.Eq.0) Then
                 ! HO-basis
                 FI1R=(two*Sqrt(xl(il))*bpi)*QHL1A
                 FI1Z=bzi*QH1LA
                 FI2D=((xh2-CNZAA)*bzi2+four*(p14-CNRAA*V2+SML2*V4)*xl(il)*bpi2 )*QHLA 
#ifndef hide_tho 
              Else 
                 ! THO-basis
                 u=xh(ih); u2=u*u;
                 FI1R=FP4(IHIL)*QHLA+FP5(IHIL)*QH1LA+FP6(IHIL)*QHL1A
                 FI1Z=FP1(IHIL)*QHLA+FP2(IHIL)*QH1LA+FP3(IHIL)*QHL1A                
                 FI2D=(FS1(IHIL)*QH1LA*QH1LA+FS2(IHIL)*QHL1A*QHL1A  & 
                      +FOUR*FS4(IHIL)*QH1LA*QHL1A  &
                      +TWO*(FS5(IHIL)*QH1LA+FS6(IHIL)*QHL1A)*QHLA  &
                      +((U2-CNZAA)*FS1(IHIL)+(p14-CNRAA*V2+SML2*V4)*FS2(IHIL)  &
                      +FS3(IHIL))*QHLA*QHLA-TWO*(FI1R*FI1R+FI1Z*FI1Z))/(TWO*QHLA)
#endif               
              Endif
              QHLA_opt(JA,ihil)=QHLA; FI2D_opt(JA,ihil)=FI2D; FI1R_opt(JA,ihil)=FI1R; FI1Z_opt(JA,ihil)=FI1Z                            
           End Do !N1
           !
        End Do !IH
     End Do !IL
  End Do !IB  
End Subroutine optHFBTHO
!=======================================================================
!
!=========================================================================
Subroutine DENSIT
  !---------------------------------------------------------------------
  ! local densities in coordinate space
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: iw,nsa,nza,nra,nla,k,i,nd,il,ih,ihil,laplus
  Integer(ipr) :: imen,ib,im,it,J,JJ,JA,JN,k0,k1,k2,ibiblo
  Real(pr)     :: s,ss,sd,yi,y,y2,sml2,cnzaa,cnraa,u,u2,v2,v4
  Real(pr)     :: anik,pnik,qhla,qh1la,qhl1a,qla,qha,fi1r,fi1z,fi2d,fidd
  Real(pr)     :: xlam,xlam2,xlamy,xlamy2,xlap,xlap2,xlapy,xlapy2,XLAMPY
  Real(pr)     :: TFIU,TFID,TFIUR,TFIDR,TFIUZ,TFIDZ,TPFIU,TPFID,TFIUD2,TFIDD2
  Real(pr)     :: PIU,PIUZ,PIUR,PIUD2,PID,PIDZ,PIDR,PIDD2
  Real(pr)     :: TW1,TW2,TW3,TW4,TW5,TW6,TW7,TW8,TW9,TW10,TW11,TW_T,WGT(nghl)
  Real(pr)     :: Takaihil,Troihil,Tdjihil,Ttauihil,Tdroihil,TSRFIihil
  Real(pr)     :: TSFIRihil,TSFIZihil,TSZFIihil,TNABLARIHIL,TNABLAZIHIL
  Real(pr), Pointer :: TAKA(:),TRO(:),TDJ(:),TTAU(:),TDRO(:)
  Real(pr), Pointer :: TSRFI(:),TSFIR(:),TSFIZ(:),TSZFI(:),TNABLAR(:),TNABLAZ(:)
  Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
  Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
  !
  Call get_CPU_time('densit',0)
  !
  !-----------------------------------------------
  ! ZERO N & P DENSITIES
  !-----------------------------------------------
  RO=ZERO; TAU=ZERO; DJ=ZERO; DRO=ZERO; AKA=ZERO; SZFI=ZERO; SFIZ=ZERO
  SRFI=ZERO; SFIR=ZERO; NABLAR=ZERO; NABLAZ=ZERO; VARMAS=ZERO
  Do it=1,2 
     !
     ! blocking
     ibiblo=bloblo(keyblo(it),it)  
     ! 
     If(it.Eq.1) Then
        EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN
     Else
        EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP
     Endif
     !
     TRO=>ro(:,it);         TTAU=>tau(:,it);       TDJ=>dj(:,it);     TDRO=>dro(:,it)
     TSZFI=>SZFI(:,it);     TSFIZ=>SFIZ(:,it);     TSRFI=>SRFI(:,it); TSFIR=>SFIR(:,it)
     TNABLAR=>NABLAR(:,it); TNABLAZ=>NABLAZ(:,it); TAKA=>aka(:,it)
     !
     ! case of zero particle number, only flush densities
     If((npr_INI(1).Eq.0).And.(it.Eq.1)) Cycle
     If((npr_INI(2).Eq.0).And.(it.Eq.2)) Cycle
     !-----------------------------------------------
     ! SCAN OVER BLOCKS
     !-----------------------------------------------
     Do ib=1,NB
        ND=ID(ib); IM=ia(ib)
        If(Parity) Then
           LAPLUS=(ib+1)/2 !Yesp
        Else
           LAPLUS=ib       !Nop
        Endif
        XLAP=LAPLUS; XLAM=XLAP-ONE; xlap2=xlap*xlap; xlam2=xlam*xlam
        !
        ! blocking
        ibiblo=bloblo(keyblo(it),it) 
        K0=0; If(ibiblo.Eq.ib) K0=blo123d(it)
        !write(*,*) 'k0=',k0,ib,ibiblo
        !
        !----------------------------------------------
        ! PAIRING WINDOW QP WAVE FUNCTIONS
        !----------------------------------------------
        k1=ka(ib,it)+1; k2=ka(ib,it)+kd(ib,it); imen=k2-k1+1
        If(IMEN.Gt.0) Then
           J=0
           Do JJ=1,nd
              Do K=K1,K2
                 J=J+1; I=KpwiPo(K)+JJ; an(J)=VqpPo(I); ank(J)=UqpPo(I)
              End Do
           End Do
           !-----------------------------------------------
           ! SCAN OVER GAUSS INTEGRATION POINTS
           !-----------------------------------------------
           Do ihil=1,nghl
              y=y_opt(ihil); xlamy=xlam*y;     xlapy=xlap*y;   XLAMPY=XLAMY+XLAPY
              y2=y*y;        xlamy2=xlam2*y2;  xlapy2=xlap2*y2                  
              Do K=1,IMEN
                 FIU(K)=ZERO; FIUZ(K)=ZERO; FIUR(K)=ZERO; FIUD2N(K)=ZERO; PFIU(K)=ZERO 
                 FID(K)=ZERO; FIDZ(K)=ZERO; FIDR(K)=ZERO; FIDD2N(K)=ZERO; PFID(K)=ZERO
              End Do
              If(K0.Ne.0) Then
                 PIU=ZERO;  PIUZ=ZERO; PIUR=ZERO; PIUD2=ZERO; PID=ZERO;  PIDZ=ZERO; PIDR=ZERO; PIDD2=ZERO
              Endif
              !-----------------------------------------------
              ! SUM OVER BASIS STATES
              !-----------------------------------------------
              Do I=1,ND
                 JA=IM+I; NSA=NS(JA); JN=(I-1)*imen 
                 QHLA=QHLA_opt(JA,ihil); FI2D=FI2D_opt(JA,ihil); FI1Z=FI1Z_opt(JA,ihil); FI1R=FI1R_opt(JA,ihil)
                 !-----------------------------------------------
                 ! QUASIPARTICLE WF IN COORDINATE SPACE
                 !-----------------------------------------------
                 If (NSA.Gt.0) Then 
                    ! SPIN Up
                    Call DAXPY(IMEN,-QHLA,ANK(JN+1),1,PFIU,1)
                    Call DAXPY(IMEN,QHLA ,AN(JN+1) ,1,FIU,1) 
                    Call DAXPY(IMEN,FI2D ,AN(JN+1) ,1,FIUD2N,1)
                    Call DAXPY(IMEN,FI1R ,AN(JN+1) ,1,FIUR,1)  
                    Call DAXPY(IMEN,FI1Z ,AN(JN+1) ,1,FIUZ,1)  
                    ! blocking
                    If(K0.Ne.0) Then
                       PNIK=ANK(JN+K0)
                       PIU=PIU+PNIK*QHLA;   PIUD2=PIUD2+PNIK*FIDD; 
                       PIUR=PIUR+PNIK*FI1R; PIUZ=PIUZ+PNIK*FI1Z
                    Endif
                 Else
                    ! SPIN Do
                    Call DAXPY(IMEN,-QHLA,ANK(JN+1),1,PFID,1)
                    Call DAXPY(IMEN,QHLA ,AN(JN+1) ,1,FID,1) 
                    Call DAXPY(IMEN,FI2D ,AN(JN+1) ,1,FIDD2N,1)
                    Call DAXPY(IMEN,FI1R ,AN(JN+1) ,1,FIDR,1)  
                    Call DAXPY(IMEN,FI1Z ,AN(JN+1) ,1,FIDZ,1)  
                    ! blocking
                    If(K0.Ne.0) Then
                       PNIK=ANK(JN+K0)
                       PID=PID+PNIK*QHLA;   PIDD2=PIDD2+PNIK*FIDD; 
                       PIDR=PIDR+PNIK*FI1R; PIDZ=PIDZ+PNIK*FI1Z
                    Endif
                 End If
              End Do ! I
              !-----------------------------------------------
              ! DENSITIES IN COORDINATE SPACE
              !-----------------------------------------------
              Takaihil=zero;    Troihil=zero;    Tdjihil=zero;   Ttauihil=zero;  Tdroihil=zero                 
              TSRFIihil=zero;   TSFIRihil=zero;  TSFIZihil=zero; TSZFIihil=zero; TNABLARIHIL=zero; TNABLAZIHIL=zero; 
              Do K=1,IMEN
                 TFIU=FIU(K); TFIUZ=FIUZ(K); TFIUR=FIUR(K); TFIUD2=FIUD2N(K); TPFIU=PFIU(K)
                 TFID=FID(K); TFIDZ=FIDZ(K); TFIDR=FIDR(K); TFIDD2=FIDD2N(K); TPFID=PFID(K) 
                 !
                 TW1=TPFIU*TFIU+TPFID*TFID;           TAKAIHIL=TAKAIHIL+TW1
                 TW2=TFIU*TFIU+TFID*TFID;             TROIHIL=TROIHIL+TW2
				 TW_T=TFIUR*TFIUR+TFIDR*TFIDR+TFIUZ*TFIUZ+TFIDZ*TFIDZ
                 TW3=TFIUR*TFIDZ-TFIDR*TFIUZ+XLAMY*TFIU*(TFIUR-TFIDZ)-XLAPY*TFID*(TFIDR+TFIUZ);         TDJIHIL=TDJIHIL+TW3
                 TW4=TW_T+XLAMY2*TFIU*TFIU+XLAPY2*TFID*TFID; TTAUIHIL=TTAUIHIL+TW4 
                 TW5=TW_T+TFIU*TFIUD2+TFID*TFIDD2;    TDROIHIL=TDROIHIL+TW5
                 TW6=TFIUR*TFID-TFIDR*TFIU;           TSRFIIHIL=TSRFIIHIL+TW6
                 TW7=TFIU*TFID*XLAMPY;                TSFIRIHIL=TSFIRIHIL+TW7
                 TW8=XLAMY*TFIU*TFIU-XLAPY*TFID*TFID; TSFIZIHIL=TSFIZIHIL+TW8
                 TW9=TFIUZ*TFID-TFIDZ*TFIU;           TSZFIIHIL=TSZFIIHIL+TW9
                 TW10=TFIUR*TFIU+TFIDR*TFID;          TNABLARIHIL=TNABLARIHIL+TW10
                 TW11=TFIUZ*TFIU+TFIDZ*TFID;          TNABLAZIHIL=TNABLAZIHIL+TW11
                 !
                 If(K.Ne.K0) Cycle
                 !
                 ! blocking
                 TAKAIHIL=TAKAIHIL-TW1;                   TW1=PIU*PIU+PID*PID                  
                 TROIHIL=TROIHIL-HALF*(TW2-TW1);  TW2=PIUR*PIDZ-PIDR*PIUZ+XLAMY*PIU*(PIUR-PIDZ)-XLAPY*PID*(PIDR+PIUZ)            
                 TDJIHIL=TDJIHIL-HALF*(TW3-TW2);  TW3=PIUR*PIUR+PIDR*PIDR+PIUZ*PIUZ+PIDZ*PIDZ+XLAMY2*PIU*PIU+XLAPY2*PID*PID     
                 TTAUIHIL=TTAUIHIL-HALF*(TW4-TW3);        TW4=PIU*PIUD2+PID*PIDD2;              
                 TDROIHIL=TDROIHIL-HALF*(TW5-TW4);        TW5=PIUR*PID-PIDR*PIU;                
                 TSRFIIHIL=TSRFIIHIL-HALF*(TW6-TW5);      TW6=PIU*PID*XLAMPY;                   
                 TSFIRIHIL=TSFIRIHIL-HALF*(TW7-TW6);      TW7=XLAMY*PIU*PIU-XLAPY*PID*PID;      
                 TSFIZIHIL=TSFIZIHIL-HALF*(TW8-TW7);      TW8=PIUZ*PID-PIDZ*PIU;                
                 TSZFIIHIL=TSZFIIHIL-HALF*(TW9-TW8);      TW9=PIUR*PIU+PIDR*PID;                
                 TNABLARIHIL=TNABLARIHIL-HALF*(TW10-TW9); TW10=PIUZ*PIU+PIDZ*PID;               
                 TNABLAZIHIL=TNABLAZIHIL-HALF*(TW11-TW10)
              End Do !K
              Taka(ihil)    = Taka(ihil)    + Takaihil;    Tro(ihil)     = Tro(ihil)     + Troihil      
              Tdj(ihil)     = Tdj(ihil)     + Tdjihil;     Ttau(ihil)    = Ttau(ihil)    + Ttauihil
              Tdro(ihil)    = Tdro(ihil)    + Tdroihil                 
              TSRFI(ihil)   = TSRFI(ihil)   + TSRFIihil;   TSFIR(ihil)   = TSFIR(ihil)   + TSFIRihil
              TSFIZ(ihil)   = TSFIZ(ihil)   + TSFIZihil;   TSZFI(ihil)   = TSZFI(ihil)   + TSZFIihil
              TNABLAR(IHIL) = TNABLAR(IHIL) + TNABLARIHIL; TNABLAZ(IHIL) = TNABLAZ(IHIL) + TNABLAZIHIL
           End Do !ihil
        End If
     End Do !IB
     s=two*Sum(tro); sd=four*Sum(tdro); drhoi(it)=sd; Sumnz(it)=Abs(s-Dble(npr(it))) 
     varmas=varmas+s; DNFactor(it)=Dble(npr(it))/s; DNFactor(3)=DNFactor(1)+DNFactor(2)
     !----------------------------------------------------
     ! REMOVES INT.WEIGHTS AND MULTIPLIES BY THE JACOBIAN
     !----------------------------------------------------
     piu=two*Dble(npr(it))/s     
     WGT=wdcori 
     Call dscal(NGHL,piu,WGT,1)
     Tro=Tro*WGT; Ttau=Ttau*WGT; Taka=Half*Taka*WGT;
     TSRFI=TSRFI*WGT; TSFIR=TSFIR*WGT; 
     TSFIZ=TSFIZ*WGT; TSZFI=TSZFI*WGT; 
     Call dscal(NGHL,two,WGT,1)
     Tdro=Tdro*WGT; Tdj=Tdj*WGT
     TNABLAR=TNABLAR*WGT; TNABLAZ=TNABLAZ*WGT
  End Do !it
  !
  Call get_CPU_time('densit',1)
  !
  !----------------------------------------------------
  ! COULOMB AND HATREE FIELDS
  !----------------------------------------------------
  Call coulom
  !Call HartreeDir
End Subroutine DENSIT
!=======================================================================
!
!=======================================================================
Subroutine field
  !---------------------------------------------------------------------
  ! calculates fields in r-space form axially symmetric densities
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: iw,it,ita,ihli
  Real(pr) :: ra,ra2,rs,rsa,rsa0
  Real(pr) :: rt,rt1,ds,da,dt,dt1,tts,tta,tt,tt1,djs,dja,djt,djt1
  Real(pr) :: rsa0A,rsa0A1,V0V1,v01a,rns,rps,rsa1,rsa12,rsa10
  Real(pr) :: rsa0An,rsa0An1,rsa0As,rsa0As1
  Real(pr) :: RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,DJ_0,DJ_1
  Real(pr) :: SZFIN,SFIZN,SRFIN,SFIRN,SZFIP,SFIZP,SRFIP,SFIRP
  Real(pr) :: SZFI_0,SFIZ_0,SRFI_0,SFIR_0,SZFI_1,SFIZ_1,SRFI_1,SFIR_1
  Real(pr) :: SNABLARN,SNABLAZN,SNABLARP,SNABLAZP
  Real(pr) :: SNABLAR_0,SNABLAZ_0,SNABLAR_1,SNABLAZ_1
  Real(pr) :: J2_0,J2_1
  Real(pr) :: cx,x
  Real(pr),Dimension(2) :: pUr,pUt,pUNr,pUNz,pUDr,pUDj,pUFIZ,pUZFI,pUFIR,pURFI
  Real(pr),Dimension(2) :: tUr,tUt,tUNr,tUNz,tUDr,tUDj,tUFIZ,tUZFI,tUFIR,tURFI
  Real(pr), Save :: ALAMBDA=0.0_pr,AEPSI=1.0_pr,CSPR=1.0_pr
  !
  Call get_CPU_time('field',0)
  !
  cx=0.0_pr
  If (icstr.Ne.0) Then
     !
     ! constraints
     !
     x=cqad*ty20*(cdef-bet)
     !
     If (icstr.Lt.0) Then
        !
        ! Agumented Lagrangian constraint
        !
        If(Abs(si).Lt.AEPSI.And.Abs(cdef-bet).Gt.0.00010_pr) Then
           If(Abs(cdef-bet).Lt.0.30_pr) Then
              AEPSI=Abs(si)
              ALAMBDA=ALAMBDA+CSPR*x
              Write(*,'(5(a,f15.8))') '    UPDATING   ALAMBDA=',ALAMBDA,'  & AEPSI=',AEPSI
           Else
              CSPR=CSPR*2.0_pr   
              Write(*,'(5(a,f15.8))') '    UPDATING CSPR=',CSPR
           Endif
        Endif
        !
        cx=ALAMBDA+CSPR*x
        !
     Else If (icstr.Gt.0) Then
        !
        ! Quadratic constraint
        !
        cx=x
        !
     End If
     Do iw=lout,lfile
        Write(iw,'(40x,a,f7.3,a,f8.3,2x,f15.8)') '  Requested cdef=',cdef,' Achieved=',bet,q2(3)/100.0_pr
     Enddo
  End If
  !
  ! fields
  Do ihli=1,nghl
     !
     RHO_0 =ro(ihli,1)+ro(ihli,2)     ; RHO_1 =ro(ihli,1)-ro(ihli,2)
     TAU_0 =tau(ihli,1)+tau(ihli,2)   ; TAU_1 =tau(ihli,1)-tau(ihli,2)
     DRHO_0=dro(ihli,1)+dro(ihli,2)   ; DRHO_1=dro(ihli,1)-dro(ihli,2)
     DJ_0  =dj(ihli,1)+dj(ihli,2)     ; DJ_1  =dj(ihli,1)-dj(ihli,2)
     SFIZ_0=SFIZ(ihli,1)+SFIZ(ihli,2) ; SFIZ_1=SFIZ(ihli,1)-SFIZ(ihli,2)
     SFIR_0=SFIR(ihli,1)+SFIR(ihli,2) ; SFIR_1=SFIR(ihli,1)-SFIR(ihli,2)
     SZFI_0=SZFI(ihli,1)+SZFI(ihli,2) ; SZFI_1=SZFI(ihli,1)-SZFI(ihli,2)
     SRFI_0=SRFI(ihli,1)+SRFI(ihli,2) ; SRFI_1=SRFI(ihli,1)-SRFI(ihli,2)
     SNABLAR_0=NABLAR(ihli,1)+NABLAR(ihli,2)
     SNABLAR_1=NABLAR(ihli,1)-NABLAR(ihli,2)
     SNABLAZ_0=NABLAZ(ihli,1)+NABLAZ(ihli,2)
     SNABLAZ_1=NABLAZ(ihli,1)-NABLAZ(ihli,2)
     !
     J2_0=SFIZ_0**2+SFIR_0**2+SZFI_0**2+SRFI_0**2
     J2_1=SFIZ_1**2+SFIR_1**2+SZFI_1**2+SRFI_1**2    
     !
     tUr=zero ; tUDr=zero ; tUNr=zero ; tUNz=zero 
     tUt=zero ; tUDj=zero ; tUFIZ=zero ; tUZFI=zero
     tUFIR=zero ; tURFI=zero ; 
     !
     Call calculate_U_parameters(RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1, &
          (SNABLAR_0**2+SNABLAZ_0**2),(SNABLAR_1**2+SNABLAZ_1**2) )
     !
     ! FUNCTIONAL       
     ! E=E+(hb0*(TAU_0+TAU_1)*HALF+hb0*(TAU_0-TAU_1)*HALF  &                         ! tau
     !+Urhotau(0,0)*RHO_0*TAU_0+Urhotau(1,0)*RHO_1*TAU_1  &                         ! rho tau
     !+Urhotau(2,0)*RHO_0*TAU_1+Urhotau(3,0)*RHO_1*TAU_0  &
     !+Urhorho(0,0)*RHO_0**2+Urhorho(1,0)*RHO_1**2  &                               ! rho^2
     !+(Urhorho(2,0)+Urhorho(3,0))*RHO_0*RHO_1  &
     !+UrhoDrho(0,0)*RHO_0*DRHO_0+UrhoDrho(1,0)*RHO_1*DRHO_1  &                     ! rho Delta rho
     !+UrhoDrho(2,0)*RHO_0*DRHO_1+UrhoDrho(3,0)*RHO_1*DRHO_0  &
     !+Unablarho(0,0)*(SNABLAR_0*SNABLAR_0+SNABLAZ_0*SNABLAZ_0)  &                  ! (nabla rho)^2
     !+Unablarho(1,0)*(SNABLAR_1*SNABLAR_1+SNABLAZ_1*SNABLAZ_1)  &
     !+(Unablarho(3,0)+Unablarho(2,0))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)  &
     !+UrhonablaJ(0,0)*RHO_0*DJ_0+UrhonablaJ(1,0)*RHO_1*DJ_1  &                     ! rho nabla J 
     !+UrhonablaJ(2,0)*RHO_0*DJ_1+UrhonablaJ(3,0)*RHO_1*DJ_0  &
     !+UJnablarho(0,0)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))  &     ! J nabla rho
     !+UJnablarho(1,0)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))  &
     !+UJnablarho(2,0)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))  &
     !+UJnablarho(3,0)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))  &
     !+UJJ(0,0)*J2_0+UJJ(1,0)*J2_1  &                                               ! JJ 
     !+(UJJ(3,0)+UJJ(2,0))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1) 
     !
     ! tUr(1)=dE/d RHO_0;       tUr(2)=dE/d RHO_1
     ! tUt(1)=dE/d TAU_0;       tUt(2)=dE/d TAU_1
     ! tUDr(1)=dE/d DeltaRHO_0; tUDr(2)=dE/d DeltaRHO_1
     ! and so on ...
     !    
     !TEST 
     !write(*,'(4(2x,g26.10))') UrhoDrho(0,0)-CrDr(0),CrDr(0),UrhoDrho(1,1)-CrDr(1),CrDr(1); pause                     
     ! Contributions in the case 'u' depends on RHO_0     
     tUr(1)=tUr(1)+two*Urhorho(0,0)*RHO_0+Urhorho(0,1)*RHO_0*RHO_0+Urhorho(1,1)*RHO_1*RHO_1  &  !! rho^2
          +(Urhorho(3,0)+Urhorho(2,0))*RHO_1+(Urhorho(3,1)+Urhorho(2,1))*RHO_0*RHO_1  
     tUr(2)=tUr(2)+two*Urhorho(1,0)*RHO_1+Urhorho(0,2)*RHO_0*RHO_0+Urhorho(1,2)*RHO_1*RHO_1  &
          +(Urhorho(3,0)+Urhorho(2,0))*RHO_0+(Urhorho(3,2)+Urhorho(2,2))*RHO_0*RHO_1  
     tUr(1)=tUr(1)+vDHartree(ihli,1)
     tUr(2)=tUr(2)+vDHartree(ihli,2)		
     !
     tUr(1)=tUr(1)+Urhotau(0,0)*TAU_0+Urhotau(0,1)*TAU_0*RHO_0+Urhotau(1,1)*TAU_1*RHO_1  &  !! rho tau
          +Urhotau(2,0)*TAU_1+Urhotau(2,1)*RHO_0*TAU_1+Urhotau(3,1)*RHO_1*TAU_0
     tUt(1)=tUt(1)+Urhotau(0,0)*RHO_0+Urhotau(3,0)*RHO_1
     tUr(2)=tUr(2)+Urhotau(1,0)*TAU_1+Urhotau(1,2)*TAU_1*RHO_1+Urhotau(0,2)*TAU_0*RHO_0  &
          +Urhotau(3,0)*TAU_0+Urhotau(3,2)*RHO_1*TAU_0+Urhotau(2,2)*RHO_0*TAU_1
     tUt(2)=tUt(2)+Urhotau(1,0)*RHO_1+Urhotau(2,0)*RHO_0
     !
     tUr(1)=tUr(1)+UrhoDrho(0,0)*DRHO_0+UrhoDrho(0,1)*RHO_0*DRHO_0+UrhoDrho(1,1)*RHO_1*DRHO_1  &  !! rho Delta rho
          +UrhoDrho(2,0)*DRHO_1+UrhoDrho(2,1)*RHO_0*DRHO_1+UrhoDrho(3,1)*RHO_1*DRHO_0
     tUDr(1)=tUDr(1)+UrhoDrho(0,0)*RHO_0+UrhoDrho(3,0)*RHO_1
     tUr(2)=tUr(2)+UrhoDrho(1,0)*DRHO_1+UrhoDrho(1,2)*RHO_1*DRHO_1+UrhoDrho(0,2)*RHO_0*DRHO_0  &
          +UrhoDrho(3,0)*DRHO_0+UrhoDrho(3,2)*RHO_1*DRHO_0+UrhoDrho(2,2)*RHO_0*DRHO_1
     tUDr(2)=tUDr(2)+UrhoDrho(1,0)*RHO_1+UrhoDrho(2,0)*RHO_0
     !
     tUr(1)=tUr(1)+Unablarho(0,1)*(SNABLAR_0**2+SNABLAZ_0**2)+Unablarho(1,1)*(SNABLAR_1**2+SNABLAZ_1**2)  &  !! (nabla rho)^2
          +(Unablarho(2,1)+Unablarho(3,1))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)
     tUNr(1)=tUNr(1)+two*Unablarho(0,0)*SNABLAR_0+(Unablarho(2,0)+Unablarho(3,0))*SNABLAR_1
     tUNz(1)=tUNz(1)+two*Unablarho(0,0)*SNABLAZ_0+(Unablarho(2,0)+Unablarho(3,0))*SNABLAZ_1
     tUr(2)=tUr(2)+Unablarho(0,2)*(SNABLAR_0**2+SNABLAZ_0**2)+Unablarho(1,2)*(SNABLAR_1**2  &
          +SNABLAZ_1**2)+(Unablarho(2,2)+Unablarho(3,2))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)
     tUNr(2)=tUNr(2)+two*Unablarho(1,0)*SNABLAR_1+(Unablarho(2,0)+Unablarho(3,0))*SNABLAR_0
     tUNz(2)=tUNz(2)+two*Unablarho(1,0)*SNABLAZ_1+(Unablarho(2,0)+Unablarho(3,0))*SNABLAZ_0
     !
     tUr(1)=tUr(1)+UrhonablaJ(0,0)*DJ_0+UrhonablaJ(0,1)*DJ_0*RHO_0+UrhonablaJ(1,1)*DJ_1*RHO_1  &  !! rho nabla J
          +UrhonablaJ(2,0)*DJ_1+UrhonablaJ(2,1)*RHO_0*DJ_1+UrhonablaJ(3,1)*RHO_1*DJ_0
     tUDj(1)=tUDj(1)+UrhonablaJ(0,0)*RHO_0+UrhonablaJ(3,0)*RHO_1
     tUr(2)=tUr(2)+UrhonablaJ(1,0)*DJ_1+UrhonablaJ(1,2)*DJ_1*RHO_1+UrhonablaJ(0,2)*DJ_0*RHO_0  &
          +UrhonablaJ(3,0)*DJ_0+UrhonablaJ(3,2)*RHO_1*DJ_0+UrhonablaJ(2,2)*RHO_0*DJ_1
     tUDj(2)=tUDj(2)+UrhonablaJ(1,0)*RHO_1+UrhonablaJ(2,0)*RHO_0
     !
     tUr(1)=tUr(1)+UJnablarho(0,1)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))      !! J nabla rho
     tUr(1)=tUr(1)+UJnablarho(1,1)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1)) 
     tUr(1)=tUr(1)+UJnablarho(2,1)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0)) 
     tUr(1)=tUr(1)+UJnablarho(3,1)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))

     tUr(2)=tUr(2)+UJnablarho(0,2)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))      
     tUr(2)=tUr(2)+UJnablarho(1,2)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1)) 
     tUr(2)=tUr(2)+UJnablarho(2,2)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0)) 
     tUr(2)=tUr(2)+UJnablarho(3,2)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))

     tUNr(1)=tUNr(1)+UJnablarho(0,0)*(SFIZ_0-SZFI_0)                                                 
     tUNr(2)=tUNr(2)+UJnablarho(1,0)*(SFIZ_1-SZFI_1)
     tUNz(1)=tUNz(1) -UJnablarho(0,0)*(SFIR_0-SRFI_0)
     tUNz(2)=tUNz(2) -UJnablarho(1,0)*(SFIR_1-SRFI_1)

     tUFIZ(1)=tUFIZ(1)+UJnablarho(0,0)*SNABLAR_0*half
     tUFIZ(2)=tUFIZ(2)+UJnablarho(1,0)*SNABLAR_1*half
     tUZFI(1)=tUZFI(1)-UJnablarho(0,0)*SNABLAR_0*half
     tUZFI(2)=tUZFI(2)-UJnablarho(1,0)*SNABLAR_1*half
     tURFI(1)=tURFI(1)+UJnablarho(0,0)*SNABLAZ_0*half
     tURFI(2)=tURFI(2)+UJnablarho(1,0)*SNABLAZ_1*half       
     tUFIR(1)=tUFIR(1)-UJnablarho(0,0)*SNABLAZ_0*half       
     tUFIR(2)=tUFIR(2)-UJnablarho(1,0)*SNABLAZ_1*half
     !
     !! J.J (Mario: not tested for N2LO)
     tUr(1)=tUr(1)+UJJ(0,1)*J2_0+UJJ(1,1)*J2_1  &
          +(UJJ(3,1)+UJJ(2,1))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
     tUr(2)=tUr(2)+UJJ(0,2)*J2_0+UJJ(1,2)*J2_1  &
          +(UJJ(3,2)+UJJ(2,2))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
     tUFIZ(1)=tUFIZ(1)+UJJ(0,0)*SFIZ_0+half*(UJJ(3,0)+UJJ(2,0))*SFIZ_1
     tUFIR(1)=tUFIR(1)+UJJ(0,0)*SFIR_0+half*(UJJ(3,0)+UJJ(2,0))*SFIR_1
     tUZFI(1)=tUZFI(1)+UJJ(0,0)*SZFI_0+half*(UJJ(3,0)+UJJ(2,0))*SZFI_1
     tURFI(1)=tURFI(1)+UJJ(0,0)*SRFI_0+half*(UJJ(3,0)+UJJ(2,0))*SRFI_1
     tUFIZ(2)=tUFIZ(2)+UJJ(1,0)*SFIZ_1+half*(UJJ(3,0)+UJJ(2,0))*SFIZ_0
     tUFIR(2)=tUFIR(2)+UJJ(1,0)*SFIR_1+half*(UJJ(3,0)+UJJ(2,0))*SFIR_0
     tUZFI(2)=tUZFI(2)+UJJ(1,0)*SZFI_1+half*(UJJ(3,0)+UJJ(2,0))*SZFI_0
     tURFI(2)=tURFI(2)+UJJ(1,0)*SRFI_1+half*(UJJ(3,0)+UJJ(2,0))*SRFI_0
     !
     tUr(1)=tUr(1)+UFnonstdr(0)                                        !! other amplitudes
     tUr(2)=tUr(2)+UFnonstdr(1)
     !
     !!  External Filed
     !!
     !!tUr(1)=tUr(1)+Vexternal(0,zero,fl(ihli),fh(ihli))
     !!tUr(2)=tUr(2)+Vexternal(1,zero,fl(ihli),fh(ihli))
     !     
     ! Contributions in the case 'u' depends on TAU_0     
     !
     tUt(1)=tUt(1)+Urhotau(0,6)*RHO_0*TAU_0  &  
          +Urhotau(1,6)*RHO_1*TAU_1+Urhotau(2,6)*RHO_0*TAU_1  &
          +Urhotau(3,6)*RHO_1*TAU_0+Urhorho(0,6)*RHO_0**2  &
          +Urhorho(1,6)*RHO_1**2+(Urhorho(2,6)+Urhorho(3,6))*RHO_0*RHO_1  &
          +UrhoDrho(0,6)*RHO_0*DRHO_0+UrhoDrho(1,6)*RHO_1*DRHO_1  &
          +UrhoDrho(2,6)*RHO_0*DRHO_1+UrhoDrho(3,6)*RHO_1*DRHO_0  &
          +Unablarho(0,6)*(SNABLAR_0*SNABLAR_0+SNABLAZ_0*SNABLAZ_0)  &
          +Unablarho(1,6)*(SNABLAR_1*SNABLAR_1+SNABLAZ_1*SNABLAZ_1)  &
          +(Unablarho(2,6)+Unablarho(3,6))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)  &
          +UrhonablaJ(0,6)*RHO_0*DJ_0+UrhonablaJ(1,6)*RHO_1*DJ_1  &  
          +UrhonablaJ(2,6)*RHO_0*DJ_1+UrhonablaJ(3,6)*RHO_1*DJ_0  &  
          +UJnablarho(0,6)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))  &
          +UJnablarho(1,6)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))  &
          +UJnablarho(2,6)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))  &
          +UJnablarho(3,6)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1)) 
     tUt(1)=tUt(1)+UJJ(0,6)*J2_0+UJJ(1,6)*J2_1  &
          +(UJJ(2,6)+UJJ(3,6))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
     !     
     ! Contributions in the case 'u' depends on DeltaRHO_0     
     !
     tUDr(1)=tUDr(1)+Urhotau(0,7)*RHO_0*TAU_0  &  
          +Urhotau(1,7)*RHO_1*TAU_1+Urhotau(2,7)*RHO_0*TAU_1  &
          +Urhotau(3,7)*RHO_1*TAU_0+Urhorho(0,7)*RHO_0**2  &
          +Urhorho(1,7)*RHO_1**2+(Urhorho(2,7)+Urhorho(3,7))*RHO_0*RHO_1  &
          +UrhoDrho(0,7)*RHO_0*DRHO_0+UrhoDrho(1,7)*RHO_1*DRHO_1  &
          +UrhoDrho(2,7)*RHO_0*DRHO_1+UrhoDrho(3,7)*RHO_1*DRHO_0  &
          +Unablarho(0,7)*(SNABLAR_0*SNABLAR_0+SNABLAZ_0*SNABLAZ_0)  &
          +Unablarho(1,7)*(SNABLAR_1*SNABLAR_1+SNABLAZ_1*SNABLAZ_1)  &
          +(Unablarho(2,7)+Unablarho(3,7))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)  &
          +UrhonablaJ(0,7)*RHO_0*DJ_0+UrhonablaJ(1,7)*RHO_1*DJ_1  &  
          +UrhonablaJ(2,7)*RHO_0*DJ_1+UrhonablaJ(3,7)*RHO_1*DJ_0  &  
          +UJnablarho(0,7)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))  &
          +UJnablarho(1,7)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))  &
          +UJnablarho(2,7)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))  &
          +UJnablarho(3,7)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1)) 
     tUDr(1)=tUDr(1)+UJJ(0,7)*J2_0+UJJ(1,7)*J2_1  &
          +(UJJ(2,7)+UJJ(3,7))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
     !
     ! proton-neutron representation
     pUr(1)  =tUr(1)+tUr(2);     pUr(2)  =tUr(1)  -tUr(2)
     pUt(1)  =tUt(1)+tUt(2)+hb0; pUt(2)  =tUt(1)  -tUt(2)+hb0
     pUDr(1) =tUDr(1)+tUDr(2);   pUDr(2) =tUDr(1) -tUDr(2)
     pUNr(1) =tUNr(1)+tUNr(2);   pUNr(2) =tUNr(1) -tUNr(2)
     pUNz(1) =tUNz(1)+tUNz(2);   pUNz(2) =tUNz(1) -tUNz(2)
     pUDj(1) =tUDj(1)+tUDj(2);   pUDj(2) =tUDj(1) -tUDj(2)
     pUFIZ(1)=tUFIZ(1)+tUFIZ(2); pUFIZ(2)=tUFIZ(1)-tUFIZ(2)
     pUZFI(1)=tUZFI(1)+tUZFI(2); pUZFI(2)=tUZFI(1)-tUZFI(2)
     pUFIR(1)=tUFIR(1)+tUFIR(2); pUFIR(2)=tUFIR(1)-tUFIR(2)
     pURFI(1)=tURFI(1)+tURFI(2); pURFI(2)=tURFI(1)-tURFI(2)
     !
     Do it=itmin,itmax   !! loop over n  & p
        ita=3-it
        ! constraining potential
        pUr(it)= pUr(it) -cx*(two*fh(ihli)**2-fl(ihli)**2)                                  
        ! coulomb
        If(it.Eq.2) Then
           If(icou.Ge.1) pUr(it)=pUr(it)+cou(ihli)                                          
           If(icou.Eq.2) pUr(it)=pUr(it)+CExPar*coex*ro(ihli,it)**p13                                     
        Endif
        ! pairing contribution to rearrangement term 
        pUr(it)=pUr(it)-CpV0(it-1) *CpV1(it-1) /rho_c*aka(ihli,it)**2 &
                       -CpV0(ita-1)*CpV1(ita-1)/rho_c*aka(ihli,ita)**2
        ! pairing contribution to delta dv(ihli,it)
        rsa0=(ro(ihli,it)+ro(ihli,ita))/rho_c
        If(it.Eq.1) Then
           dvn(ihli)=(CpV0(it-1)*(ONE-rsa0*CpV1(it-1)))*aka(ihli,it)    
        Else
           dvp(ihli)=(CpV0(it-1)*(ONE-rsa0*CpV1(it-1)))*aka(ihli,it)    
        End If
     End Do !it
     !
     vn(ihli)=pUr(1)       ; vp(ihli)=pUr(2)        !* RHO_ij  
     vhbn(ihli)=pUt(1)     ; vhbp(ihli)=pUt(2)      !* TAU_ij  
     vrn(ihli)=pUNr(1)     ; vrp(ihli)=pUNr(2)      !* NABLAr RHO__ij  
     vzn(ihli)=pUNz(1)     ; vzp(ihli)=pUNz(2)      !* NABLAz RHO__ij  
     vdn(ihli)=pUDr(1)     ; vdp(ihli)=pUDr(2)      !* DELTA RHO_ij  
     vsn(ihli)=pUDj(1)     ; vsp(ihli)=pUDj(2)      !* NABLA . J__ij  
     vSFIZn(ihli)=pUFIZ(1) ; vSFIZp(ihli)=pUFIZ(2)  !* JFIZ_ij                               
     vSZFIn(ihli)=pUZFI(1) ; vSZFIp(ihli)=pUZFI(2)  !* JZFI_ij                               
     vSFIRn(ihli)=pUFIR(1) ; vSFIRp(ihli)=pUFIR(2)  !* JFIR_ij                               
     vSRFIn(ihli)=pURFI(1) ; vSRFIp(ihli)=pURFI(2)  !* JRFI_ij   
     !
  End Do !ihli
  !
  Call get_CPU_time('field',1)
  !
End Subroutine field
!===============================================================================================
!
!=======================================================================
Subroutine gamdel
  !---------------------------------------------------------------------
  ! ph- and pp- matrices in configurational space 
  !---------------------------------------------------------------------
  Use HFBTHO
  Use pairing_HFBTHO !mario     
  Implicit None
  Integer(ipr) :: i,ih,il,ib,ibx,nd,nd2,nza,nra,nla,nsa,nsb,nsab
  Integer(ipr) :: ihil,laplus,im,JA,N1,N2,ndnd,n12,n21
  Integer(ipr) :: i1,i2,i3
  Real(pr)     :: qla,yi,y,y2,qha,qhla,xmi,u2,un,up,xxx
  Real(pr)     :: sml2,cnzaa,cnraa,SSU,SSD
  Real(pr)     :: FITW1,FITW2,FITW3,FITW4
  Real(pr)     :: fi1r,fi1z,fi2d,QHL1A,QH1LA  
  Real(pr)     :: vh,vdh,vsh,hbh,vsum
  Real(pr)     :: SRFIh,SFIRh,SFIZh,SZFIh,SNABLARh,SNABLAZh
  Real(pr)     :: xlam,xlam2,xlamy,xlamy2,xlap,xlap2,xlapy,xlapy2,XLAMPY
  Real(pr)     :: FIUN1,FIDN1,FIURN1,FIDRN1,FIUZN1,FIDZN1,FIUD2N1,FIDD2N1
  Real(pr)     :: FIUN2,FIDN2,FIURN2,FIDRN2,FIUZN2,FIDZN2,FIUD2N2,FIDD2N2
  Real(pr)     :: FIUN12,FIDN12,FIURN12,FIDRN12,FIUZN12,FIDZN12
  Real(pr)     :: vnhl,vrnhl,vznhl,vdnhl,vsnhl,vhbnhl,vSRFInhl,vSFIRnhl
  Real(pr)     :: vSFIZnhl,vSZFInhl,vphl,vrphl,vzphl,vdphl,vsphl,vhbphl
  Real(pr)     :: vSRFIphl,vSFIRphl,vSFIZphl,vSZFIphl,dvnhl,dvphl
  Integer(ipr) :: ibro
  !
  Call get_CPU_time('gamdel',0)
  !
  !----------------------------------------------
  ! START BLOCKS
  !----------------------------------------------
  brout=zero; ibro=0 
  Do ib=1,NB
     ND=ID(ib); IM=ia(ib); ibx=ib+nbx
     If(Parity) Then
        LAPLUS=(ib+1)/2 !Yesp
     Else
        LAPLUS=ib       !Nop
     Endif
     XLAP=LAPLUS; XLAM=XLAP-ONE; xlap2=xlap*xlap; xlam2=xlam*xlam
     !----------------------------------------------
     ! SUM OVER GAUSS INTEGRATION POINTS
     !----------------------------------------------
     Do ihil=1,nghl
        y=y_opt(ihil); xlamy=xlam*y;     xlapy=xlap*y;   XLAMPY=XLAMY+XLAPY
        y2=y*y;        xlamy2=xlam2*y2;  xlapy2=xlap2*y2   
        !
        vnhl=vn(ihil);         vrnhl=vrn(ihil);       vznhl=vzn(ihil);       vdnhl=vdn(ihil)
        vsnhl=vsn(ihil);       vhbnhl=vhbn(ihil);     vSRFInhl=vSRFIn(IHIL); vSFIRnhl=vSFIRn(IHIL)
        vSFIZnhl=vSFIZn(IHIL); vSZFInhl=vSZFIn(IHIL); vphl=vp(ihil);         vrphl=vrp(ihil)
        vzphl=vzp(ihil);       vdphl=vdp(ihil);       vsphl=vsp(ihil);       vhbphl=vhbp(ihil)  
        vSRFIphl=vSRFIp(IHIL); vSFIRphl=vSFIRp(IHIL); vSFIZphl=vSFIZp(IHIL); vSZFIphl=vSZFIp(IHIL)
        dvnhl=dvn(ihil);       dvphl=dvp(ihil)
        !
        Do N1=1,ND
           JA=IM+N1;               NSA=NS(JA);             SSU=Max(NSA,0);         SSD=Max(-NSA,0)
           QHLA=QHLA_opt(JA,ihil); FI1R=FI1R_opt(JA,ihil); FI1Z=FI1Z_opt(JA,ihil); FI2D=FI2D_opt(JA,ihil)                             
           FIU(N1)=QHLA*SSU;       FIUR(N1)=fi1r*SSU;      FIUZ(N1)=fi1z*SSU;      FIUD2N(N1)=(FI2D-XLAMY2*QHLA)*SSU
           FID(N1)=QHLA*SSD;       FIDR(N1)=fi1r*SSD;      FIDZ(N1)=fi1z*SSD;      FIDD2N(N1)=(FI2D-XLAPY2*QHLA)*SSD
        End Do
        !
        I=ibro
        Do N1=1,ND
           JA=IM+N1;               NSA=NS(JA)
           FIUN1=FIU(N1);          FIURN1=FIUR(N1);        FIUZN1=FIUZ(N1);        FIUD2N1=FIUD2N(N1)
           FIDN1=FID(N1);          FIDRN1=FIDR(N1);        FIDZN1=FIDZ(N1);        FIDD2N1=FIDD2N(N1)
           Do N2=1,N1
              I=I+1; i1=i+nhhdim; i2=i+nhhdim2; i3=i+nhhdim3; NSB=NS(N2+IM); NSAB=NSA+NSB
              If (NSAB.Ne.0) Then
                 If (NSB.Gt.0) Then                                                                 !spin:UpUp
                    !SRFIh=ZERO; SFIRh=ZERO; SZFIh=ZERO
                    FIUN2    = FIU(N2); FIURN2=FIUR(N2); FIUD2N2=FIUD2N(N2); FIUZN2=FIUZ(N2)
                    vh       = FIUN1*FIUN2
                    hbh      = vh*XLAMY2+FIURN1*FIURN2+FIUZN1*FIUZN2
                    vdh      = hbh+hbh+FIUN1*FIUD2N2+FIUN2*FIUD2N1
                    SNABLARh = FIURN1*FIUN2+FIURN2*FIUN1
                    SNABLAZh = FIUZN1*FIUN2+FIUZN2*FIUN1
                    vsh      = SNABLARh*XLAMY                                 
                    SFIZh    = (vh+vh)*XLAMY
                 Else                                                                               !spin:DoDo
                    !SRFIh=ZERO; SFIRh=ZERO; SZFIh=ZERO
                    FIDN2    = FID(N2); FIDRN2=FIDR(N2); FIDZN2=FIDZ(N2); FIDD2N2=FIDD2N(N2) 
                    vh       = FIDN1*FIDN2
                    hbh      = vh*XLAPY2+FIDRN1*FIDRN2+FIDZN1*FIDZN2
                    vdh      = hbh+hbh+FIDN1*FIDD2N2+FIDN2*FIDD2N1;     
                    SNABLARh = FIDRN1*FIDN2+FIDRN2*FIDN1
                    SNABLAZh = FIDZN1*FIDN2+FIDZN2*FIDN1
                    vsh      =-SNABLARh*XLAPY
                    SFIZh    =-(vh+vh)*XLAPY
                 End If
                 brout(i )=brout(i )+vSFIZnhl*SFIZh+vh*vnhl+SNABLARh*vrnhl+SNABLAZh*vznhl+vdh*vdnhl+vsh*vsnhl+hbh*vhbnhl                     
                 brout(i1)=brout(i1)+vSFIZphl*SFIZh+vh*vphl+SNABLARh*vrphl+SNABLAZh*vzphl+vdh*vdphl+vsh*vsphl+hbh*vhbphl 
				 brout(i2)=brout(i2)+vh*dvnhl                  
				 brout(i3)=brout(i3)+vh*dvphl
              Else
                 If (NSB.Gt.0) Then                                                                !spin:DoUp
                    !vh=ZERO; hbh=ZERO; vdh=ZERO; SNABLARh=ZERO; SNABLAZh=ZERO; SFIZh=ZERO
                    FIUN2 = FIU(N2); FIURN2=FIUR(N2); FIUD2N2=FIUD2N(N2); FIUZN2=FIUZ(N2)
                    FITW3 =-FIDZN1*FIUN2; FITW4=FIUZN2*FIDN1 
                    vsh   =-FIDRN1*FIUZN2+FIURN2*FIDZN1+FITW3*XLAMY-FITW4*XLAPY
                    SRFIh =-FIDRN1*FIUN2+FIURN2*FIDN1
                    SFIRh = FIDN1*FIUN2*XLAMPY
                    SZFIh = FITW3+FITW4
                 Else                                                                             !spin:UpDo
                    !vh=ZERO; hbh=ZERO; vdh=ZERO; SNABLARh=ZERO; SNABLAZh=ZERO; SFIZh=ZERO
                    FIDN2 = FID(N2); FIDRN2=FIDR(N2); FIDZN2=FIDZ(N2); FIDD2N2=FIDD2N(N2) 
                    FITW3 =-FIDZN2*FIUN1; FITW4=FIUZN1*FIDN2
                    vsh   = FIURN1*FIDZN2-FIDRN2*FIUZN1-FITW4*XLAPY+FITW3*XLAMY                                 
                    SRFIh = FIURN1*FIDN2-FIDRN2*FIUN1
                    SFIRh = FIUN1*FIDN2*XLAMPY
                    SZFIh = FITW3+FITW4
                 Endif
                 brout(i )=brout(i )+vsh*vsnhl+vSRFInhl*SRFIh+vSFIRnhl*SFIRh+vSZFInhl*SZFIh
                 brout(i1)=brout(i1)+vsh*vsphl+vSRFIphl*SRFIh+vSFIRphl*SFIRh+vSZFIphl*SZFIh
              End If
              !----------------------------------------------
              ! LN PH PART
              !----------------------------------------------
              If(kindhfb.Lt.0) Then
                 If(ihil.Eq.1) Then
                    un=zero; up=zero; 
                    If(N1.Eq.N2) Then
                       un=-ala2(1); up=-ala2(2)
                    End If
                    n12=N1+(N2-1)*ND
                    brout(i )=brout(i )+two*(ala2(1)*rk(n12,ib )+un)
                    brout(i1)=brout(i1)+two*(ala2(2)*rk(n12,ibx)+up)
                 End If
              End If
           End Do !N2                      
        End Do !N1           
     End Do !ihil
     ibro=i
  End Do !IB
  !
  Call get_CPU_time('gamdel',1)
  ! 
  ! TMR pairing matrix elements 
  If(use_TMR_pairing.Ne.0) Then   
   !if(inin.lt.0.or.(inin.gt.0.and.iiter.gt.0)) 
   Call delta
  Endif
  !
  !----------------------------------------------
  ! BROYDEN/LINEAR MIXING
  !----------------------------------------------
  !
  Call get_CPU_time('broyden',0)
  !
  Call broyden_min(nhhdim4,brout,brin,alphamix,si,iiter,nbroyden,bbroyden)
  !
  Call get_CPU_time('broyden',1)
  !
End Subroutine gamdel
!=======================================================================
!
!======================================
! lib PnProjected specIfics Start >>>>>
!======================================
!=======================================================================
Subroutine densitpj
  !---------------------------------------------------------------------
  ! calculate local densities
  ! calculate local densities in mixed canonical (\for rho C, Y) and
  ! qp (for\tilde{\rho}) representation therefore E(can) \equiv E(qp)
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  !
  Complex(pr) :: tpfiu1,tpfid1,v2ig,dig,sumsum
  Complex(pr), Allocatable :: ank1(:,:),pfiun1(:,:),pfidn1(:,:)
  Complex(pr), Allocatable :: pakapj(:),propj(:), pdjpj(:), ptaupj(:),pdropj(:)
  Complex(pr), Allocatable :: pszfipj(:),psfizpj(:),psrfipj(:),psfirpj(:)
  Complex(pr), Pointer:: ppjk(:),pcpj(:,:),prpj(:,:),pypj(:,:)
  Real(pr)    :: f,s,sd,su,sud,y,y2,sml2,cnzaa,cnraa,u,v2,tauin,xxx,yyy
  Real(pr)    :: aav,anik,anik2,qhla,qh1la,qhl1a,qla,qha,fi1r,fi1z,fi2d
  Real(pr)    :: xlam,xlam2,xlamy,xlamy2,xlap,xlap2,xlapy,xlapy2,xlampy
  Real(pr)    :: tfiu,tfid,tfiur,tfidr,tfiuz,tfidz,tfiud2,tfidd2,tpfiu2,tpfid2
  Real(pr), Allocatable :: an2(:),ank2(:),pfiun2(:),pfidn2(:)
  Integer(ipr)    :: iw,nsa,nza,nra,nla,k,i,nd,il,ih,ihil,laplus,kkymu,n12
  Integer(ipr)    :: imen,ib,m,im,ig,it,j,jj,ja,jn,ILIHLI,k1,k2,kkk,kky,mu,kkkmu
  Integer(ipr)    :: k0(2),ky(2),kk(nqx),kyk(nqx)
  !
  Allocate(ank1(nqx,ilpj),pfiun1(ndx,ilpj),pfidn1(ndx,ilpj))
  Allocate(pfiun2(ndx),pfidn2(ndx),an2(nqx),ank2(nqx))
  Allocate(pakapj(ilnghl),propj(ilnghl), pdjpj(ilnghl), ptaupj(ilnghl),pdropj(ilnghl)  &
       ,pSZFIpj(ilnghl),pSFIZpj(ilnghl),pSRFIpj(ilnghl),pSFIRpj(ilnghl))
  !
  ! Projection grid points
  ! keypj=max(1,keypj); ilpj=keypj;  ilpj2=ilpj**2 !all 
  ! when a value two*pi is used the results are precisely the same 
  ! but the accuracy for even L is slow with increasing L. 
  ! when 'pi' is used it gives regular and better convergence 
  ! with respect to both, odd and even, L.
  ! write(*,'(2x,a,i2,a,f12.8,a,f12.8)') 'point ig=',i,' phi=',yyy,' pi/2=',pi/two
  xxx=pi/Dble(ilpj) ! equivalent to xxx=two*pi/Dble(ilpj)
  Do i=1,ilpj
     yyy   =Dble(i-1)*xxx
     phypj(i)=onei*yyy
     sinphy(i)=onei*Sin(yyy)
     exp1iphy(i)=Exp(onei*yyy)
     exp1iphym(i)=Exp(-onei*yyy)
     exp2iphy(i)=Exp(two*onei*yyy)
     exp2iphym(i)=Exp(-two*onei*yyy)
  End Do
  !
  ! initialize parameters
  varmas=zero
  !
  Do it=itmin,itmax
     !
     ! zero for densities
     Do J=1,ilnghl
        pakapj(J)=zero; propj(J)=zero; pdjpj(J)=zero; ptaupj(J)=zero; pdropj(J)=zero; 
     End Do
     Do J=1,ilnghl
        pszfipj(J)=zero; psfizpj(J)=zero; psrfipj(J)=zero; psfirpj(J)=zero; 
     End Do
     !
     ! it-pointers
     prpj=> rpj(:,:,it);  pcpj=> cpj(:,:,it); 
     pypj=> ypj(:,:,it);  ppjk=> pjk(:,it); 
     !
     ! null for all pointers
     pypj=zero; prpj=zero; pcpj=zero; ppjk=one; 
     !
     ! particle-init (kkk-even: 2 x number of pairs)
     kkk=npr(it); If(kkk.Ne.2*(kkk/2)) kkk=npr(it)-1 
     ppjk(1:ilpj)=exp1iphym(1:ilpj)**kkk 
     !
     ! start blocks
     k0(it)=0; ky(it)=0
     Do ib=1,nb
        nd=id(ib); im=ia(ib) 
        If(Parity) Then
           LAPLUS=(ib+1)/2 !Yesp
        Else
           LAPLUS=ib       !Nop
        Endif
        xlap=laplus; xlap2=xlap*xlap; xlam=xlap-one; xlam2=xlam*xlam
        !
        ! charge block can quantities
        m=ib+(it-1)*nbx; k1=ka(ib,it)+1; k2=ka(ib,it)+kd(ib,it); imen=0
        If(k1.Le.k2) Then
           ! below the pwi cut-off
           imen=nd
           !lcanon(ib,it)=lc
           Do k=1,nd
              k0(it)=k0(it)+1; kk(k)=k0(it); kkk=k0(it)
              ky(it)=ky(it)+1; kyk(k)=ky(it); kky=ky(it)
              aav=vk(kkk,it)                                         ! v^2 
              Do ig=1,ilpj
                 v2ig  =exp2iphy(ig)*aav                          ! gauged v^2
                 dig   =one-aav+v2ig                          ! denominator
                 If(kkk.Ne.blocanon(it)) Then
                    ppjk(ig)=ppjk(ig)*dig                           ! y(ig,it) <<<<<
                 Endif
                 prpj(kkk,ig)=v2ig/dig                                  ! rho(mu,ig,it)
                 pcpj(kky,ig)=exp2iphy(ig)/dig                          ! c(mu,ig,it)
                 pypj(kky,ig)=exp1iphy(ig)/dig*onei*Sin(phypj(ig)/onei) ! sinphy(ig) !Y(mu,ig,it)
              End Do
           End Do
           ! At this point density (and related) are strictly equivalent in qp-and can-representation
           ! (up to 10^-14). Pairing density is not so strict (up to 10^-5) due to uv from v^2 but
           ! pairing density is taken directly in qp representation so both representations
           ! qp and can are strictly exact (up to 10^-14).  
           j=0
           Do jj=1,nd
              Do k=1,nd
                 j=j+1; n12=jj+(k-1)*nd; 
                 an2(j)=ddc(jj,kk(k),it)                   
                 ank2(j)=ak(n12,m)                 ! half \tilde{\rho} in q.p. basis
                 Do ig=1,ilpj
                    ank1(j,ig)=zero
                 End Do
                 Do mu=1,nd                          ! for half e^(-i\phy)*C(\phy)*\tilde{\rho} in q.p. basis
                    kkkmu=kk(mu); kkymu=kyk(mu)
                    Do ig=1,ilpj                     ! e^(-i\phy)*C in q.p. basis
                       ank1(j,ig)=ank1(j,ig)+ddc(jj,kkkmu,it)*ddc(k,kkkmu,it)*pcpj(kkymu,ig)*exp1iphym(ig) 
                    End Do
                 End Do
              End Do
           End Do
        Else
           ! above the pwi cut-off (NB! Attention)
           ! here imem=0 and the contribution does 
           ! not enter the densities but the Hamiltonian matrix
           ! used only in VAP regime
           ky(it)=ky(it)+1; kky=ky(it)
           Do ig=1,ilpj
              pcpj(kky,ig)=exp2iphy(ig)
              pypj(kky,ig)=exp1iphy(ig)*sinphy(ig)
           End Do
        End If
        !
        ! calculate the densities only below the PWI cutoff
        If (imen.Gt.0) Then
           ! gauss integration points
           Do il=1,ngl
              v2=half/xl(il)
              Do ih=1,ngh
                 ihil=ih+(il-1)*ngh; ilihli=(ihil-1)*ilpj
                 u=xh(ih); y=fli(ihil); y2=y*y
                 xlamy=xlam*y; xlamy2=xlam2*y2; 
                 xlapy=xlap*y; xlapy2=xlap2*y2; 
                 xlampy=xlamy+xlapy
                 !
                 ! initialize spin up/down funct
                 Do k=1,nd
                    fiu(k)=zero; fiuz(k)=zero; fiur(k)=zero; fiud2n(k)=zero; pfiun2(k)=zero; 
                    fid(k)=zero; fidz(k)=zero; fidr(k)=zero; fidd2n(k)=zero; pfidn2(k)=zero; 
                    Do ig=1,ilpj
                       pfiun1(k,ig)=zero; pfidn1(k,ig)=zero
                    End Do
                 End Do
                 !
                 ! scan over basis states
                 jn=0
                 Do i=1,nd
                    ja=i+im; nla=nl(ja); nra=nr(ja); nza=nz(ja); nsa=ns(ja); 
                    sml2=nla*nla; cnzaa=nza+nza+1; cnraa=nra+nra+nla+1
                    qha=qh(nza,ih); qla=ql(nra,nla,il); qhla=qha*qla
                    qhl1a=qha*ql1(nra,nla,il)*v2; qh1la=qh1(nza,ih)*qla
                    fi1z=fp1(ihil)*qhla+fp2(ihil)*qh1la+fp3(ihil)*qhl1a
                    fi1r=fp4(ihil)*qhla+fp5(ihil)*qh1la+fp6(ihil)*qhl1a
                    fi2d=(fs1(ihil)*qh1la**2+four*fs4(ihil)*qh1la*qhl1a  &
                         +fs2(ihil)*qhl1a**2+two*(fs5(ihil)*qh1la  &
                         +fs6(ihil)*qhl1a)*qhla+((u*u-cnzaa)*fs1(ihil)  &
                         +(p14-cnraa*v2+sml2*v2*v2)*fs2(ihil)+fs3(ihil))*qhla**2  &
                         -two*(fi1r**2+fi1z**2))/(two*qhla)
                    !
                    ! wave function(spin:up,down; grad:r,z,d2)
                    If (nsa.Gt.0)   Then
                       Do  k=1,nd
                          jn=jn+1; anik=an2(jn); anik2=ank2(jn)
                          Do ig=1,ilpj
                             pfiun1(k,ig)=pfiun1(k,ig)+ank1(jn,ig)*qhla
                          End Do
                          pfiun2(k)=pfiun2(k)+anik2*qhla
                          fiu(k)=fiu(k)+anik*qhla
                          fiur(k)=fiur(k)+anik*fi1r
                          fiuz(k)=fiuz(k)+anik*fi1z
                          fiud2n(k)=fiud2n(k)+anik*(fi2d-xlamy2*qhla)
                          !
                       End Do
                    Else
                       Do k=1,nd
                          jn=jn+1; anik=an2(jn); anik2=ank2(jn)
                          Do ig=1,ilpj
                             pfidn1(k,ig)=pfidn1(k,ig)+ank1(jn,ig)*qhla
                          End Do
                          pfidn2(k)=pfidn2(k)+anik2*qhla
                          fid(k)=fid(k)+anik*qhla
                          fidr(k)=fidr(k)+anik*fi1r
                          fidz(k)=fidz(k)+anik*fi1z
                          fidd2n(k)=fidd2n(k)+anik*(fi2d-xlapy2*qhla)
                          !
                       End Do
                    End If
                 End Do ! i
                 !
                 ! calculate densities
                 Do k=1,nd
                    kkk=kk(k)
                    tfiu=fiu(k); tfiuz=fiuz(k); tfiur=fiur(k); tfiud2=fiud2n(k); tpfiu2=pfiun2(k);  
                    tfid=fid(k); tfidz=fidz(k); tfidr=fidr(k); tfidd2=fidd2n(k); tpfid2=pfidn2(k); 
                    Do ig=1,ilpj
                       I=ig+ilihli; v2ig=prpj(kkk,ig); tpfiu1=pfiun1(k,ig); tpfid1=pfidn1(k,ig)
                       !
                       pakapj(I)=pakapj(I)+(tpfiu1*tpfiu2+tpfid1*tpfid2)
                       propj(I)=propj(I)+(tfiu**2+tfid**2)*v2ig
                       pdjpj(I)=pdjpj(I)+(tfiur*tfidz-tfidr*tfiuz+xlamy*tfiu*(tfiur-tfidz)  &
                            -xlapy*tfid*(tfidr+tfiuz))*v2ig
                       tauin=(tfiur**2+tfidr**2+tfiuz**2+tfidz**2+xlamy2*tfiu**2+xlapy2*tfid**2)
                       ptaupj(I)=ptaupj(I)+tauin*v2ig
                       pdropj(I)=pdropj(I)+(tauin+tfiu*tfiud2+tfid*tfidd2)*v2ig
                       psrfipj(I)=psrfipj(I)+(tfiur*tfid-tfidr*tfiu)*v2ig
                       psfirpj(I)=psfirpj(I)+(tfiu*tfid*xlampy)*v2ig
                       psfizpj(I)=psfizpj(I)+(xlamy*tfiu**2-xlapy*tfid**2)*v2ig
                       pszfipj(I)=pszfipj(I)+(tfiuz*tfid-tfidz*tfiu)*v2ig
                       !
                    End Do !ig
                 End Do !k
              End Do !ih
           End Do !il
        End If
     End Do !ib
     !
     ! normalized pjk
     sumsum=Sum(ppjk(1:ilpj)); ppjk(1:ilpj)=ppjk(1:ilpj)/sumsum
     !
     ! Y minus second term of Y
     Do k=1,ky(it)
        sumsum=Sum(ppjk(1:ilpj)*pypj(k,1:ilpj))
        pypj(k,1:ilpj)=pypj(k,1:ilpj)-sumsum
     End Do
     !
     ! norm of the projected/unprojected density
     s=zero; sd=zero; su=zero; sud=zero; 
     Do ihil=1,nghl
        ilihli=(ihil-1)*ilpj
        Do ig=1,ilpj
           I=ig+ilihli
           s=s+two*propj(I)*ppjk(ig); sd=sd+four*pdropj(I)*ppjk(ig)
        End Do
        I=1+ilihli
        su=su+two*propj(I); sud=sud+four*pdropj(I)
     End Do
     !
     ! print unprojected normalization
     Do iw=lout,lfile
        Write(iw,'(2(a,2(2x,D14.8)),(a,D14.8),a,i3)')  &
             '   pj/unpj  s=',s,su,'   pj/unpj sd=',sd,sud,' ala1=',ala1(it),' inner=',inner(it)
     Enddo
     varmas=varmas+su
     varmasnz(it)=su; pjmassnz(it)=s
     !
     s=Dble(npr(it))/s; dnfactor(it)=s; drhoi(it)=sd
     !
     Do ihil=1,nghl
        ilihli=(ihil-1)*ilpj
        ! wdcor moves out the int.weight and multiply by the jacobian
        f=two*wdcori(ihil)
        Do ig=1,ilpj
           I=ig+ilihli
           ropj (ihil,ig,it)=f*propj (I)
           taupj(ihil,ig,it)=f*ptaupj(I)
           dropj(ihil,ig,it)=f*pdropj(I)*two
           djpj (ihil,ig,it)=f*pdjpj (I)*two
           akapj(ihil,ig,it)=f*pakapj(I)*half
           srfipj(ihil,ig,it)=f*psrfipj(I)
           sfirpj(ihil,ig,it)=f*psfirpj(I)
           sfizpj(ihil,ig,it)=f*psfizpj(I)
           szfipj(ihil,ig,it)=f*pszfipj(I)
        End Do !ig
     End Do !ihil
     !
  End Do !it
  !
  dnfactor(3)=dnfactor(1)+dnfactor(2)
  !
  Deallocate(ank1,pfiun1,pfidn1)
  Deallocate(pakapj,propj, pdjpj, ptaupj,pdropj,pSZFIpj,pSFIZpj,pSRFIpj,pSFIRpj)
  !
  Call coulompj !complex coulomb fields
  !
End Subroutine densitpj
!=======================================================================
!
!=======================================================================
Subroutine coulompj
  !---------------------------------------------------------------------
  ! Coulom-field (direct part)
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Integer(ipr) :: i,j,k
  Real(pr) :: zd2,y1,y2,xx1,s1,vik,f,r,r1,rr2,z,z1,zd1,t
  Real(pr) :: bb,r2,r12,rrr,rz1,rz2,rrz1,rrz2,xx
  Real(pr) :: bb1=3.5156229_pr,g1=0.39894228_pr,g7=0.02635537_pr,  &
       bb2=3.0899424_pr,g2=0.01328592_pr,g8=0.01647633_pr,  &
       bb3=1.2067492_pr,g3=0.00225319_pr,g9=0.00392377_pr,  &
       bb4=0.2659732_pr,g4=0.00157565_pr,bbxx=3.750_pr,  &
       bb5=0.0360768_pr,g5=0.00916281_pr,  &
       bb6=0.0045813_pr,g6=0.02057706_pr
  If(icacoupj.Eq.0) Then
     icacoupj=1; bb=Max(bp,bz)**4; f=chargee2/Sqrt(pi);  ! f=e^2/Sqrt(pi)
     Do i=1,nghl
        r=fl(i); z=fh(i); r2=r*r
        Do k=1,i
           r1=fl(k); z1=fh(k); r12=r1*r1
           rrr=two*r*r1;    rr2=(r-r1)**2
           zd1=(z-z1)**2; zd2=(z+z1)**2
           rz1=r2+r12+zd1;  rz2=r2+r12+zd2
           rrz1=rr2+zd1;     rrz2=rr2+zd2
           !
           xx1=zero
           Do j=1,nleg
              xx=Sqrt(one-xleg(j)**2); y1=(xleg(j)/(bb*xx))**2; s1=y1*rrr
              If(s1.Le.bbxx) Then
                 t=(s1/bbxx)**2; y2=one+t*(bb1+t*(bb2+t*(bb3+t*(bb4+t*(bb5+t*bb6)))))
                 y2=y2*(Exp(-rz1*y1)+Exp(-rz2*y1))
              Else
                 t=(bbxx/s1); y2=g1+t*(g2+t*(g3+t*(-g4+t*(g5+t*(-g6+t*(g7+t*(-g8+t*g9)))))))
                 y2=y2/Sqrt(s1)*(Exp(-rrz1*y1)+Exp(-rrz2*y1))
              End If
              xx1=xx1+wleg(j)*y2/(bb*xx**3)
           End Do
           vik=f*xx1; vc(i,k)=vik*wdcor(k); vc(k,i)=vik*wdcor(i)  !wdcor=pi*wh*wl*bz*bp*bp/fd
        End Do  !k
     End Do  !i
  End If
  ! calculation of the coulomb field
  coupj=zero
  Do i=1,nghl
     Do k=1,ilpj
        coupj(:,k)=coupj(:,k)+vc(:,i)*ropj(i,k,2)
     End Do
  End Do
End Subroutine coulompj
!=======================================================================
!
!=======================================================================
Subroutine broyden_min(N,vout,vin,alpha,si,iter,M,bbroyden)
  !---------------------------------------------------------------------
  ! Modified Broyden's method: D.D.Johnson, PRB 38, 12807 (1988)
  ! Adopted from: (C) 2001 PWSCF group
  ! Input :
  !  N      dimension of arrays vin,vout
  !  vin    outpu at previous iteration
  !  vout   output at current iteration
  !  alpha  mixing factor (0 < alpha <= 1)
  !  iter   current iteration number
  !  M      number of iterations in Broyden history
  !  M=0    Linear mixing
  ! Output:
  !  si     MaxVal(|vout-vin|)
  !  vin    Broyden/Linear mixing result
  !  vout   vout-vin
  !  bbroyden='B' Broyden mixing, curvature>0
  !  bbroyden='L' Linear mixing,  curvature<0
  !---------------------------------------------------------------------
  Use HFBTHO, Only: pr,ipr,ierror_flag,ierror_info
  Implicit None
  Integer(ipr),     Intent(In)    :: N,iter,M
  Real(pr),        Intent(In)     :: alpha
  Real(pr),        Intent(Out)    :: si
  Character(1),      Intent(Out)  :: bbroyden
  Real(pr),        Intent(InOut)  :: vout(N),vin(N)  
  Integer(ipr)                    :: i,j,iter_used,ipos,inext
  Integer(ipr), Allocatable, Save :: iwork(:)
  Real(pr),    Allocatable, Save  :: beta(:,:),work(:)
  Real(pr),    Allocatable, Save  :: df(:,:),dv(:,:),curv(:)
  Real(pr),                 Save  :: w0
  Real(pr)                        :: DDOT,DNRM2,normi,gamma,curvature,sf
  !
  sf=-1.0_pr; Call DAXPY(N,sf,vin,1,vout,1)  
  si=Maxval(Abs(vout))
  ! Linear mixing
  If(M.Eq.0.Or.iter.Eq.0) Then
     bbroyden='L'; Call DAXPY(N,alpha,vout,1,vin,1)  
     !If(iter.Eq.0) Write(6,*) '  Linear mixing (alpha) : ',alpha
     Return
  End If
  ! Broyden mixing
  iter_used=Min(iter-1,M) 
  ipos=iter-1-((iter-2)/M)*M
  inext=iter-((iter-1)/M)*M
  If (iter.Eq.1) Then       
     w0=0.010_pr
     If(Allocated(df)) Deallocate(curv,df,dv,beta,work,iwork)
     Allocate(curv(N),df(N,M),dv(N,M),beta(M,M),work(M),iwork(M))
     !write(6,'(a,i3,3(2x,f18.8),a)') '   Broyden mixing (M,alpha,w0,mem) : '  &
     !     ,M,alpha,w0,(2*N*M+N)*8._pr/1.e6,' MB'
  Else
     df(:,ipos)=vout(:)-df(:,ipos); dv(:,ipos)=vin(:)-dv(:,ipos)
     Normi=1.0_pr/Sqrt((DNRM2(N,df(1,ipos),1))**2)
     Call dscal(N,Normi,df(1,ipos),1)
     Call dscal(N,Normi,dv(1,ipos),1)
  Endif
  Do i=1,iter_used
     Do j=i+1,iter_used
        beta(i,j)=DDOT(N,df(1, j),1,df(1,i),1)
     Enddo
     beta(i,i)=1.0_pr+w0*w0
  Enddo
  Call DSYTRF('U',iter_used,beta,M,iwork,work,M,i)
  If(i.Ne.0) Then
     ierror_flag=ierror_flag+1
     ierror_info(ierror_flag)='STOP: In Broyden: info at DSYTRF '
     Return  
  Endif
  Call DSYTRI('U',iter_used,beta,M,iwork,work,i)
  If(i.Ne.0) Then
     ierror_flag=ierror_flag+1
     ierror_info(ierror_flag)='STOP: In Broyden: info at DSYTRI '
     Return  
  Endif
  Do i=1,iter_used
     Do j=i+1,iter_used
        beta(j,i)=beta(i,j)
     Enddo
     work(i)=DDOT(N,df(1,i),1,vout,1)      
  Enddo
  curv=alpha*vout
  Do i=1,iter_used
     gamma=0.0_pr
     Do j=1,iter_used
        gamma=gamma+beta(j,i)*work(j)
     Enddo
     curv=curv-gamma*(dv(:,i)+alpha*df(:,i))
  Enddo
  Call dcopy(N,vout,1,df(1,inext),1)
  Call dcopy(N,vin,1,dv(1,inext),1)
  curvature=DDOT(N,vout,1,curv,1)
  If(curvature.Gt.-1.0_pr) Then
     bbroyden='B'; sf=+1.0_pr; Call DAXPY(N,sf,curv,1,vin,1)  
  Else
     bbroyden='L'; sf=alpha*0.50_pr; Call DAXPY(N,sf,vout,1,vin,1)  
  End If
End Subroutine broyden_min
!=======================================================================
!
!=======================================================================
!======================START library gauss points=====================
!=======================================================================
Subroutine gausspoints
  !---------------------------------------------------------------------
  !PROGRAM DETERMINES POINTS AND WEIGHTS FOR THE GAUSS INTEGRATION.
  !GAUSS-LEGENDRE,LAGUERRE AND HERMITTE INTEGRATION.
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Real(pr):: Work(1000),al,be,sparity
  Integer(ipr) :: N,N1,N2,N3,i,j,j1,J2,KINDI,kpts,nparity
  ! 
  al=0.0_pr; be=0.0_pr; kpts=0
  !
  !--------------------------------------------------------------------
  !------------------>> Gauss-Hermite (positive nodes) <<--------------
  !--------------------------------------------------------------------
  If(Parity) Then
     KINDI=4; N=2*ngh; N1=3; N2=N1+N; N3=N2+N  !Yesp
  Else  
     KINDI=4;  N=ngh; N1=3; N2=N1+N; N3=N2+N   !Nop
  Endif
  Call Gaussq(KINDI,N,al,be,kpts,Work(1),Work(N1),Work(N2),Work(N3))
  If(ierror_flag.Ne.0) Return	   	                 
  !WRITE(11,*) ' !'
  !WRITE(11,*) ' !GAUSS-HERMITTE INTEGRATION POINTS AND WEIGHTS'
  If(Parity) Then
     nparity=ngh  !Yesp
     sparity=2.0_pr
  Else
     nparity=0    !Nop
     sparity=1.0_pr
  Endif
  Do i=nparity+1,N 
     J=i-nparity; J1=N2+i-1; J2=N3+i-1 
     xh(J)=Work(j1); wh(J)=sparity*Dexp(Work(j1)*Work(j1)+Dlog(Work(j2)))
     !write(11,'(a,i2,a,D24.18,a,I2,a,D24.18)')  '  xh(',J,')=',xh(J),';   wh(',J,')=',wh(J)
  End Do
  !--------------------------------------------------------------------
  !---------------------------->> Gauss-Laguerre <<--------------------|
  !--------------------------------------------------------------------
  KINDI=6; N=ngl; N1=3; N2=N1+N; N3=N2+N
  Call Gaussq(KINDI,N,al,be,kpts,Work(1),Work(N1),Work(N2),Work(N3))
  If(ierror_flag.Ne.0) Return	   	                 
  !WRITE(11,*) ' !'
  !WRITE(11,*) ' !GAUSS-LAGUERRE INTEGRATION POINTS AND WEIGHTS'
  Do J=1,ngl
     J1=N2+J -1; J2=N3+J -1
     xl(J)=Work(j1); wl(J)=Dexp(Work(j1)+Dlog(Work(j2))); sxl(J)=Sqrt(xl(J))
     !write(11,'(a,i2,a,D24.18,a,I2,a,D24.18)')  '  xl(',J,')=',xl(J),';   wl(',J,')=',wl(J)
  End Do
  !--------------------------------------------------------------------
  !----------------->> Gauss-Legendre (positive nodes) <<--------------
  !--------------------------------------------------------------------
  KINDI=1; N=2*nleg; N1=3; N2=N1+N; N3= N2+N   
  Call Gaussq(KINDI,N,al,be,kpts,Work(1),Work(N1),Work(N2),Work(N3))
  If(ierror_flag.Ne.0) Return	   	                 
  !WRITE(11,*) ' !'
  !WRITE(11,*) ' !GAUSS-LEGENDRE INTEGRATION POINTS AND WEIGHTS'
  Do J=1,nleg
     J1=N2+nleg+ J -1
     J2=N3+nleg+ J -1
     xleg(J)=Work(j1); wleg(J)=Work(j2); 
     !write(11,'(a,I2,a,D24.18,a,I2,a,D24.18)') ' xleg(',J,')=',xleg(J),'; wleg(',J,')=',wleg(J)
  End Do
  !  
End Subroutine gausspoints
!=======================================================================
!  
!=======================================================================
Subroutine GAUSSQ(kindi,N,ALPHA,BETA,KPTS,ENDPTS,B,T,W)
  Use HFBTHO, Only: pr,ipr
  Implicit None
  Integer(ipr) :: N,kindi
  Real(pr):: ALPHA,BETA,MUZERO,GAM,T1
  Integer(ipr) :: j1,J2,KIND,kpts,IERR
  Real(pr):: B(N),T(N),W(N),ENDPTS(2)
  Real(pr):: GBSLVE 
  !
  Call Class(kindi,N,ALPHA,BETA,B,T,MUZERO)
  !
  If(KPTS.Eq.0) Then
     W=0.0_pr; W(1)=1._pr
     Call GBTQL2(N,T,B,W,IERR)
     W=MUZERO*W*W
     Return
  Endif
  If(KPTS.Eq.2) Then
     GAM=GBSLVE(ENDPTS(1),N,T,B)
     T1=((ENDPTS(1)-ENDPTS(2))/(GBSLVE(ENDPTS(2),N,T,B)-GAM))
     B(N-1)=Sqrt(T1)
     T(N)=ENDPTS(1)+GAM*T1
     W=0.0_pr; W(1)=1._pr
     Call GBTQL2(N,T,B,W,IERR)
     W=MUZERO*W*W
     Return
  Endif
  T(N)=GBSLVE(ENDPTS(1),N,T,B)*B(N-1)**2+ENDPTS(1)
  W=0.0_pr; W(1)=1._pr
  Call GBTQL2(N,T,B,W,IERR)
  W=MUZERO*W*W
End Subroutine GAUSSQ
!=======================================================================
!
!=======================================================================
Double Precision Function GBSLVE(SHIFT,N,A,B)
  Use HFBTHO
  Implicit None
  Integer(ipr) :: N,NM1,i
  Real(pr) :: ALPHA,SHIFT,A(N),B(N)
  ALPHA=A(1)-SHIFT
  NM1=N-1
  Do I=2,NM1
     ALPHA=A(I)-SHIFT-B(I-1)**2/ALPHA
  Enddo
  GBSLVE=1.0_pr/ALPHA
End Function GBSLVE
!=======================================================================
!
!=======================================================================
Subroutine Class(kindi,N,ALPHA,BETA,B,A,MUZERO)
  Use HFBTHO, Only: pr,ipr
  Implicit None
  Integer(ipr) :: N,kindi,i,NM1
  Real(pr) :: MUZERO,ALPHA,BETA,A(N),B(N)
  Real(pr) :: PI,ABI,DI20,AB,A2B2,FI
  Real(pr) :: DGAMMA
  Data PI / 3.1415926535897930_pr /
  NM1=N-1
  Select Case (kindi)
  Case (1)
     MUZERO=2.0_pr
     Do  I=1,NM1
        A(I)=0_pr
        ABI=I
        B(I)=ABI/Sqrt(4.0_pr*ABI*ABI-1.0_pr)
     Enddo
     A(N)=0.0_pr
  Case (2)
     MUZERO=PI
     Do  I=1,NM1
        A(I)=0.0_pr
        B(I)=0.50_pr
     Enddo
     B(1)=Sqrt(0.50_pr)
     A(N)=0.0_pr
  Case (3)
     MUZERO=PI/2.0_pr
     Do  I=1,NM1
        A(I)=0.0_pr
        B(I)=0.50_pr
     Enddo
     A(N)=0.0_pr
  Case (4)
     MUZERO=Sqrt(PI)
     Do  I=1,NM1
        A(I)=0.0_pr
        DI20=I/2.0_pr
        B(I)=Sqrt(DI20)
     Enddo
     A(N)=0.0_pr
  Case (5)
     AB=ALPHA+BETA
     ABI=2.0_pr+AB
     MUZERO=2.0_pr**(AB+1.0_pr)*DGAMMA(ALPHA+1.0_pr)*DGAMMA(BETA+1.0_pr)/DGAMMA(ABI)
     A(1)=(BETA-ALPHA)/ABI
     B(1)=Sqrt(4.0_pr*(1.0_pr+ALPHA)*(1.0_pr+BETA)/((ABI+1.0_pr)*ABI*ABI))
     A2B2=BETA*BETA-ALPHA*ALPHA
     Do  I=2,NM1
        ABI=2.0_pr*I+AB
        A(I)=A2B2/((ABI-2.0_pr)*ABI)
        FI=I
        B(I)=Sqrt(4.0_pr*FI*(FI+ALPHA)*(FI+BETA)*(FI+AB)/((ABI*ABI-1.0_pr)*ABI*ABI))
     Enddo
     ABI=2.0_pr*N+AB
     A(N)=A2B2/((ABI-2.0_pr)*ABI)
  Case (6)
     MUZERO=DGAMMA(ALPHA+1.0_pr)
     Do  I=1,NM1
        FI=I
        A(I)=2.0_pr*FI-1.0_pr+ALPHA
        B(I)=Sqrt(FI*(FI+ALPHA))
     Enddo
     A(N)=2.0_pr*N-1.0_pr+ALPHA
  Case default
  End Select
End Subroutine Class
!=======================================================================
!
!=======================================================================
Subroutine GBTQL2 ( N,D,E,Z, IERR )
  Use HFBTHO
  Implicit None
  Integer(ipr) :: N,IERR
  Real(pr) :: D(N),E(N),Z(N)
  Integer(ipr) :: I,J,K,L,M,II,MML
  Real(pr) :: MACHEP,P,G,R,S,C,F,B
  MACHEP=16.0_pr**(-14)
  IERR=0
  If (N .Eq. 1 ) Return
  E(N)=0.0_pr
  Do  L= 1,N
     J=0
     Do
        Do M=L,N
           If (M .Eq. N) Exit
           If ( Abs(E(M)) .Le. MACHEP * (Abs(D(M))+Abs(D(M+1)))) Exit
           Continue
        Enddo
        P=D(L)
        If ( M .Eq. L ) Exit
        If (J .Eq. 30) Then
           IERR=L
           Return
        Endif
        J=J+1
        G=(D(L+1)-P) /(2.0_pr * E(L))
        R =Sqrt( G*G+1.0_pr )
        G=D(M)-P+ E(L) / (G+Sign(R , G))
        S=1.0_pr
        C=1.0_pr
        P=0.0_pr
        MML=M-L
        Do  II=1, MML
           I=M-II
           F=S * E(I)
           B=C * E(I)
           If (Abs(F).Ge.Abs(G)) Then
              C=G / F
              R =Sqrt(C*C+1.0_pr )
              E(I+1)=F*R
              S=1.0_pr /R
              C=C * S
           Else
              S=F / G
              R =Sqrt(S*S+1.0_pr )
              E(I+1)=G * R
              C=1.0  / R
              S=S * C
           Endif
           G=D(I+1)-P
           R=( D(I)-G) * S+2.0_pr * C * B
           P=S * R
           D(I+1)=G+P
           G=C * R-B
           F=Z(I+1)
           Z(I+1)=S * Z(I)+C * F
           Z(I)=C * Z(I)-S * F
        Enddo
        D(L)=D(L)-P
        E(L)=G
        E(M)=0.0_pr
     End Do
  Enddo
  Do  II=2 , N
     I=II-1
     K=I
     P=D(I)
     Do  J=II,N
        If ( D(J) .Ge. P) Cycle
        K=J
        P=D(J)
     Enddo
     If (K .Eq. I ) Cycle
     D(K)=D(I)
     D(I)=P
     P=Z(I)
     Z(I)=Z(K)
     Z(K)=P
  Enddo
End Subroutine GBTQL2
!=======================================================================
!
!=======================================================================
Double Precision Function DGAMMA(Z)
  Use HFBTHO
  Implicit None
  Real(pr) :: Z
  Integer(ipr) :: K,K1
  Real(pr) :: A(18),T,P
  A(1)=1.00_pr;                A(2)=.42278433509846780_pr;
  A(3)=.41184033042636720_pr;  A(4)=.08157691925026090_pr
  A(5)=.07424901068009040_pr;  A(6)=-.00026698103334840_pr
  A(7)=.01115403602403440_pr;  A(8)=-.00285258214461970_pr
  A(9)=.00210362870245980_pr;  A(10)=-.00091848436909910_pr
  A(11)=.00048742279447680_pr; A(12)=-.00023472040189190_pr
  A(13)=.00011153395196660_pr; A(14)=-.00004787479838340_pr
  A(15)=.00001751027271790_pr; A(16)=-.00000492037509040_pr
  A(17)=.00000091991564070_pr; A(18)=-.00000008399404960_pr
  If(Z .Le. 1.00_pr ) Then
     T=Z
     P=A(18)
     Do  K1=1,17
        K=18 -K1
        P=T*P+ A(K)
     Enddo
     DGAMMA=P /(Z * (Z+1.00_pr))
     Return
  Else If (Z .Gt. 1.0 ) Then
     DGAMMA=P / Z
     Return
  End If
  If ( Z .Le. 2.00_pr ) Then
     T=Z-1.00_pr
     P=A(18)
     Do  K1=1,17
        K=18 -K1
        P=T*P+ A(K)
     Enddo
     DGAMMA=P /(Z * (Z+1.00_pr))
     Return
  Else If(Z .Gt. 2.0 ) Then 
     DGAMMA=P
     Return
  End If
  T=Z-2.00_pr
  P=A(18)
  Do  K1=1,17
     K=18 -K1
     P=T*P+ A(K)
  Enddo
  DGAMMA=P /(Z * (Z+1.00_pr))
End Function DGAMMA
!=======================================================================
!====================END library gauss points=========================
!=======================================================================
!
!=======================================================================
Subroutine expect(lpr)
  !---------------------------------------------------------------------
  ! calculates expectation values (xn is the particle number)
  ! at lpr=.true. also calculates PAV corrections
  !---------------------------------------------------------------------
  Use HFBTHO
  Use pairing_HFBTHO 
  Implicit None
  Logical :: lpr
  Integer(ipr) :: i,it,ihli,iw
  Real(pr) :: ekt(3),xn(3),q4(3),def(3),bet2(3),het4(3),econst
  Real(pr) :: z,zz,rrr,p2,p3,p4,row,r212,r222,rc
  Real(pr) :: eso,ecodi,ecoex,rn,rp,rnp1,rnp2,rt,whl,tnt,tpt,tt
  Real(pr) :: dn,dp,dt,akn,akp,akn2,akp2,adn,adp,evol,esurf,ecoul
  Real(pr) :: etens,dd1n,dd1p,rt1,tt1,dt1,djn,djp,djt,djt1
  Real(pr) :: RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,DJ_0,DJ_1,J2_0,J2_1
  Real(pr) :: SZFIN,SFIZN,SRFIN,SFIRN,SZFIP,SFIZP,SRFIP,SFIRP
  Real(pr) :: SZFI_0,SFIZ_0,SRFI_0,SFIR_0,SZFI_1,SFIZ_1,SRFI_1,SFIR_1
  Real(pr) :: SNABLARN,SNABLAZN,SNABLARP,SNABLAZP
  Real(pr) :: SNABLAR_0,SNABLAZ_0,SNABLAR_1,SNABLAZ_1
  Real(pr) :: xn1,xn2,rms1,rms2,q21,q22,q41,q42,EKIN_N,EKIN_P,ept1,ept2,del1,del2
  Real(pr) :: rsa,rsa0,rsa0A,rps,rns,rsa1,rsa10,rsa12,rsa0An,rsa0As
  Real(pr) :: ESURF_rho_DELTA_rho,ESURF_NABLA_rho_NABLA_rho,ESO_rho_NABLA_J,ESO_NABLA_rho_J
  Real(pr) :: EVOL_rho_tau,EVOL_rho_rho,EExtra,E_HARTREE_DIR,tempE_Crho0,tempREARR
  Real(pr) :: E_EXT_FIELD
  Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
  Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
  !------------------------------------------------
  ! Part called during iterations (lpr=F)
  !------------------------------------------------
  !   `
  Call densit
  If(ierror_flag.Ne.0) Return	   	                 
  !
  Call get_CPU_time('expect',0)
  !
  !------------------------------------------------
  ! zero energy variables
  !------------------------------------------------
  EKIN_N=zero;              EKIN_P=zero; 
  EVOL_rho_tau=zero;        EVOL_rho_rho=zero; 
  ESURF_rho_DELTA_rho=zero; ESURF_NABLA_rho_NABLA_rho=zero; 
  ESO_rho_NABLA_J=zero;     ESO_NABLA_rho_J=zero; E_HARTREE_DIR=zero
  ept1=zero; ept2=zero;     del1=zero; del2=zero; 
  ecodi=zero; ecoex=zero; etens=zero
  EExtra=zero ; E_EXT_FIELD = zero ;
  xn1=zero; xn2=zero; rms1=zero; rms2=zero
  q21=zero; q22=zero; q41=zero; q42=zero
  tempE_Crho0=zero; tempREARR=zero 
  !------------------------------------------------
  ! zero optimization variables
  !------------------------------------------------
  If(DO_FITT) Then
     efit_0=zero; efitV0=zero; dfitV0=zero
     efit_rhorho=zero; efit_rhorhoD=zero; 
     efit_rhotau=zero; efit_rhoDrho=zero; 
     efit_rhonablaJ=zero; efit_JJ=zero; 
  Endif
  !------------------------------------------------
  ! Integration in coordinate space
  !------------------------------------------------
  Do ihli=1,nghl
     whl=wdcor(ihli)
     !------------------------------------------------
     ! np-representation
     !------------------------------------------------
     rn=ro(ihli,1);      rp=ro(ihli,2); rnp2=rn**2+rp**2; rnp1=rn-rp
     tnt=tau(ihli,1);    tpt=tau(ihli,2); 
     dn=dro(ihli,1);     dp=dro(ihli,2); 
     djn=dj(ihli,1);     djp=dj(ihli,2); 
     akn=aka(ihli,1);    akp=aka(ihli,2)
     akn2=akn*akn;       akp2=akp*akp
     adn=akn*rn;         adp=akp*rp
     SFIZN=SFIZ(IHLI,1); SFIZP=SFIZ(IHLI,2); 
     SFIRN=SFIR(IHLI,1); SFIRP=SFIR(IHLI,2); 
     SZFIN=SZFI(IHLI,1); SZFIP=SZFI(IHLI,2); 
     SRFIN=SRFI(IHLI,1); SRFIP=SRFI(IHLI,2); 
     SNABLARN=NABLAR(IHLI,1); SNABLARP=NABLAR(IHLI,2); 
     SNABLAZN=NABLAZ(IHLI,1); SNABLAZP=NABLAZ(IHLI,2); 
     !------------------------------------------------
     ! t-representation
     !------------------------------------------------
     RHO_0=rn+rp;        RHO_1=rn-rp;   
     TAU_0=tnt+tpt;      TAU_1=tnt-tpt; 
     DRHO_0=dn+dp;       DRHO_1=dn-dp;   
     DJ_0=djn+djp;       DJ_1=djn-djp;   
     SFIZ_0=SFIZN+SFIZP; SFIZ_1=SFIZN-SFIZP; 
     SFIR_0=SFIRN+SFIRP; SFIR_1=SFIRN-SFIRP; 
     SZFI_0=SZFIN+SZFIP; SZFI_1=SZFIN-SZFIP; 
     SRFI_0=SRFIN+SRFIP; SRFI_1=SRFIN-SRFIP; 
     SNABLAR_0=SNABLARN+SNABLARP; SNABLAR_1=SNABLARN-SNABLARP; 
     SNABLAZ_0=SNABLAZN+SNABLAZP; SNABLAZ_1=SNABLAZN-SNABLAZP; 
     J2_0=SFIZ_0**2+SFIR_0**2+SZFI_0**2+SRFI_0**2
     J2_1=SFIZ_1**2+SFIR_1**2+SZFI_1**2+SRFI_1**2    
     !
     Call calculate_U_parameters(RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,  &
          (SNABLAR_0**2+SNABLAZ_0**2),(SNABLAR_1**2+SNABLAZ_1**2))
     !------------------------------------------------
     ! rms and deformations
     !------------------------------------------------
     z=fh(ihli); zz=z*z; rrr=zz+fl(ihli)**2
     p2=p32*zz-half*rrr; p3=p53*z*p2-p23*rrr*z; p4=p74*z*p3-p34*rrr*p2
     row=whl*rn; xn1=xn1+row; rms1=rms1+row*rrr; q21=q21+row*p2; q41=q41+row*p4
     row=whl*rp; xn2=xn2+row; rms2=rms2+row*rrr; q22=q22+row*p2; q42=q42+row*p4
     !------------------------------------------------
     ! PH energies
     !------------------------------------------------
     EKIN_N=EKIN_N+hb0*(TAU_0+TAU_1)*HALF*whl                                 ! kinetic, n
     EKIN_P=EKIN_P+hb0*(TAU_0-TAU_1)*HALF*whl                                 ! kinetic, p
     EVOL_rho_tau=EVOL_rho_tau+(Urhotau(0,0)*RHO_0*TAU_0  &                   ! volume rho tau
          +Urhotau(1,0)*RHO_1*TAU_1+Urhotau(2,0)*RHO_0*TAU_1  &
          +Urhotau(3,0)*RHO_1*TAU_0 )*whl
     EVOL_rho_rho=EVOL_rho_rho+(Urhorho(0,0)*RHO_0**2  &                      ! volume density dependent 
          +Urhorho(1,0)*RHO_1**2+(Urhorho(3,0)+Urhorho(2,0))*RHO_0*RHO_1)*whl        
     ESURF_rho_DELTA_rho =ESURF_rho_DELTA_rho+(UrhoDrho(0,0)*RHO_0*DRHO_0  &  ! surface: rho delta rho 
          +UrhoDrho(1,0)*RHO_1*DRHO_1+UrhoDrho(2,0)*RHO_0*DRHO_1  &
          +UrhoDrho(3,0)*RHO_1*DRHO_0 )*whl
     ESURF_NABLA_rho_NABLA_rho=ESURF_NABLA_rho_NABLA_rho  &                   ! surface: (nabla rho)**2 
          +(Unablarho(0,0)*(SNABLAR_0*SNABLAR_0+SNABLAZ_0*SNABLAZ_0)  &
          +Unablarho(1,0)*(SNABLAR_1*SNABLAR_1+SNABLAZ_1*SNABLAZ_1)  &
          +(Unablarho(3,0)+Unablarho(2,0))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1) )*whl
     ESO_rho_NABLA_J=ESO_rho_NABLA_J+(UrhonablaJ(0,0)*RHO_0*DJ_0  &           ! spin-orbit rho Nabla . J
          +UrhonablaJ(1,0)*RHO_1*DJ_1+UrhonablaJ(2,0)*RHO_0*DJ_1  &
          +UrhonablaJ(3,0)*RHO_1*DJ_0 )*whl
     ESO_NABLA_rho_J=ESO_NABLA_rho_J  &
          +(UJnablarho(0,0)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))  &  ! spin-orbit J . Nabla rho        
          +UJnablarho(1,0)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))  &
          +UJnablarho(2,0)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))  &
          +UJnablarho(3,0)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1)) )*whl        
     ETENS=ETENS+(UJJ(0,0)*J2_0+UJJ(1,0)*J2_1  &                              ! tensor J^2         
          +(UJJ(3,0)+UJJ(2,0))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1) )*whl
     EExtra=EExtra+(UEnonstdr(0)+UEnonstdr(1))*whl                            ! extra field if needed 
     E_EXT_FIELD=E_EXT_FIELD + ( Vexternal(0,zero,fl(ihli),z)*RHO_0 &         ! external field
          +Vexternal(1,zero,fl(ihli),z)*RHO_1 )*whl
     !------------------------------------------------
     ! Coulomb & Hartree
     !------------------------------------------------
     If (icou.Ge.1) ecodi=ecodi+half*cou(ihli)*rp*whl              
     If (icou.Eq.2) ecoex=ecoex-CExPar*cex*rp**p43*whl
     E_HARTREE_DIR=E_HARTREE_DIR +half*vDHartree(ihli,1)*RHO_0*whl+half*vDHartree(ihli,2)*RHO_1*whl
     ! just for printing
     tempE_Crho0=tempE_Crho0+RHO_0**2*whl
     tempREARR=tempREARR+(Cdrho(0)*RHO_0**2+Cdrho(1)*RHO_1**2)*RHO_0**sigma*whl
     !------------------------------------------------
     ! pairing energy and delta
     !------------------------------------------------
     rsa0=(RHO_0/rho_c) 
     dd1n=CpV0(0)*(ONE-rsa0*CpV1(0))*whl
     dd1p=CpV0(1)*(ONE-rsa0*CpV1(1))*whl
     ept1=ept1+dd1n*akn2; del1=del1-dd1n*adn
     ept2=ept2+dd1p*akp2; del2=del2-dd1p*adp          
     !------------------------------------------------
     ! optimization quatities
     !------------------------------------------------
     If(DO_FITT) Then
        efitV0(0)=efitV0(0)+(ONE-rsa0*CpV1(0))*akn2*whl
        efitV0(1)=efitV0(1)+(ONE-rsa0*CpV1(1))*akp2*whl
        dfitV0(0)=dfitV0(0)-(ONE-rsa0*CpV1(0))*adn*whl
        dfitV0(1)=dfitV0(1)-(ONE-rsa0*CpV1(1))*adp*whl
        !
        efit_rhotau(0)=efit_rhotau(0)+RHO_0*TAU_0*whl                                     ! rho tau 
        efit_rhotau(1)=efit_rhotau(1)+RHO_1*TAU_1*whl                                     ! rho tau 
        efit_rhorho(0)=efit_rhorho(0)+RHO_0**2*whl                                        ! rho^2   
        efit_rhorho(1)=efit_rhorho(1)+RHO_1**2*whl                                        ! rho^2 
        efit_rhorhoD(0)=efit_rhorhoD(0)+RHO_0**sigma*RHO_0**2*whl                         ! rho^2   
        efit_rhorhoD(1)=efit_rhorhoD(1)+RHO_0**sigma*RHO_1**2*whl                         ! rho^2   
        efit_rhoDrho(0)=efit_rhoDrho(0)+RHO_0*DRHO_0*whl                                  ! rho Delta rho 
        efit_rhoDrho(1)=efit_rhoDrho(1)+RHO_1*DRHO_1*whl                                  ! rho Delta rho 
        efit_rhonablaJ(0)=efit_rhonablaJ(0)+RHO_0*DJ_0*whl                                ! rho nabla J J 
        efit_rhonablaJ(1)=efit_rhonablaJ(1)+RHO_1*DJ_1*whl                                ! rho nabla J J 
        efit_JJ(0)=efit_JJ(0)+J2_0*whl                                                    ! J.J  
        efit_JJ(1)=efit_JJ(1)+J2_1*whl                                                    ! J.J 
     Endif
  End Do !ihli
  ! TMR pairing
  If(use_TMR_pairing.Ne.0) Then
	!ept1=ept1+dd1n*akn2; del1=del1-dd1n*adn
	!ept2=ept2+dd1p*akp2; del2=del2-dd1p*adp          
	ept1=e_pair(1); del1=d_pair(1)
	ept2=e_pair(2); del2=d_pair(2)
  Endif
  !------------------------------------------------
  ! after the integration
  !------------------------------------------------
  xn(1)=xn1;         xn(2)=xn2;         xn(3)=xn1+xn2; 
  rms(1)=rms1;       rms(2)=rms2 
  q2(1)=q21;         q2(2)=q22; 
  q4(1)=q41;         q4(2)=q42
  ekt(1)=EKIN_N;     ekt(2)=EKIN_P;     ekt(3)=ekt(1)+ekt(2) 
  ept(1)=ept1;       ept(2)=ept2;       ept(3)=ept(1)+ept(2) 
  del(1)=del1/xn(1); del(2)=del2/xn(2); 
  !
  EVOL=EVOL_rho_tau+EVOL_rho_rho+E_HARTREE_DIR        
  esurf=ESURF_rho_DELTA_rho+ESURF_NABLA_rho_NABLA_rho
  ESO=ESO_rho_NABLA_J+ESO_NABLA_rho_J
  ecoul=ecodi+ecoex
  etot=ekt(3)+evol+esurf+eso+ecoul+ept(3)+ETENS+EExtra+E_EXT_FIELD
  ehfb=etot
  !------------------------------------------------
  ! rms and deformations
  !------------------------------------------------
  Do it=itmin,itmax
     rms(it)=Sqrt(rms(it)/xn(it))
     q2(it)=two*q2(it)       !Qnp=<2r^2P_2(teta)>=<2z^2-x^2-y^2>
     q4(it)=ffdef4*q4(it)    !Hn=8r^4P_4(teta)=8z^4-24z^2(x^2+y^2)+3(x^2+y^2)^2
     def(it)=Sqrt(pi/5.0_pr)*q2(it)/(rms(it)**2*xn(it))
  End Do
  r212=rms(1)**2; r222=rms(2)**2
  rms(3)=Sqrt((xn(1)*r212+xn(2)*r222)/amas)
  q2(3)=q2(1)+q2(2)          ! quadrupole moment
  q4(3)=q4(1)+q4(2)          ! hexadecapole moment
  def(3)=Sqrt(pi/5.0_pr)*q2(3)/(rms(3)**2*amas) !deformation
  bet=def(3)
  !bet=ffdef6*q2(3)/(amas*r02)  ! bet=Q2*Sqrt(5 Pi)/(3A x^2);  x=r0 A^(1/3)
  !------------------------------------------------
  ! Lipkin-Nogami energy 
  !------------------------------------------------
  If(kindhfb.Lt.0) Then
     Call tracesln       
     If(ierror_flag.Ne.0) Return	   	               	 
     etot=etot+etr(3)
  End If
  !------------------------------------------------
  ! optimization quantities
  !------------------------------------------------
  If(DO_FITT) Then
     efV_0=0.0_pr
     If(kindhfb.Lt.0) Then
        efV_0(0)=ala2(1)
        efV_0(1)=ala2(2)
     End If
     dfitV0(0)=dfitV0(0)/xn(1)
     dfitV0(1)=dfitV0(1)/xn(2)
     efit_0=etot-efitV0(0)*CpV0(0)-efitV0(1)*CpV0(1)  &
          -efit_rhotau(0)*Ctau(0)-efit_rhotau(1)*Ctau(1)  &
          -efit_rhorho(0)*Crho(0)-efit_rhorho(1)*Crho(1)  &
          -efit_rhorhoD(0)*Cdrho(0)-efit_rhorhoD(1)*Cdrho(1)  &
          -efit_rhoDrho(0)*CrDr(0)-efit_rhoDrho(1)*CrDr(1)  &
          -efit_rhonablaJ(0)*CrdJ(0)-efit_rhonablaJ(1)*CrdJ(1)  &
          -efit_JJ(0)*CJ(0)-efit_JJ(1)*CJ(1)       
  Endif
  !------------------------------------------------
  ! debug
  !------------------------------------------------
  If(Print_Screen.And.IDEBUG.Gt.10) Then
     Write(*,'(4(a12,g13.6))')  &
          ' Tn=     ',ekt(1),           ' Tp=     ',ekt(2), &          
          ' EPn=    ',ept(1),           ' EPp=    ',ept(2),  &         
          ' EVOL=   ',EVOL,             ' Esurf=  ',esurf,  &           
          ' NrNr=   ',ESURF_NABLA_rho_NABLA_rho,' rDr=    ',ESURF_rho_DELTA_rho,  &
          ' Etens=  ',ETENS,            ' Eso=   ',eso,  &             
          ' rNJ=    ',ESO_rho_NABLA_J,  ' NrJ=   ',ESO_NABLA_rho_J,  &
          ' ECd=    ',ecodi,            ' ECex=  ',ecoex, &            
          ' EHd=    ',E_HARTREE_DIR,    ' Ir0^2= ',tempE_Crho0, &            			
          ' Eextra= ',EExtra,           ' Ext.Fl= ',E_EXT_FIELD, &
          ' Etot=  ',etot
     If(DO_FITT) Then
        Write(*,'(4(a12,g13.6))') 
        Write(*,'(4(a12,g13.6))')  &
             ' efrr0= ',efit_rhorho(0),     ' efrr1= ',efit_rhorho(1), &
             ' efrrD0=   ',efit_rhorhoD(0),          ' efrr1D=  ',efit_rhorhoD(1),  &
             ' efrt0= ',efit_rhotau(0),     ' efrt1= ',efit_rhotau(1), &   
             ' efrDr0=   ',efit_rhoDrho(0),          ' efrDr1=  ',efit_rhoDrho(1),  &
             ' efrDj0=',efit_rhonablaJ(0),  ' efrDj1=',efit_rhonablaJ(1), &
             ' efjj0=    ',efit_JJ(0),               ' efjj1=   ',efit_JJ(1),  &
             ' efV0_0=',efitV0(0),          ' efV0_1=',efitV0(1), &        
             ' dfV0_0=   ',dfitV0(0),                ' dfV0_1=  ',dfitV0(1),  &
             ' efV0=  ',efV_0(0),           ' efV_1= ',efV_0(1), &         
             ' ef0=      ',efit_0,                   ' etot=    ',etot 
     Endif
  Endif
  !------------------------------------------------
  ! Part called at the very end only (lpr=T)
  !------------------------------------------------
  If (lpr) Then
     !------------------------------------------------
     ! other definitions of deformations  (ffdef6=Sqrt(5.0_pr*pi)/3.0_pr)
     !------------------------------------------------
     bet2(1)=ffdef6*q2(1)/(xn(1)*r02) ! beta_n=Qn*Sqrt(5 Pi)/(3N x^2)
     bet2(2)=ffdef6*q2(2)/(xn(2)*r02) ! x=r0 A^(1/3)
     bet2(3)=ffdef6*q2(3)/(amas*r02)
     het4(1)=ffdef7*q4(1)/(xn(1)*r04)
     het4(2)=ffdef7*q4(2)/(xn(2)*r04)
     het4(3)=ffdef7*q4(3)/(amas*r04)
     rc=Sqrt(r222+0.640_pr)
     ! transitions to barn,barn^2,barn^4
     Do i=1,3
        q2(i)=q2(i)/100.0_pr; q4(i)=q4(i)/10000.0_pr
     End Do
     !------------------------------------------------
     ! STORE to unprojected buffer 'eresu'
     !------------------------------------------------
     ! ieresu=50 from module definitions
     ! ,'UEtot','Ubett','Ubetn','Ubetp',' UQt ',' UQn ',' UQp '  &
     eresu(1)=etot; eresu(2)=def(3); eresu(3)=def(1); eresu(4)=def(2); 
     eresu(5)=q2(3); eresu(6)=q2(1); eresu(7)=q2(2); 
     ! ,' Uln ',' Ulp ',' UpEn',' UpEp',' UpDn',' UpDp',' UAsn',' UAsp'  &
     eresu(8)=alast(1); eresu(9)=alast(2); eresu(10)=ept(1); eresu(11)=ept(2); 
     eresu(12)=del(1); eresu(13)=del(2); eresu(14)=ass(1); eresu(15)=ass(2); 
     ! ,' Urt ',' Urn ',' Urp ',' Urc ',' Uht ',' Uhn ',' Uhp '  &
     eresu(16)=rms(3); eresu(17)=rms(1); eresu(18)=rms(2); eresu(19)=rc; 
     eresu(20)=het4(3); eresu(21)=het4(1); eresu(22)=het4(2); 
     ! ,' Uqht',' Uqhn',' Uqhp'  &
     eresu(23)=q4(3); eresu(24)=q4(1); eresu(25)=q4(2); 
     ! ,'UKINT','UKINN','UKINP',' USO ','UCDIR',' UCEX','UDisn','UDisp'  &
     eresu(26)=ekt(3); eresu(27)=ekt(1); eresu(28)=ekt(2); eresu(29)=eso; 
     eresu(30)=ecodi; eresu(31)=ecoex; eresu(32)=Dispersion(1); eresu(33)=Dispersion(2); 
     ! ,'UV2Mn','UV2Mp'
     eresu(34)=v2min(1); eresu(35)=v2min(2); 
     !  ,'UECMT','UECMN','UECMP'
     eresu(36)=ECMHFB(3); eresu(37)=ECMHFB(1); eresu(38)=ECMHFB(2); 
     !  ,'UROTT','UROTN','UROTP'
     eresu(39)=DEROT(3); eresu(40)=DEROT(1); eresu(41)=DEROT(2); 
     !  ,'USQUJT','USQUJN','USQUJP'
     eresu(42)=SQUJ(3); eresu(43)=SQUJ(1); eresu(44)=SQUJ(2); 
     !  ,'UCRANT','UCRANN','UCRANP'
     eresu(45)=CRAN(3); eresu(46)=CRAN(1); eresu(47)=CRAN(2); 
     !  ,'UERIGT','UERIGN','UERIGP'
     eresu(48)=ERIGHFB(3); eresu(49)=ERIGHFB(1); eresu(50)=ERIGHFB(2); 
     !
     ! nucleus with wrong assymptotic
     If(iasswrong(3).Ne.0) eresu(16)=-eresu(16)
     !------------------------------------------------
     ! WRITE UNPROJECTED OUTPUT
     !------------------------------------------------
     Do iw=lout,lfile
        Write(iw,*)
        Write(iw,'(a,9x,a)')            '  NB! From expect (UNPROJECTED RESULTS)'
        Write(iw,*)
        If(iLST1.Ne.0)  &
             Write(iw,'(a,3f15.6)') '  hfb decay const. ass ',ass
        Write(iw,'(a,5f15.6)') '  pairing: CpV0,CpV1,...    ',CpV0,CpV1
        Write(iw,'(a,a)')      '  forces:   ',skyrme
        If(keyblo(1).Ne.0)  &
             Write(iw,'(a,i4,a,f10.3)')  '  Blocked neutron block    ',  &
             bloblo(keyblo(1),1)
        If(keyblo(2).Ne.0)  &
             Write(iw,'(a,i4,a,f10.3)')  '  Blocked proton  block    ',  &
             bloblo(keyblo(2),2)
        Write(iw,*)
        Write(iw,'(/,28x,a,8x,a,9x,a)') ' neutrons ','protons','total'
        Write(iw,'(a,6f15.6)')          '  Requested part.numbs.',tz,Sum(tz)
        Write(iw,'(a,6f15.6)')          '  UnPj(av) part.numbs .',xn
        Write(iw,'(a,3f15.6)')          '  b0, bz, bp ..........',b0,bz,bp
        Write(iw,*)
        Write(iw,'(a,3f15.6)') '  lambda (ala) ........',ala
        Write(iw,'(a,3f15.6)') '  Lambda (alast) ......',alast
        Write(iw,'(a,3f15.6)') '  delta(n,p), pwi .....',del,pwi
        Write(iw,'(a,3f15.6)') '  pairing energy ......',ept
        If(kindhfb.Lt.0) Then
           Write(iw,'(a,3f15.6)') '  LN lambda_2 ... ala2 ',ala2
           Write(iw,'(a,3f15.6)') '  LN energies .........',etr
           Write(iw,'(a,3f15.6)') '  delta(n,p)+ala2 .....',del+ala2
           Write(iw,'(a,3f15.6)') '  Geff(n,p) ...........',Geff
        End If
        Write(iw,*)
        Write(iw,'(a,3f15.6)') '  rms-radius ..........',rms
        Write(iw,'(a,15x,2f15.6)') '  charge-radius, r0 ...',rc,r00
        Write(iw,'(a,3f15.6)') '  deformation ......bet',bet
        Write(iw,'(a,3f15.6)') '  q2 defformation, beta',bet2
        Write(iw,'(a,3f15.6)') '  q4 defformation, heta',het4
        Write(iw,'(a,3f15.6)') '  deformation .........',def
        Write(iw,'(a,3f15.6)') '  quadrupole moment[b] ',q2
        Write(iw,'(a,3f15.6)') '  hexadecapole moment .',q4
        Write(iw,'(a,3f15.6)') '  <q20> ...............',  &
             (q2(i)*ffdef3,i=1,3)
        Write(iw,'(a,3f15.6)') '  <q40> ...............',  &
             (q4(i)*ffdef5,i=1,3)
        Write(iw,*)
        Write(iw,'(a,3f15.6)')    '  kinetic energy ......',ekt
        Write(iw,'(a,30x,f15.6)') '  volume energy .......',evol
        Write(iw,'(a,30x,f15.6)') '        rho_tau .......',EVOL_rho_tau
        Write(iw,'(a,30x,f15.6)') '        rho_rho .......',EVOL_rho_rho
        Write(iw,'(a,30x,f15.6)') '  surface energy ......',esurf
        Write(iw,'(a,30x,f15.6)') '   rho_DELTA_rho ......',ESURF_rho_DELTA_rho
        Write(iw,'(a,30x,f15.6)') '   (NABLA_rho)^2 ......',ESURF_NABLA_rho_NABLA_rho          
        Write(iw,'(a,30x,f15.6)') '  spin-orbit energy ...',eso
        Write(iw,'(a,30x,f15.6)') '        rho_NABLA_J ...',ESO_rho_NABLA_J
        Write(iw,'(a,30x,f15.6)') '        NABLA_rho_J ...',ESO_NABLA_rho_J
        Write(iw,'(a,30x,f15.6)') '  coulomb energy ......',ecodi+ecoex
        Write(iw,'(a,30x,f15.6)') '          direct ......',ecodi
        Write(iw,'(a,30x,f15.6)') '          exchange ....',ecoex
        Write(iw,'(a,30x,f15.6)') '  tensor energy .......',etens
        Write(iw,'(a,30x,f15.6)') '  direct Hartree E  ...',E_HARTREE_DIR
        Write(iw,'(a,30x,f15.6)') '  Extra E .............',EEXTRA
        Write(iw,'(a,30x,f15.6)') '  External field E ....',E_EXT_FIELD
        Write(iw,*)
        Write(iw,'(a,30x,f15.6)')    '  tEnergy: ehfb (qp)...',ehfb
        If(kindhfb.Lt.0) Then
           Write(iw,'(a,30x,f15.6)') '  tEnergy: ehfb(qp)+LN ',etot
        End If
        Write(iw,*)
        Write(iw,'(a,6f15.6)')    '  Calculated but not added corrections '
        Write(iw,'(a,6f15.6)')    '===================================='
        Write(iw,'(a,6f15.6)')    '  cmc-diagonal part ...',ekt/hb0*hbzero-ekt
        Write(iw,'(a,6f15.6)')    '  cmc-hfb .............',ECMHFB
        Write(iw,'(a,6f15.6)')    '  cranking rot corr ...',DEROT
        Write(iw,*)
        Write(iw,'(a,6f15.6)')    '  SQUJ ................',SQUJ
        Write(iw,'(a,6f15.6)')    '  CRAN x 4 ............',4.0_pr*CRAN
        Write(iw,'(a,6f15.6)')    '  Rigit Body ..........',ERIGHFB
        Write(iw,'(a,6f15.6)')
     End Do
     !------------------------------------------------
     ! START corrected Lipkin-Nogami characteristics
     !------------------------------------------------
     If(kindhfb.Lt.0) Then
        Call densitln       !density LN corrections
        If(ierror_flag.Ne.0) Return	   	               		
        Do it=itmin,itmax
           xn(it)=zero
           rms(it)=zero; q2(it)=zero; q4(it)=zero
        End Do
        !
        Do ihli=1,nghl
           whl=wdcor(ihli)
           rn=ro(ihli,1); rp=ro(ihli,2); rnp2=rn**2+rp**2
           ! rms and deformations
           z=fh(ihli); zz=z*z; rrr=zz+fl(ihli)**2
           p2=p32*zz   -half*rrr
           p3=p53*z*p2 -p23*rrr*z
           p4=p74*z*p3 -p34*rrr*p2
           row=whl*rn
           xn(1)=xn(1)+row
           rms(1)=rms(1)+row*rrr
           q2(1)=q2(1)+row*p2
           q4(1)=q4(1)+row*p4
           row=whl*rp
           xn(2)=xn(2)+row
           rms(2)=rms(2)+row*rrr
           q2(2)=q2(2)+row*p2
           q4(2)=q4(2)+row*p4
        End Do !ihli
        !------------------------------------------------
        ! rms and deformations
        !------------------------------------------------
        Do it=itmin,itmax
           rms(it)=Sqrt(rms(it)/xn(it))
           q2(it)=two*q2(it)       !Qnp=<2r^2P_2(teta)>=<2z^2-x^2-y^2>
           q4(it)=ffdef4*q4(it)  !Hn=<8r^4P_4(teta)>=<8z^4-24z^2(x^2+y^2)+3(x^2+y^2)^2>
           def(it)=Sqrt(pi/5.0_pr)*q2(it)/(rms(it)**2*xn(it))
        End Do
        r212=rms(1)**2; r222=rms(2)**2
        rms(3)=Sqrt((xn(1)*r212+xn(2)*r222)/amas)
        q2(3)=q2(1)+q2(2)  ! quadrupole moment
        q4(3)=q4(1)+q4(2)  ! hexadecapole moment
        def(3)=Sqrt(pi/5.0_pr)*q2(3)/(rms(3)**2*amas) !deformation
        ! other definitions of the same quantitsies
        bet2(1)=ffdef6*q2(1)/(xn(1)*r02) !beta_n=Qn*Sqrt(5Pi)/(3N x^2)
        bet2(2)=ffdef6*q2(2)/(xn(2)*r02) !x=r0=1.2A^(1/3)
        bet2(3)=ffdef6*q2(3)/(amas*r02)
        het4(1)=ffdef7*q4(1)/(xn(1)*r04)
        het4(2)=ffdef7*q4(2)/(xn(2)*r04)
        het4(3)=ffdef7*q4(3)/(amas*r04)
        xn(3)=xn(1)+xn(2)
        bet=def(3)
        rc=Sqrt(r222+0.640_pr)
        ! transitions to barn,barn^2,barn^4
        Do i=1,3
           q2(i)=q2(i)/100.0_pr; q4(i)=q4(i)/10000.0_pr
        End Do
        !------------------------------------------------
        ! STORE to unprojected LN buffer 'eresl'
        !------------------------------------------------
        ! ieresl=20 from module definitions
        ! ,' EHFBLN',' EHFB',' LNEt','LNbet','LNben','LNbep',' LNQt',' LNQn',' LNQp'  &
        eresl(1)=etot; eresl(2)=etot-etr(3); 
        eresl(3)=def(3); eresl(4)=def(1); eresl(5)=def(2)
        eresl(6)=q2(3); eresl(7)=q2(1); eresl(8)=q2(2); 
        ! ,'LNpEn','LNpEp','LNpDn','LNpDp',' LNrt',' LNrn',' LNrC'  &
        eresl(9)=ept(1); eresl(10)=ept(2); eresl(11)=del(1)+ala2(1); eresl(12)=del(2)+ala2(2); 
        eresl(13)=rms(3); eresl(14)=rms(1); eresl(15)=rms(2); eresl(16)=rc; 
        ! ,' LNam2n',' LNam2p',' LNe2n',' LNe2p'
        eresl(17)=ala2(1); eresl(18)=ala2(2); eresl(19)=etr(1); eresl(20)=etr(2)
        !------------------------------------------------
        ! WRITE UNPROJECTED LN OUTPUT 
        !------------------------------------------------
        Do iw=lout,lfile
           Write(iw,'(a,3f15.6)')
           Write(iw,'(a,3f15.6)') '  With Lipkin-Nogami Corrections'
           Write(iw,'(a,3f15.6)') '================================'
           Write(iw,'(a,3f15.6)') '  rms-radius ..........',rms
           Write(iw,'(a,15x,2f15.6)') '  charge-radius, r0 ...',rc,r00
           Write(iw,'(a,3f15.6)') '  deformation ......bet',bet
           Write(iw,'(a,3f15.6)') '  q2 deformation,  beta',bet2
           Write(iw,'(a,3f15.6)') '  q4 deformation,  heta',het4
           Write(iw,'(a,3f15.6)') '  deformation .........',def
           Write(iw,'(a,3f15.6)') '  quadrupole moment[b] ',q2
           Write(iw,'(a,3f15.6)') '  hexadecapole moment .',q4
           Write(iw,'(a,3f15.6)') '  <q20> ...............',  &
                (q2(i)*ffdef3,i=1,3)
           Write(iw,'(a,3f15.6)') '  <q40> ...............',  &
                (q4(i)*ffdef5,i=1,3)
           Write(iw,'(a,3f15.6)') '================================'
           Write(iw,'(a,3f15.6)')
        End Do
     End If
     !------------------------------------------------
     ! WRITE all blocking candidates
     !------------------------------------------------
     If(keyblo(3).Eq.0) Then
        Do iw=lout,lfile
           Write(iw,*)
           Do it=itmin,itmax
              If(it.Eq.1) Then
                 EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN
              Else
                 EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP
              Endif
              !
              Write(iw,*) ' ',' Blocking candidates are:'
              Write(iw,*) '  ',protn(it),' eqpmin=',eqpmin(it),' pwiblo=',pwiblo
              Do i=1,blomax(it)
                 Write(iw,'(a,i4,a,i4,a,i4,2x,i4,3(a,1x,f12.8,1x),a)') '    num=',i,  &
                      ' block=',bloblo(i,it),  &
                      ' state=',blo123(i,it),blok1k2(i,it),  &
                      ' Eqp=',EqpPo(KqpPo(blok1k2(i,it))),  &
                      ' (1-2N)E=',(one-two*uk(blok1k2(i,it),it))*EqpPo(KqpPo(blok1k2(i,it))),  &
                      ' Ovlp=',vkmax(blok1k2(i,it),it),  &
                      tb(numax(blok1k2(i,it),it))
              End Do
              Write(iw,*)
           End Do
           Write(iw,*)
        End Do
     End If
     !
     !------------------------------------------------
     ! PAV
     !------------------------------------------------
     ! Projecting on different nucleus
     If(iproj.Ne.0) Then
        npr(1)=npr1pj;       npr(2)=npr2pj
        tz(1)=Dble(npr(1)); tz(2)=Dble(npr(2))
     End If
     !Call expectpj(.True.)
  End If
  !
  Call get_CPU_time('expect',1)
  !
End Subroutine expect
!=======================================================================
!
!=======================================================================
Subroutine expectpj(lpr)
  !---------------------------------------------------------------------
  ! calculates expectation values (tz is the particle number)
  ! optimized for half gauge-angle points
  !---------------------------------------------------------------------
  Use HFBTHO
  Implicit None
  Logical :: lpr
  Integer(ipr) :: i,j,it,ihli,iw,iw1=901,iw2=902
  Complex(pr) :: SZFIN,SFIZN,SRFIN,SFIRN,SZFIP,SFIZP,SRFIP,SFIRP
  Complex(pr) :: SFIZ_0,SFIR_0,SZFI_0,SRFI_0,SFIZ_1,SFIR_1,SZFI_1,SRFI_1
  Complex(pr) :: RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,DJ_0,DJ_1
  Complex(pr) :: SNABLAR_0,SNABLAZ_0,SNABLAR_1,SNABLAZ_1
  Complex(pr) :: cekt(3),cdel(2),cept(3),cetot,etens,cq2pj(ilpj,ilpj)
  Complex(pr) :: cxn(2),crms(3),cq2(3),cq4(3),xnpj(2),rmspj(3),q2pj(3),q4pj(3)
  Complex(pr) :: evolpj,esurpj,ecdipj,ecexpj,ecoupj,ept1pj,ept2pj,epotpj  &
       ,eki1pj,eki2pj,ekinpj,etotpj,espopj,epa1pj,epa2pj,epirpj,ede1pj,ede2pj,etenspj  &
       ,eva,ev5,es5,eso,ecodi,ecoex,rn,rp,rnp1,rnp2,rt,rt2,tnt,tpt,tt  &
       ,dn,dp,dt,akn,akp,akn2,akp2,adn,adp,evol,esurf,ecoul,pijk,row,cx,dd1n,dd1p  &
       ,rt1,tt1,dt1,djn,djp,djt,djt1
  Complex(pr) :: rsa,rsa0,rsa0A,rps,rns,rsa1,rsa10,rsa12,rsa0An,rsa0As
  Real(pr) :: whl,x,xn(3),q4(3),def(3),bet2(3),het4(3),r212,r222,rc,z,zz,rrr,p2,p3,p4
  Real(pr)    :: rdelta(2),repair(3),rekin(3),revolpj,resurpj,respopj,recdipj,recexpj,retenspj
  !
  !=====================================================================
  ! IT IS NOT CLEAR HOW TO WORK WITN COMPLEX DENSITOES IN UNEDF
  ! TEMPORARY THIS PART SKIPED
  !=====================================================================
  !
  ! TEMPORARY SKIP JUMP (remove when ready)
  GO TO 200
  !
  Call densitpj ! calculates complex densities and the direct coulomb field
  If(ierror_flag.Ne.0) Return	   	                 
  !
  evolpj=zero; esurpj=zero; ecdipj=zero; ecexpj=zero; espopj=zero; 
  ept1pj=zero; ept2pj=zero; eki1pj=zero; eki2pj=zero; epj=zero; 
  etotpj=zero; epa1pj=zero; epa2pj=zero; ede1pj=zero; ede2pj=zero; 
  xnpj=zero; rmspj=zero; q2pj=zero; q4pj=zero; etenspj=zero; cq2pj=zero
  !
  Do i=1,ilpj
     Do j=1,ilpj
        pijk=pjk(i,1)*pjk(j,2)
        !
        cekt=zero;  cept=zero;  cdel=zero; 
        cxn=zero;  crms=zero;  cq2=zero; cq4=zero; 
        eva=zero;   EVOL=zero; ev5=zero; es5=zero; 
        eso=zero; ecodi=zero; ecoex=zero; etens=zero; 
        !
        Do ihli=1,nghl
           ! real
           whl=wdcor(ihli)
           !
           ! q-representation
           !
           !rn=ropj(ihli,i,1); rp=ropj(ihli,j,2); rnp2=rn**2+rp**2; rnp1=rn-rp
           !rt=rn+rp; rt2=rt*rt; rt1=rn-rp; 
           !tnt=taupj(ihli,i,1); tpt=taupj(ihli,j,2); tt=tnt+tpt; tt1=tnt-tpt
           !dn=dropj(ihli,i,1); dp=dropj(ihli,j,2); dt=dn+dp;   dt1=dn-dp
           !djn=djpj(ihli,i,1);  djp=djpj(ihli,j,2); djt=djn+djp; djt1=djn-djp
           !akn=akapj(ihli,i,1); akp=akapj(ihli,j,2)
           !akn2=akn*akn;         akp2=akp*akp
           !adn=akn*rn;          adp=akp*rp
           !SFIZN=SFIZpj(IHLI,i,1); SFIZP=SFIZpj(IHLI,j,2)
           !SFIRN=SFIRpj(IHLI,i,1); SFIRP=SFIRpj(IHLI,j,2)
           !SZFIN=SZFIpj(IHLI,i,1); SZFIP=SZFIpj(IHLI,j,2)
           !SRFIN=SRFIpj(IHLI,i,1); SRFIP=SRFIpj(IHLI,j,2)             
           !SNABLARN=NABLAR(IHLI,i,1); SNABLARP=NABLAR(IHLI,j,2); 
           !SNABLAZN=NABLAZ(IHLI,i,1); SNABLAZP=NABLAZ(IHLI,j,2);                                       
           !
           ! t-representation
           !
           !RHO_0=rn+rp;   RHO_1=rn-rp;   
           !TAU_0=tnt+tpt; TAU_1=tnt-tpt; 
           !DRHO_0=dn+dp;  DRHO_1=dn-dp;   
           !DJ_0=djn+djp;  DJ_1=djn-djp;   
           !SFIZ_0=SFIZN+SFIZP; SFIZ_1=SFIZN-SFIZP; 
           !SFIR_0=SFIRN+SFIRP; SFIR_1=SFIRN-SFIRP; 
           !SZFI_0=SZFIN+SZFIP; SZFI_1=SZFIN-SZFIP; 
           !SRFI_0=SRFIN+SRFIP; SRFI_1=SRFIN-SRFIP; 
           !SNABLAR_0=SNABLARN+SNABLARP; SNABLAR_1=SNABLARN-SNABLARP; 
           !SNABLAZ_0=SNABLAZN+SNABLAZP; SNABLAZ_1=SNABLAZN-SNABLAZP; 
           !J2_0=SFIZ_0**2+SFIR_0**2+SZFI_0**2+SRFI_0**2
           !J2_1=SFIZ_1**2+SFIR_1**2+SZFI_1**2+SRFI_1**2    
           !
           !Call calculate_U_parameters(RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1, &
           !       (SNABLAR_0**2+SNABLAZ_0**2),(SNABLAR_1**2+SNABLAZ_1**2)) 
           ! 
           z=fh(ihli); zz=z*z; rrr=zz+fl(ihli)**2
           p2=p32*zz   -half*rrr    !3/2 z*z-1/2 (z*z+r*r)=1/2(2 z*z-r*2)=1/2 Q
           p3=p53*z*p2 -p23*rrr*z
           p4=p74*z*p3 -p34*rrr*p2
           !complex
           !ig-particle number, rms and deformations
           row=whl*rn; cxn(1)=cxn(1)+row; crms(1)=crms(1)+row*rrr; cq2(1)=cq2(1)+row*p2; cq4(1)=cq4(1)+row*p4
           row=whl*rp; cxn(2)=cxn(2)+row; crms(2)=crms(2)+row*rrr; cq2(2)=cq2(2)+row*p2; cq4(2)=cq4(2)+row*p4
           !ig-energy contributions
           !
           ! ig-Pairing energy and delta
           rsa0=(rt/rho_c) 
           dd1n=CpV0(0)*(ONE-rsa0*CpV1(0))*whl
           dd1p=CpV0(1)*(ONE-rsa0*CpV1(1))*whl
           !
           cept(1)=cept(1)+dd1n*akn2; cept(2)=cept(2)+dd1p*akp2 
           cdel(1)=cdel(1)-dd1n*adn;  cdel(2)=cdel(2)-dd1p*adp  
           !
           x=hb0*whl

           !
           ! PH energy 
           !
           !             EKIN_N=EKIN_N+hb0*(TAU_0+TAU_1)*HALF*whl                                ! kinetic, n
           !             EKIN_P=EKIN_P+hb0*(TAU_0-TAU_1)*HALF*whl                                ! kinetic, p
           !             EVOL_rho_tau=EVOL_rho_tau+(Urhotau(0,0)*RHO_0*TAU_0  &  ! volume rho tau
           !+Urhotau(1,0)*RHO_1*TAU_1+Urhotau(2,0)*RHO_0*TAU_1  &
           !+Urhotau(3,0)*RHO_1*TAU_0 )*whl
           !             EVOL_rho_rho=EVOL_rho_rho+(Urhorho(0,0)*RHO_0**2  &  ! volume density dependent 
           !+Urhorho(1,0)*RHO_1**2+(Urhorho(2,0)+Urhorho(3,0))  &
           !              *RHO_0*RHO_1 )*whl
           !             ESURF_rho_DELTA_rho=ESURF_rho_DELTA_rho+(UrhoDrho(0,0)*RHO_0  &  ! surface: rho delta rho 
           !              *DRHO_0+UrhoDrho(1,0)*RHO_1*DRHO_1+UrhoDrho(2,0)*RHO_0  &
           !              *DRHO_1+UrhoDrho(3,0)*RHO_1*DRHO_0 )*whl
           !             ESURF_NABLA_rho_NABLA_rho=ESURF_NABLA_rho_NABLA_rho  &  ! surface: (nabla rho)**2 
           !+(Unablarho(0,0)*(SNABLAR_0*SNABLAR_0+SNABLAZ_0*SNABLAZ_0)  &
           !+Unablarho(1,0)*(SNABLAR_1*SNABLAR_1+SNABLAZ_1*SNABLAZ_1)  &
           !+(Unablarho(2,0)+Unablarho(3,0))*(SNABLAR_0*SNABLAR_1  &
           !+SNABLAZ_0*SNABLAZ_1) )*whl
           !             ESO_rho_NABLA_J=ESO_rho_NABLA_J+(UrhonablaJ(0,0)*RHO_0*DJ_0  &  ! spin-orbit rho Nabla . J
           !+UrhonablaJ(1,0)*RHO_1*DJ_1+UrhonablaJ(2,0)*RHO_0*DJ_1  &
           !+UrhonablaJ(3,0)*RHO_1*DJ_0 )*whl
           !             ESO_NABLA_rho_J=ESO_NABLA_rho_J+( UJnablarho(0,0)*(SNABLAR_0  &  ! spin-orbit J . Nabla rho
           !              *(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))+UJnablarho(1,0)  &
           !              *(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))  &
           !+UJnablarho(2,0)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1  &
           !              *(SFIR_0-SRFI_0))+UJnablarho(3,0)*(SNABLAR_0*(SFIZ_1-SZFI_1)  &
           !              -SNABLAZ_0*(SFIR_1-SRFI_1)) )*whl
           !             ETENS=ETENS+(UJJ(0,0)*(SFIZ_0**2+SFIR_0**2+SZFI_0**2+SRFI_0**2)  &
           !+UJJ(1,0)*(SFIZ_1**2+SFIR_1**2+SZFI_1**2+SRFI_1**2)  &
           !+(UJJ(2,0)+UJJ(3,0))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1  &
           !+SZFI_0*SZFI_1+SRFI_0*SRFI_1) )*whl

           !             cekt(1)=cekt(1)+x*tnt; cekt(2)=cekt(2)+x*tpt       !kinetic
           !             EVOL=EVOL+(U_rt_0*rt*tt+U_rt_1*rt1*tt1)*whl        ! volume
           !             EVOL=EVOL+((U_rr_00+U_rr_0D*rt**alph)*rt**2+(U_rr_10+U_rr_1D*rt**alph)*rt1**2)*whl 
           !             es5=es5+(U_rDr_0*rt*dt+U_rDr_1*rt1*dt1)*whl      ! surface  
           !             eso=eso+(U_rNj_0*rt*djt+U_rNj_1*rt1*djt1)*whl        ! spin-orbit
           !
           If (icou.Ge.1) ecodi=ecodi+half*coupj(ihli,j)*rp*whl   !Coul.dir
           If (icou.Eq.2) ecoex=ecoex-CExPar*cex*rp**p43*whl            !Coul.exc
        End Do !ihli
        !
        esurf=es5; ecoul=ecodi+ecoex
        cekt(3)=cekt(1)+cekt(2); cept(3)=cept(1)+cept(2)
        cetot=cekt(3)+evol+esurf+eso+ecoul+cept(3)+ ETENS
        cdel(1)=cdel(1)/tz(1); cdel(2)=cdel(2)/tz(2)
        !
        ! half-projected energies required for the matrix elements
        epj(i,1)=epj(i,1)+cetot*pjk(j,2)
        epj(j,2)=epj(j,2)+cetot*pjk(i,1)
        !
        ! for constraint contributions to half-projected energies  
        If (icstr.Ne.0) cq2pj(i,j)=two*(cq2(1)+cq2(2))
        !
        ! projected energies
        evolpj=evolpj+pijk*evol; esurpj=esurpj+pijk*esurf; espopj=espopj+pijk*eso
        epa1pj=epa1pj+pijk*cept(1); epa2pj=epa2pj+pijk*cept(2); epirpj=epa1pj+epa2pj
        ede1pj=ede1pj+pijk*cdel(1); ede2pj=ede2pj+pijk*cdel(2)
        ecdipj=ecdipj+pijk*ecodi;   ecexpj=ecexpj+pijk*ecoex;   ecoupj=ecdipj+ecexpj
        ept1pj=ept1pj+pijk*cept(1); ept2pj=ept2pj+pijk*cept(2); epotpj=ept1pj+ept2pj
        eki1pj=eki1pj+pijk*cekt(1); eki2pj=eki2pj+pijk*cekt(2); ekinpj=eki1pj+eki2pj
        !
        etotpj=etotpj+pijk*cetot
        etenspj=etenspj+ pijk*etens
        !
        ! unprojected hfb total energy and constraint
        If(i.Eq.1.And.j.Eq.1) Then
           rehfbcan=Real(cetot) 
        Endif
        !
        ! projected particle numbers, rms, deformations
        If(j.Eq.1) Then
           xnpj(1)=xnpj(1)+pjk(i,1)*cxn(1)
           rmspj(1)=rmspj(1)+pjk(i,1)*crms(1)
           q2pj(1)=q2pj(1)+pjk(i,1)*cq2(1)
           q4pj(1)=q4pj(1)+pjk(i,1)*cq4(1)
        Endif
        If(i.Eq.1) Then
           xnpj(2)=xnpj(2)+pjk(j,2)*cxn(2)
           rmspj(2)=rmspj(2)+pjk(j,2)*crms(2)
           q2pj(2)=q2pj(2)+pjk(j,2)*cq2(2)
           q4pj(2)=q4pj(2)+pjk(j,2)*cq4(2)
        Endif
        !
     End Do !j
  End Do !i
  !
  ! Real quantities to the end
  !
  ! Energies
  rdelta(1)=Real(ede1pj); rdelta(2)=Real(ede2pj); retotpj=Real(etotpj); 
  repair(1)=Real(epa1pj); repair(2)=Real(epa2pj); repair(3)=Real(epirpj)
  rekin(1)=Real(eki1pj);  rekin(2)=Real(eki2pj); rekin(3)=Real(ekinpj)
  revolpj=Real(evolpj);   resurpj=Real(esurpj); respopj=Real(espopj); 
  recdipj=Real(ecdipj);   recexpj=Real(ecexpj); retenspj=Real(etenspj); 
  depnp=retotpj-rehfbcan  !correlation energy due to projection
  !
  ! rms and deformations
  Do it=itmin,itmax
     xn(it)=Real(xnpj(it))
     rms(it)=Sqrt(Real(rmspj(it))/xn(it))
     q2(it)=two*Real(q2pj(it))       !Qnp=<2r^2P_2(teta)>=<2z^2-x^2-y^2>
     q4(it)=ffdef4*Real(q4pj(it))  !Hn=<8r^4P_4(teta)>=<8z^4-24z^2(x^2+y^2)+3(x^2+y^2)^2>
     def(it)=Sqrt(pi/5.0_pr)*q2(it)/(rms(it)**2*xn(it))
  End Do
  r212=rms(1)**2; r222=rms(2)**2
  rms(3)=Sqrt((xn(1)*r212+xn(2)*r222)/amas)
  q2(3)=q2(1)+q2(2)  ! quadrupole moment
  q4(3)=q4(1)+q4(2)  ! hexadecapole moment
  def(3)=Sqrt(pi/5.0_pr)*q2(3)/(rms(3)**2*amas) !deformation
  ! other definitions of the same quantities
  bet2(1)=ffdef6*q2(1)/(xn(1)*r02) !beta_n=Qn*Sqrt(5Pi)/(3N x^2)
  bet2(2)=ffdef6*q2(2)/(xn(2)*r02) !x=r0=1.2A^(1/3)
  bet2(3)=ffdef6*q2(3)/(amas*r02)
  het4(1)=ffdef7*q4(1)/(xn(1)*r04)
  het4(2)=ffdef7*q4(2)/(xn(2)*r04)
  het4(3)=ffdef7*q4(3)/(amas*r04)
  xn(3)=xn(1)+xn(2)
  bet=def(3)
  !
  !  constraint constants and contributions to half-projected energies  
  If (icstr.Ne.0) Then
     ty20=Sqrt(5.0_pr/pi)*hom/b0**2/two   
     cx=cqad*(cdef-bet)*ty20; 
     Do i=1,ilpj
        Do j=1,ilpj
           epj(i,1)=epj(i,1)+cx*cq2pj(i,j)*pjk(j,2)
           epj(j,2)=epj(j,2)+cx*cq2pj(i,j)*pjk(i,1)
        End Do
     End Do
  Endif
  !
  ! TEMPORARY SKIP LABLE (remove when done)
200 Continue    
  !
  !
  If (lpr) Then
     rc=Sqrt(r222+0.640_pr)
     ! transitions to barn,barn^2,barn^4
     Do i=1,3
        q2(i)=q2(i)/100.0_pr; q4(i)=q4(i)/10000.0_pr
     End Do
     !
     ! STORE to projected buffer 'eresj'
     ! ieresj=50 from module definitions
     ! ' si ','JININ'
     eresj(1)=si; eresj(2)=inin; 
     ! ' A','   N ','   Z ' 
     eresj(3)=npr(1)+npr(2); eresj(4)=npr(1); eresj(5)=npr(2); 
     ! ' Jln ',' Jlp '
     eresj(6)=alast(1); eresj(7)=alast(2); 
     ! ,'JEtot','Jbett','Jbetn','Jbetp',' JQt ',' JQn ',' JQp '  &
     eresj(8)=retotpj; eresj(9)=def(3); eresj(10)=def(1); eresj(11)=def(2); 
     eresj(12)=q2(3); eresj(13)=q2(1); eresj(14)=q2(2); 
     ! ' JpEn',' JpEp',' JpDn',' JpDp',' JAsn',' JAsp'  &
     eresj(15)=repair(1); eresj(16)=repair(2); 
     eresj(17)=rdelta(1); eresj(18)=rdelta(2); eresj(19)=ass(1); eresj(20)=ass(2); 
     ! ,' Jrt ',' Jrn ',' Jrp ',' Jrc ',' Jht ',' Jhn ',' Jhp '  &
     eresj(21)=rms(3); eresj(22)=rms(1); eresj(23)=rms(2); eresj(24)=rc; 
     eresj(25)=het4(3); eresj(26)=het4(1); eresj(27)=het4(2); 
     ! ,' Jqht',' Jqhn',' Jqhp'  &
     eresj(28)=q4(3); eresj(29)=q4(1); eresj(30)=q4(2); 
     ! ,' JKINt',' JKINn','JKINp',' JSO ','JCDIR',' JCEX','JDisn','JDisp'  &
     eresj(31)=rekin(3); eresj(32)=rekin(1); eresj(33)=rekin(2); eresj(34)=respopj; 
     eresj(35)=recdipj; eresj(36)=recexpj; eresj(37)=Dispersion(1); eresj(38)=Dispersion(2); 
     ! ,'JV2Mn','JV2Mp','JILST','JKIND','  JL '  &
     eresj(39)=v2min(1); eresj(40)=v2min(2)
     eresj(41)=iLST; eresj(42)=kindhfb; eresj(43)=iLpj; 
     !  ,'JECMPAV','JECMPAV','JECMPAV'
     eresj(44)=ECMPAV(3); eresj(45)=ECMPAV(1); eresj(46)=ECMPAV(2); 
     ! 'JA','JN',JZ'
     eresj(47)=Nint(xn(3)); eresj(48)=Nint(xn(1)); eresj(49)=Nint(xn(2)); 
     ! 'iter'
     eresj(50)=iiter
     ! nucleus with wrong asymptotic
     If(iasswrong(3).Ne.0) eresj(21)=-eresj(21)
     !
     ! WRITE to screen 'lout' and tape akzout.dat 'lfile'
     Do iw=lout,lfile
        Write(iw,*)
        Write(iw,'(a,9x,a,/)')            '  NB! From expectpj (PNP PAV RESULTS)'
        Write(iw,*)
        If(iLST1.Ne.0)  &
             Write(iw,'(a,6f15.6)') '  hfb decay const. ass ',ass
        Write(iw,'(a,8f15.6)') '  pairing: CpV0,CpV1,pwi... ',CpV0,CpV1,pwi
        Write(iw,'(a,a,a,i3)') '  forces:   ',skyrme,',  Gauge points:',ilpj
        If(keyblo(1).Ne.0)  &
             Write(iw,'(a,i4,a,f10.3)')  '  Blocked neutron block    ',  &
             bloblo(keyblo(1),1)
        If(keyblo(2).Ne.0)  &
             Write(iw,'(a,i4,a,f10.3)')  '  Blocked proton  block    ',  &
             bloblo(keyblo(2),2)
        Write(iw,*)
        Write(iw,'(/,28x,a,8x,a,9x,a)') ' neutrons ','protons','total'
        Write(iw,'(a,6f15.6)') '  Requested part.numbs.',tz,Sum(tz)
        Write(iw,'(a,6f15.6)') '  Projected part.numbs.',xn
        Write(iw,'(a,3f15.6)') '  Dispersion dN2 ......',Dispersion
        Write(iw,'(a,6f15.6)') '  b0, bz, bp ..........',b0,bz,bp
        Write(iw,*)
        Write(iw,'(a,6f15.6)') '  lambda (ala) ........',ala
        Write(iw,'(a,6f15.6)') '  Lambda (alast) ......',alast
        Write(iw,'(a,6f15.6)') '  delta(n,p) ..........',rdelta
        Write(iw,'(a,6f15.6)') '  pairing energy ......',repair
        Write(iw,*)
        Write(iw,'(a,6f15.6)') '  rms-radius ..........',rms
        Write(iw,'(a,15x,2f15.6)') '  charge-radius, r0 ...',rc,r00
        Write(iw,'(a,6f15.6)') '  deformation ......bet',bet
        Write(iw,'(a,6f15.6)') '  q2 deformation,  beta',bet2
        Write(iw,'(a,6f15.6)') '  q4 deformation,  heta',het4
        Write(iw,'(a,6f15.6)') '  deformation .........',def
        Write(iw,'(a,6f15.6)') '  quadrupole moment[b] ',q2
        Write(iw,'(a,6f15.6)') '  hexadecapole moment .',q4
        Write(iw,'(a,6f15.6)') '  <q20> ...............',  &
             (q2(i)*ffdef3,i=1,3)
        Write(iw,'(a,6f15.6)') '  <q40> ...............',  &
             (q4(i)*ffdef5,i=1,3)
        Write(iw,*)
        Write(iw,'(a,6f15.6)')     '  kinetic energy ......',rekin
        Write(iw,'(a,6f15.6)')     '  cmc-diagonal part ...',rekin/hb0*hbzero-rekin
        Write(iw,'(a,6f15.6)')     '  cmc-PAV .............',ECMPAV
        Write(iw,*)
        Write(iw,'(a,30x,6f15.6)') '  volume energy .......',revolpj
        Write(iw,'(a,30x,6f15.6)') '  surface energy ......',resurpj
        Write(iw,'(a,30x,6f15.6)') '  spin-orbit energy ...',respopj
        Write(iw,'(a,30x,6f15.6)') '  coulomb direct ......',recdipj
        Write(iw,'(a,30x,6f15.6)') '  coulomb exchange ....',recexpj
        Write(iw,'(a,30x,6f15.6)') '  tensor energy .......',retenspj
        Write(iw,*)
        Write(iw,'(a,30x,f15.6)')  '  Energy: ehfb(qp) ....',ehfb
        Write(iw,'(a,30x,f15.6)')  '  Energy: ehfb(can,pj).',rehfbcan
        Write(iw,'(a,30x,f15.6)')  '  ehfb(qp)-ehfb(can,pj)',ehfb-rehfbcan
        Write(iw,'(a,30x,f15.6)')  '  Epj-ehfb(can,pj) ....',depnp
        Write(iw,'(a,30x,6f15.6)') '  Energy: Epj=E(PAV) ..',retotpj
        Write(iw,*)
     End Do
     !
     ! APPEND the results to file 'thodef.dat'
     ! ieres=ieresu+ieresl+ieresj+ierebl from module definitions
     If(iappend.Ne.0) Then
        ierest=0
        ! charge buffers
        Do i=1,ieresj       !charge projected buffer
           ierest=ierest+1
           eres(ierest)=eresj(i)
        End Do
        Do i=1,ieresu       !charge unprojected buffer
           ierest=ierest+1
           eres(ierest)=eresu(i)
        End Do
        Do i=1,ieresl       !charge LN  buffer
           ierest=ierest+1
           eres(ierest)=eresl(i)
        End Do
        Do i=1,ieresbl      !charge Blocking buffer
           ierest=ierest+1
           eres(ierest)=eresbl(i)
        End Do
        If(ierest.Ne.ieres) Then
           ierror_flag=ierror_flag+1
           ierror_info(ierror_flag)='STOP: In expectpj: ierest wrong'
           Return  
        Endif
        !write(*,*) nucname,hlabels(1)
        !write(*,*) ereslbl,hlabels(2:3)
        !do i=1,ierest
        !  write(*,*) i,hlabels(i+3),eres(i)
        !enddo
        If(Print_Screen) Then
           ! recording results
100        Continue                        ! complications are due to eagle_ornl
           If(iLST1.Le.0) Then
              Open (unit=iw2,file='hodef.dat',err=100,iostat=i,position='append')
              Write(iw2,'(3(1x,a,1x),160(1x,f14.6))') nucname,ereslbl,eres(1:ierest)
              Close(iw2)
              !do i=1,ierest
              !   write(*,*) i,eres(i)
              !enddo
           Else 
              If (iasswrong(3).Eq.0) Then
                 Open (unit=iw1,file='thodef.dat',err=100,iostat=i,position='append')
                 Write(iw1,'(3(1x,a,1x),160(1x,f14.6))') nucname,ereslbl,eres(1:ierest)
                 Close(iw1)
              End If
           End If
        Endif
     End If
  End If
End Subroutine expectpj
!=======================================================================
!
!===============================================================================================
Subroutine Constraint_or_not(inin_INI0,inin0,icstr0)
  Use HFBTHO
  Implicit None
  Integer(ipr), Intent(in) :: inin_INI0
  Integer(ipr), Intent(inout) :: inin0,icstr0
  If(Abs(inin_INI0).Ge.100) Then
     icstr0=1; inin0=inin_INI0/100
  Else
     icstr0=0; inin0=inin_INI0
  End If
End Subroutine Constraint_or_not
!===============================================================================================
!
!===============================================================================================
Subroutine requested_blocked_level(ib,it)
  !------------------------------------------------------
  ! Search for the requested state to block
  !------------------------------------------------------
  Use HFBTHO ! hfb(ndx2,ndx2)
  Implicit None
  Integer(ipr), Intent(in)  :: ib,it
  Integer(ipr) :: nd,im,k,ndk,na2,nad2,iqn,k0,LAPLUS,OMEGA,n1,n2,n3
  Real(pr)     :: s1,s2,UUk,VVk 
  k0=0
  If(nkblo(it,2).Eq.0)      Return
  If(Parity) Then
     LAPLUS=(ib+1)/2 !Yesp
  Else
     LAPLUS=ib       !Nop
  Endif
  OMEGA=2*LAPLUS-1
  If(nkblo(it,1).Ne.OMEGA) Return
  nd=ID(ib); im=ia(ib); 
  Do k=1,nd
     ndk=k+nd; s1=zero
     Do na2=1,nd
	nad2=na2+nd
	UUk=hfb(na2,ndk); VVk=hfb(nad2,ndk)  
	s2=Max(s1,Abs(UUk),Abs(VVk))
	If(s2.Gt.s1) Then
           s1=s2
           iqn=na2+im  ! the position in [123] numbering
	End If
     Enddo
     ! quantum numbers: Omega,P[n1,n2,n3]=>OMEGA,tpar(npar(iqn))[nz(iqn)+2*nr(iqn)+nl(iqn),nz(iqn),nl(iqn)]
     If(nkblo(it,2).Ne.tpar(npar(iqn))) Cycle
     n3=nl(iqn);          If(nkblo(it,5).Ne.n3) Cycle
     n2=nz(iqn);          If(nkblo(it,4).Ne.n2) Cycle
     n1=n2+2*nr(iqn)+n3;	If(nkblo(it,3).Ne.n1) Cycle
     k0=iqn
     keyblo(it)=1
     bloblo(keyblo(it),it)=ib
     blo123(keyblo(it),it)=k
     Exit
  Enddo
End Subroutine requested_blocked_level
!===============================================================================================
!
!===============================================================================================
Subroutine get_CPU_time (subname,is)
  Use HFBTHO
  Implicit None
  Integer(ipr), Intent(in)   :: is
  Character*(*), Intent(in)  :: subname
  Real(pr), Save :: time1,time2
  !
  If(Print_Screen.And.IDEBUG.Gt.0) Then 
     If(is.Eq.0) Then
	   Call Cpu_time (time1)
     Else
       Call Cpu_time (time2)
       Write(*,'(a,a,a,G16.6)') '  Time in seconds -> ',Trim(subname),':',time2-time1
	 End If  
  End If
End Subroutine get_CPU_time 
!===================================================================================================================================
!#END HFBTHO_SOLVER
!===================================================================================================================================


