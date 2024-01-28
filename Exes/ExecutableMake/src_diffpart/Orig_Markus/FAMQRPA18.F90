#ifndef hide_qrpa
!===================================================================================================================================
!#START qrpa_HFBTHO MODULE
!===================================================================================================================================
Module qrpa_HFBTHO 
  Use HFBTHO 
  Implicit None 
  Character(6) :: qrpa_Version='18' 
  ! Version History 
  !=================================================================================================================================
  ! ver#19: Will add option to do K=0 calculations without simplex-y symmetry.
  !         This allows M1 calculations for K=0 case. (INCOMPLETE, UNDER CONSTRUCTION)
  ! ver#18: Enable equal filling approximation (UNDER DEBUG)
  ! ver#17: Spin-flip and M1 operators for K=1 case. Overhaul of the input file, 
  !         eta parameter removed form qrpa.inp. Some clean-up of the code. 
  !         Added various operators and functionaly coded by Kristian Petrik
  ! ver#16: Included J^2 and s.T terms
  ! ver#15: Bug fix: Icorrect term in pairing rearrangement term in fields.
  ! ver#14: A complete rewrite, almost from the scratch. Expanded to K.ne.0 cases.
  !         Densities are now computed from induced density matrix. Eta parameter is
  !         now ignored (but kept in the input file read-operation for compability).
  !         Basis states are now eigenstates of simplex-y symmetry.
  !
  !         In spherical nuclei for L=1,2,3 different K modes agree perfectly with Coulomb 
  !         included. Tests against the old K=0 FAM code agrees perfectly.
  !         Note: The use of pairing window for densities can be set by variable USE_HFB_PWI
  !
  ! ver#12: Final version with all time-odd components added
  !  Nsh=20 (HFBTHO):  
  !  Time in seconds -> hfbdiag(n): 0.626624    
  !  Time in seconds -> hfbdiag(p): 0.625412    
  !  Time in seconds -> densit: 0.583579    
  !  Time in seconds -> coulom: 0.227499E-02
  !  Time in seconds -> expect: 0.732183E-03
  !  Time in seconds -> field:  0.972986E-03
  !  Time in seconds -> gamdel: 0.960518    
  !  Time in seconds -> broyden: 0.470479E-01
  !  Time in seconds -> per HFBTHO iteration: 2.848  
  !  Nsh=20 (QRPA):  
  !  Time in seconds -> qrpa_DENSIT: 1.09663    
  !  Time in seconds -> qrpa_coulom: 0.287080E-02
  !  Time in seconds -> qrpa_field: 0.613594E-02
  !  Time in seconds -> qrpa_gamdel: 4.42919    
  !  Time in seconds -> per qrpa iteration: 6.4924  
  ! ver#11: Final version with time-odd components j^2 added
  !  Nsh=20:  
  !  Time in seconds -> qrpa_DENSIT: 1.11193           x 2
  !  Time in seconds -> qrpa_coulom: 0.288701E-02      x 2
  !  Time in seconds -> qrpa_field:  0.348902E-02      x 2
  !  Time in seconds -> qrpa_field_Todd: 0.156999E-02  x 2
  !  Time in seconds -> qrpa_gamdel: 6.57947           x 2
  !  Time in seconds -> per qrpa iteration: 16.3223  
  ! ver#10: (MK) time-odd component j^2 added
  ! ver#9:  error fixed when taking derivatives fr0m VA*VB 
  ! ver#8:  Final version with time-even terms only 
  ! ver#7:  Cleaning, strenght function ready, sum rules checked, the module looks ready. 
  !         Started [12/16/2010] first deformed qrpa Zr100 done on [01/09/2011]
  ! ver#6:  Broyden solver added. It solves the X,Y system of about 5000 nonlinear algebraic equations 
  !         for up to 10-40 iterations. Without Broyden, direct iterations diverge for small 
  !         Imag(omega) and even more often. 
  ! ver#5:  Now the main program requests qrpa at the end. The qrpa module decides whether to do it or not.
  !         Do_QRPA starts with hfbdiag to charge the solution U,V for further qrpa calculations.
  !         Option to run qrpa within the pairing window only added. The strenght function still not done. 
  ! ver#4:  Added calculating the external field q.p. matrix elements. Up to date complete
  !         solution of the qrpa problem. Only strenght function remains to be coded. 
  ! ver#2:  qrpa_GAMDEL calculates d,h(omega) and Bd,Bh(omega), substracting HFB h and d 
  !         calculated from hfbdiag. 
  !         Tested orthonormality of U,V block by block (accuracy 10^{-15}).
  !         Allocated X,Y amplitudes and pairing window enforced also to them.
  !         The changes give one-to-one the results of ptho123.f90 (accuracy 10^{-13})
  !         ###Changes in ptho123qrpa### (affect only if qrpa requested)
  !         RESU allocates RVqp,RUqp,REqp,Kqp arrays and HFBDIAG 
  !         populates them with converged HFB solution (the whole space!). Cleaned bug
  !         If(norm_to_improve) Cycle in HFBDIAG 
  !         Timing per qrpa iteration with Nsh=20 is 13.24303889274597 seconds
  ! ver#1:  All required quantyties defined as complex numbers 
  !         The changes give one-to-one the results of ptho123.f90
  !         SPEED TEST on version: 1 (Nsh=14)    
  !          Time in qrpa_densit: 0.4734477996826172 seconds 
  !          Time in densit:      0.1961390972137451 seconds x 2.41
  !          Time in qrpa_field:  0.0036630630493164 seconds
  !          Time in field:       0.0010271072387695 seconds x 3.57
  !          Time in qrpa_gamdel: 0.6993141174316406 seconds
  !          Time in gamdel:      0.4951391220092773 seconds x 1.41
  ! ver#0:  Starting setup:
  !         fields, densit, gamdel, coulom isolated in a new module qrpa_HFBTHO which calles 
  !         qrpa_calculate_U_parameters instead of calculate_U_parameters with the idea all
  !         these subroutines to carry complex (X,Y) qrpa amplitudes. 
  !         The interface to HFBTHO is the Subroutine Do_QRPA called at the main program.
  !         The changes give one-to-one the results of ptho123.f90
  !         Started [12/16/2010]
  !=================================================================================================================================
  ! time-even fields
  Complex(pr), Public, Allocatable         :: vhbn_c(:),vn_c(:),vrn_c(:),vzn_c(:),vdn_c(:),vsn_c(:),davn_c(:),dbvn_c(:)  
  Complex(pr), Public, Allocatable         :: vhbp_c(:),vp_c(:),vrp_c(:),vzp_c(:),vdp_c(:),vsp_c(:),davp_c(:),dbvp_c(:)    
  Complex(pr), Public, Allocatable, Target :: vJRRn_c(:),vJRFIn_c(:),vJRZn_c(:),vJFIRn_c(:),vJFIFIn_c(:)
  Complex(pr), Public, Allocatable, Target :: vJFIZn_c(:),vJZRn_c(:),vJZFIn_c(:),vJZZn_c(:)
  Complex(pr), Public, Allocatable, Target :: vJRRp_c(:),vJRFIp_c(:),vJRZp_c(:),vJFIRp_c(:),vJFIFIp_c(:)
  Complex(pr), Public, Allocatable, Target :: vJFIZp_c(:),vJZRp_c(:),vJZFIp_c(:),vJZZp_c(:)
  Complex(pr), Public, Allocatable         :: cou_c(:), vc_c(:,:) 
  ! time-even densities
  Complex(pr), Public, Allocatable, Target :: aka_c(:,:),bka_c(:,:),ro_c(:,:),tau_c(:,:),dro_c(:,:),dj_c(:,:) 
  Complex(pr), Public, Allocatable, Target :: JRR_c(:,:),JRFI_c(:,:),JRZ_c(:,:),JFIR_c(:,:),JFIFI_c(:,:)
  Complex(pr), Public, Allocatable, Target :: JFIZ_c(:,:),JZR_c(:,:),JZFI_c(:,:),JZZ_c(:,:)
  Complex(pr), Public, Allocatable, Target :: NABLAR_c(:,:),NABLAZ_c(:,:)
  Complex(pr), Public, Allocatable, Target :: rk_c(:,:),ak_c(:,:)   !! rk_c for LN. Not used
  ! HFB densitieds for density dependent part
  Complex(pr), Public, Allocatable, Target :: ro_c_hfb(:,:),tau_c_hfb(:,:),kap_c_hfb(:,:)
  ! time-odd coupling constants
  Real(pr)                                 :: Css(0:1),CDss(0:1),CDeltass(0:1),Cjcjc(0:1),Csnablajc(0:1),CsT(0:1),CsF(0:1)
  ! time-odd fields
  Complex(pr), Public, Allocatable         :: vSIGRn_c(:),vSIGZn_c(:),vSIGFIn_c(:)     !! s, neutrons
  Complex(pr), Public, Allocatable         :: vSIGRp_c(:),vSIGZp_c(:),vSIGFIp_c(:)     !! s, protons
  Complex(pr), Public, Allocatable         :: vSIDRn_c(:),vSIDZn_c(:),vSIDFIn_c(:)     !! \delta s, neutrons
  Complex(pr), Public, Allocatable         :: vSIDRp_c(:),vSIDZp_c(:),vSIDFIp_c(:)     !! \delta s, protons
  Complex(pr), Public, Allocatable         :: vjRn_c(:),vjZn_c(:),vjFIn_c(:)           !! j, neutrons
  Complex(pr), Public, Allocatable         :: vjRp_c(:),vjZp_c(:),vjFIp_c(:)           !! j, protons  
  Complex(pr), Public, Allocatable         :: vTRn_c(:),vTZn_c(:),vTFIn_c(:)           !! T, neutrons
  Complex(pr), Public, Allocatable         :: vTRp_c(:),vTZp_c(:),vTFIp_c(:)           !! T, protons  
  ! time-odd densities.
  Complex(pr), Public, Allocatable, Target :: sor_c(:,:),soz_c(:,:),sofi_c(:,:)        !! spin
  Complex(pr), Public, Allocatable, Target :: dsor_c(:,:),dsoz_c(:,:),dsofi_c(:,:)     !! laplacian spin
  Complex(pr), Public, Allocatable, Target :: rjr_c(:,:),rjz_c(:,:),rjfi_c(:,:)        !! current (= j)
  Complex(pr), Public, Allocatable, Target :: tsor_c(:,:),tsoz_c(:,:),tsofi_c(:,:)     !! spin kinetic (= T)
  Complex(pr), Public, Allocatable, Target :: curjr_c(:,:),curjz_c(:,:),curjfi_c(:,:)  !! curl current (= curl j)
  Complex(pr), Public, Allocatable, Target :: cursr_c(:,:),cursz_c(:,:),cursfi_c(:,:)  !! curl spin
  ! gamdel matrix elements
  Complex(pr), Public, Allocatable         :: hN_c(:),hP_c(:),dN_c(:),dP_c(:)
  Complex(pr), Public, Allocatable         :: BhN_c(:),BhP_c(:),BdN_c(:),BdP_c(:)
  Complex(pr), Public, Allocatable         :: h11N_c(:),h11P_c(:),Bh11N_c(:),Bh11P_c(:)
  Complex(pr), Public, Allocatable         :: FhNab_c(:),FhPab_c(:),BFhNab_c(:),BFhPab_c(:) !! operator F matrix elements 
  Complex(pr), Public, Allocatable         :: FhN_c(:),FhP_c(:),BFhN_c(:),BFhP_c(:)         !! F20 and F02
  Complex(pr), Public, Allocatable         :: F11N_c(:),F11P_c(:),BF11N_c(:),BF11P_c(:)     !! F11 and F11b
  ! operators to remove spurious states
  Complex(pr), Public, Allocatable         :: JyNab_c(:),JyPab_c(:),BJyNab_c(:),BJyPab_c(:) !! J_y operator
  Complex(pr), Public, Allocatable         :: JyN_c(:),JyP_c(:),BJyN_c(:),BJyP_c(:)         !! J_y^20 and J_y^02
  ! uv amplitudes
  Complex(pr), Public, Allocatable         :: GreenNplus(:),GreenPplus(:),GreenNminus(:),GreenPminus(:)  
  Complex(pr), Public, Allocatable         :: GreenNefa(:),GreenPefa(:)
  Complex(pr), Public, Allocatable, Target :: VqpN(:),UqpN(:),VqpP(:),UqpP(:) 
  ! eta density matrices for DENSIT
  Complex(pr), Public, Allocatable :: rhoetadm(:),kapetadm(:),kapetabdm(:)
  ! x,y amplitudes
  Complex(pr), Public, Allocatable, Target :: XN_c(:),YN_c(:),XP_c(:),YP_c(:)
  ! p amplitudes
  Complex(pr), Public, Allocatable, Target :: P1N_c(:),P2N_c(:),P1P_c(:),P2P_c(:)
  !
  ! f-factors for efa
  Complex(pr), Public, Allocatable         :: efa_fN_factor1(:),efa_fN_factor2(:)
  Complex(pr), Public, Allocatable         :: efa_fP_factor1(:),efa_fP_factor2(:)
  !
  ! temporary arrays
  Complex(pr), Public, Allocatable         :: qrpa_AUX(:),qrpa_BUX(:),qrpa_AUX20(:),qrpa_BUX20(:)
  Complex(pr), Public, Allocatable         :: hfb_UA(:),hfb_VA(:)
  !
  ! indexing of blocks
  integer(ipr), Public                     :: iHFB_NB, iQRPA_NB                      !! total number of HFB  and QRPA blocks
  integer(ipr), Public                     :: iHFB_NQP                               !! total number of HFB quasiparticles
  integer(ipr), Public                     :: iHFB_LEN, iQRPA_LEN                    !! total length of HFB and QRPA data arrays     
  integer(ipr), Public, Allocatable        :: iHFB_ID(:), iQRPA_IDI(:), iQRPA_IDF(:) !! dimensions of HFB and QRPA blocks
  integer(ipr), Public, Allocatable        :: iHFB_LP(:), iQRPA_LPI(:), iQRPA_LPF(:) !! Lambda+ of HFB, and Lambda+_i and Lambda+_f of QRPA blocks
  integer(ipr), Public, Allocatable        :: iHFB_OM(:), iQRPA_OMI(:), iQRPA_OMF(:) !! Omega of HFB, and Omega_i and Omega_f of QRPA blocks. Number is 2*real value. E.g. Omega=-3/2 is stored as -3
  integer(ipr), Public, Allocatable        :: iHFB_PR(:)                             !! Parity of the HFB block
  integer(ipr), Public, Allocatable        :: iHFB_DS(:), iQRPA_DS(:)                !! index indicating where the block data starts for given block
  integer(ipr), Public, Allocatable        :: iQRPA_BLI(:), iQRPA_BLF(:)             !! HFB block indeces of the qrpa block referring to the HFBTHO basis wave functions.
  integer(ipr), Public, Allocatable        :: iQRPA_BHI(:), iQRPA_BHF(:)             !! HFB block indeces of the qrpa block referring to the iHFB_??() indexing
  integer(ipr), Public, Allocatable        :: iQRPA_TRANSPOSE(:)                     !! index of QRPA block which corresponds transpose of the given qrpa block index
  Integer(ipr), Public, Allocatable        :: iHFB_QPBS(:)                           !! index of first q.p. state strating a block. Dimension = iHFB_NB
  Integer(ipr), Public, Allocatable        :: ipwi_hfb(:,:)                          !! states inside pairing window
  !
  ! qrpa_Broyden
  Character(1), Public                     :: qrpa_bbroyden
  Real(pr), Public, Allocatable            :: qrpa_broin(:),qrpa_broout(:)
  ! qrpa common variables
  Complex(pr), Public                      :: snz_qrpa(2),dsnz_qrpa(2),qrpa_omega
  Integer(ipr), Public                     :: gfpoles_imN,gfpoles_imP,iter_qrpa
  Integer(ipr)                             :: icacou_qrpa
  Real(pr), Public                         :: ImStrengthN,ImStrengthP,ReStrengthN,ReStrengthP
  Real(pr), Public                         :: ImTransStrN,ImTransStrP,ReTransStrN,ReTransStrP
  Real(pr), Public                         :: eeffN, eeffP, gleffN, gleffP, gseffN, gseffP  !! effective charges
  ! qrpa variables taken from the qrpa.inp file and deduced from those
  Integer(ipr), Public                     :: qrpa_calc_mode,qrpa_file_io                                !! calculation mode and file io
  Integer(ipr), Public                     :: T_qrpa_responce,L_qrpa_responce,K_qrpa_responce,P_qrpa_response,QRPA_Operator_type !! isospin-type, L, K, parity, and opr type
  Integer(ipr), Public                     :: max_iter_qrpa,max_qrpa_points
  Integer(ipr), Public                     :: qrpa_nbroyden
  Real(pr), Public                         :: qrpa_alphamix
  Real(pr), Public                         :: qrpa_eps
  Complex(pr), Public                      :: qrpa_eta !! unused
  Real(pr), Public                         :: qrpa_Romega,qrpa_Iomega,qrpa_Romega_step,qrpa_Iomega_step  !! line params
  Real(pr), Public                         :: qrpa_circRe,qrpa_circIm,qrpa_circRad                       !! circle params
  Real(pr), Public                         :: qrpa_HcircRe,qrpa_HcircRad                                 !! half-circle params
  !
  ! various setups
  Logical                                  :: USE_HFB_PWI                            !! If true, densities are calculated with HFB pairing window
  Logical                                  :: USE_ONLY_F20 = .false.                 !! If true, X and Y directly from F20,F02,F11. H20, H02, H11 ignored. _use for debug only_
  Logical                                  :: USE_ONLY_F11EFA = .false.               !! If true, only F11 with efa used. _use for debug only_
  Logical                                  :: USE_NORMALHO_K0 = .false.              !! If true, normal HO basis used for K=0 case
  Logical                                  :: ENABLE_N_EFA = .false.                 !! If true, activates neutron EFA. Set true in initialization, if needed
  Logical                                  :: ENABLE_P_EFA = .false.                 !! If true, activates proton EFA. Set true in initialization, if needed
  !
  !=================================================================================================================================
  !
Contains 
  !
  !=================================================================================================================================
  Subroutine qrpa_calculate_timeodd_CC
    ! from skyrme force
    Cjcjc(0)     =  -Ctau(0) ;       Cjcjc(1) =  -Ctau(1)
    Csnablajc(0) =   Crdj(0) ;   Csnablajc(1) =   Crdj(1) 
    Css(0)       = -2.0_pr*Crho(0)/3.0_pr-Crho(1)
    Css(1)       = -Crho(0)/3.0_pr
    CDss(0)      = -2.0_pr*CDrho(0)/3.0_pr-CDrho(1)
    CDss(1)      = -CDrho(0)/3.0_pr
    CDeltass(0)  = (Ctau(0)+3.0_pr*Ctau(1)  -4.0_pr*(Crdr(0)+Crdr(1)))/8.0_pr   
    CDeltass(1)  = (3.0_pr*(Ctau(0)-Ctau(1))-4.0_pr*(Crdr(0)+Crdr(1)))/24.0_pr
    CsT(0)       = -CJ(0)
    CsT(1)       = -CJ(1)
    !
  End Subroutine qrpa_calculate_timeodd_CC
  !=================================================================================================================================
  !
  !   Modified Bessel I_n(x)
  ! Returns I_n(x) * Exp(-x)
  ! Uses bessel I_0(x) and I_1(x) from HFBTHO
  !=================================================================================================================================
  function modbesselI(n,x)
    use Bessik
    implicit none
    integer(ipr), intent(in) :: n
    Real(pr), intent(in) :: x
    real(pr) :: modbesselI 
    !
    Real(pr) :: I0,I1,I2,I3, eps
    !
    eps=Spacing(1.0_pr);
    !
    I0 = besei0(x)
    I1 = besei1(x)
    I2 = I0 - 2.0_pr*I1/(x+eps)
    I3 = I1 - 4.0_pr*I2/(x+eps)
    !
    modbesselI = 0.0_pr
    select case(n)
      case(0)
        modbesselI = I0
      case(1)
        modbesselI = I1
      case(2) 
        modbesselI = I2
      case(3)
        modbesselI = I3
      case default
        STOP "modBesselI: n not in [0,3]!"
   end select
   return
  end function modbesselI
  !=================================================================================================================================
  !
  ! Direct Coulomb part. This is now done in a similar manner as in HFBTHO v200d
  !
  !=================================================================================================================================
  Subroutine qrpa_coulom
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Coulom-field (direct part), Gogny prescription, complex numbers
    !-------------------------------------------------------------------------------------------------------------------------------
    Use HFBTHO
    Use bessik
    Implicit None
    Integer(ipr), Save :: i,j,k
    Real(pr), Save :: zd2,y1,y2,xx1,s1,vik,f,r,r1,fac1,fac2,rr2,z,z1,zd1,t,  &
                      bb,r2,r12,rrr,rz1,rz2,rrz1,rrz2,xx,rk1,rip1,rkp1,alpha,&
                      beta,xxx
    Complex(pr) :: cone, czero
    cone = 1.0_pr;  czero = 0.0_pr
    !
    Call get_CPU_time('qrpa_coulom',0)
    !
    If(icacou_qrpa.Eq.0 .AND. K_QRPA_RESPONCE .EQ. 0) Then

      icacou_qrpa = 1
      vc_c = vc       ! vc_c can be taken from the HFB part, since there too K=0

    else if (icacou_qrpa .EQ. 0 .AND. K_QRPA_RESPONCE .NE. 0) Then

      icacou_qrpa = 1
      !
      ! For parity-breaking shapes, the Coulomb potential was incorrectly
      ! calculated by assuming the two intervals [0,+\infty[ and ]-infty,0]
      ! were equivalent (see also below). This bug was corrected in version
      ! 139a
      If(Parity) Then
        fac1 = (-1)**K_QRPA_RESPONCE;  fac2 = one 
      Else
        fac1 = zero; fac2 = two
      End If
      bb=50.0_pr          ! Length scale L
      beta=2.00_pr
      alpha=one/beta
      f=chargee2/Sqrt(pi) ! e^2/Sqrt(pi)

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
              !!y2=besei0(s1)           ! I0( 2 u^2 r r' ) * exp(-2 u^2 r r')
              y2=modbesselI(K_QRPA_RESPONCE,s1)
              xx1=xx1+fac2*wleg(j)*y2*(Exp(-rz1*y1) + fac1*Exp(-rz2*y1)) / xxx
           End Do
           vik=f*xx1/bb
           !
           vc_c(i,k)=vik*wdcor(k)  !wdcor=pi*wh*wl*bz*bp*bp
           vc_c(k,i)=vik*wdcor(i)  !wdcor=pi*wh*wl*bz*bp*bp
           !
        End Do  !k
      End Do  !i
    End If  !! else if (icacou_qrpa .EQ. 0 .AND. K_QRPA_RESPONCE .NE. 0)
    ! Calculation of the Coulomb field
    cou_c=zero
    Call zgemm('n','n',nghl,1,nghl,cone,vc_c,nghl,ro_c(:,2),nghl,czero,cou_c,nghl)
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('coulom',1)
  End Subroutine qrpa_coulom
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_DENSIT
    !---------------------------------------------------------------------
    ! local densities in coordinate space. 
    !---------------------------------------------------------------------
    Implicit None
    Real(pr)     :: y
    Complex(pr)  :: ss
    Complex(pr)  :: TEMP1_c,TEMP2_c,TEMP3_c,TEMP4_c,TEMP5_c,TEMP6_c,TEMP7_c,TEMP8_c,TEMP9_c,TEMP1_c0,TEMP1_c1
    Complex(pr)  :: TEMP10_c,TEMP11_c,TEMP12_c,TEMP13_c,TEMP14_c,TEMP15_c,TEMP16_c,TEMP17_c,TEMP18_c
    Complex(pr)  :: TEMP_rr,TEMP_rf,TEMP_rz,TEMP_fr,TEMP_ff,TEMP_fz,TEMP_zr,TEMP_zf,TEMP_zz
    Complex(pr)  :: FIA1FIB2
    !
    Integer(ipr) :: nsa,nla,ihil,n1,n2
    Integer(ipr) :: imena,imenb,im,ib,it,J,JA
    Integer(ipr), Allocatable :: LAM1(:), LAM2(:), SS1(:), SS2(:)
    integer(ipr) :: ovlap
    complex(pr)  :: ovlapr,ovlapphi,ovlapz
    !
    ! pointers
    Complex(pr),  Pointer :: TAKA_c(:),TBKA_c(:),TRO_c(:),TDJ_c(:),TTAU_c(:),TDRO_c(:)
    Complex(pr),  Pointer :: TJRR_c(:),TJRFI_c(:),TJRZ_c(:),TJFIR_c(:),TJFIFI_c(:),TJFIZ_c(:)
    Complex(pr),  Pointer :: TJZR_c(:), TJZFI_c(:),TJZZ_c(:)
    Complex(pr),  Pointer :: TNABLAR_c(:),TNABLAZ_c(:)
    !
    Complex(pr),  Pointer :: Psor_c(:),Psoz_c(:),Psofi_c(:)       !! spin
    Complex(pr),  Pointer :: Pdsor_c(:),Pdsoz_c(:),Pdsofi_c(:)    !! laplacian spin
    Complex(pr),  Pointer :: Prjr_c(:),Prjz_c(:),Prjfi_c(:)       !! current (=j)
    Complex(pr),  Pointer :: Ptsor_c(:),Ptsoz_c(:),Ptsofi_c(:)    !! spin kinetic
    Complex(pr),  Pointer :: Pcurjr_c(:),Pcurjz_c(:),Pcurjfi_c(:) !! curl current (= curl j)
    Complex(pr),  Pointer :: Pcursr_c(:),Pcursz_c(:),Pcursfi_c(:) !! curl s (= curl s)
    Complex(pr),  Pointer :: VqpPoA(:),UqpPoA(:),VqpPoB(:),UqpPoB(:)
    Integer(ipr), Pointer :: Kpwi_qrpa(:)
    Complex(pr),  Pointer :: Xqrpa(:),Yqrpa(:)
    Complex(pr),  Pointer :: P1qrpa(:),P2qrpa(:)
    Complex(pr),  Pointer :: Vqp(:),Uqp(:)
    ! helpers
    Integer(ipr) :: ipa,ipq,ipt,ibHFBi,ibHFBf,ibTHOi,ibTHOf,lda
    Real(pr)     :: SSU,SSD
    ! density matrix elements
    Complex(pr) :: rhme,kpme,kbme
    !
    Complex(pr)               :: iunit
    Complex(pr), Allocatable  :: FI_cA(:),FI_cB(:), FIR_cA(:),FIR_cB(:), FIZ_cA(:),FIZ_cB(:)
    Complex(pr), Allocatable  :: FID2_cA(:),FID2_cB(:)
    !
    Complex(pr), Allocatable  :: FIU_cA(:),FID_cA(:),FIUR_cA(:),FIDR_cA(:)
    Complex(pr), Allocatable  :: FIUD2N_cA(:),FIDD2N_cA(:),FIUZ_cA(:),FIDZ_cA(:)    
    Complex(pr), Allocatable  :: FIU_cB(:),FID_cB(:),FIUR_cB(:),FIDR_cB(:)
    Complex(pr), Allocatable  :: FIUD2N_cB(:),FIDD2N_cB(:),FIUZ_cB(:),FIDZ_cB(:)    
    !
    Parameter(iunit = Cmplx(0.0_pr,1.0_pr,kind=pr))    
    !
    Call get_CPU_time('qrpa_DENSIT',0)
    !
    Do it=1,2 !! flush also densities in the case of zero paricle number
       ! pointers
       If(it.Eq.1) Then
         Vqp=>VqpN; Uqp=>UqpN; Xqrpa=>XN_c; Yqrpa=>YN_c; P1qrpa=>P1N_c(:); P2qrpa=>P2N_c(:);
       Else
         Vqp=>VqpP; Uqp=>UqpP; Xqrpa=>XP_c; Yqrpa=>YP_c; P1qrpa=>P1P_c(:); P2qrpa=>P2P_c(:);
       Endif
       !
       !------------------------------------------------
       ! Calculation of density matrices.
       ! The use of HFB pairing window here is set to be 
       ! equivalent to the old K=0 FAM module. Therefore,
       ! in matrix product, only the other HFB matrix is
       ! modified
       !------------------------------------------------
       if(USE_HFB_PWI) then
         ! set aux arrays
         If(it.Eq.1) Then
           hfb_UA = UqpN ; hfb_VA = VqpN
         Else
           hfb_UA = UqpP ; hfb_VA = VqpP
         Endif
         ! remove states outside pairing window
         do ib=1,iHFB_NB
           im=0
           do ja=1,iHFB_ID(ib)
             do j=1,iHFB_ID(ib) 
               if(ja.gt.ipwi_hfb(ib,it)) then
                hfb_UA( iHFB_DS(ib)+im ) = zero;
                hfb_VA( iHFB_DS(ib)+im ) = zero;
               end if
               im=im+1
             end do !! components
           end do !! q.p. states
         end do !! blocks
         !
         ! density matrices with pairing window
         call qrpa_TriProd('n','n','c',Uqp,Xqrpa,hfb_Va,qrpa_aux)
         call qrpa_TriProd('n','t','c',hfb_Va,Yqrpa,Uqp,rhoetadm)
         rhoetadm = -rhoetadm -qrpa_aux
         !
         call qrpa_TriProd('n','t','c',Vqp,Yqrpa,hfb_Va,kapetadm)
         call qrpa_TriProd('n','n','c',hfb_Ua,Xqrpa,Uqp,qrpa_aux)
         kapetadm = -kapetadm +qrpa_aux
         call qrpa_TriProd('n','c','c',Vqp,Xqrpa,hfb_Va,kapetabdm)
         call qrpa_TriProd('n','*','c',hfb_Ua,Yqrpa,Uqp,qrpa_aux)
         kapetabdm = -kapetabdm +qrpa_aux
         !
         ! add efa contibutions
         if((it.eq.1 .and. ENABLE_N_EFA) .or. ((it.eq.2 .and. ENABLE_P_EFA)) ) then

           !! ivestigate here proper use of pairing window! Now pairing window is set for both matrices.
           !! with efa, this should not have impact since blocked states are always withing pairing window (?)
           !! with finite temperature, this may have impact...
           call qrpa_TriProd('n','n','c',hfb_Ua,P1qrpa,hfb_Ua,qrpa_aux)
           call qrpa_TriProd('n','t','c',hfb_Va,P2qrpa,hfb_Va,qrpa_bux)
           rhoetadm = rhoetadm +qrpa_aux -qrpa_bux

           call qrpa_TriProd('n','t','c',hfb_Va,P2qrpa,hfb_Ua,qrpa_aux)
           call qrpa_TriProd('n','n','c',hfb_Ua,P1qrpa,hfb_Va,qrpa_bux)
           kapetadm = kapetadm +qrpa_aux +qrpa_bux

           call qrpa_TriProd('n','*','c',hfb_Va,P2qrpa,hfb_Ua,qrpa_aux)
           call qrpa_TriProd('n','c','c',hfb_Ua,P1qrpa,hfb_Va,qrpa_bux)
           kapetabdm = kapetabdm +qrpa_aux +qrpa_bux

         end if

       else
         !
         ! density matrices without pairing window
         call qrpa_TriProd('n','n','c',Uqp,Xqrpa,Vqp,qrpa_aux)
         call qrpa_TriProd('n','t','c',Vqp,Yqrpa,Uqp,rhoetadm)
         rhoetadm = -rhoetadm -qrpa_aux
         !
         call qrpa_TriProd('n','t','c',Vqp,Yqrpa,Vqp,kapetadm)
         call qrpa_TriProd('n','n','c',Uqp,Xqrpa,Uqp,qrpa_aux)
         kapetadm = -kapetadm +qrpa_aux
         call qrpa_TriProd('n','c','c',Vqp,Xqrpa,Vqp,kapetabdm)
         call qrpa_TriProd('n','*','c',Uqp,Yqrpa,Uqp,qrpa_aux)
         kapetabdm = -kapetabdm +qrpa_aux

         if((it.eq.1 .and. ENABLE_N_EFA) .or. ((it.eq.2 .and. ENABLE_P_EFA)) ) then
           call qrpa_TriProd('n','n','c',Uqp,P1qrpa,Uqp,qrpa_aux)
           call qrpa_TriProd('n','t','c',Vqp,P2qrpa,Vqp,qrpa_bux)
           rhoetadm = rhoetadm +qrpa_aux -qrpa_bux

           call qrpa_TriProd('n','t','c',Vqp,P2qrpa,Uqp,qrpa_aux)
           call qrpa_TriProd('n','n','c',Uqp,P1qrpa,Vqp,qrpa_bux)
           kapetadm = kapetadm +qrpa_aux +qrpa_bux

           call qrpa_TriProd('n','*','c',Vqp,P2qrpa,Uqp,qrpa_aux)
           call qrpa_TriProd('n','c','c',Uqp,P1qrpa,Vqp,qrpa_bux)
           kapetabdm = kapetabdm +qrpa_aux +qrpa_bux
         end if

       end if !! USE_HFB_PWI
       !
       ! pointers to densities, time even
       TRO_c=>ro_c(:,it);         TTAU_c=>tau_c(:,it);       TDJ_c=>dj_c(:,it);     TDRO_c=>dro_c(:,it);  
       TJRR_c=>JRR_c(:,it);       TJRFI_c=>JRFI_c(:,it);     TJRZ_c=>JRZ_c(:,it);   TJFIR_c=>JFIR_c(:,it);
       TJFIFI_c=>JFIFI_c(:,it);   TJFIZ_c=>JFIZ_c(:,it);     TJZR_c=>JZR_c(:,it);   TJZFI_c=>JZFI_c(:,it);
       TJZZ_c=>JZZ_c(:,it);
       TNABLAR_c=>NABLAR_c(:,it); TNABLAZ_c=>NABLAZ_c(:,it); TAKA_c=>aka_c(:,it);   TBKA_c=>bka_c(:,it);     
       !
       ! pointers to densities, time odd
       Psor_c=>sor_c(:,it);       Psoz_c=>soz_c(:,it);       Psofi_c=>sofi_c(:,it);
       Pdsor_c=>dsor_c(:,it);     Pdsoz_c=>dsoz_c(:,it);     Pdsofi_c=>dsofi_c(:,it);
       Prjr_c=>rjr_c(:,it);       Prjz_c=>rjz_c(:,it);       Prjfi_c=>rjfi_c(:,it);
       Ptsor_c=>tsor_c(:,it);     Ptsoz_c=>tsoz_c(:,it);     Ptsofi_c=>tsofi_c(:,it);
       Pcurjr_c=>curjr_c(:,it);   Pcurjz_c=>curjz_c(:,it);   Pcurjfi_c=>curjfi_c(:,it);
       Pcursr_c=>cursr_c(:,it);   Pcursz_c=>cursz_c(:,it);   Pcursfi_c=>cursfi_c(:,it);	   
       ! ZERO THE DENSITIES
       Do ihil=1,nghl
          TRO_c(ihil)=ZERO;     TTAU_c(ihil)=ZERO;    TDJ_c(ihil)=ZERO;   TDRO_c(ihil)=ZERO;
          TJRR_c(ihil)=ZERO;    TJRFI_c(ihil)=ZERO;   TJRZ_c(ihil)=ZERO;  TJFIR_c(ihil)=ZERO;
          TJFIFI_c(ihil)=ZERO;  TJFIZ_c(ihil)=ZERO;   TJZR_c(ihil)=ZERO;  TJZFI_c(ihil)=ZERO;
          TJZZ_c(ihil)=ZERO;
          TNABLAR_c(ihil)=ZERO; TNABLAZ_c(ihil)=ZERO; TAKA_c(ihil)=ZERO;  TBKA_c(ihil)=ZERO; 
          Psor_c(ihil)=ZERO;    Psoz_c(ihil)=ZERO;    Psofi_c(ihil)=ZERO; 
          Pdsor_c(ihil)=ZERO;   Pdsoz_c(ihil)=ZERO;   Pdsofi_c(ihil)=ZERO; 
          Prjr_c(ihil)=ZERO;    Prjz_c(ihil)=ZERO;    Prjfi_c(ihil)=ZERO; 
          Ptsor_c(ihil)=ZERO;   Ptsoz_c(ihil)=ZERO;   Ptsofi_c(ihil)=ZERO;
          Pcurjr_c(ihil)=ZERO;  Pcurjz_c(ihil)=ZERO;  Pcurjfi_c(ihil)=ZERO;
          Pcursr_c(ihil)=ZERO;  Pcursz_c(ihil)=ZERO;  Pcursfi_c(ihil)=ZERO;
       End Do
       ! in the case of zero particle number, only flush densities
       If((npr_INI(1).Eq.0).And.(it.Eq.1)) Cycle
       If((npr_INI(2).Eq.0).And.(it.Eq.2)) Cycle
       !-----------------------------------------------
       ! SCAN OVER GAUSS INTEGRATION POINTS
       !-----------------------------------------------
       Do ihil=1,nghl
         !-----------------------------------------------
         ! SCAN OVER QRPA BLOCKS 
         !-----------------------------------------------
         Do ib=1,iQRPA_NB
           If(kindhfb.Lt.0) Then
             STOP "LN does not work yet!"
           End if
           ibHFBi = iQRPA_BHI(ib) ; ibHFBf = iQRPA_BHF(ib)
           ibTHOi = iQRPA_BLI(ib) ; ibTHOf = iQRPA_BLF(ib)
           ipa = iHFB_DS(ibHFBi) ; ipq = iQRPA_DS(ib) ; ipt = iQRPA_DS(iQRPA_TRANSPOSE(ib))
           !
           ! allocate arrays
           IMENA=iQRPA_IDF(ib) ; IMENB=iQRPA_IDI(ib)
           If(Allocated(FI_cA)) Deallocate(FI_cA,FIR_cA,FIZ_cA,FID2_cA)
           Allocate(FI_cA(IMENA),FIR_cA(IMENA),FIZ_cA(IMENA),FID2_cA(IMENA))
           If(Allocated(FI_cB)) Deallocate(FI_cB,FIR_cB,FIZ_cB,FID2_cB)
           Allocate(FI_cB(IMENB),FIR_cB(IMENB),FIZ_cB(IMENB),FID2_cB(IMENB))
           if(allocated(lam1)) Deallocate(lam1,ss1,lam2,ss2)
           allocate(LAM1(IMENA),SS1(IMENA),LAM2(IMENB),SS2(IMENB))
           !
           ! run trough the left basis vectors and fill the array
           Do N1=1,iQRPA_IDF(ib)
             IM=ia(iQRPA_BLF(ib)); JA=IM+N1; NSA=NS(JA); SSU=Max(NSA,0); SSD=Max(-NSA,0)
             FI_cA(N1)   = QHLA_opt(JA,ihil)
             FIR_cA(N1)  = FI1R_opt(JA,ihil)
             FIZ_cA(N1)  = FI1Z_opt(JA,ihil)
             FID2_cA(N1) = FI2D_opt(JA,ihil)
             SS1(N1) = NSA
             if(NSA.gt.0) then
               LAM1(N1) = (iQRPA_OMF(ib) -1)/2
             else
               LAM1(N1) = -(iQRPA_OMF(ib) +1)/2
             end if
             If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) LAM1(N1) = abs(LAM1(N1))    !! positive omega states have positive lambda
           End Do
           !
           ! run trough the right basis vectors and fill the array
           Do N1=1,iQRPA_IDI(ib)
             IM=ia(iQRPA_BLI(ib)); JA=IM+N1; NSA=NS(JA); SSU=Max(NSA,0); SSD=Max(-NSA,0)
             FI_cB(N1)   = QHLA_opt(JA,ihil)
             FIR_cB(N1)  = FI1R_opt(JA,ihil)
             FIZ_cB(N1)  = FI1Z_opt(JA,ihil)
             FID2_cB(N1) = FI2D_opt(JA,ihil)
             SS2(n1) = NSA
             if(NSA.gt.0) then
               LAM2(N1) = (iQRPA_OMI(ib) -1)/2
             else
               LAM2(N1) = -(iQRPA_OMI(ib) +1)/2
             end if
             If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) LAM2(N1) = abs(LAM2(N1))
           End Do

           y=y_opt(ihil);
           j=iQRPA_DS(ib)-1
           !
           Do N2=1,iQRPA_IDI(ib)
             Do N1=1,iQRPA_IDF(ib)
                !-----------------------------------------------
                ! go one step forward in density matrix arrays
                !-----------------------------------------------
                j=j+1  
                !-----------------------------------------------
                ! induced density matrix and kappa matrices
                !-----------------------------------------------
                rhme = rhoetadm(j) ; kpme = kapetadm(j) ; kbme = kapetabdm(j)  
                !-----------------------------------------------
                ! TIME-EVEN PART
                !-----------------------------------------------
                !
                ovlap = 0 ; 
                if( abs(LAM1(n1) + LAM2(n2)) .eq. abs(K_qrpa_responce) ) ovlap = 1
                if( abs(LAM1(n1) - LAM2(n2)) .eq. abs(K_qrpa_responce) ) ovlap = 1
                !
                ovlapphi = 0; ovlapr = 0; 
                ovlapz = -ovlap;
                ovlapr = (1-ovlap)*(-IUnit)
                ovlapphi = (1-ovlap)*(-1)
                FIA1FIB2 = FI_cA(N1)*FI_cB(N2)  !! \Phi_a(r,z) \Phi_b(r,z). Used to calculate many densities
                !
                ! define overlaps for normal HO basis. 
                If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) then
                   if(LAM1(n1) .eq. LAM2(n2)) then ; ovlap = 1 ; else ; ovlap = 0 ; end if
                   ovlapz = ovlap
                   ovlapr = (1-ovlap)
                   ovlapphi = (1-ovlap)*Iunit
                   if( SS1(n1).lt.0 ) ovlapz = ovlapz * (-1.0_pr)
                   if( SS2(n2).lt.0 ) ovlapphi = ovlapphi * (-1.0_pr)
                end if
                !
                !  ****************
                !   time-even part
                !  ****************
                !
                TEMP1_c = FIA1FIB2*ovlap
                !
                TEMP2_c = ( FIR_cA(N1)*FIR_cB(N2) +FIZ_cA(N1)*FIZ_cB(N2) &
                           +y*y*FIA1FIB2*LAM1(N1)*LAM2(N2) )*ovlap
                !
                TEMP3_c = ( FID2_cA(N1)*FI_cB(N2) +FI_cA(N1)*FID2_cB(N2) +two*FIR_cA(N1)*FIR_cB(N2) &
                           +two*FIZ_cA(N1)*FIZ_cB(N2) -y*y*FIA1FIB2*K_qrpa_responce**2 )*ovlap 
                !
                TEMP4_c = (-( IUnit*LAM2(N2)*y*FIR_cA(N1)*FI_cB(N2)  + IUnit*LAM1(N1)*y*FI_cA(N1)*FIR_cB(N2) )*ovlapz & 
                           +( -IUnit*LAM1(N1)*y*FI_cA(N1)*FIZ_cB(N2) + IUnit*LAM2(N2)*y*FIZ_cA(N1)*FI_cB(N2) )*ovlapr &
                           +( FIZ_cA(N1)*FIR_cB(N2)-FIR_cA(N1)*FIZ_cB(N2) )*ovlapphi )/IUnit
                !
                ! K=0 and !=0 tensor part
                TEMP_rf = -iunit*( FI_cA(N1)*FIR_cB(N2) -FIR_cA(N1)*FI_cB(N2) )*half*ovlapphi
                TEMP_fr = -iunit*( -iunit*(LAM1(N1)-LAM2(N2))*y*FIA1FIB2 )*half*ovlapr
                TEMP_fz = -iunit*( -iunit*(LAM1(N1)+LAM2(N2))*y*FIA1FIB2 )*half*ovlapz
                TEMP_zf = -iunit*( FI_cA(N1)*FIZ_cB(N2) -FIZ_cA(N1)*FI_cB(N2) )*half*ovlapphi
                !
                ! K.ne.0 tensor part
                TEMP_rr = -iunit*( FI_cA(N1)*FIR_cB(N2) -FIR_cA(N1)*FI_cB(N2) )*half*ovlapr
                TEMP_rz = -iunit*( FI_cA(N1)*FIR_cB(N2) -FIR_cA(N1)*FI_cB(N2) )*half*ovlapz
                TEMP_ff = -iunit*( -iunit*(LAM1(N1)-LAM2(N2))*y*FIA1FIB2 )*half*ovlapphi
                TEMP_zr = -iunit*( FI_cA(N1)*FIZ_cB(N2) -FIZ_cA(N1)*FI_cB(N2) )*half*ovlapr
                TEMP_zz = -iunit*( FI_cA(N1)*FIZ_cB(N2) -FIZ_cA(N1)*FI_cB(N2) )*half*ovlapz
                !
                ! these densities have (e^iKphi - e^-iKphi) structure
                ! here e^+iKphi part defines the density and sign
                if(LAM1(n1) - LAM2(n2) .gt. 0) then
                  TEMP_rz = -TEMP_rz; TEMP_zz = -TEMP_zz;
                end if
                if(LAM1(n1) + LAM2(n2) +1 .eq. -abs(K_qrpa_responce)) then
                  TEMP_rr = -TEMP_rr; TEMP_ff = -TEMP_ff; TEMP_zr = -TEMP_zr;
                end if
                ! put zero on those which don't contribute when K=0
                if(K_qrpa_responce.eq.0) then
                  TEMP_rr = zero; TEMP_rz = zero; TEMP_ff = zero; TEMP_zr = zero; TEMP_zz = zero; 
                end if
                !
                ! pairing
                TAKA_c(IHIL)=TAKA_c(IHIL)+TEMP1_c*kpme
                TBKA_c(ihil)=TBKA_c(ihil)+TEMP1_c*kbme
                !
                ! correct densities when using normal HO basis
                If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) then
                  TEMP4_c = (+( IUnit*LAM2(N2)*y*FIR_cA(N1)*FI_cB(N2) + IUnit*LAM1(N1)*y*FI_cA(N1)*FIR_cB(N2) )*ovlapz & 
                             -( IUnit*LAM1(N1)*y*FI_cA(N1)*FIZ_cB(N2) + IUnit*LAM2(N2)*y*FIZ_cA(N1)*FI_cB(N2) )*ovlapr &
                             +( FIZ_cA(N1)*FIR_cB(N2)-FIR_cA(N1)*FIZ_cB(N2) )*ovlapphi )/IUnit
                  !! correct here tensor terms too
                  if(abs(CJ(0)) .ge. spacing(1.0_pr) .or. abs(CJ(1)) .ge. spacing(1.0_pr)) stop "tensor term support missing for this case"
                end if
                !
                ! time-even matter densities
                TRO_c(IHIL)   = TRO_c(IHIL)   + TEMP1_c*rhme
                TTAU_c(IHIL)  = TTAU_c(IHIL)  + TEMP2_c*rhme
                TDRO_c(IHIL)  = TDRO_c(IHIL)  + TEMP3_c*rhme
                TDJ_c(IHIL)   = TDJ_c(IHIL)   + TEMP4_c*rhme
                TJRFI_c(ihil) = TJRFI_c(ihil) + TEMP_rf*rhme
                TJFIR_c(ihil) = TJFIR_c(ihil) + TEMP_fr*rhme
                TJFIZ_c(ihil) = TJFIZ_c(ihil) + TEMP_fz*rhme
                TJZFI_c(ihil) = TJZFI_c(ihil) + TEMP_zf*rhme
                !
                TJRR_c(ihil)  = TJRR_c(ihil)  + TEMP_rr*rhme
                TJRZ_c(ihil)  = TJRZ_c(ihil)  + TEMP_rz*rhme
                TJFIFI_c(ihil)= TJFIFI_c(ihil)+ TEMP_ff*rhme
                TJZR_c(ihil)  = TJZR_c(ihil)  + TEMP_zr*rhme
                TJZZ_c(ihil)  = TJZZ_c(ihil)  + TEMP_zz*rhme
                !  ***************
                !   time-odd part
                !  ***************
                !
                TEMP3_c = 0.0_pr; TEMP4_c = 0.0_pr; TEMP6_c = 0.0_pr; TEMP7_c = 0.0_pr; TEMP9_c = 0.0_pr; TEMP10_c= 0.0_pr; 
                TEMP11_c= 0.0_pr; TEMP12_c= 0.0_pr; TEMP13_c= 0.0_pr; TEMP14_c= 0.0_pr; TEMP15_c= 0.0_pr; TEMP16_c= 0.0_pr;
                TEMP17_c= 0.0_pr; TEMP18_c= 0.0_pr;  
                !
                TEMP1_c = -iunit*( FIR_cA(N1)*FI_cB(N2) - FI_cA(N1)*FIR_cB(N2) )*half*ovlap; !! j_r
                TEMP2_c = -iunit*( FIZ_cA(N1)*FI_cB(N2) - FI_cA(N1)*FIZ_cB(N2) )*half*ovlap; !! j_z
                !
                TEMP5_c = FIA1FIB2*ovlapphi                                                  !! s_phi
                !
                TEMP8_c = (FID2_cA(N1)*FI_cB(N2)+FI_cA(N1)*FID2_cB(N2)+two*FIR_cA(N1)*FIR_cB(N2)+two*FIZ_cA(N1)*FIZ_cB(N2) &
                          -y*y*FIA1FIB2*K_qrpa_responce**2 -y*y*FIA1FIB2 )*ovlapphi                  !! (Delta s)_phi
                !
                TEMP11_c = -iunit*( FIR_cA(N1)*FIZ_cB(N2) - FIZ_cA(N1)*FIR_cB(N2) )*two*half*ovlap   !! (curl j)_phi
                !
                TEMP13_c = -(FIZ_cA(N1)*FI_cB(N2) +FI_cA(N1)*FIZ_cB(N2))*ovlapphi                    !! (curl s)_r
                TEMP15_c = (y*FIA1FIB2 +FIR_cA(N1)*FI_cB(N2) +FI_cA(N1)*FIR_cB(N2))*ovlapphi         !! (curl s)_z
                !
                TEMP17_c = ( FIR_cA(N1)*FIR_cB(N2) +FIZ_cA(N1)*FIZ_cB(N2) &
                            -y*y*FIA1FIB2*LAM1(N1)*LAM2(N2) )*ovlapphi                               !! T_fi
                !
                ! these densities have (e^iKphi - e^-iKphi) structure
                ! here e^+iKphi part defines the density and sign
                if(LAM1(n1) + LAM2(n2) +1 .eq. abs(K_qrpa_responce)) then                
                  TEMP4_c = FIA1FIB2*ovlapr                                                 !! s_r
                  !
                  TEMP7_c = (FID2_cA(N1)*FI_cB(N2)+FI_cA(N1)*FID2_cB(N2)+two*FIR_cA(N1)*FIR_cB(N2)+two*FIZ_cA(N1)*FIZ_cB(N2) &
                            -y*y*FIA1FIB2*K_qrpa_responce**2 -y*y*FIA1FIB2 )*ovlapr         !! (Delta s)_r
                  TEMP8_c = TEMP8_c -two*IUnit*K_qrpa_responce*y*y*FIA1FIB2*ovlapr          !! (Delta s)_phi term: (2/r^2)(dv_r/dphi)
                  !
                  TEMP14_c = TEMP14_c + (FIZ_cA(N1)*FI_cB(N2)+FI_cA(N1)*FIZ_cB(N2))*ovlapr  !! (curl s)_phi
                  TEMP15_c = TEMP15_c + IUnit*K_qrpa_responce*y*FIA1FIB2*ovlapr             !! (curl s)_z
                  TEMP16_c = ( FIR_cA(N1)*FIR_cB(N2) +FIZ_cA(N1)*FIZ_cB(N2)     &
                              -y*y*FIA1FIB2*LAM1(N1)*LAM2(N2) )*ovlapr                      !! T_r
                end if
                if(LAM1(n1) + LAM2(n2) +1 .eq. -abs(K_qrpa_responce)) then                
                  TEMP4_c = -FIA1FIB2*ovlapr                                                !! s_r
                  !
                  TEMP7_c = -(FID2_cA(N1)*FI_cB(N2)+FI_cA(N1)*FID2_cB(N2)+two*FIR_cA(N1)*FIR_cB(N2)+two*FIZ_cA(N1)*FIZ_cB(N2) &
                            -y*y*FIA1FIB2*K_qrpa_responce**2 -y*y*FIA1FIB2 )*ovlapr         !! (Delta s)_r
                  TEMP8_c = TEMP8_c +two*IUnit*K_qrpa_responce*y*y*FIA1FIB2*ovlapr          !! (Delta s)_phi term: (2/r^2)(dv_r/dphi)
                  !
                  TEMP14_c = TEMP14_c - (FIZ_cA(N1)*FI_cB(N2)+FI_cA(N1)*FIZ_cB(N2))*ovlapr  !! (curl s)_phi
                  TEMP15_c = TEMP15_c - IUnit*K_qrpa_responce*y*FIA1FIB2*ovlapr             !! (curl s)_z 
                  TEMP16_c = -( FIR_cA(N1)*FIR_cB(N2) +FIZ_cA(N1)*FIZ_cB(N2)     &
                               -y*y*FIA1FIB2*LAM1(N1)*LAM2(N2) )*ovlapr                     !! T_r
                end if
                !
                TEMP7_c = TEMP7_c +two*IUnit*K_qrpa_responce*y*y*FIA1FIB2*ovlapphi !! (Delta s)_r term: (-2/r^2)(dv_phi/dphi)
                !
                ! these densities have (e^iKphi - e^-iKphi) structure
                ! here e^+iKphi part defines the density and sign
                if(LAM1(n1) - LAM2(n2) .gt. 0) then
                  TEMP3_c  = -half*y*FIA1FIB2*(LAM1(N1)+LAM2(N2))*ovlap !! j_phi
                  TEMP6_c  = -FIA1FIB2*ovlapz                           !! s_z
                  TEMP9_c  = -(FID2_cA(N1)*FI_cB(N2)+FI_cA(N1)*FID2_cB(N2)+two*FIR_cA(N1)*FIR_cB(N2)+two*FIZ_cA(N1)*FIZ_cB(N2) &
                             -y*y*FIA1FIB2*K_qrpa_responce**2 )*ovlapz  !! (Delta s)_z
                  !
                  TEMP10_c = TEMP10_c +half*y*(FIZ_cA(N1)*FI_cB(N2)+FI_cA(N1)*FIZ_cB(N2))*(LAM1(N1)+LAM2(N2))*ovlap     !! (curl j)_r
                  TEMP12_c = TEMP12_c -half*y*( +FIR_cA(N1)*FI_cB(N2) +FI_cA(N1)*FIR_cB(N2))*(LAM1(N1)+LAM2(N2))*ovlap  !! (curl j)_z
                  !
                  TEMP13_c = TEMP13_c + IUnit*K_qrpa_responce*y*FIA1FIB2*ovlapz             !! (curl s)_r
                  TEMP14_c = TEMP14_c + (FIR_cA(N1)*FI_cB(N2) +FI_cA(N1)*FIR_cB(N2))*ovlapz !! (curl s)_phi
                  TEMP18_c = -( FIR_cA(N1)*FIR_cB(N2) +FIZ_cA(N1)*FIZ_cB(N2)     &
                               +y*y*FIA1FIB2*LAM1(N1)*LAM2(N2) )*ovlapz                     !! T_z
                else 
                  TEMP3_c  = half*y*FIA1FIB2*(LAM1(N1)+LAM2(N2))*ovlap !! j_phi
                  TEMP6_c  = FIA1FIB2*ovlapz                           !! s_z
                  TEMP9_c  = (FID2_cA(N1)*FI_cB(N2)+FI_cA(N1)*FID2_cB(N2)+two*FIR_cA(N1)*FIR_cB(N2)+two*FIZ_cA(N1)*FIZ_cB(N2) &
                             -y*y*FIA1FIB2*K_qrpa_responce**2 )*ovlapz !! (Delta s)_z
                  !
                  TEMP10_c = TEMP10_c -half*y*(FIZ_cA(N1)*FI_cB(N2)+FI_cA(N1)*FIZ_cB(N2))*(LAM1(N1)+LAM2(N2))*ovlap    !! (curl j)_r
                  TEMP12_c = TEMP12_c +half*y*( +FIR_cA(N1)*FI_cB(N2) +FI_cA(N1)*FIR_cB(N2))*(LAM1(N1)+LAM2(N2))*ovlap !! (curl j)_z
                  !
                  TEMP13_c = TEMP13_c - IUnit*K_qrpa_responce*y*FIA1FIB2*ovlapz             !! (curl s)_r
                  TEMP14_c = TEMP14_c - (FIR_cA(N1)*FI_cB(N2) +FI_cA(N1)*FIR_cB(N2))*ovlapz !! (curl s)_phi
                  TEMP18_c = ( FIR_cA(N1)*FIR_cB(N2) +FIZ_cA(N1)*FIZ_cB(N2)     &
                              +y*y*FIA1FIB2*LAM1(N1)*LAM2(N2) )*ovlapz                      !! T_z
                end if
                !
                ! these densities have (e^iKphi - e^-iKphi) structure
                ! here e^+iKphi part defines the density and sign
                TEMP10_c = TEMP10_c -IUnit*K_qrpa_responce*y*(-iunit*( FIZ_cA(N1)*FI_cB(N2) - FI_cA(N1)*FIZ_cB(N2) ))*half*ovlap !! (curl j)_r
                TEMP12_c = TEMP12_c +IUnit*K_qrpa_responce*y*(-iunit*( FIR_cA(N1)*FI_cB(N2) - FI_cA(N1)*FIR_cB(N2) ))*half*ovlap !! (curl j)_z
                !
                ! correct time-odd densities when using normal HO basis
                If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) then
                 !! time-odd densities are already correct
                end if
                !
                ! time-odd matter densities
                Prjr_c(ihil)    = Prjr_c(ihil) +TEMP1_c*rhme         !! j_r
                Prjz_c(ihil)    = Prjz_c(ihil) +TEMP2_c*rhme         !! j_z
                Psofi_c(ihil)   = Psofi_c(ihil) +TEMP5_c*rhme        !! s_fi
                Pdsofi_c(ihil)  = Pdsofi_c(ihil) +TEMP8_c*rhme       !! (Delta s)_fi
                Pcurjfi_c(ihil) = Pcurjfi_c(ihil) -TEMP11_c*rhme     !! (curl j)_fi
                Pcursr_c(ihil)  = Pcursr_c(ihil) -TEMP13_c*rhme      !! (curl s)_r
                Pcursz_c(ihil)  = Pcursz_c(ihil) -TEMP15_c*rhme      !! (curl s)_z
                Ptsofi_c(ihil)  = Ptsofi_c(ihil) +TEMP17_c*rhme      !! T_fi
                if(K_qrpa_responce.ne.0) then 
                  Prjfi_c(ihil)   = Prjfi_c(ihil) +TEMP3_c*rhme      !! j_fi
                  Psor_c(ihil)    = Psor_c(ihil) +TEMP4_c*rhme       !! s_r
                  Psoz_c(ihil)    = Psoz_c(ihil) +TEMP6_c*rhme       !! s_z
                  Pdsor_c(ihil)   = Pdsor_c(ihil) +TEMP7_c*rhme      !! (Delta s)_r
                  Pdsoz_c(ihil)   = Pdsoz_c(ihil) +TEMP9_c*rhme      !! (Delta s)_z
                  Pcurjr_c(ihil)  = Pcurjr_c(ihil) -TEMP10_c*rhme    !! (curl j)_r
                  Pcurjz_c(ihil)  = Pcurjz_c(ihil) -TEMP12_c*rhme    !! (curl j)_z
                  Pcursfi_c(ihil) = Pcursfi_c(ihil) -TEMP14_c*rhme   !! (curl s)_fi
                  Ptsor_c(ihil)   = Ptsor_c(ihil) +TEMP16_c*rhme     !! T_r
                  Ptsoz_c(ihil)   = Ptsoz_c(ihil) +TEMP18_c*rhme     !! T_z
                end if
                !
             End Do
           End Do

         End Do !! qrpa blocks
         !
       End DO  !! ihil
       !----------------------------------------------------
       ! REMOVES INT.WEIGHTS AND MULTIPLIES BY THE JACOBIAN
       !----------------------------------------------------
       snz_qrpa(it)=Two*Sum(TRO_c); dsnz_qrpa(it)=Two*Sum(TDRO_c)
       Do ihil=1,nghl      
          ss=wdcori(ihil)
          if(K_QRPA_RESPONCE .ne. 0) ss = ss*half !! normalization
          !----------------------------------
          ! TIME-EVEN PART
          !----------------------------------
          TAKA_c(ihil)    = TAKA_c(ihil)*ss
          TBKA_c(ihil)    = TBKA_c(ihil)*ss
          !
          ss=ss*Two
          TJRR_c(ihil)     = TJRR_c(ihil)*ss     ! J_r,r
          TJRFI_c(ihil)    = TJRFI_c(ihil)*ss    ! J_r,fi
          TJRZ_c(ihil)     = TJRZ_c(ihil)*ss     ! J_r,z
          TJFIR_c(ihil)    = TJFIR_c(ihil)*ss    ! J_fi,r
          TJFIFI_c(ihil)   = TJFIFI_c(ihil)*ss   ! J_fi,fi
          TJFIZ_c(ihil)    = TJFIZ_c(ihil)*ss    ! J_fi,z
          TJZR_c(ihil)     = TJZR_c(ihil)*ss     ! J_z,r
          TJZFI_c(ihil)    = TJZFI_c(ihil)*ss    ! J_z,fi
          TJZZ_c(ihil)     = TJZZ_c(ihil)*ss     ! J_z,z
          !
          TRO_c(ihil)     = TRO_c(ihil)*ss
          TDRO_c(ihil)    = TDRO_c(ihil)*ss
          TTAU_c(ihil)    = TTAU_c(ihil)*ss
          TDJ_c(ihil)     = TDJ_c(ihil)*ss
          TNABLAR_c(IHIL) = TNABLAR_c(IHIL)*ss
          TNABLAZ_c(IHIL) = TNABLAZ_c(IHIL)*ss
          !----------------------------------
          ! TIME-ODD PART
          !----------------------------------		  
          Prjr_c(ihil)     = Prjr_c(ihil)*ss     ! j_r
          Prjfi_c(ihil)    = Prjfi_c(ihil)*ss    ! j_fi
          Prjz_c(ihil)     = Prjz_c(ihil)*ss     ! j_z
          Psor_c(ihil)     = Psor_c(ihil)*ss     ! s_r
          Psofi_c(ihil)    = Psofi_c(ihil)*ss    ! s_fi          
          Psoz_c(ihil)     = Psoz_c(ihil)*ss     ! s_z
          Pcursr_c(ihil)   = Pcursr_c(ihil)*ss   ! cur s
          Pcursfi_c(ihil)  = Pcursfi_c(ihil)*ss  ! cur s
          Pcursz_c(ihil)   = Pcursz_c(ihil)*ss   ! cur s
          Pdsor_c(ihil)    = Pdsor_c(ihil)*ss    ! (laplacian s)_r
          Pdsofi_c(ihil)   = Pdsofi_c(ihil)*ss   ! (laplacian s)_fi
          Pdsoz_c(ihil)    = Pdsoz_c(ihil)*ss    ! (laplacian s)_z
          Pcurjr_c(ihil)   = Pcurjr_c(ihil)*ss   ! (cur j)_r
          Pcurjfi_c(ihil)  = Pcurjfi_c(ihil)*ss  ! (cur j)_fi
          Pcurjz_c(ihil)   = Pcurjz_c(ihil)*ss   ! (cur j)_z
          Ptsor_c(ihil)    = Ptsor_c(ihil)*ss    ! T_r
          Ptsofi_c(ihil)   = Ptsofi_c(ihil)*ss   ! T_fi
          Ptsoz_c(ihil)    = Ptsoz_c(ihil)*ss    ! T_z
       Enddo
    End Do !it
    !
    Call get_CPU_time('qrpa_DENSIT',1)
    !
  End Subroutine qrpa_DENSIT
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_field
    !---------------------------------------------------------------------
    ! calculates fields in r-space form axially symmetric densities
    ! NB! includes the external field
    !---------------------------------------------------------------------
    Implicit None
    Integer(ipr)   :: iw,it,ita,ihli
    Real(pr)       :: eps
    !
    Complex(pr)    :: RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,DJ_0,DJ_1,davn0,dbvn0
    Complex(pr)    :: JRR_0,JRR_1,JRFI_0,JRFI_1,JRZ_0,JRZ_1,JFIR_0,JFIR_1,JFIFI_0,JFIFI_1
    Complex(pr)    :: JFIZ_0,JFIZ_1,JZR_0,JZR_1,JZFI_0,JZFI_1,JZZ_0,JZZ_1
    Complex(pr)    :: SNABLARN,SNABLAZN,SNABLARP,SNABLAZP
    Complex(pr)    :: SNABLAR_0,SNABLAZ_0,SNABLAR_1,SNABLAZ_1
    Complex(pr)    :: J2_0,J2_1,rsa0
    ! time-even fields
    Complex(pr)    :: pUr(2),pUt(2),pUNr(2),pUNz(2),pUDr(2),pUDj(2)
    Complex(pr)    :: tUr(2),tUt(2),tUNr(2),tUNz(2),tUDr(2),tUDj(2)
    Complex(pr)    :: pUJrr(2),pUJrfi(2),pUJrz(2),pUJfir(2),pUJfifi(2),pUJfiz(2),pUJzr(2),pUJzfi(2),pUJzz(2)
    Complex(pr)    :: tUJrr(2),tUJrfi(2),tUJrz(2),tUJfir(2),tUJfifi(2),tUJfiz(2),tUJzr(2),tUJzfi(2),tUJzz(2)
    !
    Complex(pr)    :: RHO_0_HFB, RHO_1_HFB, TAU_0_HFB, TAU_1_HFB, KAP_P_HFB, KAP_N_HFB
    !
    Complex(pr)    :: SPR_0,SPR_1,SPZ_0,SPZ_1,SPFI_0,SPFI_1
    Complex(pr)    :: DSPR_0,DSPR_1,DSPZ_0,DSPZ_1,DSPFI_0,DSPFI_1
    Complex(pr)    :: CJCr_0,CJCr_1,CJCz_0,CJCz_1,CJCfi_0,CJCfi_1
    Complex(pr)    :: Csr_0,Csr_1,Csz_0,Csz_1,Csfi_0,Csfi_1
    Complex(pr)    :: JCR_0,JCR_1,JCZ_0,JCZ_1,JCFI_0,JCFI_1
    Complex(pr)    :: Tr_0,Tr_1,Tfi_0,Tfi_1,Tz_0,Tz_1
    ! time-odd fields
    Complex(pr)    :: pSir(2),pSiz(2),pSifi(2),pSidr(2),pSidz(2),pSidfi(2),pjcr(2),pjcz(2),pjcfi(2),crpj(2),czpj(2),cfipj(2)
    Complex(pr)    :: tSir(2),tSiz(2),tSifi(2),tSidr(2),tSidz(2),tSidfi(2),tjcr(2),tjcz(2),tjcfi(2),crtj(2),cztj(2),cfitj(2)
    Complex(pr)    :: crts(2),czts(2),cfits(2),crps(2),czps(2),cfips(2)
    complex(pr)    :: tTTr(2),tTTz(2),tTTfi(2)
    complex(pr)    :: pTTr(2),pTTz(2),pTTfi(2)
    !
    Call get_CPU_time('qrpa_field',0)
    !
    ! fields
    eps=Spacing(1.0_pr);
    Do ihli=1,nghl
       !------------------------------------------------------------------------------------
       ! TIME-EVEN PART
       !------------------------------------------------------------------------------------
       !
       ! form the isoscalar and -vector densities from the arrays containin p-n densities
       RHO_0     = ro_c(ihli,1)     + ro_c(ihli,2);         RHO_1     = ro_c(ihli,1)     - ro_c(ihli,2)
       TAU_0     = tau_c(ihli,1)    + tau_c(ihli,2);        TAU_1     = tau_c(ihli,1)    - tau_c(ihli,2)
       DRHO_0    = dro_c(ihli,1)    + dro_c(ihli,2);        DRHO_1    = dro_c(ihli,1)    - dro_c(ihli,2)
       DJ_0      = dj_c(ihli,1)     + dj_c(ihli,2);         DJ_1      = dj_c(ihli,1)     - dj_c(ihli,2)
       JRR_0     = Jrr_c(ihli,1)    + Jrr_c(ihli,2);        JRR_1     = Jrr_c(ihli,1)    - Jrr_c(ihli,2)
       JRFI_0    = Jrfi_c(ihli,1)   + Jrfi_c(ihli,2);       JRFI_1    = Jrfi_c(ihli,1)   - Jrfi_c(ihli,2)
       JRZ_0     = Jrz_c(ihli,1)    + Jrz_c(ihli,2);        JRZ_1     = Jrz_c(ihli,1)    - Jrz_c(ihli,2)
       JFIR_0    = Jfir_c(ihli,1)   + Jfir_c(ihli,2);       JFIR_1    = Jfir_c(ihli,1)   - Jfir_c(ihli,2)
       JFIFI_0   = Jfifi_c(ihli,1)  + Jfifi_c(ihli,2);      JFIFI_1   = Jfifi_c(ihli,1)  - Jfifi_c(ihli,2)
       JFIZ_0    = Jfiz_c(ihli,1)   + Jfiz_c(ihli,2);       JFIZ_1    = Jfiz_c(ihli,1)   - Jfiz_c(ihli,2)
       JZR_0     = Jzr_c(ihli,1)    + Jzr_c(ihli,2);        JZR_1     = Jzr_c(ihli,1)    - Jzr_c(ihli,2)
       JZFI_0    = Jzfi_c(ihli,1)   + Jzfi_c(ihli,2);       JZFI_1    = Jzfi_c(ihli,1)   - Jzfi_c(ihli,2)
       JZZ_0     = Jzz_c(ihli,1)    + Jzz_c(ihli,2);        JZZ_1     = Jzz_c(ihli,1)    - Jzz_c(ihli,2)
       SNABLAR_0 = NABLAR_c(ihli,1) + NABLAR_c(ihli,2);     SNABLAR_1 = NABLAR_c(ihli,1) - NABLAR_c(ihli,2)
       SNABLAZ_0 = NABLAZ_c(ihli,1) + NABLAZ_c(ihli,2);     SNABLAZ_1 = NABLAZ_c(ihli,1) - NABLAZ_c(ihli,2)
       !
       RHO_0_HFB = ro_c_hfb(ihli,1)  + ro_c_hfb(ihli,2)  ;  RHO_1_HFB = ro_c_hfb(ihli,1)  - ro_c_hfb(ihli,2)     ! HFB rho_0 and rho_1 density
       TAU_0_HFB = tau_c_hfb(ihli,1) + tau_c_hfb(ihli,2) ;  TAU_1_HFB = tau_c_hfb(ihli,1) - tau_c_hfb(ihli,2)    ! HFB tau_0 and tau_1 density
       KAP_N_HFB = kap_c_hfb(ihli,1)                     ;  KAP_P_HFB = kap_c_hfb(ihli,2)                        ! HFB pairing densities
       !    
       ! At the moment basic Skyrme EDF is assumed with no extra density dependence besides the usual ones
       !
       ! initialize to zero these isospin indexed fields
       tUr=zero; tUDr=zero; tUNr=zero; tUNz=zero; tUt=zero; tUDj=zero; 
       tUJrr=zero; tUJrfi=zero; tUJrz=zero; tUJfir=zero; tUJfifi=zero;
       tUJfiz=zero; tUJzr=zero; tUJzfi=zero; tUJzz=zero;
       !
       ! rho^2 terms. Here the fields have been explicitly linearized for density dependent part
       tUr(1)=tUr(1)+two*Crho(0)*RHO_0 +(2.0_pr+sigma)*(1.0_pr+sigma)*CDrho(0)*(RHO_0_HFB**sigma)*RHO_0 &
                    +CDrho(1)*sigma*two*RHO_1*RHO_1_HFB*((RHO_0_HFB+eps)**(sigma-1.0_pr)) &
                    +CDrho(1)*sigma*(sigma-1.0_pr)*RHO_0*(RHO_1_HFB**2)*((RHO_0_HFB+eps)**(sigma-2.0_pr)) 
       tUr(2)=tUr(2)+two*Crho(1)*RHO_1 +two*CDrho(1)*(RHO_0_HFB**sigma)*RHO_1 &
                    +two*CDrho(1)*sigma*RHO_0*RHO_1_HFB*((RHO_0_HFB+eps)**(sigma-1.0_pr))
       !
       ! rho tau
       tUr(1)=tUr(1)+Ctau(0)*TAU_0
       tUr(2)=tUr(2)+Ctau(1)*TAU_1
       tUt(1)=tUt(1)+Ctau(0)*RHO_0
       tUt(2)=tUt(2)+Ctau(1)*RHO_1
       !
       ! rho Delta rho
       tUr(1)=tUr(1)+Crdr(0)*DRHO_0
       tUr(2)=tUr(2)+Crdr(1)*DRHO_1
       tUDr(1)=tUDr(1)+Crdr(0)*RHO_0
       tUDr(2)=tUDr(2)+Crdr(1)*RHO_1
       !
       ! rho nabla J
       tUr(1)=tUr(1)+Crdj(0)*DJ_0
       tUr(2)=tUr(2)+Crdj(1)*DJ_1
       tUdj(1)=tUdj(1)+Crdj(0)*RHO_0
       tUdj(2)=tUdj(2)+Crdj(1)*RHO_1
       !
       ! J^2
       tUJrr(1)   = tUJrr(1)   + CJ(0)*two*Jrr_0;     tUJrr(2) = tUJrr(2)   + CJ(1)*two*Jrr_1;
       tUJrfi(1)  = tUJrfi(1)  + CJ(0)*two*Jrfi_0;   tUJrfi(2) = tUJrfi(2)  + CJ(1)*two*Jrfi_1;
       tUJrz(1)   = tUJrz(1)   + CJ(0)*two*Jrz_0;     tUJrz(2) = tUJrz(2)   + CJ(1)*two*Jrz_1;
       tUJfir(1)  = tUJfir(1)  + CJ(0)*two*Jfir_0;   tUJfir(2) = tUJfir(2)  + CJ(1)*two*Jfir_1;
       tUJfifi(1) = tUJfifi(1) + CJ(0)*two*Jfifi_0; tUJfifi(2) = tUJfifi(2) + CJ(1)*two*Jfifi_1;
       tUJfiz(1)  = tUJfiz(1)  + CJ(0)*two*Jfiz_0;   tUJfiz(2) = tUJfiz(2)  + CJ(1)*two*Jfiz_1;
       tUJzr(1)   = tUJzr(1)   + CJ(0)*two*Jzr_0;     tUJzr(2) = tUJzr(2)   + CJ(1)*two*Jzr_1;
       tUJzfi(1)  = tUJzfi(1)  + CJ(0)*two*Jzfi_0;   tUJzfi(2) = tUJzfi(2)  + CJ(1)*two*Jzfi_1;
       tUJzz(1)   = tUJzz(1)   + CJ(0)*two*Jzz_0;     tUJzz(2) = tUJzz(2)   + CJ(1)*two*Jzz_1;
       !
       ! Proton-neutron representation
       pUr(1)  =tUr(1)   + tUr(2);       pUr(2)  =tUr(1)   - tUr(2)
       pUt(1)  =tUt(1)   + tUt(2);       pUt(2)  =tUt(1)   - tUt(2)
       pUDr(1) =tUDr(1)  + tUDr(2);      pUDr(2) =tUDr(1)  - tUDr(2)
       pUdj(1) =tUdj(1)  + tUdj(2);      pUdj(2) =tUdj(1)  - tUdj(2)
       pUNr(1) =tUNr(1)  + tUNr(2);      pUNr(2) =tUNr(1)  - tUNr(2)
       pUNz(1) =tUNz(1)  + tUNz(2);      pUNz(2) =tUNz(1)  - tUNz(2)
       !
       pUJrr(1)   = tUJrr(1)   + tUJrr(2);    pUJrr(2)   = tUJrr(1)   - tUJrr(2);
       pUJrfi(1)  = tUJrfi(1)  + tUJrfi(2);   pUJrfi(2)  = tUJrfi(1)  - tUJrfi(2);
       pUJrz(1)   = tUJrz(1)   + tUJrz(2);    pUJrz(2)   = tUJrz(1)   - tUJrz(2);
       pUJfir(1)  = tUJfir(1)  + tUJfir(2);   pUJfir(2)  = tUJfir(1)  - tUJfir(2);
       pUJfifi(1) = tUJfifi(1) + tUJfifi(2);  pUJfifi(2) = tUJfifi(1) - tUJfifi(2);
       pUJfiz(1)  = tUJfiz(1)  + tUJfiz(2);   pUJfiz(2)  = tUJfiz(1)  - tUJfiz(2);
       pUJzr(1)   = tUJzr(1)   + tUJzr(2);    pUJzr(2)   = tUJzr(1)   - tUJzr(2);
       pUJzfi(1)  = tUJzfi(1)  + tUJzfi(2);   pUJzfi(2)  = tUJzfi(1)  - tUJzfi(2);
       pUJzz(1)   = tUJzz(1)   + tUJzz(2);    pUJzz(2)   = tUJzz(1)   - tUJzz(2);
       !
       Do it=itmin,itmax   !! loop over n  & p
          ita=3-it
          ! coulomb
          If(it.Eq.2) Then 
             If(icou.Ge.1) pUr(it)=pUr(it)+cou_c(ihli) 
             If(icou.Eq.2) pUr(it)=pUr(it)+CExPar*coex*ro_c(ihli,it) / ( 3.0_pr*(eps+ro_c_hfb(ihli,2)**p23))
          Endif
          !------------------------------------------------------------------------------------
          ! PAIRING 		  
          !------------------------------------------------------------------------------------
          !
          ! density independent part
          davn0 = aka_c(ihli,it)*CpV0(it-1)
          dbvn0 = bka_c(ihli,it)*CpV0(it-1)
          !
          ! density dependent part
          davn0 = davn0 - aka_c(ihli,it)*CpV0(it-1)*CpV1(it-1)*RHO_0_HFB/rho_c
          dbvn0 = dbvn0 - bka_c(ihli,it)*CpV0(it-1)*CpV1(it-1)*RHO_0_HFB/rho_c
          !
          ! rearrangement term. density dep. proportional to rho_n+rho_p -> the same term for protons and neutrons
          pUr(it)=pUr(it)-CpV0(0)*CpV1(0)*conjg(KAP_N_HFB*aka_c(ihli,1) +KAP_N_HFB*Conjg(bka_c(ihli,1)))/rho_c &
                         -CpV0(1)*CpV1(1)*conjg(KAP_P_HFB*aka_c(ihli,2) +KAP_P_HFB*Conjg(bka_c(ihli,2)))/rho_c
          !
          If(it.Eq.1) Then 
             davn_c(ihli)=davn0 -CpV0(it-1)*CpV1(it-1)*KAP_N_HFB*RHO_0/rho_c        !! neutron fields
             dbvn_c(ihli)=dbvn0 -CpV0(it-1)*CpV1(it-1)*KAP_N_HFB*Conjg(RHO_0)/rho_c !! neutron fields
          Else
             davp_c(ihli)=davn0 -CpV0(it-1)*CpV1(it-1)*KAP_P_HFB*RHO_0/rho_c        !! proton fields
             dbvp_c(ihli)=dbvn0 -CpV0(it-1)*CpV1(it-1)*KAP_P_HFB*Conjg(RHO_0)/rho_c !! proton fields
          End If
          !------------------------------------------------------------------------------------
       End Do !it
       vn_c(ihli)     = pUr(1);   vp_c(ihli)     = pUr(2)    !* RHO_ij  
       vhbn_c(ihli)   = pUt(1);   vhbp_c(ihli)   = pUt(2)    !* TAU_ij  
       vrn_c(ihli)    = pUNr(1);  vrp_c(ihli)    = pUNr(2)   !* NABLAr RHO__ij  
       vzn_c(ihli)    = pUNz(1);  vzp_c(ihli)    = pUNz(2)   !* NABLAz RHO__ij  
       vdn_c(ihli)    = pUDr(1);  vdp_c(ihli)    = pUDr(2)   !* DELTA RHO_ij  
       vsn_c(ihli)    = pUdj(1);  vsp_c(ihli)    = pUdj(2)   !* NABLA . J__ij  
       !
       vJrrn_c(ihli)   = pUJrr(1);   vJrrp_c(ihli)   = pUJrr(2);  
       vJrfin_c(ihli)  = pUJrfi(1);  vJrfip_c(ihli)  = pUJrfi(2);  
       vJrzn_c(ihli)   = pUJrz(1);   vJrzp_c(ihli)   = pUJrz(2);  
       vJfirn_c(ihli)  = pUJfir(1);  vJfirp_c(ihli)  = pUJfir(2);  
       vJfifin_c(ihli) = pUJfifi(1); vJfifip_c(ihli) = pUJfifi(2);  
       vJfizn_c(ihli)  = pUJfiz(1);  vJfizp_c(ihli)  = pUJfiz(2);  
       vJzrn_c(ihli)   = pUJzr(1);   vJzrp_c(ihli)   = pUJzr(2);  
       vJzfin_c(ihli)  = pUJzfi(1);  vJzfip_c(ihli)  = pUJzfi(2);  
       vJzzn_c(ihli)   = pUJzz(1);   vJzzp_c(ihli)   = pUJzz(2);  
       !
       !------------------------------------------------------------------------------------
       ! * TIME-ODD PART *
       !------------------------------------------------------------------------------------
       !
       ! T-odd densities
       SPR_0   = sor_c(ihli,1)    + sor_c(ihli,2);        SPR_1 = sor_c(ihli,1)    - sor_c(ihli,2)
       SPZ_0   = soz_c(ihli,1)    + soz_c(ihli,2);        SPZ_1 = soz_c(ihli,1)    - soz_c(ihli,2)
       SPFI_0  = sofi_c(ihli,1)   + sofi_c(ihli,2);      SPFI_1 = sofi_c(ihli,1)   - sofi_c(ihli,2)
       DSPR_0  = dsor_c(ihli,1)   + dsor_c(ihli,2);      DSPR_1 = dsor_c(ihli,1)   - dsor_c(ihli,2)
       DSPZ_0  = dsoz_c(ihli,1)   + dsoz_c(ihli,2);      DSPZ_1 = dsoz_c(ihli,1)   - dsoz_c(ihli,2)
       DSPFI_0 = dsofi_c(ihli,1)  + dsofi_c(ihli,2);    DSPFI_1 = dsofi_c(ihli,1)  - dsofi_c(ihli,2)
       JCR_0   = rjr_c(ihli,1)    + rjr_c(ihli,2);        JCR_1 = rjr_c(ihli,1)    - rjr_c(ihli,2)
       JCZ_0   = rjz_c(ihli,1)    + rjz_c(ihli,2);        JCZ_1 = rjz_c(ihli,1)    - rjz_c(ihli,2)
       JCFI_0  = rjfi_c(ihli,1)   + rjfi_c(ihli,2);      JCFI_1 = rjfi_c(ihli,1)   - rjfi_c(ihli,2)
       CJCr_0  = curjr_c(ihli,1)  + curjr_c(ihli,2);     CJCr_1 = curjr_c(ihli,1)  - curjr_c(ihli,2);
       CJCz_0  = curjz_c(ihli,1)  + curjz_c(ihli,2);     CJCz_1 = curjz_c(ihli,1)  - curjz_c(ihli,2);
       CJCfi_0 = curjfi_c(ihli,1) + curjfi_c(ihli,2);   CJCfi_1 = curjfi_c(ihli,1) - curjfi_c(ihli,2);
       Csr_0   = cursr_c(ihli,1)  + cursr_c(ihli,2);      Csr_1 = cursr_c(ihli,1)  - cursr_c(ihli,2);
       Csz_0   = cursz_c(ihli,1)  + cursz_c(ihli,2);      Csz_1 = cursz_c(ihli,1)  - cursz_c(ihli,2);
       Csfi_0  = cursfi_c(ihli,1) + cursfi_c(ihli,2);    Csfi_1 = cursfi_c(ihli,1) - cursfi_c(ihli,2);
       Tr_0    = tsor_c(ihli,1)   + tsor_c(ihli,2);        Tr_1 = tsor_c(ihli,1)   - tsor_c(ihli,2);
       Tfi_0   = tsofi_c(ihli,1)  + tsofi_c(ihli,2);      Tfi_1 = tsofi_c(ihli,1)  - tsofi_c(ihli,2);
       Tz_0    = tsoz_c(ihli,1)   + tsoz_c(ihli,2);        Tz_1 = tsoz_c(ihli,1)   - tsoz_c(ihli,2);
       !
       ! put zero on this and use it for density dependent parts
       tUr=zero 
       ! 
       ! initialize isospin indexed fields to zero
       tSir=zero; tSiz=zero;  tSifi=zero; tSidr=zero; tSidz=zero; tSidfi=zero 
       tjcr=zero; tjcz=zero;  tjcfi=zero; crtj=zero;  cztj=zero;  cfitj=zero  
       crts=zero; czts=zero;  cfits=zero;
       tTTr=zero; tTTfi=zero; tTTz=zero;
       !
       ! j^2 term  (will be integreted with j_ij)
       tjcr(1)=tjcr(1)+two*Cjcjc(0)*JCR_0
       tjcr(2)=tjcr(2)+two*Cjcjc(1)*JCR_1
       tjcfi(1)=tjcfi(1)+two*Cjcjc(0)*JCFI_0
       tjcfi(2)=tjcfi(2)+two*Cjcjc(1)*JCFI_1
       tjcz(1)=tjcz(1)+two*Cjcjc(0)*JCZ_0
       tjcz(2)=tjcz(2)+two*Cjcjc(1)*JCZ_1
       !
       ! C^s s^2 term (will be integreted with s_ij)
       tSir(1) =tSir(1) +two*Css(0)*SPR_0
       tSir(2) =tSir(2) +two*Css(1)*SPR_1
       tSifi(1)=tSifi(1)+two*Css(0)*SPFI_0
       tSifi(2)=tSifi(2)+two*Css(1)*SPFI_1
       tSiz(1) =tSiz(1) +two*Css(0)*SPZ_0
       tSiz(2) =tSiz(2) +two*Css(1)*SPZ_1
       !
       ! CD^s s^2 density dependent term (will be integreted with s_ij)
       tSir(1) =tSir(1) +two*CDss(0)*SPR_0*(RHO_0_HFB**sigma)
       tSir(2) =tSir(2) +two*CDss(1)*SPR_1*(RHO_0_HFB**sigma)
       tSifi(1)=tSifi(1)+two*CDss(0)*SPFI_0*(RHO_0_HFB**sigma)
       tSifi(2)=tSifi(2)+two*CDss(1)*SPFI_1*(RHO_0_HFB**sigma)
       tSiz(1) =tSiz(1) +two*CDss(0)*SPZ_0*(RHO_0_HFB**sigma)
       tSiz(2) =tSiz(2) +two*CDss(1)*SPZ_1*(RHO_0_HFB**sigma)
       !
       ! CD^s s^2 -rearrangment term not needed (s-density zero for HFB)
       !
       ! s.\delta s term (will be integreted with s_ij) 
       tSir(1) =tSir(1) +two*CDeltass(0)*DSPR_0
       tSir(2) =tSir(2) +two*CDeltass(1)*DSPR_1
       tSifi(1)=tSifi(1)+two*CDeltass(0)*DSPFI_0
       tSifi(2)=tSifi(2)+two*CDeltass(1)*DSPFI_1
       tSiz(1) =tSiz(1) +two*CDeltass(0)*DSPZ_0
       tSiz(2) =tSiz(2) +two*CDeltass(1)*DSPZ_1
       !
       ! s dot curl j (will be integreted with j_ij and s_ij) 
       tSir(1) =tSir(1) +Csnablajc(0)*CJCr_0
       tSir(2) =tSir(2) +Csnablajc(1)*CJCr_1
       tSifi(1)=tSifi(1)+Csnablajc(0)*CJCfi_0
       tSifi(2)=tSifi(2)+Csnablajc(1)*CJCfi_1
       tSiz(1) =tSiz(1) +Csnablajc(0)*CJCz_0
       tSiz(2) =tSiz(2) +Csnablajc(1)*CJCz_1
       !
       tjcr(1) =tjcr(1) +Csnablajc(0)*Csr_0
       tjcr(2) =tjcr(2) +Csnablajc(1)*Csr_1
       tjcfi(1)=tjcfi(1)+Csnablajc(0)*Csfi_0
       tjcfi(2)=tjcfi(2)+Csnablajc(1)*Csfi_1
       tjcz(1) =tjcz(1) +Csnablajc(0)*Csz_0
       tjcz(2) =tjcz(2) +Csnablajc(1)*Csz_1
       !       
       ! s dot T term (will be integreted with s_ij and T_ij) 
       tSir(1) =tSir(1) +CsT(0)*Tr_0
       tSir(2) =tSir(2) +CsT(1)*Tr_1
       tSifi(1)=tSifi(1)+CsT(0)*Tfi_0
       tSifi(2)=tSifi(2)+CsT(1)*Tfi_1
       tSiz(1) =tSiz(1) +CsT(0)*Tz_0
       tSiz(2) =tSiz(2) +CsT(1)*Tz_1
       tTTr(1) =tTTr(1) +CsT(0)*SPR_0
       tTTr(2) =tTTr(2) +CsT(1)*SPR_1
       tTTfi(1)=tTTfi(1)+CsT(0)*SPFI_0
       tTTfi(2)=tTTfi(2)+CsT(1)*SPFI_1
       tTTz(1) =tTTz(1) +CsT(0)*SPZ_0
       tTTz(2) =tTTz(2) +CsT(1)*SPZ_1
       !
       ! add here s dot F later
       !
       ! Proton-neutron representation. 
       ! From density dependent part
       pUr(1)  =tUr(1)+tUr(2);        pUr(2)  =tUr(1)  -tUr(2)
       !
       ! delta s
       pSir(1) =tSir(1)+tSir(2);      pSir(2) =tSir(1)-tSir(2)
       pSifi(1) =tSifi(1)+tSifi(2);   pSifi(2) =tSifi(1)-tSifi(2)
       pSiz(1) =tSiz(1)+tSiz(2);      pSiz(2) =tSiz(1)-tSiz(2)
       !
       ! delta j
       pjcr(1) =tjcr(1)+tjcr(2);      pjcr(2) =tjcr(1)-tjcr(2)
       pjcz(1) =tjcz(1)+tjcz(2);      pjcz(2) =tjcz(1)-tjcz(2)
       pjcfi(1)=tjcfi(1)+tjcfi(2);    pjcfi(2)=tjcfi(1)-tjcfi(2)  
       !	   
       ! delta T
       pTTr(1)=tTTr(1)+tTTr(2);       pTTr(2)=tTTr(1)-tTTr(2);
       pTTfi(1)=tTTfi(1)+tTTfi(2);    pTTfi(2)=tTTfi(1)-tTTfi(2);
       pTTz(1)=tTTz(1)+tTTz(2);       pTTz(2)=tTTz(1)-tTTz(2);
       !
       ! Add to the time-even field arrays
       vn_c(ihli)=vn_c(ihli)+pUr(1); vp_c(ihli)=vp_c(ihli)+pUr(2)        !* RHO_ij  
       !
       ! Set to time-odd field arrays
       vSIGRn_c(ihli)  = pSir(1)  ;  vSIGRp_c(ihli)  = pSir(2)           !* s_r_ij 
       vSIGZn_c(ihli)  = pSiz(1)  ;  vSIGZp_c(ihli)  = pSiz(2)           !* s_z_ij
       vSIGFIn_c(ihli) = pSifi(1) ;  vSIGFIp_c(ihli) = pSifi(2)          !* s_fi_ij
       !
       vjRn_c(ihli)    = pjcr(1)  ;  vjRp_c(ihli)    = pjcr(2)           !* j_r_ij
       vjZn_c(ihli)    = pjcz(1)  ;  vjZp_c(ihli)    = pjcz(2)           !* j_z_ij
       vjFIn_c(ihli)   = pjcfi(1) ;  vjFIp_c(ihli)   = pjcfi(2)          !* j_fi_ij
       !
       vTRn_c(ihli)    = pTTr(1)  ;  vTRp_c(ihli)    = pTTr(2);          !* T_r_ij
       vTZn_c(ihli)    = pTTz(1)  ;  vTZp_c(ihli)    = pTTz(2);          !* T_z_ij
       vTFIn_c(ihli)   = pTTfi(1) ;  vTFIp_c(ihli)   = pTTfi(2);         !* T_fi_ij
       !
       ! add here later F -terms
       !
    End Do !ihli
    !
    Call get_CPU_time('qrpa_field',1)
    !
  End Subroutine qrpa_field
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_gamdel
    !---------------------------------------------------------------------
    ! Complex qrpa matrices BhN/P_c and BdN/P_c in configurational space 
    !---------------------------------------------------------------------
    Implicit None
    Integer(ipr)    :: i,j,ihil,nd,im,nsa,JA,N1,N2,ibq
    Integer(ipr)    :: dim1,dim2,ovlap
    complex(pr)     :: vh,vhb,vdel,vjr,vjz,vso
    complex(pr)     :: vnabr,vnabz,vsfiz,vszfi,vsfir,vsrfi,vspinfi,vspinr
    complex(pr)     :: vsob,vspinfib,vspinrb,vspinz,vspinzb,vjfi
    complex(pr)     :: vjrf,vjrfb,vjfr,vjfrb,vjfz,vjfzb,vjzf,vjzfb
    complex(pr)     :: vjrr,vjrrb,vjrz,vjrzb,vjff,vjffb,vjzr,vjzrb,vjzz,vjzzb
    complex(pr)     :: vttbase,vttr,vttfi,vttz, vttrb,vttfib,vttzb
    complex(pr)     :: ovlapr,ovlapphi,ovlapz
    !
    Integer(ipr), allocatable :: lam1(:),lam2(:),ss1(:),ss2(:)
    Real(pr), allocatable     :: FI1(:),FI2(:),FIR1(:),FIR2(:),FIZ1(:),FIZ2(:),FIDE1(:),FIDE2(:)
    !
    Real(pr) :: y !!,XLAPI,XLAMI,XLAPF,XLAMF, XLAPIY,XLAMIY,XLAPFY,XLAMFY
    !
    Complex(pr)  :: iunit,iunit2
    !
    Parameter(iunit = Cmplx(0.0_pr,1.0_pr,kind=pr), iunit2 = iunit/two)
    !
    Call get_CPU_time('qrpa_gamdel',0)
    !
    hN_c=zero; hP_c=zero; BhN_c=zero; BhP_c=zero; dN_c=zero; dP_c=zero; BdN_c=zero; BdP_c=zero
    !
    !----------------------------------------------
    ! START BLOCKS
    !----------------------------------------------
    Do ibq=1,iQRPA_NB    !! blocks
      !
      dim1 = iQRPA_IDF(ibq) ; dim2 = iQRPA_IDI(ibq)
      !
      if(allocated(FI1)) Deallocate(FI1,FI2,FIR1,FIR2,FIZ1,FIZ2,FIDE1,FIDE2)
      allocate(FI1(dim1),FI2(dim2),FIR1(dim1),FIR2(dim2),FIZ1(dim1),FIZ2(dim2),FIDE1(dim1),FIDE2(dim2))
      !
      if(allocated(lam1)) Deallocate(lam1,ss1,lam2,ss2)
      allocate(lam1(dim1),lam2(dim2),ss1(dim1),ss2(dim2))
      !
      Do ihil=1,nghl     !! mesh
        !
        ! run trough the left basis vectors and fill the array
        Do N1=1,iQRPA_IDF(ibq)
          IM=ia(iQRPA_BLF(ibq)); JA=IM+N1; NSA=NS(JA);
          if(NSA.gt.0) then
            LAM1(N1) = (iQRPA_OMF(ibq) -1)/2
          else
            LAM1(N1) = -(iQRPA_OMF(ibq) +1)/2
          end if
          If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) LAM1(N1) = abs(LAM1(N1))  !! normal HO, positive omega states have positive lambda
          FI1(N1)  = QHLA_opt(JA,ihil)
          FIR1(N1) = FI1R_opt(JA,ihil)
          FIZ1(N1) = FI1Z_opt(JA,ihil)
          FIDE1(N1)= FI2D_opt(JA,ihil) -FI1(N1)*((y_opt(ihil)*LAM1(N1))**2)
          SS1(N1) = NSA
        End Do
        !
        ! run trough the right basis vectors and fill the array
        Do N1=1,iQRPA_IDI(ibq)
          IM=ia(iQRPA_BLI(ibq)); JA=IM+N1; NSA=NS(JA);
          if(NSA.gt.0) then
            LAM2(N1) = (iQRPA_OMI(ibq) -1)/2
          else
            LAM2(N1) = -(iQRPA_OMI(ibq) +1)/2
          end if
          If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) LAM2(N1) = abs(LAM2(N1))  !! normal HO, positive omega states have positive lambda
          FI2(N1)  = QHLA_opt(JA,ihil)
          FIR2(N1) = FI1R_opt(JA,ihil)
          FIZ2(N1) = FI1Z_opt(JA,ihil)
          FIDE2(N1)=FI2D_opt(JA,ihil) -FI2(N1)*((y_opt(ihil)*LAM2(N1))**2)
          SS2(N1) = NSA
        End Do
        !
        y=y_opt(ihil);
        !
        j=iQRPA_DS(ibq)
        !
        Do N2=1,iQRPA_IDI(ibq) !! states in the block
          Do N1=1,iQRPA_IDF(ibq)   !! states in the block
             !
             ovlap = 0;
             if( abs(LAM1(n1) + LAM2(n2)) .eq. abs(K_qrpa_responce) ) ovlap = 1  !! matching lambdas
             if( abs(LAM1(n1) - LAM2(n2)) .eq. abs(K_qrpa_responce) ) ovlap = 1  !! 
             !
             ovlapphi = 0; ovlapr = 0; ovlapz = 0;
             ovlapz = -ovlap;
             ovlapr = (1-ovlap)*(-IUnit)
             ovlapphi = (1-ovlap)*(-1)
             !
             ! define overlaps for normal HO basis. 
             If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) then
                if(LAM1(n1) .eq. LAM2(n2)) then ; ovlap = 1 ; else ; ovlap = 0 ; end if
                ovlapz = ovlap
                ovlapr = (1-ovlap)
                ovlapphi = (-1)*(1-ovlap)*Iunit                    !! opposite sign here compared to DENSIT
                if( SS1(n1).lt.0 ) ovlapz = ovlapz * (-1.0_pr)
                if( SS2(n2).lt.0 ) ovlapphi = ovlapphi * (-1.0_pr)
             end if

             ! some aux. variables, time even
             vh = FI1(N1)*FI2(N2) *ovlap                                                                   !! RHO_ij
             vhb = ( FIR1(N1)*FIR2(N2) +FIZ1(N1)*FIZ2(N2) +LAM1(N1)*LAM2(N2)*y*y*FI1(N1)*FI2(N2) )*ovlap   !! TAU_ij
             vdel = ( FIDE1(N1)*FI2(N2) + FI1(N1)*FIDE2(N2))*ovlap +two*vhb                                !! Laplacian Rho_ij
             vso = (-( IUnit*LAM2(N2)*y*FIR1(N1)*FI2(N2)  + IUnit*LAM1(N1)*y*FI1(N1)*FIR2(N2) ) *ovlapz &  !! nabla.J_ij
                    -( -IUnit*LAM1(N1)*y*FI1(N1)*FIZ2(N2) + IUnit*LAM2(n2)*y*FIZ1(N1)*FI2(N2) ) *ovlapr  &
                    -( FIZ1(N1)*FIR2(N2) - FIR1(N1)*FIZ2(N2) )*ovlapphi )/IUnit
             vsob= (-( IUnit*LAM1(N1)*y*FIR2(N2)*FI1(N1)  + IUnit*LAM2(N2)*y*FI2(N2)*FIR1(N1) ) *ovlapz &   !! nabla.J_ij
                    -( -IUnit*LAM2(N2)*y*FI2(N2)*FIZ1(N1) + IUnit*LAM1(N1)*y*FIZ2(N2)*FI1(N1) ) *ovlapr  &  !! for simplex -i states
                    -( FIZ2(N2)*FIR1(N1) - FIR2(N2)*FIZ1(N1) )*ovlapphi )/IUnit
             vjrf=-iunit*( FIR1(N1)*FI2(N2) -FI1(N1)*FIR2(N2) )*half*ovlapphi             !! J_r,phi_ij
             vjfr=-iunit*( iunit*(LAM1(N1)-LAM2(N2))*y*FI1(N1)*FI2(N2) )*half*ovlapr      !! J_phi,r_ij
             vjfz=-iunit*( -iunit*(LAM1(N1)+LAM2(N2))*y*FI1(N1)*FI2(N2) )*half*ovlapz     !! J_phi,z_ij
             vjzf=-iunit*( FIZ1(N1)*FI2(N2) -FI1(N1)*FIZ2(N2) )*half*ovlapphi             !! J_z,phi_ij
             !
             ! tensor components with only K.ne.0 case
             vjrr=-iunit*( FIR1(N1)*FI2(N2) -FI1(N1)*FIR2(N2) )*half*ovlapr               !! J_r,r_ij
             vjrz=-iunit*( FIR1(N1)*FI2(N2) -FI1(N1)*FIR2(N2) )*half*ovlapz               !! J_r,r_ij
             vjff=-iunit*( iunit*(LAM1(N1)-LAM2(N2))*y*FI1(N1)*FI2(N2))*half*ovlapphi     !! J_phi,phi_ij
             vjzr=-iunit*( FIZ1(N1)*FI2(N2) -FI1(N1)*FIZ2(N2) )*half*ovlapr               !! J_z,r_ij
             vjzz=-iunit*( FIZ1(N1)*FI2(N2) -FI1(N1)*FIZ2(N2) )*half*ovlapz               !! J_z,z_ij
             !
             ! these fields have (e^+iKphi - e^-iKphi) structure, hence the sign
             if(LAM1(n1) + LAM2(n2) +1 .eq. abs(K_qrpa_responce)) then
               vjrr = -vjrr; vjff = -vjff; vjzr = -vjzr;
             end if
             !
             ! these fields have (e^+iKphi - e^-iKphi) structure, hence the sign
             if(LAM1(n1) - LAM2(n2) .gt. 0) then
               vjrz = -vjrz; vjzz = -vjzz
             end if
             !
             ! zero those that don't contribute in K=0 case
             if(K_qrpa_responce.eq.0) then
               vjrr = zero; vjff = zero; vjzr = zero; vjrz = zero; vjzz = zero;
             end if
             !
             ! matrix elements from conjugate states (simplex -i states)
             vjrfb = -vjrf ; vjfrb = -vjfr ; vjfzb =  vjfz ; vjzfb = -vjzf;
             vjrrb = -vjrr ; vjrzb =  vjrz ; vjffb = -vjff ; vjzrb = -vjzr;
             vjzzb =  vjzz ;
             !
             ! correct matrix element terms when using normal HO basis
             If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) then
               vso = (+( IUnit*LAM2(N2)*y*FIR1(N1)*FI2(N2) + IUnit*LAM1(N1)*y*FI1(N1)*FIR2(N2) ) *ovlapz &  !! nabla.J_ij
                      -( IUnit*LAM1(N1)*y*FI1(N1)*FIZ2(N2) + IUnit*LAM2(n2)*y*FIZ1(N1)*FI2(N2) ) *ovlapr &
                      -( FIZ1(N1)*FIR2(N2) - FIR1(N1)*FIZ2(N2) )*ovlapphi )/IUnit
               vsob=vso
               !! add here correction for tensor terms
             end if   
             ! * * time odd part * *
             !
             ! variation w.r.t. j
             vjr = ( FIR1(N1)*FI2(N2) - FI1(N1)*FIR2(N2) )*ovlap  !! imag. unit added later
             vjz = ( FIZ1(N1)*FI2(N2) - FI1(N1)*FIZ2(N2) )*ovlap  !! imag. unit added later
             vjfi = half*y*FI1(N1)*FI2(N2)*(LAM1(N1)+LAM2(N2))*ovlap
             !
             ! variation w.r.t. s
             vspinfi  = FI1(N1)*FI2(N2)*ovlapphi
             vspinr   = FI1(N1)*FI2(N2)*ovlapr
             vspinz   = FI1(N1)*FI2(N2)*ovlapz
             !
             ! variation w.r.t. T
             vttbase = FIR1(N1)*FIR2(N2) +FIZ1(N1)*FIZ2(N2)
             vttfi = ( vttbase -LAM1(N1)*LAM2(N2)*y*y*FI1(N1)*FI2(N2) )*ovlapphi
             vttr  = ( vttbase -LAM1(N1)*LAM2(N2)*y*y*FI1(N1)*FI2(N2) )*ovlapr
             vttz  = ( vttbase +LAM1(N1)*LAM2(N2)*y*y*FI1(N1)*FI2(N2) )*ovlapz
             !
             ! these fields have (e^+iKphi - e^-iKphi) structure, hence the sign
             if(LAM1(n1) + LAM2(n2) +1 .eq. abs(K_qrpa_responce)) then
               vspinr = -vspinr
               vttr  = -vttr
             end if
             !
             ! these fields have (e^+iKphi - e^-iKphi) structure, hence the sign
             if(LAM1(n1) - LAM2(n2) .gt. 0) then
               vspinz = -vspinz
               vjfi = -vjfi
               vttz = -vttz
             end if
             !
             ! matrix elements from conjugate states (simplex -i states)
             vspinrb  = -vspinr  ; vspinfib = -vspinfi ;  vspinzb  = vspinz
             vttrb = -vttr       ; vttfib = -vttfi     ;  vttzb = vttz
             !
             ! zero those that don't contribute in K=0 case
             if(K_qrpa_responce.eq.0) then
               vspinr = zero; vspinz = zero; vspinrb = zero; vspinzb = zero; vjfi = zero;
               vttr   = zero; vttz   = zero; vttrb   = zero; vttzb   = zero;
             end if
             !
             ! correct matrix element terms when using normal HO basis
             If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) then
               vspinfib = vspinfi
             end if
             !
             ! matrix elements with respect of +i states
             hN_c(j) = hN_c(j) +vh*vn_c(ihil)     &                    !! rho_n
                               +vhb*vhbn_c(ihil)  &                    !! tau_n
                               +vjr*vjRN_c(ihil)*iunit2 &              !! jr_n
                               +vjz*vjZN_c(ihil)*iunit2 &              !! jz_n
                               +vjfi*vjFIN_c(ihil) &                   !! jfi_n
                               +vdel*vdn_c(ihil) &                     !! Delta_n
                               +vso*vsn_c(ihil)  &                     !! nabla J
                               +vjrf*vJrfin_c(ihil) &                  !! J_r,fi_n
                               +vjfr*vJfirn_c(ihil) &                  !! J_fi,r_n
                               +vjfz*vJfizn_c(ihil) &                  !! J_fi,z_n
                               +vjzf*vJzfin_c(ihil) &                  !! J_z,fi_n
                               +vjrr*vJrrn_c(ihil) &                   !! J_r,r_n
                               +vjrz*vJrzn_c(ihil) &                   !! J_r,z_n
                               +vjff*vJfifin_c(ihil) &                 !! J_fi,fi_n
                               +vjzr*vJzrn_c(ihil) &                   !! J_z,r_n
                               +vjzz*vJzzn_c(ihil) &                   !! J_z,z_n
                               +vspinr*vSIGRn_c(ihil) &                !! spin_r
                               +vspinfi*vSIGFIn_c(ihil) &              !! spin_fi 
                               +vspinz*vSIGZn_c(ihil) &                !! spin_z
                               +vttr*vTRn_c(ihil) &                    !! T_r
                               +vttfi*vTFIn_c(ihil) &                  !! T_fi
                               +vttz*vTZn_c(ihil)                      !! T_z
             !
             hP_c(j) = hP_c(j) +vh*vp_c(ihil)     &                    !! rho_p
                               +vhb*vhbp_c(ihil)  &                    !! tau_n
                               +vjr*vjRP_c(ihil)*iunit2 &              !! jr_p
                               +vjz*vjZP_c(ihil)*iunit2 &              !! jz_p
                               +vjfi*vjFIP_c(ihil) &                   !! jfi_p
                               +vdel*vdp_c(ihil) &                     !! Delta_p
                               +vso*vsp_c(ihil) &                      !! nabla J
                               +vjrf*vJrfip_c(ihil) &                  !! J_r,fi_p
                               +vjfr*vJfirp_c(ihil) &                  !! J_fi,r_p
                               +vjfz*vJfizp_c(ihil) &                  !! J_fi,z_p
                               +vjzf*vJzfip_c(ihil) &                  !! J_z,fi_p
                               +vjrr*vJrrp_c(ihil) &                   !! J_r,r_p
                               +vjrz*vJrzp_c(ihil) &                   !! J_r,z_p
                               +vjff*vJfifip_c(ihil) &                 !! J_fi,fi_p
                               +vjzr*vJzrp_c(ihil) &                   !! J_z,r_p
                               +vjzz*vJzzp_c(ihil) &                   !! J_z,z_p
                               +vspinr*vSIGRp_c(ihil) &                !! spin_r
                               +vspinfi*vSIGFIp_c(ihil) &              !! spin_fi 
                               +vspinz*vSIGZp_c(ihil) &                !! spin_z
                               +vttr*vTRp_c(ihil) &                    !! T_r
                               +vttfi*vTFIp_c(ihil) &                  !! T_fi
                               +vttz*vTZp_c(ihil)                      !! T_z
             !
             ! matrix elements with respect of conjugate (-i) simplex-y states (note that spin-independent part is the same)
             BhN_c(j)=BhN_c(j) +vh*vn_c(ihil)     &                    !! rho_n
                               +vhb*vhbn_c(ihil)  &                    !! tau_n
                               +vjr*vjRN_c(ihil)*iunit2 &              !! jr_n
                               +vjz*vjZN_c(ihil)*iunit2 &              !! jz_n
                               +vjfi*vjFIN_c(ihil) &                   !! jfi_n
                               +vdel*vdn_c(ihil) &                     !! Delta_n
                               +vsob*vsn_c(ihil)  &                    !! nabla J
                               +vjrfb*vJrfin_c(ihil) &                 !! J_r,fi_n
                               +vjfrb*vJfirn_c(ihil) &                 !! J_fi,r_n
                               +vjfzb*vJfizn_c(ihil) &                 !! J_fi,z_n
                               +vjzfb*vJzfin_c(ihil) &                 !! J_z,fi_n   
                               +vjrrb*vJrrn_c(ihil) &                  !! J_r,r_n
                               +vjrzb*vJrzn_c(ihil) &                  !! J_r,z_n
                               +vjffb*vJfifin_c(ihil) &                !! J_fi,fi_n
                               +vjzrb*vJzrn_c(ihil) &                  !! J_z,r_n
                               +vjzzb*vJzzn_c(ihil) &                  !! J_z,z_n
                               +vspinrb*vSIGRn_c(ihil) &               !! spin_r
                               +vspinfib*vSIGFIn_c(ihil) &             !! spin_fi 
                               +vspinzb*vSIGZn_c(ihil) &               !! spin_z
                               +vttrb*vTRn_c(ihil) &                   !! T_r
                               +vttfib*vTFIn_c(ihil) &                 !! T_fi
                               +vttzb*vTZn_c(ihil)                     !! T_z

             !
             BhP_c(j)=BhP_c(j) +vh*vp_c(ihil)     &                    !! rho_p
                               +vhb*vhbp_c(ihil)  &                    !! tau_n
                               +vjr*vjRP_c(ihil)*iunit2 &              !! jr_p
                               +vjz*vjZP_c(ihil)*iunit2 &              !! jz_p
                               +vjfi*vjFIP_c(ihil) &                   !! jfi_p
                               +vdel*vdp_c(ihil) &                     !! Delta_p
                               +vsob*vsp_c(ihil) &                     !! nabla J
                               +vjrfb*vJrfip_c(ihil) &                 !! J_r,fi_p
                               +vjfrb*vJfirp_c(ihil) &                 !! J_fi,r_p
                               +vjfzb*vJfizp_c(ihil) &                 !! J_fi,z_p
                               +vjzfb*vJzfip_c(ihil) &                 !! J_z,fi_p
                               +vjrrb*vJrrp_c(ihil) &                  !! J_r,r_p
                               +vjrzb*vJrzp_c(ihil) &                  !! J_r,z_p
                               +vjffb*vJfifip_c(ihil) &                !! J_fi,fi_p
                               +vjzrb*vJzrp_c(ihil) &                  !! J_z,r_p
                               +vjzzb*vJzzp_c(ihil) &                  !! J_z,z_p
                               +vspinrb*vSIGRp_c(ihil) &               !! spin_r
                               +vspinfib*vSIGFIp_c(ihil) &             !! spin_fi 
                               +vspinzb*vSIGZp_c(ihil) &               !! spin_z
                               +vttrb*vTRp_c(ihil) &                   !! T_r
                               +vttfib*vTFIp_c(ihil) &                 !! T_fi
                               +vttzb*vTZp_c(ihil)                     !! T_z
             !
             dN_c(j)=dN_c(j)+ vh*davn_c(ihil)
             BdN_c(j)=BdN_c(j)+ vh*dbvn_c(ihil)
             dP_c(j)=dP_c(j)+ vh*davp_c(ihil)
             BdP_c(j)=BdP_c(j)+ vh*dbvp_c(ihil)
             !
             j=j+1  !! go one step forward in the matrix array
             !
          End Do !! Do N2=1,iQRPA_IDI(ib) !! states in the block
        End Do !! Do N1=1,iQRPA_IDF(ib)   !! states in the block
      End Do !! Do ihil=1,nghl
    End Do ! ibq=1,iQRPA_NB
    !
    if(USE_ONLY_F20) then
      hN_c = zero; hP_c = zero; BhN_c = zero; BhP_c = zero; 
      dN_c = zero; dP_c = zero; BdN_c = zero; BdP_c = zero; 
    end if
    !
!WRITE(*,*) "gamdel:"
!call qrpa_printProp("hN_c",hN_c)
!call qrpa_printProp("BhN_c",BhN_c)

    If(Print_Screen.And.IDEBUG.Gt.10) Then 
       Write(*,*) 
       Write(*,*) 'gamdel: maxval(hN)=',Maxval(Real(hN_c)),Maxval(Imag(hN_c))
       Write(*,*) '        maxval(hP)=',Maxval(Real(hP_c)),Maxval(Imag(hP_c))
       Write(*,*) '        maxval(dN)=',Maxval(Real(dN_c)),Maxval(Imag(dN_c))
       Write(*,*) '        maxval(dP)=',Maxval(Real(dP_c)),Maxval(Imag(dP_c))
       Write(*,*) '        maxval(BdN)=',Maxval(Real(BdN_c)),Maxval(Imag(BdN_c))
       Write(*,*) '        maxval(BdP)=',Maxval(Real(BdP_c)),Maxval(Imag(BdP_c))
       Write(*,*) 
       Write(*,*) '        minval(hN)=',Minval(Real(hN_c)),Minval(Imag(hN_c))
       Write(*,*) '        minval(hP)=',Minval(Real(hP_c)),Minval(Imag(hP_c))
       Write(*,*) '        minval(dN)=',Minval(Real(dN_c)),Minval(Imag(dN_c))
       Write(*,*) '        minval(dP)=',Minval(Real(dP_c)),Minval(Imag(dP_c))
       Write(*,*) '        minval(BdN)=',Minval(Real(BdN_c)),Minval(Imag(BdN_c))
       Write(*,*) '        minval(BdP)=',Minval(Real(BdP_c)),Minval(Imag(BdP_c))
       Write(*,*) 
    End If
    !
    Call get_CPU_time('qrpa_gamdel',1)
    !
  End Subroutine qrpa_gamdel
  !=================================================================================================================================
  !  
  !=================================================================================================================================
  Subroutine qrpa_f_gamdel
    !---------------------------------------------------------------------
    ! Fh_ab  in configurational space
    ! Fh     in qp space, both referenced in full
    !
    ! note: ovlapsy may not work properly, if operator is
    !       a product of sigma_y and phi dependent quantity
    !---------------------------------------------------------------------
    Implicit None
    Integer(ipr) :: i,j,ihil,im,nsa,JA,N1,N2,ibq,ovlap
    Real(pr)     :: vfn,vfp,zcoord,rcoord,y
    complex(pr)  :: vh,vly1,vly2, vspiny, vspinyb, vspinz, vlz, vspinzb, vlzb
    complex(pr)  :: vpz,vpx1,vpx2,vcoz,vcox
    Real(pr), allocatable :: FI1(:),FI2(:),FIR1(:),FIR2(:),FIZ1(:),FIZ2(:)
    Integer (ipr), allocatable :: LAM1(:),LAM2(:),SS1(:),SS2(:),NZ1(:),NZ2(:),NR1(:),NR2(:)
    Real(pr)     :: vfna(nghl),vfpa(nghl),vfnb(nghl),vfpb(nghl)
    complex(pr)  :: iunit
    complex(pr)  :: ovlaplz,ovlaply    !! overlap when l_z or l_y operator
    complex(pr)  :: ovlapsz,ovlapsy    !! overlap when sigma_z or sigma_y operator
    complex(pr)  :: ovlappz,ovlappx    !! overlap for momentum P_z and P_x operator
    complex(pr)  :: ovlapcoz,ovlapcox  !! overlap for coordinate Q_z operator 
    !
    Parameter(iunit = Cmplx(0.0_pr,1.0_pr,kind=pr))
    !-----------------------------------------------------------------
    ! F_ab in configurational space: FhNab_c,FhPab_c,BFhNab_c,BFhPab_c
    !-----------------------------------------------------------------
    FhNab_c=zero; FhPab_c=zero; BFhNab_c=zero; BFhPab_c=zero;
    ! fields
    select case(QRPA_Operator_type)
    case(1)
      vfna = qrpa_electric_field(0,zero,fl,fh)
      vfpa = qrpa_electric_field(1,zero,fl,fh) 
    case(2)
      vfna = qrpa_spin_field(0,zero,fl,fh)
      vfpa = qrpa_spin_field(1,zero,fl,fh) 
    case(3)
      vfna = qrpa_spin_field(0,zero,fl,fh)
      vfpa = qrpa_spin_field(1,zero,fl,fh) 
      vfnb = qrpa_l_field(0,zero,fl,fh)
      vfpb = qrpa_l_field(1,zero,fl,fh) 
    case(4)
    case(5)
    case(6)
    case(7)
    case default
      STOP "Unknown operator type"
    end select
    !
    if(allocated(FI1)) Deallocate(FI1,FI2,FIR1,FIR2,FIZ1,FIZ2)
    allocate(FI1(nqx),FI2(nqx),FIR1(nqx),FIR2(nqx),FIZ1(nqx),FIZ2(nqx))
    !
    if(allocated(lam1)) Deallocate(lam1,ss1,lam2,ss2,nz1,nz2,nr1,nr2)
    allocate(lam1(nqx),lam2(nqx),ss1(nqx),ss2(nqx),nz1(nqx),nz2(nqx),nr1(nqx),nr2(nqx))
    !
    Do ibq=1,iQRPA_NB
       Do ihil=1,nghl
         !
         ! run trough the left basis vectors and fill the array
         Do N1=1,iQRPA_IDF(ibq)
           IM=ia(iQRPA_BLF(ibq)); JA=IM+N1; NSA=NS(JA);
           FI1(N1)  = QHLA_opt(JA,ihil);
           FIR1(N1) = FI1R_opt(JA,ihil)
           FIZ1(N1) = FI1Z_opt(JA,ihil)
           SS1(N1) = NSA
           NZ1(N1) = NZ(JA)
           NR1(N1) = NR(JA)
           if(NSA.gt.0) then
             LAM1(N1) = (iQRPA_OMF(ibq) -1)/2
           else
             LAM1(N1) = -(iQRPA_OMF(ibq) +1)/2
           end if
           If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) LAM1(N1) = abs(LAM1(N1))  !! normal HO, positive omega states have positive lambda
         End Do
         !
         ! run trough the right basis vectors and fill the array
         Do N1=1,iQRPA_IDI(ibq)
           IM=ia(iQRPA_BLI(ibq)); JA=IM+N1; NSA=NS(JA);
           FI2(N1)  = QHLA_opt(JA,ihil);
           FIR2(N1) = FI1R_opt(JA,ihil)
           FIZ2(N1) = FI1Z_opt(JA,ihil)
           SS2(N1) = NSA
           NZ2(N1) = NZ(JA)
           NR2(N1) = NR(JA)
           if(NSA.gt.0) then
             LAM2(N1) = (iQRPA_OMI(ibq) -1)/2
           else
             LAM2(N1) = -(iQRPA_OMI(ibq) +1)/2
           end if
           If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) LAM2(N1) = abs(LAM2(N1))  !! normal HO, positive omega states have positive lambda
         End Do
         !
         vfn=vfna(ihil); vfp=vfpa(ihil);
         j=iQRPA_DS(ibq)
         zcoord=fh(ihil)
         rcoord=fl(ihil)
         y=y_opt(ihil)   !! y=1/r
         !
         Do N2=1,iQRPA_IDI(ibq)
           Do N1=1,iQRPA_IDF(ibq)
              !
              ovlap = 0; ovlapsz = zero; ovlapsy = zero; ovlaply = zero; ovlaplz = zero; ovlapsz = zero;
              ovlappz = zero; ovlappx = zero ; ovlapcoz = zero; ovlapcox = zero
              if( abs(LAM1(n1) + LAM2(n2)) .eq. abs(K_qrpa_responce) ) ovlap = 1
              if( abs(LAM1(n1) - LAM2(n2)) .eq. abs(K_qrpa_responce) ) ovlap = 1
              !
              if( (LAM1(n1) .eq. -LAM2(n2)) .and. (abs(K_qrpa_responce).eq.1) ) ovlapsy = -1.0_pr; !! this is also used for M1. s_y is K=1 operator
              !
              if(abs(LAM1(n1)-LAM2(n2)).eq.1) ovlaply = 1.0_pr
              if(abs(LAM1(n1)+LAM2(n2)).eq.1) ovlaply = 1.0_pr
              if(abs(K_qrpa_responce).ne.1)   ovlaply = zero;  !! l_y is K=1 operator

              ! full conditions from the angular integrals of exponentials for the LINEAR MOMENTUM and COORDINATE operator
              if( (abs(LAM1(n1) - LAM2(n2)) .eq. abs(K_qrpa_responce)) .and. (abs(K_qrpa_responce).eq.0) )   ovlappz  = 1.0_pr;              ! P_z
              if( (abs(LAM1(n1) - LAM2(n2)) .eq. abs(K_qrpa_responce)) .and. (abs(K_qrpa_responce).eq.1) )   ovlappx  = 1.0_pr;              ! P_x
              if( (abs(LAM1(n1) - LAM2(n2)) .eq. abs(K_qrpa_responce)) .and. (abs(K_qrpa_responce).eq.0) )   ovlapcoz = 1.0_pr;              ! Q_z
              if( (abs(LAM1(n1) - LAM2(n2)) .eq. abs(K_qrpa_responce)) .and. (abs(K_qrpa_responce).eq.1) )   ovlapcox = 1.0_pr;              ! Q_x
              !
              !
              !! hormal HO and s_z seems not to work properly... fix this part in future version
              If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) then
                 if( LAM1(n1).eq.LAM2(n2) ) then ; ovlaplz = 1.0_pr ; ovlapsz = 1.0 ; end if 
                 if( SS1(n1).lt.0 ) ovlapsz = ovlapsz * (-1.0_pr)
              end if

              vh=FI1(N1)*FI2(N2)*ovlap
              vspiny=FI1(N1)*FI2(N2)*half*ovlapsy  ! s = (1/2) sigma
              vspinz=FI1(N1)*FI2(N2)*half*ovlapsz  ! s = (1/2) sigma
              !
              ! terms for l_y
              vly1 =-0.25_pr*iunit*rcoord*(FI1(N1)*FIZ2(N2)-FIZ1(N1)*FI2(N2))*ovlaply &
                    +0.25_pr*iunit*zcoord*(FI1(N1)*FIR2(N2)-FIR1(N1)*FI2(N2))*ovlaply

              vly2 = 0.25_pr*iunit*( (LAM2(N2)+LAM1(N1))*zcoord*y*FI1(N1)*FI2(N2))*ovlaply
              !
              vlz = LAM1(N1)*FI1(N1)*FI2(N2)*ovlaplz
              !
              vpz  = 0.50_pr*iunit*(FI1(N1)*FIZ2(N2)-FIZ1(N1)*FI2(N2))*ovlappz
              vpx1 = 0.25_pr*iunit*(FI1(N1)*FIR2(N2)-FIR1(N1)*FI2(N2))*ovlappx
              vpx2 = 0.25_pr*iunit*FI1(N1)*y*FI2(N2)*(LAM1(N1)+LAM2(N2))*ovlappx
              !
              ! CONDITION FOR THE SIGN
              if(LAM1(N1) - LAM2(N2) .gt. 0) then
                vpx2 = -vpx2                                                     !krisa! vpx2 changes sign depending on the chosen condition for lambda and K numbers
              end if
              !
              !Q_z and Q_x momentum operator (element of Q_y vanishes)
              vcoz = (1/Dble(npr(1)+npr(2)))*FI1(N1)*zcoord*FI2(N2)*ovlapcoz
              vcox = (1/(2*Dble(npr(1)+npr(2))))*FI1(N1)*rcoord*FI2(N2)*ovlapcox
              !
              ! this term has (e^+iKphi - e^-iKphi) structure, hence the sign
              if(LAM1(n1) - LAM2(n2) .gt. 0) then
                 vly2 = -vly2
              end if
              !
              ! simplex -i states
              vspinyb = -vspiny 
              !
              ! negative Omega states
              vspinzb = -vspinz; vlzb = -vlz
              !
              select case(QRPA_Operator_type)
              case(1)
                ! electric operator
                ! simplex +i states and simplex -i state are same for electric operator
                FhNab_c(j) = FhNab_c(j)  + vh *vfn 
                FhPab_c(j) = FhPab_c(j)  + vh *vfp 
                BFhNab_c(j)= BFhNab_c(j) + vh *vfn 
                BFhPab_c(j)= BFhPab_c(j) + vh *vfp 
              case(2)
                ! spin operator
                ! simplex +i states and simplex -i state are different
                ! note: operator s_11 = i*s_y
                 FhNab_c(j)  = FhNab_c(j)  + (iunit*vspiny +vspinz) *vfna(ihil)
                 FhPab_c(j)  = FhPab_c(j)  + (iunit*vspiny +vspinz) *vfpa(ihil)
                 BFhNab_c(j) = BFhNab_c(j) + (iunit*vspinyb+vspinzb)*vfna(ihil)
                 BFhPab_c(j) = BFhPab_c(j) + (iunit*vspinyb+vspinzb)*vfpa(ihil)
              case(3)
                ! M1 operator
                ! simplex +i states and simplex -i state are different for spin part
                ! note: operator l_11 = i*l_y
                FhNab_c(j)  = FhNab_c(j)  + (iunit*vspiny +vspinz) *vfna(ihil) + (iunit*(vly1+vly2) +vlz)*vfnb(ihil)
                FhPab_c(j)  = FhPab_c(j)  + (iunit*vspiny +vspinz) *vfpa(ihil) + (iunit*(vly1+vly2) +vlz)*vfpb(ihil)
                BFhNab_c(j) = BFhNab_c(j) + (iunit*vspinyb+vspinzb)*vfna(ihil) + (iunit*(vly1+vly2) +vlzb)*vfnb(ihil)
                BFhPab_c(j) = BFhPab_c(j) + (iunit*vspinyb+vspinzb)*vfpa(ihil) + (iunit*(vly1+vly2) +vlzb)*vfpb(ihil)
              case(4)
                ! LINEAR MOMENTUM (P10,P11) OPERATOR
                ! Since P10 and P11 are not present simultaneusly, they can be added here
                FhNab_c(j)  = FhNab_c(j)  + vpz - (vpx1+vpx2) 
                FhPab_c(j)  = FhPab_c(j)  + vpz - (vpx1+vpx2)
                BFhNab_c(j) = BFhNab_c(j) + vpz - (vpx1+vpx2)
                BFhPab_c(j) = BFhPab_c(j) + vpz - (vpx1+vpx2)
              case(5)
                ! COORDINATE (Q10,Q11) OPERATOR ------------------------------------------------
                FhNab_c(j)  = FhNab_c(j)  + vcoz - vcox
                FhPab_c(j)  = FhPab_c(j)  + vcoz - vcox
                BFhNab_c(j) = BFhNab_c(j) + vcoz - vcox
                BFhPab_c(j) = BFhPab_c(j) + vcoz - vcox
              case(6)
                ! TOTAL ANGULAR MOMENTUM (J11) OPERATOR ----------------------------------------
                ! simplex +i states and simplex -i state are different for spin part !
                ! note: (tilde)L_11 = 1/Sqrt[2] * (L_(11)+L_(1-1)) = -i*L_y is the K = 1 operator (plus in the bracket is because L (axial vector) is parity=+1 operator)
                FhNab_c(j)  = FhNab_c(j)  + iunit*vspiny  + iunit*(vly1+vly2)
                FhPab_c(j)  = FhPab_c(j)  + iunit*vspiny  + iunit*(vly1+vly2)
                BFhNab_c(j) = BFhNab_c(j) + iunit*vspinyb + iunit*(vly1+vly2)
                BFhPab_c(j) = BFhPab_c(j) + iunit*vspinyb + iunit*(vly1+vly2)
              case(7)
                ! TOTAL ANGULAR MOMENTUM (Jy) OPERATOR ----------------------------------------
                ! simplex +i states and simplex -i state are different for spin part !
                ! note: we use simply J_y = S_y + L_y  (J_y is K = 1 operator)
                FhNab_c(j)  = FhNab_c(j)  + vspiny  + (vly1+vly2)
                FhPab_c(j)  = FhPab_c(j)  + vspiny  + (vly1+vly2)
                BFhNab_c(j) = BFhNab_c(j) + vspinyb + (vly1+vly2)
                BFhPab_c(j) = BFhPab_c(j) + vspinyb + (vly1+vly2)
              case default
                STOP "Unknown operator type"
              end select
              !
              j=j+1
              !
           End Do   ! N1=1,iQRPA_IDF(ibq)
         End Do     ! N2=1,iQRPA_IDI(ibq)
      End Do  !! ihil
    End Do !! ibq
    !
    !----------------------------------------------
    ! F20 and F02 in qp space: FhN_c, FhP_c, BFhN_c, BFhP_c
    ! F11 calculated if required by EFA
    !----------------------------------------------
    Call qrpa_F20(VqpN,UqpN,FhNab_c,BFhNab_c,FhN_c,BFhN_c)
    Call qrpa_F20(VqpP,UqpP,FhPab_c,BFhPab_c,FhP_c,BFhP_c)
    if(ENABLE_N_EFA) Call qrpa_F11(VqpN,UqpN,FhNab_c,BFhNab_c,F11N_c,BF11N_c)
    if(ENABLE_P_EFA) Call qrpa_F11(VqpP,UqpP,FhPab_c,BFhPab_c,F11P_c,BF11P_c)
!
!if(ENABLE_N_EFA) then
!write(*,*) "f11"
!call qrpa_printprop("F11",F11N_c)
!call qrpa_printprop("F11b",BF11N_c)
!end if

if(.false.) then
call qrpa_printprop("FhNab_c",FhNab_c)
!do ibq=99,222; write(*,'(F10.6,SP,F10.6,"i  ")',advance='NO') FhNab_c(ibq) ; end do;
call qrpa_printprop("BFhNab_c",BFhNab_c)
!do ibq=99,222; write(*,'(F10.6,SP,F10.6,"i  ")',advance='NO') BFhNab_c(ibq) ; end do;

qrpa_aux = FhNab_c-BFhNab_c ; call qrpa_printProp("FhNab_c-BFhNab_c",qrpa_aux)
!do ibq=99,222; write(*,'(F10.6,SP,F10.6,"i  ")',advance='NO') qrpa_aux(ibq) ; end do;

qrpa_aux = FhNab_c+BFhNab_c ; call qrpa_printProp("FhNab_c+BFhNab_c",qrpa_aux)
!do ibq=99,222; write(*,'(F10.6,SP,F10.6,"i  ")',advance='NO') qrpa_aux(ibq) ; end do;


!qrpa_aux = F11N_c -BF11N_c ; call qrpa_printProp("F11-BF11",qrpa_aux)
!qrpa_aux = F11N_c +BF11N_c ; call qrpa_printProp("F11+BF11",qrpa_aux)
!stop
end if

if(.false.) then
WRITE(*,*) "Eqp:" 
WRITE(*,'(10f10.6)') REqpN
WRITE(*,*) "fn:"
WRITE(*,'(10f10.6)') REAL(FhNab_c)
WRITE(*,*) "Re un:"
WRITE(*,'(10f10.6)') REAL(UqpN)
WRITE(*,*) "Im un:"
WRITE(*,'(10f10.6)') Imag(UqpN)
WRITE(*,*) "f11n"
WRITE(*,'(10f10.6)') Real(F11N_c)
WRITE(*,*) "f11nb"
WRITE(*,'(10f10.6)') Real(bF11N_c)

STOP
end if
    !
    Deallocate(FI1,FI2,FIR1,FIR2,FIZ1,FIZ2,LAM1,LAM2,ss1,ss2,nz1,nz2,nr1,nr2)
    !
  End Subroutine qrpa_f_gamdel
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpaINI_CUT_ALLOCATE
    !-----------------------------------------------------------------------------------------------------
    ! Alocate qrpa arrays and impose cut-off conditions
    !-----------------------------------------------------------------------------------------------------
    Implicit None
    Integer(ipr) :: i,j,k,iblk,iline,ip,iq,ib,ibb,nd,ndb,iw,jprot,jneut,n1,n2,jtr,nsa,lam
    Real(pr)     :: EqpCUT,PqpCUT,pn,eqpe,ela,enb
    complex(pr)  :: phase
    logical :: isactive
    !--------------------------------------------------------------------------
    ! test HFBTHO U,V for orthonormality
    !--------------------------------------------------------------------------
    Call qrpa_test_HFB_UV(RVqpN,RUqpN,1)  
    Call qrpa_test_HFB_UV(RVqpP,RUqpP,2)
    !
    ! Check here if EFA activated
    !
    !--------------------------------------------------------------------------
    ! QRPA: Allocate Arrays 
    !--------------------------------------------------------------------------
    ! Arrays depending on gauss points for time-even densities and fields
    If(Allocated(vn_c)) Deallocate(aka_c,bka_c,vhbn_c,vn_c,vrn_c,vzn_c,vdn_c,vsn_c,davn_c,dbvn_c &
         ,vhbp_c,vp_c,vrp_c,vzp_c,vdp_c,vsp_c,davp_c,dbvp_c  &
         ,ro_c,tau_c,dro_c,dj_c,NABLAR_c,NABLAZ_c,cou_c,vc_c )
    Allocate(aka_c(nghl,2),bka_c(nghl,2),ro_c(nghl,2),tau_c(nghl,2),dro_c(nghl,2),dj_c(nghl,2) & 
         ,vhbn_c(nghl),vn_c(nghl),vrn_c(nghl),vzn_c(nghl),vdn_c(nghl),vsn_c(nghl),davn_c(nghl),dbvn_c(nghl) &
         ,vhbp_c(nghl),vp_c(nghl),vrp_c(nghl),vzp_c(nghl),vdp_c(nghl),vsp_c(nghl),davp_c(nghl),dbvp_c(nghl)  &         
         ,NABLAR_c(nghl,2),NABLAZ_c(nghl,2),cou_c(nghl),vc_c(nghl,nghl) )
    ! Arrays depending on gauss points for time-even densities and fields: tensor terms
    if(Allocated(JRR_c)) Deallocate(JRR_c,JRFI_c,JRZ_c,JFIR_c,JFIFI_c,JFIZ_c,JZR_c,JZFI_c,JZZ_c &
         ,vJRRn_c,vJRFIn_c,vJRZn_c,vJFIRn_c,vJFIFIn_c,vJFIZn_c,vJZRn_c,vJZFIn_c,vJZZn_c &
         ,vJRRp_c,vJRFIp_c,vJRZp_c,vJFIRp_c,vJFIFIp_c,vJFIZp_c,vJZRp_c,vJZFIp_c,vJZZp_c )
    Allocate(JRR_c(nghl,2),JRFI_c(nghl,2),JRZ_c(nghl,2),JFIR_c(nghl,2),JFIFI_c(nghl,2),JFIZ_c(nghl,2) &
         ,JZR_c(nghl,2),JZFI_c(nghl,2),JZZ_c(nghl,2) &
         ,vJRRn_c(nghl),vJRFIn_c(nghl),vJRZn_c(nghl),vJFIRn_c(nghl),vJFIFIn_c(nghl) &
         ,vJFIZn_c(nghl),vJZRn_c(nghl),vJZFIn_c(nghl),vJZZn_c(nghl) &
         ,vJRRp_c(nghl),vJRFIp_c(nghl),vJRZp_c(nghl),vJFIRp_c(nghl),vJFIFIp_c(nghl) &
         ,vJFIZp_c(nghl),vJZRp_c(nghl),vJZFIp_c(nghl),vJZZp_c(nghl)  )
    !
    ! Arrays depending on gauss points for time-odd densities and fields
    If(Allocated(sor_c)) Deallocate(sor_c,soz_c,sofi_c,dsor_c,dsoz_c,dsofi_c, &
         rjr_c,rjz_c,rjfi_c,tsor_c,tsoz_c,tsofi_c,curjr_c,curjz_c,curjfi_c,cursr_c,cursz_c,cursfi_c, &
         vSIGRn_c,vSIGZn_c,vSIGFIn_c,vSIGRp_c,vSIGZp_c,vSIGFIp_c, &
         vSIDRn_c,vSIDZn_c,vSIDFIn_c,vSIDRp_c,vSIDZp_c,vSIDFIp_c, &
         vTRn_c,vTZn_c,vTFIn_c,vTRp_c,vTZp_c,vTFIp_c)

    Allocate(sor_c(nghl,2),soz_c(nghl,2),sofi_c(nghl,2),dsor_c(nghl,2), &
         dsoz_c(nghl,2),dsofi_c(nghl,2),rjr_c(nghl,2),rjz_c(nghl,2), &
         rjfi_c(nghl,2),tsor_c(nghl,2),tsoz_c(nghl,2),tsofi_c(nghl,2), &
         curjr_c(nghl,2),curjz_c(nghl,2),curjfi_c(nghl,2),cursr_c(nghl,2),cursz_c(nghl,2),cursfi_c(nghl,2), &
         vSIGRn_c(nghl),vSIGZn_c(nghl),vSIGFIn_c(nghl),vSIGRp_c(nghl),vSIGZp_c(nghl),vSIGFIp_c(nghl), &
         vSIDRn_c(nghl),vSIDZn_c(nghl),vSIDFIn_c(nghl),vSIDRp_c(nghl),vSIDZp_c(nghl),vSIDFIp_c(nghl), &
         vjRn_c(nghl),vjZn_c(nghl),vjFIn_c(nghl),vjRp_c(nghl),vjZp_c(nghl),vjFIp_c(nghl),             &
         vTRn_c(nghl),vTZn_c(nghl),vTFIn_c(nghl), &
         vTRp_c(nghl),vTZp_c(nghl),vTFIp_c(nghl) )

    !
    ! for density dependent case, use hfb densities
    If(Allocated(ro_c_hfb)) Deallocate(ro_c_hfb,tau_c_hfb,kap_c_hfb)
    Allocate(ro_c_hfb(nghl,2),tau_c_hfb(nghl,2),kap_c_hfb(nghl,2))
    ro_c_hfb = ro ; tau_c_hfb = tau ; kap_c_hfb = aka ;
    ! 
    ! use less cryptic variable names
    iHFB_LEN = nuv
    iHFB_NQP = nqp
    iHFB_NB = NB
    !
    ! Arrays depending on HFB block structure
    If(Allocated(iHFB_ID)) Deallocate(iHFB_ID,iHFB_DS,iHFB_LP,iHFB_OM,iHFB_QPBS,iHFB_PR)
    Allocate(iHFB_ID(iHFB_NB),iHFB_DS(iHFB_NB),iHFB_LP(iHFB_NB),iHFB_OM(iHFB_NB), &
             iHFB_QPBS(iHFB_NB),iHFB_PR(iHFB_NB))
    If(Allocated(VqpN)) Deallocate(VqpN,VqpP,UqpN,UqpP)
    Allocate(VqpN(iHFB_LEN),VqpP(iHFB_LEN),UqpN(iHFB_LEN),UqpP(iHFB_LEN)) 
    If(Allocated(hfb_UA)) Deallocate(hfb_UA,hfb_VA)
    Allocate(hfb_UA(iHFB_LEN),hfb_VA(iHFB_LEN))
    If(Allocated(ipwi_hfb)) Deallocate(ipwi_hfb)
    Allocate(ipwi_hfb(iHFB_NB,2))
    !
    ! calculate block sizes and structure for HFB
    i=1 ; j=1
    Do ib=1,NB
       ND=ID(ib)
       iHFB_DS(ib)=i
       iHFB_ID(ib)=ND
       iHFB_QPBS(ib)=j
       ipwi_hfb(ib,1)=kd(ib,1)
       ipwi_hfb(ib,2)=kd(ib,2)
       If(Parity) Then
         iHFB_LP(ib)=(ib+1)/2
         iHFB_OM(ib)=2*nl(j)+ns(j)
         iHFB_PR(ib) = 3-2*npar(j)
       else
         iHFB_LP(ib)=ib
         iHFB_OM(ib)=2*ib-1
         iHFB_PR(ib) = 0  ! parity not conserved         
       end if
       i = i+(ND*ND)
       j = j+ND
    End Do
    !
    ! parity of the operator 
    select case(QRPA_Operator_type)
    case(1) !! electric
      if( mod(abs(L_qrpa_responce),2).eq.0 ) then
        P_qrpa_response = +1
      else
        P_qrpa_response = -1
      end if
    case(2) !! spin
      if( abs(L_qrpa_responce).ne.1 ) STOP "Wrong L for spin operator"
      P_qrpa_response = +1
    case(3) !! magnetic
      if( abs(L_qrpa_responce).ne.1 ) STOP "Only L=1 case supported for magnetic operator so far..."
      P_qrpa_response = +1
    case(4) !! linear momentum
      if( abs(L_qrpa_responce).ne.1 ) STOP "Wrong L for linear momentum opr."
      P_qrpa_response = -1
    case(5) !! coordinate
      if( abs(L_qrpa_responce).ne.1 ) STOP "Wrong L for coordinate opr."
      P_qrpa_response = -1
    case(6) !! J_11
      if( abs(L_qrpa_responce).ne.1 ) STOP "Wrong L for J_11 opr."
      P_qrpa_response = +1
    case(7) !! J_11
      if( abs(L_qrpa_responce).ne.1 ) STOP "Wrong L for J_y opr."
      P_qrpa_response = +1
    case default
    end select

    If(Parity.and.(P_qrpa_response.eq.-1)) STOP "Parity conserved mode incompatible with requested operator!"  !! check if this condition can be relaxed later
    !
    ! calculate block sizes and structure for QRPA
    ! first the number of QRPA blocks
    iQRPA_NB = 0
    Do ib=1,iHFB_NB
      Do ibb=1,iHFB_NB
        isactive = .false.
        ! check Omega quantum number
        If( iHFB_OM(ib)+iHFB_OM(ibb)+2*K_qrpa_responce .eq. 0 ) isactive = .true.
        If( iHFB_OM(ib)-iHFB_OM(ibb)+2*K_qrpa_responce .eq. 0 ) isactive = .true.
        If( iHFB_OM(ib)+iHFB_OM(ibb)-2*K_qrpa_responce .eq. 0 ) isactive = .true.
        If( iHFB_OM(ib)-iHFB_OM(ibb)-2*K_qrpa_responce .eq. 0 ) isactive = .true.
        ! check parity
        if(parity .and. ((iHFB_PR(ib)*iHFB_PR(ibb)).ne.P_qrpa_response)) isactive = .false.
        ! add block
        if(isactive) iQRPA_NB = iQRPA_NB+1
      End Do
    End Do

    If(iQRPA_NB .eq. 0) STOP "iQRPA_NB .eq. 0"
    If(Allocated(iQRPA_IDI)) Deallocate(iQRPA_IDI,iQRPA_IDF,iQRPA_LPI,iQRPA_LPF,iQRPA_DS, &
                                       iQRPA_OMI,iQRPA_OMF,iQRPA_BLI,iQRPA_BLF,iQRPA_BHI, &
                                       iQRPA_BHF,iQRPA_TRANSPOSE)
    Allocate(iQRPA_IDI(iQRPA_NB),iQRPA_IDF(iQRPA_NB),iQRPA_LPI(iQRPA_NB),iQRPA_LPF(iQRPA_NB), &
             iQRPA_DS(iQRPA_NB),iQRPA_OMI(iQRPA_NB),iQRPA_OMF(iQRPA_NB),iQRPA_BLI(iQRPA_NB),  &
             iQRPA_BLF(iQRPA_NB),iQRPA_BHI(iQRPA_NB),iQRPA_BHF(iQRPA_NB),iQRPA_TRANSPOSE(iQRPA_NB))
    !
    ! QRPA block sizes etc. 
    i=1 ; j=1; iQRPA_LEN=0;
    Do ib=1,iHFB_NB
      Do ibb=1,iHFB_NB
        isactive = .false.
        ! check Omega quantum number
        If( iHFB_OM(ib)+iHFB_OM(ibb)+2*K_qrpa_responce .eq. 0 ) isactive = .true.
        If( iHFB_OM(ib)-iHFB_OM(ibb)+2*K_qrpa_responce .eq. 0 ) isactive = .true.
        If( iHFB_OM(ib)+iHFB_OM(ibb)-2*K_qrpa_responce .eq. 0 ) isactive = .true.
        If( iHFB_OM(ib)-iHFB_OM(ibb)-2*K_qrpa_responce .eq. 0 ) isactive = .true.
        ! check parity
        if(parity .and. ((iHFB_PR(ib)*iHFB_PR(ibb)).ne.P_qrpa_response)) isactive = .false.
        If(isactive) then
          iQRPA_IDI(i) = iHFB_ID(ib)      ;  iQRPA_IDF(i) = iHFB_ID(ibb)
          iQRPA_OMI(i) = iHFB_OM(ib)      ;  iQRPA_OMF(i) = iHFB_OM(ibb)
          iQRPA_LPI(i) = iHFB_LP(ib)      ;  iQRPA_LPF(i) = iHFB_LP(ibb)
          iQRPA_BLI(i) = ib               ;  iQRPA_BLF(i) = ibb
          iQRPA_BHI(i) = ib               ;  iQRPA_BHF(i) = ibb
          iQRPA_DS(i) = j
          iQRPA_LEN = iQRPA_LEN +iHFB_ID(ib)*iHFB_ID(ibb)
          j=j+iHFB_ID(ib)*iHFB_ID(ibb)
          i=i+1
        end if
      End Do
    End Do
    !
    ! figure out blocks corresponding transpose
    if((K_qrpa_responce .eq. 0).and.(P_qrpa_response.eq.+1)) then
      Do ib = 1,iQRPA_NB
        iQRPA_TRANSPOSE(ib)=ib
      end do
    else
      Do ib = 1,iQRPA_NB
        Do ibb = 1,iQRPA_NB
          if( (iQRPA_BHI(ib) .eq. iQRPA_BHF(ibb)) .and. (iQRPA_BHF(ib) .eq. iQRPA_BHI(ibb)) ) iQRPA_TRANSPOSE(ib)=ibb
        end do
      end do
    end if
    !
    ! print block structure
    if(.false.) then 
      Write(*,*) "iHFB_NB,iQRPA_NB,nuv,iHFB_LEN,iQRPA_LEN",iHFB_NB,iQRPA_NB,nuv,iHFB_LEN,iQRPA_LEN
      WRITE(*,*) "HFB blocks"
      WRITE(*,*) "iHFB_ID(ib),iHFB_DS(ib),iHFB_OM(ib),iHFB_LP(ib),iHFB_QPBS(ib),Parity"
      Do ib=1,iHFB_NB
        WRITE(*,*) iHFB_ID(ib),iHFB_DS(ib),iHFB_OM(ib),iHFB_LP(ib),iHFB_QPBS(ib),iHFB_PR(ib)
      End do
      WRITE(*,*) " "
      WRITE(*,*) "Operator parity:", P_qrpa_response
      WRITE(*,*) "QRPA blocks"
      WRITE(*,*) "ib,iQRPA_IDI(ib),iQRPA_IDF(ib),iQRPA_DS(ib),iQRPA_OMI(ib),iQRPA_OMF(ib),iQRPA_LPI(ib),iQRPA_LPF(ib),iQRPA_BLI(ib),iQRPA_BLF(ib),transpose,iQRPA_BHI(ib),iQRPA_BHF(ib)"
      Do ib=1,iQRPA_NB
        WRITE(*,*) ib,iQRPA_IDI(ib),iQRPA_IDF(ib),iQRPA_DS(ib),iQRPA_OMI(ib),iQRPA_OMF(ib),iQRPA_LPI(ib),iQRPA_LPF(ib),iQRPA_BLI(ib), &
                 iQRPA_BLF(ib),iQRPA_TRANSPOSE(ib),iQRPA_BHI(ib),iQRPA_BHF(ib)
      End do
      WRITE(*,*) "IA",IA
    end if
    !
    ! Arrays depending on QRPA block structure
    If(Allocated(XN_c)) Deallocate(XN_c,YN_c,XP_c,YP_c)
    Allocate(XN_c(iQRPA_LEN),YN_c(iQRPA_LEN),XP_c(iQRPA_LEN),YP_c(iQRPA_LEN))         
    If(Allocated(hN_c)) Deallocate(hN_c,hP_c,dN_c,dP_c)
    Allocate(hN_c(iQRPA_LEN),hP_c(iQRPA_LEN),dN_c(iQRPA_LEN),dP_c(iQRPA_LEN)) 
    If(Allocated(BhN_c)) Deallocate(BhN_c,BhP_c,BdN_c,BdP_c)
    Allocate(BhN_c(iQRPA_LEN),BhP_c(iQRPA_LEN),BdN_c(iQRPA_LEN),BdP_c(iQRPA_LEN)) 
    If(Allocated(FhN_c)) Deallocate(FhN_c,FhP_c,BFhN_c,BFhP_c,FhNab_c,FhPab_c,BFhNab_c,BFhPab_c)
    Allocate(FhN_c(iQRPA_LEN),FhP_c(iQRPA_LEN),BFhN_c(iQRPA_LEN),BFhP_c(iQRPA_LEN), &
             FhNab_c(iQRPA_LEN),FhPab_c(iQRPA_LEN),BFhNab_c(iQRPA_LEN),BFhPab_c(iQRPA_LEN)) 
    If(Allocated(GreenNplus)) Deallocate(GreenNplus,GreenPplus,GreenNminus,GreenPminus)
    Allocate(GreenNplus(iQRPA_LEN),GreenPplus(iQRPA_LEN),GreenNminus(iQRPA_LEN),GreenPminus(iQRPA_LEN))  
    If(Allocated(qrpa_AUX)) Deallocate(qrpa_AUX,qrpa_BUX,qrpa_AUX20,qrpa_BUX20)
    Allocate(qrpa_AUX(iQRPA_LEN),qrpa_BUX(iQRPA_LEN),qrpa_AUX20(iQRPA_LEN),qrpa_BUX20(iQRPA_LEN))
    If(Allocated(rhoetadm)) Deallocate(rhoetadm,kapetadm,kapetabdm)
    Allocate(rhoetadm(iQRPA_LEN),kapetadm(iQRPA_LEN),kapetabdm(iQRPA_LEN))
    If(Allocated(JyNab_c)) Deallocate(JyNab_c,JyPab_c,BJyNab_c,BJyPab_c,JyN_c,JyP_c,BJyN_c,BJyP_c)
    Allocate(JyNab_c(iQRPA_LEN),JyPab_c(iQRPA_LEN),BJyNab_c(iQRPA_LEN),BJyPab_c(iQRPA_LEN), &
             JyN_c(iQRPA_LEN),JyP_c(iQRPA_LEN),BJyN_c(iQRPA_LEN),BJyP_c(iQRPA_LEN))
    !
    ! Check blocking (later also finite temperature)
    ENABLE_N_EFA = .false. ; ENABLE_P_EFA = .false. ; !! default value
    If(keyblo(1).Ne.0) ENABLE_N_EFA = .true.
    If(keyblo(2).Ne.0) ENABLE_P_EFA = .true.
    !
    ! EFA arrays
    if(Allocated(F11N_c)) Deallocate(F11N_c,BF11N_c,h11N_c,Bh11N_c,P1N_c,P2N_c,GreenNefa, &
                                     efa_fN_factor1,efa_fN_factor2)
    if(ENABLE_N_EFA) Allocate(F11N_c(iQRPA_LEN),BF11N_c(iQRPA_LEN),h11N_c(iQRPA_LEN),Bh11N_c(iQRPA_LEN), & 
                              P1N_c(iQRPA_LEN),P2N_c(iQRPA_LEN),GreenNefa(iQRPA_LEN),efa_fN_factor1(iQRPA_LEN), &
                              efa_fN_factor2(iQRPA_LEN))
    if(Allocated(F11P_c)) Deallocate(F11P_c,BF11P_c,h11P_c,Bh11P_c,P1P_c,P2P_c,GreenPefa, &
                                     efa_fP_factor1,efa_fP_factor2)
    if(ENABLE_P_EFA) Allocate(F11P_c(iQRPA_LEN),BF11P_c(iQRPA_LEN),h11P_c(iQRPA_LEN),Bh11P_c(iQRPA_LEN), & 
                              P1P_c(iQRPA_LEN),P2P_c(iQRPA_LEN),GreenPefa(iQRPA_LEN),efa_fP_factor1(iQRPA_LEN), &
                              efa_fP_factor2(iQRPA_LEN))
    ! set efa arrays to zero, except energy denominators
    if(ENABLE_N_EFA) then; F11N_c = zero; BF11N_c = zero; h11N_c = zero; Bh11N_c = zero; P1N_c = zero; 
                            P2N_c = zero; GreenNefa = one; efa_fN_factor1 = one; efa_fN_factor2 = one; end if
    if(ENABLE_P_EFA) then; F11P_c = zero; BF11P_c = zero; h11P_c = zero; Bh11P_c = zero; P1P_c = zero; 
                            P2P_c = zero; GreenPefa = one; efa_fP_factor1 = one; efa_fP_factor2 = one; end if
    !---------------------------------------------------------------------
    ! Charging HFBTHO brin
    !---------------------------------------------------------------------
    alphamix=1.0_pr; iiter=0
    Call densit
    Call field	
    Call gamdel
    !---------------------------------------------------------------------
    ! make HFB matrices complex 
    !---------------------------------------------------------------------
    VqpN=Zero; UqpN=Zero; VqpP=Zero; UqpP=Zero
    i=0;
    Do ib=1,NB
       ND=ID(ib)
       Do N1=1,ND
         Do N2=1,ND
           i = i+1
           VqpN(i)=RVqpN(i)
           VqpP(i)=RVqpP(i)
           UqpN(i)=RUqpN(i)
           UqpP(i)=RUqpP(i)
         End Do
       End Do
    End Do
    !---------------------------------------------------------------------
    ! set U and V to simplex-y phase.
    !---------------------------------------------------------------------
    i=0;
    Do ib=1,NB
       ND=ID(ib)
       Do N1=1,ND
         Do N2=1,ND
           NSA=NS(ia(ib)+N2);
           if(NSA.gt.0) then
             LAM = (iHFB_OM(ib) -1)/2
           else
             LAM = -(iHFB_OM(ib) +1)/2
           end if
           phase = 1.0_pr
           if( mod(abs(Lam),2) .eq. 1 ) then
             phase = cmplx(0.0_pr,1.0_pr)
             if(lam.lt.0) phase = -phase
           end if
           i = i+1
           If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) phase = cmplx(1.0_pr,0.0_pr)  !! don't change phase for normal HO basis
           VqpN(i)=VqpN(i)*phase
           VqpP(i)=VqpP(i)*phase
           UqpN(i)=UqpN(i)*phase
           UqpP(i)=UqpP(i)*phase
         End Do
       End Do
    End Do
    !
    !---------------------------------------------------------------------
    ! use HFB pairing window for calculation of densities
    !---------------------------------------------------------------------
    USE_HFB_PWI = .true.
    !
    icacou_qrpa=0 !! recalculate coulomb
    !
    call qrpa_calculate_timeodd_CC
    !---------------------------------------------------------------------
    ! deallocate some HFBTHO arrrays ...
    !---------------------------------------------------------------------
    Deallocate(brout,brin,RVqpN,RUqpN,RVqpP,RUqpP)
    !---------------------------------------------------------------------
    ! calculate the external field
    !---------------------------------------------------------------------
    Call qrpa_f_gamdel
    !
  End Subroutine qrpaINI_CUT_ALLOCATE
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_N00(V,U,AN00,AN11)
    !---------------------------------------------------------------------
    ! The diagonal of N00 =  V'* V to AN00 
    ! The diagonal of N11 =  U'* U - V'* V  to AN11 
    !   AN00 =>  0 for particles,  1 for holes
    !   AN11 => -1 for particles, +1 for holes
    !
    ! Use this for the HFB matrices imported from HFBTHO
    !---------------------------------------------------------------------
    Implicit None
    Real(pr), Intent(in)  :: U(NUV),V(NUV)
    Real(pr), Intent(out) :: AN00(nqp),AN11(nqp)
    Real(pr)     :: REAUX(NUV),REBUX(NUV)
    Real(pr)     :: one,zero,onem
    Integer(ipr) :: i,ib,nd,ip,iq
    !
    one=1.0_pr; zero=0.0_pr; onem=-one
    !
    ip = 1; iq=1
    Do ib=1,NB
       ND=ID(ib)
       Call dsyrk('L','t',nd,nd,one, V(ip),nd,zero,REAUX(ip),nd) !  V'* V 
       Call dsyrk('L','t',nd,nd,one ,U(ip),nd,zero,REBUX(ip),nd) !  U'* U 
       Call dsyrk('L','t',nd,nd,onem,V(ip),nd,one ,REBUX(ip),nd) !- V'* V
       Do i=1,nd
          AN00(i+iq-1)=REAUX(ip+(i-1)*nd+i-1)
          AN11(i+iq-1)=REBUX(ip+(i-1)*nd+i-1)
       End Do
       ip = ip + nd*nd
       iq = iq + nd
    End Do
    !
  End Subroutine qrpa_N00
  !=================================================================================================================================
  !  
  !=================================================================================================================================
  Subroutine qrpa_test_HFB_UV(V,U,it)
    !---------------------------------------------------------------------
    ! calculates V' U - U' V and V' V + U' U; UV-real matrices. 
    !---------------------------------------------------------------------
    Implicit None
    Integer(ipr), Intent(in) :: it
    Real(pr),     Intent(in) :: U(NUV),V(NUV)
    Real(pr)     :: REAUX(NUV)
    Real(pr)     :: onem,sum1,qf1,qf1t,ef1,sum2,qf2,qf2t,ef2 
    Integer(ipr) :: ib,nd,ip,i,j,k 
    !
    onem=-one 
    !
    ip = 1; k=0; qf1t = zero; qf2t = zero
    Do ib=1,NB
       ND=ID(ib)       
       Call dgemm('t','n',nd,nd,nd,onem,U(ip),nd,V(ip),nd,zero,REAUX(ip),nd) !   V' U       
       Call dgemm('t','n',nd,nd,nd,one ,V(ip),nd,U(ip),nd,one ,REAUX(ip),nd) ! - U' V
       sum1 = zero
       Do i=1,ND
          Do j=1,ND
             sum1 = sum1 + REAUX(ip+(j-1)*ND+i-1)**2
          End Do
       End Do
       Call dgemm('t','n',nd,nd,nd,one ,U(ip),nd,U(ip),nd,zero,REAUX(ip),nd) !  V' V 
       Call dgemm('t','n',nd,nd,nd,one ,V(ip),nd,V(ip),nd,one ,REAUX(ip),nd) !+ U' U
       sum2 = zero
       Do i=1,ND
          Do j=1,ND
             sum2 = sum2 + REAUX(ip+(j-1)*ND+i-1)**2
          End Do
       End Do
       qf1 = sum1/Dfloat(Nd*Nd); qf1 = Sqrt(qf1)
       qf2 = sum2/Dfloat(ND);    qf2 = Sqrt(qf2)
       ef1 = dabs(qf1);          ef2 = Abs(qf2-one)
       If ((ef1.Gt.1.0d-11).Or.(ef2.Gt.1.0d-11)) Then
          Write(*,'(2x,a,2(a,d12.3,2x))') txb(ib),'Max VU-UV:',ef1,'VV+UU-1:',ef2  
          k=k+1    
       End If
       ip = ip + nd*nd
    End Do
    If(k.Eq.0.And.it.Eq.1) Write(*,'(2x,a,2(a,d12.3,2x))') ' NEUTRON HFB UV SOLUTION IS OK'
    If(k.Eq.0.And.it.Eq.2) Write(*,'(2x,a,2(a,d12.3,2x))') ' PROTON  HFB UV SOLUTION IS OK'
  End Subroutine qrpa_test_HFB_UV
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_green_function(omega) 
    !---------------------------------------------------------------------
    ! Calculates energy denomitors. With EFA calculates also the f-factors
    !
    ! GreenN_k1k2=1/(Eqp(k1)+Eqp(k2) +/- omega), all complex
    !---------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in) :: omega
    Integer(ipr) :: i,ib,idi,idf,k1,k2
    Integer(ipr), allocatable :: iHFB_EDS(:)
    Real(pr)     :: gf_AeNm,gf_AePm,gf_ReN,gf_ReP,Romega
    Complex(pr)  :: onec,gf_eN,gf_eP,ffn1,ffn2,ffp1,ffp2
    integer(ipr) :: ibln,iblp,k0n1,k0n2,k0p1,k0p2
    Complex(pr)  :: epsn,epsp
    !    
    Romega=Real(omega,kind=pr);  onec=1.0_pr
    gf_AeNm=1000.0_pr; gf_AePm=1000.0_pr
    !
    if(.not.allocated(iHFB_EDS)) allocate(iHFB_EDS(iHFB_NB))
    !
    i=0;
    do ib=1,NB
      IDI=ID(ib)
      iHFB_EDS(ib) = i
      i=i+IDI
    end do
    !
    i=0;
    !
    ! avoid division by zero with efa (in Eqp_i - Eqp_j part)
    epsn = zero; epsp = zero;
    If(ENABLE_N_EFA) epsn = spacing(1.0_pr)
    If(ENABLE_P_EFA) epsp = spacing(1.0_pr)  
    !
    Do ib=1,iQRPA_NB
       Do k1=1,iQRPA_IDI(ib)
         Do k2=1,iQRPA_IDF(ib)
             i=i+1
             ! neutrons 
             gf_eN=REqpN(iHFB_EDS(iQRPA_BHI(ib))+k1)+REqpN(iHFB_EDS(iQRPA_BHF(ib))+k2); GreenNminus(i)=onec/(gf_eN-omega); GreenNplus(i) =onec/(gf_eN+omega)
             ! protons 
             gf_eP=REqpP(iHFB_EDS(iQRPA_BHI(ib))+k1)+REqpP(iHFB_EDS(iQRPA_BHF(ib))+k2); GreenPminus(i)=onec/(gf_eP-omega); GreenPplus(i) =onec/(gf_eP+omega)
             !
             ! blocking
             ibln=bloblo(keyblo(1),1)  !! hfb block of blocked state
             iblp=bloblo(keyblo(2),2) 
             ffn1 = onec; ffn2 = zero; ffp1 = onec; ffp2 = zero;
             If(ENABLE_N_EFA .and. ibln.Eq. iQRPA_BHI(ib) ) then; K0n1=blo123d(1); else; k0n1 = 0; end if
             If(ENABLE_N_EFA .and. ibln.Eq. iQRPA_BHF(ib) ) then; K0n2=blo123d(1); else; k0n2 = 0; end if
             If(ENABLE_P_EFA .and. iblp.Eq. iQRPA_BHI(ib) ) then; K0p1=blo123d(2); else; k0p1 = 0; end if
             If(ENABLE_P_EFA .and. iblp.Eq. iQRPA_BHF(ib) ) then; K0p2=blo123d(2); else; k0p2 = 0; end if
             !
             If(ENABLE_N_EFA .and. k1.eq.k0n1) then; ffn1 = ffn1 - half; ffn2 = ffn2 -half; end if;
             If(ENABLE_N_EFA .and. k2.eq.k0n2) then; ffn1 = ffn1 - half; ffn2 = ffn2 +half; end if;
             If(ENABLE_P_EFA .and. k1.eq.k0p1) then; ffp1 = ffp1 - half; ffp2 = ffp2 -half; end if;
             If(ENABLE_P_EFA .and. k2.eq.k0p2) then; ffp1 = ffp1 - half; ffp2 = ffp2 +half; end if;
             !              
             ! proton poles
             gf_ReP=gf_eP 
             If(gf_AePm.Ge.Abs(gf_ReP-Romega)) Then
                gf_AePm=Abs(gf_ReP-Romega); gfpoles_imP=i
             Endif
             ! neutron poles
             gf_ReN=gf_eN
             If(gf_AeNm.Ge.Abs(gf_ReN-Romega)) Then
                gf_AeNm=Abs(gf_ReN-Romega); gfpoles_imN=i
             Endif
             ! 
             ! efa arrays CHECK THE SIGN OF ENERGY W.R.T F-FACTORS
             gf_eN=REqpN(iHFB_EDS(iQRPA_BHI(ib))+k1)-REqpN(iHFB_EDS(iQRPA_BHF(ib))+k2)+epsn;
             gf_eP=REqpP(iHFB_EDS(iQRPA_BHI(ib))+k1)-REqpP(iHFB_EDS(iQRPA_BHF(ib))+k2)+epsp;
             if(ENABLE_N_EFA) then
               GreenNefa(i) = onec/(gf_eN-omega);
               efa_fN_factor1(i) = ffn1;  
               efa_fN_factor2(i) = ffn2;
             end if
             if(ENABLE_P_EFA) then 
               GreenPefa(i) = onec/(gf_eP-omega);
               efa_fP_factor1(i) = ffp1;  
               efa_fP_factor2(i) = ffp2;
             end if
          Enddo
       Enddo
    Enddo

!WRITE(*,*) "green:"
!write(*,*) "efa1:"
!WRITE(*,'(20(1x,f5.2))') REAL(efa_fN_factor1)
!write(*,*) "efa2:"
!WRITE(*,'(20(1x,f5.2))') REAL(efa_fN_factor2)
!call qrpa_printprop("efa_fN_factor1",efa_fN_factor1)
!call qrpa_printprop("efa_fN_factor2",efa_fN_factor2)
!call qrpa_printprop("GreenNefa",GreenNefa)

    deallocate(iHFB_EDS)
  End Subroutine qrpa_green_function
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_xy 
    !---------------------------------------------------------------------
    ! without EFA:
    ! XN_c=-hN_c*GreenNminus; YN_c=-BhN_c*GreenNplus         
    ! XP_c=-hP_c*GreenPminus; YP_c=-BhP_c*GreenPplus
    !
    ! with EFA: X,Y,P1,P2 amplitudes, see separate documentation
    !---------------------------------------------------------------------
    Implicit None
    Integer(ipr) :: i
    Do i=1,iQRPA_LEN
      !
      ! neutrons
      if(.not.ENABLE_N_EFA) then
        XN_c(i)=-hN_c(i)*GreenNminus(i); YN_c(i)=-BhN_c(i)*GreenNplus(i)    !! normal even-even qrpa
      else
        XN_c(i)  = -hN_c(i)*GreenNminus(i)*efa_fN_factor1(i)
        YN_c(i)  = -BhN_c(i)*GreenNplus(i)*efa_fN_factor1(i)
        P1N_c(i) = -h11N_c(i)*GreenNefa(i)*efa_fN_factor2(i)
        P2N_c(i) = -Bh11N_c(i)*GreenNefa(i)*efa_fN_factor2(i)
      end if
      !
      ! protons
      if(.not.ENABLE_P_EFA) then
        XP_c(i)=-hP_c(i)*GreenPminus(i); YP_c(i)=-BhP_c(i)*GreenPplus(i)    !! normal even-even qrpa
      else
        XP_c(i)  = -hP_c(i)*GreenPminus(i)*efa_fP_factor1(i)
        YP_c(i)  = -BhP_c(i)*GreenPplus(i)*efa_fP_factor1(i)
        P1P_c(i) = -h11P_c(i)*GreenPefa(i)*efa_fP_factor2(i)
        P2P_c(i) = -Bh11P_c(i)*GreenPefa(i)*efa_fP_factor2(i)
      end if
    End Do

!if(ENABLE_N_EFA) then
!WRITE(*,*) "xy:"
!call qrpa_printprop("XN",XN_c)
!call qrpa_printprop("YN",YN_c)
!call qrpa_printprop("P1N",P1N_c)
!call qrpa_printprop("P2N",P2N_c)
!call qrpa_printprop("GreenNefa",GreenNefa)
!call qrpa_printprop("efa_fN_factor2",efa_fN_factor2)
!call qrpa_printprop("efa_fN_factor2*GreenNefa",efa_fN_factor2*GreenNefa)
!call qrpa_printprop("h11n",h11n_c)
!call qrpa_printprop("bh11n",bh11n_c)

!write(*,*) "f11:"
!WRITE(*,'(20(1x,f5.2))') REAL(F11N_c)
!write(*,*) "f11*efa2:"
!WRITE(*,'(20(1x,f5.2))') REAL(efa_fN_factor2*F11N_c)


!WRITE(*,*) "efa_fN_factor2*GreenNefa:"
!WRITE(*,*) REAL(efa_fN_factor2*GreenNefa)
!WRITE(*,*) "Re h11n"
!WRITE(*,*) REAL(h11n_c)
!WRITE(*,*) "Im h11n"
!WRITE(*,*) Imag(h11n_c)
!end if

  End Subroutine  qrpa_xy
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_TriProd(op1,op2,op3,U,X,V,Fout)
    !---------------------------------------------------------------------
    ! Computes: Fout = op1(U) . op2(X) . op3(V)
    ! Use lowercase characters for op's!
    ! U and V blocks are assumed to be square matrices and X to be general block matrix
    ! Op's as in ZGEMM, but only lowercase. New op '*' added for complex conjugation
    !---------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)  :: U(iHFB_LEN),V(iHFB_LEN),X(iQRPA_LEN)
    Complex(pr), Intent(out) :: Fout(iQRPA_LEN)
    Character, Intent(in)    :: op1,op2,op3
    Character                :: op1l,op2l,op3l
    Complex(pr)              :: czero,cone,conem
    Complex(pr), allocatable :: aux1(:),aux2(:),aux3(:),aux(:)
    Integer(ipr) :: ib,nd,ipa,ipq,ipt,ibHFBi,ibHFBf,ibTHOi,ibTHOf
    !
    czero=0.0_pr; cone=1.0_pr; conem=-cone
    !
    if ((op2 .eq. 'T') .or. (op2 .eq. 'C')) then
      STOP "Use lowercase chars"
    end if
    !
    If(allocated(aux1)) Deallocate(aux1,aux2,aux3)
    Allocate(aux1(iHFB_LEN),aux2(iQRPA_LEN),aux3(iHFB_LEN),aux(iQRPA_LEN))

    op1l = op1; op2l = op2; op3l = op3; 
    aux = czero;

    if(op1 .eq. '*') then ; aux1=Conjg(U); op1l = 'n'; else ; aux1=U; end if
    if(op2 .eq. '*') then ; aux2=Conjg(X); op2l = 'n'; else ; aux2=X; end if
    if(op3 .eq. '*') then ; aux3=Conjg(V); op3l = 'n'; else ; aux3=V; end if

    ! first X.V
    Do ib=1,iQRPA_NB
      ibHFBi = iQRPA_BHI(ib) ; ibHFBf = iQRPA_BHF(ib)
      ibTHOi = iQRPA_BLI(ib) ; ibTHOf = iQRPA_BLF(ib)
      ipa = iHFB_DS(ibHFBi) ; ipq = iQRPA_DS(ib) ; ipt = iQRPA_DS(iQRPA_TRANSPOSE(ib))
      !
      if ((op2l .eq. 't') .or. (op2l .eq. 'c')) then
        Call ZGEMM(op2l,op3l,iHFB_ID(ibHFBf),iHFB_ID(ibHFBi),iHFB_ID(ibHFBi),cone,aux2(ipt),iHFB_ID(ibHFBi),aux3(ipa),iHFB_ID(ibHFBi),czero,aux(ipq),iHFB_ID(ibHFBf)) !! XT * op3(V) or XH * op3(V)
      else
        Call ZGEMM(op2l,op3l,iHFB_ID(ibHFBf),iHFB_ID(ibHFBi),iHFB_ID(ibHFBi),cone,aux2(ipq),iHFB_ID(ibHFBf),aux3(ipa),iHFB_ID(ibHFBi),czero,aux(ipq),iHFB_ID(ibHFBf)) !! X * op3(V)
      end if
    End Do
    ! 
    ! then U.(XV)
    Do ib=1,iQRPA_NB
      ibHFBi = iQRPA_BHI(ib) ; ibHFBf = iQRPA_BHF(ib)
      ibTHOi = iQRPA_BLI(ib) ; ibTHOf = iQRPA_BLF(ib)
      ipa = iHFB_DS(ibHFBf) ; ipq = iQRPA_DS(ib) ; ipt = iQRPA_DS(iQRPA_TRANSPOSE(ib))
      Call ZGEMM(op1l,'n',iHFB_ID(ibHFBf),iHFB_ID(ibHFBi),iHFB_ID(ibHFBf),cone,aux1(ipa),iHFB_ID(ibHFBf),aux(ipq),iHFB_ID(ibHFBf),czero,Fout(ipq),iHFB_ID(ibHFBf)) !! op1(U) * (XV) to Fout
    End Do
    Deallocate(aux1,aux2,aux3,aux)
  End Subroutine qrpa_TriProd
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_DiProdDiag(op1,A,B,Fout)
    !---------------------------------------------------------------------
    ! Computes diagonal blocks of Fout = op1(A) . B, where A and B have
    ! the same block structure as QRPA matrices.
    ! Op as in ZGEMM, but only lowercase. New op '*' added for complex conjugation
    ! NOTE: Off-diagonal blocks are NOT calculated!
    !---------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)  :: A(iQRPA_LEN),B(iQRPA_LEN)
    Complex(pr), Intent(out) :: Fout(iHFB_LEN)
    Character, Intent(in)    :: op1
    Complex(pr)              :: czero,cone,conem
    Complex(pr), allocatable :: aux1(:),aux2(:)
    Character                :: op1l
    Integer(ipr) :: ib,ic,ict,ihdata,iadata,ibdata
    !
    czero=0.0_pr; cone=1.0_pr; conem=-cone
    !
    If(allocated(aux1)) Deallocate(aux1,aux2)
    Allocate(aux1(iHFB_LEN),aux2(iQRPA_LEN))
    aux1=czero;
    !
    op1l = op1
    if(op1 .eq. '*') then ; aux2=Conjg(A); op1l = 'n'; else ; aux2=A; end if
    !
    Do ib=1,iHFB_NB
      ihdata = iHFB_DS(ib)
      Do ic=1,iQRPA_NB
        if(iQRPA_BHF(ic).eq.ib) then
          ict = iQRPA_TRANSPOSE(ic)
          iadata = iQRPA_DS(ic)
          ibdata = iQRPA_DS(ict)
          !
          if ((op1l .eq. 't') .or. (op1l .eq. 'c')) then
            Call ZGEMM(op1l,'n',iHFB_ID(ib),iQRPA_IDF(ic),iQRPA_IDI(ic),cone,aux2(ibdata),iQRPA_IDI(ic),B(ibdata),iQRPA_IDI(ic),cone,aux1(ihdata),iHFB_ID(ib)) 
          else
            Call ZGEMM(op1l,'n',iHFB_ID(ib),iQRPA_IDF(ic),iQRPA_IDI(ic),cone,aux2(iadata),iHFB_ID(ib),B(ibdata),iQRPA_IDI(ic),cone,aux1(ihdata),iHFB_ID(ib)) 
          end if
          !
        End If
      End Do
    End Do
    Fout = aux1
    Deallocate(aux1,aux2)
  End Subroutine qrpa_DiProdDiag
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_GiveTransposed(Fin,Fout)
    !
    ! Gives transpose of QRPA-block type of matrix. Output overwritten in Fout.
    ! Not used, but could be usedful for debug.
    !
    Complex(pr), Intent(in)  :: Fin(iQRPA_LEN)
    Complex(pr), Intent(out) :: Fout(iQRPA_LEN)
    Integer(ipr) :: i,j,ib,ibt,iadata,ibdata
    !
    Do ib=1,iQRPA_NB
      ibt = iQRPA_TRANSPOSE(ib)
      iadata = iQRPA_DS(ib)
      ibdata = iQRPA_DS(ibt)
      Do i=1,iQRPA_IDI(ib)
        Do j=1,iQRPA_IDF(ib)
          Fout(ibdata-1 +i +(j-1)*iQRPA_IDI(IB) ) = Fin(iadata-1 +j +(i-1)*iQRPA_IDF(IB) )
        End Do
      End Do
    End Do
  End Subroutine qrpa_GiveTransposed
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Function qrpa_HFBtrace(Fin)
    !
    ! gives trace of a HFB-block type of matrix
    !
    Complex(pr), Intent(in)  :: Fin(iHFB_LEN)
    Complex(pr) :: qrpa_HFBtrace
    Complex(pr) :: trtmp,ciptrace
    integer(ipr) :: ib,nd,ip
    !
    trtmp = 0.0_pr
    Do ib = 1,iHFB_NB
      nd = iHFB_ID(ib)
      ip = iHFB_DS(ib)
      Call qrpa_ctrace(Fin(ip),nd,ciptrace)
      trtmp = trtmp + ciptrace
    End Do
    qrpa_HFBtrace = trtmp
    !
  End Function qrpa_HFBtrace
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_H20(V,U,H,D,BH,BD)
    !
    ! Calculates H20 and H02 matrix elements. Output to H and BH, which are overwritten
    !
    Implicit None
    Complex(pr), Intent(in)     :: U(iHFB_LEN),V(iHFB_LEN),D(iQRPA_LEN),BD(iQRPA_LEN)
    Complex(pr), Intent(inout)  :: H(iQRPA_LEN),BH(iQRPA_LEN)
    Complex(pr), Allocatable    :: temp(:)
    !
    if(allocated(temp)) Deallocate(temp)
    Allocate(temp(iQRPA_LEN))
    !
    ! V and U are complex
    call qrpa_TriProd('c','n','n',U,H,V,qrpa_AUX)
    call qrpa_TriProd('c','t','n',V,BH,U,qrpa_AUX20)
    call qrpa_TriProd('c','c','n',V,BD,V,qrpa_BUX)
    call qrpa_TriProd('c','n','n',U,D,U,qrpa_BUX20)

    temp = -qrpa_AUX -qrpa_AUX20 -qrpa_BUX +qrpa_BUX20

    call qrpa_TriProd('t','t','*',U,H,V,qrpa_AUX)
    call qrpa_TriProd('t','n','*',V,BH,U,qrpa_AUX20)
    call qrpa_TriProd('t','t','*',V,D,V,qrpa_BUX)
    call qrpa_TriProd('t','*','*',U,BD,U,qrpa_BUX20)

    H = temp
    BH = -qrpa_AUX -qrpa_AUX20 -qrpa_BUX +qrpa_BUX20
    Deallocate(temp)    
  End Subroutine qrpa_H20
  !=================================================================================================================================
  !   
  !=================================================================================================================================
  subroutine qrpa_H11(V,U,H,D,BH,BD,H11,BH11)
    !
    ! Calculates H11_1 and H11_2 matrix elements.
    !
    Implicit None
    Complex(pr), Intent(in)     :: U(iHFB_LEN),V(iHFB_LEN),D(iQRPA_LEN),BD(iQRPA_LEN),H(iQRPA_LEN),BH(iQRPA_LEN)
    Complex(pr), Intent(out)    :: H11(iQRPA_LEN),BH11(iQRPA_LEN)
    !
    ! V and U are complex
    call qrpa_TriProd('c','n','n',U,H,U,qrpa_AUX)
    call qrpa_TriProd('c','t','n',V,BH,V,qrpa_AUX20)
    call qrpa_TriProd('c','n','n',U,D,V,qrpa_BUX)
    call qrpa_TriProd('c','t','n',V,BD,U,qrpa_BUX20)

    H11 = qrpa_AUX -qrpa_AUX20 +qrpa_BUX -qrpa_BUX20

    call qrpa_TriProd('t','n','*',U,BH,U,qrpa_AUX)
    call qrpa_TriProd('t','t','*',V,H,V,qrpa_AUX20)
    call qrpa_TriProd('t','t','*',U,D,V,qrpa_BUX)
    call qrpa_TriProd('t','n','*',V,BD,U,qrpa_BUX20)

    BH11 = qrpa_AUX -qrpa_AUX20 +qrpa_BUX -qrpa_BUX20
    !
  end subroutine qrpa_H11
  !=================================================================================================================================
  !   
  !=================================================================================================================================
  Subroutine qrpa_strength(X,Y,P1,P2,F20,F02,F11,F11b,efa,ImStrength,ReStrength)
    !
    ! Computes imag and real parts of qrpa strength function: sum (F20*X +F02*Y)
    ! Note: transition strength can be obtained from imaginary part by multiplying with -1/pi
    !
    Implicit None
    Complex(pr), Intent(in)  :: X(iQRPA_LEN),Y(iQRPA_LEN),F20(iQRPA_LEN),F02(iQRPA_LEN)
    Complex(pr), Intent(in)  :: P1(iQRPA_LEN),P2(iQRPA_LEN),F11(iQRPA_LEN),F11b(iQRPA_LEN)
    logical, Intent(in)      :: efa
    Real(pr), Intent(out)    :: ImStrength,ReStrength
    Complex(pr)              :: ctrace
    Complex(pr), Allocatable :: temp1(:)
    !
    if(Allocated(temp1)) Deallocate(temp1)
    Allocate(temp1(iQRPA_LEN))
    !
    temp1 = 0.0_pr ; ImStrength = 0.0_pr; ReStrength = 0.0_pr
    temp1 = conjg(F20)*X +conjg(F02)*Y
!write(*,*) "F20 part terms: ",sum(conjg(F20)*X),sum(conjg(F02)*Y)

    if(efa) temp1 = temp1 + conjg(F11)*P1 + conjg(F11b)*P2
    !
!if(efa) write(*,*) "F11 part terms: ",sum(conjg(F11)*P1),sum(conjg(F11b)*P2)
!if(efa) call qrpa_printProp("P1",P1)
!if(efa) call qrpa_printProp("P2",P2)
!if(efa) then ; qrpa_aux = P1-P2 ; call qrpa_printProp("P1-P2",qrpa_aux) ; end if
!if(efa) then ; qrpa_aux = P1+P2 ; call qrpa_printProp("P1+P2",qrpa_aux) ; end if
!if(efa) then ; qrpa_aux = F11-F11b ; call qrpa_printProp("F11-F11b",qrpa_aux) ; end if
!if(efa) then ; qrpa_aux = F11+F11b ; call qrpa_printProp("F11+F11b",qrpa_aux) ; end if



    ImStrength = Imag( sum(temp1) )
    ReStrength = Real( sum(temp1) ,kind=pr)
    Deallocate(temp1)
  End Subroutine qrpa_strength
  !=================================================================================================================================
  ! 
  !=================================================================================================================================
  Subroutine qrpa_F20(V,U,F1,F2,FF1,FF2)
    !---------------------------------------------------------------------
    ! Computes: 
    ! FF1 = F20 = -u^\dag f1 V - v^\dag f2^T u
    ! FF2 = F02 = -u^T f1^T V^* - v^T f2 u^*
    !---------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)  :: U(iHFB_LEN),V(iHFB_LEN),F1(iQRPA_LEN),F2(iQRPA_LEN)
    Complex(pr), Intent(out) :: FF1(iQRPA_LEN),FF2(iQRPA_LEN)
    Complex(pr), allocatable :: aux1(:),aux2(:)
    Integer(ipr) :: ib,nd,ipa,ipq,ipt,ibHFBi,ibHFBf,ibTHOi,ibTHOf
    !
    If(allocated(aux1)) Deallocate(aux1,aux2)
    Allocate(aux1(iQRPA_LEN),aux2(iQRPA_LEN))
    !
    call qrpa_TriProd('c','n','n',U,F1,V,aux1)
    call qrpa_TriProd('c','t','n',V,F2,U,aux2)
    FF1 = -aux1-aux2
    !
    call qrpa_TriProd('t','*','*',U,F1,V,aux1)
    call qrpa_TriProd('t','c','*',V,F2,U,aux2)
    FF2 = -aux1-aux2
    !
    Deallocate(aux1,aux2)
  End Subroutine qrpa_F20
  !=================================================================================================================================
  ! 
  !=================================================================================================================================
  Subroutine qrpa_F11(V,U,F1,F2,FF1,FF2)
    !---------------------------------------------------------------------
    ! Computes: 
    ! FF1 = F11  = u^\dag f1 u - v^\dag f2^\dag v
    ! FF2 = F11b = u^T f2^* u^*  - v^T f1^T v^*
    !---------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)  :: U(iHFB_LEN),V(iHFB_LEN),F1(iQRPA_LEN),F2(iQRPA_LEN)
    Complex(pr), Intent(out) :: FF1(iQRPA_LEN),FF2(iQRPA_LEN)
    Complex(pr), allocatable :: aux1(:),aux2(:)
    Integer(ipr) :: ib,nd,ipa,ipq,ipt,ibHFBi,ibHFBf,ibTHOi,ibTHOf
    !
    If(allocated(aux1)) Deallocate(aux1,aux2)
    Allocate(aux1(iQRPA_LEN),aux2(iQRPA_LEN))
    !
    call qrpa_TriProd('c','n','n',U,F1,U,aux1)
    call qrpa_TriProd('c','t','n',V,F2,V,aux2)
    FF1 = aux1-aux2
    !
    call qrpa_TriProd('t','n','*',U,F2,U,aux1)
    call qrpa_TriProd('t','t','*',V,F1,V,aux2)
    FF2 = aux1-aux2
!
! MOD HERE !
!call qrpa_TriProd('c','n','n',U,F1,U,aux1)
!call qrpa_TriProd('c','t','n',V,F2,V,aux2)
!FF1=aux1
!FF2=aux2


    !
    Deallocate(aux1,aux2)
  End Subroutine qrpa_F11
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine Calculate_QRPA_H20
    !---------------------------------------------------------------------
    ! calculates H20 using:               Output: 
    !       - HFB matrices                 - H20, H02
    !       - QRPA amplitudes              - H11, H11b (if needed)
    !
    !---------------------------------------------------------------------
    Implicit None
    integer(ipr) :: i !! remove later
    !---------------------------------------------------------------------
    ! calculates h,d,Bd for neutrons and protons
    !---------------------------------------------------------------------
    !
    Call qrpa_DENSIT
    !
    Call qrpa_coulom
    !
    Call qrpa_field 
    !
    ! matrices h1,h2,d,Bd w.r.t basis state indices
    Call qrpa_gamdel
    !
    ! H11, w.r.t quasiparticel indices
    if(ENABLE_N_EFA) Call qrpa_H11(VqpN,UqpN,hN_c,dN_c,BhN_c,BdN_c,h11N_c,Bh11N_c)
    if(ENABLE_P_EFA) Call qrpa_H11(VqpP,UqpP,hP_c,dP_c,BhP_c,BdP_c,h11P_c,Bh11P_c)
    !---------------------------------------------------------------------
    ! H20,BH20 as h,Bh. Note: subroutines below overwrite h_c,d_c,Bh_c,Bd_c
    !---------------------------------------------------------------------
    Call qrpa_H20(VqpN,UqpN,hN_c,dN_c,BhN_c,BdN_c)
    Call qrpa_H20(VqpP,UqpP,hP_c,dP_c,BhP_c,BdP_c)
    !
    !---------------------------------------------------------------------
    If(Print_Screen.And.IDEBUG.Gt.10) Then 
       Write(*,*) 
       Write(*,*) 'maxval( HN)=',Maxval(Real( hN_c)),Maxval(Imag( hN_c))
       Write(*,*) 'maxval( HP)=',Maxval(Real( hP_c)),Maxval(Imag( hP_c))
       Write(*,*) 'maxval(BHN)=',Maxval(Real(BhN_c)),Maxval(Imag(BhN_c))
       Write(*,*) 'maxval(BHP)=',Maxval(Real(BhP_c)),Maxval(Imag(BhP_c))
       Write(*,*) 
       Write(*,*) 'minval( HN)=',Minval(Real( hN_c)),Minval(Imag( hN_c))
       Write(*,*) 'minval( HP)=',Minval(Real( hP_c)),Minval(Imag( hP_c))
       Write(*,*) 'minval(BHN)=',Minval(Real(BhN_c)),Minval(Imag(BhN_c))
       Write(*,*) 'minval(BHP)=',Minval(Real(BhP_c)),Minval(Imag(BhP_c))
       Write(*,*) 
       Write(*,*) 'maxval(FhN)=',Maxval(Real(FhN_c))
       Write(*,*) 'maxval(FhP)=',Maxval(Real(FhP_c))
       Write(*,*) 'minval(FhN)=',Minval(Real(FhN_c))
       Write(*,*) 'minval(FhP)=',Minval(Real(FhP_c))
       Write(*,*) 
       Write(*,*) 'FhN_c(1)=',FhN_c(1:5)
       Write(*,*) ' HN_c(1)=', hN_c(1:5)
       Write(*,*) 'BHN_c(1)=',BhN_c(1:5)
       Write(*,*) 
    End If
    !---------------------------------------------------------------------
    ! add external fields
    !---------------------------------------------------------------------
    !
    ! neutrons
    hN_c=hN_c+FhN_c; BhN_c=BhN_c+BFhN_c
    if(ENABLE_N_EFA) then ; h11N_c = h11N_c+F11N_c;  Bh11N_c = Bh11N_c+BF11N_c; end if
    if(ENABLE_N_EFA.and.USE_ONLY_F11EFA) then ; h11N_c = F11N_c;  Bh11N_c = BF11N_c; end if
    !
    ! protons
    hP_c=hP_c+FhP_c; BhP_c=BhP_c+BFhP_c
    if(ENABLE_P_EFA) then ; h11P_c = h11P_c+F11P_c;  Bh11P_c = Bh11P_c+BF11P_c; end if
    if(ENABLE_P_EFA.and.USE_ONLY_F11EFA) then ; h11P_c = F11P_c;  Bh11P_c = BF11P_c; end if
    !
  End Subroutine Calculate_QRPA_H20
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_ctrace(A,N,ctrace)
    !---------------------------------------------------------------------
    ! block trace - expand A(ip) to A(N,N) 
    !---------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)  :: A(N,N)
    Complex(pr), Intent(out) :: ctrace
    Integer(ipr) :: N,I 
    ctrace = 0.0_pr
    Do I=1,N
       ctrace = ctrace + A(I,I)
    End Do
  End Subroutine qrpa_ctrace
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_broyden(niter)
    !-------------------------------------
    ! Broyden mixing, could be optimized
    !-------------------------------------
    Implicit None
    Integer(ipr), Intent(in)  :: niter
    Integer(ipr) :: id,nd,i,j,k1,k2 

    ! determine broyden array dimension
    id=8*iQRPA_LEN
    if(ENABLE_N_EFA) id = id +4*iQRPA_LEN
    if(ENABLE_P_EFA) id = id +4*iQRPA_LEN
    !
    If(niter.Eq.0) Then
       ! zero iteration
       qrpa_bbroyden='i'; si=1.0_pr;
       If(Allocated(qrpa_broin)) Deallocate(qrpa_broin,qrpa_broout)    
       Allocate(qrpa_broin(id),qrpa_broout(id)); qrpa_broin=0.0_pr
    End If
    !
    ! to real variable
    Do i=1,iQRPA_LEN
      qrpa_broout(i+0*iQRPA_LEN)=Real(XN_c(i),kind=pr)
      qrpa_broout(i+1*iQRPA_LEN)=Imag(XN_c(i))
      qrpa_broout(i+2*iQRPA_LEN)=Real(XP_c(i),kind=pr)
      qrpa_broout(i+3*iQRPA_LEN)=Imag(XP_c(i)) 
      qrpa_broout(i+4*iQRPA_LEN)=Real(YN_c(i),kind=pr)
      qrpa_broout(i+5*iQRPA_LEN)=Imag(YN_c(i))
      qrpa_broout(i+6*iQRPA_LEN)=Real(YP_c(i),kind=pr)
      qrpa_broout(i+7*iQRPA_LEN)=Imag(YP_c(i))
      if(ENABLE_N_EFA) then
        qrpa_broout(i+8*iQRPA_LEN)=Real(P1N_c(i),kind=pr)
        qrpa_broout(i+9*iQRPA_LEN)=Imag(P1N_c(i))
        qrpa_broout(i+10*iQRPA_LEN)=Real(P2N_c(i),kind=pr)
        qrpa_broout(i+11*iQRPA_LEN)=Imag(P2N_c(i))
      end if
      if(ENABLE_P_EFA .and. .not.ENABLE_N_EFA) then
        qrpa_broout(i+8*iQRPA_LEN)=Real(P1P_c(i),kind=pr)
        qrpa_broout(i+9*iQRPA_LEN)=Imag(P1P_c(i))
        qrpa_broout(i+10*iQRPA_LEN)=Real(P2P_c(i),kind=pr)
        qrpa_broout(i+11*iQRPA_LEN)=Imag(P2P_c(i))
      end if
      if(ENABLE_P_EFA .and. ENABLE_N_EFA) then
        qrpa_broout(i+12*iQRPA_LEN)=Real(P1P_c(i),kind=pr)
        qrpa_broout(i+13*iQRPA_LEN)=Imag(P1P_c(i))
        qrpa_broout(i+14*iQRPA_LEN)=Real(P2P_c(i),kind=pr)
        qrpa_broout(i+15*iQRPA_LEN)=Imag(P2P_c(i))
      end if
    End Do
    !
    ! Broyden
    Call qrpa_broyden_method(id,qrpa_broout,qrpa_broin,qrpa_alphamix,si,niter,qrpa_nbroyden,qrpa_bbroyden)
    If(niter.Eq.0) si=1.0_pr
    !
    Do i=1,iQRPA_LEN
      XN_c(i)=Cmplx(qrpa_broin(i+0*iQRPA_LEN),qrpa_broin(i+1*iQRPA_LEN),kind=pr)
      XP_c(i)=Cmplx(qrpa_broin(i+2*iQRPA_LEN),qrpa_broin(i+3*iQRPA_LEN),kind=pr)
      YN_c(i)=Cmplx(qrpa_broin(i+4*iQRPA_LEN),qrpa_broin(i+5*iQRPA_LEN),kind=pr)
      YP_c(i)=Cmplx(qrpa_broin(i+6*iQRPA_LEN),qrpa_broin(i+7*iQRPA_LEN),kind=pr)
      if(ENABLE_N_EFA) then
        P1N_c(i)=Cmplx(qrpa_broin(i+8*iQRPA_LEN),qrpa_broin(i+9*iQRPA_LEN),kind=pr)
        P2N_c(i)=Cmplx(qrpa_broin(i+10*iQRPA_LEN),qrpa_broin(i+11*iQRPA_LEN),kind=pr)
      end if
      if(ENABLE_P_EFA .and. .not.ENABLE_N_EFA) then
        P1P_c(i)=Cmplx(qrpa_broin(i+8*iQRPA_LEN),qrpa_broin(i+9*iQRPA_LEN),kind=pr)
        P2P_c(i)=Cmplx(qrpa_broin(i+10*iQRPA_LEN),qrpa_broin(i+11*iQRPA_LEN),kind=pr)
      end if
      if(ENABLE_P_EFA .and. ENABLE_N_EFA) then
        P1P_c(i)=Cmplx(qrpa_broin(i+12*iQRPA_LEN),qrpa_broin(i+13*iQRPA_LEN),kind=pr)
        P2P_c(i)=Cmplx(qrpa_broin(i+14*iQRPA_LEN),qrpa_broin(i+15*iQRPA_LEN),kind=pr)
      end if
    End Do
  End Subroutine qrpa_broyden
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_broyden_method(N,vout,vin,alpha,si,iter,M,bbroyden)
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
    !         (M=0 Linear mixing)
    ! Output:
    !  si     MaxVal(Abs(vout-vin))
    !  vin    Broyden/Linear mixing result
    !  vout   vout-vin
    !  bbroyden='B' Broyden mixing, ='L' Linear mixing
    !---------------------------------------------------------------------
    Implicit None
    Integer(ipr),    Intent(In)     :: N,iter,M
    Real(pr),        Intent(In)     :: alpha
    Real(pr),        Intent(Out)    :: si
    Character(1),    Intent(Out)    :: bbroyden
    Real(pr),        Intent(InOut)  :: vout(N),vin(N)  
    Integer(ipr)                    :: i,j,iter_used,ipos,inext
    Integer(ipr), Allocatable,Save  :: iwork(:)
    Real(pr),    Allocatable, Save  :: beta(:,:),work(:)
    Real(pr),    Allocatable, Save  :: df(:,:),dv(:,:),curv(:)
    Real(pr),                 Save  :: w0
    Real(pr)                        :: DDOT,DNRM2,normi,gamma,curvature,sf
    !
    sf=-1.0_pr; Call DAXPY(N,sf,vin,1,vout,1)         ! vout = vout - vin
    si=Maxval(Abs(vout))    
    !---------------------------------------------------------------------
    ! Linear mixing
    !---------------------------------------------------------------------
    If(M.Eq.0.Or.iter.Eq.0) Then
       bbroyden='L'; Call DAXPY(N,alpha,vout,1,vin,1) ! vin = vin + alpha*vout         
       If(Print_Screen.And.IDEBUG.Gt.10) Then 
          If(iter.Eq.0) Write(6,*) '  Linear mixing (alpha) : ',alpha
       Endif
       Return
    End If
    !---------------------------------------------------------------------
    ! Broyden mixing
    !---------------------------------------------------------------------
    bbroyden='B'
    iter_used=Min(iter-1,M); ipos=iter-1-((iter-2)/M)*M; inext=iter-((iter-1)/M)*M
    If (iter.Eq.1) Then       
       w0=0.010_pr
       If(Allocated(curv)) Deallocate(curv,df,dv,beta,work,iwork)
       Allocate(curv(N),df(N,M),dv(N,M),beta(M,M),work(M),iwork(M))
       If(Print_Screen.And.IDEBUG.Gt.10) Then 	   
          Write(6,'(a,i3,3(2x,f18.8),a)') '   Broyden mixing (M,alpha,w0,mem) : '  &
               ,M,alpha,w0,(2*N*M+N)*8._pr/1.e6,' MB'
       End If
    Else
       df(:,ipos)=vout(:)-df(:,ipos) 
       dv(:,ipos)= vin(:)-dv(:,ipos)
       Normi=1.0_pr/DNRM2(N,df(1,ipos),1)
       Call dscal(N,Normi,df(1,ipos),1)
       Call dscal(N,Normi,dv(1,ipos),1)
    Endif
    Do i=1,iter_used
       Do j=i+1,iter_used
          beta(i,j)=DDOT(N,df(1, j),1,df(1,i),1)
       Enddo
       beta(i,i)= w0*w0  + 1.0_pr 
    Enddo
    Call DSYTRF('U',iter_used,beta,M,iwork,work,M,i)
    If(i.Ne.0) Stop '  In Broyden: info at DSYTRF '
    Call DSYTRI('U',iter_used,beta,M,iwork,work,i)
    If(i.Ne.0) Stop '  In Broyden: info at DSYTRI '
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
    Call dcopy(N,vin ,1,dv(1,inext),1)
    sf=+1.0_pr; Call DAXPY(N,sf,curv,1,vin,1)      
  End Subroutine qrpa_broyden_method
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Elemental Function qrpa_electric_field(t,x,y,z)
    !------------------------------------------------------------------
    ! For electric transitions operator
    ! called in qrpa_f_gamdel: qrpa_f_gamdel(t,fl(ihli),fh(ihli)) 
    ! t=0:neutrons, 1:protons
    ! Memo: fl(ihli)^2+fh(ihli)^2=x^2+y^2+z^2
    ! Note: Here in HFBTHO implementation y is the coordinate to radial direction and x = 0
    !------------------------------------------------------------------
    Implicit None
    Integer(ipr), Intent(in) :: t
    Real(pr),     Intent(in) :: x,y,z
    Real(pr)                 :: t_phase,qrpa_electric_field,qrpa_operator_strenght,qrpa_operator_field
    !
    qrpa_operator_strenght = 0.0_pr
    qrpa_operator_field = 0.0_pr
    !
    select case(L_qrpa_responce)
    case(0)
    !----------------------------------------------------    
    ! Monopole operator F_00
    !----------------------------------------------------    
      qrpa_operator_field    = (x**2+y**2+z**2) 
    !
    case(1)
    !----------------------------------------------------    
    ! Dipole operator F_1K
    !----------------------------------------------------    
      select case(K_qrpa_responce)
        case(-1)
          qrpa_operator_field    = 0.5_pr*Sqrt(1.5_pr/Pi)*y  !! here y <--> x
        case(0)
          qrpa_operator_field    = 0.5_pr*Sqrt(3.0_pr/Pi)*z
        case(1)
          qrpa_operator_field    = -0.5_pr*Sqrt(1.5_pr/Pi)*y  !! here y <--> x
        case default
      end select
    !
    case(2)
    !----------------------------------------------------    
    ! Quadrupole operator F_2K
    !----------------------------------------------------    
      select case(K_qrpa_responce)
        case(0)
          qrpa_operator_field    =  0.25_pr*Sqrt(5.0_pr/Pi)*(2.0_pr*z**2-(x**2+y**2))
        case(1)
          qrpa_operator_field    =  0.50_pr*Sqrt(7.5_pr/Pi)*z*y
        case(-1)
          qrpa_operator_field    =  0.50_pr*Sqrt(7.5_pr/Pi)*z*y
        case(2)
          qrpa_operator_field    = -0.25_pr*Sqrt(7.5_pr/Pi)*y*y
        case(-2)
          qrpa_operator_field    = -0.25_pr*Sqrt(7.5_pr/Pi)*y*y
        case default
      end select
    !----------------------------------------------------    
    ! Octupole operator F_3K
    !----------------------------------------------------    
    case(3)
       select case(K_qrpa_responce)
        case(0)
           qrpa_operator_field    = 0.25_pr*Sqrt(7.0_pr/Pi)*(2.0_pr*z**2 -3.0_pr*(x**2+y**2))*z
        case(1)
           qrpa_operator_field    = -0.125_pr*Sqrt(21.0_pr/Pi)*(4.0_pr*z**2 - y**2 )*y
        case(2)
           qrpa_operator_field    = 0.25_pr*Sqrt(105.0_pr/(2.0_pr*Pi))*(x**2+y**2)*z
        case(3)
           qrpa_operator_field    = -0.125_pr*Sqrt(35.0_pr/Pi)*y*y*y
         case(-3)
           qrpa_operator_field    = -0.125_pr*Sqrt(35.0_pr/Pi)*y*y*y
         case default
       end select
     end select
    !----------------------------------------------------    
    ! T of the operator. Note: convention to the old FAM
    ! module has been changed. The operator is no longer
    ! multiplied by charge e
    ! Note: with EFA, npr(3) is not the total A, but 
    ! npr(1)+npr(2) is correct
    !----------------------------------------------------    
    select case(T_qrpa_responce)
    case(0)
      ! isoscalar F = e (Z/A) * opr_i
      If(t.Eq.0) Then
      ! neutrons
        qrpa_operator_strenght = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Else
      ! Protons
        qrpa_operator_strenght = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Endif 
    case(1)
      ! isovector
      If(t.Eq.0) Then
      ! neutrons
        qrpa_operator_strenght = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Else
      ! Protons
        qrpa_operator_strenght = - (Dble(npr(1))/Dble(npr(1)+npr(2)))
      Endif 
    case(2)
      ! bare charge
      If(t.Eq.0) Then
      ! neutrons
        qrpa_operator_strenght = 0.0_pr
      Else
      ! Protons
        qrpa_operator_strenght = +1.0_pr
      Endif 
    case default
    end select
    !
    if( K_qrpa_responce .ne. 0 ) qrpa_operator_field = qrpa_operator_field *(1.0_pr/Sqrt(2.0_pr))  !! simplex-y
    qrpa_electric_field    = qrpa_operator_strenght*qrpa_operator_field
  End Function qrpa_electric_field
  !=================================================================================================================================
  ! 
  !=================================================================================================================================
  Elemental Function qrpa_spin_field(t,x,y,z)
    !
    ! constant external field for spin operator. Depends only on 
    ! particle type and T_qrpa_responce
    Implicit None
    Integer(ipr), Intent(in) :: t
    Real(pr),     Intent(in) :: x,y,z
    Real(pr)                 :: qrpa_spin_field
    !
    qrpa_spin_field = 0.0_pr
    select case(T_qrpa_responce)
    case(0)
      ! isoscalar type
      If(t.Eq.0) Then
      ! neutrons
        qrpa_spin_field = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Else
      ! Protons
        qrpa_spin_field = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Endif 
    case(1)
      ! isovector type
      If(t.Eq.0) Then
      ! neutrons
        qrpa_spin_field = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Else
      ! Protons
        qrpa_spin_field = - (Dble(npr(1))/Dble(npr(1)+npr(2)))
      Endif
    case(2)
      ! with bare g-factors
      If(t.Eq.0) Then
      ! neutrons
        qrpa_spin_field = -3.826_pr
      Else
      ! Protons
        qrpa_spin_field = +5.586_pr
      Endif  
    case default
    end select
    !
    qrpa_spin_field = qrpa_spin_field*Sqrt(0.75_pr/Pi)
  End Function qrpa_spin_field
  !=================================================================================================================================
  ! 
  !=================================================================================================================================
  Elemental Function qrpa_l_field(t,x,y,z)
    !
    ! constant external field for l operator. Depends only on 
    ! particle type and T_qrpa_responce
    Implicit None
    Integer(ipr), Intent(in) :: t
    Real(pr),     Intent(in) :: x,y,z
    Real(pr)                 :: qrpa_l_field
    !
    qrpa_l_field = 0.0_pr
    select case(T_qrpa_responce)
    case(0)
      ! isoscalar type
      If(t.Eq.0) Then
      ! neutrons
        qrpa_l_field = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Else
      ! Protons
        qrpa_l_field = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Endif 
    case(1)
      ! isovector type
      If(t.Eq.0) Then
      ! neutrons
        qrpa_l_field = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Else
      ! Protons
        qrpa_l_field = - (Dble(npr(1))/Dble(npr(1)+npr(2)))
      Endif
    case(2)
      ! with bare g-factors
      If(t.Eq.0) Then
      ! neutrons
        qrpa_l_field = 0.0_pr
      Else
      ! Protons
        qrpa_l_field = +1.0_pr
      Endif  
    case default
    end select
    !
    qrpa_l_field = qrpa_l_field*Sqrt(0.75_pr/Pi)*2.0_pr/(real(abs(L_qrpa_responce),kind=pr)+1.0_pr)
    !
  End Function qrpa_l_field
  !=================================================================================================================================
  ! 
  !=================================================================================================================================
  Subroutine qrpa_HFBTHO_solver
    !----------------------------------------------------------------------------------------------------
    ! qrpa solver. Input parameters:
    !   integer:   max_iter_qrpa,qrpa_nbroyden
    !   real:      qrpa_eps,qrpa_alphamix
    !   complex:   qrpa_eta,qrpa_omega
    !----------------------------------------------------------------------------------------------------
    Implicit None
    Logical      :: file_exists
    Integer(ipr) :: i,j,k1,k2,kk1,kk2,imen,ih,il,ib,ibx,nd,nza,nra,nla,nsa,iw       
    Real(pr)     :: snorm_xN,snorm_yN,snorm_N,snorm_xP,snorm_yP,snorm_P
    Real(pr)     :: xymix,DZNRM2,xynorm_N,xynorm_P
    Real(pr)     :: Rbox,eFermiN
    Real(pr)     :: time1,time2
    Real(pr)     :: strN,strP
    Character(LEN=100) :: fname
    Character(LEN=6)   :: tname,kname,mname,opname
    !----------------------------------------------------------------------------------------------------
    ! QRPA INITIALIZATION
    !----------------------------------------------------------------------------------------------------
    If(iqrpa.Eq.0) Then
       iqrpa=1
       Call qrpaINI_CUT_ALLOCATE      

       Do iw=lout,lfile
         Write(iw,*) '          < Starting FAM calculation >'
         Write(iw,*) '  FAM QRPA version: ',qrpa_Version
         Write(iw,*) '  '
         Write(iw,*) '  Time-odd coupling constants:'
         Write(iw,*) '  -----------------------------------'
         Write(iw,'("   Css(0)=  ",g26.18,"; Css(1)=  ",g26.18)') Css
         Write(iw,'("   CDss(0)= ",g26.18,"; CDss(1)= ",g26.18)') CDss
         Write(iw,'("   Cjj(0)=  ",g26.18,"; Cjj(1)=  ",g26.18)') Cjcjc
         Write(iw,'("   CsNj(0)= ",g26.18,"; CsNj(1)= ",g26.18)') Csnablajc
         Write(iw,'("   CsDs(0)= ",g26.18,"; CsDs(1)= ",g26.18)') CDeltass
         Write(iw,'("   CsT(0)=  ",g26.18,"; CsT(1)=  ",g26.18)') CsT
         Write(iw,*) '  '
         select case(QRPA_Operator_type)
         case(1)
           Write(iw,*) '  Electric transition operator used'
         case(2)
           Write(iw,*) '  Spin transition operator used'
         case(3)
           Write(iw,*) '  Magnetic transition operator used'
         case(4)
           Write(iw,*) '  Linear momentum operator used'
         case(5)
           Write(iw,*) '  Coordinate operator used'
         case(6)
           Write(iw,*) '  Total angular momentum J11 operator used'
         case(7)
           Write(iw,*) '  Total angular momentum J_y operator used'
         case default
           STOP "Unknown operator type"
         end select

         If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) Write(iw,*) '  Using normal harmonic oscillator basis.'
         if(USE_ONLY_F20) then
          Write(iw,*) '  '
          Write(iw,*) '  >>  WARNING! WARNING!  H20, H02 and all induced fields are set to zero!  WARNING! WARNING!  <<'
          Write(iw,*) '  '
         end if 

         if(USE_ONLY_F11EFA.and.(ENABLE_N_EFA.or.ENABLE_P_EFA)) then
          Write(iw,*) '  '
          Write(iw,*) '  >>  WARNING! WARNING!  Using only F11 part with EFA!  WARNING! WARNING!  <<'
          Write(iw,*) '  '
         end if 

         Write(iw,*) '  '

         if(ENABLE_N_EFA) Write(iw,*) '  Enabled neutron equal filling approximation due to q.p. blocking '
         if(ENABLE_P_EFA) Write(iw,*) '  Enabled proton equal filling approximation due to q.p. blocking '
         if(ENABLE_N_EFA.or.ENABLE_P_EFA) Write(iw,*) '  '

         Write(iw,*) '  QRPA CALCULATIONS '
         Write(iw,*) '  DIMENSION OF U,V iHFB_LEN:',iHFB_LEN
         Write(iw,*) '  DIMENSION OF Eqp iHFB_NQP:',iHFB_NQP
         Write(iw,*) '  MAX ITERATIONS           :',max_iter_qrpa   
         Write(iw,*) '  BROYDEN HISTORY          :',qrpa_nbroyden
         Write(iw,*) '  MIXING PARAMETER         :',qrpa_alphamix
         Write(iw,*) '  QRPA STOP CRITERIA       :',qrpa_eps   
       End Do
    Endif
    !----------------------------------------------------------------------------------------------------
    ! Set the file name for results
    !----------------------------------------------------------------------------------------------------
    select case(T_qrpa_responce)
    case(0)
      tname = 'is'
    case(1)
      tname = 'iv'
    case(2)
      tname = 'bare'
    case default
      STOP "Unknown T_qrpa_responce"
    end select
    !
    select case(QRPA_Operator_type)
    case(1)
      opname = 'E'
    case(2)
      opname = 'S'
    case(3)
      opname = 'M'
    case(4)
      opname = 'P'
    case(5)
      opname = 'C'
    case(6)
      opname = 'J11'
    case(7)
      opname = 'J_y'
    case default
      STOP "Unknown operator type"
    end select

    write (kname, "(I1,I1)") L_qrpa_responce,K_qrpa_responce
    write (mname, "(I6)") iam_mpi
    fname = trim('qrpa_results_' // trim(tname) // '_' //  trim(opname) //  trim(kname) // '_' // trim(adjustl(mname)) // '.dat')
    !----------------------------------------------------------------------------------------------------
    ! QRPA ITERATIONS 
    !----------------------------------------------------------------------------------------------------
    Do iter_qrpa=0,max_iter_qrpa
       Call Cpu_time(time1)                                          ! qrpa time per iteration
       If(iter_qrpa.Eq.0) Then
          Call qrpa_green_function(qrpa_omega)                      ! calculate green functions		
          hN_c=zero;  hP_c=zero;  BhN_c=zero;  BhP_c=zero           ! start with X=Y=0
          si=1.0_pr                                                 ! defect on zero iteration 
          Do iw=lout,lfile                                          ! print header information
             Write(iw,*) 
             Write(iw,*) '  OMEGA               :',qrpa_omega
             Write(iw,*) '  NEUTRON POLE(-)     :',one/GreenNminus(gfpoles_imN)
             Write(iw,*) '  PROTON  POLE(-)     :',one/GreenPminus(gfpoles_imP)
             Write(iw,*) '  T,L,K RESP.OPERATOR :',T_qrpa_responce,L_qrpa_responce,K_qrpa_responce             
             select case(QRPA_Operator_type)
             case(1)
               if(L_qrpa_responce.eq.0) then
                 Write(iw,*) '  Transition strength in units of [e^2 MeV^{-1} fm^{4}]'
               else
                 Write(iw,*) '  Transition strength in units of [e^2 MeV^{-1} fm^{2L}]'
               end if
             case(2,3)
                 if(L_qrpa_responce.eq.1) Write(iw,*) '  Transition strength in units of [(\mu_N/c)^2 MeV^{-1} ]'
             case(4)
                 if(L_qrpa_responce.eq.1) Write(iw,*) '  Transition strength in units of [ fm^{-2} MeV^{-1} ]'        !krisa! 4 is P operator
             case(5)
                 if(L_qrpa_responce.eq.1) Write(iw,*) '  Transition strength in units of [ fm^2 MeV^{-1} ]'           !krisa! 5 is Q operator
             case(6,7)
                 if(L_qrpa_responce.eq.1) Write(iw,*) '  Transition strength in units of [ MeV^{-1} ]'                !krisa! 6 is J11 and 7 is J_y
             end select
             Write(iw,'(a,f7.3,3(a,i3),(a,g14.8),(a,2g14.8),a)') &
                  '  |QRPA> iterations... b0=',b0,', Nsh=',n00,', N=',npr(1),', Z=',npr(2), &
                  ', Eta=',Real(qrpa_eta,kind=pr),', Omega=',qrpa_omega
             Write(iw,'(2x,154(''-''))')
             Write(iw,'(a,a)') '    i        si    mix        reSn         reSp         reSt         imSn         imSp         imSt', &
                               '           TrStrn       TrStrp       TrStrt          Time'
             Write(iw,'(2x,154(''-''))')
          End Do
       Endif
       Call qrpa_xy                                                        ! calculate X,Y from h,Bh
       Call qrpa_broyden(iter_qrpa)                                        ! use Broyden
       Call Calculate_QRPA_H20                                             ! calculate h,Bh from X,Y 

       Call qrpa_strength(XN_c,YN_c,P1N_c,P2N_c,FhN_c,BFhN_c,F11N_c,BF11N_c, &
                          ENABLE_N_EFA,ImStrengthN,ReStrengthN)                   ! qrpa strenght function neutrons, imag and real parts
       Call qrpa_strength(XP_c,YP_c,P1P_c,P2P_c,FhP_c,BFhP_c,F11P_c,BF11P_c, &
                          ENABLE_P_EFA,ImStrengthP,ReStrengthP)                   ! qrpa strenght function protons, imag and real parts
       !
       ! transition strength
       ReTransStrN = -ReStrengthN/pi
       ReTransStrP = -ReStrengthP/pi
       ImTransStrN = -ImStrengthN/pi
       ImTransStrP = -ImStrengthP/pi
       !		
       Call Cpu_time(time2); time2=time2-time1                             ! time per iteration 
       If(iter_qrpa.Lt.500.Or.si.Lt.qrpa_eps*10.0_pr) Then
          Do iw=lout,lfile
             Write(iw,'(i4,a,1x,f12.8,f5.2," | ",6(1x,g12.6)," | ",3(1x,g12.6)," | ",f7.2)') &
                  iter_qrpa,qrpa_bbroyden,si,qrpa_alphamix,&
                  ReStrengthN,ReStrengthP,ReStrengthN+ReStrengthP,ImStrengthN,ImStrengthP,ImStrengthN+ImStrengthP,&
                  ImTransStrN,ImTransStrP,ImTransStrN+ImTransStrP,time2     !krisa! new variables introduced, snz_qrpa(1:2) variables were removed
          End Do
       Endif
       If(si.Lt.qrpa_eps) Exit   ! Exit if converges
    End Do ! iter_qrpa 
    !----------------------------------------------------------------------------------------------------
    ! APPEND RESULTS TO qrpa_results.dat FILE 
    !----------------------------------------------------------------------------------------------------
    Inquire(FILE=fname, EXIST=file_exists)   
    Open (101,file=fname,position='append')
    If(.Not.file_exists) Write(101,'(20a)') '!  N 1 Z 2 ITER 3  ', &                                !krisa! new strength variable names introduced
         ' si 4      ReOmega 5      ImOmega 6      reStrN 7       reStrP  8      ', &
         'reStrT 9       imStrN 10      imStrP 11      imStrT 12      TrnStrN 13     TrnStrP 14     TrnStrT 15'
    Write(101,'(1x,3(1x,I3),1x,f9.6,200(f15.6))') npr(1),npr(2),iter_qrpa,si,qrpa_omega,               &
         ReStrengthN,ReStrengthP,(ReStrengthN+ReStrengthP),ImStrengthN,ImStrengthP,(ImStrengthN+ImStrengthP),    &
         ImTransStrN,ImTransStrP,(ImTransStrN+ImTransStrP)
    Close(101)
    !
    !----------------------------------------------------------------------------------------------------
    ! Write amplitudes on file
    !----------------------------------------------------------------------------------------------------
    if(qrpa_file_io.eq.1) call write_amplitudes
    !
    if(0.eq.1) then
     call print_density(ro_c(:,1),'rhon') !! neutrons rho
     call print_density(ro_c(:,2),'rhop') !! protons rho
     call print_density(sofi_c(:,1),'sfin') !! neutrons spin_phi
     call print_density(sofi_c(:,2),'sfip') !! protons spin_phi
     call print_density(rjr_c(:,1),'jrn') !! neutrons j_r
     call print_density(rjr_c(:,2),'jrp') !! protons j_r
     call print_density(rjz_c(:,1),'jzn') !! neutrons j_z
     call print_density(rjz_c(:,2),'jzp') !! protons j_z
     call print_density(aka_c(:,1),'kan') !! neutrons kap_a
     call print_density(bka_c(:,1),'kbn') !! neutrons \bar{kap}_a
     call print_density_hfb(ro(:,1)+ro(:,2),'rho') !! hfb total rho
    end if
    
  End Subroutine qrpa_HFBTHO_solver
 !=================================================================================================================================
  !
  ! Subroutine which writes the FAM-QRPA amplitudes on disk in
  ! binary format. The file content is following:
  !
  ! T_qrpa_responce,L_qrpa_responce,K_qrpa_responce,P_qrpa_response,QRPA_Operator_type,iQRPA_LEN  : 6 integers
  ! qrpa_omega                                                                                    : 1 complex
  ! ImStrengthN,ImStrengthP,ReStrengthN,ReStrengthP                                               : 4 reals
  ! XN_c(i),YN_c(i),XP_c(i),YP_c(i)                                                               : 4*iQRPA_LEN complex (i=1 to iQRPA_LEN)
  !
  !=================================================================================================================================
  Subroutine write_amplitudes
    Implicit None
    Integer(ipr) :: i
    open(20,file="amplitudes.dat",form='unformatted',status='replace',action='write')

    Write(20) T_qrpa_responce,L_qrpa_responce,K_qrpa_responce,P_qrpa_response,QRPA_Operator_type,iQRPA_LEN
    Write(20) ENABLE_N_EFA,ENABLE_P_EFA
    Write(20) qrpa_omega
    Write(20) ImStrengthN,ImStrengthP,ReStrengthN,ReStrengthP
    Do i = 1,iQRPA_LEN
      Write(20) XN_c(i),YN_c(i),XP_c(i),YP_c(i)
    End Do

!! ADD HERE ALSO EFA AMPLITUDES LATER

    close(20)
  End subroutine
  !=================================================================================================================================
  !
  ! Subroutine which reads the canonical conjugate operator Q
  ! not tested yet.
  !
  !=================================================================================================================================
  Subroutine read_operator(QN,QP)
    Implicit None
    Complex(pr), Intent(out) :: QN(iQRPA_LEN),QP(iQRPA_LEN)
    Complex(pr)              :: omegaRead,XNread,YNread,XPread,YPread
    Real(pr)                 :: ImStN,ImStP,ReStN,ReStP
    Complex(pr)              :: iunit
    Integer(ipr)             :: i,Tread,Lread,Kread,Pread,OPRread,LenRead
    logical                  :: ENABLE_N_EFAin,ENABLE_P_EFAin
    Parameter(iunit = Cmplx(0.0_pr,1.0_pr,kind=pr))

    open(20,file="amplitudes.dat",form='unformatted',status='old',action='read')

    Read(20) Tread,Lread,Kread,Pread,OPRread,LenRead
    Read(20) ENABLE_N_EFAin,ENABLE_P_EFAin
    Read(20) omegaRead
    if(ENABLE_N_EFAin.or.ENABLE_P_EFAin) STOP "EFA NOT YET SUPPORTED WITH READ AMPLITUDES"
    if(Kread.ne.K_qrpa_responce) stop "Wrong K in the read operator!"
    if(Pread.ne.P_qrpa_response) stop "Wrong parity in the read operator!"
    if(LenRead.ne.iQRPA_LEN)     stop "Incorrect length in read operator!"
    if(abs(omegaRead).gt.spacing(1.0_pr)) stop "Solution on the file has omega > 0!"

    Read(20) ImStN,ImStP,ReStN,ReStP
    Do i = 1,iQRPA_LEN
      Read(20) XNread,YNread,XPread,YPread
      QN(i) = -0.5_pr*iunit*(XNread+conjg(YNread))/cmplx(ReStN+ReStP,ImStN+ImStP,kind=pr)
      QP(i) = -0.5_pr*iunit*(XPread+conjg(YPread))/cmplx(ReStN+ReStP,ImStN+ImStP,kind=pr)
    End Do

    close(20)
  End subroutine
  !=================================================================================================================================
  !  
  ! Helper subroutine to print qrpa densitiens to file
  !
  !=================================================================================================================================
  Subroutine print_density(densi,fnamepart)
    Implicit None
    Complex(pr), Intent(in)       :: densi(nghl)
    Character(LEN=*), intent(in) :: fnamepart
    Integer(ipr)      :: ihli
    integer(ipr)      :: outformat
    Character(LEN=100) :: fnamer,fnamei
    !
    outformat = 1 !! Table
    fnamer = trim('qrpa_density_' // trim(fnamepart) // '_re.dat')
    fnamei = trim('qrpa_density_' // trim(fnamepart) // '_im.dat')

    Open (123,file=fnamer,status='replace')
    Write(123,*) "! Real part of density"

    if(outformat.eq.1) write(123,*) "! list of r_p, z, Re(density)"
    Do ihli=1,nghl
      if(outformat.eq.1) write(123,'(f11.6, " ", f11.6, " ", f18.12)') fl(ihli),fh(ihli),Real(densi(ihli))
    End do
    Close(123)

    Open (123,file=fnamei,status='replace')
    Write(123,*) "! Imaginary part of density"

    if(outformat.eq.1) write(123,*) "! list of r_p, z, Im(density)"
    Do ihli=1,nghl
      if(outformat.eq.1) write(123,'(f11.6, " ", f11.6, " ", f18.12)') fl(ihli),fh(ihli),Imag(densi(ihli))
    End do
    Close(123)
    !
  End Subroutine print_density
  !=================================================================================================================================
  !  
  ! Helper subroutine to print hfb densitiens to file
  !
  !=================================================================================================================================
  Subroutine print_density_hfb(densi,fnamepart)
    Implicit None
    real(pr), Intent(in)          :: densi(nghl)
    Character(LEN=*), intent(in) :: fnamepart
    Integer(ipr)      :: ihli
    integer(ipr)      :: outformat
    Character(LEN=100) :: fname
    !
    outformat = 1 !! Table
    fname = trim('hfb_density_' // trim(fnamepart) // '_re.dat')

    Open (123,file=fname,status='replace')
    Write(123,*) "! Real part of density"

    if(outformat.eq.1) write(123,*) "! list of r_p, z, density"
    Do ihli=1,nghl
      if(outformat.eq.1) write(123,'(f11.6, " ", f11.6, " ", f18.12)') fl(ihli),fh(ihli),densi(ihli)
    End do
    Close(123)
  End Subroutine print_density_hfb
  !=================================================================================================================================
  !
  ! Helper routines for QRPA matrix handling
  !
  !=================================================================================================================================
  function qrpa_isHermitian(A,epsin)
    !-----------------------------------------------------------------------
    ! Returns true if A is Hermitian. Tolerance set by spacing, if not
    ! given. Used for debug. Note: If A=zero, returns true
    !-----------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)        :: A(:)
    Real(pr), Intent(in), optional :: epsin
    Complex(pr), Allocatable    :: temp(:)
    Real(pr) :: eps
    logical  :: qrpa_isHermitian
    integer  :: i
 
    eps = spacing(1.0_pr)
    if(present(epsin)) eps = epsin
    qrpa_isHermitian = .true.
    Allocate(temp(iQRPA_LEN))
    call qrpa_GiveTransposed(A,temp)
    temp = (A - conjg(temp))
    Do i=1,iQRPA_LEN
      if( abs(temp(i)).gt.eps ) then
         qrpa_isHermitian = .false.
         return
      end if
    End do
    Deallocate(temp)
  End function qrpa_isHermitian
  !=================================================================================================================================
  !
  !=================================================================================================================================
  function qrpa_isSymmetric(A,epsin)
    !-----------------------------------------------------------------------
    ! Returns true if A is symmetric. Tolerance set by spacing, if not
    ! given. Used for debug. Note: If A=zero, returns true
    !-----------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)        :: A(:)
    Real(pr), Intent(in), optional :: epsin
    Complex(pr), Allocatable    :: temp(:)
    Real(pr) :: eps
    logical  :: qrpa_isSymmetric
    integer  :: i
 
    eps = spacing(1.0_pr)
    if(present(epsin)) eps = epsin
    qrpa_isSymmetric = .true.
    Allocate(temp(iQRPA_LEN))
    call qrpa_GiveTransposed(A,temp)
    temp = (A - temp)
    Do i=1,iQRPA_LEN
      if( abs(temp(i)).gt.eps ) then
         qrpa_isSymmetric = .false.
         Deallocate(temp)
         return
      end if
    End do
    Deallocate(temp)
  End function qrpa_isSymmetric
  !=================================================================================================================================
  !
  !=================================================================================================================================
  function qrpa_isAntiSymmetric(A,epsin)
    !-----------------------------------------------------------------------
    ! Returns true if A is antisymmetric. Tolerance set by spacing, if not
    ! given. Used for debug. Note: If A=zero, returns true
    !-----------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)        :: A(:)
    Real(pr), Intent(in), optional :: epsin
    Complex(pr), Allocatable    :: temp(:)
    Real(pr) :: eps
    logical  :: qrpa_isAntiSymmetric
    integer  :: i
 
    eps = spacing(1.0_pr)
    if(present(epsin)) eps = epsin
    qrpa_isAntiSymmetric = .true.
    Allocate(temp(iQRPA_LEN))
    call qrpa_GiveTransposed(A,temp)
    temp = (A + temp)
    Do i=1,iQRPA_LEN
      if( abs(temp(i)).gt.eps ) then
         qrpa_isAntiSymmetric = .false.
         Deallocate(temp)
         return
      end if
    End do
    Deallocate(temp)
  End function qrpa_isAntiSymmetric
  !=================================================================================================================================
  !
  !=================================================================================================================================
  function qrpa_isReal(A,epsin)
    !-----------------------------------------------------------------------
    ! Returns true if A is real. Tolerance set by spacing, if not
    ! given. Used for debug. Note: If A=zero, returns true
    !-----------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)        :: A(:)
    Real(pr), Intent(in), optional :: epsin
    Complex(pr), Allocatable    :: temp(:)
    Real(pr) :: eps
    logical  :: qrpa_isReal
    integer  :: i
 
    eps = spacing(1.0_pr)
    if(present(epsin)) eps = epsin
    qrpa_isReal = .true.
    Allocate(temp(iQRPA_LEN))
    temp = Imag(A)
    Do i=1,iQRPA_LEN
      if( abs(temp(i)).gt.eps ) then
         qrpa_isReal = .false.
         Deallocate(temp)
         return
      end if
    End do
    Deallocate(temp)
  End function qrpa_isReal
  !=================================================================================================================================
  !
  !=================================================================================================================================
  function qrpa_isImag(A,epsin)
    !-----------------------------------------------------------------------
    ! Returns true if A is purely imaginary. Tolerance set by spacing, if not
    ! given. Used for debug. Note: If A=zero, returns true
    !-----------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)        :: A(:)
    Real(pr), Intent(in), optional :: epsin
    Complex(pr), Allocatable    :: temp(:)
    Real(pr) :: eps
    logical  :: qrpa_isImag
    integer  :: i
 
    eps = spacing(1.0_pr)
    if(present(epsin)) eps = epsin
    qrpa_isImag = .true.
    Allocate(temp(iQRPA_LEN))
    temp = Real(A)
    Do i=1,iQRPA_LEN
      if( abs(temp(i)).gt.eps ) then
         qrpa_isImag = .false.
         Deallocate(temp)
         return
      end if
    End do
    Deallocate(temp)
  End function qrpa_isImag
  !=================================================================================================================================
  !
  !=================================================================================================================================
  function qrpa_isZero(A,epsin)
    !-----------------------------------------------------------------------
    ! Returns true if A is real. Tolerance set by spacing, if not
    ! given. Used for debug
    !-----------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)        :: A(:)
    Real(pr), Intent(in), optional :: epsin
    Complex(pr), Allocatable    :: temp(:)
    Real(pr) :: eps
    logical  :: qrpa_isZero
    integer  :: i
 
    eps = spacing(1.0_pr)
    if(present(epsin)) eps = epsin
    qrpa_isZero = .true.
    Do i=1,iQRPA_LEN
      if( abs(A(i)).gt.eps ) then
        qrpa_isZero = .false.
        return
      end if
    End do
  End function qrpa_isZero
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_printProp(s,A)
    !
    ! prints qrpa matrix A properties, for debug. Put variable name, or such, for s
    !
    Implicit None
    Character(LEN=*), intent(in) :: s
    Complex(pr), Intent(in)      :: A(:)
    write(*,*) trim(s)," is Hermitian: ", qrpa_isHermitian(A,0.0001_pr), &
                       " symmetric: "   , qrpa_isSymmetric(A,0.0001_pr), &
                       " antisymmetric: ",qrpa_isAntiSymmetric(A,0.0001_pr), &
                       " real: "        , qrpa_isReal(A,0.0001_pr), &
                       " imaginary: "   , qrpa_isImag(A,0.0001_pr), &
                       " zero: "        , qrpa_isZero(A,0.0001_pr)
  end subroutine qrpa_printProp
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_SYM(A)
    !-----------------------------------------------------------------------
    ! symmetrize QRPA matrix so that A^T = A. This routine is only for debug
    !-----------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(inout)  :: A(:)
    Complex(pr), Allocatable    :: temp(:)

    Allocate(temp(iQRPA_LEN))
    call qrpa_GiveTransposed(A,temp)
    A = half * (A + temp)    
    Deallocate(temp)
  End Subroutine qrpa_SYM
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_ADJ(A)
    !-------------------------------------------------------------------------
    ! Hermitice QRPA matrix so that A^\dag = A. This routine is only for debug
    !-------------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(inout)  :: A(:)
    Complex(pr), Allocatable    :: temp(:)

    Allocate(temp(iQRPA_LEN))
    call qrpa_GiveTransposed(A,temp)
    A = half * (A + conjg(temp))
    Deallocate(temp)
  End Subroutine qrpa_ADJ
  !  
End Module qrpa_HFBTHO
!===================================================================================================================================
!#END qrpa_HFBTHO MODULE
!===================================================================================================================================
#endif
