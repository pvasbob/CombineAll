Module RBM
   ! precision
   Implicit None
   Integer, Parameter, Public :: ipr=Kind(1)     
   Integer, Parameter, Public :: pr =Kind(1.0d0) 
   Complex(pr), Parameter :: iunit = (0.0_pr, 1.0_pr)
   Real(pr) :: pi, Norm

   Integer :: nop ! number of the operators to be computed
   Integer :: ntrain ! number of the training points read from file

   Integer :: nrbm ! number of the training points used in RBM


   ! FAM training data read from binary file
   Integer :: nuv ! dimension of the two-quasiparticle states (including both neutrons and protons in the like-particle case)
   Complex(pr), Dimension(:),   Allocatable :: twoEqp     ! two-quasipartciel energies
   Complex(pr), Dimension(:),   Allocatable :: omegatrain ! energies of the training points
   Complex(pr), Dimension(:,:), Allocatable :: Xtrain, Ytrain, dH20, dH02 ! X, Y, deltaH20, deltaH02 at the training points
   Complex(pr), Dimension(:,:), Allocatable :: F20, F02 ! one-body operators whose strengths to be computed

   ! arrays that depends on ntrain
   Complex(pr), Dimension(:),   Allocatable :: SFtrain, TFtrain, norm_eigen
   Complex(pr), Dimension(:,:), Allocatable :: NormKernel, HamiltonianKernel, NormKernelHalf, NormKernelRegularized, &
   & unitmat, NormKernelHalfInv, u_norm, u_norminv
   Real(pr), Allocatable :: realtemp(:)

   Logical :: SMALLESTFIRST
   Integer ::  colldim, ierr
    Integer, Allocatable :: collidx(:), SortedOrder(:), SortedOrder_Hcoll(:)
   
   ! arrays that depends on colldim
   Complex(pr), Dimension(:),   Allocatable :: RBMenergy(:), sqrt_norm(:), RBMstrength(:), RBMstrength1(:), RBMstrength2(:)
   Complex(pr), Dimension(:,:), Allocatable ::  H_coll(:,:), g_coll(:,:),g_collinv(:,:), ugn(:,:), ugn2(:,:), X_QRPA_RBM(:,:), Y_QRPA_RBM(:,:)
   Real(pr), Allocatable :: RBMstrengthfromXY(:,:)

   ! emulator
   Complex(pr) :: qrpa_omega, SFomega 
   Real(pr) :: sumrule_emulator(-4:4)

   !---------------------------------------------------------------------
   ! FAM Namelists
   !---------------------------------------------------------------------
   Real(pr) :: normcut
   Integer(ipr) :: number_of_training_points, number_of_operators, nmax_emulatorrun, VARIATION
   Complex(pr) :: qrpa_omega_emulatorrun_start, qrpa_omega_emulatorrun_step
   Logical :: STRENGTH_EMULATOR_RUN, mirror_points(1:4)

   Character(50) :: fam_training_inputfile,rbm_outputfile,emulator_outputfile

   Namelist /FAM_RBM_INPUT/ normcut, number_of_training_points, number_of_operators, fam_training_inputfile, &
   & rbm_outputfile, emulator_outputfile, VARIATION, STRENGTH_EMULATOR_RUN, qrpa_omega_emulatorrun_start, qrpa_omega_emulatorrun_step, nmax_emulatorrun, mirror_points
   !

Contains
   !=======================================================================
   ! Read NAMELIST from rbm_NAMELIST.dat
   !=======================================================================
   Subroutine read_rbmfam_NAMELIST  
      Implicit None
      Integer(ipr) :: i,ios,lnamelist=16
      !------------------------------------
      ! Namelist (default values)
      !------------------------------------   
      normcut=1.0d-14
      number_of_training_points=0
      number_of_operators=0
      fam_training_inputfile="temp"
      rbm_outputfile="temp"
      emulator_outputfile="temp"
      VARIATION=1
      STRENGTH_EMULATOR_RUN=.True.
      qrpa_omega_emulatorrun_start=0.0_pr
      qrpa_omega_emulatorrun_step=0.1_pr
      nmax_emulatorrun=1
      mirror_points(1:4)=.True.      
      !------------------------------------
      ! Namelist (handling)
      !------------------------------------   
      Open(unit=lnamelist,file='rbmfam_NAMELIST.dat',delim='apostrophe')
      !
      Read(unit=lnamelist,nml=FAM_RBM_INPUT, iostat=ios)
      If (ios.Ne.0) Print *, "Namelist FAM_RBM_INPUT error, ios = ", ios
      Close(lnamelist)     
      !
      nop=number_of_operators
      ntrain=number_of_training_points
      ! dimension of the RBM
      nrbm=0
      Do i = 1, 4
         If(mirror_points(i)) nrbm =nrbm + ntrain
      End Do
!      Print *, "ntrain = ", ntrain, "nrbm = ", nrbm
      Return
   End Subroutine read_rbmfam_NAMELIST
   !----------------------------------------------------------------------------------------------------
   ! Read FAM data
   !----------------------------------------------------------------------------------------------------
   Subroutine read_fam_training
      Implicit None
      Integer :: i, trainingdata=15
      Integer :: qqi, qqj
      Character(50) :: TempInput
      Logical :: file_exists

      Do qqj = 0, number_of_training_points-1 
         Write(TempInput, '(A,I3.3)') trim(fam_training_inputfile), qqj
         Inquire(File = TempInput, exist = file_exists)
         If(file_exists) then 
            Write(*,*) trim(TempInput), "existsqq"
            !Open(unit=trainingdata,file=fam_training_inputfile,action='read',form='unformatted')
            !Open(unit=trainingdata,file=fam_training_inputfile,action='read')
             Open(unit=trainingdata,file=TempInput,action='read')
             Print *, "**** reading FAM training data ****"
              
             If(qqj .eq. 0) Then                  

               Read(trainingdata, *) nuv
               write(*,*) 'qq, nuv', nuv
               ! Allocations
               Allocate( twoEqp(nuv), F20(nuv,nop), F02(nuv,nop))      
               Allocate( omegatrain(nrbm), Xtrain(nuv,nrbm), Ytrain(nuv,nrbm), dH20(nuv,nrbm), dH02(nuv,nrbm))
               twoEqp = 0.0_pr; F20 = 0.0_pr; F02 = 0.0_pr
               omegatrain = 0.0_pr; Xtrain = 0.0_pr; Ytrain = 0.0_pr; dH20 = 0.0_pr; dH02 = 0.0_pr
      
               Read(trainingdata, *) twoEqp(:)
               write(*,*) 'twoEqp: ', size(twoEqp) 
      
               do qqi = 1, size(twoEqp)
                  write(*,*) qqi, twoEqp(qqi)
               end do
      
               write(*,*) 'twoEqp: ', twoEqp
               Do i = 1, nop
                  Read(trainingdata, *) F20(:,i), F02(:,i)
                  write(*,*) 'F20: ', size(F20)
                  write(*,*) 'F20: ', F20
                  Print *, "operator i = ", i, "read"
               End Do
            End if ! 
            !  Do i = 1, ntrain
                Read(trainingdata, *) omegatrain(qqj+1), Xtrain(:,qqj+1), Ytrain(:,qqj+1), dH20(:,qqj+1), dH02(:, qqj+1)
               !  Print *, "training i  = ", i, "read.  omega = ", omegatrain(i)
            !  End Do
   
             Close(trainingdata)
         Endif
      Enddo

      Return
   End Subroutine read_fam_training
   !----------------------------------------------------------------------------------------------------
   ! Diagonalize a general complex matrix
   !----------------------------------------------------------------------------------------------------
   Subroutine DiagGenComplexMat(dim, mat, pmat, eigenvalue, pmatinv, info)   
      ! mat = pmat * eigenvalue * pmatinv
      Implicit None
      Integer, Intent(in) :: dim
      Complex(pr), Intent(in) :: mat(dim,dim)
      Complex(pr), Intent(out) :: pmat(dim,dim), eigenvalue(dim), pmatinv(dim,dim)
      Integer, Intent(out) :: info
   
      !zgesvd/zgesv variables
      Character :: jobvl, jobvr
      Integer :: m, n, lda, ldb, nrhs, ldvl, ldvr, lwork
      Integer, Allocatable :: ipiv(:)
      Complex(pr), Allocatable :: a(:,:), w(:), vl(:,:), vr(:,:), work(:), b(:,:)
      Real(pr), Allocatable :: rwork(:)
   
      Complex(pr) :: tempmat(dim,dim)
      Integer :: i, j, k
       
      jobvl = 'n'
      jobvr = 'v'
      m = dim
      n = dim
      lda = dim
      ldvl = dim
      ldvr = dim
      lwork = 4*n+10 ! >= 2*N
      Allocate(a(lda,n), w(n),vl(ldvl,n),vr(ldvr,n), work(lwork),rwork(2*n))
      a(:,:) = mat(:,:)
      w = 0.0_pr; vl = 0.0_pr; vr = 0.0_pr; work = 0.0_pr; rwork = 0.0_pr  
!      print *, "zgeev"
      Call zgeev(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
!      print *, "zgeev done"
      If (info .ne. 0) Then
         Print *,   "zgeev error: info = ", info
         Write(0,*) "zgeev error: info = ", info
         Stop "DiagGenComplexMat: diagonalizatioin failed"
      End If
   
      pmat(:,:) = vr(:,:) 
      eigenvalue(:) = w(:)
   
      ldb = dim
      nrhs = dim
      Allocate(b(ldb,nrhs), ipiv(n))
      b = 0.0_pr
      Do i = 1, dim
         b(i,i) = 1.0_pr
      End Do
      a(:,:) = vr(:,:)
      ipiv = 0
!      Print *, "zgesv"
      Call zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
!      Print *, "zgesv done"
      If (info .ne. 0) Then
         Print *,   "zgesv error: info = ", info
         Write(0,*) "zgesv error: info = ", info
         Stop "DiagGenComplexMat: linear equation failed"
      End If
      pmatinv(:,:) = b(:,:)
   
      Deallocate(a,w,vl,vr,work,rwork,b,ipiv)
   
      tempmat = 0.0d0
      Do i = 1, dim
         Do j = 1, dim
            Do k = 1, dim 
               tempmat(i,j) = tempmat(i,j) + pmat(i,k) * eigenvalue(k) * pmatinv(k,j)
            End Do
         End Do
      End Do
   
      Print *, "diff: A - PEP^{-1}", Maxval( Abs( tempmat(:,:) - mat(:,:)))
   
      Return
   
   End Subroutine DiagGenComplexMat
   !----------------------------------------------------------------------------------------------------
   ! Sort vector
   !----------------------------------------------------------------------------------------------------
   Subroutine SORT2(dim, vector, SortedOrder, SMALLESTFIRST)
      Implicit None
  
      Integer,                    Intent(in)  :: dim
      Real(pr), Dimension(1:dim), Intent(in)  :: vector
      Integer,  Dimension(1:dim), Intent(out) :: SortedOrder
      Logical,                    Intent(in)  :: SMALLESTFIRST
      
      Integer :: i,SortedOrderTemp
      Logical :: COMPLETED
      
      ! Initialization of SortedOrder Array                                                                                 
      Do i = 1, dim
         SortedOrder(i) = i
      End Do
      
      Sort_loop :Do         
         COMPLETED = .True.
         Do i = 1, dim -1

            If(SMALLESTFIRST) Then               
               If( Vector(SortedOrder(i)) > Vector(SortedOrder(i+1)) ) Then
                  SortedOrderTemp  = SortedOrder(i)
                  SortedOrder(i)   = SortedOrder(i+1)
                  SortedOrder(i+1) = SortedOrderTemp
                  COMPLETED = .False.
               End If
            Else
               If( Vector(SortedOrder(i)) < Vector(SortedOrder(i+1)) ) Then
                  SortedOrderTemp  = SortedOrder(i)
                  SortedOrder(i)   = SortedOrder(i+1)
                  SortedOrder(i+1) = SortedOrderTemp
                  COMPLETED = .False.
               End If               
            End If            
         End Do

         If (COMPLETED) Exit Sort_loop
                  
      End Do Sort_loop
      
      ! check                                                                                                                                                     
      Do i = 1, dim
         If(SortedOrder(i) > dim .OR. SortedOrder(i) < 0 ) Stop "SORT2 Error" 
      End Do
      Return      
   End Subroutine SORT2
   !----------------------------------------------------------------------------------------------------
   ! FAM strength emulator
   !----------------------------------------------------------------------------------------------------
   Subroutine FAMEmulator(dim, omega_c, energy_c, strengthsquared, SFout)
      Implicit None
      Integer, Intent(in) :: dim
      Complex(pr), Intent(in) :: omega_c
      Complex(pr), Intent(in) :: energy_c(dim), strengthsquared(dim)
      Complex(pr), Intent(out) :: SFout      
      Integer :: i

      SFout = 0.0_pr
      Do i = 1, dim
         SFout = SFout - strengthsquared(i)/ (energy_c(i) - omega_c)
      End Do

      Return
   End Subroutine FamEmulator
End Module RBM



Program RBM_FAM
   Use RBM
   Implicit None

   Integer :: i, j, k, l, i_rbm, outp=20, emu=21
   Complex(pr), Dimension(:,:), Allocatable :: tempmat, tempmat2, tempmat3
   pi = Acos(-1.0_pr)

   !----------------------------------------------------------------------------------------------------
   ! Read Namelist
   !----------------------------------------------------------------------------------------------------
   Call read_rbmfam_NAMELIST

   Print *, "**** Parameters on training inputfile ****"
   Print *, "fam_training_inputfile       = ", fam_training_inputfile
   Print *, "number_of_training_points    = ", number_of_training_points
   Print *, "number_of_operators          = ", number_of_operators
   Write(*,*)
   Print *, "**** Parameters on RBM calculation ****"
   Print *, "mirror_points                = ", mirror_points
   Print *, "VARIATION                    = ", VARIATION
   Print *, "normcut                      = ", normcut
   Print *, "rbm_outputfile               = ", rbm_outputfile
   Write(*,*)
   Print *, "**** Parameters on FAM Emulator Run ****"
   Print *, "STRENGTH_EMULATOR_RUN        = ", STRENGTH_EMULATOR_RUN
   Print *, "qrpa_omega_emulatorrun_start = ", qrpa_omega_emulatorrun_start
   Print *, "qrpa_omega_emulatrorun_step  = ", qrpa_omega_emulatorrun_step
   Print *, "nmax_emulatorrun             = ", nmax_emulatorrun
   Print *, "emulator_outputfile          = ", emulator_outputfile
   !

   !----------------------------------------------------------------------------------------------------
   ! Read FAM training Data from file
   !----------------------------------------------------------------------------------------------------
   Call read_fam_training

   !----------------------------------------------------------------------------------------------------
   ! Complete data
   !----------------------------------------------------------------------------------------------------
   !  dH^{20}_{\mu\nu} + X_{\mu\nu}(omega)*(E_\mu + E_\nu)
   !  dH^{02}_{\mu\nu} + Y_{\mu\nu}(omega)*(E_\mu + E_\nu)
   Do i = 1, ntrain
      Do k = 1, nuv
         dH20(k,i) = dH20(k,i) + Xtrain(k,i)*twoEqp(k) 
         dH02(k,i) = dH02(k,i) + Ytrain(k,i)*twoEqp(k)
      End Do
   End Do

   i_rbm = 0
   If(mirror_points(1)) Then ! training is assumed to be performed in the 1st quadrant
      Do i = 1, ntrain
         i_rbm = i_rbm + 1
      End Do
   End IF
   If(mirror_points(2)) Then ! includes points in the 2nd quadrant
      Do i = 1, ntrain
         i_rbm = i_rbm + 1
         Xtrain(:,i_rbm) = Ytrain(:,i)
         Ytrain(:,i_rbm) = Xtrain(:,i)
         dH20(:,i_rbm)   = dH02(:,i)
         dH02(:,i_rbm)   = dH20(:,i)
         omegatrain(i_rbm) = -Dble(omegatrain(i)) + iunit*Aimag(omegatrain(i))
      End Do      
   End If
   If(mirror_points(4)) Then! includes points in the 4th quadrant
      Do i = 1, ntrain
         i_rbm = i_rbm + 1
         Xtrain(:,i_rbm) = Conjg(Xtrain(:,i))
         Ytrain(:,i_rbm) = Conjg(Ytrain(:,i))
         dH20(:,i_rbm)   = Conjg(dH20(:,i))
         dH02(:,i_rbm)   = Conjg(dH02(:,i))
         omegatrain(i_rbm) = Conjg(omegatrain(i)) 
      End Do      
   End If
   If(mirror_points(3)) Then ! includes points in the 3rd quadrant
      Do i = 1, ntrain
         i_rbm = i_rbm + 1
         Xtrain(:,i_rbm) = -Conjg(Ytrain(:,i))
         Ytrain(:,i_rbm) = -Conjg(Xtrain(:,i))
         dH20(:,i_rbm)   = -Conjg(dH02(:,i))
         dH02(:,i_rbm)   = -Conjg(dH20(:,i))
         omegatrain(i_rbm) = -omegatrain(i) 
      End Do      
   End If

   If(i_rbm .ne. nrbm) Stop "dimension error"

   !----------------------------------------------------------------------------------------------------
   ! Strength function at the training 
   !----------------------------------------------------------------------------------------------------
   Allocate(SFtrain(nrbm), TFtrain(nrbm))
   Do i = 1, nrbm
      SFtrain(i) = Dot_Product(      F20(:,1), Xtrain(:,i)) + Dot_Product(      F02(:,1),  Ytrain(:,i)) 
      TFtrain(i) = Dot_Product(Conjg(F20(:,1)),Xtrain(:,i)) + Dot_Product(Conjg(F02(:,1)), Ytrain(:,i))
   End Do
   Print *, "Strength functions at the training energies (SF, TF)"
   Write(*,*) "i    SF(Re, Im),   TF(Re, Im) omega(Re, Im)"
   Do i = 1, nrbm
      Write(*,'(I5,6ES20.10)') i, SFtrain(i), TFtrain(i), omegatrain(i)
   End Do
   !----------------------------------------------------------------------------------------------------
   ! Allocation of matrices
   !----------------------------------------------------------------------------------------------------
   Allocate(NormKernel(nrbm,nrbm), HamiltonianKernel(nrbm,nrbm), u_norm(nrbm,nrbm), norm_eigen(nrbm), u_norminv(nrbm,nrbm))
   Allocate(collidx(nrbm),sqrt_norm(nrbm))

   Allocate(NormKernelHalf(nrbm,nrbm), NormKernelHalfInv(nrbm,nrbm))
   Allocate(tempmat(nrbm,nrbm), NormKernelRegularized(nrbm,nrbm), unitmat(nrbm,nrbm), tempmat2(nrbm,nrbm), tempmat3(nrbm,nrbm), realtemp(nrbm))
   NormKernel = 0.0d0;  HamiltonianKernel = 0.0d0;   NormKernelHalf = 0.0d0; NormKernelHalfInv = 0.0d0
   u_norm = 0.0d0; norm_eigen = 0.0d0; collidx = 0; sqrt_norm = 0.0d0; u_norminv = 0.0_pr
   tempmat = 0.0d0; NormKernelRegularized = 0.0d0
   unitmat = 0.0d0; tempmat2 = 0.0d0; tempmat3 = 0.0d0

   ! unit matrix
   unitmat = 0.0d0
   Do i = 1, nrbm
      unitmat(i,i) = 1.0_pr
   End Do

   !----------------------------------------------------------------------------------------------------
   ! Kernel Calculation
   !----------------------------------------------------------------------------------------------------
   If( VARIATION == 1 ) Then ! variation with complex conjugate (standard)
      Do i = 1, nrbm
      DO j = 1, nrbm
         ! Norm Kernel
         NormKernel(i,j) = Dot_Product(Xtrain(:,i),Xtrain(:,j)) - Dot_Product(Ytrain(:,i),Ytrain(:,j)) 
         ! Hamiltonian Kernel
         HamiltonianKernel(i,j) = Dot_Product(Xtrain(:,i),dH20(:,j)) + Dot_Product(Ytrain(:,i),dH02(:,j))
      End Do
      End Do
   Else If( VARIATION == 2 ) Then ! variation without complex conjugate
      Do j = 1, nrbm
      Do i = 1, nrbm
         ! Norm Kernel
         NormKernel(i,j) = Dot_Product(Conjg(Xtrain(:,i)),Xtrain(:,j)) - Dot_Product(Conjg(Ytrain(:,i)),Ytrain(:,j)) 
         ! Hamiltonian Kernel'
         HamiltonianKernel(i,j) = Dot_Product(Conjg(Xtrain(:,i)),dH20(:,j)) + Dot_Product(Conjg(Ytrain(:,i)),dH02(:,j)) 
      End Do
      End Do
   Else 
      Stop "VARIATION should be 1 or 2"
   End If

   !----------------------------------------------------------------------------------------------------
   ! Remove tiny numerical error from Norm Kernel
   !----------------------------------------------------------------------------------------------------   
   Do i = 1, nrbm
      Do j = 1, nrbm
         If( Abs(Dble(NormKernel(i,j)))  < 1.0d-15) NormKernel(i,j) = iunit * Aimag(NormKernel(i,j))
         If( Abs(Aimag(NormKernel(i,j))) < 1.0d-15) NormKernel(i,j) = Dble(NormKernel(i,j))
      End Do
   End Do

   !----------------------------------------------------------------------------------------------------
   ! Diagonalization of Norm Kernel
   !----------------------------------------------------------------------------------------------------      
   Print *, "**** Diagonalizing Norm Kernel ****"
   Call DiagGenComplexMat(nrbm, NormKernel, u_norm, norm_eigen, u_norminv, ierr)
   !
   If (ierr .ne. 0) Then
      Print *,   "ZGEEV ERROR: INFO = ", ierr
      Write(0,*) "ZGEEV ERROR: INFO = ", ierr
      Stop "norm kernel diagonalizatioin failed"
    End If
   !----------------------------------------------------------------------------------------------------
   ! Sort norm eigenvalues
   !----------------------------------------------------------------------------------------------------      
   SMALLESTFIRST = .false.
   If(Allocated(SortedOrder)) Deallocate(SortedOrder)
   Allocate(SortedOrder(1:nrbm))
   realtemp(1:nrbm) = Dble(norm_eigen(1:nrbm))
   Call SORT2(nrbm, realtemp(1:nrbm), SortedOrder(1:nrbm), SMALLESTFIRST)
   
   If(VARIATION == 1) Then 
      !norm kernel should be Hermitian. Remove small imaginary part
      Do i = 1, nrbm
         If( Abs(Aimag(norm_eigen(i)))<1.0d-15) norm_eigen(i) = Dble(norm_eigen(i))
      End Do
   End If

   Print *, "Norm eigenvalues (Re, Im)"
   Do i = 1, nrbm
      Print *, i, Dble(norm_eigen(SortedOrder(i))), Aimag(norm_eigen(SortedOrder(i)))
   End Do
   Print *, "----------------------------"
   Print *, "normcut = ", normcut

   !----------------------------------------------------------------------------------------------------
   ! Norm eigenvalue cutoff
   !----------------------------------------------------------------------------------------------------      
   j = 0
   Print *, "   i      Sqrt(norm) (Re, Im) "
   Do i = 1, nrbm
      If( Abs(norm_eigen(SortedOrder(i)))< normcut) Then
         Cycle
      Else 
         j = j + 1
         collidx(j) = SortedOrder(i)
      End If

      sqrt_norm(j) = norm_eigen(SortedOrder(i))**0.5_pr
      If(VARIATION==1 .and. Dble(norm_eigen(SortedOrder(i)))<0.0_pr .and. Abs(Dble(sqrt_norm(j)))<1.0d-15) sqrt_norm(j)=Aimag(sqrt_norm(j))*iunit

      Write(*,'(I5,2ES20.10)') i, sqrt_norm(j)
   End Do
   colldim = j
   print *, "colldim = ", colldim  ! dimension excluding eigenvalues below cutoff
   If( colldim == 0 ) Stop "colldim = 0"

   !----------------------------------------------------------------------------------------------------
   ! Norm kernel excluding small eigenvalues
   !----------------------------------------------------------------------------------------------------      
   NormKernelRegularized = 0.0d0
   Do j = 1, nrbm
      Do i = 1, nrbm
         Do k = 1, colldim
               NormKernelRegularized(i,j) = NormKernelRegularized(i,j) + u_norm(i,collidx(k)) * norm_eigen(collidx(k)) * u_norminv(collidx(k),j)
         End Do
      End Do
   End Do

   ! check
   tempmat = 0.0d0
   Do j = 1, nrbm
      Do i = 1, nrbm
         Do k = 1, nrbm
               tempmat(i,j) = tempmat(i,j) + u_norm(i,k) * norm_eigen(k) * u_norminv(k,j)
         End Do
      End Do
   End Do
   Print *, "norm kernel diagonalization check"
   Print *, "max diff N - u n u^{-1} : ", maxval(Abs( NormKernel(:,:) - tempmat(:,:)))
   Print *, "max diff N - Nreg       : ", maxval(Abs( NormKernel(:,:) - NormKernelRegularized(:,:)))
   Print *, "norm kernel diagonalization check completed"


   !----------------------------------------------------------------------------------------------------
   ! N^{1/2}
   !----------------------------------------------------------------------------------------------------      
   NormKernelHalf = 0.0d0; NormKernelHalfInv = 0.0d0
   Do i = 1, nrbm
      Do j = 1, nrbm
         Do k = 1, colldim
            NormKernelHalf(i,j)    = NormKernelHalf(i,j)    + u_norm(i,collidx(k)) * sqrt_norm(k)       * u_norminv(collidx(k),j) 
            NormKernelHalfInv(i,j) = NormKernelHalfInv(i,j) + u_norm(i,collidx(k)) * sqrt_norm(k)**(-1) * u_norminv(collidx(k),j) 
         End Do
      End Do
   End Do
   tempmat(:,:) = Matmul(NormKernelHalf(:,:), NormKernelHalf(:,:)) 

   ! check
   Print *, "squareroot of norm kernel check"
   Print *, "max diff N - N^{1/2} N^{1/2}: ", maxval(Abs( NormKernel(:,:) - tempmat(:,:)))
   Print *, "squareroot of norm kernel check completed"

   !----------------------------------------------------------------------------------------------------
   ! Allocation involving colldim
   !----------------------------------------------------------------------------------------------------      
   Allocate( H_coll(colldim,colldim), g_coll(colldim,colldim), RBMenergy(colldim), SortedOrder_Hcoll(colldim), g_collinv(colldim,colldim),&
   & RBMstrength(colldim), ugn(nrbm,colldim), ugn2(colldim,nrbm), RBMstrength1(colldim), RBMstrength2(colldim), &
   & X_QRPA_RBM(nuv,colldim), Y_QRPA_RBM(nuv,colldim),RBMstrengthfromXY(colldim,nop))
   H_coll = 0.0d0; g_coll = 0.0d0; RBMenergy = 0.0d0; SortedOrder_Hcoll = 0; g_collinv = 0.0d0
   RBMstrength = 0.0d0; ugn = 0.0d0; ugn2 = 0.0d0; RBMstrength1 = 0.0d0; RBMstrength2 = 0.0d0
   X_QRPA_RBM = 0.0_pr; Y_QRPA_RBM = 0.0_pr; RBMstrengthfromXY = 0.0_pr

   !----------------------------------------------------------------------------------------------------
   ! Hcoll
   !----------------------------------------------------------------------------------------------------      
   H_coll = 0.0_pr
   Do j = 1, colldim
      Do i = 1, colldim
         Do l = 1, nrbm
            Do k = 1, nrbm
               H_coll(i,j) = H_coll(i,j) + u_norminv(collidx(i),k) * HamiltonianKernel(k,l) * (u_norm(l,collidx(j))) / (sqrt_norm(i)*sqrt_norm(j))
            End Do
         End Do
      End Do
   End Do

   !----------------------------------------------------------------------------------------------------
   ! Collective Hamiltonian diagonalization
   !----------------------------------------------------------------------------------------------------     
   Print *, "**** Diagonalizing collective Hamiltonian ****" 
   Call DiagGenComplexMat(colldim, H_coll, g_coll, RBMenergy, g_collinv, ierr)
   If (ierr .ne. 0) Then
      Print *,   "Hcoll diagonalization ERROR: INFO = ", ierr
      Write(0,*) "Hcoll diagonalization ERROR: INFO = ", ierr
      Stop "Hcoll diagonalizatioin failed"
   End If
   ! ----------------------------------------------------------
   SMALLESTFIRST=.true.
   If(Allocated(realtemp)) Deallocate(realtemp)
   Allocate(realtemp(colldim))
   realtemp(1:colldim) = Dble(RBMenergy(1:colldim))
   Call SORT2(colldim, realtemp, SortedOrder_Hcoll, SMALLESTFIRST)

   Print *, "QRPA energies"
   Print *, " i    Re E       Im E"
   Do k = 1, colldim
      Write(*,'(I5,2ES20.10)') k, RBMenergy(SortedOrder_Hcoll(k))
   End Do
   
   If(Allocated(tempmat)) Deallocate(tempmat)
   Allocate(tempmat(colldim,colldim))
   tempmat = 0.0d0
   Do i = 1, colldim
      Do j = 1, colldim
         Do k = 1, colldim
            tempmat(i,j) = tempmat(i,j) + g_coll(i,SortedOrder_Hcoll(k)) * RBMenergy(SortedOrder_Hcoll(k)) * g_collinv(SortedOrder_Hcoll(k),j)
         End Do
      End Do
   End Do
   Print *, "check Hcoll diagonalization"
   Print *, "max diff Hcoll - g Omega g^{-1}: ", maxval( abs( tempmat(1:colldim,1:colldim) - H_coll(1:colldim,1:colldim)))


   !----------------------------------------------------------------------------------------------------
   ! Caclulate strength
   !----------------------------------------------------------------------------------------------------      
   ugn = 0.0d0; ugn2 = 0.0d0
   Do i = 1, nrbm
      Do l = 1, colldim
         Do k = 1, colldim
            ugn(i,l)  = ugn(i,l)  + u_norm(i,collidx(k)) * g_coll(k,l) / sqrt_norm(k)       !  \tilde{u}
            ugn2(l,i) = ugn2(l,i) + g_collinv(l,k) * u_norminv(collidx(k),i) / sqrt_norm(k) !  \tilde{\tilde{u}}
         End Do
      End Do
   End Do

   RBMstrength1 = 0.0d0; 
   RBMstrength1(1:colldim) = Matmul( SFtrain(1:nrbm),     ugn(1:nrbm, 1:colldim))  ! S \tilde{u}

   RBMstrength2 = 0.0d0
   If(VARIATION==1) Then
      RBMstrength2(1:colldim) = Matmul( ugn2(1:colldim,1:nrbm), Conjg(SFtrain(1:nrbm))) ! \tilde{\tilde{u}} S*
   Else If(VARIATION == 2) Then
      RBMstrength2(1:colldim) = Matmul( ugn2(1:colldim,1:nrbm), TFtrain(1:nrbm))
   End If

   ! QRPA strength  |<\tilde{i}|F|0>|^2
   Do i = 1, colldim
      RBMstrength(i) = RBMstrength1(i)*RBMstrength2(i)
   End Do   

   !----------------------------------------------------------------------------------------------------
   ! Calculate QRPA XY amplitudes
   !----------------------------------------------------------------------------------------------------      

   X_QRPA_RBM(1:nuv,1:colldim) = Matmul( Xtrain(1:nuv,1:nrbm), ugn(1:nrbm,1:colldim)) !+ Matmul(Conjg(YN_train_c(1:nuv,1:nrbm/2)), ugn(nrbm/2+1:nrbm,i))
   Y_QRPA_RBM(1:nuv,1:colldim) = Matmul( Ytrain(1:nuv,1:nrbm), ugn(1:nrbm,1:colldim)) !+ Matmul(Conjg(XN_train_c(1:nuv,1:nrbm/2)), ugn(nrbm/2+1:nrbm,i))

   ! check normalization
   
   Do i = 1, colldim
      Norm = Dble(Dot_Product(X_QRPA_RBM(:,i), X_QRPA_RBM(:,i)) - Dot_Product(Y_QRPA_RBM(:,i), Y_QRPA_RBM(:,i)))   
!      Print *, "QRPA emulator vector normalization i = ", i, Norm

      If(Norm>0.0d0) Then
         X_QRPA_RBM(:,i) = X_QRPA_RBM(:,i) / Sqrt(Norm)
         Y_QRPA_RBM(:,i) = Y_QRPA_RBM(:,i) / Sqrt(Norm)
      Else
         X_QRPA_RBM(:,i) = X_QRPA_RBM(:,i) / Sqrt(-Norm) * iunit
         Y_QRPA_RBM(:,i) = Y_QRPA_RBM(:,i) / Sqrt(-Norm) * iunit
      End If

      Norm = Dble(Dot_Product(X_QRPA_RBM(:,i), X_QRPA_RBM(:,i)) - Dot_Product(Y_QRPA_RBM(:,i), Y_QRPA_RBM(:,i)))
!      Print *, "QRPA emulator vector normalized: i = ", i, Norm

   End Do

   Do j = 1, nop
      Do i = 1, colldim
         ! |<i|Fj|0>|^2
         RBMstrengthfromXY(i,j) = Abs( Dot_Product( X_QRPA_RBM(:,i), F20(:,j)) + Dot_Product( Y_QRPA_RBM(:,i), F02(:,j)))**2
      End Do
   End Do

   Open(unit=outp, file=rbm_outputfile, action='write')
   Write(outp,*) "Results of the RBM"
   Write(outp,*) "1idx&
               &           2ReEnergy           3ImEnergy        4ReStr(SuuS)&
               &        5ImStr(SuuS)     6StrfromXY(op1)     7StrfromXY(op2) ..."
   Do i = 1, colldim
      Write(outp,'(I5,50ES20.10)') i,  RBMenergy(SortedOrder_Hcoll(i)), RBMstrength(SortedOrder_Hcoll(i)), &
      & (RBMstrengthfromXY(SortedOrder_Hcoll(i),j), j = 1, nop) 
   End Do
   Close(outp)
   Print *, "**** RBM results wrote on a file ****"

   If(STRENGTH_EMULATOR_RUN) Then
      Print *, "**** Strength Emulator Run Start ****"
      ! emulator run
      qrpa_omega = qrpa_omega_emulatorrun_start
      sumrule_emulator = 0.0d0
      Open(emu,file=emulator_outputfile,action='write')

      Write(emu,*) '#1qrpa_omega(Re) 2qrpa_omega(Im) 3SFomega(Re) 4SFomega(Im) 5-13Sumrule'
      Do l = 1, nmax_emulatorrun

         Call FAMEmulator(colldim, qrpa_omega, RBMenergy, RBMstrength, SFomega)
         If( Abs(Dble(qrpa_omega))>1.0d-15 ) Then
            Do j = -4, 4
               sumrule_emulator(j) = sumrule_emulator(j) - Dble(qrpa_omega)**j * Aimag(SFomega)/pi * Dble(qrpa_omega_emulatorrun_step)
            End Do
         End If
         Write(emu,'(50ES20.10)') qrpa_omega, SFomega, (sumrule_emulator(j), j = -4, 4)

         qrpa_omega = qrpa_omega + qrpa_omega_emulatorrun_step

      End Do
      Close(emu)

      Print *, "Sum rules in the energy region from ", Dble(qrpa_omega_emulatorrun_start), " to ", Dble(qrpa_omega-qrpa_omega_emulatorrun_step)
      Do k = -4, 4
         Write(*,'(A,I2,A,ES20.10)') "m_", k, " = ", sumrule_emulator(k)
      End Do


   End If

End Program RBM_FAM
