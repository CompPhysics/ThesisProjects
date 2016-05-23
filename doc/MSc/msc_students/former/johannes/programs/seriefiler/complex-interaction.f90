!             Program block complex-interaction.f90
!
!             Author:   Gaute Hagen and Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, Universities of Bergen and  Oslo, Norway
!             E-MAIL:   gaute.hagen@ift.uib.no, morten.hjorth-jensen@fys.uio.no
!             LANGUAGE: F90/F95  
!             LAST UPGRADE : February 2005
!
!             This file sets up a complex two-body interaction in the CoM 
!             system in momentum space using complex scaling.
!             It allows for the use of the free interaction and a similarity
!             transformed interaction as well.
!
!       setup the  free interaction 
!

SUBROUTINE g_channel(i)
  USE relcm_gmatrix
  USE partial_waves
  USE constants
  USE wave_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i  ! loop variable over NN channels
  INTEGER :: ncoup, pq_confs, p_confs, mspace
  COMPLEX(DPC), ALLOCATABLE :: vkk(:,:), heff(:,:)

  
  ncoup=1
  IF ( orb_lrel_max(i) == jang_rel(i) ) ncoup = 1
  IF ( orb_lrel_max(i) /= jang_rel(i) ) ncoup = 2
  ALLOCATE(vkk (ncoup*n_k,ncoup*n_k))
  !     setup the NN interaction
  CALL potential_interface(ncoup,vkk,i)
  
  
  !     construct free V-matrix in k-space
  WRITE(6,'(12H Channel Nr:,I3,7H l_min:,I3,7H l_max:,I3,3H J:,I3,3H S:,I3,4H Tz:,I3)') &
     i, orb_lrel_min(i), orb_lrel_max(i), jang_rel(i), spin_rel(i), iso(i)
  SELECT CASE(type_of_interaction)
  CASE('free')
     vfree(1:ncoup*n_k,1:ncoup*n_k,i) = vkk(1:ncoup*n_k,1:ncoup*n_k)
  CASE('renorm')
     ALLOCATE(heff (ncoup*n_k1,ncoup*n_k1))
     CALL g_mtx(ncoup,vkk,heff,i)
     vfree(1:ncoup*n_k1,1:ncoup*n_k1,i) = heff(1:ncoup*n_k1,1:ncoup*n_k1)
     DEALLOCATE(heff)
  END SELECT
  DEALLOCATE(vkk)

END SUBROUTINE g_channel

!
!       Compute the T(G)-mtx for the actual channel 
!
SUBROUTINE g_mtx(ncoup,vzz,heff,ichan)
  USE wave_functions
  USE constants
  USE partial_waves
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ncoup, ichan
  COMPLEX(DPC), DIMENSION(ncoup*n_k,ncoup*n_k), INTENT(IN) :: vzz
  COMPLEX(DPC), DIMENSION(ncoup*n_k1,ncoup*n_k1), INTENT(INOUT) :: heff

  INTEGER :: i, j, ntot, i1, np, nq, j1,j2, tisoz, no_zero 
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:) :: ham, cvec, temp
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:) ::  ceig
  COMPLEX(DPC) :: sum, input_energy


  ! dimension of vectors and matrices
  ntot = ncoup*n_k
  ALLOCATE( ham(ntot, ntot), cvec(ntot,ntot), ceig(ntot), temp(ntot,ntot)) 
  ! setup hamiltonian to be diagonalized
  ham = DCMPLX(0.0_DP,0.0_DP) 
  cvec = dcmplx(0.0_dp, 0.0_dp)
  ceig = dcmplx(0.0_dp, 0.0_dp)
  

  CALL complex_hamiltonian(ham,vzz, ntot, ichan)
  temp = ham
  no_zero = 0 
  do i = 1, ntot
     do j = 1, ntot
        if ( abs(ham(i,j)) == 0 ) no_zero = no_zero + 1
     end do
  end do
  
  if ( no_zero > ntot*ntot/2 ) then
     call lapack_diag(temp, cvec, ceig, ntot )
  else
     CALL diag_exact( temp, cvec, ceig, ntot )
  end if
  
  ! construct renormalized nn-interaction via the Lee-Suzuki sim. transform
  !
  ! size of model space     : np = n_k1*ncoup
  ! size of complement space: nq = n_k2*ncoup
  np = n_k1*ncoup; nq = n_k2*ncoup
  CALL effective_int( np, nq, ntot, ncoup, cvec, ceig, heff )   
  !
  ! subtract diagonal elements to obtain Veff, 
  ! and compare exact with effective interactions: 
  !
  tisoz = iso(ichan) 
  DO i = 1, np
     i1 = i
     IF ( i > n_k1 ) i1 = i-n_k1
     heff(i,i) = heff(i,i) - ( krel(i1) * krel(i1) )/p_mass(tisoz)  
     DO j = 1, np
        j1 = j
        IF ( j > n_k1 ) j1 = j-n_k1
        heff(i,j) = heff(i,j)/SQRT( wkrel(i1)*wkrel(j1) )/krel(i1)/krel(j1)
     ENDDO
  ENDDO

!  do i = 1, np 
!     write(6,*) abs(krel(i)), dble(heff(i,i)), aimag(heff(i,i)) 
!  end do

  DEALLOCATE(ham, temp, cvec );   DEALLOCATE(ceig)

END SUBROUTINE g_mtx

!
! complex scaled hamiltonian in momentum representation 
!
SUBROUTINE complex_hamiltonian(h,vzz,ntot,ichan)
  USE wave_functions
  USE constants
  USE partial_waves
  IMPLICIT NONE
  REAL(DP) :: delta
  INTEGER, INTENT(IN) :: ntot, ichan
  COMPLEX(DPC), DIMENSION(ntot,ntot), INTENT(INOUT) :: h
  COMPLEX(DPC), DIMENSION(ntot,ntot), INTENT(IN) :: vzz
  INTEGER :: i, j, i1, i2, tisoz
  COMPLEX(DPC) :: wi

  tisoz = iso(ichan) 
  h = 0.
  DO i = 1, ntot 
     i1 = i
     IF ( i > n_k ) i1 = i-n_k
     h(i,i) = ( krel(i1) * krel(i1) )/p_mass(tisoz)  + & 
          (krel(i1)) * (krel(i1)) * wkrel(i1) * Vzz(i,i)
     DO j = 1, ntot 
        i2 = j
        IF ( j > n_k ) i2 = j-n_k
        IF (i /= j ) THEN
           h(i,j) = SQRT( wkrel(i2) * wkrel(i1) ) *  krel(i1) * krel(i2) * Vzz(i,j) 
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE complex_hamiltonian


! setting up potential vkk for complex arguments z,z', used in   
! diagonalization of hamiltonian 
!
!          This routine contains old-fashioned common blocks
!          acting as interface between the effective interaction program
!          and the potential routines of the Bonn type, Nijmegen or 
!          Argonne groups. 
!          vkk is in units of MeV^-2,  
!          rgkk are mesh points in rel coordinates, 
!          units of fm^-1
!          in partial wave basis, v(6) 
!          means
!          v(1)= singlet uncoupled
!          V(2)= triplet uncoupled
!          V(3)= triplet coupled, < J+1 | V | J+1>
!          V(4)= triplet coupled, < J-1 | V | J-1>
!          V(5)= triplet coupled, < J+1 | V | J-1>
!          V(6)= triplet coupled, < J-1 | V | J+1>
!
SUBROUTINE potential_interface(ncoup,vzz,ichan)
  USE wave_functions
  USE constants
  USE partial_waves
  IMPLICIT NONE
  INTEGER :: i,j, ncoup, k, n, jxx, ichan, inn, isospin_tz, spin, ix, iy
  !COMPLEX(DPC) :: v, ymev, xmev
  REAL(DP) :: v, ymev, xmev
  CHARACTER (LEN=4) :: label
  COMMON /cnn/ inn
  COMMON /cpot/ v(6),xmev,ymev
  COMMON /cstate/ jxx, heform, sing, trip, coup, endep, label
  LOGICAL :: sing, trip, coup, heform, endep
  COMPLEX(DPC), DIMENSION(ncoup*n_k,ncoup*n_k), INTENT(INOUT) :: vzz

  heform=.FALSE.
  jxx=jang_rel(ichan)
  sing=spin_rel(ichan) == 0
  trip=spin_rel(ichan) == 1
  coup=(ncoup == 2) 
  isospin_tz=iso(ichan)
  spin = spin_rel(ichan)

  SELECT CASE ( isospin_tz)
  CASE (-1)
     inn = 1    !  pp case = 1
  CASE (0) 
     inn = 2    !  pn case = 2  if all inn=2, no CSB or ISB
  CASE ( 1)
     inn = 3    !  nn case = 3  if pp == nn, only ISB
  END SELECT
  DO i=1,n_k
     xmev = dble( krel(i) )
     ix = i
     IF (ncoup == 2)   k=i+n_k
     DO j=1,n_k
        ymev= dble( krel(j) )
        iy = j
        IF (ncoup == 2) n=j+n_k
        
        SELECT CASE (type_of_pot)
!        CASE('CD-bonn')
!           CALL cdbonn
!        CASE('Idaho-A')
!           CALL idaho
!        CASE('Idaho-B')
!           CALL idaho
!        CASE ( 'reid93')
!           CALL reid93
!        CASE ( 'argonnev18')
!           CALL argonp
!        CASE ( 'nij-I-II')
!           CALL nijm2
        CASE ( 'n3lo')
 !          write(6,*) 'n3lo'
          CALL n3lo
        END SELECT
      
        !  SELECT CASE (type_of_pot)
        !  CASE('cmplx_cdonn')
        ! CALL cmplx_cdbonn
        !  END SELECT
        IF (sing ) THEN
           vzz(j,i)= dcmplx( v(1) )
        ELSEIF ((trip).AND.(.NOT.coup )  ) THEN 
           vzz(j,i) = dcmplx( v(2) )
        ELSEIF (coup ) THEN
           vzz(j,i)= dcmplx( v(4) )
           vzz(j,k)= dcmplx( v(5) )
           vzz(n,i)= dcmplx( v(6) )
           vzz(n,k)= dcmplx( v(3) )
        ENDIF
     ENDDO
  ENDDO

END  SUBROUTINE potential_interface

!
! calculate effective interaction for gamow shell model calcs
!
subroutine effective_int( np, nq, ntot, ncoup, cvec, ceig, heff )   
  use constants
  USE wave_functions

  implicit none
  complex*16, dimension(ntot,ntot), intent(in) :: cvec
  complex*16, dimension(ntot,ntot) :: cvec_temp
  complex*16, dimension(ntot), intent(in) :: ceig
  integer, intent(in) ::  np, nq, ntot, ncoup
  integer :: npp, k1, i1,i,k2,nqq,j
  complex*16 :: cvec_pp(np,np), cvec_qp(nq,np), cvec_model(nq,np) 
  complex*16 :: ceig_p(np), ceig_model(np)
  double precision,  dimension(ntot) ::  temp
  integer ::  model(np), orb, orbit_no(ntot)
  double precision :: cvec_max, overlap(ntot)
  complex*16 :: e_a, e_b
  complex*16, dimension(np,np), intent(inout) :: heff
  real*8 :: a1,a2,b1,b2

  ! calculate P-space overlap of all eigenvectors 
  ! loop over all eigen vectors
  do orb = 1, ntot

     overlap(orb) = 0.
     do i = 1, n_k1
        if ( ncoup == 1 )then
           overlap(orb) = overlap(orb) + abs( cvec(i,orb) )**2
        elseif( ncoup == 2 ) then
           overlap(orb) = overlap(orb) + abs( cvec(i,orb) )**2 + abs( cvec(n_k + i,orb) )**2
        end if
     end do
     orbit_no(orb) = orb
  end do

  ! sort all overlaps and corresponding orbit numbers 
  call eigenvalue_sort(overlap,orbit_no, ntot)


  cvec_pp = 0.; cvec_qp = 0.
  ! Set up eigenvalues and eigenvectors for model space and excluded space
  if ( ncoup == 1 ) then
     DO i=1, np
        ! loop over all model space coefficients of exact eigenvectors |k>
        ! 
        k1 = orbit_no(i)

        DO j = 1, n_k

           IF ( j <= n_k1 ) THEN
              cvec_pp(j,i) = cvec(j,k1)
              ceig_p(i) = ceig(k1)
           ELSEIF ( j > n_k1 ) THEN
              cvec_qp(j-n_k1,i) = cvec(j,k1)
           ENDIF
        ENDDO
     ENDDO
  end if

  if ( ncoup == 2 ) then
     DO i=1, np
        ! loop over all model space coefficients of exact eigenvectors |k>
        ! 
        k1 = orbit_no(i)
        ceig_p(i) = ceig(k1)
        DO j = 1, n_k1
           cvec_pp(j,i) = cvec(j,k1)
           cvec_pp(j+n_k1,i) = cvec(j+n_k,k1)
        ENDDO
     ENDDO

     DO i=1, np
        k1 = orbit_no(i)
        do j = 1, n_k2
           cvec_qp(j,i) = cvec(n_k1+j,k1)
           cvec_qp(j+n_k2,i) = cvec(n_k+n_k1+j,k1)
        end do
     end DO
  end if


  heff = 0.


  call lee_suzuki2( cvec_pp, cvec_qp, ceig_p, np, nq, heff )
  !  ! subtract diagonal terms


end subroutine effective_int

!
! Lee Suzuki similarity transformation 
!
subroutine lee_suzuki( cvec_pp, cvec_qp, ceig, np, nq, heff )
  !use gms_constants
  !use ang_mom_functions
  !use gms_arrays
  !use shell_model_ham

  implicit none
  integer, intent(in) :: np, nq
  complex*16, dimension(np,np),intent(in) :: cvec_pp
  complex*16, dimension(nq,np),intent(in)  :: cvec_qp
  complex*16, dimension(np), intent(in) :: ceig
  complex*16, dimension(np,np) :: temp, omega2, sim1, sim2, omega_2inv, u, u_inv, & 
       cvec_pp_inv
  complex*16, dimension(nq,np) :: omega
  complex*16, dimension(np) :: ceig_p, omega2_eig, eigval2
  complex*16, dimension(np,np) :: eigen_vec, vl, heff_rhs
  complex*16, dimension(np,np), intent(out) ::  heff
  double precision, dimension(2*np) :: rwork
  complex*16, dimension(10000) :: work1
  complex*16 :: d, sum1, temp1(np,np), temp2(np,nq), norm, determinant
  integer :: i_p,j_p, j_q, ii, jj, k, i_q , a_p, a_pp 
  INTEGER :: i, lda, ldb, ldvl, LDVR, info, lwork, ILO , IHI
  CHARACTER*1 :: jobvl, JOBVR, BALANC, SENSE
  DOUBLE PRECISION, DIMENSION(np) :: SCALE, RCONDE, RCONDV
  DOUBLE PRECISION :: ABNRM
  integer :: ipiv(np) ,j
  real*8 :: a1, a2, b1, b2

  BALANC = 'N'
  JOBVL = 'N'
  JOBVR = 'V'
  SENSE = 'N'
  LDA = np
  LDVL = 1
  LDVR = np
  LWORK = 10000

  open(13,file='specter.dat')

  eigen_vec = 0. 

  cvec_pp_inv = cvec_pp
  call cmplxmatinv(cvec_pp_inv, np, d)
  !write(6,*) cvec_pp_inv

  omega = matmul( cvec_qp, cvec_pp_inv )


  ! non-hermitian effective interaction
  ! setup 2-p effective interaction in P-space

  !write(6,*) 'step1'
  heff = 0.
  do i_p = 1, np
     do j_p = 1, np 

        heff(i_p,j_p) = sum( cvec_pp(i_p,:)*ceig(:)*cvec_pp(j_p,:) ) 

        sum1 = 0.
        do k = 1, np
           do i_q = 1, nq
              sum1 = sum1 + cvec_pp(i_p,k) * ceig(k) * cvec_qp(i_q,k) * omega(i_q,j_p) 

           end do
        end do
        heff(i_p,j_p) = heff(i_p,j_p) + sum1

        !write(6,*) heff(i_p,j_p)
     end do
  end do
  !write(6,*) 'step2'


  !write(6,*) cvec_pp_inv
  !                                      
  ! *****  organizing matrix (P(1 + omega * omega)P)
  !
  temp = 0.; omega2 = 0.
  omega2 = matmul(transpose(cvec_pp_inv),cvec_pp_inv)
  temp = omega2 


  !
  ! calculate sqrt and inverse sqrt of matrix (P(1 + omega * omega)
  !call sqrtmat_newton( omega2, U, U_inv ,np ) 
  call sqrtmat_db( omega2, U, U_inv ,np ) 
  do i = 1, np
     do j = 1, np
        !if ( abs(temp1(i,j)- u(i,j)) < 1.e-6 ) cycle 
        ! write(6,*) u(i,j) !, temp1(i,j) 
     end do
  end do
  heff_rhs =  matmul( U, matmul( heff,u_inv ) ) 

  ! check if heff is symmetrized:
  heff = 0.
  heff = heff_rhs
  do i_p = 1, np
     do j_p = 1, np
        ! make heff manifestly symmetric

        if ( i_p /= j_p ) then
           if ( abs(heff(i_p,j_p)- heff(j_p,i_p)) < 1.e-6 ) cycle 
           write(13,*) 'sym test', heff(i_p,j_p), heff(j_p,i_p)
        end if
     end do
  end do

  ! diagonalize 2p-effective shell model hamiltonian
  !CALL ZGEEV( JOBVL, JOBVR, np, heff_rhs, LDA, ceig_p, VL, LDVL, eigen_vec, LDVR, &
  !     WORK1, LWORK, RWORK, INFO )
  !write(6,*) 'step3'
  call diag_hpp(heff_rhs, eigen_vec, ceig_p, np, 0 )
  !write(6,*) 'step4'
  call eigenvalue_sort_cmplx(ceig_p, np)
  call eigenvalue_sort_cmplx(ceig, np)
  ! compare spectrum from exact and P-space diagonalization
  write(13,*)
  write(13,*) 'Compare model space two-body spectrum with exact spectrum:' 
  do i_p = 1, np
     a1 = real( ceig_p(i_p))
     a2 = aimag( ceig_p(i_p))
     b1 = real( ceig(i_p) )
     b2 = aimag( ceig(i_p) ) 
    
     write(6,*) a1, b1
  end do



end subroutine lee_suzuki



subroutine diag_hpp(h_pp, cvec, ceig, n, option )

  implicit none
  integer, intent(in) :: n
  complex*16, dimension(n,n), intent(in) :: h_pp
  integer :: i,j, i1, kvec, lwork, option, il ,iu, info
  real*8 :: Avec(n*(n+1))
  real*8,allocatable,dimension(:) :: work 
  real*8 :: thresh
  complex*16 :: cvl, cvu
  logical :: flag(n)  
  complex*16, dimension(n), intent(out) :: ceig
  complex*16, dimension(n,n), intent(out) :: cvec
  real*8 vl,vu

  i1 = 0
  do i =  1, n
     do j =  1, i
        i1 = i1 + 1

        avec(i1) = dble(h_pp(j,i))
        avec(i1+n*(n+1)/2) = aimag(h_pp(j,i))

     end do
  end do

  lwork =  300*n
  thresh = 30.0
  kvec = n

  !option = 4
  info = 1

  allocate( work(lwork) )

  !write(6,*) 'begynner diag' 
  call cs(n,avec,ceig,kvec,cvec,lwork,work,thresh, & 
       option,il,iu,cvl,cvu,flag,info)


  deallocate( work )
  !write(12,*) 'three-body resonant state energy:'
  !write(6,*) 'ferdig diag' 

end subroutine diag_hpp


!
! eigenvalue sort
! sort cmplx vector real(a(1))<real(a(2)) < ... < real(a(n))
!
SUBROUTINE eigenvalue_sort(A,B, n)
  IMPLICIT NONE

  INTEGER :: i, j, n
  double precision, DIMENSION(n), INTENT(INOUT) :: A
  integer, dimension(n), intent(inout) :: B
  double precision :: temp, temp1
  double precision, DIMENSION(n) :: temp2
  integer :: orb

  DO i = 1, n
     DO j = 1, n
        IF ( abs( A(i) )  > abs( A(j) ) ) THEN

           temp = A(i)
           A(i) = A(j) 
           A(j) = temp

           orb = B(i) 
           B(i) = B(j)
           B(j) = orb

        END IF
     END DO
  END DO

END SUBROUTINE eigenvalue_sort

!
! Lee Suzuki similarity transformation 
!

SUBROUTINE lee_suzuki2( cvec_pp, cvec_qp, ceig, np, nq, heff )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: np, nq
  COMPLEX*16, DIMENSION(np,np), INTENT(IN) :: cvec_pp
  COMPLEX*16, DIMENSION(nq,np), INTENT(IN)  :: cvec_qp
  COMPLEX*16, DIMENSION(np), INTENT(IN) :: ceig
  COMPLEX*16, DIMENSION(np,np) :: temp, omega2, sim1, sim2, omega_2inv, u, u_inv, & 
       cvec_pp_inv
  COMPLEX*16, DIMENSION(nq,np) :: omega
  COMPLEX*16, DIMENSION(np) :: ceig_p, omega2_eig, eigval2
  COMPLEX*16, DIMENSION(np,np) :: eigen_vec, vl, heff_rhs
  COMPLEX*16, DIMENSION(np,np), INTENT(INOUT) ::  heff
  DOUBLE PRECISION, DIMENSION(2*np) :: rwork
  COMPLEX*16, DIMENSION(10000) :: work1
  COMPLEX*16 :: d, sum1, temp1(np,np), temp2(np,nq), norm, determinant
  INTEGER :: i_p,j_p, j_q, ii, jj, k, i_q , a_p, a_pp 
  INTEGER :: i, lda, ldb, ldvl, ldvr, info, lwork, ilo , ihi
  CHARACTER*1 :: jobvl, jobvr, balanc, sense
  DOUBLE PRECISION, DIMENSION(np) :: scale, rconde, rcondv
  DOUBLE PRECISION :: abnrm
  INTEGER :: ipiv(np) ,j
  DOUBLE PRECISION :: a1, a2, b1, b2

  balanc = 'n';  jobvl = 'n' ;  jobvr = 'v';  sense = 'n';  lda = np
  ldvl = 1;  ldvr = np;  lwork = 10000
  eigen_vec = 0. 
  temp1 = TRANSPOSE(cvec_pp) 
  CALL zgeev( jobvl, jobvr, np, temp1, lda, omega2_eig, vl, ldvl, eigen_vec, ldvr, &
       work1, lwork, rwork, info )
  determinant = PRODUCT(omega2_eig(:))
  !  write(6,*) 'check determinant', determinant
  ! the P->Q space transformation matrix, omega 
  cvec_pp_inv = cvec_pp
  call cmplxmatinv(cvec_pp_inv, np, d)
  DO i_p = 1, np
     DO i_q = 1, nq
        omega(i_q,i_p) = SUM( cvec_pp_inv(:,i_p)*cvec_qp(i_q,:) )
     ENDDO
  ENDDO
  ! non-hermitian effective interaction
  ! setup 2-p effective interaction in P-space
  heff = 0.
  DO i_p = 1, np
     DO j_p = 1, np 
        heff(i_p,j_p) = SUM( cvec_pp(i_p,:)*ceig(:)*cvec_pp(j_p,:) ) 
        sum1 = 0.
        DO k = 1, np
           DO i_q = 1, nq
              sum1 = sum1 + cvec_pp(i_p,k) * ceig(k) * cvec_qp(i_q,k) * omega(i_q,j_p) 
           ENDDO
        ENDDO
        heff(i_p,j_p) = heff(i_p,j_p) + sum1
     ENDDO
  ENDDO
  ! organizing the matrix (P(1 + omega * omega)P)
  omega2 = MATMUL(TRANSPOSE(cvec_pp_inv),cvec_pp_inv)
  ! calculate sqrt and inverse sqrt of matrix (P(1 + omega * omega)P)
  CALL sqrtmat( omega2, U, U_inv ,np ) 
  heff_rhs =  matmul( U, matmul( heff,u_inv ) ) 
  ! check if heff is symmetrized:
  heff = 0.
  heff = heff_rhs
  DO i_p = 1, np
     DO j_p = 1, np
        ! make heff manifestly symmetric
        IF ( i_p /= j_p ) THEN
           IF ( ABS(heff(i_p,j_p)- heff(j_p,i_p)) < 1.E-6 ) CYCLE 
           !           WRITE(6,*) 'sym test', heff(i_p,j_p), heff(j_p,i_p)
        ENDIF
     ENDDO
  ENDDO
  ! diagonalize 2p-effective shell model hamiltonian
  CALL zgeev( jobvl, jobvr, np, heff_rhs, lda, ceig_p, vl, ldvl, eigen_vec, ldvr, &
       work1, lwork, rwork, info )
  call eigenvalue_sort_cmplx(ceig_p, np)
  call eigenvalue_sort_cmplx(ceig, np)

  !   compare spectrum from exact and P-space diagonalization
  !  WRITE(6,*) 'Compare model space two-body spectrum with exact spectrum:' 
  DO i_p = 1, np
     a1 = REAL( ceig_p(i_p))
     a2 = AIMAG( ceig_p(i_p))
     b1 = REAL( ceig(i_p) )
     b2 = AIMAG( ceig(i_p) ) 
     if ( abs( ceig_p(i_p) - ceig(i_p) ) > 0.001 ) & 
          write(6,*) 'eigenvalues wrong', A1, a2, B1, b2
  ENDDO

END SUBROUTINE lee_suzuki2



! 
! routine for calculating the principal square root ( and inverse ) 
! of a general matrix AA, i.e. solving the matrix equation: AA - A = 0. 
! The routine is built on the property : 
!
!            | 0  AA |      | 0    A | 
! sign(B) = (|       |) =  (|        |)
!            | 1   0 |      | A^-1 0 |
! 
! the sign of the matrix AA is calculated using Newtons iteration method,
! see Higham et. al. ref. Numerical Algorithms 15 (1997) 227-242 
! 

SUBROUTINE sqrtmat( aa, a, a_inv ,n ) 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: aa
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: a, a_inv
  COMPLEX*16, DIMENSION(2*n,2*n) :: x, x0, x1, x2
  COMPLEX*16, DIMENSION(n,n) :: i_mat, temp
  INTEGER :: i,j,k
  COMPLEX*16 :: d

  ! setup identity matrix
  i_mat = 0.
  DO i = 1, n
     i_mat(i,i) = 1.d0
  ENDDO

  !setup starting matrix x0 
  x0 = 0.
  DO i = n+1, 2*n
     DO j = 1, n
        x0(j,i) = aa(j,i-n)
     ENDDO
  ENDDO

  DO i = 1, n
     DO j = n+1, 2*n
        x0(j,i) = i_mat(j-n,i)
     ENDDO
  ENDDO
  k = 0
  temp = 0.
  x1 = 0.; x2 = 0.; x = 0.
  DO WHILE( MAXVAL(ABS(temp-aa)) > 1.E-14 .AND.  k < 1000 )
     x1 = x0 
     x2 = x0 
     CALL cmplxmatinv(x2,2*n,d)
     x = 0.5d0 * ( x1 + x2 ) 
     x0 = x 
     k = k + 1
     DO i = 1, n
        DO j =1, n
           a(i,j) = x(i,j+n)
           a_inv(i,j) = x(i+n,j) 
        ENDDO
     ENDDO
     temp = MATMUL( a,a )
  ENDDO
  !WRITE(6,*) 'number of iterations in sqrt_mat:', k

END SUBROUTINE sqrtmat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!fast diagonalizing of complex symmetric matrices

subroutine diag_exact( h, cvec, ceig, n )
  use constants
  implicit none
  integer, intent(in) :: n
  complex(dpc), dimension(n,n), intent(in) :: h
  integer :: k, np
  integer :: p, i,j, i1, kvec, lwork, option, il ,iu, info
  real(dp) :: A(n*(n+1)), work(300*n)
  real(dp) :: thresh
  complex(dpc) :: cvl, cvu
  logical :: flag(n)  
  complex(dpc) :: ceig(n), sum1
  complex(dpc) :: cvec(n,n)
  real(dp) vl,vu

  lwork = 300*n
  thresh = 30.0
  kvec = n

  option = 4
  info = 1


  i1 = 0
  do i =  1, n
     do j =  1, i
        i1 = i1 + 1

        a(i1) = dble(h(j,i))
        a(i1+n*(n+1)/2) = aimag(h(j,i))

     end do
  end do


  call cs(n,a,ceig,kvec,cvec,lwork,work,thresh, option,il,iu,cvl,cvu,flag,info)


end subroutine diag_exact


! jobvl, jobvr, np, temp1, lda, omega2_eig, vl, ldvl, eigen_vec, ldvr, &
!       work1, lwork, rwork, info
SUBROUTINE lapack_diag(h, cvec, ceig, n )
  
  use constants
  implicit none
  integer, intent(in) :: n
  complex*16, dimension(n,n), intent(in) :: h
  COMPLEX*16, DIMENSION(n,n), intent(out) :: cvec
  COMPLEX*16, DIMENSION(n,n) ::  vl
  COMPLEX*16, DIMENSION(n), intent(out) :: ceig
  DOUBLE PRECISION, DIMENSION(2*n) :: rwork
  COMPLEX*16, DIMENSION(3168) :: work1
  INTEGER :: lda, ldb, ldvl, ldvr, info, lwork, ilo , ihi
  CHARACTER*1 :: jobvl, jobvr
  DOUBLE PRECISION, DIMENSION(n) :: scale, rconde, rcondv
  complex*16 :: norm
  integer :: i

  jobvl = 'N' ;  jobvr = 'V';  lda = n
  ldvl = 1;  ldvr = n;  lwork = 3168
  ceig = 0.; cvec = 0.
  CALL zgeev( jobvl, jobvr, n, h, lda, ceig, vl, ldvl, cvec, ldvr, &
       work1, lwork, rwork, info )
  ! write(6,*) 'info', info, work1(1)
  ! berggren normalization
  do i = 1, n
     norm = sum( cvec(:,i)*cvec(:,i) )
     cvec(:,i) = cvec(:,i)/sqrt(norm)
  end do
end SUBROUTINE lapack_diag
