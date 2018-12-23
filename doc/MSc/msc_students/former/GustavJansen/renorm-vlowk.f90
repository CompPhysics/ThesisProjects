!             Program block renorm-vlowk.f90
!
!             Authors:  Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
!             E-MAIL:   morten.hjorth-jensen@fys.uio.no
!             LANGUAGE: F90/F95  
!             LAST UPGRADE : March 2007
!
!                        
!
!                    Begin setup of similarity transformation
!
SUBROUTINE  setup_vlowk
  USE constants
  USE relcm_gmatrix
  USE single_particle_orbits
  USE wave_functions
  USE configurations
  USE partial_waves
  IMPLICIT NONE
  INTEGER :: loop, number_orbits, max_conf

  ! Free possibly used memory
  IF(ALLOCATED(ra)) DEALLOCATE(ra)
  IF(ALLOCATED(wra)) DEALLOCATE(wra)
  IF(ALLOCATED(krel)) DEALLOCATE(krel)
  IF(ALLOCATED(wkrel)) DEALLOCATE(wkrel)
  IF(ALLOCATED(rnlr)) DEALLOCATE(rnlr)
  IF(ALLOCATED(coulomb_relcom)) DEALLOCATE(coulomb_relcom)
  IF (ASSOCIATED(rel_conf%nconfs_rel)) DEALLOCATE(rel_conf%nconfs_rel)
  IF (ASSOCIATED(rel_conf%nconfs_relmodel)) &
    DEALLOCATE(rel_conf%nconfs_relmodel)
  IF (ASSOCIATED(rel_conf%rel_ab)) DEALLOCATE(rel_conf%rel_ab)
  IF (ALLOCATED(v_com_rel)) DEALLOCATE(v_com_rel)
  IF (ASSOCIATED(relcm_conf%nconfs_relcm)) &
    DEALLOCATE(relcm_conf%nconfs_relcm)
  IF (ASSOCIATED(relcm_conf%relcm_ab)) &
    DEALLOCATE(relcm_conf%relcm_ab)

  !     reserve space in memory for mesh point and h.o. wave function  arrays
  ALLOCATE ( ra (n_rel), wra (n_rel));   ALLOCATE ( krel (n_rel), wkrel (n_rel))
  ALLOCATE ( rnlr (n_rel, 0:lmax, 0:nmax) )
  ALLOCATE(coulomb_relcom(0:nmax, 0:nmax, 0:lmax, 6))
  !     set up Coulomb interaction in relative coordinates using Brody-Jacob-Moshinsky
  coulomb_relcom = 0.0_dp  
   IF ( coulomb_included =='yes') CALL coulomb_integral
  !     set up quantum numbers nl and NL relative and cm system
  CALL find_spdata_relcm(number_orbits)
  relcm_sp_data%max_value=number_orbits
  CALL allocate_relcm_array(relcm_sp_data, number_orbits) 
  CALL make_configurations_relcm
  !     set up mesh points for relative coordinates
  CALL vlowk_mesh                   
  !     setup h.o. wave functions, only relative coordinates needed
  CALL vlowkho_wfunction                

!
!     Note: first we find all possible configurations in the relative system
!     only based on h.o. quantum number nl for all partial waves
!     and allocate array for the given partial wave
!     all configs are set up in the call to setup_configurations_rel
!     First we search for the number of configurations in each
!     channel
!

  !WRITE(*,*) "Allocating rel_conf%nconfs_rel: ", no_channels

  ALLOCATE ( rel_conf%nconfs_rel(no_channels) )
  ALLOCATE ( rel_conf%nconfs_relmodel(no_channels) )
  CALL number_configurations_rel(rel_conf)
  ! Find maximum number of configurations
  max_conf=0
  DO loop = 1, no_channels
    IF ( max_conf <  rel_conf%nconfs_rel(loop) ) THEN
         max_conf= rel_conf%nconfs_rel(loop)
    ENDIF
  ENDDO

  ! Take care of the number of subchannels. Max 4 - as spesified in
  ! partial_waves
  ! rel_conf are only the possible combinations of n,l for a specific channel.
  ! This is the same for all subchaZZnnels. 
  ALLOCATE( rel_conf%rel_ab(no_channels, max_conf))
  ALLOCATE ( v_com_rel(no_channels,max_sub*max_conf,max_sub*max_conf))
  v_com_rel = 0.
  CALL setup_configurations_rel(rel_conf)

  ! 
  !     Then we find all possible configurations in the relative and cm system
  !     based on h.o. quantum number nlNL for all partial waves
  !     and allocate array for the given partial wave
  !     all configs are set up in the call to setup_configurations_cm
  !     First we search for the number of configurations in each
  !     channel
  !
  ALLOCATE ( relcm_conf%nconfs_relcm(no_channels) )
  CALL number_configurations_relcm(relcm_conf)
  max_conf=0
  DO loop = 1, no_channels
    IF ( max_conf <  relcm_conf%nconfs_relcm(loop) ) THEN
        max_conf= relcm_conf%nconfs_relcm(loop)
    ENDIF
  ENDDO
  ALLOCATE( relcm_conf%relcm_ab(no_channels, max_conf+max_conf))
  CALL setup_configurations_relcm(relcm_conf)

END SUBROUTINE setup_vlowk
!
!     Obtain the bare interaction in a harmonic oscillator 
!     basis plus the kinetic energy and the Coulomb part. It contains
!     also the CoM correction to the relative coordinates. The latter depends
!     on the mass number of the nucleus
! 
SUBROUTINE vlowk_channel(i, v, nv, vint)
  USE wave_functions
  USE relcm_gmatrix
  USE partial_waves
  USE configurations
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i, nv
  INTEGER :: ncoup, n,np, bra, ket, k1, k2, a, c, subab, subcd
  INTEGER :: la, lb, jang, lim1, lim2, lim3, lim4, nsub, nconf, ida, idb, lookup
  REAL(KIND = 8) :: e_coulomb, vsum, mab, mcd
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: heff(:,:), vzz(:,:)
  REAL(KIND = 8), INTENT(IN) :: v(nv,nv)
  REAL(KIND = 8),  INTENT(INOUT):: &
       vint(rel_conf%nconfs_rel(i)*no_sub(i),rel_conf%nconfs_rel(i)*no_sub(i))
  REAL(KIND(1.0D0)), ALLOCATABLE :: rnl_ab(:), rnl_cd(:)


    !WRITE(*,*) "vlowk_channel called"
    !WRITE(*,*) "n_rel, nv, i:", n_rel, nv, i
  ALLOCATE(rnl_ab(n_rel), rnl_cd(n_rel))
  jang=jang_rel(i)
  ncoup=1
  IF ( orb_lrel_max(i) == jang_rel(i) ) ncoup = 1
  IF ( (orb_lrel_max(i) /= jang_rel(i)) .AND. jang /= 0 ) ncoup = 2
  IF ( (jang == 0) .AND. (orb_lrel_max(i) == 1)) ncoup = 2
  nsub = no_sub(i)
  IF (nv /= nsub*n_rel*ncoup) THEN
    WRITE(*,*) "Error in matrix size: ", i, nv, nsub*n_rel*ncoup
    RETURN
  ENDIF

  ALLOCATE(vzz (ncoup*n_rel*nsub,ncoup*n_rel*nsub))
  vzz = v
  WRITE(6,'(12H Channel Nr:,I3,10H strange:,I3,7H l_min:,I3,7H l_max:,I3,3H J:,I3,3H S:,I3,4H Tz:,I3)') &
       i, strange(i),orb_lrel_min(i), orb_lrel_max(i), jang_rel(i), spin_rel(i), iso(i)
  ALLOCATE(heff (ncoup*n_k1*nsub,ncoup*n_k1*nsub))
  heff =0.0D0
  ! get similarity transformed interaction
  CALL vlowk_mtx(ncoup,nsub,vzz,heff,i)
  !CALL vlowk_nochange(ncoup*nsub, vzz, heff)

  !     make now transformation to h.o. basis in the rel and cm system
  !     loop over all cm and rel coordinate configurations

  ! *** Implement loop over subchannels and use relative mass in HO function
  nconf = rel_conf%nconfs_rel(i)
  DO subab = 1, nsub
  mab = mass(i,subab)
  DO subcd = 1, nsub
  mcd = mass(i,subcd)
        IF (subab == subcd) THEN
            ida = sub_id(i,subab)/100
            idb = sub_id(i,subab) - ida*100
            lookup = coulomb_lookup(ida,idb)
        ENDIF
  DO bra =1, rel_conf%nconfs_rel(i)
     k1=bra
     a=rel_conf%rel_ab(i,k1)
     n=relcm_sp_data%nrel(a)
     la=relcm_sp_data%lrel(a)
     IF ((n+n+la) > nlmax) CYCLE
     DO ket =1,  rel_conf%nconfs_rel(i) 
        k2=ket
        c=rel_conf%rel_ab(i,k2)
        np=relcm_sp_data%nrel(c)        
        lb=relcm_sp_data%lrel(c)
        !  No dependence of the bare interaction upon the CoM momenta
        !  The Hamiltonian is also diagonal in L and N 
        vsum = 0.
        IF ((np+np+lb) > nlmax) CYCLE
        IF ( ncoup == 1) THEN
           lim1=1; lim2=n_k1 ; lim3=1 ; lim4=n_k1
        ELSEIF ( ncoup == 2 ) THEN
           IF ( (la == lb).AND. ( jang > la) ) THEN
              lim1=1; lim2=n_k1
              lim3=1 ; lim4=n_k1
           ELSEIF ( (la == lb).AND. ( jang < la) ) THEN
              lim1=1+n_k1*nsub; lim2=n_k1+n_k1*nsub
              lim3=1+n_k1*nsub ; lim4=n_k1+n_k1*nsub
           ELSEIF ( la >  lb ) THEN
              lim1=1+n_k1*nsub; lim2=n_k1+n_k1*nsub
              lim3=1 ; lim4=n_k1
           ELSEIF ( la <  lb ) THEN
              lim1=1; lim2=n_k1
              lim3=1+n_k1*nsub ; lim4=n_k1+n_k1*nsub
           ENDIF
        ENDIF
        lim1 = lim1 + (subab-1)*n_k1
        lim2 = lim2 + (subab-1)*n_k1
        lim3 = lim3 + (subcd-1)*n_k1
        lim4 = lim4 + (subcd-1)*n_k1
        !WRITE(*,*) lim1, lim2, lim3, lim4

        CALL vlowk_rnl(n,la,mab,rnl_ab)
        CALL vlowk_rnl(np,lb,mcd,rnl_cd)
        CALL vlowk_hosc(rnl_cd, rnl_ab,REAL(heff(lim1:lim2,lim3:lim4)),vsum)
        !       Here we set up the Coulomb part in an oscillator basis. 
        e_coulomb = 0.
        IF ((la == lb) .AND. (lookup /= 0)) THEN
           e_coulomb = coulomb_relcom(n,np,la,lookup)
        ENDIF
        !  only potential energy : coulomb + V_NN
        vint(ket+(subcd-1)*nconf,bra+(subab-1)*nconf)=vsum+e_coulomb
     ENDDO
  ENDDO
  ENDDO
  ENDDO
  DEALLOCATE(heff); DEALLOCATE(vzz)

END SUBROUTINE vlowk_channel

SUBROUTINE vlowk_transform_nl(i, heff, nh, vrel)
  USE wave_functions
  USE relcm_gmatrix
  USE partial_waves
  USE configurations
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i, nh
  INTEGER :: ncoup, n,np, bra, ket, k1, k2, a, c, subab, subcd
  INTEGER :: la, lb, jang, lim1, lim2, lim3, lim4, nsub, nconf 
  REAL(KIND = 8) :: e_coulomb, vsum, mab, mcd
  COMPLEX*16, INTENT(IN) :: heff(nh,nh)
  REAL(KIND = 8),  INTENT(INOUT):: &
       vrel(rel_conf%nconfs_rel(i)*no_sub(i),rel_conf%nconfs_rel(i)*no_sub(i))
  REAL(KIND(1.0D0)), ALLOCATABLE :: rnl_ab(:), rnl_cd(:)

  jang=jang_rel(i)
  ncoup=1
  IF ( orb_lrel_max(i) == jang_rel(i) ) ncoup = 1
  IF ( (orb_lrel_max(i) /= jang_rel(i)) .AND. jang /= 0 ) ncoup = 2

  nsub = no_sub(i)
  nconf = rel_conf%nconfs_rel(i)
  DO subab = 1, nsub
  mab = mass(i,subab)
  DO subcd = 1, nsub
  mcd = mass(i,subcd)
  DO bra =1, rel_conf%nconfs_rel(i)
     k1=bra
     a=rel_conf%rel_ab(i,k1)
     n=relcm_sp_data%nrel(a)
     la=relcm_sp_data%lrel(a)
     IF ((n+n+la) > nlmax) CYCLE
     DO ket =1,  rel_conf%nconfs_rel(i) 
        k2=ket
        c=rel_conf%rel_ab(i,k2)
        np=relcm_sp_data%nrel(c)        
        lb=relcm_sp_data%lrel(c)
        !  No dependence of the bare interaction upon the CoM momenta
        !  The Hamiltonian is also diagonal in L and N 
        vsum = 0.
        IF ((np+np+lb) > nlmax) CYCLE
        IF ( ncoup == 1) THEN
           lim1=1; lim2=n_k1 ; lim3=1 ; lim4=n_k1
        ELSEIF ( ncoup == 2 ) THEN
           IF ( (la == lb).AND. ( jang > la) ) THEN
              lim1=1; lim2=n_k1
              lim3=1 ; lim4=n_k1
           ELSEIF ( (la == lb).AND. ( jang < la) ) THEN
              lim1=1+n_k1*nsub; lim2=n_k1*nsub+n_k1*nsub
              lim3=1+n_k1*nsub ; lim4=n_k1*nsub+n_k1*nsub
           ELSEIF ( la >  lb ) THEN
              lim1=1+n_k1*nsub; lim2=n_k1*nsub+n_k1*nsub
              lim3=1 ; lim4=n_k1
           ELSEIF ( la <  lb ) THEN
              lim1=1; lim2=n_k1
              lim3=1+n_k1*nsub ; lim4=n_k1*nsub+n_k1*nsub
           ENDIF
        ENDIF
        lim1 = lim1 + (subab-1)*n_k1
        lim2 = lim2 + (subab-1)*n_k1
        lim3 = lim3 + (subcd-1)*n_k1
        lim4 = lim4 + (subcd-1)*n_k1
        CALL vlowk_rnl(n,la,mab,rnl_ab)
        CALL vlowk_rnl(np,lb,mcd,rnl_cd)
        CALL vlowk_hosc(rnl_cd, rnl_ab,REAL(heff(lim1:lim2,lim3:lim4)),vsum)
        !       Here we set up the Coulomb part in an oscillator basis. 
        e_coulomb = 0.
        ! What to do here when bra is idp, but ket isn't or viceversa.
        IF ( ( la == lb ).AND.( iso(i) == -1 ))  THEN
           !e_coulomb = coulomb_relcom(n,np,la) 
        ENDIF
        !  only potential energy : coulomb + V_NN
        vrel(ket+(subcd-1)*nconf,bra+(subab-1)*nconf)=vsum+e_coulomb
     ENDDO
  ENDDO
  ENDDO
  ENDDO

END SUBROUTINE vlowk_transform_nl
SUBROUTINE vlowk_nochange(coup,vzz,heff)
  USE wave_functions
  USE partial_waves
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: coup
  COMPLEX*16, DIMENSION(coup*n_k,coup*n_k), INTENT(IN) :: vzz
  COMPLEX*16, DIMENSION(coup*n_k1,coup*n_k1), INTENT(INOUT) :: heff
  INTEGER :: i, j, k, l

  DO k = 0, (coup-1)
    DO l= 0, (coup-1)
      DO i = 1, n_k1
        DO j = 1, n_k1
          heff(i+k*n_k1, j+l*n_k1) = vzz(i+k*n_k, j+l*n_k)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE vlowk_nochange

!
!       Compute the T(G)-mtx for the actual channel 
!
SUBROUTINE vlowk_mtx(ncoup,nsub,vzz,heff,ichan)
  USE wave_functions
  USE partial_waves
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ncoup, ichan, nsub
  COMPLEX*16, DIMENSION(ncoup*n_k*nsub,ncoup*n_k*nsub), INTENT(IN) :: vzz
  COMPLEX*16, DIMENSION(ncoup*n_k1*nsub,ncoup*n_k1*nsub), INTENT(INOUT) :: heff
  INTEGER :: i, j, ntot, i1, i2, np, nq, j1, tisoz 
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: ham, cvec, temp
  COMPLEX*16, ALLOCATABLE, DIMENSION(:) ::  ceig
  REAL(KIND(1.0D0)) :: imass(nsub)

    !WRITE(*,*) "vlowk_mtx called with channel:", ichan
  imass = mass(ichan,:)
  ! dimension of vectors and matrices
  ntot = ncoup*n_k*nsub
  ALLOCATE( ham(ntot, ntot), cvec(ntot,ntot), ceig(ntot), temp(ntot,ntot)) 
  ! setup hamiltonian to be diagonalized
  ham = DCMPLX(0.0D0,0.0D0) 
  cvec = dcmplx(0.0d0, 0.0d0)
  ceig = dcmplx(0.0d0, 0.0d0)

  CALL complex_hamiltonian(ham,vzz, ntot, nsub, ichan)
  temp = ham
  CALL vlowkdiag_exact( temp, ntot, cvec, ceig )
  !
  ! construct renormalized nn-interaction via the Lee-Suzuki sim. transform
  !
  ! size of model space     : np = n_k1*ncoup*nsub
  ! size of complement space: nq = n_k2*ncoup*nsub
  np = n_k1*ncoup*nsub; nq = n_k2*ncoup*nsub
  CALL effective_int( np, nq, ntot, ncoup, nsub, cvec, ceig, heff )   
  !
  ! subtract diagonal elements to obtain Veff, 
  ! and compare exact with effective interactions: 
  !
  tisoz = iso(ichan) 
  i2 = 0
  DO i = 1, np
     i1 = i
     IF ( i > n_k1 ) THEN
        i1 = MOD(i,n_k1)
        IF (i1 == 0) i1 = n_k1
     ENDIF
     IF (i1 == 1) THEN
        i2 = MOD(i2 + 1,nsub)
        IF (i2 == 0) i2 = nsub
     ENDIF
     heff(i,i) = heff(i,i) - ( krel(i1) * krel(i1) )/(2*imass(i2))
     DO j = 1, np
        j1 = j
        IF ( j > n_k1 ) THEN
            j1 = MOD(j,n_k1)
            IF(j1 == 0) j1 = n_k1
        ENDIF
        heff(i,j) = heff(i,j)/SQRT( wkrel(i1)*wkrel(j1) )/krel(i1)/krel(j1)
     ENDDO
  ENDDO
  DEALLOCATE(ham, temp, cvec );   DEALLOCATE(ceig)

END SUBROUTINE vlowk_mtx

SUBROUTINE complex_hamiltonian(h,vzz,ntot,nsub,channel)
    USE wave_functions
    USE constants
    USE partial_waves
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ntot, channel, nsub
    COMPLEX*16, DIMENSION(ntot,ntot), INTENT(INOUT) :: h
    COMPLEX*16, DIMENSION(ntot,ntot), INTENT(IN) :: vzz
    REAL*8 :: imass(nsub)
    INTEGER :: i, j, i1, i2, j1

    imass = mass(channel,:)
    h = 0.
    i2 = 0
    DO i = 1, ntot 
        i1 = i
        IF ( i > n_k ) THEN
            i1 = MOD(i,n_k)
            if (i1 == 0 ) i1 = n_k
        ENDIF
        IF (i1 == 1) THEN 
            i2 = MOD(i2 + 1, nsub)
            IF (i2 == 0) i2 = nsub
        ENDIF
        h(i,i) = ( krel(i1) * krel(i1) )/(2*imass(i2))  + & 
            (krel(i1)) * (krel(i1)) * wkrel(i1) * vzz(i,i) 
        DO j = 1, ntot 
            j1 = j
            IF ( j > n_k ) THEN
                j1 = MOD(j,n_k)
                if (j1 == 0 ) j1 = n_k
            ENDIF
            IF (i /= j ) THEN
                h(i,j) =  SQRT(wkrel(j1) * wkrel(i1)) * krel(j1) * krel(i1) * vzz(i,j)   
            ENDIF
        ENDDO
    ENDDO

END SUBROUTINE complex_hamiltonian

!
! calculate effective interaction for gamow shell model calcs
!
SUBROUTINE effective_int( np, nq, ntot, ncoup, nsub, cvec, ceig, heff )   
  USE constants
  USE wave_functions
  IMPLICIT NONE
  INTEGER, INTENT(in) ::  np, nq, ntot, ncoup, nsub
  INTEGER :: k1,i,j,k, coup
  INTEGER ::  orb
  COMPLEX*16, DIMENSION(ntot,ntot), INTENT(in) :: cvec
  COMPLEX*16, DIMENSION(ntot), INTENT(in) :: ceig
  COMPLEX*16, DIMENSION(np,np), INTENT(out) :: heff
  COMPLEX*16 , ALLOCATABLE :: cvec_pp(:,:), cvec_qp(:,:)
  COMPLEX*16 , ALLOCATABLE:: ceig_p(:)
  REAL(kind = 8) , ALLOCATABLE:: overlap(:)
  INTEGER, ALLOCATABLE :: orbit_no(:)

    ALLOCATE(cvec_pp(np,np), cvec_qp(nq,np))
    ALLOCATE(ceig_p(np))
    ALLOCATE(overlap(ntot))
    ALLOCATE(orbit_no(ntot))
  !WRITE(*,*) "effective_int called..."
    coup = ncoup*nsub
  !WRITE(*,*) "coupled: ", coup
    
  ! calculate P-space overlap of all eigenvectors 
  ! loop over all eigen vectors

    DO orb = 1, ntot
        overlap(orb) = 0.
        DO i = 1, n_k1
            DO j = 0, (coup-1)
                overlap(orb) = overlap(orb) + ABS( cvec(i+j*n_k,orb) )**2
            ENDDO
        ENDDO
        orbit_no(orb) = orb
    ENDDO

  ! sort all overlaps and corresponding orbit numbers 
  CALL eigenvalue_sort(overlap,orbit_no, ntot)
  !WRITE(*,*) "eigenvalue_sort returned..."
  cvec_pp = 0.; cvec_qp = 0.

    DO i = 1, np
        k1 = orbit_no(i)
        ceig_p(i) = ceig(k1)
        DO j = 1, n_k
            DO k = 0, (coup-1)
                IF (j <= n_k1) THEN
                    cvec_pp(j+k*n_k1,i) = cvec(j+k*n_k,k1)
                ELSEIF (j > n_k1) THEN
                    cvec_qp(j+(k-1)*n_k1,i) = cvec(j+k*n_k,k1)
                ENDIF
            ENDDO
        ENDDO
    ENDDO

  heff = 0.
  !WRITE(*,*) "Before calling lee_suzuki2...", np, nq
  !IF (coup/=8) RETURN

  CALL lee_suzuki2( cvec_pp, cvec_qp, ceig_p, np, nq, heff )
  !WRITE(*,*) "After calling lee_suzuki2..."
  DEALLOCATE(cvec_pp, cvec_qp, ceig_p, overlap, orbit_no)

END SUBROUTINE effective_int

!
! eigenvalue sort
! sort cmplx vector real(a(1))<real(a(2)) < ... < real(a(n))
!
SUBROUTINE eigenvalue_sort(A,B, n)
  IMPLICIT NONE
  INTEGER :: i, j, n
  REAL(kind = 8), DIMENSION(n), INTENT(INOUT) :: A
  INTEGER, DIMENSION(n), INTENT(inout) :: B
  REAL(kind = 8) :: temp
  INTEGER :: orb

  !WRITE(*,*) "eigenvalue_sort called..."
  DO i = 1, n
     DO j = 1, n
        IF ( ABS( A(i) )  > ABS( A(j) ) ) THEN
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
  COMPLEX*16, DIMENSION(np,np), INTENT(INOUT) ::  heff
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: temp, omega2, sim1, sim2, &
    omega_2inv, u, u_inv, cvec_pp_inv, eigen_vec, vl, heff_rhs, omega, temp1
  COMPLEX*16, ALLOCATABLE, DIMENSION(:):: ceig_p, omega2_eig
  REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: rwork
  COMPLEX*16, DIMENSION(10000) :: work1
  COMPLEX*16 :: d, sum1, determinant
  INTEGER :: i_p,j_p, k, i_q 
  INTEGER :: lda, ldvl, ldvr, info, lwork
  CHARACTER*1 :: jobvl, jobvr, balanc, sense
  REAL(KIND = 8) :: a1, a2, b1, b2

  !WRITE(*,*) "lee_suzuki2 - np, nq:", np, nq
  ALLOCATE(rwork(2*np))
  ALLOCATE(ceig_p(np), omega2_eig(np))
  ALLOCATE(temp(np,np), omega2(np,np), sim1(np,np), sim2(np,np))
  ALLOCATE(omega_2inv(np,np), u(np,np), u_inv(np,np), cvec_pp_inv(np,np))
  ALLOCATE(eigen_vec(np,np), vl(np,np), heff_rhs(np,np))
  ALLOCATE(omega(nq,np), temp1(np,np))

  balanc = 'n';  jobvl = 'n' ;  jobvr = 'v';  sense = 'n';  lda = np
  ldvl = 1;  ldvr = np;  lwork = 10000
  eigen_vec = 0. 
  temp1 = TRANSPOSE(cvec_pp) 
  CALL zgeev( jobvl, jobvr, np, temp1, lda, omega2_eig, vl, ldvl, eigen_vec, ldvr, &
       work1, lwork, rwork, info )
  determinant = PRODUCT(omega2_eig(:))
  !write(*,*) 'check determinant', determinant
  ! the P->Q space transformation matrix, omega 
  cvec_pp_inv = cvec_pp

  !write(*,*) "cmplxmatinv called..."
  CALL cmplxmatinv(cvec_pp_inv, np, d)
  !write(*,*) "cmplxmatinv returned..."
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
  !write(*,*) "Before calling sqrtmat..."
  CALL sqrtmat( omega2, U, U_inv ,np ) 
  !write(*,*) "sqrtmat returned..."
  heff_rhs =  MATMUL( U, MATMUL( heff,u_inv ) ) 
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
  !write(*,*) "Before calling zgeev..."
  CALL zgeev( jobvl, jobvr, np, heff_rhs, lda, ceig_p, vl, ldvl, eigen_vec, ldvr, &
       work1, lwork, rwork, info )
  !write(*,*) "zgeev returned..."
  !   compare spectrum from exact and P-space diagonalization
  WRITE(6,*) 'Compare model space two-body spectrum with exact spectrum:' 
  DO i_p = 1, np
     a1 = REAL( ceig_p(i_p))
     a2 = AIMAG( ceig_p(i_p))
     b1 = REAL( ceig(i_p) )
     b2 = AIMAG( ceig(i_p) ) 
     WRITE(6,*) A1, B1
  ENDDO

  DEALLOCATE(rwork)
  DEALLOCATE(ceig_p, omega2_eig)
  DEALLOCATE(temp, omega2, sim1, sim2)
  DEALLOCATE(omega_2inv, u, u_inv, cvec_pp_inv)
  DEALLOCATE(eigen_vec, vl, heff_rhs)
  DEALLOCATE(omega, temp1)

END SUBROUTINE lee_suzuki2
!
! eigenvalue sort
! sort cmplx vector real(a(1))<real(a(2)) < ... < real(a(n))
!
SUBROUTINE eigenvalue_sort_cmplx(A, n)
  IMPLICIT NONE
  INTEGER :: i, j, n
  COMPLEX*16, DIMENSION(n), INTENT(INOUT) :: A
  COMPLEX*16 :: temp

  DO i = 1, n
     DO j = 1, n
        IF ( ABS( A(i) )  > ABS( A(j) ) ) THEN
           temp = A(i)
           A(i) = A(j) 
           A(j) = temp
        END IF
     END DO
  END DO

END SUBROUTINE eigenvalue_sort_cmplx

SUBROUTINE vlowkdiag_exact( h, n, cvec, ceig )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: h
  INTEGER :: i,j, i1, kvec, lwork, option, il ,iu, info
  REAL(kind = 8) :: A(n*(n+1)), work(300*n)
  REAL(kind = 8) :: thresh
  COMPLEX*16 :: cvl, cvu
  LOGICAL :: flag(n)  
  COMPLEX*16 , INTENT(OUT):: ceig(n)
  COMPLEX*16 , INTENT(OUT) :: cvec(n,n)

  !WRITE(*,*) "vlowkdiag_exact called..."
  lwork = 300*n
  thresh = 30.0
  kvec = n
  option = 4
  info = 1
  i1 = 0
  DO i =  1, n
     DO j =  1, i
        i1 = i1 + 1
        a(i1) = DBLE(h(j,i))
        a(i1+n*(n+1)/2) = AIMAG(h(j,i))
     END DO
  END DO
  CALL cs(n,a,ceig,kvec,cvec,lwork,work,thresh, option,il,iu,cvl,cvu,flag,info)

END SUBROUTINE vlowkdiag_exact

!
!                 Set up h.o. wf for rel system
!
SUBROUTINE vlowkho_wfunction
  USE constants
  USE wave_functions
  IMPLICIT NONE
  INTEGER :: n, l, i
  REAL(KIND = 8)  :: cx(0:200), factor, z_rel, xp, ph, oscl_r, sum_rel, contrib

  oscl_r=oscl*SQRT(2.)        ! Oscillator parameter for relative
  DO n=0,nmax
     ph=(-1.D0)**n
     DO l=0,lmax
        sum_rel=0.
        factor = 0.5D0*((n+l+2)*LOG(2.D0)+fac(n)-dfac(2*n+2*l+1)-0.5D0*LOG(pi))
        factor = EXP(factor)
        DO i=1,n_rel
           z_rel= ra(i)*ra(i)*oscl_r*oscl_r
           CALL laguerre_general( n, l+0.5D0, z_rel, cx )
           xp = EXP(-z_rel*0.5)*((ra(i)*oscl_r)**l)*cx(n)
           rnlr(i,l,n) = xp*(wra(i)*(ra(i)**2))*ph*factor*(oscl_r**(1.5D0))     ! rel wf
           contrib = xp*ph*factor*(oscl_r**(1.5D0)) 
           sum_rel=sum_rel+ wra(i)*(contrib*ra(i))**2
        ENDDO
        WRITE(6,*) 'Norm rel ho wf n,l : ', n, l, sum_rel
     ENDDO
  ENDDO

END SUBROUTINE vlowkho_wfunction

!
!                 Set up h.o. wf for rel system
!
SUBROUTINE vlowk_rnl(n,l,m, rnlm)
    USE constants
    USE wave_functions
    IMPLICIT NONE
    INTEGER :: n, l, i
    REAL(KIND = 8)  :: cx(0:200), factor, z_rel, xp, ph, oscl_r, sum_rel, contrib,m
    REAL(KIND(1.0D0)), INTENT(INOUT):: rnlm(n_rel)

    oscl_r =  hbarc/SQRT(m*hbar_omega)
    ph=(-1.D0)**n
    sum_rel=0.
    factor = 0.5D0*((n+l+2)*LOG(2.D0)+fac(n)-dfac(2*n+2*l+1)-0.5D0*LOG(pi))
    factor = EXP(factor)
    DO i=1,n_rel
        z_rel= ra(i)*ra(i)*oscl_r*oscl_r
        CALL laguerre_general( n, l+0.5D0, z_rel, cx )
        xp = EXP(-z_rel*0.5)*((ra(i)*oscl_r)**l)*cx(n)
        rnlm(i) = xp*(wra(i)*(ra(i)**2))*ph*factor*(oscl_r**(1.5D0))     ! rel wf
        contrib = xp*ph*factor*(oscl_r**(1.5D0)) 
        sum_rel=sum_rel+ wra(i)*(contrib*ra(i))**2
    ENDDO
    !WRITE(6,*) 'Norm rel ho wf n,l,m: ', n, l, m, sum_rel

END SUBROUTINE vlowk_rnl
!
!          Obtain the bare potential in oscillator basis
!          Final V is in units of MeV or eV
!          Time consuming part. Parallel version available.
!
SUBROUTINE vlowk_hosc(wave_bra, wave_ket, a,vsum)
  USE wave_functions
  USE constants
  USE relcm_gmatrix
  IMPLICIT NONE
  REAL(KIND = 8), INTENT(INOUT) :: vsum
  INTEGER :: i, j
  REAL(KIND = 8), DIMENSION(n_k1,n_k1), INTENT(IN) :: a
  REAL(KIND = 8) :: sum1, sum2, hbarc3
  REAL(KIND = 8), DIMENSION(n_rel), INTENT(IN) :: wave_bra, wave_ket 

  hbarc3=hbarc**3
  sum1=0.
  DO i=1, n_k1
     sum2=0.
     DO j=1, n_k1
        sum2=sum2+wave_ket(j)*a(j,i)
     ENDDO
     sum1=sum1+sum2*wave_bra(i)
  ENDDO
  vsum = sum1*hbarc**3

END SUBROUTINE vlowk_hosc

! *************** Added from gmatrix

SUBROUTINE mosh_transf(it,ij,str,gmatrix_configs,lab_to_relcoeff, &
     lab_to_relconf,lab_to_relnumber,max_coeff)
  USE single_particle_orbits
  USE configurations
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN) :: gmatrix_configs
  INTEGER, INTENT(IN) :: it, ij, max_coeff, str
  INTEGER :: k1,k2,i,j,n1,l1,j1d,n2,l2,j2d,k, number, id1, id2
  REAL(DP), DIMENSION(gmatrix_configs%number_confs,max_coeff), & 
       INTENT(INOUT) :: lab_to_relcoeff
  INTEGER, DIMENSION(gmatrix_configs%number_confs,max_coeff), &
       INTENT(INOUT) :: lab_to_relconf 
  INTEGER, DIMENSION(gmatrix_configs%number_confs), &
       INTENT(INOUT) :: lab_to_relnumber 
  REAL(DP), DIMENSION(max_coeff) :: coeb
  INTEGER, DIMENSION(max_coeff) :: mqnb

  DO k=1,gmatrix_configs%number_confs              
     k2=k*2
     k1=k2-1
     i=gmatrix_configs%config_ab(k1)
     j=gmatrix_configs%config_ab(k2)
     n1=all_orbit%nn(i) 
     l1=all_orbit%ll(i)
     j1d=all_orbit%jj(i)
     id1=all_orbit%ptype(i)
     n2=all_orbit%nn(j)
     l2=all_orbit%ll(j)
     j2d=all_orbit%jj(j)
     id2=all_orbit%ptype(j)
     CALL transf_mbs(n1,l1,j1d,id1,n2,l2,j2d,id2,it,ij,str,mqnb,coeb,max_coeff,number)
     lab_to_relcoeff(k,:)=coeb(:)
     lab_to_relconf(k,:)=mqnb(:)
     lab_to_relnumber(k)=number
  ENDDO

END SUBROUTINE  mosh_transf
!
!           Explicit evaluation of transf coeff
!
SUBROUTINE transf_mbs(n1,l1,j1d,id1,n2,l2,j2d,id2,itz,ij,str,mqnb,coeb,max_coeff,number)
  USE single_particle_orbits
  USE partial_waves
  USE ang_mom_functions
  USE configurations
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n1,l1,j1d,id1,n2,l2,j2d,id2,itz,ij,str,max_coeff
  INTEGER, DIMENSION(max_coeff), INTENT(OUT) :: mqnb
  REAL(DP), DIMENSION(max_coeff), INTENT(OUT) :: coeb
  INTEGER :: ld, l, n, lcd, lambd, isign, lambdd,  a, b, k1, k2, &
       l2d, l1d, nemx, ijd, l_min, l_max, &
       jlsd, nch, isd, jls, ispin,nc, lc, nconfs
  INTEGER, INTENT(OUT) :: number
  REAL(DP) :: bk, w6, w9, s , co, cst, delta,mass_ab
  LOGICAL triag

  mqnb=0 ; coeb=0.
  l1d=l1*2 ;   l2d=l2*2 ;   ijd=ij*2 ;   nemx=2*n1+l1+2*n2+l2
  mass_ab = sp_mass(id1)/sp_mass(id2)
  number=0
  !     loop over possible partial wave channels which specify J, l, S and Tz
  DO nch=1,no_channels
     IF ( itz /= iso(nch) ) CYCLE
     IF ( str /= strange(nch) ) CYCLE
     l_min=orb_lrel_min(nch)
     l_max=orb_lrel_max(nch)
     ispin=spin_rel(nch) ; isd=2*ispin
     jls=jang_rel(nch) ; jlsd=2*jls
     !     loop over rel and cm configs specifying n, l, N and L
     DO nconfs =1, relcm_conf%nconfs_relcm(nch)
        k2=nconfs*2
        k1=k2-1
        b=relcm_conf%relcm_ab(nch,k2)
        a=relcm_conf%relcm_ab(nch,k1)
        n=relcm_sp_data%nrel(a)
        l=relcm_sp_data%lrel(a)
        IF ( (l /= l_max) .AND. ( l /= l_min ) ) CYCLE
        nc=relcm_sp_data%nrel(b)
        lc=relcm_sp_data%lrel(b)
        !     pauli test for identical particles in partial waves
        IF( (id1==id2) .AND. ( (-1)**(l+ispin+itz/2)>0) ) CYCLE
        !IF((ABS(itz) == 1).AND. ((-1)**(l+ispin+1) > 0)) CYCLE
        !     get right sign from isospin clebsch-gordan coeff for the partial waves
        !     since we have for Tz=0 always a coupling order | (pn) J >
        !     this means that the T=0 contrib goes like  -(1-(-)^(l+spin))/sqrt(2)
        !     and l+s have to be odd       
        cst = 1.0_dp
        IF ( id1==id2 ) THEN
           cst = SQRT(2.D0)/SQRT( 1.d0+ delta(n1,n2)*delta(l1,l2)*delta(j1d,j2d) )
        !ELSEIF ( (itz == 0) ) THEN
        ELSEIF ( id1 /= id2) THEN
           IF ( (MOD(l+ispin,2) /= 0 ) )  cst= -1.0_dp
        ENDIF
        IF(triag(jls,ij,lc)) CYCLE
        IF((-1)**(l1+l2+l+lc) < 0) CYCLE
        IF((2*n+l+2*nc+lc) /= nemx) CYCLE
        IF((2*n+l+2*nc+lc) > nlmax) CYCLE
        IF(triag(jls,l,ispin)) CYCLE
        ld=l*2 ; lcd=lc*2
        co=0.D0
        DO lambd = ABS(l1-l2), l1+l2
           lambdd=lambd*2
           IF(triag(lambd,l,lc)) CYCLE
           IF(triag(ij,lambd,ispin)) CYCLE
           bk= gmosh(n,l,nc,lc,n1,l1,n2,l2,lambd,mass_ab)
           IF( ABS (bk) == 0.D0) CYCLE
           w9=snj(l1d,1,j1d,l2d,1,j2d,lambdd,isd,ijd)
           IF ( ABS( w9) == 0.D0) CYCLE
           w9=w9*SQRT((j1d+1.)*(j2d+1.)*(lambdd+1.)*(isd+1.))
           isign=(lcd+ld+ijd+isd)/2
           s=1.-2.*MOD(isign,2)
           w6=s*sjs(lcd,ld,lambdd,isd,ijd,jlsd)
           w6=w6*SQRT(FLOAT((lambdd+1)*(jlsd+1)))
           co=co+cst*bk*w6*w9*(-1.)**(l+lambd+jls+ij)
        ENDDO
        IF(DABS(co) >= 1.D-5)THEN
           number=number+1
           mqnb(number)=nconfs*1000+nch
           coeb(number)=co
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE transf_mbs
!
!           Finds total number of transformation  coeffs 
!           for given J, parity and isospin projection
!
SUBROUTINE find_ncoeffs(it,ij,str,this,max_coeff)
  USE single_particle_orbits
  USE configurations
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN) :: this
  INTEGER, INTENT(IN) :: it, ij, str
  INTEGER :: k1,k2,i,j,n1,l1,n2,l2,number, k, id1, id2
  INTEGER, INTENT(OUT) :: max_coeff
  max_coeff=0
  DO k=1, this%number_confs              
     k2=k*2
     k1=k2-1
     i=this%config_ab(k1)
     j=this%config_ab(k2)
     n1=all_orbit%nn(i) 
     l1=all_orbit%ll(i)
     n2=all_orbit%nn(j)
     l2=all_orbit%ll(j)
     id1=all_orbit%ptype(i)
     id2=all_orbit%ptype(j)
     CALL count_trans_mbs(id1,id2,n1,l1,n2,l2,it,ij,str,number)
     IF ( max_coeff < number ) max_coeff = number
  ENDDO

END SUBROUTINE find_ncoeffs
!
!           counts max possible number  of transf coeff
!
SUBROUTINE count_trans_mbs(id1, id2, n1,l1,n2,l2,itz,ij,str,number)
  USE single_particle_orbits
  USE partial_waves
  USE configurations
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: number
  INTEGER, INTENT(IN) :: n1,l1,n2,l2,itz,ij, id1, id2, str
  INTEGER :: l, n, nemx, l_min, l_max, nconfs , nch, jls, & 
       ispin, nc, lc, a, b,  k1, k2 
  LOGICAL triag

  !WRITE(*,*) 'Called count_trans_mbs: id1, n1, l1, id2, n2, l2, itz, ij, str'
  !WRITE(*,'(9I4)') id1, n1, l1, id2, n2, l2, itz, ij, str

  nemx=2*n1+l1+2*n2+l2
  number=0
  DO nch=1,no_channels
     !WRITE(*,*) "count_trans_mbs, nch:", nch
     IF ( itz /= iso(nch) ) CYCLE
     IF ( str /= strange(nch) ) CYCLE
     l_min=orb_lrel_min(nch)
     l_max=orb_lrel_max(nch)
     ispin=spin_rel(nch)
     jls=jang_rel(nch)
     DO nconfs =1, relcm_conf%nconfs_relcm(nch)
        !WRITE(*,*) "Starting conf: ", nconfs
        k2=nconfs*2
        k1=k2-1
        !WRITE(*,*) "Got k1 and k2...", k1, k2
        b=relcm_conf%relcm_ab(nch,k2)
        a=relcm_conf%relcm_ab(nch,k1)
        !WRITE(*,*) "Got a and b...", a, b
        n=relcm_sp_data%nrel(a)
        l=relcm_sp_data%lrel(a)
        !WRITE(*,*) "Got n and l...", n, l
        IF ( (l /= l_max) .AND. ( l /= l_min ) ) CYCLE
        nc=relcm_sp_data%nrel(b)
        lc=relcm_sp_data%lrel(b)
        !WRITE(*,*) "Got nc and lc...", nc, lc
        IF((id1==id2).AND.((-1)**(l+ispin+itz/2) > 0)) CYCLE
        IF(triag(jls,ij,lc)) CYCLE
        IF((-1)**(l1+l2+l+lc) < 0) CYCLE
        IF((2*n+l+2*nc+lc) /= nemx) CYCLE
        IF((2*n+l+2*nc+lc) > nlmax) CYCLE
        IF(triag(jls,l,ispin)) CYCLE
        number=number+1
        !WRITE(*,*) "Done with conf, number: ", nconfs, number
     ENDDO
     !WRITE(*,*) "Done with ch: ", nch
  ENDDO

END SUBROUTINE count_trans_mbs
!
!           Set up the effective interaction in the lab-frame using the
!           mtx elements in the rel-CoM frame
!
SUBROUTINE final_vlowk_labsystem
  USE single_particle_orbits  
  USE configurations
  USE constants
  USE relcm_gmatrix
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: gmatrix_configs 
  REAL(KIND = 8), ALLOCATABLE :: lab_to_relcoeff(:,:), reduced_p1(:,:)
  INTEGER, ALLOCATABLE :: lab_to_relconf(:,:)
  INTEGER, ALLOCATABLE :: lab_to_relnumber(:)
  INTEGER ::  p_parity, ang_mom, isospin_z, max_coeff, pq_confs, istrange
  COMPLEX*16, ALLOCATABLE::  gna(:,:)
  REAL(KIND=8), ALLOCATABLE:: twobody_com(:,:), twobody_r2(:,:), &
    twobody_p2(:,:)

  ! momentum space representation
  ! *** This is for all combinations of particles
  ALLOCATE( reduced_p1(all_orbit%total_orbits, all_orbit%total_orbits) )
  CALL reduced_matel_p1_osc(reduced_p1)

    !WRITE(6,*) "Strange  isospin  parity  ang_mom"
DO istrange = strange_min, strange_max
    itzmin = istrange-2; itzmax = -1*itzmin
  !     loop over isospin projection
  DO isospin_z=itzmin,itzmax,2
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=j_lab_min,j_lab_max

           !     find all possible configurations, large space and model space 
           CALL  number_gmatrix_confs&
                (ang_mom,p_parity,isospin_z,istrange,gmatrix_configs)
           IF (gmatrix_configs%number_confs <= 0 ) CYCLE
           pq_confs=gmatrix_configs%number_confs
           !WRITE(*,*) istrange, isospin_z, p_parity, ang_mom
           !WRITE(*,*) "Allocating gmatrix_configs%config_ab..." , pq_confs
           ALLOCATE(gmatrix_configs%config_ab(pq_confs+pq_confs) )
           !WRITE(*,*) "Calling setup_gmatrix_configurations..."
           CALL setup_gmatrix_configurations &
                (ang_mom,p_parity,isospin_z,istrange,gmatrix_configs)
           !     find max possible number of transformation coeffs

           !WRITE(*,*) "Calling find_ncoeffs..."
           CALL find_ncoeffs(isospin_z,ang_mom,istrange,gmatrix_configs,max_coeff)
           !     allocate space for various arrays needed in transform rel-cm -> lab 
           ALLOCATE(lab_to_relconf(pq_confs,max_coeff), &
                lab_to_relcoeff(pq_confs,max_coeff))
           ALLOCATE(lab_to_relnumber(pq_confs))
           ALLOCATE(twobody_p2(pq_confs, pq_confs))
           ALLOCATE(twobody_r2(pq_confs, pq_confs))
           ALLOCATE(twobody_com(pq_confs, pq_confs))
           !     setup transformation coefficients for oscillator basis
           !     transformations from the c.m. frame to the lab frame
           !WRITE(*,*) "Calling mosh_transf..."
           CALL mosh_transf(isospin_z,ang_mom,istrange,gmatrix_configs, &
                lab_to_relcoeff, lab_to_relconf, lab_to_relnumber, &
                max_coeff)
           ALLOCATE(gna(pq_confs, pq_confs))
           gna=0.0D0; twobody_p2 = 0.0D0; twobody_r2 = 0.0D0
           twobody_com = 0.0D0
           !     Performs HO transformation from rel and cm coordinates to
           !     lab system.
           !WRITE(*,*) "Calling vlowk_free..."
           CALL vlowk_free(gna,lab_to_relcoeff, lab_to_relconf,&
                lab_to_relnumber,max_coeff,gmatrix_configs)

           ! find pipj elements
           !WRITE(*,*) "Calling pipj_osc_basis..."
           CALL pipj_osc_basis(ang_mom, gmatrix_configs, reduced_p1, twobody_p2 )

           !WRITE(*,*) "Calling vlowk_print..."
           CALL vlowk_print(istrange, isospin_z,p_parity,ang_mom,gna,twobody_com,twobody_r2,twobody_p2,gmatrix_configs)
           !      free space
           DEALLOCATE(gmatrix_configs%config_ab)
           DEALLOCATE(twobody_com)
           DEALLOCATE(twobody_r2,twobody_p2)
           DEALLOCATE(lab_to_relcoeff, lab_to_relconf )
           DEALLOCATE (lab_to_relnumber); DEALLOCATE(gna)

           !WRITE(*,*) "Done deallocating..."
        ENDDO
     ENDDO
  ENDDO
  ENDDO

  DEALLOCATE( reduced_p1 )


END SUBROUTINE final_vlowk_labsystem

SUBROUTINE vlowk_free(gfree,lab_to_relcoeff, lab_to_relconf,&
     lab_to_relnumber,max_coeff,gmatrix_configs)
  USE configurations
  USE constants
  USE single_particle_orbits
  USE partial_waves
  USE relcm_gmatrix
  IMPLICIT NONE
  TYPE (configuration_descriptor) , INTENT(IN)  :: gmatrix_configs
  INTEGER :: k2,k1,i,j,nlc1,nlc2, conf1, conf2, lcm1, lcm2, chan1, chan2, as, ls, ns,&
       max_coeff, ncm1, ncm2, lr1, lr2, nr1, nr2, bra, ket, state, a,&
       b,c,d,ida, idb, idc, idd, subab, subcd, isub, nconf, id1, id2
  COMPLEX*16, DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs), INTENT(INOUT) :: gfree
  REAL(KIND=8) :: dcoem
  REAL(KIND=8), DIMENSION(gmatrix_configs%number_confs,max_coeff), &
       INTENT(IN)  :: lab_to_relcoeff
  INTEGER, DIMENSION(gmatrix_configs%number_confs,max_coeff), &
       INTENT(IN)  :: lab_to_relconf
  INTEGER, DIMENSION(gmatrix_configs%number_confs), &
       INTENT(IN)  :: lab_to_relnumber

  gfree = 0.D0
  DO k2=1,gmatrix_configs%number_confs
     a=gmatrix_configs%config_ab(2*k2)
     b=gmatrix_configs%config_ab(2*k2-1)
     ida=all_orbit%ptype(a)
     idb=all_orbit%ptype(b)
     DO k1=1,gmatrix_configs%number_confs
        c=gmatrix_configs%config_ab(2*k1-1)
        d=gmatrix_configs%config_ab(2*k1)
        idc=all_orbit%ptype(c)
        idd=all_orbit%ptype(d)
        DO j=1,lab_to_relnumber(k2)
           nlc2=lab_to_relconf(k2,j)
           conf2=nlc2/1000;  chan2=nlc2-conf2*1000
           lcm2=relcm_sp_data%lrel(relcm_conf%relcm_ab(chan2,conf2+conf2))
           ncm2=relcm_sp_data%nrel(relcm_conf%relcm_ab(chan2,conf2+conf2))
           DO i=1,lab_to_relnumber(k1)
              nlc1=lab_to_relconf(k1,i)
              conf1=nlc1/1000;  chan1=nlc1-conf1*1000
              IF ( chan1 /= chan2 ) CYCLE
              lcm1=relcm_sp_data%lrel(relcm_conf%relcm_ab(chan1,conf1+conf1))
              IF ( lcm1 /= lcm2 ) CYCLE   ! cm orbital mom must be equal
              ncm1=relcm_sp_data%nrel(relcm_conf%relcm_ab(chan1,conf1+conf1))
              IF ( ncm1 /= ncm2 ) CYCLE   ! cm N must be equal
              lr1=relcm_sp_data%lrel(relcm_conf%relcm_ab(chan1,conf1+conf1-1))
              nr1=relcm_sp_data%nrel(relcm_conf%relcm_ab(chan1,conf1+conf1-1))
              lr2=relcm_sp_data%lrel(relcm_conf%relcm_ab(chan2,conf2+conf2-1))
              nr2=relcm_sp_data%nrel(relcm_conf%relcm_ab(chan2,conf2+conf2-1))

              ! Find no elements pr subchannel
              nconf = rel_conf%nconfs_rel(chan1)
              ! Find subchannel
              subab = 0; subcd = 0
              DO isub = 1,no_sub(chan1)
                id1 = sub_id(chan1,isub)/100
                id2 = sub_id(chan1,isub)-100*id1
                IF ( (id1 == ida) .AND. (id2 == idb)) subab = isub
                IF ( (id2 == ida) .AND. (id1 == idb)) subab = isub
                IF ( (id1 == idc) .AND. (id2 == idd)) subcd = isub
                IF ( (id2 == idc) .AND. (id1 == idd)) subcd = isub
              ENDDO

              IF (subab == 0) THEN
                WRITE(*,*) "Something is wrong. Couldn't find subchannel.."
                WRITE(*,*) ida,idb,chan1
                DO isub = 1,no_sub(chan1)
                    WRITE(*,*) "nsub: ", no_sub(chan1)
                    WRITE(*,*) isub, sub_id(chan1,isub)
                ENDDO
              ENDIF

              bra = 0; ket = 0
              DO state =1, rel_conf%nconfs_rel(chan1)
                 as=rel_conf%rel_ab(chan1,state)
                 ns=relcm_sp_data%nrel(as)
                 ls=relcm_sp_data%lrel(as)
                 IF ( (ns == nr1) .AND. (ls == lr1)) ket = state
                 IF ( (ns == nr2) .AND. (ls == lr2)) bra = state
              ENDDO
              bra = bra + nconf*(subab-1)
              ket = ket + nconf*(subcd-1)
              dcoem=lab_to_relcoeff(k1,i)*lab_to_relcoeff(k2,j)

              gfree(k1,k2)= gfree(k1,k2)+dcoem*v_com_rel(chan1,bra,ket) 
              !WRITE(*,*) v_com_rel(chan1,bra,ket), chan1, bra, ket

           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE vlowk_free

SUBROUTINE vlowk_print(str,it,ip,ij,gna,twobody_com,twobody_r2, twobody_p2,gmatrix_configs)
  USE constants
  USE configurations
  USE single_particle_orbits
  USE relcm_gmatrix
  USE wave_functions
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN)  :: gmatrix_configs
  INTEGER :: i,j, ijd, it, ip, ij, ia, ib, ic,id, str
  COMPLEX*16, DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs), INTENT(IN)  :: gna
  REAL(KIND=8), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs), INTENT(IN) :: twobody_com,twobody_r2, twobody_p2
  ijd=ij+ij
  DO i=1,gmatrix_configs%number_confs
     ia= gmatrix_configs%config_ab(i*2-1)
     ib= gmatrix_configs%config_ab(i*2)
     DO j=i,gmatrix_configs%number_confs
        ic= gmatrix_configs%config_ab(j+j-1)
        id= gmatrix_configs%config_ab(j+j)
        WRITE(7,'(8I4,5X,4(5X,E12.6))') str,it, ip, ijd, ia, ib, ic, id, REAL(gna(j,i)), &
             twobody_com(j,i), twobody_r2(j,i), twobody_p2(j,i)       
     ENDDO
  ENDDO

END SUBROUTINE vlowk_print

SUBROUTINE coulomb_integral
  USE constants
  USE wave_functions
  USE relcm_gmatrix
  IMPLICIT NONE
  INTEGER :: n1, n2, lr, i, nrel_max, sub, id1, id2, q
  REAL(DP) :: oscl_r, int_sum, xr, xp, z, factor1, factor2, m
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: rr, wrr
  REAL(DP) :: cx(0:200)

  ! Change this to calculate 6 different integrals. One for each subchannel
  nrel_max = 400
  ALLOCATE ( rr(nrel_max ), wrr(nrel_max ))
  CALL gauss_legendre(0.d0, 20.d0, rr, wrr, nrel_max )
  DO sub = 1,6
    id1 = coulomb_channel(sub)/100
    id2 = coulomb_channel(sub) - id1*100
    q = sp_charge(id1)*sp_charge(id2) ! 1 or -1
    m = sp_mass(id1) * sp_mass(id2)/(sp_mass(id1) + sp_mass(id2)) ! Rel mass
    oscl_r = hbarc/SQRT(m*hbar_omega)
  DO lr = 0, lmax, 1
     DO n1 = 0, nmax, 1
        factor1 = 0.5D0*((n1+lr+2)*LOG(2.D0)+fac(n1)-dfac(2*n1+2*lr+1)-0.5D0*LOG(pi))
        factor1 = EXP(factor1)
        DO n2 = 0, nmax,1
           IF ( n2 == n1) THEN
              factor2 = factor1
           ELSE
              factor2 = 0.5D0*((n2+lr+2)*LOG(2.D0)+fac(n2)-dfac(2*n2+2*lr+1)-0.5D0*LOG(pi))
              factor2 = EXP(factor2)
           ENDIF
           int_sum = 0.D0
           DO i=1,nrel_max 
              z= rr(i)/oscl_r
              CALL laguerre_general( n1, lr+0.5D0, z*z, cx )
              xp = cx(n1)*EXP(-z*z*0.5)*(z**lr)
              IF ( n1 == n2) THEN 
                 xr = xp
              ELSE
                 CALL laguerre_general( n2, lr+0.5D0, z*z, cx )
                 xr = cx(n2)*EXP(-z*z*0.5)*(z**lr)
              ENDIF
              int_sum=int_sum+wrr(i)*rr(i)*xp*xr
           ENDDO
           !  Coulomb energy in MeV
           ! Should be changed to q_1q_2/4/pi/epislon_0 for this subchannel
           coulomb_relcom(n2, n1, lr, sub) = int_sum*factor1*factor2*q*1.439965183D0/(oscl_r**3)
        ENDDO
     ENDDO
  ENDDO
  ENDDO
  DEALLOCATE ( rr, wrr)

END SUBROUTINE coulomb_integral

