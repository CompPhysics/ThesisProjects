!             Program block bhf-matrix.f
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
!             E-MAIL:   morten.hjorth-jensen@fys.uio.no
!             LANGUAGE: F90/F95  
!
!                    Begin bhf-matrix code
!
SUBROUTINE  setup_hfmatrix
  USE constants
  USE single_particle_orbits
  USE wave_functions
  USE configurations
  IMPLICIT NONE
  REAL(DP), ALLOCATABLE, DIMENSION(:,:)  :: coeffs  

  !     reserve space in memory for various arrays
  ALLOCATE ( ra (n_rel), wra (n_rel));  ALLOCATE ( rgkk (n_rel), wgkk (n_rel))
  ALLOCATE ( hol (n_rel, 0:lmax, 0:nmax, 6) )
  !     set up mesh points in lab frame
  CALL rel_mesh                   
  !     setup ho wave functions
  CALL ho_wfunction                
  !    Expansion coefficients for sp wave functions. 
  ALLOCATE(coeffs(all_orbit%total_orbits,all_orbit%total_orbits))
  ALLOCATE ( bhf_hol (n_rel,all_orbit%total_orbits))
  ALLOCATE ( wave_function (n_rel,all_orbit%total_orbits))
  coeffs = 0.0_dp; wave_function = 0.0_dp; bhf_hol = 0.0_dp
  !  perform the HF calculation
  IF (hf_iterations /= 0) THEN
    CALL brueckner_hartree_fock(coeffs,all_orbit%total_orbits)
    !  update the g-matrix
    CALL setupg_bhf(coeffs,all_orbit%total_orbits)
  ENDIF
  SELECT CASE (type_of_interaction)
  CASE ('open-diagrams')
     CALL onebody_contribution
  END SELECT
  DEALLOCATE(coeffs) ;    DEALLOCATE ( bhf_hol,wave_function )
  DEALLOCATE ( rgkk, wgkk, ra, wra) ; DEALLOCATE ( hol)

END SUBROUTINE setup_hfmatrix
!
!           Set up the BHF G-mtx in the lab-frame 
!
SUBROUTINE setupg_bhf(coeffs,ncoeffs)
  USE single_particle_orbits
  USE configurations
  USE constants
  USE hoosc_gmatrix
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: gmatrix_configs
  REAL(DP), ALLOCATABLE :: bhf_coeff(:,:)
  INTEGER, INTENT(IN) :: ncoeffs
  REAL(DP), DIMENSION(ncoeffs,ncoeffs), INTENT(IN)  :: coeffs  
  INTEGER ::  p_parity, ang_mom, isospin_z, n_confs, ie, istrange
  REAL(DP), ALLOCATABLE :: gna(:,:,:),  temp(:,:,:)

  !     loop over isospin projection
  DO istrange = strange_min, strange_max
    itzmin = istrange -2; itzmax = -1*itzmin
    DO isospin_z=itzmin,itzmax,2 
      !     loop over parity values, here positive parity is 0, negative 1
      DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=j_lab_min,j_lab_max
           !     find all possible configurations 
           CALL number_gmatrix_confs &
                (ang_mom,p_parity,isospin_z,istrange, gmatrix_configs)
           IF (gmatrix_configs%number_confs <= 0 ) CYCLE
           n_confs=gmatrix_configs%number_confs
           !WRITE(*,*) istrange, isospin_z, p_parity, ang_mom, n_confs
           ALLOCATE(gmatrix_configs%config_ab(n_confs+n_confs) )
           CALL setup_gmatrix_configurations &
                (ang_mom,p_parity,isospin_z,istrange, gmatrix_configs)
           ALLOCATE(temp(n_confs, n_confs,n_startenergy_g+3))
           ALLOCATE(gna(n_confs, n_confs,n_startenergy_g+3))
           ALLOCATE(bhf_coeff(n_confs, n_confs))
           gna=0.0_dp; bhf_coeff = 0.0_dp; temp = 0.0_dp
           CALL fetch_matrix(ang_mom,gna,gmatrix_configs)
           CALL bhf_coefficients(ang_mom,n_confs,ncoeffs,gmatrix_configs,bhf_coeff,coeffs)
           DO ie=1,n_startenergy_g+3
              temp(:,:,ie) = MATMUL(gna(:,:,ie),TRANSPOSE(bhf_coeff(:,:)))
              gna(:,:,ie) = MATMUL(bhf_coeff(:,:),temp(:,:,ie))
           ENDDO
           !     update table of matrix elements 
           CALL update(istrange,isospin_z,p_parity,ang_mom,gna,gmatrix_configs)
           !     free space in heap
           DEALLOCATE(gmatrix_configs%config_ab)
           DEALLOCATE(bhf_coeff)
           DEALLOCATE(gna, temp)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE setupg_bhf
!
!                    Update out G-mtx
!
SUBROUTINE update(istr, it,ip,ij,gna,gmatrix_configs)
  USE configurations
  USE constants
  USE single_particle_orbits
  USE relcm_gmatrix
  USE wave_functions
  USE stored_bare_interaction
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN)  :: gmatrix_configs
  INTEGER :: i,j, ijd, it, ip, ij, ia, ib, ic,id,ie, istr
  REAL(DP), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs,n_startenergy_g+3), INTENT(IN)  :: gna

  
  !WRITE(*,*) "Update called..."
  ijd=ij+ij
  DO i=1,gmatrix_configs%number_confs
     ia=gmatrix_configs%config_ab(i*2-1)
     ib=gmatrix_configs%config_ab(i*2)
     DO j=i,gmatrix_configs%number_confs
        ic=gmatrix_configs%config_ab(j+j-1)
        id=gmatrix_configs%config_ab(j+j)
        CALL replace_g(ia,ib, ic, id, istr, it, ip, ij ,gna(j,i,:))
        WRITE(9,'(8I4,5X,10(5X,E12.6))') istr, it, ip, ijd, ia, ib, ic, id, (gna(j,i,ie),ie=1,n_startenergy_g+3)
     ENDDO
  ENDDO

END SUBROUTINE update
!
!                    Get just the G-mtx to be used in various HF iterations
!
SUBROUTINE fetch_matrix(ij,gna,gmatrix_configs)
  USE configurations
  USE constants
  USE single_particle_orbits
  USE relcm_gmatrix
  USE wave_functions
  USE stored_bare_interaction
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN)  :: gmatrix_configs
  INTEGER :: i,j, ij, ia, ib, ic,id
  REAL(DP) :: norm, dij
  REAL(DP), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs,n_startenergy_g+3), INTENT(INOUT)  :: gna
  REAL(DP), DIMENSION(n_startenergy_g+3) :: ans

  
  !WRITE(*,*) "Fetch_matrix called..."
  DO i=1,gmatrix_configs%number_confs
     ia=gmatrix_configs%config_ab(i*2-1)
     ib=gmatrix_configs%config_ab(i*2)
     DO j=i,gmatrix_configs%number_confs
        ic=gmatrix_configs%config_ab(j+j-1)
        id=gmatrix_configs%config_ab(j+j)
        CALL pphhmtx(ia,ib,ic,id,ij,ans)
        norm=1.0_dp/dij(ia,ib)/dij(ic,id)
        gna(j,i,:)=ans(:)*norm
        gna(i,j,:) = ans(:)*norm
     ENDDO
  ENDDO

END SUBROUTINE fetch_matrix

SUBROUTINE bhf_coefficients(ang_mom,n_confs,ncoeffs,gmatrix_configs,bhf_coeff,coeffs)
  USE single_particle_orbits
  USE configurations
  USE constants
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN)  :: gmatrix_configs
  INTEGER, INTENT(IN) :: ncoeffs, n_confs, ang_mom
  REAL(DP), DIMENSION(ncoeffs,ncoeffs), INTENT(IN)  :: coeffs  
  REAL(DP), DIMENSION(n_confs,n_confs), INTENT(INOUT)  :: bhf_coeff
  REAL(DP) :: dij, fnorm, fnorm_ab, fnorm_pq
  REAL(DP) :: e_dir, e_exc, e
  INTEGER :: ket , bra, a,b, la, ja, lb, jb, p,q, lq, lp, jp, jq, &
    ida, idb, idq, idp

  
  !WRITE(*,*) "bhf_coefficients called..."
  fnorm = 0.d0; fnorm_ab = 1.0d0; fnorm_pq = 1.0d0
  DO bra = 1, n_confs
     a = gmatrix_configs%config_ab(bra+bra-1)
     b = gmatrix_configs%config_ab(bra+bra)
     la=all_orbit%ll(a); jb=all_orbit%jj(b)
     ja=all_orbit%jj(a); lb=all_orbit%ll(b)
     ida = all_orbit%ptype(a); idb = all_orbit%ptype(b)
     DO ket = 1, n_confs
        p = gmatrix_configs%config_ab(ket+ket-1)
        q = gmatrix_configs%config_ab(ket+ket)
        lp = all_orbit%ll(p); jp=all_orbit%jj(p)
        lq = all_orbit%ll(q); jq=all_orbit%jj(q)
        idp = all_orbit%ptype(p); idq = all_orbit%ptype(q)
        IF (ida == idb) fnorm_ab = 1.d0/ dij(a,b)
        IF (idp == idq) fnorm_pq = 1.d0/ dij(q,p)
        IF (( ida /= idb ) .AND. (idp /= idq) )THEN
           fnorm = 1.d0
        ELSE
           fnorm = fnorm_ab*fnorm_pq
        ENDIF
        e_dir=0.0_dp; e_exc=0.D0; e = 0.0_dp
        ! direct term
        IF ( (la == lp).AND.( lb == lq ).AND.( ja == jp ).AND.( jb == jq )) THEN
            IF (ida == idp) e_dir = coeffs(p,a)
            IF (idb == idq) e_dir = e_dir*coeffs(q,b)
        ENDIF
        ! exchange term
        IF ( (la == lq).AND.( lb == lp ).AND.( ja == jq ).AND.( jb == jp ) ) THEN
            IF (ida == idq) e_exc = coeffs(q,a)
            IF (idb == idp) e_exc = e_exc*coeffs(p,b)
        ENDIF
        e= (e_dir-e_exc*((-1.0_dp)**((2*ang_mom-jp-jq)/2)))*fnorm
        bhf_coeff(bra,ket) = e
     ENDDO
  ENDDO

END SUBROUTINE bhf_coefficients
!
!                 Set up h.o. wf for rel cm system and lab frame
!                 It computes also a complex
!
SUBROUTINE ho_wfunction
  USE constants
  USE wave_functions
  IMPLICIT NONE
  INTEGER :: n, l, i, p
  REAL(DP) :: ph, sum_rel
  REAL(DP) :: osclp
  REAL(DP)  :: cx(0:200), factor, z_lab,xp

  DO p = 1,6
  osclp = hbarc/SQRT(sp_mass(p)*hbar_omega)
  DO n=0,nmax
     ph=(-1.D0)**n
     DO l=0,lmax
        factor = 0.5D0*((n+l+2)*LOG(2.D0)+fac(n)-dfac(2*n+2*l+1)-0.5D0*LOG(pi))
        factor = EXP(factor)
        sum_rel=0.0_dp
        DO i=1,n_rel
           !  real ho wave function
           z_lab= ra(i)*ra(i)*osclp*osclp; cx = 0.0_dp
           CALL laguerre_general( n, l+0.5D0, z_lab, cx )
           xp = EXP(-z_lab*0.5D0)*((ra(i)*osclp)**l)*cx(n)
           hol(i,l,n, p) = xp*ph*factor*(osclp**(1.5D0)) ! lab wf
           sum_rel=sum_rel+ wra(i)*(hol(i,l,n, p)*ra(i))**2
        ENDDO
        WRITE(6,'(22HNorm cm ho wf n,l, id:,3I3,2X,F12.7)') n, l, p, sum_rel 
     ENDDO
  ENDDO
  ENDDO

END SUBROUTINE ho_wfunction
!
!     Brueckner-Hartree-Hock self-consistency
!
SUBROUTINE brueckner_hartree_fock(coeffs, ncoeffs)
  USE single_particle_orbits
  USE configurations
  USE constants
  USE wave_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ncoeffs
  REAL(DP), DIMENSION(ncoeffs,ncoeffs), INTENT(INOUT)  :: coeffs  
  INTEGER :: a, c, ket , bra, nrot, hf_iter, max_hfiter
  INTEGER :: la, lc, na, nc, k, i, id
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: hf_vect
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: hf_mtx
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: hf_eigen
  REAL(DP), DIMENSION(n_rel) ::  sum_wf
  REAL(DP) :: e_kin, hf, sigma, factor
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: int_factor, kin_energy

  ALLOCATE (int_factor(n_rel), kin_energy(n_rel) )
  int_factor(:)=wra(:)*ra(:)*ra(:)
  kin_energy(:)=hbarc*hbarc*ra(:)*ra(:)/2

  !     initialize bhf harmonic oscillator sp wave function
  bhf_hol=0.0_dp; max_hfiter = 100
  ALLOCATE (hf_mtx(all_orbit%total_orbits,all_orbit%total_orbits))
  ALLOCATE (hf_vect(all_orbit%total_orbits,all_orbit%total_orbits))
  ALLOCATE (hf_eigen(all_orbit%total_orbits))
  hf_vect = 0.0_dp
  !  for first iteration coeffs has only diagonal values equal 1
  DO bra = 1, all_orbit%total_orbits
     coeffs(bra,bra) = 1.0_dp
  ENDDO
  hf_iter = 0; sigma = 1.0_dp
  DO WHILE ((hf_iter <= max_hfiter).AND.(ABS(sigma) > 1E-8))
     !     set up bhf matrix to be diagonalized
     hf_mtx = 0.0_dp; hf_vect = 0.0_dp
     DO bra = 1, all_orbit%total_orbits
        a = bra
        IF (all_orbit%included(a) == 'exclude') CYCLE
        DO ket= bra, all_orbit%total_orbits
           c = ket
           IF (all_orbit%included(c) == 'exclude') CYCLE
           IF(all_orbit%ptype(a) /= all_orbit%ptype(c) ) CYCLE
           IF(all_orbit%ll(a) /= all_orbit%ll(c)) CYCLE
           IF ( all_orbit%jj(a) /= all_orbit%jj(c)) CYCLE
           id = all_orbit%ptype(a)
           hf = 0.0_dp               
           CALL  diagram_HF(a, c, coeffs, hf, ncoeffs)
           !  compute the kinetic energy or unperturbed one-body H
           e_kin = 0.0_dp
           SELECT CASE (type_of_renormv)
           CASE ('vlowk')
              DO k=1,n_rel
                 e_kin=e_kin+hol(k,all_orbit%ll(a),all_orbit%nn(a), id)* &
                      hol(k,all_orbit%ll(c),all_orbit%nn(c), id)* &
                      int_factor(k)*kin_energy(k)
              ENDDO
           END SELECT
           factor = (1.0/sp_mass(id) - 1.0/total_mass)
           e_kin = e_kin*factor
           hf_mtx(bra,ket) = e_kin+hf
           hf_mtx(ket,bra) = hf_mtx(bra,ket)
           WRITE(6,'(2I4,3X,3F16.8)') a, c, hf, e_kin, e_kin+hf
        ENDDO
     ENDDO
     !     obtain the BHF coefficients 
     CALL matrix_diag(hf_mtx,all_orbit%total_orbits, all_orbit%total_orbits,&
          hf_eigen, hf_vect,nrot)
     !     set up the new bhf harmonic  oscillator wave function
     !     Note memory stride
     sigma = 0.0_dp
     DO bra = 1, all_orbit%total_orbits
        a = bra
        la = all_orbit%ll(a)
        na = all_orbit%nn(a)
        sum_wf = 0.0_dp
        DO ket= 1, all_orbit%total_orbits
           c = ket
           lc = all_orbit%ll(c)
           nc = all_orbit%nn(c)
           id = all_orbit%ptype(c)
           sum_wf(:)=sum_wf(:)+hf_vect(ket,bra)*hol(:,lc,nc, id)
           coeffs(a,c) = hf_vect(bra,ket)
        ENDDO
        bhf_hol(:,a)=sum_wf(:)
        ! set up of new single-particle energy 
        !WRITE(6,*) a, "\n", bhf_hol(:,a)
        sigma = sigma +ABS(all_orbit%e(a) - hf_eigen(bra))
        all_orbit%e(a) = hf_eigen(bra)                    
        !WRITE(6,*) bra, hf_eigen(bra)
     ENDDO
     sigma = sigma/all_orbit%total_orbits
     WRITE(6,*) 'iteration nr and sigma', hf_iter, sigma
     DO i=1, all_orbit%total_orbits
        IF (all_orbit%orbit_status(i) /= 'hole') CYCLE         
        WRITE(6,'(7HNumber:,6(I3,2X),2X,E20.10)') i, all_orbit%nn(i), all_orbit%ll(i), &
             all_orbit%jj(i), &
             all_orbit%nshell(i), all_orbit%itzp(i), all_orbit%e(i)
     ENDDO
     hf_iter = hf_iter+1
  ENDDO
  DO i=1, all_orbit%total_orbits
     all_orbit%evalence(i) = all_orbit%e(i)
     all_orbit%e(i) = all_orbit%e_original(i)
     WRITE(8,'(7HNumber:,8(I4,2X),2X,E12.6,2X,E12.6,2X,A10,2X,A10,2X,A10)') i, &
        all_orbit%ptype(i), &
        all_orbit%strange(i), all_orbit%nn(i), all_orbit%ll(i), all_orbit%jj(i), &
        all_orbit%itzp(i), all_orbit%nshell(i), all_orbit%e(i), all_orbit%evalence(i), &
          all_orbit%orbit_status(i), all_orbit%model_space(i),&
          all_orbit%included(i)
  ENDDO
  DEALLOCATE ( hf_mtx)
  DEALLOCATE ( hf_eigen)
  DEALLOCATE ( hf_vect)
  !     new harmonic oscillator wave function
  wave_function=bhf_hol
  DEALLOCATE (int_factor, kin_energy )

END SUBROUTINE brueckner_hartree_fock
!
!    The Hartree-Fock diagram
!
SUBROUTINE diagram_HF(a,c,coeffs, onebody_diagram_HF, ncoeffs)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c, ncoeffs
  INTEGER :: j_min, j_max, jph, h, h1, h2
  REAL(DP), DIMENSION(ncoeffs,ncoeffs), INTENT(IN)  :: coeffs  
  REAL(DP) :: ang_mom_factor
  REAL(DP), INTENT(INOUT) :: onebody_diagram_HF
  REAL(DP), DIMENSION(n_startenergy_g+3) :: ans

  onebody_diagram_HF=0.0_dp
  DO h=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE         
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     DO h1=1, all_orbit%total_orbits
        IF  (all_orbit%jj(h) /= all_orbit%jj(h1)) CYCLE
        IF  (all_orbit%ll(h) /= all_orbit%ll(h1)) CYCLE
        IF  (all_orbit%ptype(h) /= all_orbit%ptype(h1)) CYCLE
        DO h2=1, all_orbit%total_orbits
           IF  (all_orbit%jj(h) /= all_orbit%jj(h2)) CYCLE
           IF  (all_orbit%ll(h) /= all_orbit%ll(h2)) CYCLE
           IF  (all_orbit%ptype(h) /= all_orbit%ptype(h2)) CYCLE
           SELECT CASE (type_of_renormv)
           CASE ('vlowk')
              DO jph=j_min,j_max
                 ang_mom_factor=(2.*jph+1.)/(all_orbit%jj(a)+1.)
                 CALL pphhmtx(a,h1,c,h2,jph,ans); IF ( ans(1) == 0.0_dp ) CYCLE
                 onebody_diagram_HF=onebody_diagram_HF+ans(1)*ang_mom_factor*coeffs(h1,h)*coeffs(h2,h)
                 IF (coeffs(h1,h)*coeffs(h2,h) > 1E-6) THEN
                 !write(6,'(5I4,2X,F4.2X,2x4F16.8)') a,h1,c,h2,jph, ang_mom_factor,&
                 !   coeffs(h1,h), coeffs(h2,h), ans(1), onebody_diagram_HF
                 ENDIF
              ENDDO
           END SELECT
        ENDDO
     ENDDO
  ENDDO

END  SUBROUTINE diagram_HF

SUBROUTINE diagram_HF_2a(a,c,coeffs,onebody_diagram_HF, ncoeffs)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c, ncoeffs
  INTEGER :: j_min, j_max, jph, h, h1, h2, p1, p2, p11, p12, p21, p22
  REAL(DP), DIMENSION(ncoeffs,ncoeffs), INTENT(IN)  :: coeffs  
  REAL(DP) :: ang_mom_factor, e_factor
  REAL(DP), INTENT(INOUT) :: onebody_diagram_HF
  REAL(DP), DIMENSION(n_startenergy_g+3) :: ans1, ans2

  onebody_diagram_HF=0.0_dp
  DO h=1, all_orbit%total_orbits
   IF (all_orbit%orbit_status(h) /= 'hole') CYCLE         
   j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
   j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
   DO p1=1, all_orbit%total_orbits
    IF(all_orbit%orbit_status(p1) /= 'particle') CYCLE
    DO p2=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p2) /= 'particle') CYCLE
     DO h1=1, all_orbit%total_orbits
      IF  (all_orbit%jj(h) /= all_orbit%jj(h1)) CYCLE
      IF  (all_orbit%ll(h) /= all_orbit%ll(h1)) CYCLE
      IF  (all_orbit%ptype(h) /= all_orbit%ptype(h1)) CYCLE
      DO h2=1, all_orbit%total_orbits
       IF  (all_orbit%jj(h) /= all_orbit%jj(h2)) CYCLE
       IF  (all_orbit%ll(h) /= all_orbit%ll(h2)) CYCLE
       IF  (all_orbit%ptype(h) /= all_orbit%ptype(h2)) CYCLE
       DO p11=1, all_orbit%total_orbits
        IF  (all_orbit%jj(p1) /= all_orbit%jj(p11)) CYCLE
        IF  (all_orbit%ll(p1) /= all_orbit%ll(p11)) CYCLE
        IF  (all_orbit%ptype(p1) /= all_orbit%ptype(p11)) CYCLE
        DO p12=1, all_orbit%total_orbits
         IF  (all_orbit%jj(p1) /= all_orbit%jj(p12)) CYCLE
         IF  (all_orbit%ll(p1) /= all_orbit%ll(p12)) CYCLE
         IF  (all_orbit%ptype(p1) /= all_orbit%ptype(p12)) CYCLE
         DO p21=1, all_orbit%total_orbits
          IF  (all_orbit%jj(p2) /= all_orbit%jj(p21)) CYCLE
          IF  (all_orbit%ll(p2) /= all_orbit%ll(p21)) CYCLE
          IF  (all_orbit%ptype(p2) /= all_orbit%ptype(p21)) CYCLE
          DO p22=1, all_orbit%total_orbits
           IF  (all_orbit%jj(p2) /= all_orbit%jj(p22)) CYCLE
           IF  (all_orbit%ll(p2) /= all_orbit%ll(p22)) CYCLE
           IF  (all_orbit%ptype(p2) /= all_orbit%ptype(p22)) CYCLE
              e_factor = all_orbit%e(c) + all_orbit%e(h) - &
              all_orbit%e(p1) - all_orbit%e(p2)
              IF (ABS(e_factor) < 1E-05) THEN
                WRITE(*,*) "Denominator too small in diagram_HF_2a."
                CYCLE
              ELSE
                e_factor = 1.0/e_factor
              ENDIF
              SELECT CASE (type_of_renormv)
              CASE ('vlowk')
                DO jph=j_min,j_max
                  ang_mom_factor=(2.*jph+1.)/(4*(all_orbit%jj(a)+1.))
                  CALL pphhmtx(a,h1,p11,p21,jph,ans1); IF ( ans1(1) == 0.0_dp ) CYCLE
                  CALL pphhmtx(p21,p22,c,h2,jph,ans2); IF ( ans1(1) == 0.0_dp ) CYCLE
                  onebody_diagram_HF=onebody_diagram_HF+ans1(1)*ans2(1)*ang_mom_factor*&
                    coeffs(h1,h)*coeffs(h2,h)*coeffs(p11,p1)*coeffs(p12,p1)*&
                    coeffs(p21,p2)*coeffs(p22,p2)*e_factor
                ENDDO
              END SELECT
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO

END  SUBROUTINE diagram_HF_2a

SUBROUTINE diagram_HF_2b(a,c,coeffs,onebody_diagram_HF, ncoeffs)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c,  ncoeffs
  INTEGER :: j_min, j_max, jph, p, p1, p2, h1, h2, h11, h12, h21, h22
  REAL(DP), DIMENSION(ncoeffs,ncoeffs), INTENT(IN)  :: coeffs  
  REAL(DP) :: ang_mom_factor, e_factor
  REAL(DP), INTENT(INOUT) :: onebody_diagram_HF
  REAL(DP), DIMENSION(n_startenergy_g+3) :: ans1, ans2

  onebody_diagram_HF=0.0_dp
  DO p=1, all_orbit%total_orbits
   IF (all_orbit%orbit_status(p) /= 'particle') CYCLE         
   j_min=ABS((all_orbit%jj(a)-all_orbit%jj(p))/2)
   j_max=(all_orbit%jj(a)+all_orbit%jj(p))/2
   DO h1=1, all_orbit%total_orbits
    IF(all_orbit%orbit_status(h1) /= 'hole') CYCLE
    DO h2=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h2) /= 'hole') CYCLE
     DO p1=1, all_orbit%total_orbits
      IF  (all_orbit%jj(p) /= all_orbit%jj(p1)) CYCLE
      IF  (all_orbit%ll(p) /= all_orbit%ll(p1)) CYCLE
      IF  (all_orbit%ptype(p) /= all_orbit%ptype(p1)) CYCLE
      DO p2=1, all_orbit%total_orbits
       IF  (all_orbit%jj(p) /= all_orbit%jj(p2)) CYCLE
       IF  (all_orbit%ll(p) /= all_orbit%ll(p2)) CYCLE
       IF  (all_orbit%ptype(p) /= all_orbit%ptype(p2)) CYCLE
       DO h11=1, all_orbit%total_orbits
        IF  (all_orbit%jj(h1) /= all_orbit%jj(h11)) CYCLE
        IF  (all_orbit%ll(h1) /= all_orbit%ll(h11)) CYCLE
        IF  (all_orbit%ptype(h1) /= all_orbit%ptype(h11)) CYCLE
        DO h12=1, all_orbit%total_orbits
         IF  (all_orbit%jj(h1) /= all_orbit%jj(h12)) CYCLE
         IF  (all_orbit%ll(h1) /= all_orbit%ll(h12)) CYCLE
         IF  (all_orbit%ptype(h1) /= all_orbit%ptype(h12)) CYCLE
         DO h21=1, all_orbit%total_orbits
          IF  (all_orbit%jj(h2) /= all_orbit%jj(h21)) CYCLE
          IF  (all_orbit%ll(h2) /= all_orbit%ll(h21)) CYCLE
          IF  (all_orbit%ptype(h2) /= all_orbit%ptype(h21)) CYCLE
          DO h22=1, all_orbit%total_orbits
           IF  (all_orbit%jj(h2) /= all_orbit%jj(h22)) CYCLE
           IF  (all_orbit%ll(h2) /= all_orbit%ll(h22)) CYCLE
           IF  (all_orbit%ptype(h2) /= all_orbit%ptype(h22)) CYCLE
              e_factor = all_orbit%e(h1) + all_orbit%e(h2) - &
              all_orbit%e(p) - all_orbit%e(a)
              IF (ABS(e_factor) < 1E-05) THEN
                WRITE(*,*) "Denominator too small in diagram_HF_2b."
                CYCLE
              ELSE
                e_factor = 1.0/e_factor
              ENDIF
              SELECT CASE (type_of_renormv)
              CASE ('vlowk')
                DO jph=j_min,j_max
                  ang_mom_factor=(2.*jph+1.)/(4*(all_orbit%jj(a)+1.))
                  CALL pphhmtx(a,p1,h11,h21,jph,ans1); IF ( ans1(1) == 0.0_dp ) CYCLE
                  CALL pphhmtx(h21,h22,c,p2,jph,ans2); IF ( ans1(1) == 0.0_dp ) CYCLE
                  onebody_diagram_HF=onebody_diagram_HF+ans1(1)*ans2(1)*ang_mom_factor*&
                    coeffs(p1,p)*coeffs(p2,p)*coeffs(h11,h1)*coeffs(h12,h1)*&
                    coeffs(h21,h2)*coeffs(h22,h2)*e_factor
                ENDDO
              END SELECT
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  onebody_diagram_HF = -onebody_diagram_HF

END  SUBROUTINE diagram_HF_2b


!
!     This function returns the matrix for V
!
SUBROUTINE pphhmtx(ja,jb,jc,jd,jt,gmtpn)
  USE stored_bare_interaction
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  LOGICAL TRIAG
  INTEGER, INTENT(IN) :: ja,jb,jc,jd,jt
  DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n_startenergy_g+3) :: gmtpn

  !WRITE(*,*) "pphhmatrix called..."
  gmtpn=0.
  IF(2*all_orbit%nn(ja)+all_orbit%ll(ja) +2*all_orbit%nn(jb)+all_orbit%ll(jb) > nlmax ) RETURN
  IF(2*all_orbit%nn(jc)+all_orbit%ll(jc) +2*all_orbit%nn(jd)+all_orbit%ll(jd) > nlmax ) RETURN
  IF((all_orbit%strange(ja) + all_orbit%strange(jb)) /= &
    (all_orbit%strange(jc)+all_orbit%strange(jd))) RETURN
  IF((all_orbit%itzp(ja)+all_orbit%itzp(jb)) /= &
       (all_orbit%itzp(jc)+all_orbit%itzp(jd))) RETURN
  IF((-1)**(all_orbit%ll(ja)+all_orbit%ll(jb)) /=  &
       (-1)**(all_orbit%ll(jc)+all_orbit%ll(jd))) RETURN
  IF((ja == jb).AND.(MOD(jt,2)/=0)) RETURN      
  IF((jc == jd).AND.(MOD(jt,2)/=0)) RETURN       
  IF(triag(all_orbit%jj(ja),all_orbit%jj(jb),2*jt)) RETURN
  IF(triag(all_orbit%jj(jc),all_orbit%jj(jd),2*jt)) RETURN
  CALL mtx_elements(ja,jb,jc,jd,jt,gmtpn)
  !WRITE(6,*) "pphhmtx: ", gmtpn

END SUBROUTINE pphhmtx

