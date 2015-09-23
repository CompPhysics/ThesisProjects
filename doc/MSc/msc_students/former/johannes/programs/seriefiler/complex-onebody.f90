!             Program block complex-onebody.f
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
!             E-MAIL:   morten.hjorth-jensen@fys.uio.no
!             LANGUAGE: F90/F95  
!             LAST UPGRADE : February 2005
!
!             This file sets up one-body diagrams to second order 
!             for the sp self-consistency part. The present version includes
!             only the Hartree-Fock piece. The 2p-1h and 2h-1p piece will be
!             added upon later. Only steps 1 and 2 can be performed now.
!
!             The approach is as follows: 
!             step 1: compute the HF part with HO oscillator wfs for the
!             occupied states using a similarity transformed interaction
!             as starting point. The latter is computed in the CoM system.
!             step 2: generate new complex scaled sp orbits in momentum space
!             by diagonalizing the one-body problem. This serves simply to 
!             generate a new basis which is orthogonal. Iterate till self-consistency. 
!             step 3: This basis is then used to compute the 2p-1h and 2h-1p 
!             diagrams and the self-consistency procedure is repeated.
!             step 4: This basis can now be used to generate an effective 
!             two-body complex scaled interaction for shell-model studies.

!
!             The Hartree-Fock diagram in units of MeV^-2. Momenta are in units
!             of MeV.
!
!
!     Hartree-Hock self-consistency
!
SUBROUTINE hartree_fock(hf_iteration, sigma,selfconsistency )
  USE single_particle_orbits
  USE configurations
  USE constants
  USE wave_functions
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: sp_configs
  logical, intent(in) :: selfconsistency
  INTEGER :: a, c, ij, ipar, itz, ket , bra, number_of_iterations
  INTEGER :: la, lc, ja, jc, tza, tzc, iph
  REAL(DP) :: diaghf, pa, pc
  double precision, intent(out) :: sigma
  integer, intent(in) :: hf_iteration
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:) :: hf_mtx, hf_vect
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:) :: hf_eigen, e_temp

  !     loop over single-particle isospin projection (twice the value)
  WRITE(7,1000)
  DO itz= -1, 1, 2
     !     loop over single-particle parities
     DO ipar = 0, 1
        !     loop over single-particle angular momentum (twice the value)
        DO ij = 1, 10, 2
           !     find sp configurations with given j, parity and isospin
           CALL number_sp_confs(ij,ipar,itz,sp_configs)
           IF ( sp_configs%number_confs == 0 ) CYCLE
           ALLOCATE (sp_configs%config_ab(sp_configs%number_confs))
           CALL setup_sp_configurations(ij,ipar,itz,sp_configs)
           ALLOCATE &
                (hf_mtx(sp_configs%number_confs,sp_configs%number_confs) )
           ALLOCATE &
                (hf_vect(sp_configs%number_confs,sp_configs%number_confs) )
           ALLOCATE &
                (hf_eigen(sp_configs%number_confs) )
           ALLOCATE &
                (e_temp(sp_configs%number_confs) )
           WRITE(6,*) 'Number of sp configurations for Tz, Parity and J'
           WRITE(6,*) sp_configs%number_confs, itz, ipar, ij
           !     set up hf matrix to be diagonalized
           hf_mtx = 0.0_dp; hf_vect = 0.0_dp; hf_eigen = 0.0_dp!; sigma = 1E-4
           !number_of_iterations = 0
           !DO WHILE( (number_of_iterations < 1) .AND. (ABS(sigma) > 1E-5) )         
           e_temp = hf_eigen
           DO bra = 1, sp_configs%number_confs
              a = sp_configs%config_ab(bra)
              IF (all_orbit%model_space(a) == 'outside') CYCLE
              la = all_orbit%ll(a) ; ja = all_orbit%jj(a)
              pa = all_orbit%p(a)  ; tza = all_orbit%itzp(a)
              DO ket= bra, sp_configs%number_confs
                 c = sp_configs%config_ab(ket)
                 IF (all_orbit%model_space(c) == 'outside') CYCLE
                 lc = all_orbit%ll(c) ; jc = all_orbit%jj(c)
                 pc = all_orbit%p(c)  ; tzc = all_orbit%itzp(c)         
                 IF(iph(la) /= iph(lc)) CYCLE
                 IF ( ja /= jc) CYCLE; IF(tza /= tzc ) CYCLE
                 IF( ABS(pa-pc) >  2*hole_cutoff) CYCLE; diaghf = 0.0_dp
                 CALL hf(diaghf,pa,la,ja,tza,pc,lc,jc,tzc)
                 WRITE(7,1001) la, ja, tza,pa, lc, jc, tzc,pc, diaghf
                 
                 ! Interaction in units of MeV and preparing for diagonalization
                 hf_mtx(bra,ket) = pa*pc*diaghf*SQRT(all_orbit%w_p(a)*all_orbit%w_p(c))  
                 hf_mtx(ket,bra)=hf_mtx(bra,ket)
                 ! Add kinetic energy
                 IF ( bra == ket) THEN
                    hf_mtx(bra,bra)=hf_mtx(bra,bra)+pa*pa*0.5_dp/p_mass(tza)
                 ENDIF
              ENDDO
              WRITE(6,*) la, ja, tza,pa, lc, jc, tzc,pc, diaghf
           ENDDO
           ! diagonalize and get new eigenvalues and energies
           CALL diag_exact( hf_mtx, hf_vect, hf_eigen, sp_configs%number_confs)
           ! write out eigenvalues
           DO bra = 1, sp_configs%number_confs
              WRITE(6,*) bra, hf_eigen(bra)
           ENDDO
           sigma = 0.0_dp
           DO bra = 1, sp_configs%number_confs
              sigma = sigma +ABS( hf_eigen(bra)-e_temp(bra))
           ENDDO
           sigma = sigma/(sp_configs%number_confs*sp_configs%number_confs)
           WRITE(6,*) 'Sigma for this iteration', sigma, number_of_iterations
           number_of_iterations = number_of_iterations +1 
           
           DEALLOCATE ( hf_mtx); DEALLOCATE ( hf_eigen)
           DEALLOCATE ( hf_vect); DEALLOCATE ( e_temp)
        ENDDO
        
     ENDDO
  ENDDO

1000 FORMAT(//2X,4H--->,2X,36H Momenta in MeV HF diagram in MeV^-2,2X,4H<---//)
1001 FORMAT(2X,12HA: l, j, tz=,3I3,2X,3HPA:,E12.6,2X,12HC: l, j, tz=,3I3,2X,3HPC:,E12.6,2X,3HHF=,E12.6)

END SUBROUTINE hartree_fock

!
!                  SUBROUTINE to calculate the Hartree-Fock contrib.    
!                  Dimensionality = MeV^-2, as the interaction      
!

SUBROUTINE hf(hfe,pa,la,ja,tza,pc,lc,jc,tzc)
  USE constants
  USE single_particle_orbits
  USE wave_functions  
  IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: hfe
  REAL(DP), INTENT(IN) :: pa, pc
  INTEGER, INTENT(IN) :: ja, jc, la, lc, tza, tzc
  INTEGER :: h,nh, lh, jh, tzh, j_min, j_max, jph, itz, iq
  REAL(DP) :: vmtx, qmin, qmax, sum_hf, ph, zh, ptest
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: q_points, weight_q

  hfe=0.0_dp
  qmin = 0.0_dp; qmax = hole_cutoff
  IF (ABS(pa-pc) - hole_cutoff > 0.) qmin = (ABS(pa-pc) - hole_cutoff)
  ALLOCATE ( q_points(int_points), weight_q(int_points) )
  CALL gauss_legendre(qmin,qmax,q_points,weight_q,int_points)
  ptest = all_orbit%p(1)
  DO h=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE
     lh = all_orbit%ll(h) ; jh = all_orbit%jj(h)
     nh = all_orbit%nn(h)  ; tzh = all_orbit%itzp(h)
     IF (all_orbit%p(h) /= ptest) CYCLE   ! only one hole state, independent of p
     j_min=ABS((ja-all_orbit%jj(h))/2)
     j_max=(ja+all_orbit%jj(h))/2
     IF ( tza+tzh /= tzc+tzh) CYCLE; itz = (tza+tzh)/2
     DO jph=j_min,j_max
        sum_hf = 0.0_dp
        ! now integrate over the HO wf for the  ket side
        DO iq=1, int_points
           ph = q_points(iq)
           IF (ABS(pc-ph)-pa > hole_cutoff ) CYCLE 
           IF (pa-ph-pc > hole_cutoff ) CYCLE 
           zh = ph*oscillator_length
           !  Get <pa la ja nh lh jh| v |pc lc jc ph lh jh> 
           CALL leftside_transf(vmtx,pa,la,ja,nh,lh,jh,pc,lc,jc,ph,lh,jh,itz,jph)
           sum_hf=sum_hf+weight_q(iq)*dble( rnl(nh,lh,dcmplx(zh)) ) &
                         *ph*ph*SQRT(oscillator_length**3)*vmtx
        END DO
        hfe=hfe+sum_hf*(2*jph+1.0_dp)
     ENDDO
  ENDDO
  hfe=hfe/(ja+1.0_dp)
  DEALLOCATE ( weight_q, q_points)

END SUBROUTINE hf

!
!     Calculates the matrix elements
!     <pa la ja pb lb jb | V | pc lc jc pd ld jd  >
!     using as input the matrix elements
!     <pa la ja nb lb jb | V | pc lc jc pd ld jd  >
!     Dimensionality = MeV^-2, as the interaction      
!

SUBROUTINE leftside_transf(sum_v,pa,la,ja,nb,lb,jb,pc,lc,jc,pd,ld,jd,itz,jph)
  USE constants
  USE single_particle_orbits
    USE wave_functions
  IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: sum_v
  REAL(DP), INTENT(IN) :: pa, pc, pd
  INTEGER, INTENT(IN) :: ja, jc, la, lc, nb, lb, jb, ld, jd, jph, itz
  INTEGER :: iq
  REAL(DP) :: vmtx, qmin, qmax, pb, zb
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: q_points, weight_q

  ! Finding integration limits
  qmin = 0.0_dp; 
  IF (ABS(pc-pd) - pa > 0.0_dp ) qmin = ABS(pc-pd) - pa 
  IF (pa-pc-pd > qmin ) qmin = pa-pc-pd
  qmax = hole_cutoff
  IF (pa+pc+pd < qmax) qmax = pa+pc+pd
  ALLOCATE ( q_points(int_points), weight_q(int_points) )
  CALL gauss_legendre(qmin,qmax,q_points,weight_q,int_points)
  sum_v = 0.0_dp
  ! now integrate over the HO for the  bra side
  DO iq=1, int_points
     pb = q_points(iq)
     IF ( (ABS(pa-pb) > pc+pd) .OR. (ABS(pc-pd) > pa+pb) ) CYCLE
     zb = pb*oscillator_length
     !  Get <pa la ja pb lb jb| v |pc lc jc pd ld jd> 
     CALL v_labmomentum_space(vmtx,la,ja,pa,lb,jb,pb,lc,jc,pc,ld,jd,pd,itz,jph)
     sum_v=sum_v+weight_q(iq)*dble( rnl(nb,lb,dcmplx(zb)) ) *SQRT(oscillator_length**3)*pb*pb*vmtx
  END DO
  DEALLOCATE ( weight_q, q_points)

END SUBROUTINE leftside_transf

!
!     Calculates the matrix elements
!     <pa la ja pb lb jb J T_z | V |  pc lc jc pd ld jd J T_z>
!     Present version for identical masses only.
!     Dimensionality is energy^-5, with energy in units of MeV

SUBROUTINE v_labmomentum_space(v_pot,la,ja,pa,lb,jb,pb,lc,jc,pc,   &
     ld,jd,pd,isospin_tz,j_lab)
  USE configurations
  USE constants
  USE single_particle_orbits
  USE partial_waves
  USE relcm_gmatrix
  USE wave_functions
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: isospin_tz, j_lab, ja, jb, jc, jd, la, lb, lc, ld
  INTEGER :: llmin, llmax, lambda_min, lambda_max, ncoup, channel
  INTEGER :: spin, l, jsmin, jsmax, ll, js, lambda_ket
  INTEGER :: lp1, lp2, lam1, lam2, lambda_bra, iq, iph, lp, is1, is2
  REAL(DP), INTENT(IN) :: pa, pb, pc, pd
  REAL(DP), INTENT(INOUT) :: v_pot
  REAL(DP) :: p_cd, p_ab, q_ab, q_cd, q_rel, k_rel, qmin, qmax, kinetic_ket
  REAL(DP) :: sum6j, sum9j, lsjj_coeff_bra, lsjj_coeff_ket, int_sum
  REAL(DP) :: vb_bra, vb_ket, kinetic_bra, pos, p_com, vv
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: q_points, weight_q
  LOGICAL triag

  ALLOCATE ( q_points(int_points), weight_q(int_points) )
  v_pot = 0.0D0; kinetic_ket =  0.5*(pc*pc+pd*pd)
  kinetic_bra =  0.5*(pa*pa+pb*pb)
  ! Tests of momenta relation needed for the computation of vector brackets
  ! Max values of two-body momenta
  p_cd = pc+pd; p_ab = pa+pb
  ! Min value of two-body momenta
  q_cd=ABS(pc-pd); q_ab=ABS(pa-pb)
  ! Now we set up the integration limits for the relative momentum q_rel
  ! of the ket state | q l P L S, Js T_z> performing tests based on max and min
  ! momenta
  IF(p_ab >=  p_cd) qmin = 0.5D0*q_cd
  IF(p_cd > p_ab) qmin = SQRT(kinetic_ket-.25D0*p_ab*p_ab)
  IF(q_cd >= q_ab) qmax = 0.5D0*p_cd
  IF(q_ab > q_cd) qmax = SQRT(kinetic_ket-0.25D0*q_ab*q_ab)

  ! adjust max relative momentum (qmax) for renormalized interaction
  ! since V_lowk(k,k') = 0 for k,k' > kcut = abs(krel(n_k1)) 
  SELECT CASE(type_of_interaction)
  CASE('free')
     qmax = qmax
  CASE('renorm')
     if ( qmax > abs( krel(n_k1) ) ) then
        qmax = abs( krel(n_k1) )
     else
        qmax = qmax 
     end if
  END SELECT

  
  CALL gauss_legendre(qmin,qmax,q_points,weight_q,int_points)
  ! Start performing the transformation for CoM and relative coordinates
  ! to lab frame, ket side first
  is1 = 1; is2=1
  DO spin= 0, 1  ! S = 0 or 1 for nucleons only
     DO l=0, jmax  ! relative orbital momentum
        jsmin=ABS(l-spin); jsmax=l+spin
        IF ( (ABS(isospin_tz) == 1).AND.(MOD(l+spin+ABS(isospin_tz),2)==0)) CYCLE
        DO js=jsmin,jsmax    ! Relative angular momentum
           IF(js > jmax) CYCLE 
           llmin=ABS(j_lab-js); llmax=j_lab+js
           DO ll=llmin,llmax ! CoM orbital momentum
              ! parity test   
              IF(MOD(l+ll+lc+ld,2) /= 0) CYCLE    ! parity check ket side
              lambda_min=ABS(lc-ld);  lambda_max=lc+ld
              DO lambda_ket = lambda_min,lambda_max   
                 IF(triag(lambda_ket,l,ll)) CYCLE
                 IF(triag(lambda_ket,j_lab,spin)) CYCLE
                 lsjj_coeff_ket = sum6j(ll,l,lambda_ket,spin,j_lab,js)* &
                      sum9j(jc,jd,lc,ld,is1,is2,j_lab,lambda_ket,spin)
                 IF (ABS(lsjj_coeff_ket) < epsilon) CYCLE
                 ! Transformation for CoM and relative coordinates
                 ! to lab frame, bra side now
                 lp1=ABS(js-spin); lp2=js+spin
                 lam1=ABS(j_lab-spin); lam2=j_lab+spin
                 ! The integral over momenta is the last loop
                 DO iq=1, int_points
                    q_rel = q_points(iq)
                    pos = 4.D0*(kinetic_ket-q_rel*q_rel)
                    IF(pos <=  0.D0) CYCLE; p_com = SQRT(pos)
                    IF(p_com  > p_ab) CYCLE; IF(p_com  < q_ab) CYCLE
                    IF(p_com  > p_cd) CYCLE; IF(p_com  < q_cd) CYCLE
                    pos  = kinetic_bra - 0.25D0*p_com*p_com
                    IF ( pos < 0.D0 ) CYCLE; k_rel = SQRT(pos)  
                    vb_ket = vector_trcoefficients(q_rel,p_com,pc,pd,l,ll,&
                             lambda_ket,lc,ld,isospin_tz)/(pc*pd)
                    IF (ABS(vb_ket) < epsilon) CYCLE
                    int_sum = 0.0_dp
                    DO lp=lp1,lp2
                       IF ((ABS(isospin_tz) == 1).AND.(MOD(lp+spin+ABS(isospin_tz),2)==0)) CYCLE
                       IF(lp > lmax) CYCLE
                       IF(MOD(lp+ll+la+lb,2) /= 0) CYCLE         !  parity check bra side
                       !  get the interaction in the CoM coordinates
                       IF ( ( lp == l) .AND. ( lp == js)) THEN
                          ncoup = 1
                       ELSEIF ( ( lp /= js) .OR. ( l /= js)) THEN
                          ncoup = 2
                       ENDIF
                       CALL get_channel(lp,l,js,spin,isospin_tz,ncoup,channel)
                       IF ( channel == 0) CYCLE
                       CALL g_interpolate(channel, ncoup, lp,l,js,k_rel,q_rel,vv)
!                       CALL v_bare(lp,k_rel,l,q_rel,spin,js,isospin_tz,vv)
                       vv = vv*iph(ABS((l-lp))/2)
                       DO lambda_bra=lam1,lam2
                          IF(triag(lambda_bra,la,lb)) CYCLE           ! triangular rels.
                          IF(triag(lambda_bra,lp,ll)) CYCLE
                          lsjj_coeff_bra = sum9j(ja,jb,la,lb,is1,is2,j_lab,lambda_bra,spin)* &
                               sum6j(ll,lp,lambda_bra,spin,j_lab,js)
                          IF (ABS(lsjj_coeff_bra) < epsilon) CYCLE
                          vb_bra = vector_trcoefficients &
                           (k_rel,p_com,pa,pb,lp,ll,lambda_bra,la,lb,isospin_tz)/(pa*pb)
                          IF (ABS(vb_bra) < epsilon) CYCLE
                          int_sum = int_sum+lsjj_coeff_bra*vv*vb_bra
                       ENDDO
                    ENDDO
                    v_pot=v_pot+int_sum*vb_ket*lsjj_coeff_ket*weight_q(iq)*q_rel/p_com
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  DEALLOCATE ( weight_q, q_points)

END SUBROUTINE v_labmomentum_space


REAL(DP) FUNCTION sum9j(j1,j2,l1,l2,is1,is2,jt,lam,is)
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: j1,j2,l1,l2,is1,is2,jt,lam,is
  REAL(DP) :: s1, s2, sls

  s1=SQRT((1.*j1+1.)*(1.*j2+1.)*(2.*is+1.))
  sls=snj(2*l1,2*l2,2*lam,is1,is2,2*is,j1,j2,2*jt)     
  s2=2.*lam+1.
  sum9j=sls*s1*s2

END FUNCTION sum9j



REAL(DP) FUNCTION sum6j(lc,l,lam,is,jt,js)
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: lc,l,lam,is,jt,js
  REAL(DP) :: s3, sign, sx

  s3=SQRT(2.*js+1.)
  sign = 1.- 2.*MOD(lam+js+lc+is,2)
  sx=sjs(2*lc,2*l,2*lam,2*is,2*jt,2*js)
  sum6j=sign*s3*sx

END FUNCTION sum6j


SUBROUTINE g_interpolate(channel,ncoup,la,lb,jang,k_rel,q_rel,vv)
  USE wave_functions
  USE constants
  USE relcm_gmatrix
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: k_rel, q_rel
  REAL(DP), INTENT(INOUT) :: vv
  REAL(DP)  :: v_interpol
  INTEGER, INTENT(IN) :: ncoup, channel, la, lb, jang
  INTEGER :: lim1, lim2, lim3, lim4, dim

  !  set in test here so that n_k = total in case of free and n_k1 in case of not free
  SELECT CASE(type_of_interaction)
  CASE('free')
     dim = n_k
  CASE('renorm')
     dim = n_k1
  END SELECT
  vv = 0.0_dp; v_interpol = 0.0_dp
  IF ( ncoup == 1) THEN
     lim1=1; lim2=dim ; lim3=1 ; lim4=dim
  ELSEIF ( ncoup == 2 ) THEN
     IF ( (la == lb).AND. ( jang > la) ) THEN
        lim1=1; lim2=dim ; lim3=1 ; lim4=dim
     ELSEIF ( (la == lb).AND. ( jang < la) ) THEN
        lim1=1+dim; lim2=dim+dim ; lim3=1+dim ; lim4=dim+dim
     ELSEIF ( la >  lb ) THEN
        lim1=1+dim; lim2=dim+dim ; lim3=1 ; lim4=dim
     ELSEIF ( la <  lb ) THEN
        lim1=1; lim2=dim ; lim3=1+dim ; lim4=dim+dim
     ENDIF
  ENDIF
  CALL lagrange_2dim(k_rel,q_rel,REAL(vfree(lim1:lim2,lim3:lim4,channel)), &
                     v_interpol,dim,REAL(krel))
  vv = v_interpol

END SUBROUTINE g_interpolate



SUBROUTINE  get_channel(lp,l,js,spin,isospin_tz,ncoup,channel)
  USE constants
  USE partial_waves
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: isospin_tz, spin, js, lp, l, ncoup
  INTEGER, INTENT(INOUT) :: channel
  INTEGER :: loop, la, lb
  
  channel = 0
  IF ( ncoup == 1) THEN
     la = js; lb = js
  ENDIF
  IF ( ncoup == 2) THEN
     IF ( js > 0) THEN
        la = js-spin; lb = js+spin
     ELSEIF( js == 0 ) THEN
        la = js+spin; lb = js+spin
     ENDIF
  ENDIF
  DO loop=1, no_channels
     IF ((orb_lrel_max(loop) == lb) .AND. (orb_lrel_min(loop) == la) .AND. &
         ( spin_rel(loop) == spin) .AND. &
         (jang_rel(loop) == js) .AND. (iso(loop) == isospin_tz ) ) THEN 
         channel = loop
         EXIT
     ENDIF
  ENDDO

END SUBROUTINE  get_channel
