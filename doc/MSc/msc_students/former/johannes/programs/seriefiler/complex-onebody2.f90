
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
SUBROUTINE hartree_fock(iteration_number, sigma,selfconsistency )
  USE single_particle_orbits
  USE configurations
  USE constants
  use hf_constants
  USE wave_functions
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: sp_configs
  INTEGER :: a, c, ij, ipar, itz, ket , bra, number_of_iterations, number_of_holes
  INTEGER :: la, lc, ja, jc, tza, tzc, iph
  REAL(DP) ::  pa, pc, kk, energy_check,norm_check
  complex(dpc) :: diaghf, klab, norm_sum,zh
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:) :: hf_mtx, hf_vect
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:) :: hf_eigen, e_temp
  integer :: i, number_mtxel
  integer, intent(in) :: iteration_number
  logical, intent(in) :: selfconsistency
  double precision, intent(out) :: sigma
  
  
  !call interpol_1dim( 0,0,1,-1,dcmplx(10.,0.d0), diaghf)
  !write(6,*) rnl(0,0,zh)*dSQRT(oscillator_length**3),dSQRT(oscillator_length**3)
  !write(6,*) n_lab, int_points
  !  do kk= 30.d0, 100.d0, 10.d0
  !     klab = dcmplx(kk,0.d0) 
  
  !     call interpol_1dim2( 0,0,1,-1, klab,diaghf)
  !     call interpol_1dim2( 0,0,1,-1,klab, diaghf)
  !     zh =   klab*oscillator_length
  !     write(6,*) kk,dble(diaghf), rnl(0,0,zh)*dSQRT(oscillator_length**3) 
  !  end do
  !  return
  !  
  !     loop over single-particle isospin projection (twice the value)
  !WRITE(7,1000)
  

  !
  ! for hf_iteration > 1 no longer woods-saxon basis
  ! setup new limits for interpolation routine. 
  ! 
  !if ( iteration_number > 1 ) then
  k_interpol_min = k_mesh(1)
  k_interpol_max = k_mesh(n_lab) 
  !end if
  
  
  
  sigma = 0.0_dp
  v_bhf = 0.d0
  number_mtxel = 0
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
           hf_mtx = 0.0_dp; hf_vect = 0.0_dp; hf_eigen = 0.0_dp
                         
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
                 !if ( itz == 1 ) cycle
                 
                 CALL hf(diaghf,pa,la,ja,tza,pc,lc,jc,tzc, iteration_number, number_mtxel)
                 
                 v_bhf( bra, ket, la,(ja+1)/2 ,(tza+1)/2 ) = dble( diaghf ) 
                 if (  selfconsistency ) then 
                    WRITE(10,*) bra, la, ja, tza,pa, ket, lc, jc, tzc,pc, real(diaghf)
                 END if
                 
                 ! Interaction in units of MeV and preparing for diagonalization
                 hf_mtx(bra,ket) = phase**3 * pa*pc*diaghf*SQRT(all_orbit%w_p(a)*all_orbit%w_p(c))  
                 !hf_mtx(bra,ket) = phase * diaghf*SQRT(all_orbit%w_p(a)*all_orbit%w_p(c))  
                 hf_mtx(ket,bra)=hf_mtx(bra,ket)
                 ! Add kinetic energy
                 IF ( bra == ket) THEN
                    write(6,*) pa, dble(diaghf),aimag(diaghf)
                    hf_mtx(bra,bra)=hf_mtx(bra,bra) + phase**2*pa*pa*0.5_dp/p_mass(tza)
                 ENDIF
              ENDDO
           ENDDO
           ! diagonalize and get new eigenvalues and energies
           CALL diag_exact( hf_mtx, hf_vect, hf_eigen, sp_configs%number_confs)
           
           
           !    set up the new bhf wave function
           bhf_hol_temp(:,la,(ij+1)/2,(tza+1)/2) = hf_vect(:,1)
           
           energy_check = dot_product(bhf_hol_temp(:,la,(ij+1)/2,(tza+1)/2), & 
                matmul( hf_mtx, bhf_hol_temp(:,la,(ij+1)/2,(tza+1)/2) ))
           
           !    re-normalize bhf wave function; sum( bhf**2*k_mesh**2*w_mesh = 1 ) 
           bhf_hol_temp(:,la,(ij+1)/2,(tza+1)/2) = bhf_hol_temp(:,la,(ij+1)/2,(tza+1)/2)/k_mesh(:)/sqrt(w_mesh(:)) 
           norm_check   = dot_product( k_mesh(:)**2*w_mesh(:)*bhf_hol_temp(:,la,(ij+1)/2,(tza+1)/2), & 
                bhf_hol_temp(:,la,(ij+1)/2,(tza+1)/2))
           
           
           if ( (.not.  selfconsistency ) ) then 
              write(6,*) 'sp quantum numbers:',itz, ij,la 
              write(6,*) 'check energy and norm of bhf_hol for this iteration:', iteration_number, energy_check, norm_check 
           END if
              
           if (  selfconsistency ) then 
              write(6,*) 'sp quantum numbers:',itz, ij,la 
              write(6,*) 'check energy and norm of bhf_hol:', energy_check, norm_check 
              DO bra = 1, n_lab
                 write(6,*) 'Single-particle energy:',hf_eigen(bra)
              ENDDO
           end if
           
           ! calculate variance of bhf wave function
           DO bra = 1, n_lab
              sigma = sigma + ABS( k_mesh(bra)*sqrt(w_mesh(bra))* & 
                   (bhf_hol_temp(bra,la,(ij+1)/2,(tza+1)/2) - bhf_hol(bra,la,(ij+1)/2,(tza+1)/2) ) )
           ENDDO
           
           
           DEALLOCATE ( hf_mtx); DEALLOCATE ( hf_eigen)
           DEALLOCATE ( hf_vect); DEALLOCATE ( e_temp)
           
        ENDDO
     ENDDO
  ENDDO
  
  ! allocate space for momentum space interaction
  !if ( iteration_number == 1 ) allocate( vlab_mom(number_mtxel) )
  
  write(6,*) 'number of matrix elements in this iteration',number_mtxel, hf_iterations
  

  ! deallocate matrix for hole functions
  !call deallocate_hf_holes


  !1000 FORMAT(//2X,4H--->,2X,36H Momenta in MeV HF diagram in MeV^-2,2X,4H<---//)
  !1001 FORMAT(2X,12HA: l, j, tz=,3I3,2X,3HPA:,E12.6,2X,12HC: l, j, tz=,3I3,2X,3HPC:,E12.6,2X,3HHF=,E12.6)

END SUBROUTINE hartree_fock

!
!                  SUBROUTINE to calculate the Hartree-Fock contrib.    
!                  Dimensionality = MeV^-2, as the interaction      
!

SUBROUTINE hf(hfe,pa,la,ja,tza,pc,lc,jc,tzc, iteration_number, number_mtxel)
  USE constants
  use hf_constants
  USE single_particle_orbits
  USE wave_functions  
  IMPLICIT NONE
  complex(DPC), INTENT(INOUT) :: hfe
  REAL(DP), INTENT(IN) :: pa, pc
  INTEGER, INTENT(IN) :: ja, jc, la, lc, tza, tzc, iteration_number
  integer, intent(out) :: number_mtxel
  INTEGER :: h,nh, lh, jh, tzh, j_min, j_max, jph, itz, iq
  REAL(DP) :: qmin, qmax,  ph, ptest
  complex(dpc) :: zh
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: q_points, weight_q
  complex(dpc) :: vmtx,sum_hf, hole_func

  hfe=0.0_dp
  qmin = 0.0_dp; qmax = hole_cutoff
  IF (ABS(pa-pc) - hole_cutoff > 0.) qmin = (ABS(pa-pc) - hole_cutoff)
  ALLOCATE ( q_points(int_points), weight_q(int_points) )
  
  ! check that integration points are within the klab  limits 
  ! for the bhf wave function -> no extrapolation!
  ! if not -> new integration limits.
  if ( qmin < k_interpol_min ) qmin = k_interpol_min
  if ( qmax > k_interpol_max ) qmax = k_interpol_max
  
  q_points(:) = 0.5d0*( qmax+qmin ) + 0.5d0*( qmax-qmin )*x_mesh(:)
  weight_q(:) = 0.5d0* ( qmax-qmin ) * wx_mesh(:)
    
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
           zh = phase * ph*oscillator_length
           
           
           ! in first iteration use woods-saxon basis 
           ! after that we interpolate the new bhf wave function on the integration points
           !if ( iteration_number == 1 ) then
           !   call interpol_1dim( nh,lh,jh,tzh, phase*ph, hole_func)
           !elseif (iteration_number > 1 ) then
           call interpol_1dim2( nh,lh,jh,tzh, phase*ph, hole_func)
           !end if
           !  Get <pa la ja nh lh jh| v |pc lc jc ph lh jh> 
           CALL leftside_transf(vmtx,pa,la,ja,nh,lh,jh,tzh,pc,lc,jc,ph,lh,jh,itz,jph,& 
                iteration_number, number_mtxel)
           !sum_hf=sum_hf+phase**3*weight_q(iq)*rnl(nh,lh,zh) &
           !     *ph*ph*SQRT(oscillator_length**3)*vmtx
           
           sum_hf=sum_hf+phase**3*weight_q(iq)*ph*ph*hole_func*vmtx
        
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

SUBROUTINE leftside_transf(sum_v,pa,la,ja,nb,lb,jb,tzb,pc,lc,jc,pd,ld,jd,itz,jph, iteration_number, number_mtxel)
  USE constants
  use hf_constants
  USE single_particle_orbits
  USE wave_functions
  
  IMPLICIT NONE
  complex(DPC), INTENT(INOUT) :: sum_v
  REAL(DP), INTENT(IN) :: pa, pc, pd
  INTEGER, INTENT(IN) :: ja, jc, la, lc, nb, lb, jb, tzb,ld, jd, jph, itz, iteration_number
  integer, intent(out) :: number_mtxel
  INTEGER :: iq
  REAL(DP) :: qmin, qmax, pb
  complex(dpc) :: vmtx, zb, hole_func
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: q_points, weight_q

  ! Finding integration limits
  qmin = 0.0_dp; 
  IF (ABS(pc-pd) - pa > 0.0_dp ) qmin = ABS(pc-pd) - pa 
  IF (pa-pc-pd > qmin ) qmin = pa-pc-pd
  qmax = hole_cutoff
  IF (pa+pc+pd < qmax) qmax = pa+pc+pd
  
  if ( qmin < k_interpol_min ) qmin = k_interpol_min
  if ( qmax > k_interpol_max ) qmax = k_interpol_max
 
  ALLOCATE ( q_points(int_points), weight_q(int_points) )
  q_points(:) = 0.5d0*( qmax+qmin ) + 0.5d0*( qmax-qmin )*x_mesh(:)
  weight_q(:) = 0.5d0* ( qmax-qmin ) * wx_mesh(:)
 
  !CALL gauss_legendre(qmin,qmax,q_points,weight_q,int_points)
  sum_v = 0.0_dp
  ! now integrate over the HO for the  bra side
  DO iq=1, int_points
     pb = q_points(iq)
     IF ( (ABS(pa-pb) > pc+pd) .OR. (ABS(pc-pd) > pa+pb) ) CYCLE
     zb = phase * pb*oscillator_length


     !if ( iteration_number == 1 ) then
     !   call interpol_1dim( nb,lb,jb,tzb, phase*pb, hole_func)
     !elseif ( iteration_number > 1 ) then
     call interpol_1dim2( nb,lb,jb,tzb, phase*pb, hole_func)
     !end if
         
     
     number_mtxel = number_mtxel +1 
     
     !  Get <pa la ja pb lb jb| v |pc lc jc pd ld jd> 
     if ( iteration_number == 1 ) then
        CALL v_labmomentum_space(vmtx,la,ja,pa,lb,jb,pb,lc,jc,pc,ld,jd,pd,itz,jph)
        vlab_mom_re(number_mtxel) = real( vmtx )
        vlab_mom_im(number_mtxel) = real(aimag(vmtx))
     !elseif ( iteration_number == 2 ) then
     !   CALL v_labmomentum_space(vmtx,la,ja,pa,lb,jb,pb,lc,jc,pc,ld,jd,pd,itz,jph)
     !   vlab_mom(number_mtxel) = real( vmtx )  
     !elseif ( iteration_number > 2 ) then
     !   vmtx = dcmplx( vlab_mom(number_mtxel) ) 
     !end if
     elseif ( number_mtxel > 1 ) then
        vmtx = dcmplx( vlab_mom_re(number_mtxel), 0.d0 ) 
     end if

     !sum_v=sum_v+phase**3*weight_q(iq)*rnl(nb,lb,zb)*SQRT(oscillator_length**3)*pb*pb*vmtx
     sum_v=sum_v+phase**3*weight_q(iq)*hole_func*pb*pb*vmtx
  END DO
  DEALLOCATE ( weight_q, q_points)

END SUBROUTINE leftside_transf

!
!     Calculates the matrix elements
!     <pa la ja pb lb jb J T_z | V |  pc lc jc pd ld jd J T_z>
!     Present version for identical masses only.
!     Dimensionality is energy^-5, with energy in units of MeV
!     Note that all momenta are in absolute values so all
!     momenta should be scaled by exp(-i*theta)
!
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
  COMPLEX(DPC) :: vv_cmplx
  complex(DPC), INTENT(INOUT) :: v_pot
  REAL(DP) :: p_cd, p_ab, q_ab, q_cd, q_rel, k_rel, qmin, qmax, kinetic_ket
  REAL(DP) :: sum6j, sum9j, lsjj_coeff_bra, lsjj_coeff_ket
  complex(dpc) :: int_sum
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

  !q_points(:) = 0.5d0*( qmax+qmin ) + 0.5d0*( qmax-qmin )*x_mesh(:)
  !weight_q(:) = 0.5d0* ( qmax-qmin ) * wx_mesh(:)
 
  q_points(:) = (0.5d0*(qmax-qmin)*(t_mesh(:)+2.d0*qmax/(qmax-qmin) -1.d0))
  weight_q(:) = (0.5d0*(qmax-qmin)*wt_mesh(:))* dsqrt( 1.d0 - t_mesh(:)*t_mesh(:) )
  !q_rel = phase* dcmplx(0.5d0*(qmax-qmin)*(t(iq)+2.d0*qmax/(qmax-qmin) -1.d0))
  !dq_rel = phase*dcmplx(0.5d0*(qmax-qmin)*w_t(iq))
  !cheb_weight = dsqrt( 1.d0 - t(iq)*t(iq) )
                    

  !CALL gauss_legendre(qmin,qmax,q_points,weight_q,int_points)
  ! Start performing the transformation for CoM and relative coordinates
  ! to lab frame, ket side first
  is1 = 1; is2=1
  DO spin= 0, 1  ! S = 0 or 1 for nucleons only
     DO l=0, lmax  ! relative orbital momentum, (jmax changed to lmax)
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
                    IF(pos <=  0.D0) CYCLE; p_com = SQrt(pos)
                    IF(p_com  > p_ab) CYCLE; IF(p_com  < q_ab) CYCLE
                    IF(p_com  > p_cd) CYCLE; IF(p_com  < q_cd) CYCLE
                    pos  = kinetic_bra - 0.25D0*p_com*p_com
                    IF ( pos < 0.D0 ) CYCLE; k_rel = SQRT(pos)  
                    vb_ket = vector_trcoefficients(q_rel,p_com,pc,pd,l,ll,&
                             lambda_ket,lc,ld,isospin_tz) /(pc*pd)
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
                       
                       !if ( l > 1 ) cycle
                       !if ( lp /= l ) cycle
                       !call malfliet_tjon( vv_cmplx, l, phase*k_rel, phase*q_rel ) 
                       
                       CALL g_interpolate(channel, ncoup, lp,l,js,k_rel,q_rel,vv_cmplx)
                       !CALL v_bare(lp,k_rel,l,q_rel,spin,js,isospin_tz,vv_cmplx)
                       vv_cmplx = vv_cmplx*iph(ABS((l-lp))/2)
                       DO lambda_bra=lam1,lam2
                          IF(triag(lambda_bra,la,lb)) CYCLE           ! triangular rels.
                          IF(triag(lambda_bra,lp,ll)) CYCLE
                          lsjj_coeff_bra = sum9j(ja,jb,la,lb,is1,is2,j_lab,lambda_bra,spin)* &
                               sum6j(ll,lp,lambda_bra,spin,j_lab,js)
                          IF (ABS(lsjj_coeff_bra) < epsilon) CYCLE
                          vb_bra = vector_trcoefficients &
                           (k_rel,p_com,pa,pb,lp,ll,lambda_bra,la,lb,isospin_tz) /(pa*pb)
                          IF (ABS(vb_bra) < epsilon) CYCLE
                          int_sum = int_sum+lsjj_coeff_bra*vv_cmplx*vb_bra
                       ENDDO
                    ENDDO
                    !int_sum = int_sum+cheb_weight*dq_rel*q_rel*vb_bra*vb_ket*vv/p_com
                    v_pot=v_pot+ int_sum*vb_ket*lsjj_coeff_ket*weight_q(iq)*q_rel/p_com/phase**3
                 ENDDO
              ENDDO
              
           ENDDO
        ENDDO
     ENDDO
  end DO
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
  ! REAL(DP), INTENT(INOUT) :: vv
  COMPLEX(DPC), INTENT(INOUT) :: vv
  REAL(DP)  :: v_interpol_re, v_interpol_im
  INTEGER, INTENT(IN) :: ncoup, channel, la, lb, jang
  INTEGER :: lim1, lim2, lim3, lim4, dim

  !  set in test here so that n_k = total in case of free and n_k1 in case of not free
  SELECT CASE(type_of_interaction)
  CASE('free')
     dim = n_k
  CASE('renorm')
     dim = n_k1
  END SELECT
  vv = 0.0_dp; v_interpol_re = 0.0_dp; v_interpol_im = 0.0_dp
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
  
  if ( phase == 1.d0 ) then
     CALL lagrange_2dim(k_rel,q_rel,REAL(vfree(lim1:lim2,lim3:lim4,channel)), &
          v_interpol_re,dim,abs(krel))
  
  else
     CALL lagrange_2dim(k_rel,q_rel,REAL(vfree(lim1:lim2,lim3:lim4,channel)), &
          v_interpol_re,dim,abs(krel))
  
     CALL lagrange_2dim(k_rel,q_rel,aimag(vfree(lim1:lim2,lim3:lim4,channel)), &
          v_interpol_im,dim,abs(krel))

  end if
  
  vv = dcmplx( v_interpol_re, v_interpol_im )

!  if ( channel == 1 .and. k_rel == q_rel ) then
!     write(6,*) k_rel, q_rel, dble(vv),aimag(vv)
!  end if

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

subroutine malfliet_tjon( mt_pot, l, k, q ) 
  implicit none
  complex*16, intent(in) :: k, q
  complex*16, intent(out) :: mt_pot
  integer, intent(in) :: l
  double precision :: V_r, V_a, mu_a, mu_r
  double precision :: a4, b4, a2, b2, a6, b6 
  complex*16 :: k4,k2, x1, x2, x, y, y1, y2, legendre1, legendre2

  ! Mt-potential parameters
  V_r =  7.291d0 
  V_a =  -2.!-2.6047d0 !5.42      
  mu_r = 613.69d0  
  mu_a = 305.86d0  


  x = ( k * k + q * q + mu_r*mu_r ) / ( 2.d0 * k * q  )
  x1 = k * k + q * q + 2.d0 * k * q + mu_r*mu_r
  x2 = k * k + q * q - 2.d0 * k * q + mu_r*mu_r
  
  y = ( k * k + q * q + mu_a*mu_a ) / ( 2.d0 * k * q  )
  y1 = k * k + q * q + 2.d0 * k * q + mu_a*mu_a
  y2 = k * k + q * q - 2.d0 * k * q + mu_a*mu_a
  
  IF ( l == 0 ) THEN
     legendre1 = 0.5d0 * log(x1/x2)
     legendre2 = 0.5d0 * log(y1/y2)
  ELSEIF ( l == 1 ) THEN
     legendre1 = ( 0.5d0 * x * log(x1/x2) - 1.d0 )
     legendre2 = ( 0.5d0 * y * log(y1/y2) - 1.d0 )
  END IF
  


  mt_pot = (V_r * ( 1.d0 / ( 2.d0 * k * q ) ) * &
       legendre1  + &
       V_a * ( 1.d0 /( 2.d0 * k * q ) ) * &
       legendre2  )

end subroutine malfliet_tjon



!
!          This routine contains old-fashioned common blocks
!          acting as interface between the effective interaction program
!          and the potential routines of the Bonn type, Nijmegen or 
!          Argonne groups. 
!          v is in units of MeV^-2,  
!          In partial wave basis, v(6) 
!          means
!          v(1)= singlet uncoupled
!          V(2)= triplet uncoupled
!          V(3)= triplet coupled, < J+1 | V | J+1>
!          V(4)= triplet coupled, < J-1 | V | J-1>
!          V(5)= triplet coupled, < J+1 | V | J-1>
!          V(6)= triplet coupled, < J-1 | V | J+1>
!
!          v_bare(lp, k_rel,l,q_rel,spin,js,isospin_tz,vv)
SUBROUTINE v_bare(lp,k_rel,l,q_rel,spin,js,isospin_tz,vv)
  USE constants
  IMPLICIT NONE
  !double precision, INTENT(INOUT) :: vv
  complex*16, INTENT(INOUT) :: vv
  INTEGER, INTENT(IN) :: lp, l, spin, js, isospin_tz
  real(dp), INTENT(INOUT) :: k_rel, q_rel
  INTEGER :: jxx, inn, n1, ix, iy
  complex*16  :: v,xmev,ymev
  !double precision  :: v,xmev,ymev
  CHARACTER (LEN=4) :: label
  COMMON /cnn/ inn
  COMMON /cpot/ v(6),xmev,ymev
  !common /cpts/   q(97),c,n1,ix,iy
  COMMON /cstate/ jxx, heform, sing, trip, coup, endep, label
  LOGICAL :: sing, trip, coup, heform, endep


  
  vv = dcmplx(0.d0,0.d0)
  heform=.FALSE.
  jxx=js
  sing= spin == 0; trip = spin == 1
  coup = .FALSE.
  IF ( (lp /= js).OR.(l /= js)) coup = .TRUE.
 ! write(6,*) l,lp,js,spin,trip, coup
  SELECT CASE ( isospin_tz)
  CASE (-1)
     inn = 2    !  pp case = 1
  CASE (0) 
     inn = 2    !  pn case = 2  if all inn=2, no CSB or ISB
  CASE ( 1)
     inn = 2    !  nn case = 3  if pp == nn, only ISB
  END SELECT
  SELECT CASE (type_of_pot)
  CASE('Idaho-A')
     inn = 1
  CASE('Idaho-B')
     inn = 2
  END SELECT
  xmev= phase*(k_rel)
  ymev = phase*(q_rel) 
  !xmev=abs(k_rel)
  !ymev=abs(q_rel) 
  
!  SELECT CASE (type_of_pot)
!  CASE('CD-bonn')
     CALL cmplx_cdbonn
!  CASE('Idaho-A')
!     CALL idaho
!  CASE('Idaho-B')
!     CALL idaho
!  CASE ( 'reid93')
!     CALL reid93
!  CASE ( 'argonnev18')
!     CALL argonp
!  CASE ( 'nij-I-II')
!     CALL nijm2
!  END SELECT
  IF (sing ) THEN
     vv=v(1)
     !write(6,*) vv
  ELSEIF ((trip).AND.(.NOT.coup )  ) THEN 
     vv = v(2)
  ELSEIF (coup ) THEN
     IF ( ( lp == js - spin).AND.(l == js - spin)) THEN 
        vv=v(4)
     ELSEIF ( ( lp == js + spin).AND.(l == js - spin)) THEN 
        vv=v(5)
     ELSEIF ( ( lp == js - spin).AND.(l == js + spin)) THEN 
        vv=v(6)
     ELSEIF ( ( lp == js + spin).AND.(l == js + spin)) THEN 
        vv= v(3)
     ENDIF
  ENDIF

END SUBROUTINE v_bare

!
! setup harmonic oscillator functions for 1st iteration
! 
SUBROUTINE HO_OSC_FUNC
  use hf_constants
  use constants
  use wave_functions
  implicit none
  integer ::  nh, lh, i,jh,tzh
  COMPLEX*16 :: zh, norm_sum
  double precision :: ph, qmin, qmax
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: q_points, weight_q
  
  DO nh=0,nh_max
     
     ph= (-1.D0)**nh
     DO lh=0,lh_max
        jh_min = abs(2*lh-1); jh_max = 2*lh+1
        do jh = jh_min, jh_max, 2
           do tzh = itzh_min, itzh_max, 2
              
              !norm_sum = 0.
              DO i=1,n_lab
                 
                 zh = phase * k_mesh(i) * oscillator_length
                 if ( nh == 0 ) then
                    !write(6,*) lh,jh,tzh
                    bhf_hol(i, lh, (jh+1)/2, (tzh+1)/2 ) =  ph*rnl(nh,lh,zh)*dSQRT(oscillator_length**3) 
                 end if
                 !norm_sum = norm_sum + weight_q(i)* q_points(i)**2*hf_holes(nh, lh, jh, tzh, i )**2 
                 !rnl(nh,lh,zh)**2*(oscillator_length**3) 
                 
              ENDDO
              !WRITE(6,*) 'Norm rel ho wf n,l : ', lh,jh,tzh,norm_sum,nh,lh,oscillator_length 
           end do
        ENDDO
     ENDDO
  end DO
  
  k_interpol_min = k_mesh(1)
  k_interpol_max = k_mesh(n_lab)
!  write(6,*) k_interpol_min, k_interpol_max

!  do i = 1, n_lab
!     zh = phase * q_points(i) * oscillator_length
!     write(6,*) dble(q_points(i)),bhf_hol(i,0,(1+1)/2,(-1+1)/2), rnl(0,0,zh)*dSQRT(oscillator_length**3) 
!  end do

END SUBROUTINE HO_OSC_FUNC

!
! klab = the interpolation point
! wavefunc = the interpolated functionat klab
! n,l,j,tz = quantum numbers of hole state
! 
subroutine interpol_1dim( n,l,j,tz,klab, wavefunc) 
  use constants
  USE wave_functions 
  use hf_constants
  
  implicit none
  integer, intent(in) :: n,l,j,tz
  complex*16 , intent(in) :: klab
  complex*16, intent(out) :: wavefunc
  !double precision, dimension(n_lab) :: x, y
  double precision, dimension(100) :: x, y
  double precision :: s, qmin, qmax, xval
  double precision, allocatable,dimension(:) :: q_points
  double precision, allocatable, dimension(:) :: lambda, C, wrk
  !double precision, dimension(n_lab) :: xx, yy
  double precision, dimension(100) :: xx, yy
  integer :: i, lck, lwrk, ifail
  
  !allocate( q_points(n_lab) )
  !qmin = 0.0_dp; qmax = hole_cutoff
  !q_points(:) = 0.5d0*( qmax+qmin ) + 0.5d0*( qmax-qmin )*x_mesh(:)
  
  !  write(6,*) klab, q_points
  if ( abs(klab) - k_ws(1) < 0. ) write(6,*) 'wrong in 1dim interpolation'
  if ( abs(klab) - k_ws(100) > 0. ) write(6,*) 'wrong in 1dim interpolation'
  
  x = k_ws
  
  do i = 1, 100
    
     !x(i) = 0.09*i
     !y(i) = exp(-x(i)) 
     y(i) = dble( ws_basis(l,j,i) )
     !y(i) = dble( hf_holes(n, l, j, tz, i ) )
     !if ( y(i) == 0. ) then
     !write(6,*) l,j, x(i), y(i),ws_basis(l,j,i) 
     !end if
     
  end do
  
    
  !write(6,*) x, klab
  xx = (x); yy = (y)
  !write(6,*) xx
  ! LCK > n_lab+4
  LCK = 100 + 4
  lwrk = 6*100 + 16 
  ifail = 0
  allocate( lambda(lck), c(lck), wrk(lwrk) )
  
  ! input: 
  ! constraint M > 4
  ! x(1) <...< x(M)  
  ! y(1)...y(M)
  ! LCK >= M + 4
  ! LWRK >= 6*M + 16
  ! E01BAF(M, x, y, lambda, C, LCK, WRK, LWRK, IFAIL)
  ! See documentation for further details. 
  call E01BAF(100, xx, yy, lambda, C, LCK, WRK, LWRK, IFAIL)
  
  xval = dble(klab)
  call E02BBF(100+4, dble(lambda), dble(C), xval, S, IFAIL)
  
  wavefunc = dcmplx(s)
  
  deallocate( lambda, c, wrk ) 

end subroutine interpol_1dim

!
!
!
subroutine interpol_1dim2( n,l,j,tz,klab, wavefunc) 
  use constants
  USE wave_functions 
  use hf_constants
  
  implicit none
  integer, intent(in) :: n,l,j,tz
  complex*16 , intent(in) :: klab
  complex*16, intent(out) :: wavefunc
  double precision, dimension(n_lab) :: x, y
  double precision :: s, qmin, qmax, xval
  double precision, allocatable, dimension(:) :: lambda, C, wrk
  double precision, dimension(n_lab) :: xx, yy
  integer :: i, lck, lwrk, ifail
  

  if ( abs(klab) < k_interpol_min  ) write(6,*) 'wrong in 1dim interpolation'
  if ( abs(klab) > k_interpol_max ) write(6,*) 'wrong in 1dim interpolation'
  
!  x = dble( k_mesh )
!  y = dble(  bhf_hol(:, l, (j+1)/2 ,(tz+1)/2 ) )
  !x = k_mesh
  do i = 1, n_lab
    
     !x(i) = 0.09*i
     !y(i) = exp(-x(i)) 
     x(i) = dble( k_mesh(i) )
     y(i) = dble( bhf_hol(i, l, (j+1)/2 ,(tz+1)/2 ) )
     !y(i) = dble( hf_holes(n, l, j, tz, i ) )
     !if ( y(i) == 0. ) then
     !write(6,*) l,j, x(i), y(i)
     !end if
     
  end do
  
  xx = (x); yy = (y)
  !write(6,*) xx
  ! LCK > n_lab+4
  LCK = n_lab + 4
  lwrk = 6*n_lab + 16 
  ifail = 0
  allocate( lambda(lck), c(lck), wrk(lwrk) )

  call E01BAF(n_lab, xx, yy, lambda, C, LCK, WRK, LWRK, IFAIL)
  
  xval = dble(klab)
  call E02BBF(n_lab+4, dble(lambda), dble(C), xval, S, IFAIL)
  
  wavefunc = dcmplx(s)
  
  deallocate( lambda, c, wrk ) 

end subroutine interpol_1dim2
!
!
! 
subroutine hartree_fock_final(n_osc, nlab_large)
  USE configurations
  USE constants
  use hf_constants
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n_osc, nlab_large
  INTEGER :: a, c, ij, ipar, itz, ket , bra
  INTEGER :: na,nc, la, lc, ja, jc, tza, tzc, iph
  REAL(DP) ::  vpot, pa, wpa, pc, wpc, kmin, kmax, fxy_re, pa_temp, pc_temp, ph
  double precision, allocatable, dimension(:) :: k_input, klab_large, wlab_large
  complex(dpc) :: diaghf, klab, norm_sum,zh
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:) :: hf_mtx, hf_vect
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:) :: hf_eigen, e_temp
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: vlab_interp, vlab_input
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:,:) :: ho_osc_func
  double precision, allocatable, dimension(:,:) :: vint
  integer :: i, j,  number, nstep, nsmall, la_temp, ja_temp, tza_temp
  integer :: ntot, i1,i2,l1,j1,tz1
  double precision, allocatable, dimension(:) ::X,Y,XX, WORK,W, D
  double precision :: val,vall
  integer:: ifail, IG1
  integer :: nh, lh, jh,tzh
  complex(dpc) :: kinetic, pot, cmplx_phase, cmplx_oscl, cmplx_theta

  open(10,file='hf6_pot12_12.dat') 
  ! count number of data lines in woods-saxon basis
  number=0
3 READ(10,*,end=4)
  number=number+1
  GOTO 3
4 CONTINUE
  REWIND 10
  write(6,*) 'qbox number' , number

  nsmall = 0 
  jh_min = 1; jh_max = 1
  lh_min = 0; lh_max = 0
  do i = 1, number
     read(10,*) bra, la, ja, tza,pa, ket, lc, jc, tzc,pc, vpot
     if ( bra > nsmall ) nsmall = bra
     if ( la > lh_max ) lh_max = la
     if ( ja > jh_max ) jh_max = ja
  end do
  rewind 10
  
  IG1 = nsmall
  ALLOCATE ( vlab_input(nsmall,nsmall, 0:lh_max, (jh_min+1)/2:(jh_max+1)/2 ,(-1+1)/2:(1+1)/2) )
  ALLOCATE ( k_input(nsmall) ) 
  allocate( x(nsmall), y(nsmall) )
  allocate( XX(nsmall), WORK(nsmall),W(nsmall), D(nsmall) )
  allocate( vint(nsmall, nsmall) )

  write(6,*) jh_min, jh_max,lh_min, lh_max
  !write(6,*) 'min/max s.p. quantum number: j,l,itz', jh_min, jh_max, lh_min,lh_max, itzh_min,itzh_max
  do i = 1, number
     read(10,*) bra, la, ja, tza,pa, ket, lc, jc, tzc,pc, vpot
     i1 = bra; i2 = ket; l1 =la;j1=ja;tz1=tza
     vlab_input(i1,i2,l1,(j1+1)/2,(tz1+1)/2) = dcmplx( vpot )
     if ( bra /= ket ) then
        vlab_input(i2,i1,l1,(j1+1)/2,(tz1+1)/2) = vlab_input(i1,i2,l1,(j1+1)/2,(tz1+1)/2) 
     end if
     k_input(bra) = pa
  end do
  kmin = k_input(1)
  kmax = k_input(nsmall)
  
  
  !nlab_large = 100
  allocate ( klab_large(nlab_large), wlab_large(nlab_large) )
  ALLOCATE ( vlab_interp(nlab_large, nlab_large, 0:lh_max, & 
       (jh_min+1)/2:(jh_max+1)/2 ,(-1+1)/2:(1+1)/2) )
 
  
  call gauss_legendre(kmin,kmax,klab_large,wlab_large,nlab_large )
  
  nh_max = n_osc
  ALLOCATE ( ho_osc_func(nlab_large,0:lh_max, 0:nh_max ))
  
  !
  ! construct complex harmonic oscillator basis
  ! see article by Gyarmati and Kruppa: phys.rev.c volume34,page 34 (1986)
  !
  cmplx_theta = dcmplx(0.d0,0.1d0)
  cmplx_phase = exp( cmplx_theta )  
  ! set oscillator length/energy 
  oscillator_energy = 7.
  oscillator_length = 1./SQRT(p_mass(0)*oscillator_energy)
  cmplx_oscl = oscillator_length * cmplx_phase
  DO nh=0,nh_max
     ph= (-1.D0)**nh
     DO lh=0,lh_max
        norm_sum = 0.
        DO i=1,nlab_large
           
           
           zh =  klab_large(i) * cmplx_oscl
           ho_osc_func(i, lh, nh ) =  ph*rnl(nh,lh,zh)* SQRT(cmplx_oscl**3) 
           norm_sum = norm_sum + wlab_large(i)* klab_large(i)**2*ho_osc_func(i, lh, nh)**2
                    
        ENDDO
        WRITE(6,*) 'Norm rel ho wf n,l : ', nh, lh,norm_sum 
     end do
  ENDDO
  
  
  !
  ! set up hartree-fock potential on large grid using cubic spline interpolation 
  ! 
  ! loop over single-particle isospin
  DO tza = -1, 1, 2
     ! loop over single-particle orbital angular momentum
     do la = lh_min, lh_max, 1 
        !     loop over single-particle angular momentum (twice the value)
        DO ja = jh_min, jh_max, 2 
           
           if ( ja > ( 2*la + 1) ) cycle
           if ( ja < abs( 2*la - 1 ) ) cycle
           
           vint = 0.
           vint(:,:) = real( vlab_input(:,:,la,(ja+1)/2,(tza+1)/2 ) )
           DO bra = 1, nlab_large
              pa = klab_large(bra); wpa = wlab_large(bra)
              DO ket= bra, nlab_large
                 pc = klab_large(ket); wpc = wlab_large(ket)

                 call E01ACF((pa), (pc), (k_input), (k_input), vint, & 
                      VAL, VALL, IFAIL, XX, WORK,W, D, IG1, nsmall,nsmall)

                 vlab_interp(bra,ket,la,(ja+1)/2,(tza+1)/2) = dcmplx( vall) !dcmplx(fxy_re)
                 
                 if ( bra /= ket ) then
                    vlab_interp(ket,bra,la,(ja+1)/2,(tza+1)/2) = vlab_interp(bra,ket,la,(ja+1)/2,(tza+1)/2)
                 end if
                 
                 
              end DO
           end DO
        end DO
     end do
  end DO

  !
  ! diagonalize hartree-fock potential using plane-wave basis on large grid
  ! 
  write(6,*) 'hartree-fock spectrum calculated using large plane-wave basis:'
  ! loop over single-particle isospin
  DO tza = -1, 1, 2
     ! loop over single-particle orbital angular momentum
     do la = lh_min, lh_max, 1 
        !     loop over single-particle angular momentum (twice the value)
        DO ja = jh_min, jh_max, 2 
           
           if ( ja > ( 2*la + 1) ) cycle
           if ( ja < abs( 2*la - 1 ) ) cycle
           
           ALLOCATE(hf_mtx(nlab_large,nlab_large) )
           ALLOCATE(hf_vect(nlab_large,nlab_large) )
           ALLOCATE(hf_eigen(nlab_large) )
                     
           !     set up hf matrix to be diagonalized
           hf_mtx = 0.0_dp; hf_vect = 0.0_dp; hf_eigen = 0.0_dp
           DO bra = 1, nlab_large
              pa = klab_large(bra); wpa = wlab_large(bra)
              DO ket= bra, nlab_large
                 pc = klab_large(ket); wpc = wlab_large(ket)
                 
                 
                 diaghf = dble( vlab_interp(bra,ket,la,(ja+1)/2,(tza+1)/2) )
                 !v_bhf( bra, ket, la,(ja+1)/2 ,(tza+1)/2 )
                 ! Interaction in units of MeV and preparing for diagonalization
                 hf_mtx(bra,ket) = phase**3 * pa*pc*diaghf*SQRT(wpa*wpc)
                 hf_mtx(ket,bra)=hf_mtx(bra,ket)
                 ! Add kinetic energy
                 IF ( bra == ket) THEN
                    hf_mtx(bra,bra)=hf_mtx(bra,bra) + phase**2*pa*pa*0.5_dp/p_mass(tza)
                 ENDIF
              ENDDO
           ENDDO
           ! diagonalize and get new eigenvalues and energies
           CALL diag_exact( hf_mtx, hf_vect, hf_eigen, nlab_large)
           
           
           do i = 1, nlab_large
              
              if ( dble( hf_eigen(i) )< 0.d0 ) then
                 write(6,*) 'sp quantum numbers:tza,ja,la and energy',tza, ja, la,real( hf_eigen(i) )  
              end if
           end do

           DEALLOCATE ( hf_mtx); DEALLOCATE ( hf_eigen)
           DEALLOCATE ( hf_vect)
           
        ENDDO
     ENDDO
  ENDDO
  

  !
  ! diagonalize hartree-fock potential using plane-wave basis on large grid
  ! 
  write(6,*) 'hartree-fock spectrum calculated using large harmonic oscillator basis:'
  ! loop over single-particle isospin
  DO tza = -1, 1, 2
     ! loop over single-particle orbital angular momentum
     do la = lh_min, lh_max, 1 
        lc = la
        !     loop over single-particle angular momentum (twice the value)
        DO ja = jh_min, jh_max, 2 
           
           if ( ja > ( 2*la + 1) ) cycle
           if ( ja < abs( 2*la - 1 ) ) cycle
           
           ALLOCATE(hf_mtx(n_osc+1,n_osc+1) )
           ALLOCATE(hf_vect(n_osc+1,n_osc+1) )
           ALLOCATE(hf_eigen(n_osc+1) )
           
           !     set up hf matrix to be diagonalized
           hf_mtx = 0.0_dp; hf_vect = 0.0_dp; hf_eigen = 0.0_dp
           DO na = 0, n_osc
              bra = na + 1
              
              DO nc= na, n_osc
                 ket = nc + 1
                 
                 ! kinetic energy 
                 if ( na == nc ) then
                    kinetic = 0.
                    do i = 1, nlab_large 
                       pa = klab_large(i); wpa = wlab_large(i)
                       kinetic = kinetic + wpa * pa**2 * & 
                            ho_osc_func(i,la,na)*ho_osc_func(i,lc,nc)*pa*pa*0.5_dp/p_mass(tza)
                    end do
                    hf_mtx(bra,ket) = kinetic
                 elseif ( na == nc - 1 ) then
                    kinetic = 0.
                    do i = 1, nlab_large 
                       pa = klab_large(i); wpa = wlab_large(i)
                       kinetic = kinetic + wpa * pa**2 * & 
                            ho_osc_func(i,la,na)*ho_osc_func(i,lc,nc)*pa*pa*0.5_dp/p_mass(tza)
                    end do
                    hf_mtx(bra,ket) = kinetic
                 elseif ( na == nc + 1 ) then
                    kinetic = 0.
                    do i = 1, nlab_large
                       pa = klab_large(i); wpa = wlab_large(i)
                       kinetic = kinetic + wpa * pa**2 * & 
                            ho_osc_func(i,la,na)*ho_osc_func(i,lc,nc)*pa*pa*0.5_dp/p_mass(tza)
                    end do
                    hf_mtx(bra,ket) = kinetic
                 end if
                 
                 pot = 0.
                 do i = 1, nlab_large
                    pa = klab_large(i); wpa = wlab_large(i)
                    do j = 1, nlab_large
                       pc = klab_large(j); wpc = wlab_large(j)
                       pot = pot + wpa*wpc*pa**2*pc**2* ho_osc_func(i,la,na)*ho_osc_func(j,lc,nc)*& 
                            vlab_interp(i,j,la,(ja+1)/2,(tza+1)/2 )
                    end do
                 end do
                 
                 ! Interaction in units of MeV and preparing for diagonalization
                 hf_mtx(bra,ket) =hf_mtx(bra,ket)+ pot
                 
                 if ( bra /= ket ) then
                    hf_mtx(ket,bra)=hf_mtx(bra,ket)
                 end if

              ENDDO
           ENDDO
           ! diagonalize and get new eigenvalues and energies
           CALL diag_exact( hf_mtx, hf_vect, hf_eigen, n_osc+1)
           
           ! check for bound states 
           do i = 1, n_osc
              if ( dble( hf_eigen(i) )< 0.d0 .and. abs( aimag( hf_eigen(i)) ) < 0.5 ) then
                 write(6,*) 'sp quantum numbers:tza,ja,la and energy',tza, ja, la,real( hf_eigen(i) ), aimag(hf_eigen(i))  
              end if
           end do
           

!           do i = 1, n_osc
!              !if ( dble( hf_eigen(i) )< 0.d0 .and. abs( aimag( hf_eigen(i)) ) < 0.5 ) then
!              if ( ja == 5 .and. la == 2 ) then
!                 !write(6,*) 'sp quantum numbers:tza,ja,la and energy',tza, ja, la,real( hf_eigen(i) ), cmplx(hf_eigen(i))  
!                 write(6,*) real( hf_eigen(i) ), aimag(hf_eigen(i))  
!              end if
!           end do
           
                   
           DEALLOCATE ( hf_mtx); DEALLOCATE ( hf_eigen)
           DEALLOCATE ( hf_vect)
           
        ENDDO
     ENDDO
  ENDDO
  
  
END SUBROUTINE hartree_fock_final
!
! routine for counting total number of matrix elements used in
! self consistent evaluation of hartree-fock potential
!
SUBROUTINE count_number_mtxel(number_mtxel)
  USE single_particle_orbits
  USE configurations
  USE constants
  use hf_constants
  USE wave_functions
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: sp_configs
  INTEGER :: a, c, ij, ipar, itz, ket , bra
  INTEGER :: la, lc, ja, jc, tza, tzc, iph
  REAL(DP) ::  pa, pc, pb,pd, qmin_left, qmax_left, qmin_right, qmax_right, ptest,ph
  real(dp), ALLOCATABLE, DIMENSION(:) :: q_points_right, q_points_left
  integer :: i,  h,jh,lh,tzh,nh,jph,j_min, j_max,iq1,iq2
  integer, intent(out) :: number_mtxel
  
  k_interpol_min = k_mesh(1)
  k_interpol_max = k_mesh(n_lab) 

  ALLOCATE ( q_points_right(int_points), q_points_left(int_points) )
  
  number_mtxel = 0
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
                 IF( ABS(pa-pc) >  2*hole_cutoff) CYCLE
                 
                 ! bra side
                 qmin_left = 0.0_dp; qmax_left = hole_cutoff
                 IF (ABS(pa-pc) - hole_cutoff > 0.) qmin_left = (ABS(pa-pc) - hole_cutoff)
                 
                 
                 ! check that integration points are within the klab  limits 
                 ! for the bhf wave function -> no extrapolation!
                 ! if not -> new integration limits.
                 if ( qmin_left < k_interpol_min ) qmin_left = k_interpol_min
                 if ( qmax_left > k_interpol_max ) qmax_left = k_interpol_max
                 
                 q_points_left(:) = 0.5d0*( qmax_left+qmin_left ) + & 
                      0.5d0*( qmax_left-qmin_left )*x_mesh(:)
                 
                 ptest = all_orbit%p(1)
                 DO h=1, all_orbit%total_orbits
                    IF (all_orbit%orbit_status(h) /= 'hole') CYCLE
                    lh = all_orbit%ll(h) ; jh = all_orbit%jj(h)
                    nh = all_orbit%nn(h)  ; tzh = all_orbit%itzp(h)
                    IF (all_orbit%p(h) /= ptest) CYCLE   ! only one hole state, independent of p
                    j_min=ABS((ja-all_orbit%jj(h))/2)
                    j_max=(ja+all_orbit%jj(h))/2
                    IF ( tza+tzh /= tzc+tzh) CYCLE
                    DO jph=j_min,j_max
                      
                       ! now integrate over the HO wf for the  ket side
                       DO iq1=1, int_points
                          ph = q_points_left(iq1)
                          IF (ABS(pc-ph)-pa > hole_cutoff ) CYCLE 
                          IF (pa-ph-pc > hole_cutoff ) CYCLE 
                          pd = ph
                          ! leftside transformation
                          qmin_right= 0.0_dp; 
                          IF (ABS(pc-pd) - pa > 0.0_dp ) qmin_right= ABS(pc-pd) - pa 
                          IF (pa-pc-pd > qmin_right ) qmin_right= pa-pc-pd
                          qmax_right = hole_cutoff
                          IF (pa+pc+pd < qmax_right) qmax_right = pa+pc+pd
  
                          if ( qmin_right < k_interpol_min ) qmin_right = k_interpol_min
                          if ( qmax_right > k_interpol_max ) qmax_right = k_interpol_max
                          
                          q_points_right(:) = 0.5d0*( qmax_right+qmin_right ) + & 
                               0.5d0*( qmax_right-qmin_right )*x_mesh(:)
                          
                          ! now integrate over the HO for the  bra side
                          DO iq2=1, int_points
                             pb = q_points_right(iq2)
                             IF ( (ABS(pa-pb) > pc+pd) .OR. (ABS(pc-pd) > pa+pb) ) CYCLE
                             
                             number_mtxel = number_mtxel +1 
     
                          END DO
                       END DO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  DEALLOCATE ( q_points_right, q_points_left)
  write(6,*) 'number of matrix elements',number_mtxel
  
END SUBROUTINE count_number_mtxel
