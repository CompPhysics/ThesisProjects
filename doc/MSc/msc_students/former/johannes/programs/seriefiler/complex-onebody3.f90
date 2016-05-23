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
  INTEGER :: la, lc, ja, jc, tza, tzc, iph, r1, iq
  REAL(DP) ::  pa, pc, kk, check_sign
  real(8) :: re_e, im_e
  complex(dpc) :: diaghf,norm_check, int_sum, klab, kr, bessel_p1, norm_sum,zh, energy_check
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:) :: hf_mtx, hf_vect
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:) :: hf_eigen, e_temp, temp
  complex(dpc), dimension(n_lab) :: z_mesh, wz_mesh
  integer :: i, no_zeros, number_mtxel
  integer, intent(in) :: iteration_number
  logical, intent(in) :: selfconsistency
  double precision, intent(out) :: sigma
  integer :: n,l,NN,LL, n1,l1,n2,l2

  !n = 1; l = 1; NN = 0; LL = 0
  !lambda = 1
  !n1 = 1; n2 = 0; l1 = 1; l2= 0
  !call check_vector_bras( n,l,NN,LL, n1,l1,n2,l2, lambda )
    
  !return
  
  ! complex mesh setup
  z_mesh = phase*k_mesh; wz_mesh = phase*w_mesh 

  !
  ! for hf_iteration > 1 no longer woods-saxon basis
  ! setup new limits for interpolation routine. 
  ! 
  !if ( iteration_number > 1 ) then
  !k_interpol_min = k_mesh(1)
  !k_interpol_max = k_mesh(n_lab) 
  !end if
  
  !
  ! here we calculate hf-potentials self consistently.
  ! Note: this is done for real momenta only!
  !
  sigma = 0.0_dp
  v_bhf = 0.d0
  number_mtxel = 0
  no_zeros = 0
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
           ALLOCATE (e_temp(sp_configs%number_confs), temp(sp_configs%number_confs) )
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
                 
                 
                 !if ( bra /= ket ) cycle
                 
                 CALL hf(diaghf,pa,la,ja,tza,pc,lc,jc,tzc, iteration_number, number_mtxel, no_zeros)
                 
                 v_bhf( bra, ket, la,(ja+1)/2 ,(tza+1)/2 ) = ( diaghf ) 
                 if (  selfconsistency ) then 
                    WRITE(10,*) bra, la, ja, tza,pa, ket, lc, jc, tzc,pc, & 
                         real(diaghf), aimag(diaghf)
                 END if
                 
                 
                 ! Interaction in units of MeV and preparing for diagonalization
                 hf_mtx(bra,ket) = diaghf*SQRT(all_orbit%w_p(a)*all_orbit%w_p(c))*pa*pc  
                 hf_mtx(ket,bra)=hf_mtx(bra,ket)
                 ! Add kinetic energy
                 
                 IF ( bra == ket) THEN
                    write(6,*) la,ja,tza, pa, dble(diaghf),aimag(diaghf)
                    hf_mtx(bra,bra)=hf_mtx(bra,bra) + pa*pa*0.5_dp/p_mass(tza)*(1.d0-1.d0/mass_closed_shell_core)
                 ENDIF
              ENDDO
           ENDDO
           ! diagonalize and get new eigenvalues and energies
           CALL diag_exact( hf_mtx, hf_vect, hf_eigen, sp_configs%number_confs)
                      
           
           !    set up the new bhf wave function
           bhf_hol_temp(:,la,(ij+1)/2,(tza+1)/2) = hf_vect(:,1)
           
           
           
           temp = matmul( hf_mtx, bhf_hol_temp(:,la,(ij+1)/2,(tza+1)/2) )
           energy_check = 0.
           do i = 1, n_lab
              energy_check = energy_check + bhf_hol_temp(i,la,(ij+1)/2,(tza+1)/2) * temp(i)
           end do
           
           !    re-normalize bhf wave function; sum( bhf**2*k_mesh**2*w_mesh = 1 ) 
           bhf_hol_temp(:,la,(ij+1)/2,(tza+1)/2) = bhf_hol_temp(:,la,(ij+1)/2,(tza+1)/2)/z_mesh(:)/sqrt(wz_mesh(:))
           
           !
           ! calculate r-space functions
           !
!          do r1 = 1, nrad
!             
!             int_sum = 0.0
!             do iq = 1, n_lab
!                
!                kr = r_lab(i)*k_mesh(iq)/hbarc
!                call CMPLX_BESSELJ(la,kr,BESSEL_p1)
!                   
!                int_sum = int_sum + sqrt(2./pi)*k_mesh(iq)**2*w_mesh(iq)* & 
!                     bhf_hol_temp(iq,la,(ij+1)/2,(tza+1)/2)*dble( bessel_p1 )*sqrt(1./hbarc**3)
!                
!             end do
!             bhf_hol_rspace_temp(r1,la,(ij+1)/2,(tza+1)/2) = int_sum
!          end do
           

           check_sign =  sum(  dble(bhf_hol_temp(:,la,(ij+1)/2,(tza+1)/2))  )
           if ( check_sign < 0. ) then
              write(6,*) 'check sign', check_sign
              bhf_hol_temp(:,la,(ij+1)/2,(tza+1)/2) = -bhf_hol_temp(:,la,(ij+1)/2,(tza+1)/2)
           end if
           
           write(6,*) 'norm diff', sum(  z_mesh(:)**2*wz_mesh(:)* & 
                bhf_hol(:,la,(ij+1)/2,(tza+1)/2)*bhf_hol_temp(:,la,(ij+1)/2,(tza+1)/2)  ) 
           do i = 1, n_lab
          !    write(6,*) abs(z_mesh(i)), dble(bhf_hol_temp(i,la,(ij+1)/2,(tza+1)/2)),  &
          !         aimag( bhf_hol_temp(i,la,(ij+1)/2,(tza+1)/2)  ) , & 
          !         dble(bhf_hol(i,la,(ij+1)/2,(tza+1)/2) ), aimag( bhf_hol(i,la,(ij+1)/2,(tza+1)/2) )
          !    
              write(6,*) dble(hf_eigen(i)), aimag(hf_eigen(i))
           end do

           
           norm_check   = 0.d0 

           do i = 1, n_lab
              norm_check = norm_check + z_mesh(i)**2*wz_mesh(i)* & 
                   bhf_hol_temp(i,la,(ij+1)/2,(tza+1)/2) **2 
           end do
           re_e = real(energy_check)
           im_e = aimag(energy_check)
           if ( (.not.  selfconsistency ) ) then 
              write(6,*) 'sp quantum numbers:',itz, ij,la 
              write(6,*) 'check energy and norm of bhf_hol for this iteration:', iteration_number, re_e,im_e,&
abs(norm_check) 
           END if


           if (  selfconsistency ) then 
              write(6,*) 'sp quantum numbers:',itz, ij,la 
              write(6,*) 'check energy and norm of bhf_hol:', energy_check, norm_check 
              DO bra = 1, n_lab
                 write(6,*) DBLE(hf_eigen(bra)), AIMAG(hf_eigen(bra))
                 !write(6,*) 'Single-particle energy:',hf_eigen(bra)
              ENDDO
           end if
           
           ! calculate variance of bhf wave function
           DO bra = 1, n_lab
              sigma = sigma + ABS( z_mesh(bra)*sqrt(wz_mesh(bra))* & 
                   (bhf_hol_temp(bra,la,(ij+1)/2,(tza+1)/2) - bhf_hol(bra,la,(ij+1)/2,(tza+1)/2) ) )
           ENDDO
           
           
           DEALLOCATE ( hf_mtx); DEALLOCATE ( hf_eigen)
           DEALLOCATE ( hf_vect); DEALLOCATE ( e_temp, temp)
           
        ENDDO
     ENDDO
  ENDDO
  
  ! allocate space for momentum space interaction
  !if ( iteration_number == 1 ) allocate( vlab_mom(number_mtxel) )
  
  write(6,*) 'number of matrix elements in this iteration',number_mtxel, no_zeros, no_zeros/number_mtxel, hf_iterations
  

  ! deallocate matrix for hole functions
  !call deallocate_hf_holes


  !1000 FORMAT(//2X,4H--->,2X,36H Momenta in MeV HF diagram in MeV^-2,2X,4H<---//)
  !1001 FORMAT(2X,12HA: l, j, tz=,3I3,2X,3HPA:,E12.6,2X,12HC: l, j, tz=,3I3,2X,3HPC:,E12.6,2X,3HHF=,E12.6)

END SUBROUTINE hartree_fock

!
!                  SUBROUTINE to calculate the Hartree-Fock contrib.    
!                  Dimensionality = MeV^-2, as the interaction      
!

SUBROUTINE hf(hfe,pa,la,ja,tza,pc,lc,jc,tzc, iteration_number, number_mtxel, no_zeros)
  USE constants
  use hf_constants
  USE single_particle_orbits
  USE wave_functions  
  IMPLICIT NONE
  complex(DPC), INTENT(INOUT) :: hfe
  REAL(DP), INTENT(IN) :: pa, pc
  INTEGER, INTENT(IN) :: ja, jc, la, lc, tza, tzc, iteration_number
  integer, intent(out) :: number_mtxel, no_zeros
  INTEGER :: h,nh, lh, jh, tzh, j_min, j_max, jph, itz, iq
  REAL(DP) :: qmin, qmax,  ph, ptest
  complex(dpc) :: zh, kr, bessel_p1
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: q_points, weight_q, temp
  complex(dpc) :: vmtx,sum_hf, hole_func, int_sum
  integer ::  no_r, i

  hfe=0.0_dp
  qmin = 0.0_dp
  qmax = 2.d0*Lambda - pc 
  
  
  !IF (ABS(pa-pc) - hole_cutoff > 0.) qmin = (ABS(pa-pc) - hole_cutoff)
  ALLOCATE ( q_points(int_points), weight_q(int_points), temp(int_points) )
    
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
        !temp = 0.
        ! now integrate over the HO wf for the  ket side
        DO iq=1, int_points
           
           ph = q_points(iq)
           
           !IF (ABS(pc-ph)-pa > hole_cutoff ) CYCLE 
           !IF (pa-ph-pc > hole_cutoff ) CYCLE 
           !zh = phase * ph*oscillator_length
           
           
           ! in first iteration use woods-saxon basis 
           ! after that we interpolate the new bhf wave function on the integration points
           !if ( iteration_number == 1 ) then
           !   call interpol_1dim( nh,lh,jh,tzh, phase*ph, hole_func)
           !elseif (iteration_number > 1 ) then
           
           call interpol_1dim2( nh,lh,jh,tzh, phase*ph, hole_func)
           
           !  Get <pa la ja nh lh jh| v |pc lc jc ph lh jh> 
           CALL leftside_transf(vmtx,pa,la,ja,nh,lh,jh,tzh,pc,lc,jc,ph,lh,jh,itz,jph,& 
                iteration_number, number_mtxel, no_zeros)
           !sum_hf=sum_hf+phase**3*weight_q(iq)*rnl(nh,lh,zh) &
           !     *ph*SQRT(oscillator_length**3)*vmtx
           
           !temp(iq) = weight_q(iq)*ph*ph*vmtx
           
           sum_hf=sum_hf + weight_q(iq)*ph*ph*hole_func*vmtx
           
           !sum_hf=sum_hf+phase**2*weight_q(iq)*ph*ph*hole_func*vmtx

           
        END DO
        
        
        ! 
        ! radial integration
        !
   !     do i = 1, nrad
   !        int_sum = 0.
   !                   
   !        hole_func = bhf_hol_rspace(i,lh,(jh+1)/2,(tzh+1)/2) 
   !        
   !        !write(6,*) nh, lh, jh, tzh,jph, r_lab(i), (hole_func)
   !        do iq = 1, int_points
   !           ph = q_points(iq)
   !           kr = r_lab(i)*ph/hbarc
   !           call CMPLX_BESSELJ(lh,kr,BESSEL_p1)
   !                 
   !           int_sum = int_sum + sqrt(2./pi)*temp(iq)*dble( bessel_p1 )*sqrt(1./hbarc**3)
   !
!           end do
 !          sum_hf=sum_hf + wr_lab(i)*r_lab(i)**2*hole_func*int_sum
           
  !         !write(6,*) pa,pc, r_lab(i), int_sum, hole_func
  !      end do



        hfe=hfe+sum_hf*(2*jph+1.0_dp)
     ENDDO
  ENDDO
  
  hfe=hfe/(ja+1.0_dp)

  
 

  DEALLOCATE ( weight_q, q_points, temp )

END SUBROUTINE hf

!
!     Calculates the matrix elements
!     <pa la ja pb lb jb | V | pc lc jc pd ld jd  >
!     using as input the matrix elements
!     <pa la ja nb lb jb | V | pc lc jc pd ld jd  >
!     Dimensionality = MeV^-2, as the interaction      
!

SUBROUTINE leftside_transf(sum_v,pa,la,ja,nb,lb,jb,tzb,pc,lc,jc,pd,ld,jd,itz,jph, iteration_number, number_mtxel,& 
no_zeros)
  USE constants
  use hf_constants
  USE single_particle_orbits
  USE wave_functions
  
  IMPLICIT NONE
  complex(DPC), INTENT(INOUT) :: sum_v
  REAL(DP), INTENT(IN) :: pa, pc, pd
  INTEGER, INTENT(IN) :: ja, jc, la, lc, nb, lb, jb, tzb,ld, jd, jph, itz, iteration_number
  integer, intent(out) :: number_mtxel, no_zeros
  INTEGER :: iq, no_r,i
  REAL(DP) :: qmin, qmax, pb
  complex(dpc) :: vmtx, zb, hole_func, kr, bessel_p1, int_sum
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: q_points, weight_q, temp
  
  ! Finding integration limits
  
  !IF ( (ABS(pa-pb) > pc+pd) .OR. (ABS(pc-pd) > pa+pb) ) CYCLE
  ! is taken care of here:
  qmin = max( pa-pc-pd, ABS(pc-pd) - pa, k_interpol_min ) 
  qmax = min( 2.d0*lambda-pa, pa+pc+pd, k_interpol_max )
  
  !qmin = max( pa-pc-pd, ABS(pc-pd) - pa )
  !qmax = min( 2.d0*lambda-pa, pa+pc+pd )
  
  if ( qmax < qmin ) write(6,*) qmin, qmax


  ALLOCATE ( q_points(int_points), weight_q(int_points), temp(int_points) )
    
  q_points(:) = 0.5d0*( qmax+qmin ) + 0.5d0*( qmax-qmin )*x_mesh(:)
  weight_q(:) = 0.5d0* ( qmax-qmin ) * wx_mesh(:)
 
  !CALL gauss_legendre(qmin,qmax,q_points,weight_q,int_points)
  sum_v = 0.0_dp
  ! now integrate over the HO for the  bra side
  DO iq=1, int_points
     pb = q_points(iq)
     
     zb = phase * pb*oscillator_length

     
     call interpol_1dim2( nb,lb,jb,tzb, phase*pb, hole_func)
          
     number_mtxel = number_mtxel +1 
     
     !  Get <pa la ja pb lb jb| v |pc lc jc pd ld jd> 
     if ( iteration_number == 1 ) then
        CALL v_labmomentum_space(vmtx,la,ja,pa,lb,jb,pb,lc,jc,pc,ld,jd,pd,itz,jph)
        vlab_mom_re(number_mtxel) = real(vmtx)
        !vlab_mom_im(number_mtxel) = real(aimag(vmtx))
     elseif ( iteration_number > 1 ) then
        vmtx =  dcmplx(vlab_mom_re(number_mtxel), 0.d0 )
        
     end if
     
     if ( abs(vmtx) == 0.d0 ) no_zeros = no_zeros + 1
     
     !temp(iq) = weight_q(iq)*pb*pb*vmtx
     !sum_v=sum_v+phase**3*weight_q(iq)*rnl(nb,lb,zb)*SQRT(oscillator_length**3)*pb*vmtx
     sum_v=sum_v + weight_q(iq)*hole_func*pb*pb*vmtx
     
  END DO
  
  
  ! 
  ! radial integration
  !
!  do i = 1, nrad
!     int_sum = 0.
!     hole_func = bhf_hol_rspace(i,lb,(jb+1)/2,(tzb+1)/2) 
!     do iq = 1, int_points
!        pb = q_points(iq)
!        kr = r_lab(i)*pb/hbarc
!        call CMPLX_BESSELJ(lb,kr,BESSEL_p1)
!                    
!        int_sum = int_sum + sqrt(2./pi)*temp(iq)*dble( bessel_p1 )
!        
!     end do
!     sum_v=sum_v + wr_lab(i)*r_lab(i)**2*hole_func*int_sum*sqrt(1./hbarc**3)!

     !write(6,*) pa,pc,pd, r_lab(i), int_sum
!  end do


  

  DEALLOCATE ( weight_q, q_points, temp )

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
  !IF(p_ab >=  p_cd) 
  qmin = 0.5D0*q_cd
  IF(p_cd > p_ab) qmin = SQRT(kinetic_ket-.25D0*p_ab*p_ab)
  !IF(q_cd >= q_ab) 
  qmax = 0.5D0*p_cd
  IF(q_ab > q_cd) qmax = SQRT(kinetic_ket-0.25D0*q_ab*q_ab)
  ! adjust max relative momentum (qmax) for renormalized interaction
  ! since V_lowk(k,k') = 0 for k,k' > kcut = abs(krel(n_k1)) 
  qmax = min( qmax, abs(krel(n_k1) ) )
!  SELECT CASE(type_of_interaction)
!  CASE('free')
!     qmax = qmax
!  CASE('renorm')
!     if ( qmax > abs( krel(n_k1) ) ) then
!        qmax = abs( krel(n_k1) )
!     else
!        qmax = qmax 
!     end if
!  END SELECT

  !q_points(:) = 0.5d0*( qmax+qmin ) + 0.5d0*( qmax-qmin )*x_mesh(:)
  !weight_q(:) = 0.5d0* ( qmax-qmin ) * wx_mesh(:)
  !Under gir ikke lenger NaN 
  q_points(:) = (0.5d0*(qmax-qmin)*(t_mesh(:)+2.d0*qmax/(qmax-qmin) -1.d0))
  weight_q(:) = (0.5d0*(qmax-qmin)*wt_mesh(:))* dsqrt( 1.d0 - t_mesh(:)*t_mesh(:) )
 ! write(*,*) 'hei, en test for debugging' !HER ER DET EN SEGMENT FEIL
  !q_rel = phase* dcmplx(0.5d0*(qmax-qmin)*(t(iq)+2.d0*qmax/(qmax-qmin) -1.d0))
  !dq_rel = phase*dcmplx(0.5d0*(qmax-qmin)*w_t(iq))
  !cheb_weight = dsqrt( 1.d0 - t(iq)*t(iq) )
                    
  !CALL gauss_legendre(qmin,qmax,q_points,weight_q,int_points)
  ! Start performing the transformation for CoM and relative coordinates
  ! to lab frame, ket side first
  is1 = 1; is2=1
  DO spin= 0, 1  ! S = 0 or 1 for nucleons only
     DO js = jmin, jmax
        
        DO l = abs(js-spin), js+spin
           !DO l=0, lmax  ! relative orbital momentum, (jmax changed to lmax)
           !jsmin=ABS(l-spin); jsmax=l+spin
        
           IF ( (ABS(isospin_tz) == 1).AND.(MOD(l+spin+ABS(isospin_tz),2)==0)) CYCLE
           !   DO js=jsmin,jsmax    ! Relative angular momentum
           !      IF(js > jmax) CYCLE 
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
               !  write(*,*) lsjj_coeff_ket
                 ! The integral over momenta is the last loop
                 DO iq=1, int_points
                    
                    q_rel = q_points(iq)
                    pos = 4.D0*(kinetic_ket-q_rel*q_rel)
                    IF(pos <=  0.D0) CYCLE 
                    p_com = SQrt(pos)
                    IF(p_com  > p_ab) CYCLE; IF(p_com  < q_ab) CYCLE
                    IF(p_com  > p_cd) CYCLE; IF(p_com  < q_cd) CYCLE
                    pos  = kinetic_bra - 0.25D0*p_com*p_com
                    IF ( pos < 0.D0 ) CYCLE; k_rel = SQRT(pos)  
                  !  write(*,*) kinetic_bra, p_com, pc+pd, DABS(pc-pd) 
                    vb_ket = vector_trcoefficients(q_rel,p_com,pc,pd,l,ll,&
                             lambda_ket,lc,ld,isospin_tz)/(pc*pd)
                  !  write(*,*) 'vb_ket', vb_ket
                   ! TIL HIT ER DET BRA, INGEN UENDELIGHETER I vb_ket 
                    IF (ABS(vb_ket) < epsilon) CYCLE !ok
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
                       
                       
                       !vv_cmplx = 0.0
                       !if ( l > 0 ) cycle
                       !if ( lp /= l ) cycle
                       !                       if ( channel /= 1 ) cycle
!                       if ( js > 3 ) cycle
                       
                       !call malfliet_tjon( vv_cmplx, l, phase*k_rel, phase*q_rel ) 
                       !call gaussian( vv_cmplx, l, js, phase*k_rel, phase*q_rel ) 
 
!                       if ( l == 0 .and. lp == l ) then
 !                      call yama( vv_cmplx, k_rel, q_rel ) 
!                       end if
                       
                       

                       CALL g_interpolate(channel, ncoup, lp,l,js,k_rel,q_rel,vv_cmplx)
                       !CALL v_bare(lp,k_rel,l,q_rel,spin,js,isospin_tz,vv_cmplx)
                       !write(*,*) 'vv_cmplx', vv_cmplx
                       !infinity i vv_cmplx
                       
                       vv_cmplx = vv_cmplx*iph(ABS((l-lp))/2) 

                       DO lambda_bra=lam1,lam2
                          IF(triag(lambda_bra,la,lb)) CYCLE           ! triangular rels.
                          IF(triag(lambda_bra,lp,ll)) CYCLE
                          lsjj_coeff_bra = sum9j(ja,jb,la,lb,is1,is2,j_lab,lambda_bra,spin)* &
                               sum6j(ll,lp,lambda_bra,spin,j_lab,js)
                        !  write(*,*) 'lsjj',lsjj_coeff_bra ! Ingen uendeligheter
                          IF (ABS(lsjj_coeff_bra) < epsilon) CYCLE
                          vb_bra = vector_trcoefficients &
                           (k_rel,p_com,pa,pb,lp,ll,lambda_bra,la,lb,isospin_tz)/(pa*pb)
              !            write(*,*) 'vb_bra', real(vb_bra) ! ingen uendeligheter her
                          IF (ABS(vb_bra) < epsilon) CYCLE
                     !   write(*,*) 'HEI DETTE ER EN TEST' 
                          int_sum = int_sum+lsjj_coeff_bra*vv_cmplx*vb_bra
                   !     write(*,*)  real(vv_cmplx), real(vb_bra)
                       ENDDO
                    ENDDO
               !     if ( abs(int_sum) < epsilon ) cycle 
                 !  write(*,*)'int_sum', int_sum    
                    v_pot=v_pot + int_sum*vb_ket*lsjj_coeff_ket*weight_q(iq)*q_rel/p_com
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
  
 ! write(*,*) 'n_k, n_k1, n_k2' , n_k, n_k1, n_k2
  !  set in test here so that n_k = total in case of free and n_k1 in case of not free
 ! SELECT CASE(type_of_interaction)
 ! CASE('free')
     dim = n_k
 ! CASE('renorm')
 !    dim = n_k1
 ! END SELECT
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

  CALL lagrange_2dim(k_rel,q_rel,REAL(vfree(lim1:lim2,lim3:lim4,channel)), &
       v_interpol_re,dim,abs(krel))
    !  write(*,*) 'vfree', REAL(vfree(lim1:lim2,lim3:lim4,channel))
  
  vv = dcmplx( v_interpol_re, 0.d0 ) 
  
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

!
! yamaguchi potential
!
subroutine yama( vpot, k, q ) 
  use constants
  
  implicit none
  double precision, intent(in) :: k, q
  complex*16, intent(out) :: vpot
  double precision :: v0, beta

  beta = hbarc
  v0 = 10.d0 * hbarc
  ! 3.*(2./pi)*beta**3/M_np
  vpot = -v0/ (k**2 + beta**2)/ (q**2 + beta**2)

  
end subroutine yama


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
  V_a =  -2.79d0 !-2.6047d0 !5.42      
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


subroutine gaussian(v_gauss, l_ang, j_ang, k, q ) 
  use constants
  IMPLICIT NONE
  integer, intent(in) :: l_ang, j_ang
  complex*16, intent(out) :: v_gauss
  COMPLEX*16, INTENT(IN) :: k, q
  COMPLEX*16 :: x,  x1, x2, GAUSS, GAUSS1
  DOUBLE PRECISION :: a, v0, V1,B, Vc, Vsl, ac, asl, sl 
  !double precision, PARAMETER :: pi = 3.141592741012573d0
  
   if ( l_ang == 1 ) then
     vc =  -47.4/2.9  !47.32d0
     ac = 2.305d0 /hbarc 
     vsl = -5.86d0*2.d0/2.1 !11.71d0
     asl = 2.31d0 /hbarc
     
     !vc = 1.
     ! ac = 2.3
     !vc = -9.
     !ac =  1.d0 
     !!vsl = 0.d0
     !asl = 2.31d0

     !Vc = -47.4d0
     ! ac = 2.3d0
     !Vsl = -25.49d0
     !asl = 1.72d0
  elseif ( l_ang == 0 ) then
     vc = -47.4d0/2.9 !48.2d0
     ac =  2.305d0/hbarc   !2.33d0 
     vsl = 0.d0
  end if
     
  if ( j_ang == 1 ) sl = -1.d0
  if ( j_ang == 3 ) sl = 0.5d0
  
  x = k * k + q * q + 2.d0 * k * q 
  x1 =  k * k + q * q - 2.d0 * k * q
  ! x2 = k * q * a * a * 0.5d0
  
  
  if ( l_ang == 0 ) THEN
     
     GAUSS =  (SQRT(PI) * ac / ( 4.* k * q ) ) * ( EXP( -  x1 * ac * ac/ 4.d0 ) - &
          exp( -  x * ac * ac/ 4.d0 ) )
     
  elseif ( L_ang == 1 ) then
     
     GAUSS =  (sqrt(PI) * ac / ( 4.d0 * k  * q ) ) * ( EXP( -  x1 * ac * ac/ 4.d0 ) * &
          ( 1.d0 - 1.d0/ (k * q * ac * ac * 0.5d0)  ) +  exp( -  x * ac * ac/ 4.d0 ) * & 
          ( 1.d0 + 1.d0/ (k * q * ac * ac * 0.5d0) ) )
     GAUSS1 =   (sqrt(PI) * asl / ( 4.d0 * k  * q ) ) * ( EXP( -  x1 * asl * asl/ 4.d0 ) * &
          ( 1.d0 - 1.d0/ ( k * q * asl * asl * 0.5d0 ) ) + exp( - x * asl * asl/ 4.d0 ) * & 
          ( 1.d0 + 1.d0/ ( k * q * asl * asl * 0.5d0 ) ) )
  end if
  
  v_gauss =  Vc*GAUSS  + sl*Vsl*GAUSS1

ENd subroutine gaussian


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
  !complex*16  :: v,xmev,ymev
  double precision  :: v,xmev,ymev
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
  xmev= dble(k_rel)
  ymev =dble(q_rel) 
  !xmev=abs(k_rel)
  !ymev=abs(q_rel) 
  
  SELECT CASE (type_of_pot)
  CASE('CD-bonn')
!     CALL cdbonn
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
  CASE ( 'n3lo' ) 
     !write(6,*) 'n3lo'
     call n3lo
  END SELECT
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
  COMPLEX*16 :: zh, norm_sum, hole_func
  double precision :: ph, qmin, qmax, alpha
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: q_points, weight_q

  alpha = 1./oscillator_length/ hbarc

  DO nh=0, 0 !nh_max
     
     ph= (-1.D0)**nh
     DO lh=0,lh_max
        jh_min = abs(2*lh-1); jh_max = 2*lh+1
        do jh = jh_min, jh_max, 2
           do tzh = itzh_min, itzh_max, 2
              
              norm_sum = 0.
              DO i=1,n_lab
                 
                 zh = phase * k_mesh(i) * oscillator_length
                 if ( nh == 0 ) then
                    !write(6,*) lh,jh,tzh
                    bhf_hol(i, lh, (jh+1)/2, (tzh+1)/2 ) =  ph*rnl(nh,lh,zh)*dSQRT(oscillator_length**3) 
                   
                 end if
                 norm_sum = norm_sum + k_mesh(i)**2* w_mesh(i)*bhf_hol(i, lh, (jh+1)/2, (tzh+1)/2 )**2
                 !rnl(nh,lh,zh)**2*(oscillator_length**3) 
                 !write(6,*) k_mesh(i), w_mesh(i),ph*rnl(nh,lh,zh)*dSQRT(oscillator_length**3) 
              ENDDO
              WRITE(6,*) 'Norm rel ho wf n,l : ', lh,jh,tzh,norm_sum,nh,lh,alpha

!              norm_sum = 0.
!              do i = 1, nrad
!                 if ( nh == 0 ) then
!                    zh = dcmplx( r_lab(i)*alpha )
                    !bhf_hol_rspace(i, lh, (jh+1)/2, (tzh+1)/2 ) =  ph*rnl(nh,lh,zh)*dSQRT(alpha**3)
                    
!                    
!                    call interpol_1dim( nh,lh,jh, tzh, r_lab(i), hole_func)
!                    bhf_hol_rspace(i, lh, (jh+1)/2, (tzh+1)/2 ) =  hole_func
!                    norm_sum = norm_sum +  r_lab(i)**2*wr_lab(i)* hole_func**2
!                    
!                 end if
!              end do
!              
!              WRITE(6,*) 'Norm rel ho wf n,l : ', lh,jh,tzh,norm_sum,nh,lh,alpha
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
subroutine interpol_1dim( n,l,j,tz,rlab, wavefunc) 
  use constants
  USE wave_functions 
  use hf_constants
  
  implicit none
  integer, intent(in) :: n,l,j,tz
  double precision , intent(in) :: rlab
  complex*16, intent(out) :: wavefunc
  !double precision, dimension(n_lab) :: x, y
  double precision, dimension(size(k_ws)) :: x
  double precision :: re_s,im_s, qmin, qmax, xval
  double precision, allocatable,dimension(:) :: q_points
  double precision, allocatable, dimension(:) :: lam, C, wrk
  !double precision, dimension(n_lab) :: xx, yy
  double precision, dimension(size(k_ws)) :: xx, y, yy, re_y, im_y
  integer :: i, lck, lwrk, ifail
  
  !allocate( q_points(n_lab) )
  !qmin = 0.0_dp; qmax = hole_cutoff
  !q_points(:) = 0.5d0*( qmax+qmin ) + 0.5d0*( qmax-qmin )*x_mesh(:)
  
  !  write(6,*) klab, q_points
  if ( abs(rlab) - k_ws(1) < 0. ) write(6,*) 'wrong in 1dim interpolation'
  if ( abs(rlab) - k_ws(size(k_ws)) > 0. ) write(6,*) 'wrong in 1dim interpolation'
  
  x = k_ws
  
  do i = 1, size(k_ws)
     
     re_y(i) = dble( ws_basis(l,j,i) )
     im_y(i) = aimag( ws_basis(l,j,i) )
          
  end do
  
    
  !write(6,*) x, klab
  xx = (x); y = (re_y); yy = im_y
  
  !write(6,*) xx
  ! LCK > n_lab+4
  LCK = size(k_ws) + 4
  lwrk = 6*size(k_ws) + 16 
  ifail = 0
  allocate( lam(lck), c(lck), wrk(lwrk) )
  
  ! input: 
  ! constraint M > 4
  ! x(1) <...< x(M)  
  ! y(1)...y(M)
  ! LCK >= M + 4
  ! LWRK >= 6*M + 16
  ! E01BAF(M, x, y, lambda, C, LCK, WRK, LWRK, IFAIL)
  ! See documentation for further details. 
  xval = dble(rlab)
  
  call E01BAF(size(k_ws), x, y, lam, C, LCK, WRK, LWRK, IFAIL)
  call E02BBF(size(k_ws)+4, dble(lam), dble(C), xval, re_S, IFAIL)
  

  call E01BAF(size(k_ws), x, yy, lam, C, LCK, WRK, LWRK, IFAIL)
  call E02BBF(size(k_ws)+4, dble(lam), dble(C), xval, IM_S, IFAIL)
  
  wavefunc = dcmplx(re_s,im_s)
  
  deallocate( lam, c, wrk ) 

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
  double precision, dimension(n_lab) :: x, re_y, im_y
  double precision :: re_s, im_s, qmin, qmax, xval, check_sum
  double precision, allocatable, dimension(:) :: lam, C, wrk
  double precision, dimension(n_lab) :: xx, re_yy, im_yy
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
     re_y(i) = dble( bhf_hol(i, l, (j+1)/2 ,(tz+1)/2 ) )
     !im_y(i) = aimag( bhf_hol(i, l, (j+1)/2 ,(tz+1)/2 ) )
     !y(i) = dble( hf_holes(n, l, j, tz, i ) )
     !if ( y(i) == 0. ) then
     !write(6,*) l,j, x(i), y(i)
     !end if
     
  end do
  
  xx = (x); re_yy = re_y
  !; im_yy = im_y
  re_s = 0.d0
  !; im_s = 0.d0
  ! LCK > n_lab+4
  LCK = n_lab + 4
  lwrk = 6*n_lab + 16 
  ifail = 0

  allocate( lam(lck), c(lck), wrk(lwrk) )

  call E01BAF(n_lab, xx, re_yy, lam, C, LCK, WRK, LWRK, IFAIL)
  
  xval = abs(klab)
  call E02BBF(n_lab+4, dble(lam), dble(C), xval, re_S, IFAIL)
  
 ! if ( aimag(phase) /= 0.d0 ) then
 !    lam= 0.d0; c = 0.d0; wrk = 0.d0
 !    call E01BAF(n_lab, xx, im_yy, lam, C, LCK, WRK, LWRK, IFAIL)
 ! 
 !    xval = abs(klab)
 !    call E02BBF(n_lab+4, dble(lam), dble(C), xval, im_S, IFAIL)
 ! end if

  wavefunc = dcmplx(re_s,0.d0)
  

  deallocate( lam, c, wrk ) 

end subroutine interpol_1dim2

!
!
!
subroutine interp1dim_rspace( n,l,j,tz, rlab, wavefunc) 
  use constants
  USE wave_functions 
  use hf_constants
  
  implicit none
  integer, intent(in) :: n,l,j,tz
  double precision, intent(in) :: rlab
  complex*16, intent(out) :: wavefunc
  double precision, dimension(nrad) :: x, re_y, im_y
  double precision :: re_s, im_s, qmin, qmax, xval, check_sum
  double precision, allocatable, dimension(:) :: lam, C, wrk
  double precision, dimension(nrad) :: xx, re_yy, im_yy
  integer :: i, lck, lwrk, ifail
  

  write(6,*) rlab
  if ( abs(rlab) < r_lab(1)  ) write(6,*) 'wrong in 1dim interpolation'
  if ( abs(rlab) > r_lab(nrad) ) write(6,*) 'wrong in 1dim interpolation'
  
  do i = 1, nrad
    
     x(i) =  r_lab(i)
     re_y(i) = dble( bhf_hol_rspace(i, l, (j+1)/2 ,(tz+1)/2 ) )
     
  end do
  
  xx = (x); re_yy = re_y
  
  re_s = 0.d0
  !; im_s = 0.d0
  ! LCK > n_lab+4
  LCK = nrad + 4
  lwrk = 6*nrad + 16 
  ifail = 0

  allocate( lam(lck), c(lck), wrk(lwrk) )

  call E01BAF(nrad, xx, re_yy, lam, C, LCK, WRK, LWRK, IFAIL)
  
  xval = abs(rlab)
  call E02BBF(nrad+4, dble(lam), dble(C), xval, re_S, IFAIL)
  
  wavefunc = dcmplx(re_s,0.d0)
  

  deallocate( lam, c, wrk ) 

end subroutine interp1dim_rspace
!
!
! 
subroutine hartree_fock_final(n_osc, nlab_large)
  USE configurations
  USE constants
  use hf_constants
  use setup_cmplx_mesh
  use wave_functions
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n_osc, nlab_large
  INTEGER :: a, c, ij, ipar, itz, ket , bra
  INTEGER :: na,nc, la, lc, ja, jc, tza, tzc, iph
  REAL(DP) ::  vpot_re, vpot_im, pa_real, wpa_real, pc_real, wpc_real, kmin, kmax, fxy_re, pa_temp, pc_temp, ph
  double precision, allocatable, dimension(:) :: k_input, klab_large, wlab_large
  complex(dpc) :: diaghf, klab, norm_sum,zh, pa, pc, wpa, wpc
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:) :: hf_mtx, hf_vect
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:) :: hf_eigen, e_temp
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: vlab_interp, vlab_input, vlab_rspace
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:,:) :: ho_osc_func
  double precision, allocatable, dimension(:,:) :: vint_re, vint_im
  integer :: i, j,  number, nstep, nsmall, la_temp, ja_temp, tza_temp
  integer :: ntot, i1,i2,l1,j1,tz1
  double precision, allocatable, dimension(:) ::X,Y,XX, WORK,W, D, rad, wrad
  double precision :: val_re, val_im,vall, vall_re, vall_im, b_trans, angle
  integer:: ifail, IG1
  integer :: nh, lh, jh,tzh, nrad_lab, r1, r2, nk1, nk2, nk_tot
  complex(dpc) :: kinetic, pot, cmplx_phase, cmplx_oscl, cmplx_theta, bessel_p1, bessel_p2, kr, qr, int_sum
  real*8 :: e_re, e_im
  character(len=100) :: contour

  !
  ! in data file: hf_pot_'nucleus'_'jrel_max'_'n_lab'_'int_points'_'cutoff'.dat
  !
  open(10, file='hfpot_n3lo_016_4_16_14_20.dat')
  !open(10,file='test.dat')
  !open(10,file='hf_pot_O16_6_16_14_21.dat')
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
     
     read(10,*) bra, la, ja, tza,pa_real, ket, lc, jc, tzc,pc_real, vpot_re, vpot_im
     
     if ( bra > nsmall ) nsmall = bra
     if ( la > lh_max ) lh_max = la
     if ( ja > jh_max ) jh_max = ja
  end do
  rewind 10

  write(6,*) lh_max, jh_max 
  
  

  IG1 = nsmall
  ALLOCATE ( vlab_input(nsmall,nsmall, 0:lh_max, (jh_min+1)/2:(jh_max+1)/2 ,(-1+1)/2:(1+1)/2) )
  ALLOCATE ( k_input(nsmall) ) 
  allocate( x(nsmall), y(nsmall) )
  allocate( XX(nsmall), WORK(nsmall),W(nsmall), D(nsmall) )
  allocate( vint_re(nsmall, nsmall), vint_im(nsmall, nsmall) )

  write(6,*) jh_min, jh_max,lh_min, lh_max
  write(6,*) 'min/max s.p. quantum number: j,l,itz', jh_min, jh_max, lh_min,lh_max, itzh_min,itzh_max
  do i = 1, number
     read(10,*) bra, la, ja, tza,pa_real, ket, lc, jc, tzc,pc_real, vpot_re,  vpot_im
     i1 = bra; i2 = ket; l1 =la;j1=ja;tz1=tza
     
     vlab_input(i1,i2,l1,(j1+1)/2,(tz1+1)/2) = dcmplx( vpot_re, 0.d0 )
     
     if ( bra /= ket ) then
        vlab_input(i2,i1,l1,(j1+1)/2,(tz1+1)/2) = vlab_input(i1,i2,l1,(j1+1)/2,(tz1+1)/2) 
     end if
     k_input(bra) = pa_real
  end do
  kmin = k_input(1)
  kmax = k_input(nsmall)

  
  allocate ( klab_large(nlab_large), wlab_large(nlab_large) )
  ALLOCATE ( vlab_interp(nlab_large, nlab_large, 0:lh_max, & 
       (jh_min+1)/2:(jh_max+1)/2 ,(-1+1)/2:(1+1)/2) )
 

  call gauss_legendre(kmin,kmax,klab_large,wlab_large,nlab_large )
  
  
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
           
           vint_re = 0.
           vint_im = 0.
           vint_re(:,:) = real( vlab_input(:,:,la,(ja+1)/2,(tza+1)/2 ) )
           vint_im(:,:) = aimag( vlab_input(:,:,la,(ja+1)/2,(tza+1)/2 ) )
           
           DO bra = 1, nlab_large
              pa_real = klab_large(bra); wpa_real = wlab_large(bra)
              DO ket= bra, nlab_large
                 pc_real = klab_large(ket); wpc_real = wlab_large(ket)

                 call E01ACF((pa_real), (pc_real), (k_input), (k_input), vint_re, & 
                      VAL_re, VALL_re, IFAIL, XX, WORK,W, D, IG1, nsmall,nsmall)
                 if ( phase /= 1.d0 ) then
                    call E01ACF((pa_real), (pc_real), (k_input), (k_input), vint_im, & 
                      VAL_im, VALL_im, IFAIL, XX, WORK,W, D, IG1, nsmall,nsmall)
                 end if
                 
                 vlab_interp(bra,ket,la,(ja+1)/2,(tza+1)/2) = dcmplx( vall_re, vall_im ) !dcmplx(fxy_re)
                 
                 if ( bra /= ket ) then
                    vlab_interp(ket,bra,la,(ja+1)/2,(tza+1)/2) = vlab_interp(bra,ket,la,(ja+1)/2,(tza+1)/2)
                 end if
                 
                 
              end DO
              !write(6,*) pa, vall_re
           end DO
        end DO
     end do
  end DO

  !la = 0; ja = 1; tza = 1
  !do bra = 1, nlab_large
  !   write(6,*) klab_large(bra),dble( vlab_interp(bra,2,la,(ja+1)/2,(tza+1)/2) )
  !end do
  !return
  
  nrad_lab = 40
  allocate(rad(nrad_lab), wrad(nrad_lab) )
  call gauss_legendre(0.d0,8.d0,rad,wrad,nrad_lab )
  
  ALLOCATE ( vlab_rspace(nrad_lab,nrad_lab, 0:lh_max, (jh_min+1)/2:(jh_max+1)/2 ,(-1+1)/2:(1+1)/2) )
  !
  ! here we set up r-space hf-potentials
  ! 
  ! loop over single-particle isospin
  vlab_rspace = 0.d0
  DO tza = -1, 1, 2
     ! loop over single-particle orbital angular momentum
     do la = lh_min, lh_max, 1 
        !     loop over single-particle angular momentum (twice the value)
        DO ja = jh_min, jh_max, 2 
           
           if ( ja > ( 2*la + 1) ) cycle
           if ( ja < abs( 2*la - 1 ) ) cycle
           
           if ( la /= 2 ) cycle
           if ( ja /= 5 ) cycle
           !if ( tza /= -1 ) cycle
           !           write(6,*) la,ja,tza
           do r1 = 1, nrad_lab
              !r2 = r1 
              do r2 = 1, r1 !nrad
                 int_sum = 0.0
                 DO bra = 1, nlab_large
                    pa = klab_large(bra); wpa = wlab_large(bra)
                    kr = rad(r1)*pa/hbarc
                    call CMPLX_BESSELJ(la,kr,BESSEL_p1)
                    
                    DO ket= 1, nlab_large
                       pc = klab_large(ket); wpc = wlab_large(ket)
                       qr = rad(r2)*pc/hbarc
                       call CMPLX_BESSELJ(la,qr,BESSEL_p2)
                       
                       int_sum  = int_sum + (2./pi)*wpa*pa**2*wpc*pc**2* & 
                            ( vlab_interp(bra,ket,la,(ja+1)/2,(tza+1)/2) ) * bessel_p1*bessel_p2/hbarc**3
                       
                    end DO
                 end DO
                 vlab_rspace(r1,r2, la, (ja+1)/2 ,(tza+1)/2)  = int_sum

                 vlab_rspace(r2,r1, la, (ja+1)/2 ,(tza+1)/2) = & 
                      vlab_rspace(r1,r2, la, (ja+1)/2 ,(tza+1)/2)
                              
                 if ( r1 == r2 .and. la == 0 ) then
                    write(6,*) la,ja,tza,rad(r1), dble(int_sum)
                 end if
                 
              end do

           end do
        end DO
     end do
  end DO

  int_sum = 0.0
  do r1 = 1, nrad_lab
     do r2 = 1, nrad_lab
        
        int_sum = int_sum + vlab_rspace(r1,r2, 0, (1+1)/2 ,(-1+1)/2)* rad(r1)**2*rad(r2)**2*wrad(r1)*wrad(r2)
     end do
  end do
  write(6,*) 'jaj', int_sum


  !
  ! here we set up integration contour in the complex k-plane 
  !  
  contour = 'triangle'
  nk1 = 4
  nk2 = 15
  nk_tot = 2*nk1 + nk2 
  allocate( k_cmplx(nk_tot), wk_cmplx(nk_tot) )
  angle = pi/5.d0
  b_trans = 0.6
  call GENERAL_MESH( nk1, nk2, angle, b_trans, contour )

  DEALLOCATE ( vlab_interp )
  ALLOCATE ( vlab_interp(nk_tot, nk_tot, 0:lh_max, & 
       (jh_min+1)/2:(jh_max+1)/2 ,(-1+1)/2:(1+1)/2) )
 
  
  !
  ! 
  ! here we set up potential in the complex k-plane 
  ! 
  ! loop over single-particle isospin
  DO tza = -1, 1, 2
     ! loop over single-particle orbital angular momentum
     do la = lh_min, lh_max, 1 
        !     loop over single-particle angular momentum (twice the value)
        DO ja = jh_min, jh_max, 2 
           
           if ( ja > ( 2*la + 1) ) cycle
           if ( ja < abs( 2*la - 1 ) ) cycle
           
           if ( la /= 2 ) cycle
           if ( ja /= 5 ) cycle
           !if ( tza /= 1 ) cycle
          
           !if ( la /= 2 .and. ja /= 3 .and. tza /= 1 ) cycle
           pa = 0.0
           DO bra = 1, nk_tot
              pa = k_cmplx(bra); wpa = wk_cmplx(bra)
                            
              pc = 0.0
              DO ket= 1, bra !nlab_large
                 pc = k_cmplx(ket); wpc = wk_cmplx(ket)
                 
                 int_sum = 0.0
                 do r1 = 1, nrad_lab
                    kr = rad(r1)*pa/hbarc
                    call CMPLX_BESSELJ(la,kr,BESSEL_p1)
                    
                    
                    do r2 = 1, nrad_lab
                       qr = rad(r2)*pc/hbarc
                       call CMPLX_BESSELJ(la,qr,BESSEL_p2)
                       
                       int_sum  = int_sum + (2.d0/pi)*wrad(r1)*rad(r1)**2*wrad(r2)*rad(r2)**2 * & 
                            vlab_rspace(r1,r2, la, (ja+1)/2 ,(tza+1)/2) * bessel_p1*bessel_p2/hbarc**3
                       !int_sum  = int_sum + (2.d0/pi)*wrad(r1)*rad(r1)**2 * & 
                       !     vlab_rspace(r1,r2, la, (ja+1)/2 ,(tza+1)/2) * bessel_p1*bessel_p2/hbarc**3
                     
                    end DO
                 end DO
                 
                 vlab_interp(bra,ket,la,(ja+1)/2,(tza+1)/2) = int_sum
                 vlab_interp(ket,bra,la,(ja+1)/2,(tza+1)/2) =vlab_interp(bra,ket,la,(ja+1)/2,(tza+1)/2) 
                 if ( bra == ket ) then
                    write(6,*) la,ja,tza, dble(pa),aimag(pa), dble(int_sum),aimag(int_sum)
                 end if
                 
              end do

           end do
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
           
 
           if ( la /= 2 ) cycle
           if ( ja /= 5 ) cycle
           !if ( tza /= 1 ) cycle
          

           ALLOCATE(hf_mtx(nk_tot,nk_tot) )
           ALLOCATE(hf_vect(nk_tot,nk_tot) )
           ALLOCATE(hf_eigen(nk_tot) )
                     
           !     set up hf matrix to be diagonalized
           hf_mtx = 0.0_dp; hf_vect = 0.0_dp; hf_eigen = 0.0_dp
           DO bra = 1, nk_tot
              pa = k_cmplx(bra); wpa = wk_cmplx(bra)
              DO ket= bra, nk_tot
                 pc = k_cmplx(ket); wpc = wk_cmplx(ket)
                 
                 
                 diaghf = ( vlab_interp(bra,ket,la,(ja+1)/2,(tza+1)/2) )
                 !v_bhf( bra, ket, la,(ja+1)/2 ,(tza+1)/2 )
                 ! Interaction in units of MeV and preparing for diagonalization
                 hf_mtx(bra,ket) = pa*pc*diaghf*SQRT(wpa*wpc)
                 hf_mtx(ket,bra)=hf_mtx(bra,ket)
                 ! Add kinetic energy
                 IF ( bra == ket) THEN
                    hf_mtx(bra,bra)=hf_mtx(bra,bra) + pa*pa*0.5_dp/p_mass(tza)*(1.d0-1.d0/mass_closed_shell_core)
                 ENDIF
              ENDDO
           ENDDO
           ! diagonalize and get new eigenvalues and energies
           CALL diag_exact( hf_mtx, hf_vect, hf_eigen, nk_tot)
           
           write(6,*) 'sp quantum numbers:tza,ja,la and energy',tza, ja, la
           do i = 1, nk_tot
              
           !   if ( dble( hf_eigen(i) )< 0.d0 ) then
                 e_re = real( hf_eigen(i) ); e_im = real(aimag( hf_eigen(i) )) 
           !      write(6,*) 'sp quantum numbers:tza,ja,la and energy',tza, ja, la,e_re, e_im
                
           !   end if
                 write(6,*) e_re, e_im
           end do

           DEALLOCATE ( hf_mtx); DEALLOCATE ( hf_eigen)
           DEALLOCATE ( hf_vect)
           
        ENDDO
     ENDDO
  ENDDO
  
  return

  !
  ! bare glemm oscillator basis og det nedenfor her...
  !

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
        !WRITE(6,*) 'Norm rel ho wf n,l : ', nh, lh,norm_sum 
     end do
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
                 write(6,*) 'sp quantum numbers:tza,ja,la and energy',tza, ja, la,real( hf_eigen(i) ), &
 aimag(hf_eigen(i))  
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
!
!
!
subroutine check_vector_bras( n,l,NN,LL, n1,l1,n2,l2, lam )
  USE constants
  USE wave_functions
  use ang_mom_functions
  
  IMPLICIT NONE

  integer, intent(in) :: n,l,NN,LL, n1,l1,n2,l2, lam
  integer :: lab_dim, rel_dim, i, j, i1, j1
  double precision, allocatable, dimension(:) :: k_lab, wk_lab, k_rel,wk_rel 
  double precision :: pcom, k, qmin, qmax, k1, k2, alpha, alpha_k, alpha_P, int_sum, vb_bra
  complex*16 :: zk1, zk2, zk,zp
  double precision :: osc1, osc2, osc_k, osc_p, wk, wk1, wk2, int_sum2, q1min, q1max, d

  lab_dim = 50
  rel_dim = 50
  allocate( k_rel(rel_dim), wk_rel(rel_dim) )
  allocate( k_lab(lab_dim), wk_lab(lab_dim) )


  ! oscillator length in lab
  alpha = 1./hbarc !oscillator_length
  ! oscillator length in rel
  alpha_k = dsqrt(2.d0)*alpha 
  alpha_P = dsqrt(0.5d0)*alpha

  write(6,*) alpha, alpha_k, alpha_P
  qmin = 0.d0
  qmax = 10.d0
  CALL gauss_legendre(qmin*hbarc,qmax*hbarc,k_lab,wk_lab,lab_dim)
  
  write(6,*)n,l,NN,LL, n1,l1,n2,l2, lam 

  int_sum = 0.d0
  do i = 1, lab_dim
     k1 = k_lab(i)
     wk1 = wk_lab(i)
     zk1 = k1*alpha
     
     osc1=sqrt( alpha**3 )*dble( rnl(n1,l1,zk1) )
     
     !
     
     do j = 1, lab_dim
        k2 = k_lab(j)
        wk2 = wk_lab(j)
        zk2 = k2*alpha
        osc2=sqrt( alpha**3 )*dble( rnl(n2,l2,zk2) )
      

        q1min = 0.5*abs(k1-k2)
        q1max = 0.5*abs(k1+k2)
        
        
        CALL gauss_legendre(q1min,q1max,k_rel,wk_rel,rel_dim)
        
        do i1 = 1, rel_dim
           k = k_rel(i1) !k_rel(i1)
           wk = wk_rel(i1) !wk_rel(i1)
           zk = k_rel(i1)* alpha_k !dsqrt(2.d0)*oscillator_length
           osc_k=sqrt( alpha_k**3 )*dble( rnl(n,l,zk) )
           
           if ( 0.5*k1**2+0.5*k2**2-k*k < 0.d0 ) cycle
           
           Pcom = 2.d0*sqrt(0.5*k1**2+0.5*k2**2-k*k) 
           zP= Pcom*alpha_P
           
           osc_P=sqrt( alpha_P**3 )*dble( rnl(NN,LL,zP) )
           
           vb_bra = vector_trcoefficients(k,pcom,k1,k2,l,ll,lam,l1,l2,1) 
           
           int_sum = int_sum + 2.d0*k1*wk1*k2*wk2*k*wk*osc1*osc2*osc_k *osc_P* vb_bra
           
        end do

     end do
  end do

  int_sum2 = gmosh(n,l,nn,ll,n1,l1,n2,l2,lam,1.d0)
  
  write(6,*) 'moshinsky bracket', int_sum2
  write(6,*) 'moshinsky bracket', int_sum
  
  call check_vector_bras_rspace( n,l,NN,LL, n1,l1,n2,l2, lam )
  
  deallocate( k_rel, wk_rel, k_lab, wk_lab )
  
end subroutine check_vector_bras

!
!
!
subroutine check_vector_bras_rspace( n,l,NN,LL, n1,l1,n2,l2, lam )
  USE constants
  USE wave_functions
  use ang_mom_functions
  
  IMPLICIT NONE

  integer, intent(in) :: n,l,NN,LL, n1,l1,n2,l2, lam
  integer :: lab_dim, rel_dim, i, j, i1, j1
  double precision, allocatable, dimension(:) :: rlab, wrlab, r_rel,wr_rel 
  double precision :: rcom, r, rmin, rmax, r1, r2, alpha, alpha_r, alpha_cm, int_sum, vb_bra
  complex*16 :: zr1, zr2, zr,zcm
  double precision :: osc1, osc2, osc_r, osc_CM, wr, wr1, wr2, int_sum2, q1min, q1max, d

  lab_dim = 70
  rel_dim = 70
  allocate( r_rel(rel_dim), wr_rel(rel_dim) )
  allocate( rlab(lab_dim), wrlab(lab_dim) )


  ! oscillator length in lab
  alpha = 1./oscillator_length/hbarc
  ! oscillator length in rel
  alpha_r = dsqrt(2.d0)*alpha 
  alpha_cm = dsqrt(0.5d0)*alpha

  write(6,*) alpha, alpha_r, alpha_cm, oscillator_length*hbarc
  rmin = 0.d0
  rmax = 10.d0
  CALL gauss_legendre(rmin,rmax,rlab,wrlab,lab_dim)
  
  write(6,*)n,l,NN,LL, n1,l1,n2,l2, lam 

  int_sum = 0.d0
  do i = 1, lab_dim
     r1 = rlab(i)
     wr1 = wrlab(i)
     zr1 = r1*alpha
     
     osc1=sqrt( alpha**3 )*dble( rnl(n1,l1,zr1) )
     
          
     !
     
     do j = 1, lab_dim
        r2 = rlab(j)
        wr2 = wrlab(j)
        zr2 = r2*alpha
        osc2=sqrt( alpha**3 )*dble( rnl(n2,l2,zr2) )
      
        
        
        q1min =  0.5*abs(r1-r2)
        q1max =  0.5*abs(r1+r2)
        
        
        CALL gauss_legendre(q1min,q1max,r_rel,wr_rel,rel_dim)
      
        do i1 = 1, rel_dim
           r = r_rel(i1) 
           wr = wr_rel(i1) 
           zr = r_rel(i1) *alpha_r
           osc_r=sqrt( alpha_r**3 )*dble( rnl(n,l,zr) )
           
           if ( 0.5*r1**2+0.5*r2**2-r*r < 0.d0 ) cycle
           
           rcom = 2.d0*sqrt(0.5*r1**2+0.5*r2**2-r*r) 
           zcm= rcom*alpha_cm
           
           osc_cm=sqrt( alpha_cm**3 )*dble( rnl(NN,LL,zcm) )
           
           vb_bra = vector_trcoefficients(r,rcom,r1,r2,l,ll,lam,l1,l2,1)
           
           !           int_sum = int_sum + r1**2*wr1*osc1**2 * r2**2*wr2*osc2**2 * wr *r**2 * osc_r**2 
           !           int_sum = int_sum +  wr *r**2 * osc_r**2 
           int_sum = int_sum + 2.d0*r1*wr1*r2*wr2*wr*r*osc1*osc2*osc_r *osc_cm* vb_bra
           

        end do
        
     end do
  end do

  int_sum2 = gmosh(n,l,nn,ll,n1,l1,n2,l2,lam,1.d0)
  
  write(6,*) 'moshinsky bracket', int_sum2
  write(6,*) 'moshinsky bracket', int_sum
  
  deallocate( r_rel, wr_rel, rlab, wrlab )
  
end subroutine check_vector_bras_rspace
