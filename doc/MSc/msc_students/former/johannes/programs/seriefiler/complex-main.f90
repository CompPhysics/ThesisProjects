
!             Program block complex-main.f90    
!
!             Authors:  Gaute Hagen and Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO 
!             E-MAIL:   mhjensen@fys.uio.no
!             LANGUAGE: F90 
!             LAST UPGRADE : February 2005
!             Program to calculate effective interactions using 
!             a similarity transformed 
!             interactions  for finite nuclei employing plane waves only
!             The program runs in the proton-neutron formalism. 
!
!     Main program starts here
!

PROGRAM MAIN
  USE ang_mom_functions
  USE constants
  USE wave_functions
  USE partial_waves
  USE configurations
  USE single_particle_orbits
  USE relcm_gmatrix
  use hf_constants
  
  IMPLICIT NONE
  INTEGER :: loop, number_orbits, max_conf
  integer :: i, iteration_number
  double precision :: sigma
  logical :: selfconsistency
  INTEGER :: number_mtxel, n_osc, nlab_large

  
  
  !     Read in all relevant data
  CALL effective_interaction_data


  !     set up mesh points in rel frame for integration for computation of
  !     bare interaction or similarity transformed interaction
  ALLOCATE ( krel (n_k), wkrel (n_k))
  CALL rel_mesh
  ! cutoff in vlowk 
  Lambda = abs( krel(n_k1) ) 

  ! number of rspace points
  nrad = 30
  
  !     set up mesh points in lab frame
  ALLOCATE ( k_mesh (n_lab), w_mesh (n_lab))
  ALLOCATE ( r_lab (nrad), wr_lab (nrad))
  
  ! find max value of lab momentum based on cutoff in vlwok
  hole_cutoff = min( 2.d0*lambda, 5.d0*hbarc )  
  
  
  
  ! setup mesh points for 
  CALL lab_mesh
  
  !call rspace_lab_mesh
  
  !     read sp data and set up their variables plus com and rel momenta 
  CALL read_sp_data

  
  SELECT CASE (type_of_calculation)
  CASE ('hf_setup')
     ! calculate hartree-fock potential self consistently using plane-wave basis
     call hartree_fock_self_consistent
  CASE ('hf_final') 
     ! calculate hartree-fock single particle spectrum using oscillator/plane-wave basis
     ! self consistent hartree-fock potential is read from disk
     n_osc = 20; nlab_large = 40
     call  hartree_fock_final(n_osc, nlab_large)
  END SELECT

  DEALLOCATE ( krel, wkrel)

END PROGRAM MAIN
!
! calculate hartree-fock potential self consistently using plane-wave basis
! 
subroutine  hartree_fock_self_consistent
  USE ang_mom_functions
  USE constants
  USE wave_functions
  USE partial_waves
  USE configurations
  USE single_particle_orbits
  USE relcm_gmatrix
  use hf_constants
  
  IMPLICIT NONE
  INTEGER :: loop, number_orbits, max_conf
  integer :: i, iteration_number
  double precision :: sigma
  logical :: selfconsistency
  INTEGER :: number_mtxel

  
  
  
  !     Factorials for 3j, 6j and 9j symbols and for the moshinsky/vector
  !     bracket  transf coeffs  
  CALL commons_to_angmom  
  !     find all partial waves with given jmin jmax
  CALL setup_channels     
  !     Setup of effective similarity transformed interaction
  !     loop over channels in rel&cm system  
  !     allocate space for the free NN interaction, allowing for coupled channels
  SELECT CASE(type_of_interaction)
     CASE('free')
       ALLOCATE(vfree(2*n_k,2*n_k,no_channels))
     CASE('renorm')
       ALLOCATE(vfree(2*n_k1,2*n_k1,no_channels))
  END SELECT
  vfree = 0.0_dp
  DO loop=1, no_channels    
     CALL g_channel(loop)
  ENDDO

  !   Perform single-particle self-consistency 
  
  allocate( x_mesh(int_points), wx_mesh(int_points) )
  ALLOCATE( t_mesh (int_points), wt_mesh (int_points))
  call hf_mesh
  call cheb_mesh
  
  jh_min = 1
  jh_max = 1
  lh_min = 0
  lh_max = 0
  itzh_min = 0
  itzh_max = 0
  do i= 1, all_orbit%total_orbits
     if ( all_orbit%orbit_status(i) == 'outside' ) cycle
     
     if ( all_orbit%jj(i) < jh_min ) then
        jh_min = all_orbit%jj(i)
     end if
     if ( all_orbit%jj(i) > jh_max ) then
        jh_max = all_orbit%jj(i)
     end if
     if ( all_orbit%ll(i) < lh_min ) then
        lh_min = all_orbit%ll(i)
     end if
     if ( all_orbit%ll(i) > lh_max ) then
        lh_max = all_orbit%ll(i)
     end if
     
     if ( all_orbit%itzp(i) < itzh_min ) then
        itzh_min = all_orbit%itzp(i)
     end if
     if ( all_orbit%itzp(i) > itzh_max ) then
        itzh_max = all_orbit%itzp(i)
     end if
     
  end do
  
  
  write(6,*) 'min/max s.p. quantum number: j,l,itz', jh_min, jh_max, lh_min,lh_max, itzh_min,itzh_max
  
  
  ALLOCATE ( v_bhf(n_lab, n_lab, 0:lh_max, & 
       (jh_min+1)/2:(jh_max+1)/2 ,(-1+1)/2:(1+1)/2) )
  ALLOCATE ( bhf_hol(n_lab, 0:lh_max, & 
       (jh_min+1)/2:(jh_max+1)/2 ,(-1+1)/2:(1+1)/2) )
  ALLOCATE ( bhf_hol_temp(n_lab, 0:lh_max, & 
       (jh_min+1)/2:(jh_max+1)/2 ,(-1+1)/2:(1+1)/2) )
  ALLOCATE ( bhf_hol_rspace_temp(nrad, 0:lh_max, & 
       (jh_min+1)/2:(jh_max+1)/2 ,(-1+1)/2:(1+1)/2) )
  ALLOCATE ( bhf_hol_rspace(nrad, 0:lh_max, & 
       (jh_min+1)/2:(jh_max+1)/2 ,(-1+1)/2:(1+1)/2) )
  
  ! allocate space for interaction in lab coords
  
  !setup woods-saxon hole functions
  !call setup_woods_saxon

  !stop
  
  selfconsistency = .FALSE. 
  iteration_number = 1
  
  ! count number of matrix elements <(pa,la,ja),(pb,lb,jb)|H|(pc,lc,jc),(pd,ld,jd)>
  ! needed in Hartree-Fock calculation.
  call count_number_mtxel(number_mtxel)
  allocate( vlab_mom_re(number_mtxel), vlab_mom_im(number_mtxel) )
  !     initialize bhf harmonic oscillator sp wave function
  bhf_hol=0.

  ! setup harmonic osc. functions for 1st iteration
  call HO_OSC_FUNC
  write(6,*) 'number of matrix elements calculated', number_mtxel
  write(6,*) 'new limits k_interpol_min/max', k_interpol_min, k_interpol_max
  
  
  sigma = 10.; iteration_number = 1
  DO WHILE ( iteration_number < hf_iterations .and. (ABS(sigma) > 1.E-5) ) 
     
     !     start BHF calculation
     WRITE(6,*) 'Sigma for this iteration', sigma, iteration_number
     v_bhf = 0. 
     CALL hartree_fock(iteration_number, sigma,selfconsistency )
          
     bhf_hol = bhf_hol_temp 
     !bhf_hol_rspace = bhf_hol_rspace_temp 

     iteration_number = iteration_number+1
  ENDDO
  
 
  selfconsistency = .TRUE. 
  !
  CALL hartree_fock(iteration_number, sigma,selfconsistency )
  

  
  DEALLOCATE ( bhf_hol, bhf_hol_temp )
  DEALLOCATE ( vlab_mom_re,vlab_mom_im )
  !CALL hartree_fock
  DEALLOCATE ( k_mesh, w_mesh)
  deallocate( x_mesh, wx_mesh, t_mesh, wt_mesh )
  DEALLOCATE(vfree)

END subroutine hartree_fock_self_consistent
!
!
!
subroutine setup_woods_saxon
  USE constants
  USE wave_functions
  use hf_constants

  implicit none
  integer :: i, la, ja, ntot
  complex*16, allocatable, dimension(:) :: wavefunc, ws_energy
  double precision, allocatable ::  re_z(:), im_z(:)
  integer :: n0,n2,n1, number, ja_min, ja_max, la_min, la_max, no_orbits
  complex*16 :: energy

  open(7,file='woods_saxon_basis.dat')
  ! count number of data lines in woods-saxon basis
  number=0
3 READ(7,*,end=4)
  number=number+1
  GOTO 3
4 CONTINUE
  REWIND 7
  !write(6,*) 'qbox number' , number
  
  ntot = number
  allocate( wavefunc(ntot), re_z(ntot), im_z(ntot) )
  
  !find lmin/lmax and jmin/jmax in woods-saxon basis
  la_min = 0; la_max = 0
  ja_min = 1; ja_max = 1 
  do i = 1, ntot
     read(7,*)  la, ja
     
     if ( la > la_max ) la_max = la
     if ( ja > ja_max ) ja_max = ja
     
  end do
  write(6,*) la_min, la_max, ja_min,ja_max
  
  no_orbits = 0
  do la = la_min, la_max
     do ja = ja_max, ja_min, -2
        
        if ( ja < abs( 2*la-1 ) ) cycle
        if ( ja > (2*la+1) ) cycle
        no_orbits = no_orbits + 1
     end do
  end do
  
  write(6,*) no_orbits
  
  ! only s1/2 and p3/2 orbits
  no_orbits = 2
  allocate( ws_basis(la_min:la_max, ja_min:ja_max, ntot/no_orbits), k_ws(ntot/no_orbits) )
  allocate( ws_energy( no_orbits) )
  
  rewind 7
  
  n0 = 0; n1 = 0; n2 = 0
  do i = 1, ntot
     read(7,*)  la,ja, energy, re_z(i), wavefunc(i)  
     if ( la == 0 .and. ja == 1 ) then
        n0 = n0 + 1
        ws_basis(la,ja,n0) = (wavefunc(i))
        k_ws(n0) = re_z(i)
        ws_energy(1) = energy
        !write(6,*) k_ws(n0), dble(wavefunc(i))
     elseif ( la ==1 .and. ja == 3 ) then
        n1 = n1 + 1
        ws_basis(la,ja,n1) = (wavefunc(i))
        ws_energy(2) = energy
        !write(6,*) dble(wavefunc(i))
     elseif ( la ==1 .and. ja == 1 ) then
        n2 = n2 + 1
        ws_basis(la,ja,n2) = (wavefunc(i))
        ws_energy(3) = energy
     end if
  end do
  
  k_interpol_min = k_ws(1)
  k_interpol_max = k_ws(ntot/no_orbits)
  
  no_orbits = 0
  do la = la_min, la_max
     do ja = ja_max, ja_min, -2
        
        if ( ja < abs( 2*la-1 ) ) cycle
        if ( ja > (2*la+1) ) cycle
        
        if ( ja == 1 .and. la == 1 ) cycle
        no_orbits = no_orbits + 1
        write(6,*) 'woods-saxon data: energy, la,ja', ws_energy(no_orbits), la, ja  
     end do
  end do
  
  write(6,*) 'k_interpol_min/max for woods-saxon basis', k_interpol_min, k_interpol_max
  
end subroutine setup_woods_saxon

!
!     This subroutine reads in variables needed to specify the
!     effective interaction
!
SUBROUTINE effective_interaction_data
  USE constants
  USE wave_functions 
  IMPLICIT NONE
  CHARACTER (LEN=100) :: results_file
  REAL(DP) :: pi_fraction

  !     Read in file names for effective interaction, q-box and g-matrix
  READ(5,*)
  READ(5,*) results_file
  OPEN (UNIT=10,FILE=results_file)
  READ(5,*);READ(5,*) type_of_calculation
  READ(5,*); READ(5,*) itzmin, itzmax, j_lab_min, j_lab_max
  WRITE(6,*) ' min & max Tz and J: ', itzmin, itzmax, j_lab_min, j_lab_max
  READ(5,*); READ(5,*) mass_closed_shell_core
  !     set up the oscillator length and energy
  WRITE(6,*) 'Mass of closed shell core', mass_closed_shell_core
  oscillator_energy=45.0_dp*DFLOAT(mass_closed_shell_core)**(-1./3.)- &  
       25.0_dp*DFLOAT(mass_closed_shell_core)**(-2./3.)  
  oscillator_length = 1./SQRT(p_mass(0)*oscillator_energy)   !  osc length in MeV^-1
  WRITE(6,*) 'Oscillator energy and length ', oscillator_energy, oscillator_length
  !hole_cutoff = 1.0_dp/oscillator_length
  ! the number of integration mesh points are fixed by the routine setup_hole_cutoff
  !WRITE(6,*) 'Trial number of  integration points for intermediate states', int_points
  !     limit on lab momentum in fm^-1
  READ(5,*) ; READ(5,*) k_cutoff, k_max
  write(6,*) 'k_cutoff and k_max for models space', k_cutoff, k_max
  READ(5,*)
  READ(5,*) jmin, jmax            ! min and max value of J(L) for NN partial waves
  lmax = jmax+1 ! applies only to nuclear case with S_max = 1
  WRITE(6,*) 'Min and max value of partial wave ang. mom'
  WRITE(6,*) jmin, jmax
  READ(5,*)
  READ(5,*) n_lab1,n_lab2,n_lab3 ! nr. of mesh points in lab system
  n_lab = n_lab1+n_lab2+n_lab3
  READ(5,*) 
  READ(5,*) int_points ! nr. of mesh points for integration of hole functions
  WRITE(6,*) 'number of lab mesh points and integration points ', n_lab, int_points
  READ(5,*);READ(5,*) hf_iterations !  maximum number of iterations in self consisten solution of hf-pot
  ! nr of mesh points along line L1 and line L2
  READ(5,*);READ(5,*) n_k1,  n_k2
  n_k = n_k1+n_k2
  WRITE(6,*) 'number of CoM system points along L1, L2 and total', n_k1, n_k2, n_k
  ! complex scaling angle theta = pi/pi_fraction
  READ(5,*); READ(5,*) pi_fraction
  ! rotation angle in coordinate space pi > theta > 0
  ! corresponding to a rotation of the contour in momentum space by -theta
  theta = 0.0_DP
  IF ( pi_fraction /= 0.0_DP) theta = pi/pi_fraction
  phase = EXP(-DCMPLX(0.D0,1.D0)*theta)
  WRITE(6,*) 'Complex scaling phase', phase
  READ(5,*)
  READ(5,*) type_of_pot
  WRITE(6,*) 'Potential:', type_of_pot
  READ(5,*)
  READ(5,*) type_of_interaction
  WRITE(6,*) 'Potential:', type_of_interaction

END SUBROUTINE effective_interaction_data
!
!     Reads in and allocates sp data for the nuclear case
! 
SUBROUTINE read_sp_data
  USE single_particle_orbits
  USE constants
  use hf_constants
  USE wave_functions 
  IMPLICIT NONE
  INTEGER :: i, j, i1, l1, j1, t1, n1, neutron_orbs_file, proton_orbs_file, &
       n_orbits, p_orbits, count, number_nholes, number_pholes
  CHARACTER(LEN= 10) :: model, status
  REAL(DP), DIMENSION(n_lab) :: u, s

  READ(5,*) ; READ(5,*) n_orbits, number_nholes
  READ(5,*) ; READ(5,*) p_orbits, number_pholes
  neutron_data%total_orbits = n_lab*n_orbits
  proton_data%total_orbits = n_lab*p_orbits
  !     Setup all possible orbit information
  all_orbit%total_orbits=neutron_data%total_orbits +    &
       proton_data%total_orbits
  !     Allocate space in heap for all single-particle data
  CALL allocate_sp_array(neutron_data,neutron_data%total_orbits) 
  CALL allocate_sp_array(proton_data,proton_data%total_orbits) 
  CALL allocate_sp_array(all_orbit,all_orbit%total_orbits) 
  !     Read neutron single-particle data
  READ(5,*);       READ(5,*);       READ(5,*)

  write(6,*) n_orbits, number_nholes
  write(6,*) p_orbits, number_pholes
  count = 0
  DO j=1, n_orbits
     READ(5,*) n1, l1, j1, t1, status, model
     
     DO i=1, n_lab
        count = count +1
        neutron_data%nn(count)=n1
        neutron_data%p(count)= k_mesh(i)    ! momenta in units of MeV
        neutron_data%w_p(count)= w_mesh(i)  ! weights in units of MeV
        neutron_data%ll(count)=l1
        neutron_data%jj(count)=j1
        neutron_data%itzp(count)=t1
        neutron_data%orbit_status(count)=status
        neutron_data%model_space(count)=model
        neutron_data%e(count)=k_mesh(i)*k_mesh(i)*0.5/p_mass(t1)

!        write(6,*) neutron_data%nn(count), neutron_data%ll(count), neutron_data%jj(count), neutron_data%orbit_status(count)
     ENDDO
  ENDDO
  !     Read proton single-particle data
  READ(5,*);       READ(5,*);       READ(5,*)
  count = 0
  DO j=1,  p_orbits
     READ(5,*) n1, l1, j1, t1, status, model
     
     DO i=1, n_lab
        count = count + 1
        proton_data%nn(count)=n1
        proton_data%ll(count)=l1
        proton_data%jj(count)=j1
        proton_data%itzp(count)=t1 
        proton_data%orbit_status(count)=status
        proton_data%model_space(count)=model
        proton_data%p(count)=k_mesh(i)      ! momenta in units of MeV
        proton_data%w_p(count)=w_mesh(i)      ! weights in units of MeV
        proton_data%e(count)=k_mesh(i)*k_mesh(i)*0.5/p_mass(t1)
!        write(6,*) proton_data%nn(count), proton_data%ll(count), proton_data%jj(count), proton_data%orbit_status(count)
     ENDDO
  ENDDO
  !     Neutrons are in the internal structure always even numbers 
  DO i=1, neutron_data%total_orbits
     all_orbit%p(i*2)=neutron_data%p(i)
     all_orbit%w_p(i*2)=neutron_data%w_p(i)
     all_orbit%ll(i*2)=neutron_data%ll(i)
     all_orbit%nn(i*2)=neutron_data%nn(i)
     all_orbit%jj(i*2)=neutron_data%jj(i)
     all_orbit%e(i*2)=neutron_data%e(i)
     all_orbit%itzp(i*2)=neutron_data%itzp(i)
     all_orbit%model_space(i*2)=neutron_data%model_space(i)
     all_orbit%orbit_status(i*2)=neutron_data%orbit_status(i)
     
     !all_orbit%p(i)=neutron_data%p(i)
     !all_orbit%w_p(i)=neutron_data%w_p(i)
     !all_orbit%ll(i)=neutron_data%ll(i)
     !all_orbit%nn(i)=neutron_data%nn(i)
     !all_orbit%jj(i)=neutron_data%jj(i)
     !all_orbit%e(i)=neutron_data%e(i)
     !all_orbit%itzp(i)=neutron_data%itzp(i)
     !all_orbit%model_space(i)=neutron_data%model_space(i)
     !all_orbit%orbit_status(i)=neutron_data%orbit_status(i)
     

  ENDDO
  !     protons are in the internal structure always odd numbers
  DO i=1, proton_data%total_orbits
     !i1 = i + neutron_data%total_orbits
     all_orbit%p(i*2-1)=proton_data%p(i)
     all_orbit%w_p(i*2-1)=proton_data%w_p(i)
     all_orbit%nn(i*2-1)=proton_data%nn(i)
     all_orbit%ll(i*2-1)=proton_data%ll(i)
     all_orbit%jj(i*2-1)=proton_data%jj(i)
     all_orbit%e(i*2-1)=proton_data%e(i)
     all_orbit%itzp(i*2-1)=proton_data%itzp(i)
     all_orbit%model_space(i*2-1)=proton_data%model_space(i)
     all_orbit%orbit_status(i*2-1)=proton_data%orbit_status(i)
     
     !all_orbit%p(i1)=proton_data%p(i)
     !all_orbit%w_p(i1)=proton_data%w_p(i)
     !all_orbit%nn(i1)=proton_data%nn(i)
     !all_orbit%ll(i1)=proton_data%ll(i)
     !all_orbit%jj(i1)=proton_data%jj(i)
     !all_orbit%e(i1)=proton_data%e(i)
     !all_orbit%itzp(i1)=proton_data%itzp(i)
     !all_orbit%model_space(i1)=proton_data%model_space(i)
     !all_orbit%orbit_status(i1)=proton_data%orbit_status(i)

  ENDDO
  
  DO i=1, all_orbit%total_orbits
     WRITE(6,*) all_orbit%nn(i),all_orbit%ll(i), &
          all_orbit%jj(i), all_orbit%itzp(i), all_orbit%orbit_status(i)
  ENDDO
 
 ! DO i=1, all_orbit%total_orbits
 !    WRITE(6,'(4I3,2X,A10,2X,A10,2X,E12.6,2X,E12.6)') all_orbit%nn(i),all_orbit%ll(i), &
 !         all_orbit%jj(i), all_orbit%itzp(i), &
 !         all_orbit%orbit_status(i), all_orbit%model_space(i), &
 !         all_orbit%p(i), all_orbit%e(i)
 ! ENDDO
  !   max sp kinetic energy given energy of last orbit, largest momentum 
  max_kinetic = 0.5_dp*(all_orbit%p(all_orbit%total_orbits)**2) ! not divided by mass
  max_momentum = all_orbit%p(all_orbit%total_orbits)
  ! read in nmax, lmin: lmax, jmin: jmax, tz_min: tz_max for holes 
  !READ(5,*);READ(5,*) nh_max, lh_min, lh_max, jh_min, jh_max, tzh_min, tzh_max
  

END SUBROUTINE read_sp_data
      
!
!                  SUBROUTINE to fix the hole cutoff depending
!                  on the oscillator parameter.    
!

SUBROUTINE setup_of_holecutoff
  USE constants
  USE single_particle_orbits
  USE wave_functions  
  IMPLICIT NONE
  INTEGER :: h,nh, lh, iq, number_of_iterations
  REAL(DP) :: qmin, qmax, ph,  ptest, sigma, norm
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: q_points, weight_q
  complex(DPC), allocatable, dimension(:) ::sum_norm
  complex(DPC) :: zh, sum_hf

  sigma = 1.0_dp; norm = 1.0_dp
  number_of_iterations = 0
  ALLOCATE ( sum_norm(all_orbit%total_orbits) )
  DO WHILE( (number_of_iterations < 10) .AND. (ABS(sigma) > 1E-5) )         
     qmin = 0.0_dp; qmax = hole_cutoff
     ALLOCATE ( q_points(int_points), weight_q(int_points) )
     CALL gauss_legendre(qmin,qmax,q_points,weight_q,int_points)
!     ptest = all_orbit%p(1)
     sum_norm  = 0.0_dp
     DO h=1, all_orbit%total_orbits
        IF (all_orbit%orbit_status(h) /= 'hole') CYCLE
        lh = all_orbit%ll(h) 
        nh = all_orbit%nn(h) 
!        IF (all_orbit%p(h) /= ptest) CYCLE   ! only one hole state, independent of p
        sum_hf = 0.
        DO iq=1, int_points
           ph = q_points(iq)
           zh = phase*ph*oscillator_length
           sum_hf = sum_hf+phase**3*weight_q(iq)*ph*ph*rnl(nh,lh,zh)*rnl(nh,lh,zh)&
                *oscillator_length**3
        ENDDO
        sum_norm(h) = sum_hf
        write(6,'(21H Norm of hole state #,I3,2X,E12.6)') h, dble(sum_norm(h))
     ENDDO
     sigma = 0.0_dp
     DO h = 1, all_orbit%total_orbits
        sigma = sigma +ABS( sum_norm(h)-norm)
     ENDDO
     sigma = sigma/all_orbit%total_orbits
     number_of_iterations = number_of_iterations +1 
     WRITE(6,*) 'Sigma for this iteration', sigma, number_of_iterations
     DEALLOCATE ( weight_q, q_points)
     IF (ABS(sigma) > 1E-5) THEN 
        hole_cutoff = 2.0_dp*hole_cutoff          
        int_points = int_points*2
     ENDIF
  ENDDO
  WRITE(6,'(43H New hole cutoff and # integration points: ,E12.6,2X,I3)') &
           hole_cutoff, int_points
  DEALLOCATE (sum_norm)
END SUBROUTINE setup_of_holecutoff











