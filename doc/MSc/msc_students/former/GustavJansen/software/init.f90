SUBROUTINE read_data
  USE constants
  USE single_particle_orbits  
  USE stored_bare_interaction
  IMPLICIT NONE
  INTEGER :: i, n_particles

  CALL read_sp_file

  itzmin = -4; itzmax=4; j_lab_min=0
  j_lab_max = 0; nmax = 0; lmax = 0
  DO i = 1, all_orbit%total_orbits
    IF (all_orbit%orbit_status(i) == 'hole') THEN
       ! Find total mass of nucleus
       n_particles = all_orbit%jj(i) +1
       total_mass = total_mass + n_particles*sp_mass(all_orbit%ptype(i))
    ENDIF
    ! Find max j in lab, n and l in rel/cm
    IF ( all_orbit%jj(i) > j_lab_max) j_lab_max = all_orbit%jj(i)
    IF ( all_orbit%ll(i) > lmax) lmax = all_orbit%ll(i)
    IF ( all_orbit%nn(i) > nmax) nmax = all_orbit%nn(i)
  ENDDO

  CALL read_data_file
  IF ( nlmax > nlmax_model) THEN
    nlmax = 2*nmax
  ENDIF
  strange_min = 0; strange_max = 0
  IF (yn_elements > 0) strange_min = -1
  IF (yy_elements > 0) strange_min = -2

  ! Now print info to new files for the Hartree-Fock output,
  ! single-particle data and interaction data
  IF ( hf_iterations > 0 ) THEN
    CALL write_hf_sp_file
    CALL write_hf_data_file
  ENDIF

END SUBROUTINE read_data

SUBROUTINE read_sp_file
  USE constants
  USE single_particle_orbits  
  CHARACTER :: text
  INTEGER :: i,j

  READ(12,*)
  READ(12,'(A64,I11)') text, mass_nucleus
  READ(12,'(A54,E13.6,2X,E14.6)')  text, oscl, hbar_omega
  !Skip oscillator lengths for all particles
  READ(12,*);READ(12,*);READ(12,*);READ(12,*);READ(12,*);READ(12,*)
  READ(12,'(A43,2I12)')text, jmin, jmax

  READ(12,*); READ(12,*)
  READ(12,'(A56,I10)') text, nlmax
  READ(12,'(A56,I10)') text, nlmax_model
  READ(12,'(A39,I15)') text, all_orbit%total_orbits 
  READ(12,*) 
  CALL allocate_sp_array(all_orbit,all_orbit%total_orbits) 
  DO j = 1, all_orbit%total_orbits
     READ(12,*) text, i, all_orbit%ptype(j), all_orbit%strange(j), &
        all_orbit%nn(j), all_orbit%ll(j), &
        all_orbit%jj(j), all_orbit%itzp(j), all_orbit%nshell(j), all_orbit%e(j),&
        all_orbit%evalence(j), all_orbit%orbit_status(j), all_orbit%model_space(j),&
        all_orbit%included(j)
        all_orbit%e_original(j) = all_orbit%e(j)
  ENDDO
END SUBROUTINE read_sp_file

SUBROUTINE write_hf_sp_file
  USE constants
  USE single_particle_orbits  
  ! Now print info to new files for the Hartree-Fock output, single-particle data and interaction data
  WRITE(8,'(68H   ----> Oscillator parameters, Model space and single-particle data)')
  WRITE(8,'(65HMass number A of chosen nucleus (important for CoM corrections): ,I10)') mass_nucleus 
  WRITE(8,'(55HOscillator length and energy - proton/neutron average: ,E12.6,2X,E12.6)')  oscl, hbar_omega
  WRITE(8,*) 'Min and max value of partial wave ang. mom', jmin, jmax
  WRITE(8,*) 'Max value of l in lab frame: ', lmax
  WRITE(8,*) 'Max value of n in lab frame: ', nmax
  WRITE(8,*) 'Max value of 2*n + l+ cm 2*N +L for large space:', nlmax
  WRITE(8,*) 'Max value of 2*n + l+ cm 2*N +L for model space:', nlmax_model
  WRITE(8,*) 'Total number of single-particle orbits', all_orbit%total_orbits 
  WRITE(8,'(131HLegend:        id    str    n     l     2j   tz    2n+l  HO-energy     evalence     particle/hole  inside/outside   include/exclude)') 

END SUBROUTINE write_hf_sp_file

SUBROUTINE read_data_file
  USE constants
  USE stored_bare_interaction
  CHARACTER :: text
  INTEGER :: i
  READ(7,*)
  READ(7,'(A34, A20)') text, type_of_nn_pot
  READ(7,'(A34, A20)') text, type_of_yn_pot
  READ(7,'(A34, A20)') text, type_of_yy_pot
  READ(7,'(A20,A20)')text, type_of_renormv
  READ(7,'(A38,I4)') text, n_startenergy_g
  IF ( n_startenergy_g == 0) n_startenergy_g = 1
  ALLOCATE(e_start_g(n_startenergy_g))
  READ(7,*)  (e_start_g(i), i=1,n_startenergy_g )
  READ(7,'(A38,I15)') text, number_twobody_elements
  ! read NN elements
  READ(7,'(A33,4I15)') text, nn_elements, number_pp_elements, &
    number_pn_elements, number_nn_elements
  ! read YN elements
  READ(7,'(A35,5I15)') text, yn_elements, number_Spp_elements, &
    number_Lp_elements, number_Ln_elements, number_Smn_elements
  ! read YY elements
  READ(7,'(A35,6I15)') text, yy_elements, number_SpSp_elements, &
    number_LSp_elements, number_LS0_elements, number_LSm_elements, &
    number_SmSm_elements
  READ(7,*); READ(7,*)
  ! Then read the matrix elements
  CALL read_gmatrix_yy
  CALL read_gmatrix_yn
  CALL read_gmatrix_nn
END SUBROUTINE read_data_file

SUBROUTINE write_hf_data_file
  USE constants
  USE stored_bare_interaction
  INTEGER :: i
  WRITE(9,'(34H   ----> Interaction part         )')
  WRITE(9,'(34HNucleon-Nucleon interaction model:,A20)') type_of_nn_pot
  WRITE(9,'(34HHyperon-Nucleon interaction model:,A20)') type_of_yn_pot
  WRITE(9,'(34HHyperon-Hyperon interaction model:,A20)') type_of_yy_pot
  WRITE(9,'(20HType of calculation:,A20)') type_of_renormv
  WRITE(9,'(38HNumber and value of starting energies:,I4)') n_startenergy_g
  WRITE(9,'(10(1X,E12.6) )') (e_start_g(i), i=1,n_startenergy_g )
  WRITE(9,'(38HTotal number of twobody matx elements:,I15 )') number_twobody_elements
  WRITE(9,'(33HTotal number of NN matx elements:,4I15 )')nn_elements, number_pp_elements, &
    number_pn_elements, number_nn_elements
  WRITE(9,'(35HTotal number of YN matrix elements:,5I15 )') yn_elements, &
    number_Spp_elements, number_Lp_elements, number_Ln_elements, &
    number_Smn_elements
  WRITE(9,'(35HTotal number of YY matrix elements:,6I15 )') yy_elements, &
    number_SpSp_elements, number_LSp_elements, number_LS0_elements, &
    number_LSm_elements, number_SmSm_elements
  WRITE(9,'(91HMatrix elements with the following legend, NOTE no hbar_omega/A for Hcom, p_ip_j and r_ir_j)')
  WRITE(9,'(108H  Str Tz Par  2J   a   b   c   d          <ab|V|cd>         <ab|Hcom|cd>     <ab|r_ir_j|cd>   <ab|p_ip_j|cd>)')
END SUBROUTINE write_hf_data_file
!
!                  SUBROUTINE to fix the cutoff in mesh points
!                  for the harmonic oscillator wave functions
!                  which depend on the oscillator parameter.    
!
SUBROUTINE setup_ho_cutoff
  USE single_particle_orbits
  USE constants
  USE wave_functions  
  IMPLICIT NONE
  INTEGER :: h,nh, lh, iq, number_of_iterations, int_points
  REAL(DP) :: sigma, norm, oscl_r
  REAL(DP) :: qmin, qmax, sum_hf, sum_norm(0:lmax)
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: q_points, weight_q
  REAL(DP) :: z_rel, factor, xp, cx(0:200), contrib

  sigma = 1.0_dp; norm = 1.0_dp; int_points = 10
  number_of_iterations = 0
  nh = nmax; oscl_r=oscl*SQRT(2.)
  !  for momentum space if we set z= 0.5*(oscl_r*k)**2 > 60-70, then
  !  EXP(-60) < 10E-27. Making it twice as large ensures that we account
  !  for extensions due to the Laguerre polynoms which depend on z**(2n).
  cutoff = 2*SQRT(60.0_dp)/oscl
  qmin = 0.0D0; qmax = cutoff
  DO WHILE( (number_of_iterations < 20) .AND. (ABS(sigma) > 1E-4) )         
     ALLOCATE ( q_points(int_points), weight_q(int_points) )
     CALL gauss_legendre(qmin,qmax,q_points,weight_q,int_points)
     sum_norm  = 0.0D0
     DO h=0, lmax
        lh = h
        sum_hf = 0.
        factor = 0.5D0*((nh+lh+2)*LOG(2.D0)+fac(nh)-dfac(2*nh+2*lh+1)-0.5D0*LOG(pi))
        factor = EXP(factor)
        DO iq=1, int_points
           z_rel= q_points(iq)*q_points(iq)*oscl_r*oscl_r
           CALL laguerre_general( nh, lh+0.5D0, z_rel, cx )
           xp = EXP(-z_rel*0.5)*((q_points(iq)*oscl_r)**lh)*cx(nh)
           contrib = xp*factor*(oscl_r**(1.5D0)) 
           sum_hf=sum_hf+ weight_q(iq)*(contrib*q_points(iq))**2
        ENDDO
        sum_norm(h) = sum_hf
     ENDDO
     sigma = 0.D0
     DO h = 0, lmax
        sigma = sigma +ABS( sum_norm(h)-norm)
     ENDDO
     sigma = sigma/(lmax+1)
     number_of_iterations = number_of_iterations +1 
     !     WRITE(6,*) 'Sigma for this iteration', sigma, number_of_iterations
     DEALLOCATE ( weight_q, q_points)
     IF (ABS(sigma) > 1E-4) THEN 
        int_points = int_points+10
     ENDIF
  ENDDO
  n_rel = int_points
  n_cm_mom = int_points
  WRITE(6,'(41H New HO cutoff and # integration points: ,E12.6,2X,I5)') &
       cutoff, n_rel

END SUBROUTINE setup_ho_cutoff

