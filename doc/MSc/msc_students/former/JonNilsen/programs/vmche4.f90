
! Unadorned VMC code for liquid helium


MODULE variational_mc
  DOUBLE PRECISION, PRIVATE :: density, correlation_coef
  INTEGER, PRIVATE :: dimension, number_particles, number_measurements, &
       number_cycles
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), PRIVATE :: r
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sq_r
  DOUBLE PRECISION, PRIVATE :: length, length_2, length_sq, c, length_4, ll4
  DOUBLE PRECISION, PRIVATE, PARAMETER :: a = 9.27647
  DOUBLE PRECISION, PRIVATE, PARAMETER :: b = 40.88

CONTAINS

  ! Monte Carlo sampling with the Metropolis algorithm  

  SUBROUTINE mc_sampling
    IMPLICIT NONE
    INTEGER ::  cycles, variate, accept, dim, i, j, l
    INTEGER :: idum
    DOUBLE PRECISION :: total_energy, total_energy2, r2, s2, variance,&
         average_energy, ran1, sum, delta_prob, local_energy, factor, &
         average_energy2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: r_old, r_new

    length = (number_particles/density)**(1/3.)
    length_2 = length*0.5; length_sq = length*length
    length_4 = length*0.25; ll4 = length_4*length
    idum=-1
    c = correlation_coef**5
    ! allocate matrices which contain the position of the particles  
    ALLOCATE(r(dimension, number_particles))
    ALLOCATE(r_new(dimension)); ALLOCATE(r_old(dimension))
    ALLOCATE(sq_r(dimension))
    sq_r = 0.; r_new = 0.; r_old = 0.; r = 0.
    ! setup initial position and variables
    DO i = 1,  number_particles
       DO j=1, dimension
          r(j,i) = length*ran1(idum)
       ENDDO
    ENDDO
    ! loop over possible measurements

    total_energy = 0; total_energy2 = 0.
    variance = 0.
    factor = 1./number_measurements/number_particles
    DO variate = 1, number_measurements
       average_energy = 0.; accept = 0; average_energy2 = 0.
       ! loop over monte carlo cycles 
       DO cycles = 1, number_cycles
          ! loop over particles
          DO i = 1,  number_particles
             DO j=1, dimension
                r_old(j)=r(j,i)
                r_new(j)=r_old(j) + length_4*(ran1(idum)-0.5)
                IF ( r_new(j) < 0. ) r_new(j) = r_new(j) + length
                IF ( r_new(j) > length ) r_new(j) = r_new(j) - length
             ENDDO
             sum = 0.
             DO l = 1,  number_particles
                IF ( l == i) CYCLE
                CALL sphere(r_old,l,r2)
                CALL sphere(r_new,l,s2) 
                IF ( r2 > ll4 )   r2 = ll4
                IF ( s2 > ll4 )   s2 = ll4
                sum = sum + 1./(r2*r2*sqrt(r2))-1./(s2*s2*sqrt(s2))
             ENDDO
             sum = sum*c
             delta_prob = EXP(sum)
             IF ( ran1(idum) <= delta_prob ) THEN
                r(:,i) = r_new(:)
                accept = accept +1 
             ENDIF
          ENDDO
          CALL energy_contrib(local_energy)
          average_energy = average_energy + local_energy
          average_energy2 = average_energy2 + local_energy*local_energy
       ENDDO    ! end of loop over variational  steps 
       WRITE(6,*) 'Measurement number: ', variate, &
       ' Accepted moves per measurement in %: ', &
         (100.*accept)/number_cycles/number_particles
       average_energy = average_energy/(number_cycles)
       average_energy2 = average_energy2/(number_cycles)
       total_energy = total_energy +  average_energy
       total_energy2 = total_energy2 +average_energy2
    ENDDO
!    total_energy = total_energy*factor
!    total_energy2 = total_energy2*factor/number_particles
    variance = (total_energy2-total_energy*total_energy)*factor
    WRITE(6,*) 'Density', density
    WRITE(6,*) 'Total number of particle', number_particles
    WRITE(6,*) 'Total energy per particle', total_energy*factor
    WRITE(6,*) 'Variance', variance
    WRITE(6,*) 'Standard deviation', SQRT(variance/(number_measurements-1.))
    DEALLOCATE ( r_old, r_new);  DEALLOCATE ( r ) 
  END SUBROUTINE  mc_sampling  ! end mc_sampling function  

  SUBROUTINE sphere(position,l,radius2)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: position
    DOUBLE PRECISION, INTENT(INOUT) :: radius2
    INTEGER, INTENT(IN):: l

    radius2 = 0.
    sq_r(:) = ABS(position(:)-r(:,l))
    WHERE (sq_r > length_2 ) sq_r = length - sq_r  
    radius2 = SUM(sq_r*sq_r)

  END SUBROUTINE sphere

  ! Function to calculate the local energy 

  SUBROUTINE energy_contrib(local_energy)
  IMPLICIT NONE
  INTEGER :: i, j , k1, k2, k3
  DOUBLE PRECISION :: r2, r6, sum, z1, z2, z3, kinetic
  DOUBLE PRECISION, INTENT(OUT) :: local_energy
  
  kinetic = a*c
  local_energy = 0.
  DO i = 1, number_particles-1
     DO j = i+1, number_particles
        CALL sphere(r(:,j),i,r2)      
        ! correlation energy cutoff
        IF ( r2 <= ll4 ) THEN
           local_energy = local_energy+kinetic/(r2*r2*r2*SQRT(r2))
        ENDIF
        r2 = 0.; sum = 0.
        DO k1 = -1,  0, 1
           z1 = (sq_r(1)+k1*length)**2
           DO k2 = -1,  0, 1
              z2 = (sq_r(2)+k2*length)**2
              DO k3 = -1,  0, 1
                 z3 = (sq_r(3)+k3*length)**2
                 r2 = z1+z2+z3
                 ! potential cutoff
                 IF ( r2 > length_sq) CYCLE
                 r6 = 1./(r2*r2*r2)
                 sum = sum +(r6-1)*r6
              ENDDO
           ENDDO
        ENDDO
        local_energy = local_energy+sum*b

     ENDDO
  ENDDO

END SUBROUTINE energy_contrib


SUBROUTINE initialise
  IMPLICIT NONE
  WRITE(*,*)'number of particles = '
  READ(*,*) number_particles
  WRITE(*,*)'dimensionality = '
  READ(*,*) dimension
  WRITE(*,*)'maximum number of independent measurements, min 2 = '
  READ(*,*) number_measurements
  WRITE(*,*)'# Correlation coefficients= '
  READ(*,*) correlation_coef
  WRITE(*,*)'# MC steps= '
  READ(*,*)number_cycles
  WRITE(*,*)'# Total density= '
  READ(*,*)density
END SUBROUTINE initialise

END MODULE variational_mc

! Begin of main program   

PROGRAM vmc_helium
USE variational_mc
IMPLICIT NONE

!   Read in data 
    CALL initialise    
!   Do the mc sampling  
    CALL mc_sampling
END PROGRAM vmc_helium





