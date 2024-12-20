
PROGRAM twoholes
USE constants
USE configurations
USE ang_mom_functions
USE single_particle_orbits
USE wave_functions
USE partial_waves
USE relcm_gmatrix
IMPLICIT NONE

INTEGER :: t_z, parity, j_total, icount, lmax_lab
INTEGER ::  p, l, j, j_min, j_max, l1, l2, j1, j2,p2 
INTEGER :: j1_min, j2_min, j1_max, j2_max
REAL(DP) ::  cutoff, p1
INTEGER, DIMENSION(100000) :: la, ja, pa
INTEGER ::  number_meshpoint
REAL(DP), DIMENSION(10) :: mom,weight
COMPLEX(DPC) :: v_pot,energy
REAL(DP), DIMENSION(30) :: vekt_rel

TYPE (configuration_descriptor)  :: this
TYPE(single_particle_descript)   :: this_array

call commons_to_angmom
call setup_channels ! må kalle på disse to for å hente noen arrayer
int_points=10 !sjekk hva denne skal være
n_k=60         
n_k1=30       
n_k2=30
p1=0.
cutoff=2.5
k_cutoff=20
number_meshpoint=10
phase=1.
call gauss_legendre(p1,cutoff,mom,weight,number_meshpoint)
!write(*,*) mom(4)
!sjekk om dette er riktige dimensjoner
allocate(krel(n_k1),wkrel(n_k1),t_mesh(n_k1),wt_mesh(n_k1))
call rel_mesh

lmax_lab = 6
icount=0
DO p=1, number_meshpoint
 DO l=0, lmax_lab
!  IF(mom(p)**2/2. .le. cutoff) CYCLE
  !To make it easier the total angular momentum is set to integers,  
  ! I multiply it with 2
  j_min= 2*l-1; j_max= 2*l+1
  IF(j_min < 0 ) j_min=1
  DO j=j_max,j_min,-2
   icount = icount+1
 !write(*,*) j, icount
   ja(icount)=j
   la(icount)=l
   pa(icount)=p
  END DO
 END DO
END DO
 
all_orbit%total_orbits= 2*icount
write(*,*), 'total orbits = ', all_orbit%total_orbits
call allocate_sp_array(all_orbit,all_orbit%total_orbits)
call cheb_to_set(n_k1, t_mesh,wt_mesh)
!Det virker som t_mesh og wt_mesh er de som gir feil
!write(*,*) 't_mesh(3)o g wt_mesh(3)', t_mesh(3), wt_mesh(3)

energy=0.

DO j_total = 0, 2
 DO t_z= -1,1
  DO parity = -1, 1, 2
   call large_number_confs(j_total,parity,t_z,this)
    DO  p=1, number_meshpoint !loop over hole 1
     DO l1=0,6 !Er vel egentlig for stor l når j_total ikke er større enn 2
      j1_min= 2*l1-1; j1_max= 2*l1+1
      IF(j1_min < 0 ) j1_min=1
      DO j1=j1_max,j1_min,-2
       DO p2=1,number_meshpoint !loop over hole 2
        DO l2=0,6
         j2_min= 2*l2-1; j2_max= 2*l2+1
         IF(j2_min < 0 ) j2_min=1
         DO j2=j2_max,j2_min,-2
       
          jmin=ABS(j2-j1); jmax=j1+j2 
          call v_labmomentum_space(v_pot, l1,j1,mom(p),l2,j2, mom(p2), &
          l1, j1, mom(p), l2, j2, mom(p2), t_z, j_total)
          energy=energy+(2*j_total+1)*v_pot
        !  write(*,*) 'energy', energy
 
         
         END DO
        END DO
       END DO
      END DO
     END DO
    END DO
  END DO
 END DO
END DO
write(*,*) v_pot
!call large_number_confs(j_total,parity,t_z,this)
!deallocate(krel,wkrel,t_mesh,wt_mesh)


END PROGRAM twoholes
