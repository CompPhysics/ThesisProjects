
PROGRAM SECONDORDER
USE constants
USE configurations
USE ang_mom_functions
USE single_particle_orbits
USE wave_functions
USE partial_waves
USE relcm_gmatrix
IMPLICIT NONE

INTEGER :: t_z, parity, j_total, icount, lmax_lab
INTEGER ::  p, l, j, j_min, j_max, l1, l2, j1, j2, h1, h2, p1, p2
INTEGER :: j1_min, j2_min, j1_max, j2_max, jc_min,jc_max, jd_min, jd_max
INTEGER :: pc, pd
REAL(DP) ::  cutoff, p1, fermi
INTEGER, DIMENSION(100000) :: la, ja, pa, wa, pp, wp
INTEGER ::  number_meshpoint
REAL(DP), DIMENSION(10) :: mom, weight, mom_part, weight_part
COMPLEX(DPC) :: v_pot,energy
!REAL(DP), DIMENSION(30) :: vekt_rel

TYPE (configuration_descriptor)  :: this
!TYPE(single_particle_descript)   :: this_array

call commons_to_angmom
call setup_channels ! må kalle på disse to for å hente noen arrayer
int_points=10 !sjekk hva denne skal være
n_k=60         
n_k1=30       
n_k2=30
p1=0.
fermi = 1.3
cutoff=2.5
k_cutoff=20
number_meshpoint=10
phase=1.
call gauss_legendre(p1,fermi,mom,weight,number_meshpoint)
call gauss_legendre(fermi,cutoff,mom_part,weight_part,number_meshpoint)
!sjekk om dette er riktige dimensjoner
allocate(krel(n_k),wkrel(n_k),t_mesh(n_k),wt_mesh(n_k))
call rel_mesh

!write(*,*) 'n_k, n_k1, n_k2', n_k, n_k1, n_k2

lmax_lab = 2
icount=0
DO p=1, number_meshpoint
 DO l=0, lmax_lab
!  IF(mom(p)*2. .gt. cutoff) CYCLE
  !To make it easier the total angular momentum is set to integers,  
  ! I multiply it with 2
  j_min= 2*l-1; j_max= 2*l+1
  IF(j_min < 0 ) j_min=1
  DO j=j_max,j_min,-2
   icount = icount+1
   ja(icount)=j
   la(icount)=l
   pa(icount)=mom(p)
   wa(icount)=weight(p)
   pp(icount)=mom_part(p)
   wp(icount)=weight_part(p)
  END DO
 END DO
END DO
 
all_orbit%total_orbits= 2*icount
write(*,*), 'total orbits = ', all_orbit%total_orbits
call allocate_sp_array(all_orbit,all_orbit%total_orbits)
call cheb_to_set(n_k1, t_mesh,wt_mesh)

write(*,*) 'no_channels', no_channels
allocate(vfree(2*n_k,2*n_k,no_channels))
energy=0.

vfree(:,:,:)=1.


 DO h1=1, all_orbit%total_orbits
  DO h2=1, all_orbit%total_orbits
   DO p1=1,all_orbit%total_orbits
    Do p2=1,all_orbit%total_orbits

      j_max=(ja(h1)+ja(h2))
      j_min=ABS(ja(h1)-ja(h2))     

      DO j_total=j_min,j_max,2    
       DO t_z= -1,1

  call v_labmomentum_space(v_pot,la(p1),ja(p1),pp(p1),la(p2),ja(p2), pp(p2), &
          la(h1), ja(h1), pa(h1), la(h2), ja(h2), pa(h2), t_z, j_total)

          energy=energy+((2*j_total+1)*v_pot*mom_part(pc)**2*mom_part(pd)**2&
*mom(p)**2*mom(p2)**2*weight(p)*weight(p2)*weight_part(pc)*weight_part(pd))**2&
 /(0.5*(mom_part(pc)**2+mom_part(pd)**2-mom(p)**2-mom(p2)**2 ) )
          write(*,*) 'energy', energy

      END DO   !slutt paa t_z
    END DO  !slutt p2
   END DO !slutt p1
  END DO !slutt h2
 END DO ! slutt h1
write(*,*) 'v_pot', v_pot
!call large_number_confs(j_total,parity,t_z,this)
!deallocate(krel,wkrel,t_mesh,wt_mesh)


END PROGRAM SECONDORDER
