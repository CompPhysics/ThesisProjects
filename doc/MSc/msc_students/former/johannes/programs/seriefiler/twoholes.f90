PROGRAM twoholes
USE constants
USE configurations
USE ang_mom_functions
USE single_particle_orbits
USE wave_functions
USE partial_waves
USE relcm_gmatrix
IMPLICIT NONE

INTEGER :: t_z, parity, j_total, icount, lmax_lab,h1,h2
INTEGER ::  p, l, j, j_min, j_max 
REAL(DP) ::  cutoff, p1, fermi, relpot, relint 
INTEGER, DIMENSION(100000) :: la, ja 
REAL(DP), DIMENSION(100000) :: pa, wa  
REAL(DP), DIMENSION(120,120) ::  vkk0
INTEGER ::  number_meshpoint
REAL(DP), DIMENSION(10) :: mom,weight1, momh, weight2 
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
cutoff=2.5
k_cutoff=20
fermi=1.0
number_meshpoint=10
phase=1.
call gauss_legendre(p1,cutoff,mom,weight1,number_meshpoint) !cutoff istedet for fermi??
!sjekk om dette er riktige dimensjoner
allocate(krel(n_k),wkrel(n_k),t_mesh(n_k),wt_mesh(n_k))
call rel_mesh
no_channels=1
allocate(vfree(2*n_k,2*n_k,no_channels))
call gauss_legendre(p1,fermi,momh,weight2, number_meshpoint) !cutoff istedet for fermi??
!
!momh=momh*hc
!weight2=weight2*hc
!call calcvfree(vkk0,2*n_k,2*n_k,mom,number_meshpoint)
!write(*,*) vkk0


!Her settes potensialet opp
Do p=1,10
 DO h1=1,10
 vkk0(p,h1)=  -3.7368/(momh(p)*momh(h1))*log((0.49+(momh(p)+momh(h1))**2) &
/(0.49+(momh(p)-momh(h1))**2))-589.5/(momh(p)*momh(h1))*log((7.84+(momh(p)+momh(h1))**2)&
/(7.84+(momh(p)-momh(h1))**2))+1323.3265/(momh(p)*momh(h1))*log((24.01+(momh(p)&
+momh(h1))**2)&
/(24.01+(momh(p)-momh(h1))**2))
 END DO
END DO


!Algoritma for aa finne antall orbitaler
lmax_lab = 2
icount=0
DO p=1, number_meshpoint
 DO l=0, lmax_lab
  !passer på at to partikkel systemet ikke overstiger cutoff
  if(2.*mom(p) .gt. cutoff) CYCLE

  !To make it easier the total angular momentum is set to integers,  
  ! I multiply it with 2
  j_min= 2*l-1; j_max= 2*l+1

  IF(j_min < 0 ) j_min=1

  DO j=j_max,j_min,-2
   icount = icount+1
   ja(icount)=j !angulaer moment 
   la(icount)=l !orbital moment
   pa(icount)=mom(p)  !bev mengde
   wa(icount)=weight1(p) !vektene til gaussisk kvad
  END DO
 END DO
END DO
 
all_orbit%total_orbits= 2*icount !Hvorfor 2*icount ant orb pga prot og noy
write(*,*), 'total orbits = ', all_orbit%total_orbits
write(*,*) 'jmin', jmin, 'jmax', jmax

! Allokerer noen arrayer som brukes
call allocate_sp_array(all_orbit,all_orbit%total_orbits)
call cheb_to_set(n_k1, t_mesh,wt_mesh) !nk_1 istedet for 10000, t_mesh, 
!og til sist wt_mesh

write(*,*) 'no_channels', no_channels
!allocate(vfree(2*n_k,2*n_k,no_channels))
energy=0.
!t_z=-1
vfree(:,:,1)= 1.!  vkk0(:,:)!1. Gir verdier til potensialet vfree


relint=0.

! Her begynner jeg aa regne ut diagrammet
    DO h1=1,all_orbit%total_orbits/2
     DO h2=1,all_orbit%total_orbits/2
      if(pa(h1) .gt. fermi .or. pa(h2) .gt. fermi) CYCLE

        j_max= (ja(h1)+ja(h2))
        j_min=ABS(ja(h1)-ja(h2))

        DO j_total=j_min,j_max,2

          DO t_z= -1,1
         v_pot=1.

! Bruker vfree for aa omgjore vv i rel. koord til lab koord.
        call v_labmomentum_space(v_pot, la(h1),ja(h1),pa(h1),la(h2),ja(h2),&
       pa(h2), la(h1), ja(h1),pa(h1), la(h2), ja(h2), pa(h2), t_z, j_total)
       
         if(v_pot .eq. 0) CYCLE
         ! Regner ut energien
         energy=energy+(2*j_total+1)*v_pot*wa(h1)*wa(h2)*pa(h1)**2*pa(h2)**2
      
 !       write(*,*) energy
 !        relint=relint +(2*j_total+1)*vfree(h1,h2,1)*pa(h1)**2*pa(h2)**2*wa(h1)*wa(h2)   
       END DO
      END DO
    END DO
    END DO 
   write(*,*) 'sum av energi', energy/(pi**2)
!   write(*,*) relint

relint=0.
DO p=1,10
do h1=1,10
!if (mom(p) .gt. fermi .or. mom(h1) .gt. fermi) cycle
!if(mom(h1) .gt. fermi .or. mom(p) .gt. fermi) CYCLE
 relpot=  vfree(p,h1,1)*momh(p)**2*momh(h1)**2*weight2(p)*weight2(h1)
 relint=relint+relpot/(pi**2)
end do 
end do
write(*,*) 'integrasjon over potensialet', relint!, relint**(1./3)
!   call large_number_confs(j_total,parity,t_z,this)
call large_number_confs(j_total,parity,t_z,this)
deallocate(krel,wkrel,t_mesh,wt_mesh)


END PROGRAM twoholes


