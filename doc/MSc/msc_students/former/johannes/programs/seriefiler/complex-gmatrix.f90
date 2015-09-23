 
!             Program block heff-gmatrix.f
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
!             E-MAIL:   morten.hjorth-jensen@fys.uio.no
!             LANGUAGE: F90/F95  
!             LAST UPGRADE : Oct 2004, Similarity transformed effective
!             two-body interaction for plane waves in lab system.
!             The similarity transformation is done in the lab system.
!
!

SUBROUTINE setup_vsim
  USE single_particle_orbits
  USE configurations
  USE constants
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: gmatrix_configs 
  !  TYPE (configuration_descriptor) :: 
  INTEGER ::  p_parity, j_lab, isospin_z, p_confs, pq_confs, mspace, ia, ib, i
  COMPLEX(DPC), ALLOCATABLE :: cvec_pp(:,:)
  COMPLEX(DPC), ALLOCATABLE :: cvec_qp(:,:)
  COMPLEX(DPC), ALLOCATABLE::  ceig_p(:) 
  COMPLEX(DPC), ALLOCATABLE::  heff(:,:), hamiltonian(:,:)
  COMPLEX(DPC), ALLOCATABLE::  kinetic_energy(:,:), gfree(:,:)

  
  !     loop over isospin projection
  DO isospin_z=  itzmin,itzmax 
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity= 0,1           
        !     loop over angular momenta
        DO j_lab=  j_lab_min,j_lab_max
           !     find all possible configurations, large space and model space            
           CALL large_number_confs &
                (j_lab,p_parity,isospin_z,gmatrix_configs)
           IF (gmatrix_configs%number_confs <= 0 ) CYCLE
           pq_confs=gmatrix_configs%number_confs
           ALLOCATE(gmatrix_configs%config_ab(pq_confs+pq_confs) )
           CALL number_configurations_model &
                (j_lab,p_parity,isospin_z,gmatrix_configs)
           p_confs = gmatrix_configs%model_space_confs 
           IF ( p_confs <= 0 ) CYCLE
           CALL setup_configurations &
                (j_lab,p_parity,isospin_z,gmatrix_configs)
           ALLOCATE(hamiltonian(pq_confs, pq_confs))
           ALLOCATE(kinetic_energy(pq_confs, pq_confs))
           ALLOCATE(gfree(pq_confs, pq_confs))
           hamiltonian = 0; kinetic_energy = 0.; gfree = 0.
           ! set up kinetic energy
           DO  mspace = 1, pq_confs
              ia= gmatrix_configs%config_ab(mspace*2-1)
              ib= gmatrix_configs%config_ab(mspace*2)
              !kinetic_energy(mspace,mspace) = 0.5d0*(all_orbit%p(ia)**2+all_orbit%p(ib)**2)/p_mass
              !write(6,*)dble( kinetic_energy(mspace,mspace)), aimag( kinetic_energy(mspace,mspace))
              !
           ENDDO
           !     setup the interaction to diagonalize in the big space
           !     Performs transformation from rel and cm coordinates to
           !     lab system. 
  
           write(6,*) gmatrix_configs%number_confs, all_orbit%total_orbits
           !CALL vmtx_free(gfree, gmatrix_configs,j_lab,isospin_z)
           hamiltonian = kinetic_energy+gfree
           
           
           ALLOCATE(cvec_pp(p_confs, p_confs))
           ALLOCATE(cvec_qp(pq_confs-p_confs,p_confs))
           ALLOCATE(ceig_p(p_confs)); ALLOCATE(heff(p_confs, p_confs))
           write(6,*) pq_confs
           call diag_exact( hamiltonian,cvec_pp, ceig_p, p_confs)
           do i = 1, p_confs 
           !   write(6,*) dble(ceig_p(i)), aimag(ceig_p(i))
           end do
           
           
!           heff = 0.
           !write(6,*) p_confs
           !      Get model space eigenvalues, eigenvectors of 
           !      model space and excluded space and
           !      model space configurations which match the corresponding ones 
           !      of the large space
!           CALL eigenvalues_large_maxvector(cvec_pp,cvec_qp,ceig_p,hamiltonian,pq_confs,&
!                p_confs,gmatrix_configs)
           !      setup 2p-effective interaction in P-space using           
           !      the Lee-Suzuki similarity transformation
!           CALL lee_suzuki( cvec_pp, cvec_qp, ceig_p, p_confs, pq_confs-p_confs, heff )
           !      subtract kinetic energy, obtain Veffective and print it
!           DO  mspace = 1, p_confs
!              heff(mspace,mspace) = heff(mspace,mspace)-kinetic_energy(mspace,mspace)
!           ENDDO
           !  print the interaction
!           CALL gprint(isospin_z,p_parity,j_lab,heff,gmatrix_configs)
           !      free space
           DEALLOCATE(cvec_pp); DEALLOCATE(cvec_qp); DEALLOCATE(ceig_p)
           DEALLOCATE(gmatrix_configs%config_ab)
           DEALLOCATE(kinetic_energy,gfree); DEALLOCATE(heff)
           DEALLOCATE(hamiltonian)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE setup_vsim
