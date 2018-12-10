SUBROUTINE onebody_contribution
  USE single_particle_orbits
  USE onebody_diagrams
  USE constants
  USE wave_functions
  IMPLICIT NONE
  REAL(DP), DIMENSION(n_startenergy_veff) :: onebody_diagram_1, &
       onebody_diagram_2,onebody_diagram_3
  REAL(DP) :: sum_kin, factor, one_a, tot_pot, tot_kin
  REAL(DP), DIMENSION(all_orbit%total_orbits,all_orbit%total_orbits) :: kinenergy
  INTEGER :: a, c, iph, k, id
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: int_factor, kin_energy

  tot_pot = 0.0d0
  tot_kin = 0.0d0
  !WRITE(*,*) n_startenergy_veff
  !WRITE(*,*) "onebody_contribution called..."
  ALLOCATE (int_factor(n_rel), kin_energy(n_rel) )
  int_factor(:)=wra(:)*ra(:)*ra(:)
  kin_energy(:)=hbarc*hbarc*ra(:)*ra(:)/2

  ALLOCATE (  one_body_terms(all_orbit%total_orbits, &
       all_orbit%total_orbits,n_startenergy_veff) )
       one_body_terms = 0.0d0
  !WRITE(*,*) "after alloc: ", one_body_terms(1,1,:)
  ALLOCATE (  one_body_folded &
       (all_orbit%total_orbits, &
       all_orbit%total_orbits) )
  kinenergy = 0.0_dp
  WRITE(6,*) 'a, c,  onebody part, kinetic energy and  total sp energy'
  DO a = 1, all_orbit%total_orbits
     IF (all_orbit%included(a) == 'exclude') CYCLE         
     IF (all_orbit%model_space(a) == 'outside') CYCLE         
     DO c= a, all_orbit%total_orbits
        !write(*,*) "a,c", a, c
        !        CALL one_body_diagram_phase(a,c,phase)
        IF (all_orbit%included(c) == 'exclude') CYCLE         
        IF (all_orbit%model_space(c) == 'outside') CYCLE         
        IF(all_orbit%ptype(a) /= all_orbit%ptype(c) ) CYCLE
        IF(iph(all_orbit%ll(a)) /= iph(all_orbit%ll(c))) CYCLE
        IF ( all_orbit%jj(a) /= all_orbit%jj(c)) CYCLE
        id = all_orbit%ptype(a)
        !write(*,*) "a,c", a, c, all_orbit%ptype(a)
        sum_kin=0.0_dp
        IF ( hf_iterations == 0 ) THEN
           DO k=1,n_rel
              sum_kin=sum_kin+hol(k,all_orbit%ll(a),all_orbit%nn(a),id)*hol(k,all_orbit%ll(c),all_orbit%nn(c),id)  &
                   *int_factor(k)*kin_energy(k)
           ENDDO
        ELSE
           DO k=1,n_rel
              sum_kin=sum_kin+wave_function(k,a)*wave_function(k,c) *&
              int_factor(k)*kin_energy(k)
           ENDDO
        ENDIF
        factor = (1.0/sp_mass(id) - 1.0/total_mass)
        !write(*,*) "Factor: ", factor
        sum_kin = sum_kin*factor
        kinenergy(c,a) = sum_kin
        onebody_diagram_1 = 0.0_dp
        SELECT CASE ( order_of_interaction )
        CASE ('first')
              !write(*,*) onebody_diagram_1
              CALL diagram_1(a,c,onebody_diagram_1)
        !write(*,*) "a,c after", a, c
              !write(*,*) onebody_diagram_1
              !write(*,*) one_body_terms(a,c,:)
              one_body_terms(a,c,:)=onebody_diagram_1(:)
        !write(*,*) "a,c after", a, c
        CASE ('second')
           CALL diagram_1(a,c,onebody_diagram_1)
           CALL diagram_2(a,c,onebody_diagram_2)
           CALL diagram_3(a,c,onebody_diagram_3)
           one_body_terms(a,c,:)= &
                (onebody_diagram_1(:)+ &
                onebody_diagram_2(:)+ &
                onebody_diagram_3(:))
        END SELECT
        one_body_terms(c,a,:)=one_body_terms(a,c,:)
        one_a = one_body_terms(a,c,n_startenergy_veff/2+1)
        WRITE(6,'(2I4,2X,3F12.6)') a , c, one_body_terms(a,c,n_startenergy_veff/2+1), &
                                   kinenergy(c,a),one_body_terms(a,c,n_startenergy_veff/2+1)+kinenergy(c,a)
        IF ((a==c) .AND. (all_orbit%orbit_status(a) == 'hole') ) THEN
            
            tot_pot = tot_pot + one_body_terms(a,c,n_startenergy_veff/2+1)*(all_orbit%jj(a) +1)
            tot_kin = tot_kin + kinenergy(c,a)*(all_orbit%jj(a) +1)
        ENDIF

     ENDDO
  ENDDO
  WRITE(6,*) "Total kinetic energy (hole): ", tot_kin
  WRITE(6,*) "Total potential energy (hole): ", tot_pot
  WRITE(6,*) "Total energy (hole): ", tot_kin+tot_pot

END SUBROUTINE onebody_contribution

SUBROUTINE diagram_1(a,c,onebody_diagram_1)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: j_min, j_max, jph, h, i
  REAL(DP), DIMENSION(n_startenergy_veff) :: w 
  REAL(DP) :: val, ang_mom_factor
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_1
  REAL(DP), DIMENSION(n_startenergy_g+3) :: ans

  !WRITE(*,*) "diagram_1 called: a,c", a, c
  onebody_diagram_1=0.0_dp
  !WRITE(*,*) e_start_g
  DO h=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE         
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     w=starting_energy*0.5_dp+all_orbit%e(h)+wcn
     DO jph=j_min,j_max
        ang_mom_factor=(2.*jph+1.)/(all_orbit%jj(a)+1.)
        CALL pphhmtx(a,h,c,h,jph,ans) ; IF ( ans(1) == 0. ) CYCLE
        !write(6,*) a,h,c,h,jph,ans(1), onebody_diagram_1(1)
        DO i=1, n_startenergy_veff
          
           !WRITE(*,*) "i, w", i, w(i)
           call interpolate(w(i),e_start_g,ans,val)
           !write(*,*) "i, val", i, val
           onebody_diagram_1(i)=onebody_diagram_1(i)+val*ang_mom_factor
        ENDDO
        write(6,'(5I4,2X,F6.2,3X,3F16.8)') a,h,c,h,jph,ang_mom_factor, ans(1), onebody_diagram_1(1),&
                ans(4)
     ENDDO
  ENDDO
END  SUBROUTINE diagram_1
!
!     2p1h diagram
!
SUBROUTINE diagram_2(a,c,onebody_diagram_2)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, h, p1, p2, i
  REAL(DP) :: val1, val2, factr
  REAL(DP), DIMENSION(n_startenergy_veff) :: w , de 
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_2
  REAL(DP), DIMENSION(n_startenergy_g+3) :: ans1, ans2

  onebody_diagram_2=0.
  DO h=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE         
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/(all_orbit%jj(a)+1.)
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) == 'hole') CYCLE         
           DO p2=1, all_orbit%total_orbits
              IF (all_orbit%orbit_status(p2) == 'hole') CYCLE         
              IF( ( all_orbit%model_space(h) == 'inside' ).AND.( all_orbit%model_space(p1) == 'inside' ) & 
                   .AND.( all_orbit%model_space(p2) == 'inside' ) ) CYCLE
              de=starting_energy*0.5_dp+all_orbit%e(h)-all_orbit%e(p1)-all_orbit%e(p2)+wcn
              CALL pphhmtx(a,h,p1,p2,jtot,ans1) ; IF ( ans1(1) == 0. ) CYCLE
              CALL pphhmtx(p1,p2,c,h,jtot,ans2) ; IF (ans2(1) == 0. ) CYCLE
              w=starting_energy*0.5+all_orbit%e(h)+wcn
              DO  i=1, n_startenergy_veff
                 CALL interpolate(w(i),e_start_g,ans1,val1)
                 CALL interpolate(w(i),e_start_g,ans2,val2)
                 onebody_diagram_2(i)=onebody_diagram_2(i)+ &
                      factr*0.5*val1*val2/de(i)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_2
!
!     2h1p diagram
!
SUBROUTINE diagram_3(a,c,onebody_diagram_3)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, p, i, h1, h2
  REAL(DP) :: val1, val2, factr
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2 , de 
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_3
  REAL(DP), DIMENSION(n_startenergy_g+3) :: ans1, ans2

  onebody_diagram_3=0.
  DO p=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p) == 'hole') CYCLE         
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(p))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(p))/2
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/(all_orbit%jj(a)+1.)
        DO h1=1,all_orbit%total_orbits
           IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE         
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE         
              IF( ( all_orbit%model_space(p) == 'inside' ).AND.( all_orbit%model_space(h1) == 'inside' ) & 
                   .AND.( all_orbit%model_space(h2) == 'inside' ) ) CYCLE
              de=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p)-starting_energy*0.5_dp+wcn
              CALL pphhmtx(h1,h2,c,p,jtot,ans1) ; IF ( ans1(1) == 0. ) CYCLE
              CALL pphhmtx(a,p,h1,h2,jtot,ans2) ; IF (ans2(1) == 0. ) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+ wcn 
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              DO i=1, n_startenergy_veff
                 CALL  interpolate(w1(i),e_start_g,ans1,val1)
                 CALL  interpolate(w2(i),e_start_g,ans2,val2)
                 onebody_diagram_3(i)=onebody_diagram_3(i)- & 
                      factr*0.5d0*val1*val2/de(i)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_3

