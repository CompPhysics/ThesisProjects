
!
! setting up interaction in lab-system
!
subroutine vmtx_free(gfree, gmatrix_configs,j_lab, isospin_z)
  use constants
  use configurations
  use wave_functions
  use ang_mom_functions
  use single_particle_orbits
  implicit none
  TYPE (configuration_descriptor) , INTENT(IN)  :: gmatrix_configs
  integer, intent(in) :: j_lab, isospin_z
  complex*16 :: int,pa, pb,pc,pd, wa, wb, wc, wd, v_pot, vlab, opt_pot
  integer :: i, j, a, b, c, d, no_mom
  INTEGER :: ja, jb, jc, jd, la, lb, lc, ld
  complex*16, DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs), INTENT(INOUT) :: gfree
    
  gfree = 0.; 
  ! loops over configurations
  DO j=1,gmatrix_configs%number_confs
     c= gmatrix_configs%config_ab(j+j-1)
     d= gmatrix_configs%config_ab(j+j)
     lc = all_orbit%ll(c) ; jc = all_orbit%jj(c)
     pc = all_orbit%p(c)  ; wc = all_orbit%w_p(c) 
     ld = all_orbit%ll(d); jd = all_orbit%jj(d)
     pd = all_orbit%p(d); wd = all_orbit%w_p(d) 
     DO i=1, gmatrix_configs%number_confs
        a= gmatrix_configs%config_ab(i*2-1)
        b= gmatrix_configs%config_ab(i*2)
        la = all_orbit%ll(a) ;ja = all_orbit%jj(a)
        pa = all_orbit%p(a) ; wa = all_orbit%w_p(a) 
        lb = all_orbit%ll(b); jb = all_orbit%jj(b)
        pb = all_orbit%p(b); wb = all_orbit%w_p(b) 
        v_pot = 0.0D0 

        !if ( a /= b ) cycle
        !if ( c /= d ) cycle
        !if ( a /= c ) cycle

        !if ( abs(pb) /= abs(pd) ) cycle
        !if ( abs(pa) /= abs(pc) ) cycle
        !       get < (ab)JT_z | V_NN | (cd)JT_z > 
        !       Test on momenta relations
        IF ( (ABS(pa-pb) > abs(pc+pd)) .OR. (ABS(pc-pd) > abs(pa+pb)) ) CYCLE
        
        CALL v_labmomentum_space(v_pot,la,ja,pa,lb,jb,pb,lc,jc,pc,   &
             ld,jd,pd,isospin_z,j_lab)
        !       Then multiply with integration factor
        !       (p_a^2dp_a p_b^2dp_b p_c^2dp_c p_d^2dp_d)^1/2
        !v_pot = 1.
        gfree(i,j) = v_pot !*SQRT(wa*wb*wc*wd)
        !opt_pot = opt_pot + wb*pb**2*v_pot
        
        if ( j /= i ) then
           gfree(j,i) = gfree(i,j) 
        end if
        write(6,*) dble(pa),dble(pc), dble(v_pot)
     ENDDO
  ENDDO
end subroutine vmtx_free
