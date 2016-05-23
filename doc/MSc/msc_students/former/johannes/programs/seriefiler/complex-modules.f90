!             Program block complex-modules.f90   
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
!             E-MAIL:   morten.hjorth-jensen@fys.uio.no
!             LANGUAGE: F90/F95
!             LAST UPGRADE : February  2005
!             This program block contains the definiton of
!             all modules used by the various program blocks.
!
!    This module contains all constants and declarations 
!    of variables read in by the function read_data. These
!    variables are used by many functions.

MODULE constants
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
  ! min and max isospin projection
  INTEGER, PUBLIC :: itzmin, itzmax
  ! min and max total two-body angular momentum in lab frame
  INTEGER, PUBLIC :: j_lab_min, j_lab_max
  ! min and max total two-body angular momentum in Rel-CoM frame
  INTEGER, PUBLIC :: jmin, jmax
  ! lmax is max relative l
  INTEGER, PUBLIC :: lmax
  ! number of integration points for lab, relative and CoM momenta
  INTEGER , PUBLIC :: n_k, int_points, n_lab,n_lab1,n_lab2,n_lab3, n_k1, n_k2, nrad
  ! A of closed shell core
  INTEGER, PUBLIC :: mass_closed_shell_core
  ! max number of iteration in hartree-fock potential
  INTEGER, PUBLIC :: hf_iterations
  ! the nn potential, cd-bonn, idaho-a, idaho-b etc
  CHARACTER (LEN=100), PUBLIC  :: type_of_pot, type_of_interaction, type_of_calculation
  COMPLEX(DPC), PUBLIC :: phase
  REAL(DP), PUBLIC :: theta, k_cutoff, k_max, lambda
  REAL(DP), PUBLIC :: max_kinetic, max_momentum
  REAL(DP), PUBLIC :: oscillator_length, hole_cutoff, oscillator_energy
  REAL(DP), PUBLIC, DIMENSION(-1:1), PARAMETER:: p_mass = &
            (/938.27231_dp, 938.91897_dp, 939.56563_dp/)
  REAL(DP), PARAMETER, PUBLIC :: hbarc = 197.327053_dp
  REAL(DP), PUBLIC, PARAMETER :: pi = 3.141592741012573_dp
  REAL(DP), PUBLIC, PARAMETER :: pi_2 = 1.570796370506287_dp
  REAL(DP), PUBLIC, PARAMETER :: pi_4 = 0.7853981852531433_dp
  REAL(DP), PUBLIC, PARAMETER :: epsilon = 1.E-15
  DOUBLE PRECISION, PARAMETER :: M_p = 938.27231d0
  DOUBLE PRECISION, PARAMETER :: M_n = 939.56563d0
  DOUBLE PRECISION, PARAMETER :: M_np = 938.926d0
  DOUBLE PRECISION, PARAMETER :: redm_he5 = M_n*(2.*M_n + 2.*M_p)/&
       (3.*M_n + 2.*M_p )
 

END MODULE constants

MODULE HF_constants
  use constants
  INTEGER :: nh_max, lh_min, lh_max, jh_min, jh_max, itzh_min, itzh_max
  COMPLEX(DPC), ALLOCATABLE, PUBLIC :: v_bhf(:,:,:,:,:)
  COMPLEX(DPC), ALLOCATABLE, PUBLIC :: bhf_hol(:,:,:,:)
  COMPLEX(DPC), ALLOCATABLE, PUBLIC :: bhf_hol_temp(:,:,:,:)
  COMPLEX(DPC), ALLOCATABLE, PUBLIC :: bhf_hol_rspace_temp(:,:,:,:)
  COMPLEX(DPC), ALLOCATABLE, PUBLIC :: bhf_hol_rspace(:,:,:,:)
  complex(dpc),allocatable,  public :: ws_basis(:,:,:)
  double precision,allocatable,  public :: k_ws(:)
  double precision :: k_interpol_min, k_interpol_max
  real, allocatable, public :: vlab_mom_re(:),vlab_mom_im(:)
  
END MODULE HF_constants
  
!     Definition of single particle data

MODULE single_particle_orbits
  USE constants
  TYPE, PUBLIC :: single_particle_descript
     INTEGER :: total_orbits
     INTEGER, DIMENSION(:), POINTER :: nn, ll, jj, itzp, nshell
     CHARACTER (LEN=10), DIMENSION(:), POINTER :: orbit_status, model_space
     REAL(DP), DIMENSION(:), POINTER :: e, p, w_p
  END TYPE single_particle_descript
  TYPE, PUBLIC :: rel_cm_data
     INTEGER, DIMENSION(:), POINTER :: lrel
     INTEGER :: max_value
  END TYPE rel_cm_data
  TYPE (rel_cm_data), PUBLIC :: relcm_sp_data
  TYPE (single_particle_descript), PUBLIC :: all_orbit, neutron_data, &
       proton_data
CONTAINS
  SUBROUTINE allocate_sp_array(this_array,n)
    TYPE (single_particle_descript), INTENT(INOUT) :: this_array
    INTEGER , INTENT(IN) :: n
    IF (ASSOCIATED (this_array%nn) ) DEALLOCATE(this_array%nn)
    ALLOCATE(this_array%nn(n))
    IF (ASSOCIATED (this_array%ll) ) DEALLOCATE(this_array%ll)
    ALLOCATE(this_array%ll(n))
    IF (ASSOCIATED (this_array%jj) ) DEALLOCATE(this_array%jj)
    ALLOCATE(this_array%jj(n))
    IF (ASSOCIATED (this_array%itzp) ) DEALLOCATE(this_array%itzp)
    ALLOCATE(this_array%itzp(n))
    IF (ASSOCIATED (this_array%e) ) DEALLOCATE(this_array%e)
    ALLOCATE(this_array%e(n))
    IF (ASSOCIATED (this_array%p) ) DEALLOCATE(this_array%p)
    ALLOCATE(this_array%p(n))
    IF (ASSOCIATED (this_array%w_p) ) DEALLOCATE(this_array%w_p)
    ALLOCATE(this_array%w_p(n))
    IF (ASSOCIATED (this_array%nshell) ) DEALLOCATE(this_array%nshell)
    ALLOCATE(this_array%nshell(n))
    IF (ASSOCIATED (this_array%orbit_status) ) DEALLOCATE(this_array%orbit_status)
    ALLOCATE(this_array%orbit_status(n))
    IF (ASSOCIATED (this_array%model_space) ) DEALLOCATE(this_array%model_space)
    ALLOCATE(this_array%model_space(n))
    !           blank all characters and zero all other values
    DO i= 1, n
       this_array%model_space(i)= ' '
       this_array%orbit_status(i)= ' '
       this_array%e(i)=0.
       this_array%p(i)=0.
       this_array%w_p(i)=0.
       this_array%nn(i)=0
       this_array%ll(i)=0
       this_array%jj(i)=0
       this_array%nshell(i)=0
       this_array%itzp(i)=0
    ENDDO

  END SUBROUTINE allocate_sp_array

  SUBROUTINE deallocate_sp_array(this_array)
    TYPE (single_particle_descript), INTENT(INOUT) :: this_array
    DEALLOCATE(this_array%nn) ; DEALLOCATE(this_array%ll)
    DEALLOCATE(this_array%jj) ;DEALLOCATE(this_array%itzp)
    DEALLOCATE(this_array%e) ;DEALLOCATE(this_array%nshell)
    DEALLOCATE(this_array%orbit_status); DEALLOCATE(this_array%model_space)
    DEALLOCATE(this_array%p);     DEALLOCATE(this_array%w_p)
  END SUBROUTINE deallocate_sp_array

  SUBROUTINE allocate_relcm_array(this_array,n)
    TYPE (rel_cm_data), INTENT(INOUT) :: this_array
    INTEGER , INTENT(IN) :: n
    IF (ASSOCIATED (this_array%lrel) ) DEALLOCATE(this_array%lrel)
    ALLOCATE(this_array%lrel(n))
    !           blank all values
    DO i= 1, n
       this_array%lrel(i)=0
    ENDDO

  END SUBROUTINE allocate_relcm_array

  SUBROUTINE deallocate_relcm_array(this_array)
    TYPE (rel_cm_data), INTENT(INOUT) :: this_array
    DEALLOCATE(this_array%lrel)

  END SUBROUTINE deallocate_relcm_array


END MODULE single_particle_orbits


!     the bare interaction-matrix in the rel and c.m. frame
MODULE relcm_gmatrix
  USE constants
  COMPLEX(DPC), ALLOCATABLE, PUBLIC :: vfree(:,:,:)

END MODULE relcm_gmatrix


!     setup all partial wave  combinations from jmin to jmax
MODULE partial_waves
  USE constants
  INTEGER, PUBLIC :: no_channels, nchans
  PARAMETER ( nchans=1000)
  INTEGER,PUBLIC :: orb_lrel_min(nchans),spin_rel(nchans), &
       jang_rel(nchans),iso(nchans), orb_lrel_max(nchans)
  CONTAINS

  !
  !          Setup partial waves.
    !          The max value of the partial wave J is def
  !          in the input file and transferred as
  !          jmin and jmax
  !
  SUBROUTINE setup_channels
    IMPLICIT NONE
    INTEGER :: j_ang, lorb_min, lorb_max, it, i_spin,l_orb, loop
    LOGICAL :: triag
    no_channels=0

    DO j_ang=jmin,jmax
       DO i_spin=0,1
          lorb_min=j_ang ; lorb_max=j_ang+i_spin
          DO l_orb=lorb_min,lorb_max
             DO it=-1, 1                 
                IF(.NOT.triag(j_ang,l_orb,i_spin)) THEN
                   IF ( ABS( it ) == 1) THEN
                      !     Pauli principle test for identical particles, pp or nn
                      IF(MOD(l_orb+1+i_spin,2) /= 0) THEN
                         IF ( l_orb == j_ang ) THEN
                            no_channels=no_channels+1
                            orb_lrel_max(no_channels)=l_orb
                            orb_lrel_min(no_channels)=l_orb
                         ELSEIF( l_orb /= j_ang ) THEN
                            IF ( j_ang > 0 ) THEN
                               no_channels=no_channels+1
                               orb_lrel_max(no_channels)=j_ang+i_spin
                               orb_lrel_min(no_channels)=j_ang-i_spin
                            ELSEIF ( j_ang == 0) THEN
                               no_channels=no_channels+1 
                               orb_lrel_max(no_channels)=j_ang+i_spin
                               orb_lrel_min(no_channels)=j_ang+i_spin
                            ENDIF
                         ENDIF
                         spin_rel(no_channels)=i_spin
                         jang_rel(no_channels)=j_ang
                         iso(no_channels)=it
                      ENDIF
                   ELSEIF (it == 0) THEN
                      !     For T_z=0 all possible waves are included, no restrictions
                      IF ( l_orb == j_ang ) THEN
                         no_channels=no_channels+1
                         orb_lrel_max(no_channels)=l_orb
                         orb_lrel_min(no_channels)=l_orb
                      ELSEIF( l_orb /= j_ang ) THEN
                         IF ( j_ang > 0 ) THEN
                            no_channels=no_channels+1
                            orb_lrel_max(no_channels)=j_ang+i_spin
                            orb_lrel_min(no_channels)=j_ang-i_spin
                         ELSEIF ( j_ang == 0) THEN
                            no_channels=no_channels+1
                            orb_lrel_max(no_channels)=j_ang+i_spin
                            orb_lrel_min(no_channels)=j_ang+i_spin
                         ENDIF
                      ENDIF
                      spin_rel(no_channels)=i_spin
                      jang_rel(no_channels)=j_ang
                      iso(no_channels)=it
                   ENDIF
                ENDIF

             ENDDO
          ENDDO
       ENDDO
    ENDDO
!   Write out all channels
    DO loop = 1, no_channels
       WRITE(6,'(12H Channel Nr:,I3,7H l_min:,I3,7H l_max:,I3,3H J:,I3,3H S:,I3,4H Tz:,I3)') &
            loop, orb_lrel_min(loop), orb_lrel_max(loop), jang_rel(loop), spin_rel(loop), iso(loop)
    ENDDO

  END SUBROUTINE setup_channels

END MODULE partial_waves

!           
!     This module contains the angular momentun functions
!     and transformation coefficients when going from
!     lab system  <--> cm system
!
MODULE ang_mom_functions
  USE constants
  REAL(DP), PRIVATE :: f_mb(50),g_mb(50),w_mb(50)
  INTEGER, PRIVATE :: kh(200)
  REAL(DP), PRIVATE :: q(50,50), cn(0:51,0:51) !PRIVATE vs PUBLIC

CONTAINS
  !
  !     factorials for 3j,6j and 9j symbols            
  !     for moshinsky trans brackets and for           
  !     vector brackets                                
  !
  SUBROUTINE commons_to_angmom
    IMPLICIT NONE
    INTEGER :: l, k, i, j
    REAL(DP) :: a , sq_pi, fj, tfj, fk
    !    3j, 6j and 9j symbols
    kh=1
    kh(100) =0
    DO l=1,50
       q(l,1)=1.0_dp
       q(l,l)=1.0_dp
       kh(l+l+100)=0
    ENDDO
    DO l=2,49
       DO k=2,l
          q(l+1,k)=q(l,k-1)+q(l,k)
       ENDDO
    ENDDO
    !    Moshinsky brackets
    f_mb(1)=0.0_dp
    g_mb(1)=LOG(0.5_dp)
    w_mb(1)=0.0_dp
    DO i=2,50
       a=i-1
       f_mb(i)=f_mb(i-1)+LOG(a)
       g_mb(i)=g_mb(i-1)+LOG(a+0.5_dp)
       w_mb(i)=LOG(a+a+1.)
    ENDDO
    !    spherical harmonics
    cn=0.
    sq_pi=1./SQRT(2.0_dp*pi)
    DO j=0,51
       cn(0,j)=SQRT(0.5_dp*(2.0_dp*j+1.))
    ENDDO
    DO j=1,51
       tfj=2.*j
       cn(j,j)=cn(j-1,j-1)*SQRT((tfj+1.)/tfj)
    ENDDO
    DO j=0,51
       fj=FLOAT(j)
       DO k=1,j-1
          fk=FLOAT(k)
          cn(k,j)=cn(k-1,j)*SQRT((fj+fk)*(fj-fk+1.))*0.5_dp/fk
       ENDDO
    ENDDO
    cn=cn*sq_pi

  END SUBROUTINE commons_to_angmom
  !
  !     calculates 3j-symbols           
  !
  REAL(DP) FUNCTION tjs(j_a,j_b,j_c,m_a,m_b,m_c)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j_a,j_b,j_c,m_a,m_b,m_c
    INTEGER :: ja, jb, jc, mb, ma, mc, la, lb, lc, lt, ld, ja2, jb2, &
         jc2, i, k0, k1, k, ip
    REAL(DP) :: eps, x, fn, p

    eps=1.0d-2
    tjs=0.
    ja=(j_a+m_a)/2+1
    ma=(j_a-m_a)/2+1
    jb=(j_b+m_b)/2+1
    mb=(j_b-m_b)/2+1
    jc=(j_c+m_c)/2+1
    mc=(j_c-m_c)/2+1
    la=(j_b+j_c-j_a)/2+1
    lb=(j_c+j_a-j_b)/2+1
    lc=(j_a+j_b-j_c)/2+1
    lt=(j_a+j_b+j_c)/2+1
    ld=MIN(ja,jb,jc,ma,mb,mc,la,lb,lc)
    IF(((m_a+m_b+m_c) <= 0).AND.(ld > 0)) THEN
       ja2=j_a+m_a
       jb2=j_b+m_b
       jc2=j_c+m_c
       i=ja2+jb2+jc2-ja2/2*2-jb2/2*2-jc2/2*2
       IF(i == 0) then
          fn=q(ja+ma-1,lc)*q(jb+mb-1,lc)/(q(lt,jc+mc-1)*q(lt+1,2) &
               *q(ja+ma-1,ja)*q(jb+mb-1,jb)*q(jc+mc-1,jc))
          k0=MAX(0,lc-ja,lc-mb)+1
          k1=MIN(lc,ma,jb)
          x=0.
          DO k=k0,k1
             x=-x-q(lc,k)*q(lb,ma-k+1)*q(la,jb-k+1)
          ENDDO
          ip=k1+lb+jc
          p=1-2*(ip-ip/2*2)
          tjs=p*x*SQRT(fn)
       ENDIF
    ENDIF

  END FUNCTION tjs
  !
  !     calculates 6j-symbols           
  !
  REAL(DP) FUNCTION sjs(j_a,j_b,j_c,l_a,l_b,l_c)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j_a,j_b,j_c,l_a,l_b,l_c
    INTEGER :: ja,jb,jc,la,lb,lc,i,mt,ma,mb,mc,na,nb,nc,ka,&
         kb,kc,l,l0,l1
    REAL(DP) :: x, fs, fss

    sjs=0.0d0
    ja=j_a + 1
    jb=j_b + 1
    jc=j_c + 1
    la=l_a + 1
    lb=l_b + 1
    lc=l_c + 1
    i=kh(ja+jb-jc+99)+kh(jb+jc-ja+99)+kh(jc+ja-jb+99)+kh(ja+lb-lc+99) &
         +kh(lb+lc-ja+99)+kh(lc+ja-lb+99)+kh(la+jb-lc+99)+kh(jb+lc-la+99) &
         +kh(lc+la-jb+99)+kh(la+lb-jc+99)+kh(lb+jc-la+99)+kh(jc+la-lb+99)
    IF(i <= 0) THEN
       mt=(j_a+j_b+j_c)/2 + 2
       ma=(j_a+l_b+l_c)/2+ 2
       mb=(l_a+j_b+l_c)/2+ 2
       mc=(l_a+l_b+j_c)/2+ 2
       na=mt-ja
       nb=mt-jb
       nc=mt-jc
       ka=ma-lc
       kb=mb-lc
       kc=mc-jc
       fss=q(mt,ja+1)*q(ja,nc)/(q(ma,ja+1)*q(ja,ka)*q(mb,la+1)* &
            q(la,kb)*q(mc,la+1)*q(la,kc))
       fs=SQRT(fss)/(l_a + 1.)
       l0=MAX(mt,ma,mb,mc)+1
       l1=MIN(ma+na,mb+nb,mc+nc)
       x=0.
       DO l=l0,l1
          x=-x+q(l-1,mt)*q(na,l-ma)*q(nb,l-mb)*q(nc,l-mc)
       ENDDO
       sjs=-(1+2*(l1/2*2-l1))*fs*x
    ENDIF

  END FUNCTION sjs
  !
  !     calculates ninej-symbols
  !      
  REAL(DP) FUNCTION snj (ia,ib,ie,ic,id,if,ig,ih,it)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ia,ib,ie,ic,id,if,ig,ih,it
    INTEGER :: ja,jb,je,jc,jd,jf,jg,jh,jt,i,la,ld,ma,mc,na,nb,le,lf,&
         lg,me,mf,mg,ne,nf,ng,lx,mx,nx,jsi,jsf, js,is,lb, lc, &
         mb, ly, my,ny,l,l0,m0,n0,l1,m1,n1,m,n,ihx
    REAL(DP) :: x, fn, fd, ps, fs, u, y, z, ud, p

    snj=0.
    ja=ia+1
    jb=ib+1
    jc=ic+1
    jd=id+1
    je=ie+1
    jf=IF+1
    jg=ig+1
    jh=ih+1
    jt=it+1
    i=kh(ja+jb-je+99)+kh(jb+je-ja+99)+kh(je+ja-jb+99)+kh(jc+jd-jf+99) &
         +kh(jd+jf-jc+99)+kh(jf+jc-jd+99)+kh(jg+jh-jt+99)+kh(jh+jt-jg+99) &
         +kh(jt+jg-jh+99)+kh(ja+jc-jg+99)+kh(jc+jg-ja+99)+kh(jg+ja-jc+99) &
         +kh(jb+jd-jh+99)+kh(jd+jh-jb+99)+kh(jh+jb-jd+99)+kh(je+jf-jt+99) &
         +kh(jf+jt-je+99)+kh(jt+je-jf+99)
    IF(i <= 0) THEN
       la=(ie+IF+it)/2+2
       ld=(ig+ih+it)/2+2
       ma=(ia+ic+ig)/2+2
       mc=(IF+ic+id)/2+2
       na=(ib+id+ih)/2+2
       nb=(ib+ie+ia)/2+2
       le=(ie+IF-it)/2+1
       lf=(IF+it-ie)/2+1
       lg=(it+ie-IF)/2+1
       me=(ia+ic-ig)/2+1
       mf=(ic+ig-ia)/2+1
       mg=(ig+ia-ic)/2+1
       ne=(ib+id-ih)/2+1
       nf=(id+ih-ib)/2+1
       ng=(ih+ib-id)/2+1
       lx=(it+ig-ih)/2+1
       mx=(ic+id-IF)/2+1
       nx=(ib+ie-ia)/2+1
       fn=q(la,jt+1)*q(jt,lg)*q(ma,jc+1)*q(jc,mf)*q(na,jb+1)*q(jb,ne)
       fd=q(ld,jt+1)*q(jt,lx)*q(mc,jc+1)*q(jc,mx)*q(nb,jb+1)*q(jb,nx)
       jsi=MAX(ABS(je-jh),ABS(jg-jf),ABS(ja-jd))+1
       jsf=MIN(je+jh,jg+jf,ja+jd)-1
       ps=-1-2*(jsi/2*2-jsi)
       fs=ps*SQRT(fn/fd)/FLOAT((ig+1)*(ie+1))
       u=0.
       DO js=jsi,jsf,2
          is=js-1
          lb=(ie+ih+is)/2+2
          lc=(ig+IF+is)/2+2
          mb=(ia+id+is)/2+2
          ly=(ie+ih-is)/2+1
          my=(ig+IF-is)/2+1
          ny=(ia-id+is)/2+1
          ud=q(lb,je+1)*q(je,ly)*q(lc,jg+1)*q(jg,my)*q(mb,js+1)*q(js,ny)
          l0=MAX(la,lb,lc,ld)+1
          m0=MAX(ma,mb,mc,lc)+1
          n0=MAX(na,nb,mb,lb)+1
          l1=MIN(le+ld,lf+lb,lg+lc)
          m1=MIN(me+lc,mf+mb,mg+mc)
          n1=MIN(ne+lb,nf+nb,ng+mb)
          x=0.
          DO l=l0,l1
             x=-x-q(l-1,la)*q(le,l-ld)*q(lf,l-lb)*q(lg,l-lc)
          ENDDO
          y=0.
          DO m=m0,m1
             y=-y-q(m-1,ma)*q(me,m-lc)*q(mf,m-mb)*q(mg,m-mc)
          ENDDO
          z=0.
          DO n=n0,n1
             z=-z-q(n-1,na)*q(ne,n-lb)*q(nf,n-nb)*q(ng,n-mb)
          ENDDO
          ihx=l1+m1+n1
          p=1+2*(ihx/2*2-ihx)
          u=u+p*x*y*z/ud
       ENDDO
       snj=u*fs
    ENDIF

  END FUNCTION snj

  !
  !     This routine calculates the moshinsky vector bracket      
  !     Note that D=mass1/mass2                                   
  !     Ref  m.sotona and m.gmitro  comp.phys.comm 3(1972)53      
  !
  REAL(DP) FUNCTION gmosh &
       (n,l,nc,lc,n1,l1,n2,l2,lr,d)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n,l,nc,lc,n1,l1,n2,l2,lr
    REAL(DP), INTENT(IN) :: d
    INTEGER :: ip,ixf,ix, iyi, iyf, j1f,j2,k1i,k1f,m1f,iy,m2f,k2, &
         m2,m2i,m1,j1,k2f,k2i,k1
    REAL(DP) :: dl, d1l, bb, ba, anorm, y, p, bc, cfac, bm , &
         sm, s, sxy, bxy

    gmosh=0.
    IF(n+n+nc+nc+l+lc-n1-n1-n2-n2-l1-l2 /= 0 ) RETURN
    IF(l+lc-lr < 0 ) RETURN
    IF(l1+l2-lr < 0 ) RETURN
    IF(ABS(l-lc)-lr > 0 ) RETURN
    IF(ABS(l1-l2)-lr > 0 ) RETURN
    DL=LOG(D)
    D1L=LOG(D+1.)
    bb=f_mb(n1+1)+f_mb(n2+1)+f_mb(n+1)-f_mb(nc+1)+ &
         g_mb(n1+l1+1)+g_mb(n2+l2+1) &
         -g_mb(n+l+1)-g_mb(nc+lc+1)
    ba=w_mb(l1+1)+w_mb(l2+1)+w_mb(lc+1)+w_mb(l+1)+ &
         f_mb(l1+l2-lr+1)+f_mb(l+lc+lr+2) &
         +f_mb(l+lc-lr+1)+f_mb(lc+lr-l+1)+ &
         f_mb(lr+l-lc+1)-f_mb(l1+l2+lr+2) &
         -f_mb(l1+lr-l2+1)-f_mb(l2+lr-l1+1)-DBLE(l)*d1l
    ip=lr+n+n1+n2
    p=1+2*(ip/2*2-ip)
    anorm=p*EXP(0.5D0*(bb+ba))
    y=0.
    j1f=l+1
    DO j1=1,j1f
       j2=l+2-j1
       k1i=ABS(l1-j1+1)+1
       k1f=l1+j1
       DO k1=k1i,k1f,2
          m1f=n1-(j1+k1-l1)/2+2
          IF(m1f-1 < 0 )  CYCLE
          k2i=MAX(ABS(l2-j2+1),ABS(lc-k1+1))+1
          k2f=MIN(l2+j2,lc+k1)
          IF(k2i-k2f > 0 ) CYCLE
          DO k2=k2i,k2f,2
             m2f=n2-(j2+k2-l2)/2+2
             IF(m2f-1 < 0 )  CYCLE
             ip=j2-1+(l1+k1+j1+l2+k2+j2)/2
             p=1+2*(ip/2*2-ip)
             bc=0.5D0*(DBLE(k1+j2-2)*dl-DBLE(k1+k2-2)*d1l) &
                  +f_mb(k1+l1-j1+1)+f_mb(k1+k2-lc-1)+ &
                  f_mb(k2+l2-j2+1)-f_mb(k1+l1+j1)-f_mb(k1+k2+lc)- &
                  f_mb(k2+l2+j2)+w_mb(k1)+w_mb(k2)+f_mb((k1+l1+j1)/2)+ &
                  f_mb((k1+k2+lc)/2)+f_mb((k2+l2+j2)/2)- &
                  f_mb((k1+l1-j1)/2+1)-f_mb((l1+j1-k1)/2+1)- &
                  f_mb((j1+k1-l1)/2)-f_mb((k1+k2-lc)/2)- &
                  f_mb((k2+lc-k1)/2+1)-f_mb((lc+k1-k2)/2+1) &
                  -f_mb((k2+l2-j2)/2+1)-f_mb((l2+j2-k2)/2+1)- &
                  f_mb((j2+k2-l2)/2)
             cfac=p*EXP(bc)
             sxy=0.
             ixf=MIN(k1+k1,k1+k2-lc)-1
             DO ix=1,ixf
                iyi=MAX(1,ix+j1+l2-k1-lr)
                iyf=MIN(l2+l2+1,l1+l2-lr+1,l2+lc+ix-k1-j2+2)
                IF(iyi-iyf > 0 ) CYCLE
                DO iy=iyi,iyf
                   ip=ix+iy
                   p=1+2*(ip/2*2-ip)
                   bxy=f_mb(k1+k1-ix)+f_mb(l2+l2-iy+2)+ &
                        f_mb(k2+lc-k1+ix)+f_mb(l1+lr-l2+iy) &
                        -f_mb(ix)-f_mb(iy)-f_mb(k1+k2-lc-ix)- &
                        f_mb(l1+l2-lr-iy+2)-f_mb(k1-l2+lr-j1+iy-ix+1)- &
                        f_mb(l2-k1+lc-j2+ix-iy+3)
                   sxy=sxy+p*EXP(bxy)
                ENDDO
             ENDDO
             s=cfac*sxy
             sm=0.
             DO m1=1,m1f
                m2i=MAX(1,nc-m1-(k1+k2-lc)/2+3)
                IF(m2i-m2f > 0 ) CYCLE
                DO m2=m2i,m2f
                   ip=m1+m2
                   p=1+2*(ip/2*2-ip)
                   bm=DBLE(m1-1)*DL-DBLE(m1+m2-2)*d1l+g_mb(1) &
                        +g_mb(m1+m2+(k1+k2+lc)/2-2)-g_mb(k1+m1-1)- &
                        g_mb(k2+m2-1)+f_mb(m1+m2+(k1+k2-lc)/2-2)- &
                        f_mb(m1)-f_mb(m2)-f_mb(n1-m1-(j1+k1-l1)/2+3)- &
                        f_mb(n2-m2-(j2+k2-l2)/2+3) &
                        -f_mb(m1+m2-nc+(k1+k2-lc)/2-2)
                   sm=sm+p*EXP(bm)
                ENDDO
             ENDDO
             y=y+s*sm
          ENDDO
       ENDDO
    ENDDO
    gmosh=anorm*y

  END FUNCTION  gmosh

  !
  !     This routine calculates the vector bracket      
  !     allowing for a transformation from rel and com coordinates
  !     to the lab frame for plane waves. Present version assumes 
  !     identical masses only. Dimensionless.
  !

  REAL(DP) FUNCTION vector_trcoefficients(ak,akk,ak1,ak2,l,ll,lam,l1,l2,tisoz)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: ak,akk,ak1,ak2
    INTEGER , INTENT(IN) :: l,ll,lam,l1,l2,tisoz
    INTEGER :: ic1, ic2, m, mm, mu, mmu, mi, mp, ix, ixx, md
    REAL(DP) :: sak, dak, tmass, xm, xm1, xm2, xm3, x, xx, &
         xa, xb, xc, aiii, xd, psa, sl2, sll, sggn,&
         sumin, sign, sl, sl1, dr, dr1, dr2

    vector_trcoefficients=0.0_dp
    sak=ak1+ak2
    dak=DABS(ak1-ak2)
!    IF (akk > sak) RETURN
!    IF (akk < dak) RETURN
    !    tmass=p_mass(it1)+p_mass(it2)
    tmass=p_mass(tisoz)+p_mass(tisoz)
    !    xm=p_mass(it1)/tmass
    xm=p_mass(tisoz)/tmass
    !    xm1=(p_mass(it2)/tmass)**2
    xm1=(p_mass(tisoz)/tmass)**2
    xm2=xm**2
    !    xm3=p_mass(it1)*p_mass(it2)/(tmass**2)
    xm3=p_mass(tisoz)*p_mass(tisoz)/(tmass**2)
    x=(ak1*ak1-ak*ak-akk*akk*xm2)/(2.D0*xm*ak*akk)
    xx=x*x
    IF (xx > 1.0_dp) RETURN
    ic1=l1+l2+l+ll
    ic2=ic1-2*(ic1/2)
    IF (ic2 /= 0) RETURN
    aiii=0.0_dp
    xa=(ak1*ak1+ak*ak-akk*akk*xm2)/(2.d0*ak*ak1)
    xb=(ak1*ak1+akk*akk*xm2-ak*ak)/(2.D0*xm*ak1*akk)
    xc=(ak1*ak1*xm1+ak2*ak2*xm2-ak*ak)/(2.d0*ak1*ak2*xm3)
    md=0
    xd=1.0_dp
    sl1=spherical_harmonics(md,l1,xd)
    IF (sl1 == 0.0_dp) RETURN
    DO m=-l,l
       sl=spherical_harmonics(m,l,xa)
       IF (sl == 0.0_dp)  CYCLE
       DO mm=-ll,ll
          mu=m+mm
          IF (IABS(mu) > l2) CYCLE
          mmu=-mu
          dr1=tjs(l+l,ll+ll,lam+lam,m+m,mm+mm,mmu+mmu)
          dr2=tjs(l1+l1,l2+l2,lam+lam,md+md,mu+mu,mmu+mmu)
          dr=dr1*dr2
          sign=1.
          mi=m-l-l1+l2+ll
          mp=mi-2*(mi/2)
          IF(mp /= 0) sign=-1.0_dp
          sll= spherical_harmonics(mm,ll,xb)
          IF (sll == 0.0_dp) CYCLE
          sl2= spherical_harmonics(mu,l2,xc)
          IF (sl2 == 0.0_dp) CYCLE
          psa=sl*sll*sl1*sl2
          sumin=0.
          sumin=dr*psa*sign
          aiii=aiii+dr*psa*sign
       ENDDO
    ENDDO
    vector_trcoefficients=16.0d0*pi*pi*aiii
    sggn=1.
    ix=(l1+l2-l-ll)/2
    ixx=ix-2*(ix/2)
    IF (ixx /= 0) sggn=-1.
    vector_trcoefficients=vector_trcoefficients*sggn

  END FUNCTION  vector_trcoefficients

  !  Spherical harmonics from Num. Recipes 

  REAL(DP) FUNCTION spherical_harmonics(m1,l,x)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m1, l
    REAL(DP), INTENT(IN) ::  x
    REAL(DP), DIMENSION(0:51) :: y
    INTEGER :: iphase, m, j
    REAL(DP) :: fj, z, fac, div, sum, a, b, c
    spherical_harmonics=0.
    m=IABS(m1)
    IF(m.LT.0) m=-m1
    y(0)=1.
    IF(l.EQ.0) THEN
       sum=y(0)
    ELSE
       a=m-l
       b=l+m+1
       c=m+1
       z=0.5_dp-x*0.5_dp
       DO j=1,l-m+1
          fj=j-1
          y(j)=y(j-1)*(a+fj)*(b+fj)*z
          div=(c+fj)*(fj+1.)
          y(j)=y(j)/div
       ENDDO
       IF(m > 0) then
          fac=(1.-x*x)**m
          fac=SQRT(fac)
       ELSE
          fac=1.
       ENDIF
       sum=0.
       DO j=0,l-m
          sum=sum+y(j)
       ENDDO
       iphase=m
       IF(m1.LT.0) then
          iphase=0
       ENDIF
       sum=sum*fac*((-1)**iphase)
    ENDIF
    spherical_harmonics=cn(m,l)*sum
  END FUNCTION spherical_harmonics

END MODULE ang_mom_functions

!
!     arrays containing mesh points and harmonic oscillator wave functions
!     In addition, routines for various wave functions are also included
!

MODULE wave_functions
  USE constants
  COMPLEX(DPC), ALLOCATABLE, PUBLIC :: krel(:),wkrel(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: k_mesh(:), w_mesh(:), x_mesh(:), wx_mesh(:), t_mesh(:), wt_mesh(:), r_lab(:), wr_lab(:)
CONTAINS

  !
  !     H.O. functions using Kummers function   
  !
  !DOUBLE PRECISION FUNCTION rnl(n,l,z)
  COMPLEX*16 FUNCTION rnl(n,l,z)
  IMPLICIT NONE
  INTEGER :: lll, nn
  INTEGER, INTENT(IN) :: l, n
  DOUBLE PRECISION :: dl, gamfaa, pi, dfll, gamfab, dfnn !, hypkum !,y
  !DOUBLE PRECISION, INTENT(IN) :: z
  COMPLEX*16, INTENT(IN) :: z
  COMPLEX*16 :: y
 
  pi =  3.141592653589793d0
    rnl = dcmplx(0.) ; y= dcmplx( 0.5d0*z*z )
    IF(abs(y) > 60.d0) RETURN
    dl = l
    IF((ABS(z) < 1.0d-6) .AND. (l == 0)) rnl = dcmplx( 1.0d0,0.d0 )
    IF( ABS(z) > 1.0d-6) rnl = (z**l) * EXP(-y) * hypkum(n,dl+1.5,z*z)

    gamfaa = 0.5 * SQRT(pi)
    IF(l /= 0) THEN
       DO lll = 1, l
          dfll = lll - 1
          gamfaa = gamfaa * (dfll + 1.5)
       ENDDO
    ENDIF
    gamfab = gamfaa
    IF(n /= 0) THEN
       dfll = dl + 0.5
       DO nn = 1, n
          dfnn = nn
          gamfab = gamfab * ((dfnn + dfll) / dfnn)
       ENDDO
    ENDIF
    rnl = rnl * (SQRT(2.0 * gamfab) / gamfaa) 
 
  END FUNCTION rnl
  !
  !     Kummers function, Abramowitz & Stegun   
  !     exp. 13.1.2. a(there) equals (-n)       
  !  
  !DOUBLE PRECISION FUNCTION hypkum(n,b,z)
  COMPLEX*16 FUNCTION hypkum(n,b,z)
    IMPLICIT NONE
    INTEGER :: nmax, nf
    INTEGER, INTENT(IN)  :: n
    DOUBLE PRECISION :: af, bf, dfnf, xadd
    DOUBLE PRECISION, INTENT(IN) :: b 
    COMPLEX*16, INTENT(IN) :: z
    COMPLEX*16 :: zf, term, sum

    IF(n < 0) WRITE (6,*)' error exit in hypkum ',  n,b,z
    hypkum = dcmplx(1.0d0, 0.d0)
    IF(n == 0) RETURN
    nmax = n ; af = - n ; bf = b ; zf = z ; sum = 1.0 ; term = 1.0
    DO nf = 1, nmax
       dfnf = nf
       xadd = dfnf - 1.0
       term = term * ((af + xadd) / (bf + xadd)) * (zf / dfnf)
       IF(ABS(term) <  1.0d-12) EXIT
       sum = sum + term
    ENDDO
    hypkum = sum

  END FUNCTION hypkum


  !  This function sets up the recursive relation
  !  for the associated Legendre polynomials

  REAL(DP) FUNCTION legendre_polynomials(l, m, x)
    IMPLICIT NONE
    REAL(DP) ::  fact,pll,pmm,pmmp1,somx2
    REAL(DP), INTENT(IN)  :: x
    INTEGER ::  i,ll
    INTEGER, INTENT(IN) :: l, m

    !  check whether m, l and x are ok

    IF((M < 0).OR.(M > L).OR.(ABS(X) > 1.)) THEN
       WRITE(6,*) 'bad arguments', m, l, x; RETURN
    ENDIF

    !  calculate now pmm as starting point for iterations

    pmm=1.0
    IF (m > 0) THEN
       somx2=SQRT((1.0-x)*(1.0+x))
       fact=1.0_dp;
       DO i=1, m
          pmm = -fact*somx2*pmm
          fact = fact+2.0_dp
       ENDDO
    ENDIF

    !  if l == m we do not need to use recursion relation

    IF (l == m) THEN
       legendre_polynomials=pmm

       !  recursive relation for associated Legendre polynomials

    ELSE
       pmmp1=x*(2*m+1)*pmm

       !  analytical formula for the case l == m+1

       IF (l == (m+1)) THEN
          legendre_polynomials=pmmp1
       ELSE
          DO ll=m+2, l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          ENDDO
          legendre_polynomials= pll
       ENDIF
    ENDIF

  END FUNCTION legendre_polynomials
    !
  !
  !      This routine calculates gauss-legendre mesh points and weights      
  !      input:                                                              
  !      x1   : lower limit of the integration interval                      
  !      x2   : upper limit ---------- "" -------------                      
  !      n    : the desired number of mesh points                            
  !      output :                                                            
  !      x     : gauss-legendre mesh points on the interval (x1,x2)          
  !      w     : the corresponding weights                                   
  !      From  : Numerical recipes
  !      F90 version : M. Hjorth-Jensen
  !
  SUBROUTINE gauss_legendre(x1,x2,x,w,n)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j, m
    REAL(DP), INTENT(IN) :: x1, x2
    REAL(DP), INTENT(INOUT) :: x, w
    REAL(DP) :: eps
    DIMENSION :: x(n), w(n)
    PARAMETER (eps=3.D-14)
    REAL(DP) :: p1,p2,p3,pp,xl,xm,z,z1

    m=(n+1)/2
    xm=0.5_dp*(x2+x1)
    xl=0.5_dp*(x2-x1)
    DO i=1,m
       z1=0.
       z=COS(pi*(i-.25_dp)/(n+.5_dp))
       DO WHILE ( ABS(z-z1) > EPS)
          p1=1.0_dp
          p2=0.0_dp
          DO j=1,n
             p3=p2
             p2=p1
             p1=((2.0_dp*j-1.)*z*p2-(j-1.0_dp)*p3)/j
          ENDDO
          pp=n*(z*p1-p2)/(z*z-1.)
          z1=z
          z=z-p1/pp
       ENDDO
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.0_dp*xl/((1.0_dp-z*z)*pp*pp)
       w(n+1-i)=w(i)
    ENDDO

  END SUBROUTINE gauss_legendre


  !
  ! Lab frame meshes, used for setting up single-particle configurations
  ! Units are fm^-1, converted to MeV in heff-main.f90
  !

  SUBROUTINE lab_mesh
    IMPLICIT NONE
    REAL(DP) :: u1(n_lab1),s1(n_lab1), & 
         u2(n_lab2),s2(n_lab2), u3(n_lab3),s3(n_lab3) 
    real(DP) :: ymax
    integer :: i, nn1,nn2,nn3
    !   mesh points in units of fm^-1
    !ymax = tanh(k_cutoff)
    !CALL gauss_legendre(0.0_dp,ymax,u,s,n_lab)
    !CALL gauss_legendre(0.0_dp,k_cutoff*hbarc,u,s,n_lab)
    CALL gauss_legendre(0.d0,1.d0,u1,s1,n_lab1)
    CALL gauss_legendre(1.d0,300.d0,u2,s2,n_lab2)
    CALL gauss_legendre(300.d0,hole_cutoff,u3,s3,n_lab3)
    ! mesh points in units of MeV
    
    nn1 = n_lab1
    nn2 = n_lab1 + n_lab2
    nn3 = n_lab1 + n_lab2 + n_lab3
    do i= 1, n_lab
       if ( i <= nn1 ) then
          k_mesh(i)=u1(i)  
          w_mesh(i)=s1(i)
       elseif ( i > nn1 .and. i <= nn2 )then
          k_mesh(i)=u2(i-nn1) 
          w_mesh(i)=s2(i-nn1)
       elseif ( i > nn2 .and. i <= nn3 ) then
          k_mesh(i)=u3(i-nn2)  
          w_mesh(i)=s3(i-nn2)
       end if
    end do
    
    !write(6,*) k_mesh(1), w_mesh(1),sum( w_mesh(:)*exp(-0.01*k_mesh(:)) )
    
    !k_mesh=0.5d0*log((1.d0+u)/(1.d0-u)) * hbarc
    !w_mesh=s/(1-u*u) * hbarc

  END SUBROUTINE lab_mesh

  SUBROUTINE rspace_lab_mesh
    IMPLICIT NONE
    REAL(DP) :: u1(nrad),s1(nrad)
    integer :: i, nn1,nn2,nn3
    !   mesh points in units of fm^-1
    CALL gauss_legendre(0.d0,16.d0,u1,s1,nrad)
    
    ! mesh points in units of MeV
    
    r_lab=u1  
    wr_lab=s1
    
  END SUBROUTINE rspace_lab_mesh
!
!             Set up of k mesh and weights for vfree obtained
!             through a Lee-Suzuki similarity transformation
!

  SUBROUTINE rel_mesh
    IMPLICIT NONE
    INTEGER :: i
    REAL(DP) :: u1(n_k1), s1(n_k1),u2(n_k2), s2(n_k2)

    !   mesh points in units of MeV
    CALL gauss_legendre(0.D0,k_cutoff*hbarc,u1,s1,n_k1)
    CALL gauss_legendre(k_cutoff*hbarc, k_max*hbarc ,u2,s2,n_k2)
    DO i= 1, n_k
       IF ( i <= n_k1 ) THEN
          krel(i) = phase*DCMPLX( u1(i) )
          wkrel(i)= phase*DCMPLX( s1(i) )
       ELSEIF ( i > n_k1 ) THEN
          krel(i) = phase*DCMPLX( u2(i-n_k1) )
          wkrel(i)= phase*DCMPLX( s2(i-n_k1) )
       ENDIF
    ENDDO

  END SUBROUTINE rel_mesh


  !
  ! setup mesh used in Hartree-Fock integrations, using mapping (a,b) -> (-1,1)
  ! 
  subroutine hf_mesh
    implicit none
    
    CALL gauss_legendre(-1.D0,1.d0,x_mesh,wx_mesh,int_points)
    write(6,*) int_points
  end subroutine hf_mesh

  SUBROUTINE cheb_mesh
    IMPLICIT NONE
    ! set up integration points on interval [-1,1]
    ! for gauss-chebyshev quadrature
    CALL cheb_to_set ( int_points, t_mesh, wt_mesh )

  END SUBROUTINE cheb_mesh
  

END MODULE wave_functions

!     module which defines configuration type, general structure
!     either for lab frame case or relcm system
MODULE configurations
  USE wave_functions
  USE constants  
  TYPE configuration_descriptor        
     INTEGER :: number_confs, model_space_confs
     INTEGER, DIMENSION(:), POINTER :: config_ab
  END TYPE configuration_descriptor
  TYPE, PUBLIC ::  configuration_relcm        
     INTEGER, DIMENSION(:), POINTER :: nconfs_relcm
     INTEGER, DIMENSION(:,:), POINTER :: relcm_ab
  END TYPE configuration_relcm
  TYPE (configuration_relcm), PUBLIC :: relcm_conf

CONTAINS

  !                       
  !     setting up all configurations for given J, Tz and parity for
  !     the large space. Note energy cutoff, triangular cutoff.

  SUBROUTINE large_number_confs(ij,ipar,itz,this)
    USE single_particle_orbits
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, la, lb, ja, jb, a, b, b_end, nconfs, &
         itza, itzb
    LOGICAL triag
    REAL(DP) :: pa, pb
    nconfs=0
    DO a=1, all_orbit%total_orbits
       pa = all_orbit%p(a)
       la=all_orbit%ll(a)
       ja=all_orbit%jj(a)
       itza=all_orbit%itzp(a)
       IF ( itz /= 0 ) THEN
          b_end = a
       ELSE
          b_end = all_orbit%total_orbits
       ENDIF
       DO b=1, b_end
          pb = all_orbit%p(b)
          lb=all_orbit%ll(b)
          itzb=all_orbit%itzp(b)
          jb=all_orbit%jj(b)
          IF ( pa*pa+pb*pb > 2.0_dp*max_kinetic) CYCLE
          IF ( itzb > itza ) CYCLE
          IF (itzb+itza /= itz*2 ) CYCLE
          IF ( (-1)**(la+lb+ipar) < 0) CYCLE
          IF ( triag( 2*ij, ja, jb) ) CYCLE
          IF ( (a == b) .AND. (MOD(ij,2) /= 0) ) CYCLE
          nconfs=nconfs+1
       ENDDO
    ENDDO
    this%number_confs=nconfs

  END SUBROUTINE large_number_confs

  !
  !     setting up all configurations for given J, Tz and parity for
  !     the effective interaction model space. Note energy cutoff, triangular cutoff.
  !

  SUBROUTINE number_configurations_model(ij,ipar,itz,this)
    USE single_particle_orbits
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, la, lb, ja, jb, a, b, b_end, nconfs, &
         itza, itzb
    LOGICAL triag
    REAL(DP) :: pa, pb
    nconfs=0
    DO a=1, all_orbit%total_orbits
       IF (all_orbit%model_space(a) == 'outside') CYCLE
       la=all_orbit%ll(a)
       ja=all_orbit%jj(a)
       pa = all_orbit%p(a)
       itza=all_orbit%itzp(a)
       IF ( itz /= 0 ) THEN
          b_end = a
       ELSE
          b_end = all_orbit%total_orbits
       ENDIF
       DO b=1, b_end
          IF (all_orbit%model_space(b) == 'outside') CYCLE
          lb=all_orbit%ll(b)
          itzb=all_orbit%itzp(b)
          pb = all_orbit%p(b)
          jb=all_orbit%jj(b)
          IF ( pa*pa+pb*pb > 2.0_dp*max_kinetic) CYCLE
          IF ( itzb > itza ) CYCLE
          IF (itzb+itza /= itz*2 ) CYCLE
          IF ( (-1)**(la+lb+ipar) < 0) CYCLE
          IF ( triag( 2*ij, ja, jb) ) CYCLE
          IF ( (a == b) .AND. (MOD(ij,2) /= 0) ) CYCLE
          nconfs=nconfs+1
       ENDDO
    ENDDO
    this%model_space_confs=nconfs
  END SUBROUTINE number_configurations_model

  !
  !  Here we set up both the total number of configurations 
  !  and those pertaining to the model space. Note triangular cutoff in energy.
  !

  SUBROUTINE setup_configurations(ij,ipar,itz,this)
    USE single_particle_orbits
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, la, lb, ja, jb, a, b, b_end, &
         k1, k2, itza, itzb, j, n_model, total_confs 
    LOGICAL triag
    REAL(DP) :: pa, pb

    total_confs=0
    !   First we set up all model space configurations and sort them in ascending order 
    !   for printout of the effective interaction
    n_model=0
    DO a = 1, all_orbit%total_orbits
       IF (all_orbit%model_space(a) == 'outside') CYCLE
       la=all_orbit%ll(a)
       ja=all_orbit%jj(a)
       pa = all_orbit%p(a)
       itza=all_orbit%itzp(a)
       IF ( itz /= 0 ) THEN
          b_end = a
       ELSE
          b_end = all_orbit%total_orbits
       ENDIF
       DO b = 1, b_end
          IF (all_orbit%model_space(b) == 'outside') CYCLE
          lb=all_orbit%ll(b)
          jb=all_orbit%jj(b)
          pb = all_orbit%p(b)
          itzb=all_orbit%itzp(b)
          IF ( itzb > itza ) CYCLE
          IF (itzb+itza /= itz*2 ) CYCLE
          IF ( (-1)**(la+lb+ipar) < 0) CYCLE
          IF ( triag( 2*ij, ja, jb) ) CYCLE
          IF ( (a == b) .AND. (MOD(ij,2) /= 0) ) CYCLE
          IF ( pa*pa+pb*pb > 2.0_dp*max_kinetic) CYCLE
          n_model=n_model+1
          k2=n_model*2
          k1=k2-1
          this%config_ab(k1)=b
          this%config_ab(k2)=a
       ENDDO
    ENDDO
    IF ( n_model /= this%model_space_confs ) THEN
       WRITE(6,*) ' Error in number of model space' ; STOP
    ENDIF
    !   sort only as ascending series the model space configurations
    CALL sort_configs(this,n_model)
    !   Now we set up the remaining configurations, not ordered
    total_confs = n_model
    DO a = 1, all_orbit%total_orbits
       la=all_orbit%ll(a)
       ja=all_orbit%jj(a)
       pa = all_orbit%p(a)
       itza=all_orbit%itzp(a)
       IF ( itz /= 0 ) THEN
          b_end = a
       ELSE
          b_end = all_orbit%total_orbits
       ENDIF
       DO b = 1, b_end
          lb=all_orbit%ll(b)
          jb=all_orbit%jj(b)
          pb = all_orbit%p(b)
          itzb=all_orbit%itzp(b)
          IF ( pa*pa+pb*pb > 2.0_dp*max_kinetic) CYCLE
          IF ( itzb > itza ) CYCLE
          IF (itzb+itza /= itz*2 ) CYCLE
          IF ( (-1)**(la+lb+ipar) < 0) CYCLE
          IF ( triag( 2*ij, ja, jb) ) CYCLE
          IF ( (a == b) .AND. (MOD(ij,2) /= 0) ) CYCLE
          IF(((all_orbit%model_space(a) == 'inside').AND.(all_orbit%model_space(b) == 'outside')).OR. &
               ((all_orbit%model_space(b) == 'inside').AND.(all_orbit%model_space(a) == 'outside')).OR. &
               ((all_orbit%model_space(a) == 'outside').AND.(all_orbit%model_space(b) == 'outside')) )THEN
             total_confs=total_confs+1
             k2=total_confs*2
             k1=k2-1
             this%config_ab(k1)=b
             this%config_ab(k2)=a
          ENDIF
       ENDDO
    ENDDO
    IF (total_confs /= this%number_confs ) THEN
       WRITE(6,*) ' Error in configuration allocation ' ; STOP
    ENDIF
    SELECT CASE (itz )
    CASE ( -1)
       WRITE(6,*) total_confs,' Large space proton-proton configurations for J',ij
       WRITE(6,*) n_model,' Model space proton-proton configurations for J',ij
       !       DO j=1, this%number_confs
       !          WRITE(6,*) this%config_ab(j*2-1)*100+this%config_ab(j*2)
       !       ENDDO
    CASE (1)
       WRITE(6,*) total_confs,' Large space neutron-neutron configurations for J ',ij
       WRITE(6,*) n_model,' Model space neutron-neutron configurations for J ',ij
       !       DO j=1, this%number_confs
       !          WRITE(6,*) this%config_ab(j*2-1)*100+this%config_ab(j*2)
       !       ENDDO
    CASE (0)
       WRITE(6,*) total_confs,' Large space proton-neutron configurations for J ',ij
       WRITE(6,*) n_model,' Model space proton-neutron configurations for J ',ij
       !       DO j=1, this%number_confs
       !          WRITE(6,*) this%config_ab(j*2-1)*100+this%config_ab(j*2)
       !       ENDDO
    END SELECT

  END SUBROUTINE setup_configurations

  !
  !     setting up all configurations for given J, Tz and parity for
  !     the single particle model space
  !

  SUBROUTINE number_sp_confs(ij,ipar,itz,this)
    USE single_particle_orbits
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, a, nconfs

    nconfs=0
    DO a=1, all_orbit%total_orbits
       IF ( MOD(all_orbit%ll(a),2) /= ipar ) CYCLE
       IF ( all_orbit%jj(a) /= ij) CYCLE
       IF ( all_orbit%itzp(a) /= itz ) CYCLE
       IF (all_orbit%model_space(a) == 'outside') CYCLE
       nconfs=nconfs+1
    ENDDO
    this%number_confs=nconfs

  END SUBROUTINE number_sp_confs
  !
  !     Allocates space for the sp configurations and sets them up
  !
  SUBROUTINE setup_sp_configurations(ij,ipar,itz,this)
    USE single_particle_orbits
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: a, ij, ipar, itz, nconfs
    nconfs=0
    DO a = 1, all_orbit%total_orbits
       IF ( MOD(all_orbit%ll(a),2) /= ipar ) CYCLE
       IF ( all_orbit%jj(a) /= ij) CYCLE
       IF ( all_orbit%itzp(a) /= itz ) CYCLE
       IF (all_orbit%model_space(a) == 'outside') CYCLE
       nconfs=nconfs+1
       this%config_ab(nconfs)=a
    ENDDO
    IF ( nconfs /= this%number_confs ) THEN
       WRITE(6,*) ' Error in configuration allocation ' ; STOP
    ENDIF

  END SUBROUTINE setup_sp_configurations

  !
  !           Sort the two-particle configurations as
  !           ascending series
  !

  SUBROUTINE sort_configs(this,n_confs)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: i, j, tempa, tempb, n_confs

    DO i = 1, n_confs
       DO j =i, n_confs
          IF( (this%config_ab(2*i-1)*100 +this%config_ab(2*i)) >  &
               (this%config_ab(2*j-1)*100 +this%config_ab(2*j)) ) THEN 
             tempa=this%config_ab(2*i-1)
             tempb=this%config_ab(2*i)
             this%config_ab(2*i-1)=this%config_ab(2*j-1)
             this%config_ab(2*i)=this%config_ab(2*j)
             this%config_ab(2*j-1)=tempa
             this%config_ab(2*j)=tempb                             
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE sort_configs


  !
  !     setting up all configurations for given J, Tz, spin
  !     and parity for the rel and cm system specified by the
  !     number channel. First find the number of configurations.
  !     Note : plane wave basis
  !
  SUBROUTINE number_configurations_relcm(this)
    USE partial_waves
    USE constants
    USE single_particle_orbits
    IMPLICIT NONE
    TYPE (configuration_relcm), INTENT(INOUT) :: this
    INTEGER :: l_rel, l_cm, a, b, nconfs, channel, l_min, l_max

    DO channel=1, no_channels
       l_min=orb_lrel_min(channel)
       l_max=orb_lrel_max(channel)
       nconfs=0
       DO a=1, relcm_sp_data%max_value
          l_rel=relcm_sp_data%lrel(a)
          IF ( (l_rel /= l_max) .AND. ( l_rel /= l_min) ) CYCLE
          DO b=1, relcm_sp_data%max_value
             l_cm=relcm_sp_data%lrel(b)
             nconfs=nconfs+1
          ENDDO
       ENDDO
       this%nconfs_relcm(channel)=nconfs              
    ENDDO

  END SUBROUTINE number_configurations_relcm
  !
  !     Now we find the possible configurations
  !     Note : Harmonic oscillator basis
  !
  SUBROUTINE setup_configurations_relcm(this)
    USE partial_waves
    USE constants
    USE single_particle_orbits
    IMPLICIT NONE
    TYPE (configuration_relcm), INTENT(INOUT) :: this
    INTEGER :: l_rel, l_cm, a, b, channel, l_min, l_max, k1, k2, nconfs

    DO channel=1, no_channels
       l_min=orb_lrel_min(channel)
       l_max=orb_lrel_max(channel)
       nconfs=0
       DO a=1, relcm_sp_data%max_value
          l_rel=relcm_sp_data%lrel(a)
          IF ( (l_rel /= l_max) .AND. ( l_rel /= l_min ) ) CYCLE
          DO b=1, relcm_sp_data%max_value
             l_cm=relcm_sp_data%lrel(b)
             nconfs=nconfs+1
             k2=nconfs*2
             k1=k2-1
             this%relcm_ab(channel,k1)=a
             this%relcm_ab(channel,k2)=b
          ENDDO
       ENDDO
       IF ( nconfs /= this%nconfs_relcm(channel) ) THEN
          WRITE(6,*) ' Error in configuration allocation ' ; STOP
       ENDIF
       WRITE(6,*) nconfs,' Total configurations for channel', channel
    ENDDO

  END SUBROUTINE setup_configurations_relcm

  !
  !     Finding all possible values for rel and cm variables l, L
  !     Note : only rel and com orbital momenta
  !
  SUBROUTINE find_spdata_relcm(a)
    USE constants
    IMPLICIT NONE
    INTEGER , INTENT(OUT) :: a
    INTEGER :: l

    a=0
    DO l=0, lmax
       a=a+1
    ENDDO

  END SUBROUTINE find_spdata_relcm
  !
  !     setting up all quantum numbers for relative and
  !     c.m coordinates in the block relcm_sp_data
  !     Note : only rel and com orbital momenta
  !
  SUBROUTINE make_configurations_relcm
    USE constants
    USE single_particle_orbits
    IMPLICIT NONE
    INTEGER :: a, l

    a=0
    DO l=0, lmax
       a=a+1
       relcm_sp_data%lrel(a)=l
    ENDDO

  END SUBROUTINE make_configurations_relcm

END MODULE configurations

! complex hyperbolic functions not defined 
! as intrinsic procedures in fortran
MODULE complexfunctions
  USE constants 
  IMPLICIT NONE
  COMPLEX(DPC) :: z 
  REAL(DP) :: re_z
CONTAINS
  REAL(DP) FUNCTION  ATANH(re_z)
    REAL(DP), INTENT(IN) :: re_z
    ATANH = 0.5d0*LOG((1+re_z)/(1-re_z))
  END FUNCTION ATANH
  
  REAL(DP) FUNCTION  sec(re_z)
    REAL(DP), INTENT(IN) :: re_z
    sec = 1/cosh(re_z)
  END FUNCTION SEC
  
  COMPLEX(DPC) FUNCTION cmplxatanh(z)
    COMPLEX(DPC), INTENT(IN) :: Z
    CMPLXATANH = 0.5d0*log((1.d0+z)/(1.d0-z))
  END FUNCTION CMPLXATANH
  
  COMPLEX(DPC) FUNCTION cmplxsinh(z) 
    COMPLEX(DPC), INTENT(IN) :: Z
    CMPLXSINH = -dcmplx(0.,1.) * sin(dcmplx(0.,1.)*z)
  END FUNCTION CMPLXSINH
  
  COMPLEX(DPC) FUNCTION cmplxcosh(z) 
    COMPLEX(DPC), INTENT(IN) :: Z
    CMPLXCOSH = cos(dcmplx(0.,1.)*z)
  END FUNCTION CMPLXCOSH
  
  COMPLEX(DPC) FUNCTION cmplxcsch(z)
    COMPLEX(DPC), INTENT(IN) :: Z
    ! CSH = 1/SINH
    cmplxCSCH = dcmplx(0.,1.)/ sin(dcmplx(0.,1.)*z)
  END FUNCTION CMPLXCSCH

  COMPLEX(DPC) FUNCTION cmplxtanh(z)
    COMPLEX(DPC), INTENT(IN) :: Z
    CMPLXTANH = -dcmplx(0.,1.)*sin(dcmplx(0.,1.)*z)/cos(dcmplx(0.,1.)*z)
  END FUNCTION CMPLXTANH
  
END MODULE complexfunctions

module setup_cmplx_mesh
  complex*16, ALLOCATABLE, PUBLIC :: k_cmplx(:),wk_cmplx(:)
end module setup_cmplx_mesh
