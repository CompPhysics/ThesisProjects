
!     Definition of single particle data

MODULE single_particle_orbits
  USE constants
  TYPE, PUBLIC :: single_particle_descript
     INTEGER :: total_orbits
     INTEGER, DIMENSION(:), POINTER :: nn, ll, jj, itzp, nshell,mvalue, &
     ptype, strange
     CHARACTER (LEN=10), DIMENSION(:), POINTER :: orbit_status, model_space,&
        included
     REAL(DP), DIMENSION(:), POINTER :: e, evalence, e_original
  END TYPE single_particle_descript
  TYPE, PUBLIC :: rel_cm_data
     INTEGER, DIMENSION(:), POINTER :: nrel, lrel
     INTEGER :: max_value
  END TYPE rel_cm_data
  TYPE (rel_cm_data), PUBLIC :: relcm_sp_data
  TYPE (single_particle_descript), PUBLIC :: all_orbit, neutron_data, &
       proton_data, lambda_data, sigmap_data, sigma0_data, sigmam_data, &
       mscheme_basis
CONTAINS
  SUBROUTINE allocate_sp_array(this_array,n)
    TYPE (single_particle_descript), INTENT(INOUT) :: this_array
    INTEGER , INTENT(IN) :: n
    INTEGER :: i
    IF (ASSOCIATED (this_array%nn) ) DEALLOCATE(this_array%nn)
    ALLOCATE(this_array%nn(n))
    IF (ASSOCIATED (this_array%ll) ) DEALLOCATE(this_array%ll)
    ALLOCATE(this_array%ll(n))
    IF (ASSOCIATED (this_array%jj) ) DEALLOCATE(this_array%jj)
    ALLOCATE(this_array%jj(n))
    IF (ASSOCIATED (this_array%itzp) ) DEALLOCATE(this_array%itzp)
    ALLOCATE(this_array%itzp(n))
    IF (ASSOCIATED (this_array%mvalue) ) DEALLOCATE(this_array%mvalue)
    ALLOCATE(this_array%mvalue(n))
    IF (ASSOCIATED (this_array%ptype) ) DEALLOCATE(this_array%ptype)
    ALLOCATE(this_array%ptype(n))
    IF (ASSOCIATED (this_array%strange) ) DEALLOCATE(this_array%strange)
    ALLOCATE(this_array%strange(n))
    IF (ASSOCIATED (this_array%e) ) DEALLOCATE(this_array%e)
    ALLOCATE(this_array%e(n))
    IF (ASSOCIATED (this_array%evalence) ) DEALLOCATE(this_array%evalence)
    ALLOCATE(this_array%evalence(n))
    IF (ASSOCIATED (this_array%e_original) ) DEALLOCATE(this_array%e_original)
    ALLOCATE(this_array%e_original(n))
    IF (ASSOCIATED (this_array%nshell) ) DEALLOCATE(this_array%nshell)
    ALLOCATE(this_array%nshell(n))
    IF (ASSOCIATED (this_array%orbit_status) ) DEALLOCATE(this_array%orbit_status)
    ALLOCATE(this_array%orbit_status(n))
    IF (ASSOCIATED (this_array%model_space) ) DEALLOCATE(this_array%model_space)
    ALLOCATE(this_array%model_space(n))
    IF (ASSOCIATED (this_array%included) ) DEALLOCATE(this_array%included)
    ALLOCATE(this_array%included(n))

    !           blank all characters and zero all other values
    DO i= 1, n
       this_array%model_space(i)= ' '
       this_array%orbit_status(i)= ' '
       this_array%included(i)= ' '
       this_array%e(i)=0.0_dp
       this_array%evalence(i)=0.0_dp
       this_array%e_original(i)=0.0_dp
       this_array%nn(i)=0
       this_array%ll(i)=0
       this_array%jj(i)=0
       this_array%nshell(i)=0
       this_array%itzp(i)=0
       this_array%mvalue(i)=0
       this_array%ptype(i)=0
       this_array%strange(i)=0
    ENDDO

  END SUBROUTINE allocate_sp_array

  SUBROUTINE deallocate_sp_array(this_array)
    TYPE (single_particle_descript), INTENT(INOUT) :: this_array
    DEALLOCATE(this_array%nn) ; DEALLOCATE(this_array%ll)
    DEALLOCATE(this_array%jj) ;DEALLOCATE(this_array%itzp)
    DEALLOCATE(this_array%evalence) ;DEALLOCATE(this_array%mvalue)
    DEALLOCATE(this_array%e) ;DEALLOCATE(this_array%nshell)
    DEALLOCATE(this_array%orbit_status); DEALLOCATE(this_array%model_space)
    DEALLOCATE(this_array%ptype); DEALLOCATE(this_array%strange)
    DEALLOCATE(this_array%e_original);DEALLOCATE(this_array%included)
  END SUBROUTINE deallocate_sp_array

  SUBROUTINE allocate_relcm_array(this_array,n)
    TYPE (rel_cm_data), INTENT(INOUT) :: this_array
    INTEGER , INTENT(IN) :: n
    INTEGER :: i
    IF (ASSOCIATED (this_array%nrel) ) DEALLOCATE(this_array%nrel)
    ALLOCATE(this_array%nrel(n))
    IF (ASSOCIATED (this_array%lrel) ) DEALLOCATE(this_array%lrel)
    ALLOCATE(this_array%lrel(n))
    !           blank all characters and zero all other values
    DO i= 1, n
       this_array%nrel(i)=0
       this_array%lrel(i)=0
    ENDDO

  END SUBROUTINE allocate_relcm_array

  SUBROUTINE deallocate_relcm_array(this_array)
    TYPE (rel_cm_data), INTENT(INOUT) :: this_array
    DEALLOCATE(this_array%nrel) ; DEALLOCATE(this_array%lrel)

  END SUBROUTINE deallocate_relcm_array

END MODULE single_particle_orbits

!     In this module we set up all one-body contributions. Saves CPU
!     time when setting up two- or three-body contributions
MODULE onebody_diagrams
  USE constants
  REAL(DP), ALLOCATABLE, PUBLIC :: one_body_terms(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: one_body_folded(:,:)
END MODULE onebody_diagrams

!     This module sets up the mpi start and end points
!     for where to search for Relative and cm configs
!     In addition we include the variables size_mpi
!     and my_rank_mpi given in the main function through
!     the calls to
!      CALL MPI_INIT ( err_mpi )
!      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, size_mpi, err_mpi )
!      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank_mpi, err_mpi )
!      which initialize MPI

MODULE jobs_mpi     
  INTEGER, PUBLIC :: length_config_all, length_config_mpi
  INTEGER, PUBLIC :: size_mpi, my_rank_mpi
  TYPE, PUBLIC :: jobs_mpi_descriptor    
     !        MPI modification
     INTEGER, DIMENSION(:), POINTER :: loop_start, loop_end
     INTEGER, DIMENSION(:,:), POINTER :: loop_start_all, loop_end_all
  END TYPE jobs_mpi_descriptor
  TYPE (jobs_mpi_descriptor), PUBLIC :: jobs

END MODULE jobs_mpi

!     the g-matrix in the rel and c.m. frame
MODULE relcm_gmatrix
  USE constants
  COMPLEX(DPC), ALLOCATABLE, PUBLIC :: v_com_rel(:,:,:)
  COMPLEX(DPC), ALLOCATABLE, PUBLIC :: gtf_tot_all(:,:,:), &
       gtf_tot_rank(:,:,:), &
       gtf_buffer(:,:,:)
END MODULE relcm_gmatrix

MODULE hoosc_gmatrix
  USE constants
  INTEGER, PUBLIC :: no_elements
  REAL(DP), ALLOCATABLE, PUBLIC :: gf(:,:), vcom(:), vr1r2(:), vp1p2(:)

END MODULE hoosc_gmatrix

!     setup all partial wave combinations from jmin to jmax
!     defined on input, as default use jmin = 0 and jmax = 10
MODULE partial_waves
  USE constants
  INTEGER, PUBLIC, PARAMETER ::  max_sub = 4
  INTEGER, PUBLIC :: no_channels
  INTEGER,PUBLIC, ALLOCATABLE:: orb_lrel_min(:),spin_rel(:), &
       jang_rel(:),iso(:), orb_lrel_max(:), strange(:), &
       no_sub(:), sub_id(:,:)
  REAL(KIND(1.0D0)), PUBLIC, ALLOCATABLE :: mass(:,:)
END MODULE partial_waves

!           
!     This module contains the angular momentun functions
!     and transformation coefficients when going from
!     lab system  <--> cm system
!
MODULE ang_mom_functions
  USE constants
  REAL(DP), PRIVATE :: f_mb(100),g_mb(100),w_mb(100), sfact(100), dfact(100)
  INTEGER, PRIVATE :: kh(200)
  REAL(DP), PRIVATE :: q(100,100), cn(0:101,0:101)

CONTAINS
  !
  !     factorials for 3j,6j and 9j symbols            
  !     for moshinsky trans brackets and for           
  !     vector brackets                                
  !
  SUBROUTINE commons_to_angmom
    IMPLICIT NONE
    INTEGER :: l, k, i, j
    REAL(DP) :: a , sq_pi, fj, tfj, fk, s
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
    ! Legendre functions of the second kind
    sfact(1) = 1.0_dp
    dfact(1) = 1.0_dp
    DO i = 1,30
       s = i
       sfact(i+1) = sfact(i)*s
       dfact(i+1) = (s+s+1.0d0)*dfact(i)
    ENDDO

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
       IF(i == 0) THEN
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
  REAL(DP) FUNCTION snj (ia,ib,ie,ic,id,IF,ig,ih,it)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ia,ib,ie,ic,id,IF,ig,ih,it
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

  !     calculates the vector bracket                           
  !     <ak,l,akk,ll,lam | ak1,l1,ak2,l2,lam>                   
  !     Ref. Kung et al, PRC 19 (1979) 1063                     
  !     Momenta in units of MeV, vector bracket in Units of MeV^-4

  REAL(DP) FUNCTION vector_trcoefficients &
       (ak,akk,ak1,ak2,l,ll,lam,l1,it1,l2,it2)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: ak,akk,ak1,ak2
    INTEGER , INTENT(IN) :: l,ll,lam,l1,it1,l2,it2
    INTEGER :: ic1, ic2, m, mm, mu, mmu, mi, mp, ix, ixx, md
    REAL(DP) :: sak, dak, tmass, xm, xm1, xm2, xm3, x, xx, &
         xa, xb, xc, aiii, xd, psa, sl2, sll, sggn,&
         sumin, sign, sl, sl1, dr, dr1, dr2

    vector_trcoefficients=0.
    sak=ak1+ak2
    dak=DABS(ak1-ak2)
    IF (akk > sak) RETURN
    IF (akk < dak) RETURN
    !    tmass=p_mass(it1)+p_mass(it2)
    !    xm=p_mass(it1)/tmass
    !    xm1=(p_mass(it2)/tmass)**2
    xm2=xm**2
    !    xm3=p_mass(it1)*p_mass(it2)/(tmass**2)
    x=(ak1*ak1-ak*ak-akk*akk*xm2)/(2.0_dp*xm*ak*akk)
    xx=x*x
    IF (xx > 1.) RETURN
    ic1=l1+l2+l+ll
    ic2=ic1-2*(ic1/2)
    IF (ic2 /= 0) RETURN
    aiii=0.
    xa=(ak1*ak1+ak*ak-akk*akk*xm2)/(2.0_dp*ak*ak1)
    xb=(ak1*ak1+akk*akk*xm2-ak*ak)/(2.0_dp*xm*ak1*akk)
    xc=(ak1*ak1*xm1+ak2*ak2*xm2-ak*ak)/(2.0_dp*ak1*ak2*xm3)
    md=0
    xd=1.
    sl1=spherical_harmonics(md,l1,xd)
    IF (sl1 == 0.0_dp) RETURN
    DO m=-l,l
       sl=spherical_harmonics(m,l,xa)
       IF (sl == 0.)  CYCLE
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
          IF (sll == 0) CYCLE
          sl2= spherical_harmonics(mu,l2,xc)
          IF (sl2 == 0.0_dp) CYCLE
          psa=sl*sll*sl1*sl2
          sumin=0.
          sumin=dr*psa*sign
          aiii=aiii+dr*psa*sign
       ENDDO
    ENDDO
    vector_trcoefficients=16.0_dp*pi*pi*aiii
    sggn=1.
    ix=(l1+l2-l-ll)/2
    ixx=ix-2*(ix/2)
    IF (ixx /= 0) sggn=-1.
    vector_trcoefficients=vector_trcoefficients/(ak*akk*ak1*ak2)*sggn

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
       IF(m > 0) THEN
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
       IF(m1.LT.0) THEN
          iphase=0
       ENDIF
       sum=sum*fac*((-1)**iphase)
    ENDIF
    spherical_harmonics=cn(m,l)*sum

  END FUNCTION spherical_harmonics

  !
  !      Legendre function of the second kind. See        
  !      Haftel and Tabakin Nucl. Phys. A145 (1970) p. 1  
  !      for further references                           
  !             

  REAL(DP) FUNCTION qlx(l,x)
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: l
    REAL(DP), INTENT(IN) :: x
    INTEGER :: k, n
    REAL(DP) :: y, yd, a, b, c, gk

    qlx = 0.0_dp
    IF(x-100.0_dp >= 0.0_dp) THEN           !  For large X values.
       y = 0.0_dp
       DO n = 1, 5
          yd = sfact(n)*dfact(l+n)*((2.0d0*x*x)**(n-1))
          y = y + sfact(l+n+n-1)/yd
       ENDDO
       qlx = y/(x**(l+1))
    ELSEIF((x-100.0_dp < 1.0_dp).AND.(l >= 0) ) THEN
       a = 0.5_dp*LOG((x+1.0_dp)/(x-1.0_dp))
       IF(l == 0) THEN
	  qlx = a
       ELSEIF(l > 0) THEN
	  b = a
	  a = b*x-1.0_dp
	  IF(l-1 == 0) THEN
             qlx = a
	  ELSEIF(L-1  > 0) THEN
             DO k =2, l
                c = b
                b = a
                gk = 1.0_dp/k
                a = (2.0_dp-gk)*x*b - (1.0_dp - gk)*c
             ENDDO
             qlx = a
	  ENDIF
       ENDIF
    ENDIF

  END FUNCTION qlx

END MODULE ang_mom_functions
!
!     Modules specific to the g-matrix calculation and effective operators
!     only
!     arrays containing mesh points and harmonic oscillator wave functions
!     In addition, routines for various wave functions are also included
!
MODULE wave_functions  
  USE constants
  INTEGER , PUBLIC :: n_k1, n_k, n_k2
  REAL(DP), ALLOCATABLE, PUBLIC :: ra(:),wra(:), krel(:), wkrel(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: rgkk(:),wgkk(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: rca(:),wrca(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: hol(:,:,:,:)
  ! ho onebody energies used for RG flow equations
  REAL(DP), PUBLIC, ALLOCATABLE :: ho_onebodyenergy(:)
  REAL(DP), PUBLIC :: k_cutoff, k_max
  COMPLEX(DPC), ALLOCATABLE, PUBLIC :: chol(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: bhf_hol(:,:)
  ! For hf wavefunctions
  REAL(DP), ALLOCATABLE, PUBLIC :: wave_function(:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: rnlc(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: rnlr(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: propagator(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: osc_propagator(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: coulomb_relcom(:,:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC:: wave_cm(:)
CONTAINS

  !
  !     H.O. functions using Kummers function   
  !
  REAL(DP) FUNCTION rnl(n,l,z)
    IMPLICIT NONE
    INTEGER :: lll, nn
    INTEGER, INTENT(IN) :: l, n
    REAL(DP) :: y, dl, gamfaa, dfll, gamfab, dfnn
    REAL(DP), INTENT(IN) :: z

    rnl=0. ; y=0.5_dp*z*z
    IF(y > 60.0_dp) RETURN
    dl = l
    IF((ABS(z) < 1.0d-6) .AND. (l == 0)) rnl = 1.0_dp
    IF( ABS(z) > 1.0d-6) rnl = (z**l) * EXP(-y) * hypkum(n,dl+1.5_dp,z*z)
    gamfaa = 0.5_dp * SQRT(pi)
    IF(l /= 0) THEN
       DO lll = 1, l
          dfll = lll - 1
          gamfaa = gamfaa * (dfll + 1.5_dp)
       ENDDO
    ENDIF
    gamfab = gamfaa
    IF(n /= 0) THEN
       dfll = dl + 0.5_dp
       DO nn = 1, n
          dfnn = nn
          gamfab = gamfab * ((dfnn + dfll) / dfnn)
       ENDDO
    ENDIF
    rnl = rnl * (SQRT(2.0_dp * gamfab) / gamfaa)

  END FUNCTION rnl
  !
  !     Kummers function, Abramowitz & Stegun   
  !     exp. 13.1.2. a(there) equals (-n)       
  !  
  REAL(DP) FUNCTION hypkum(n,b,z)
    IMPLICIT NONE
    INTEGER :: nmax, nf
    INTEGER, INTENT(IN)  :: n
    REAL(DP) :: af, bf, zf, term, dfnf, xadd, sum
    REAL(DP), INTENT(IN) :: b, z

    IF(n < 0) WRITE (6,*)' error exit in hypkum ',  n,b,z
    hypkum = 1.0_dp
    IF(n == 0) RETURN
    nmax = n ; af = - n ; bf = b ; zf = z ; sum = 1.0 ; term = 1.0_dp
    DO nf = 1, nmax
       dfnf = nf
       xadd = dfnf - 1.0_dp
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
  !             Set up of c.m. mesh and weights
  !
  SUBROUTINE cm_mesh
    IMPLICIT NONE
    REAL(DP), DIMENSION(n_cm_mom) :: u,s

    CALL gauss_legendre(0.0_dp,cutoff,u,s,n_cm_mom)
    rca=u
    wrca=s

  END SUBROUTINE cm_mesh
  !
  !             Set up of relative mesh and weights
  !
  SUBROUTINE nocore_mesh
    IMPLICIT NONE
    REAL(KIND = 8) :: u,s,c
    PARAMETER (c=0.75)
    DIMENSION u(n_rel), s(n_rel)
    !   set up points for interpolation with HO basis, limited extent
    !   of the wave function
    CALL gauss_legendre(0.D0,cutoff,u,s,n_rel)
    ra=u
    wra=s

  END SUBROUTINE nocore_mesh
  !
  !             Set up of k mesh and weights for vfree obtained
  !             through a Lee-Suzuki similarity transformation
  !
  SUBROUTINE vlowk_mesh
    IMPLICIT NONE
    INTEGER :: i
    REAL(KIND=8) :: u1(n_k1), s1(n_k1),u2(n_k2), s2(n_k2)

    !   mesh points in units of fm-1
    CALL gauss_legendre(0.D0,k_cutoff,u1,s1,n_k1)
    CALL gauss_legendre(k_cutoff, k_max ,u2,s2,n_k2)
    DO i= 1, n_k
       IF ( i <= n_k1 ) THEN
          ra(i) = u1(i)
          wra(i)= s1(i)
       ELSEIF ( i > n_k1 ) THEN
          ra(i) = u2(i-n_k1)
          wra(i)= s2(i-n_k1)
       ENDIF
    ENDDO
    krel = ra*hbarc
    wkrel = wra*hbarc
    WRITE(6,*) 'points for ho oscillator basis:',  n_rel
    WRITE(6,'(5D16.8)') ra
    WRITE(6,'(5D16.8)') wra

  END SUBROUTINE vlowk_mesh
  !
  !             Set up of relative mesh and weights
  !
  SUBROUTINE rel_mesh
    IMPLICIT NONE
    INTEGER :: i
    REAL(DP) :: u,s,xx,c
    PARAMETER (c=0.75_dp)
    DIMENSION u(n_rel), s(n_rel)
    !   set up mesh points for G-mat calc with x \in [0,\infty]
    !   mapping of Gauss-Legendre mesh points u \in [-1,1]
    CALL gauss_legendre(-1.0_dp,1.0_dp,u,s,n_rel)
    DO i=1,n_rel
       xx=pi_4*(u(i)+1.0_dp); rgkk(i)=DTAN(xx)*c
       wgkk(i)=pi_4*c/DCOS(xx)**2*s(i)
    ENDDO
    WRITE(6,*) 'points for gkk matrix calculation:',  n_rel
    WRITE(6,*) "q:"
    WRITE(6,'(5D16.8)') rgkk
    WRITE(6,*) "w:"
    WRITE(6,'(5D16.8)') wgkk
!   set up points for interpolation with HO basis, limited extent
    CALL gauss_legendre(0.0_dp,cutoff,u,s,n_rel)
    ra=u
    wra=s
!    ra=rgkk
!    wra=wgkk
    WRITE(6,*) 'points for ho oscillator basis:',  n_rel
    WRITE(6,*) "q:"
    WRITE(6,'(5D16.8)') ra
    WRITE(6,*) "w:"
    WRITE(6,'(5D16.8)') wra

  END SUBROUTINE rel_mesh

  !
  !             Set up of k mesh and weights for vfree obtained
  !             through the solution of RG flow equations
  !
  SUBROUTINE vkrgk_mesh               
    IMPLICIT NONE
    INTEGER :: i
    REAL(KIND=8) :: u(n_k), s(n_k)

    !   mesh points in units of fm-1
    CALL gauss_legendre(0.D0, k_max ,u,s,n_k)
    DO i= 1, n_k
       ra(i) = u(i)
       wra(i)= s(i)
    ENDDO
    krel = ra*hbarc
    wkrel = wra*hbarc
    WRITE(6,*) 'points for ho oscillator basis:',  n_rel
    WRITE(6,'(5D16.8)') ra
    WRITE(6,'(5D16.8)') wra

  END SUBROUTINE vkrgk_mesh


  SUBROUTINE laguerre_general( n, alpha, x, cx )
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL (dp ) ::  alpha
    REAL ( dp ) :: cx(0:n)
    INTEGER :: i
    REAL ( dp ), INTENT(IN) ::  x

    IF ( alpha <= -1.0D+00 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'LAGUERRE_GENERAL - Fatal error!'
       WRITE ( *, '(a,g14.6)' ) '  The input value of ALPHA is ', alpha
       WRITE ( *, '(a)' ) '  but ALPHA must be greater than -1.'
       STOP
    END IF
    IF ( n < 0 ) THEN
       RETURN
    END IF
    cx(0) = 1.0D+00
    IF ( n == 0 ) THEN
       RETURN
    END IF
    cx(1) = 1.0D+00 + alpha - x
    DO i = 2, n
       cx(i) = ( ( REAL ( 2 * i - 1, kind = 8 ) + alpha - x ) * cx(i-1)   &
            + ( REAL (   - i + 1, kind = 8 ) - alpha     ) * cx(i-2) ) &
            / REAL (     i,     kind = 8 )
    END DO

  END SUBROUTINE laguerre_general

  SUBROUTINE laguerre_complex( n, alpha, x, cx )
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL (dp ) ::  alpha
    COMPLEX ( DPC ) :: cx(0:n)
    INTEGER :: i
    COMPLEX ( DPC ), INTENT(IN) ::  x

    IF ( alpha <= -1.0D+00 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'LAGUERRE_GENERAL - Fatal error!'
       WRITE ( *, '(a,g14.6)' ) '  The input value of ALPHA is ', alpha
       WRITE ( *, '(a)' ) '  but ALPHA must be greater than -1.'
       STOP
    END IF
    IF ( n < 0 ) THEN
       RETURN
    END IF
    cx(0) = 1.0D+00
    IF ( n == 0 ) THEN
       RETURN
    END IF
    cx(1) = 1.0D+00 + alpha - x
    DO i = 2, n
       cx(i) = ( ( REAL ( 2 * i - 1, kind = 8 ) + alpha - x ) * cx(i-1)   &
            + ( REAL (   - i + 1, kind = 8 ) - alpha     ) * cx(i-2) ) &
            / REAL (     i,     kind = 8 )
    END DO

  END SUBROUTINE laguerre_complex

  DOUBLE PRECISION FUNCTION  fac(m)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: m
    INTEGER :: i

    fac = 0.0D0
    IF(m == 0) RETURN
    DO i=1,m
       fac=fac+LOG(FLOAT(i))
    ENDDO

  END FUNCTION  fac

  DOUBLE PRECISION FUNCTION  dfac(m)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: m
    INTEGER :: i

    IF (MOD(m,2).NE.1) STOP 'wrong argument to dfac'
    dfac = 0.0D0
    IF (m == 1)RETURN
    DO i=3,m,2
       dfac=dfac+LOG(FLOAT(i))
    ENDDO
  END FUNCTION  dfac

END MODULE wave_functions

!     module which defines configuration type, general structure
!     either for lab frame case or relcm system
MODULE configurations
  USE constants
  USE single_particle_orbits
  USE partial_waves
  USE ang_mom_functions
  TYPE configuration_descriptor        
     INTEGER :: number_confs
     INTEGER, DIMENSION(:), POINTER :: config_ab
  END TYPE configuration_descriptor
  TYPE, PUBLIC ::  configuration_relcm        
     INTEGER, DIMENSION(:), POINTER :: nconfs_relcm
     INTEGER, DIMENSION(:,:), POINTER :: relcm_ab, rel_ab 
     INTEGER, DIMENSION(:), POINTER :: nconfs_rel, nconfs_relmodel
  END TYPE configuration_relcm
  TYPE (configuration_relcm), PUBLIC :: relcm_conf, rel_conf

CONTAINS
  !                       
  !     setting up all configurations for given J, 2Tz, parity and
  !     strangeness for the g-matrix model space

  SUBROUTINE number_gmatrix_confs(ij,ipar,itz,istr, this)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, na, la, nb, lb, ja, jb, a, b, nconfs, &
         itza, itzb, istr, istra, istrb, ida, idb
    LOGICAL triag
    nconfs=0
    DO a=1, all_orbit%total_orbits
       na = all_orbit%nn(a)
       la=all_orbit%ll(a)
       ja=all_orbit%jj(a)
       itza=all_orbit%itzp(a)
       istra=all_orbit%strange(a)
       ida = all_orbit%ptype(a)
       DO b=1, all_orbit%total_orbits
          nb = all_orbit%nn(b)
          lb=all_orbit%ll(b)
          itzb=all_orbit%itzp(b)
          jb=all_orbit%jj(b)
          istrb=all_orbit%strange(b)
          idb = all_orbit%ptype(b)
          IF ((all_orbit%model_space(a) == 'outside').AND.  &
               (all_orbit%model_space(b) == 'outside') ) CYCLE
          IF ((square_calculation).AND.(all_orbit%model_space(a) == 'inside').AND.  &
               (all_orbit%model_space(b) == 'outside') ) CYCLE
          IF ((square_calculation).AND.(all_orbit%model_space(a) == 'outside').AND.  &
               (all_orbit%model_space(b) == 'inside') ) CYCLE
          IF ((ida == idb) .AND.(b>a)) CYCLE ! identical particles duplicate
          IF (idb > ida) CYCLE ! Combination allready counted
          IF ( na+na+la+nb+nb+lb > nlmax) CYCLE
          IF (istra+istrb /= istr) CYCLE
          IF (itzb+itza /= itz ) CYCLE
          IF ( (-1)**(la+lb+ipar) < 0) CYCLE
          IF ( triag( 2*ij, ja, jb) ) CYCLE
          IF ( (a == b) .AND. (MOD(ij,2) /= 0) ) CYCLE
          nconfs=nconfs+1
       ENDDO
    ENDDO
    this%number_confs=nconfs

  END SUBROUTINE number_gmatrix_confs
  !
  !
  !
  SUBROUTINE setup_gmatrix_configurations(ij,ipar,itz,istr, this)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, na, la, nb, lb, ja, jb, a, b, nconfs, &
         k1, k2, itza, itzb, istr, istra, istrb, ida, idb
    LOGICAL triag
    nconfs=0
    DO a = 1, all_orbit%total_orbits
       la=all_orbit%ll(a)
       ja=all_orbit%jj(a)
       na = all_orbit%nn(a)
       itza=all_orbit%itzp(a)
       istra=all_orbit%strange(a)
       ida = all_orbit%ptype(a)
       DO b = 1, all_orbit%total_orbits
          lb=all_orbit%ll(b)
          jb=all_orbit%jj(b)
          nb = all_orbit%nn(b)
          itzb=all_orbit%itzp(b)
          istrb=all_orbit%strange(b)
          idb = all_orbit%ptype(b)
          IF ( na+na+la+nb+nb+lb > nlmax) CYCLE
          IF ((all_orbit%model_space(a) == 'outside').AND.  &
               (all_orbit%model_space(b) == 'outside') ) CYCLE
          IF ((square_calculation).AND.(all_orbit%model_space(a) == 'inside').AND.  &
               (all_orbit%model_space(b) == 'outside') ) CYCLE
          IF ((square_calculation).AND.(all_orbit%model_space(a) == 'outside').AND.  &
               (all_orbit%model_space(b) == 'inside') ) CYCLE
          IF ( idb > ida ) CYCLE
          IF ((ida == idb) .AND.(b>a)) CYCLE ! identical particles duplicate
          IF (itzb+itza /= itz ) CYCLE
          IF (istra+istrb /= istr) CYCLE
          IF ( (-1)**(la+lb+ipar) < 0) CYCLE
          IF ( triag( 2*ij, ja, jb) ) CYCLE
          IF ( (a == b) .AND. (MOD(ij,2) /= 0) ) CYCLE
          nconfs=nconfs+1
          k2=nconfs*2
          k1=k2-1
          this%config_ab(k1)=b
          this%config_ab(k2)=a
       ENDDO
    ENDDO
    IF ( nconfs /= this%number_confs ) THEN
       WRITE(6,*) ' Error in configuration allocation ' ; STOP
    ENDIF
    CALL sort_configs(this,nconfs)
    SELECT CASE (istr)
    CASE (0)
        SELECT CASE (itz )
        CASE ( 2)
           WRITE(6,*) nconfs,' proton-proton configurations for J',ij
        CASE (-2)
           WRITE(6,*) nconfs,' neutron-neutron configurations for J ',ij
        CASE (0)
           WRITE(6,*) nconfs,' proton-neutron configurations for J ',ij
        END SELECT
    CASE (-1)
        SELECT CASE (itz)
        CASE (-3)
           WRITE(6,*) nconfs,' Sigma--neutron configurations for J',ij
        CASE (-1)
           WRITE(6,*) nconfs,' Lambda-neutron configurations for J',ij
        CASE (1)
           WRITE(6,*) nconfs,' Lambda-proton configurations for J',ij
        CASE (3)
           WRITE(6,*) nconfs,' Sigma+-proton configurations for J',ij
        END SELECT
    CASE (-2)
        SELECT CASE (itz)
        CASE (-4)
           WRITE(6,*) nconfs,' Sigma--Sigma- configurations for J',ij
        CASE (-2)
           WRITE(6,*) nconfs,' Lambda-Sigma- configurations for J',ij
        CASE (0)
           WRITE(6,*) nconfs,' Lambda-Sigma0 configurations for J',ij
        CASE (2)
           WRITE(6,*) nconfs,' Lambda-Sigma+ configurations for J',ij
        CASE (4)
           WRITE(6,*) nconfs,' Sigma+-Sigma+ configurations for J',ij
        END SELECT
    END SELECT

  END SUBROUTINE setup_gmatrix_configurations
  !                       
  !     setting up all configurations for given J, Tz and parity for
  !     the nocore triangular  model space
  !
  SUBROUTINE number_nocore_confs(ij,ipar,itz,this)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, na, la, nb, lb, ja, jb, a, b, b_end, nconfs, &
         itza, itzb
    LOGICAL triag
    nconfs=0
    DO a=1, all_orbit%total_orbits
       na = all_orbit%nn(a)
       la=all_orbit%ll(a)
       ja=all_orbit%jj(a)
       itza=all_orbit%itzp(a)
       IF ( itz /= 0 ) THEN
          b_end = a
       ELSE
          b_end = all_orbit%total_orbits
       ENDIF
       DO b=1, b_end
          nb = all_orbit%nn(b)
          lb=all_orbit%ll(b)
          itzb=all_orbit%itzp(b)
          jb=all_orbit%jj(b)
          IF ( na+na+la+nb+nb+lb > nlmax_model) CYCLE
          IF ( itzb > itza ) CYCLE
          IF (itzb+itza /= itz*2 ) CYCLE
          IF ( (-1)**(la+lb+ipar) < 0) CYCLE
          IF ( triag( 2*ij, ja, jb) ) CYCLE
          IF ( (a == b) .AND. (MOD(ij,2) /= 0) ) CYCLE
          nconfs=nconfs+1
       ENDDO
    ENDDO
    this%number_confs=nconfs

  END SUBROUTINE number_nocore_confs
  !
  !
  !
  SUBROUTINE setup_nocore_configurations(ij,ipar,itz,this)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, na, la, nb, lb, ja, jb, a, b, b_end, nconfs, &
         k1, k2, itza, itzb
    LOGICAL triag
    nconfs=0
    DO a = 1, all_orbit%total_orbits
       la=all_orbit%ll(a)
       ja=all_orbit%jj(a)
       na = all_orbit%nn(a)
       itza=all_orbit%itzp(a)
       IF ( itz /= 0 ) THEN
          b_end = a
       ELSE
          b_end = all_orbit%total_orbits
       ENDIF
       DO b = 1, b_end
          lb=all_orbit%ll(b)
          jb=all_orbit%jj(b)
          nb = all_orbit%nn(b)
          itzb=all_orbit%itzp(b)
          IF ( na+na+la+nb+nb+lb > nlmax_model) CYCLE
          IF ( itzb > itza ) CYCLE
          IF (itzb+itza /= itz*2 ) CYCLE
          IF ( (-1)**(la+lb+ipar) < 0) CYCLE
          IF ( triag( 2*ij, ja, jb) ) CYCLE
          IF ( (a == b) .AND. (MOD(ij,2) /= 0) ) CYCLE
          nconfs=nconfs+1
          k2=nconfs*2
          k1=k2-1
          this%config_ab(k1)=b
          this%config_ab(k2)=a
       ENDDO
    ENDDO
    IF ( nconfs /= this%number_confs ) THEN
       WRITE(6,*) ' Error in configuration allocation ' ; STOP
    ENDIF
    CALL sort_configs(this,nconfs)
    SELECT CASE (itz )
    CASE ( -1)
       WRITE(6,*) nconfs,' proton-proton configurations for J',ij
    CASE (1)
       WRITE(6,*) nconfs,' neutron-neutron configurations for J ',ij
    CASE (0)
       WRITE(6,*) nconfs,' proton-neutron configurations for J ',ij
    END SELECT

  END SUBROUTINE setup_nocore_configurations


  !
  !     setting up all configurations for given J, Tz and parity for
  !     the single particle model space
  !
  SUBROUTINE number_sp_confs(ij,ipar,itz,this)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, a, nconfs

    nconfs=0
    DO a=1, all_orbit%total_orbits
       IF ( MOD(all_orbit%ll(a),2) /= ipar ) CYCLE
       IF ( all_orbit%jj(a) /= ij) CYCLE
       IF ( all_orbit%itzp(a) /= itz ) CYCLE
       nconfs=nconfs+1
    ENDDO
    this%number_confs=nconfs

  END SUBROUTINE number_sp_confs
  !
  !     Allocates space for the sp configurations and sets them up
  !
  SUBROUTINE setup_sp_configurations(ij,ipar,itz,this)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: a, ij, ipar, itz, nconfs
    nconfs=0
    DO a = 1, all_orbit%total_orbits
       IF ( MOD(all_orbit%ll(a),2) /= ipar ) CYCLE
       IF ( all_orbit%jj(a) /= ij) CYCLE
       IF ( all_orbit%itzp(a) /= itz ) CYCLE
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
  !     and parity for the rel coordinates only by the
  !     number channel. First find the number of configurations.
  !
  SUBROUTINE number_configurations_rel(this)
    IMPLICIT NONE
    TYPE (configuration_relcm), INTENT(INOUT) :: this
    INTEGER :: l_rel, n_rel, a, nconfs, channel, l_min, l_max, nconfs_model

    DO channel=1, no_channels
       l_min=orb_lrel_min(channel)
       l_max=orb_lrel_max(channel)
       nconfs=0; nconfs_model = 0
       DO a=1, relcm_sp_data%max_value
          n_rel=relcm_sp_data%nrel(a)
          l_rel=relcm_sp_data%lrel(a)
          IF ( (l_rel /= l_max) .AND. ( l_rel /= l_min) ) CYCLE
          IF (2*n_rel+l_rel > nlmax) CYCLE
          IF ( 2*n_rel+l_rel <= nlmax_model) THEN
             nconfs_model = nconfs_model +1
          ENDIF
          nconfs=nconfs+1
       ENDDO
       this%nconfs_rel(channel)=nconfs              
       this%nconfs_relmodel(channel)=nconfs_model              
    ENDDO

  END SUBROUTINE number_configurations_rel
  !
  !     Now we find the possible configurations
  !
  SUBROUTINE setup_configurations_rel(this)
    IMPLICIT NONE
    TYPE (configuration_relcm), INTENT(INOUT) :: this
    INTEGER :: l_rel, n_rel, a, channel, l_min, l_max, nconfs

    DO channel=1, no_channels
       l_min=orb_lrel_min(channel)
       l_max=orb_lrel_max(channel)
       nconfs=0
       DO a=1, relcm_sp_data%max_value
          n_rel=relcm_sp_data%nrel(a)
          l_rel=relcm_sp_data%lrel(a)
          IF ( (l_rel /= l_max) .AND. ( l_rel /= l_min ) ) CYCLE
          IF ( 2*n_rel+l_rel > nlmax) CYCLE
          nconfs=nconfs+1
          this%rel_ab(channel,nconfs)=a
       ENDDO
       IF ( nconfs /= this%nconfs_rel(channel) ) THEN
          WRITE(6,*) ' Error in configuration allocation ' ; STOP
       ENDIF
    ENDDO

  END SUBROUTINE setup_configurations_rel
  !
  !     setting up all configurations for given J, Tz, spin
  !     and parity for the rel and cm system specified by the
  !     number channel. First find the number of configurations.
  !
  SUBROUTINE number_configurations_relcm(this)
    IMPLICIT NONE
    TYPE (configuration_relcm), INTENT(INOUT) :: this
    INTEGER :: l_rel, l_cm, n_rel, n_cm, a, b, &
         nconfs, channel, l_min, l_max

    DO channel=1, no_channels
       l_min=orb_lrel_min(channel)
       l_max=orb_lrel_max(channel)
       nconfs=0
       DO a=1, relcm_sp_data%max_value
          n_rel=relcm_sp_data%nrel(a)
          l_rel=relcm_sp_data%lrel(a)
          IF ( (l_rel /= l_max) .AND. ( l_rel /= l_min) ) CYCLE
          DO b=1, relcm_sp_data%max_value
             n_cm=relcm_sp_data%nrel(b)
             l_cm=relcm_sp_data%lrel(b)
             IF ( 2*n_cm+l_cm+2*n_rel+l_rel > nlmax) CYCLE
             nconfs=nconfs+1
          ENDDO
       ENDDO
       this%nconfs_relcm(channel)=nconfs              
    ENDDO

  END SUBROUTINE number_configurations_relcm
  !
  !     Now we find the possible configurations
  !
  SUBROUTINE setup_configurations_relcm(this)
    IMPLICIT NONE
    TYPE (configuration_relcm), INTENT(INOUT) :: this
    INTEGER :: l_rel, l_cm, n_rel, n_cm, a, b, &
         channel, l_min, l_max, k1, k2, nconfs

    DO channel=1, no_channels
       l_min=orb_lrel_min(channel)
       l_max=orb_lrel_max(channel)
       nconfs=0
       DO a=1, relcm_sp_data%max_value
          n_rel=relcm_sp_data%nrel(a)
          l_rel=relcm_sp_data%lrel(a)
          IF ( (l_rel /= l_max) .AND. ( l_rel /= l_min ) ) CYCLE
          DO b=1, relcm_sp_data%max_value
             n_cm=relcm_sp_data%nrel(b)
             l_cm=relcm_sp_data%lrel(b)
             IF ( 2*n_cm+l_cm+2*n_rel+l_rel > nlmax) CYCLE
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
    ENDDO

  END SUBROUTINE setup_configurations_relcm

  !
  !     Finding all possible values for rel and cm variables n, N, l, L
  !
  SUBROUTINE find_spdata_relcm(a)
    IMPLICIT NONE
    INTEGER , INTENT(OUT) :: a
    INTEGER :: n, l

    a=0
    DO n=0, nmax
       DO l=0, lmax
          a=a+1
       ENDDO
    ENDDO

  END SUBROUTINE find_spdata_relcm
  !
  !     setting up all quantum numbers for relative and
  !     c.m coordinates in the block relcm_sp_data
  !
  SUBROUTINE make_configurations_relcm
    IMPLICIT NONE
    INTEGER :: a, n, l

    a=0
    DO n=0, nmax
       DO l=0, lmax
          a=a+1
          relcm_sp_data%nrel(a)=n
          relcm_sp_data%lrel(a)=l
       ENDDO
    ENDDO

  END SUBROUTINE make_configurations_relcm

END MODULE configurations
!
!     This module contains data on the g-matrix or screened interaction
!     The data are read in by the function read_g_matrix
!
MODULE stored_bare_interaction
  USE constants
  USE single_particle_orbits
  TYPE, PUBLIC :: jcoupled_gmatrix_storage
     SEQUENCE
     !     INTEGER :: amount
     REAL(DP),         DIMENSION(:,:), ALLOCATABLE :: mtxel
     CHARACTER(LEN=10), DIMENSION(:),   ALLOCATABLE :: mtxel_label
  END TYPE jcoupled_gmatrix_storage
  ! NUmber of elements of different types
  INTEGER, PUBLIC :: number_twobody_elements, nn_elements, yn_elements,&
  yy_elements
  INTEGER, PUBLIC :: number_pp_elements, number_pn_elements, number_nn_elements
  INTEGER, PUBLIC :: number_Spp_elements, number_Lp_elements, &
    number_Ln_elements, number_Smn_elements
  INTEGER, PUBLIC :: number_SpSp_elements, number_LSp_elements, &
    number_LS0_elements, number_LSm_elements, number_SmSm_elements

  TYPE(jcoupled_gmatrix_storage), TARGET, PRIVATE :: gmatrix_pppp, gmatrix_pnpn, gmatrix_nnnn, &
    gmatrix_Spp, gmatrix_Lp, gmatrix_Ln, gmatrix_Smn, gmatrix_SpSp, &
    gmatrix_LSp, gmatrix_LS0, gmatrix_LSm, gmatrix_SmSm

CONTAINS

  SUBROUTINE read_gmatrix_nn
    INTEGER :: amount
    INTEGER, TARGET :: amount_pppp, amount_pnpn, amount_nnnn
    INTEGER, POINTER :: mtxel_number
    INTEGER :: strange, t_z, jtot, PARITY, a, b, c, d, i, ii
    TYPE(jcoupled_gmatrix_storage), POINTER :: gmatrix_temp
    REAL(DP), DIMENSION(1:n_startenergy_g+3) :: gmtxel
    REAL(DP), DIMENSION(1:n_startenergy_g+3) :: gmat

    !WRITE(*,*) "read_gmatrix_nn called..."
    amount = nn_elements
    amount_pppp = number_pp_elements
    amount_pnpn = number_pn_elements
    amount_nnnn = number_nn_elements
    IF ((amount_pppp + amount_pnpn + amount_nnnn) /= amount) THEN
       WRITE(6,*) 'Error in  amount of NN G-matrix elements:', &
            amount_pppp + amount_pnpn + amount_nnnn, amount
    ENDIF
    WRITE (6,*) 'amount of J-coupled NN G-matrix elements:', &
         amount_pppp, amount_pnpn, amount_nnnn, amount
    IF( amount == 0) RETURN
    !construct storage for the J-coupled gmatrix
    CALL construct_gmatrix_storage(gmatrix_pppp, amount_pppp)
    CALL construct_gmatrix_storage(gmatrix_pnpn, amount_pnpn)
    CALL construct_gmatrix_storage(gmatrix_nnnn, amount_nnnn)
    amount_pppp = 0; amount_pnpn = 0; amount_nnnn = 0
    DO ii = 1, amount
       READ(7,*) strange, t_z, PARITY, jtot, a, b, c, d, (gmat(i), i=1,n_startenergy_g+3)
       gmtxel = gmat
       SELECT CASE (type_of_renormv)
       CASE ('vlowk')
          gmtxel(1:n_startenergy_g) =  gmtxel(1:n_startenergy_g)- &
             gmtxel(n_startenergy_g+3)/total_mass
       END SELECT
       NULLIFY (gmatrix_temp, mtxel_number) 
       SELECT CASE (t_z)
       CASE(2)
          gmatrix_temp => gmatrix_pppp; mtxel_number => amount_pppp
       CASE( 0)
          gmatrix_temp => gmatrix_pnpn; mtxel_number => amount_pnpn
       CASE( -2)
          gmatrix_temp => gmatrix_nnnn; mtxel_number => amount_nnnn
       CASE DEFAULT
          WRITE (6,*) 'unknown t_z value in associate_jcoupled_gmatrix_storage.'; STOP
       END SELECT ! t_z 
       jtot = jtot / 2
       mtxel_number = mtxel_number + 1
       gmatrix_temp%mtxel(mtxel_number, :) = gmtxel(:)
       gmatrix_temp%mtxel_label(mtxel_number) = &
            construct_jcoupled_gmatrix_label(PARITY, jtot, a, b, c, d)
       !check if the strings are in ascending order in ASCII sense.
       !necessary for the search of m-scheme G-matrix element to operate correctly 
       IF (mtxel_number /= 1) THEN
          IF (gmatrix_temp%mtxel_label(mtxel_number) <= &                
               gmatrix_temp%mtxel_label(mtxel_number-1)) THEN
             WRITE (6,*) 'Ascending ordering failed for mtxel_label in read_gmatrix_new', t_z, mtxel_number
             STOP
          ENDIF
       ENDIF
    ENDDO
    IF ((amount_pppp + amount_pnpn + amount_nnnn) /= amount) THEN
       WRITE(6,*) 'read_gmatrix_new: read wrong amount of G-matrix elements.'
       WRITE(6,*) 'Expected:', amount, ' read:', amount_pppp + amount_pnpn + amount_nnnn
    ENDIF

  END SUBROUTINE read_gmatrix_nn

  SUBROUTINE read_gmatrix_yn
    INTEGER :: amount
    INTEGER, TARGET :: amount_Spp, amount_Lp, amount_Ln, amount_Smn
    INTEGER, POINTER :: mtxel_number
    INTEGER :: strange, t_z, jtot, PARITY, a, b, c, d, i, ii
    TYPE(jcoupled_gmatrix_storage), POINTER :: gmatrix_temp
    REAL(DP), DIMENSION(1:n_startenergy_g+3) :: gmtxel
    REAL(DP), DIMENSION(1:n_startenergy_g+3) :: gmat

    amount = yn_elements
    amount_Spp = number_Spp_elements
    amount_Lp = number_Lp_elements
    amount_Ln = number_Ln_elements
    amount_Smn = number_Smn_elements
    IF ((amount_Spp + amount_Lp + amount_Ln + amount_Smn) /= amount) THEN
       WRITE(6,*) 'Error in  amount of YN G-matrix elements:', &
            amount_Spp + amount_Lp + amount_Ln + amount_Smn, amount
    ENDIF
    WRITE (6,*) 'amount of J-coupled YN G-matrix elements:', &
         amount_Spp, amount_Lp, amount_Ln, amount_Smn, amount
    IF( amount == 0) RETURN
    !construct storage for the J-coupled gmatrix
    CALL construct_gmatrix_storage(gmatrix_Spp, amount_Spp)
    CALL construct_gmatrix_storage(gmatrix_Lp, amount_Lp)
    CALL construct_gmatrix_storage(gmatrix_Ln, amount_Ln)
    CALL construct_gmatrix_storage(gmatrix_Smn, amount_Smn)
    amount_Spp = 0; amount_Lp = 0; amount_Ln = 0; amount_Smn = 0
    DO ii = 1, amount
       READ(7,*) strange, t_z, PARITY, jtot, a, b, c, d, (gmat(i), i=1,n_startenergy_g+3)
       gmtxel = gmat
       SELECT CASE (type_of_renormv)
       CASE ('vlowk')
          gmtxel(1:n_startenergy_g) =  gmtxel(1:n_startenergy_g)- &
             gmtxel(n_startenergy_g+3)/total_mass
       END SELECT
       NULLIFY (gmatrix_temp, mtxel_number) 
       SELECT CASE (t_z)
       CASE(3)
          gmatrix_temp => gmatrix_Spp; mtxel_number => amount_Spp
       CASE( 1)
          gmatrix_temp => gmatrix_Lp; mtxel_number => amount_Lp
       CASE( -1)
          gmatrix_temp => gmatrix_Ln; mtxel_number => amount_Ln
       CASE( -3)
          gmatrix_temp => gmatrix_Smn; mtxel_number => amount_Smn
       CASE DEFAULT
          WRITE (6,*) 'unknown t_z value in associate_jcoupled_gmatrix_storage.'; STOP
       END SELECT ! t_z 
       jtot = jtot / 2
       mtxel_number = mtxel_number + 1
       gmatrix_temp%mtxel(mtxel_number, :) = gmtxel(:)
       gmatrix_temp%mtxel_label(mtxel_number) = &
            construct_jcoupled_gmatrix_label(PARITY, jtot, a, b, c, d)
       !check if the strings are in ascending order in ASCII sense.
       !necessary for the search of m-scheme G-matrix element to operate correctly 
       IF (mtxel_number /= 1) THEN
          IF (gmatrix_temp%mtxel_label(mtxel_number) <= &                
               gmatrix_temp%mtxel_label(mtxel_number-1)) THEN
             WRITE (6,*) 'Ascending ordering failed for mtxel_label in read_gmatrix_new', strange, t_z, mtxel_number
             !          STOP
          ENDIF
       ENDIF
    ENDDO
    IF ((amount_Spp + amount_Lp + amount_Ln + amount_Smn) /= amount) THEN
       WRITE(6,*) 'read_gmatrix_new: read wrong amount of YN G-matrix elements.'
       WRITE(6,*) 'Expected:', amount, ' read:', amount_Spp + amount_Lp + &
            amount_Ln + amount_Smn
    ENDIF

  END SUBROUTINE read_gmatrix_yn

  SUBROUTINE read_gmatrix_yy
    INTEGER :: amount
    INTEGER, TARGET :: amount_SpSp, amount_LSp, amount_LS0, amount_LSm, &
        amount_SmSm
    INTEGER, POINTER :: mtxel_number
    INTEGER :: strange, t_z, jtot, PARITY, a, b, c, d, i, ii
    TYPE(jcoupled_gmatrix_storage), POINTER :: gmatrix_temp
    REAL(DP), DIMENSION(1:n_startenergy_g+3) :: gmtxel
    REAL(DP), DIMENSION(1:n_startenergy_g+3) :: gmat

    amount = yy_elements
    amount_SpSp = number_SpSp_elements
    amount_LSp = number_LSp_elements
    amount_LS0 = number_LS0_elements
    amount_LSm = number_LSm_elements
    amount_SmSm = number_SmSm_elements
    IF ((amount_SpSp + amount_LSp + amount_LS0 + amount_LSm+ amount_SmSm) &
            /= amount) THEN
       WRITE(6,*) 'Error in  amount of YY G-matrix elements:', &
            amount_SpSp + amount_LSp + amount_LS0 + amount_LSm+ amount_SmSm, amount
    ENDIF
    WRITE (6,*) 'amount of J-coupled YY G-matrix elements:', &
         amount_SpSp, amount_LSp, amount_LS0, amount_LSm, amount_SmSm, amount
    IF( amount == 0) RETURN
    !construct storage for the J-coupled gmatrix
    CALL construct_gmatrix_storage(gmatrix_SpSp, amount_SpSp)
    CALL construct_gmatrix_storage(gmatrix_LSp, amount_LSp)
    CALL construct_gmatrix_storage(gmatrix_LS0, amount_LS0)
    CALL construct_gmatrix_storage(gmatrix_LSm, amount_LSm)
    CALL construct_gmatrix_storage(gmatrix_SmSm, amount_SmSm)
    amount_SpSp=0;amount_LSp=0;amount_LS0=0;amount_LSm=0;amount_SmSm=0
    DO ii = 1, amount
       READ(7,*) strange, t_z, PARITY, jtot, a, b, c, d, (gmat(i), i=1,n_startenergy_g+3)
       gmtxel = gmat
       SELECT CASE (type_of_renormv)
       CASE ('vlowk')
          gmtxel(1:n_startenergy_g) =  gmtxel(1:n_startenergy_g)- &
            gmtxel(n_startenergy_g+3)/total_mass
       END SELECT
       NULLIFY (gmatrix_temp, mtxel_number) 
       SELECT CASE (t_z)
       CASE(4)
          gmatrix_temp => gmatrix_SpSp; mtxel_number => amount_SpSp
       CASE( 2)
          gmatrix_temp => gmatrix_LSp; mtxel_number => amount_LSp
       CASE( 0)
          gmatrix_temp => gmatrix_LS0; mtxel_number => amount_LS0
       CASE( -2)
          gmatrix_temp => gmatrix_LSm; mtxel_number => amount_LSm
       CASE( -4)
          gmatrix_temp => gmatrix_SmSm; mtxel_number => amount_SmSm
       CASE DEFAULT
          WRITE (6,*) 'unknown t_z value in associate_jcoupled_gmatrix_storage.'; STOP
       END SELECT ! t_z 
       jtot = jtot / 2
       mtxel_number = mtxel_number + 1
       gmatrix_temp%mtxel(mtxel_number, :) = gmtxel(:)
       gmatrix_temp%mtxel_label(mtxel_number) = &
            construct_jcoupled_gmatrix_label(PARITY, jtot, a, b, c, d)
       !check if the strings are in ascending order in ASCII sense.
       !necessary for the search of m-scheme G-matrix element to operate correctly 
       IF (mtxel_number /= 1) THEN
          IF (gmatrix_temp%mtxel_label(mtxel_number) <= &                
               gmatrix_temp%mtxel_label(mtxel_number-1)) THEN
             WRITE (6,*) 'Ascending ordering failed for mtxel_label in read_gmatrix_new', strange, t_z, mtxel_number
             !          STOP
          ENDIF
       ENDIF
    ENDDO
    IF ((amount_SpSp + amount_LSp + amount_LS0 + amount_LSm+ amount_SmSm) /= amount) THEN
       WRITE(6,*) 'read_gmatrix_new: read wrong amount of YY G-matrix elements.'
       WRITE(6,*) 'Expected:', amount, ' read:', amount_SpSp + amount_LSp + &
         amount_LS0 + amount_LSm+ amount_SmSm
    ENDIF

  END SUBROUTINE read_gmatrix_yy


  SUBROUTINE construct_gmatrix_storage (gmatrix_temp, amount)
    INTEGER, INTENT(IN) :: amount
    TYPE (jcoupled_gmatrix_storage), INTENT(INOUT) :: gmatrix_temp

    !check if arrays has not been already allocated
    IF (ALLOCATED (gmatrix_temp%mtxel) ) DEALLOCATE(gmatrix_temp%mtxel)
    IF (ALLOCATED (gmatrix_temp%mtxel_label) ) DEALLOCATE(gmatrix_temp%mtxel_label)
    ALLOCATE (gmatrix_temp%mtxel(1:amount, 1:n_startenergy_g+3), &
         gmatrix_temp%mtxel_label(1:amount))

  END SUBROUTINE construct_gmatrix_storage


  FUNCTION construct_jcoupled_gmatrix_label(PARITY, jtot, a, b, c, d)
    USE constants
    USE configurations
    IMPLICIT NONE
    CHARACTER(LEN=10) :: construct_jcoupled_gmatrix_label
    INTEGER, INTENT (IN) :: PARITY, jtot, a, b, c, d
    INTEGER, PARAMETER :: ascii_limit = 255
    INTEGER :: a_i, b_i, c_i, d_i, a_f, b_f, c_f, d_f

    IF (jtot < 0)   STOP 'negative jtot in construct_jcoupled_gmatrix_label'
    IF (PARITY < 0) STOP 'negative parity in construct_jcoupled_gmatrix_label'
    IF (a < 0)      STOP 'negative a in construct_jcoupled_gmatrix_label'
    IF (b < 0)      STOP 'negative b in construct_jcoupled_gmatrix_label'
    IF (c < 0)      STOP 'negative c in construct_jcoupled_gmatrix_label'
    IF (d < 0)      STOP 'negative d in construct_jcoupled_gmatrix_label'
    a_i = a / ascii_limit;  a_f = a - a_i * ascii_limit
    b_i = b / ascii_limit;  b_f = b - a_i * ascii_limit
    c_i = c / ascii_limit;  c_f = c - a_i * ascii_limit
    d_i = d / ascii_limit;  d_f = d - a_i * ascii_limit
    construct_jcoupled_gmatrix_label = ACHAR(PARITY)//ACHAR(jtot)//&
         ACHAR(a_i)//ACHAR(a_f)//ACHAR(b_i)//ACHAR(b_f)// &
         ACHAR(c_i)//ACHAR(c_f)//ACHAR(d_i)//ACHAR(d_f)

  END FUNCTION construct_jcoupled_gmatrix_label
  !      
  !   finds j-coupled G-matrix element given \pi, J, and abcd orbit identifiers
  !   
  SUBROUTINE mtx_elements (a, b, c, d, jtot, gmtxel)
    IMPLICIT NONE
    LOGICAL triag
    INTEGER, INTENT (IN) :: a, b, c, d, jtot
    REAL(DP), DIMENSION (n_startenergy_g+3), INTENT (INOUT) :: gmtxel
    REAL(DP) :: fnorm, dij, fnorm_ab, fnorm_cd
    INTEGER  :: isospin_z, PARITY, local_a, local_b, local_c, local_d, &
        tot_strange
    INTEGER  :: tz_a, tz_b, tz_c, tz_d, stra, strb, strc, strd, &
        ida, idb, idc, idd
    INTEGER  :: low, high, mid, checker, phase_ab, phase_cd, iph
    CHARACTER(LEN=10) :: search_pattern
    TYPE (jcoupled_gmatrix_storage), POINTER :: gmatrix_temp

    gmtxel(:) = 0.0_dp

    tz_a = all_orbit%itzp(a) ; tz_b = all_orbit%itzp(b)
    tz_c = all_orbit%itzp(c) ; tz_d = all_orbit%itzp(d)
    stra = all_orbit%strange(a); strb = all_orbit%strange(b)
    strc = all_orbit%strange(c); strd = all_orbit%strange(d)
    ida = all_orbit%ptype(a); idb = all_orbit%ptype(b)
    idc = all_orbit%ptype(c); idd = all_orbit%ptype(d)

    IF ( (stra + strb) /= (strc + strd) ) RETURN
    IF ( (tz_a + tz_b) /= (tz_c + tz_d) ) RETURN
    IF(iph(all_orbit%ll(a) + all_orbit%ll(b)) /=  &
         iph(all_orbit%ll(c) + all_orbit%ll(d))) RETURN
    IF((a == b).AND.(MOD(jtot,2)/=0)) RETURN       
    IF((c == d).AND.(MOD(jtot,2)/=0)) RETURN       
    IF(triag(all_orbit%jj(a),all_orbit%jj(b),2*jtot)) RETURN
    IF(triag(all_orbit%jj(c),all_orbit%jj(d),2*jtot)) RETURN
    PARITY = MOD((all_orbit%ll(a) + all_orbit%ll(b)),2)

    !set values of local variables
    isospin_z = tz_a + tz_b
    phase_ab = 0; phase_cd = 0
    local_a = a; local_b = b; local_c = c; local_d = d
    tot_strange = stra + strb

    !select isospin configuration and swap local orbit indexes if necessary
    !For ip a<b, and ida < idb for not
    IF ((ida == idb) .AND. (local_a > local_b)) THEN 
       CALL SWAP (local_a, local_b); phase_ab = (all_orbit%jj(a) + all_orbit%jj(b))/2 - jtot + 1
    ENDIF
    IF ((idc == idd) .AND. (local_c > local_d)) THEN 
       CALL SWAP (local_c, local_d); phase_cd = (all_orbit%jj(c) + all_orbit%jj(d))/2 - jtot + 1
    ENDIF
    IF (ida > idb) THEN
       CALL SWAP (local_a, local_b); phase_ab = (all_orbit%jj(a) + all_orbit%jj(b))/2 - jtot + 1
    ENDIF
    IF (idc > idd) THEN
       CALL SWAP (local_c, local_d); phase_cd = (all_orbit%jj(c) + all_orbit%jj(d))/2 - jtot + 1
    ENDIF

    NULLIFY (gmatrix_temp)
    SELECT CASE (tot_strange)
    CASE (0)
      SELECT CASE (isospin_z)
      CASE (2) 
        gmatrix_temp => gmatrix_pppp
      CASE (0)
        !maintain pn-ordering
        gmatrix_temp => gmatrix_pnpn
      CASE (-2)
        gmatrix_temp => gmatrix_nnnn
      CASE DEFAULT
        WRITE (6,*) 'jcoupled_gmatrix_element:'
        WRITE (6,*) 'isospin_z value is not recognized:', isospin_z; STOP
      END SELECT
    CASE (-1)
      SELECT CASE (isospin_z)
      CASE (3)
        gmatrix_temp => gmatrix_Spp
      CASE (1)
        gmatrix_temp => gmatrix_Lp
      CASE (-1)
        gmatrix_temp => gmatrix_Ln
      CASE (-3)
        gmatrix_temp => gmatrix_Smn
      END SELECT
    CASE (-2)
      SELECT CASE (isospin_z)
      CASE (4)
        gmatrix_temp => gmatrix_SpSp
      CASE (2)
        gmatrix_temp => gmatrix_LSp
      CASE (0) 
        gmatrix_temp => gmatrix_LS0
      CASE (-2)
        gmatrix_temp => gmatrix_LSm
      CASE (-4)
        gmatrix_temp => gmatrix_SmSm
      END SELECT
    END SELECT

    !initialize value of CHARACTER search pattern
    !swap bra <-> ket if necessary
    IF (ACHAR(local_a)//ACHAR(local_b) <= ACHAR(local_c)//ACHAR(local_d)) THEN
       search_pattern = construct_jcoupled_gmatrix_label(PARITY, jtot, local_a, local_b, local_c, local_d)
    ELSE
       search_pattern = construct_jcoupled_gmatrix_label(PARITY, jtot, local_c, local_d, local_a, local_b)
    ENDIF

    !set values of search variables
    low = 1; high = SIZE(gmatrix_temp%mtxel_label); checker = 0
    !find matrix element address in gmatrix_temp_mscheme%mtxel(:,:) array
    DO WHILE(checker <= SIZE(gmatrix_temp%mtxel_label))
       mid = (low + high)/2
       checker = checker + 1
       IF(search_pattern < gmatrix_temp%mtxel_label(mid)) THEN
          high = mid - 1
       ELSEIF(search_pattern > gmatrix_temp%mtxel_label(mid)) THEN
          low = mid + 1
       ELSEIF(search_pattern == gmatrix_temp%mtxel_label(mid)) THEN
          !found matrix element address, set the return value
          fnorm_ab = 1.d0
          fnorm_cd = 1.d0
          IF ( ida == idb ) THEN
             fnorm_ab = dij(local_a,local_b)
          ENDIF
          IF ( idc == idd ) THEN
             fnorm_cd = dij(local_c, local_d)
          ENDIF
          IF ((ida /= idb) .AND. (idc /= idd)) THEN
             fnorm = 1.d0
          ELSE
             fnorm = fnorm_ab * fnorm_cd
          ENDIF
          !WRITE(6,*) "fnorm: ", fnorm, gmatrix_temp%mtxel(mid,1:n_startenergy_g+3)
          gmtxel(1:n_startenergy_g+3) = gmatrix_temp%mtxel(mid,1:n_startenergy_g+3) * iph( phase_ab + phase_cd ) * fnorm
          !WRITE(6,*) fnorm * gmatrix_temp%mtxel(mid,1:n_startenergy_g+3)
          RETURN
       ENDIF
    ENDDO
    WRITE (*,*) "Could not find matrix element...", a, b, c, d, &
        PARITY, jtot, isospin_z, tot_strange

  END SUBROUTINE mtx_elements
  !      
  !   finds j-coupled G-matrix element given \pi, J, and abcd orbit identifiers
  !   
  SUBROUTINE replace_g(a, b, c, d, istr,it,ip,jtot, gmtxel)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: a, b, c, d, it,ip,jtot, istr
    REAL(DP), DIMENSION (n_startenergy_g+3), INTENT (IN) :: gmtxel
    INTEGER  :: local_a, local_b, local_c, local_d
    INTEGER  :: low, high, mid, checker
    CHARACTER(LEN=10) :: search_pattern
    TYPE (jcoupled_gmatrix_storage), POINTER :: gmatrix_temp

    local_a = a; local_b = b; local_c = c; local_d = d

    NULLIFY (gmatrix_temp)
    SELECT CASE (istr)
    CASE (0)
      SELECT CASE (it)
      CASE (2) 
        gmatrix_temp => gmatrix_pppp
      CASE (0)
        !maintain pn-ordering
        gmatrix_temp => gmatrix_pnpn
      CASE (-2)
        gmatrix_temp => gmatrix_nnnn
      CASE DEFAULT
        WRITE (6,*) 'jcoupled_gmatrix_element:'
        WRITE (6,*) 'isospin_z value is not recognized:', it; STOP
      END SELECT
    CASE (-1)
      SELECT CASE (it)
      CASE (3)
        gmatrix_temp => gmatrix_Spp
      CASE (1)
        gmatrix_temp => gmatrix_Lp
      CASE (-1)
        gmatrix_temp => gmatrix_Ln
      CASE (-3)
        gmatrix_temp => gmatrix_Smn
      END SELECT
    CASE (-2)
      SELECT CASE (it)
      CASE (4)
        gmatrix_temp => gmatrix_SpSp
      CASE (2)
        gmatrix_temp => gmatrix_LSp
      CASE (0) 
        gmatrix_temp => gmatrix_LS0
      CASE (-2)
        gmatrix_temp => gmatrix_LSm
      CASE (-4)
        gmatrix_temp => gmatrix_SmSm
      END SELECT
    END SELECT
    !initialize value of CHARACTER search pattern
    !swap bra <-> ket if necessary
    IF (ACHAR(local_a)//ACHAR(local_b) <= ACHAR(local_c)//ACHAR(local_d)) THEN
       search_pattern = construct_jcoupled_gmatrix_label(ip, jtot, local_a, local_b, local_c, local_d)
    ELSE
       search_pattern = construct_jcoupled_gmatrix_label(ip, jtot, local_c, local_d, local_a, local_b)
    ENDIF
    !set values of search variables
    low = 1; high = SIZE(gmatrix_temp%mtxel_label); checker = 0
    !find matrix element address in gmatrix_temp_mscheme%mtxel(:,:) array
    DO WHILE(checker <= SIZE(gmatrix_temp%mtxel_label))
       mid = (low + high)/2
       checker = checker + 1
       IF(search_pattern < gmatrix_temp%mtxel_label(mid)) THEN
          high = mid - 1
       ELSEIF(search_pattern > gmatrix_temp%mtxel_label(mid)) THEN
          low = mid + 1
       ELSEIF(search_pattern == gmatrix_temp%mtxel_label(mid)) THEN
          !found matrix element address, set the return value
          gmatrix_temp%mtxel(mid,:)  = 0.0_dp
          gmatrix_temp%mtxel(mid,1:n_startenergy_g+3)  = gmtxel(1:n_startenergy_g+3)
          RETURN
       ENDIF
    ENDDO

  END SUBROUTINE replace_g

END MODULE  stored_bare_interaction
!
!     Function to calculate phase factors    
!
INTEGER FUNCTION iph(n)
  IMPLICIT NONE
  INTEGER :: n
  iph=(-1)**n
END FUNCTION iph
!
!     Function to check # osc. excitations   
! 
LOGICAL FUNCTION dencheck(i)
  USE constants
  IMPLICIT NONE
  INTEGER :: i 
  DENCHECK = ((ABS(i) == 0).OR.(ABS(i) > number_homega_exct))

END FUNCTION dencheck
!
!     Function to check triangular relations     
!
LOGICAL FUNCTION triag(i,j,k)
  IMPLICIT NONE
  INTEGER :: i, j, k
  triag = ((i-j-k)*(i-ABS(j-k)) > 0)

END FUNCTION triag
!
!      Function to calculate norm of g-mat    
!
REAL(KIND=8) FUNCTION dij(ja,jb)
  USE constants
  IMPLICIT NONE
  INTEGER :: ja, jb
  IF(ja == jb ) THEN
     dij=SQRT(2.0_dp)
  ELSE
     dij=1.0_dp
  ENDIF

END FUNCTION dij
!
!      Function to calculate norm of g-mat
!
REAL(DP) FUNCTION delta(ja,jb)
  USE constants
  IMPLICIT NONE
  INTEGER :: ja, jb
  IF(ja == jb ) THEN
     delta=1.d0
  ELSE
     delta=0.d0
  ENDIF

END FUNCTION delta
!
!     swaps values of 2 integers
!
SUBROUTINE SWAP(a, b)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: a, b
  INTEGER :: c

  c = a; a = b; b = c

END SUBROUTINE SWAP
!    Block which contains matrix operations
!
!     Eigenvalues and eigenvectors from Jacobi method
!     a is matrix to be diagonalized, n and np the dimension
!     of the system, d the eigenvalues and v the eigenvectors
!
SUBROUTINE matrix_diag(a,n,np,d,v,nrot)
  USE constants
  IMPLICIT NONE
  INTEGER :: n,np,nrot,n_max
  REAL(DP) ::  a(np,np),d(np),v(np,np)
  PARAMETER (n_max=2000)
  INTEGER :: i,ip,iq,j
  REAL(DP) ::  c,g,h,s,sm,t,tau,theta,tresh,b(n_max),z(n_max)

  DO ip=1,n
     DO iq=1,n
        v(ip,iq)=0.
     ENDDO
     v(ip,ip)=1.
  ENDDO
  DO ip=1,n
     b(ip)=a(ip,ip)
     d(ip)=b(ip)
     z(ip)=0.
  ENDDO
  nrot=0
  DO i=1,50
     sm=0.
     DO ip=1,n-1
        DO iq=ip+1,n
           sm=sm+abs(a(ip,iq))
        ENDDO
     ENDDO
     IF (sm == 0.) RETURN
     IF (i < 4) THEN
        tresh=0.2*sm/n**2
     ELSE
        tresh=0.
     ENDIF
     DO ip=1,n-1
        DO iq=ip+1,n
           g=100.*ABS(a(ip,iq))
           IF ((i > 4).and.(ABS(d(ip))+ &
                g == ABS(d(ip))).and.(ABS(d(iq))+g == ABS(d(iq)))) THEN
              a(ip,iq)=0.
           ELSEIF (ABS (a(ip,iq)) > tresh) THEN
              h=d(iq)-d(ip)
              IF (ABS(h)+g == ABS(h)) THEN
                 t=a(ip,iq)/h
              ELSE
                 theta=0.5*h/a(ip,iq)
                 t=1./(ABS(theta)+SQRT(1.+theta**2))
                 IF (theta < 0.)t=-t
              ENDIF
              c=1./SQRT(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              DO j=1,ip-1
                 g=a(j,ip)
                 h=a(j,iq)
                 a(j,ip)=g-s*(h+g*tau)
                 a(j,iq)=h+s*(g-h*tau)
              ENDDO
              DO j=ip+1,iq-1
                 g=a(ip,j)
                 h=a(j,iq)
                 a(ip,j)=g-s*(h+g*tau)
                 a(j,iq)=h+s*(g-h*tau)
              ENDDO
              DO j=iq+1,n
                 g=a(ip,j)
                 h=a(iq,j)
                 a(ip,j)=g-s*(h+g*tau)
                 a(iq,j)=h+s*(g-h*tau)
              ENDDO
              DO j=1,n
                 g=v(j,ip)
                 h=v(j,iq)
                 v(j,ip)=g-s*(h+g*tau)
                 v(j,iq)=h+s*(g-h*tau)
              ENDDO
              nrot=nrot+1
           ENDIF
        ENDDO
     ENDDO
     DO ip=1,n
        b(ip)=b(ip)+z(ip)
        d(ip)=b(ip)
        z(ip)=0.
     ENDDO
  ENDDO

END SUBROUTINE matrix_diag
!
!            Routines to do mtx inversion, from Numerical
!            Recepies, Teukolsky et al. Routines included
!            below are MATINV, LUDCMP and LUBKSB. See chap 2
!            of Numerical Recipes for further details
!            Recoded in FORTRAN 90 by M. Hjorth-Jensen
!
SUBROUTINE matinv(a,n)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  INTEGER :: i, j
  REAL(DP), DIMENSION(n,n), INTENT(INOUT)  :: a
  REAL(DP), ALLOCATABLE :: y(:,:)
  REAL(DP) :: d
  INTEGER, ALLOCATABLE :: indx(:)

  ALLOCATE (y( n, n))  ; ALLOCATE ( indx (n))
  y=0.
  !     setup identity matrix
  DO i=1,n
     y(i,i)=1.
  ENDDO
  !     LU decompose the matrix just once
  CALL  lu_decompose(a,n,indx,d)

  !     Find inverse by columns
  DO j=1,n
     CALL lu_linear_equation(a,n,indx,y(:,j))
  ENDDO
  !     The original matrix a was destroyed, now we equate it with the inverse y 
  a=y

  DEALLOCATE ( y ); DEALLOCATE ( indx )

END SUBROUTINE matinv

!     Given an NxN matrix A(N,N), this routine replaces it by the LU 
!     decomposed one, where the matrix elements are stored in the same 
!     matrix A. The array indx is  an output vector which records the row
!     permutation effected by the partial pivoting. d is the determinant
!
SUBROUTINE lu_decompose(a,n,indx,d)
  USE constants
  IMPLICIT NONE
  INTEGER :: n, i, j, k, imax
  REAL(DP) :: sum , tiny, aamax, dum, d
  REAL(DP), DIMENSION(n,n) :: a
  INTEGER, DIMENSION(n) :: indx
  REAL(DP), ALLOCATABLE :: vv(:)

  tiny=1.0e-20
  ALLOCATE ( vv(n) )
  D=1.
  DO i=1,n
     aamax=0.
     DO j=1,n
        IF (ABS(a(i,j)) > aamax) aamax=ABS(a(i,j))
     ENDDO
     !     Zero is the largest element
     IF (aamax == 0.) STOP 'Singular matrix.'
     !     No nonzero largest element
     vv(i)=1./aamax
  ENDDO
  !     loop over columns
  DO j=1,n
     !     solves equation 2.3.12 except for i=j of Numerical Recipes
     IF (j > 1) THEN
        DO i=1,j-1
           sum=a(i,j)
           IF (i > 1)THEN
              DO k=1,i-1
                 sum=sum-a(i,k)*a(k,j)
              ENDDO
              a(i,j)=sum
           ENDIF
        ENDDO
     ENDIF
     !    start searching for largest pivot element
     aamax=0.
     DO i=j,n
        sum=a(i,j)
        IF (j > 1)THEN
           DO k=1,j-1
              sum=sum-a(i,k)*a(k,j)
           ENDDO
           a(i,j)=sum
        ENDIF
        dum=vv(i)*ABS(sum)
        IF (dum >= aamax) THEN
           imax=i
           aamax=dum
        ENDIF
     ENDDO
     !    interchange of rows
     IF (j /= imax)THEN
        DO k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        ENDDO
        !    change of parity for determinant
        d=-d
        vv(imax)=vv(j)
     ENDIF
     indx(j)=imax
     IF(j /= n) THEN
        IF(a(j,j) == 0.) a(j,j)=tiny
        dum=1./a(j,j)
        DO i=j+1,n
           a(i,j)=a(i,j)*dum
        ENDDO
     ENDIF
     !    set up determinant
     d=d*a(j,j)
  ENDDO
  IF(a(n,n) == 0.)  a(n,n)=tiny
  DEALLOCATE ( vv)

END SUBROUTINE lu_decompose

!     Solves set of linear equations Ax=b, A is input as an LU decompomsed
!     matrix and indx keeps track of the permutations of the rows. b is input
!     as the right-hand side vector b and returns the solution x. A, n and indx
!     are not modified by this routine. This function takes into that b can contain
!     many zeros and is therefore suitable for matrix inversion


SUBROUTINE lu_linear_equation(a,n,indx,b)
  USE constants
  IMPLICIT NONE
  INTEGER :: n, ii, ll, i, j
  REAL(DP) :: sum 
  REAL(DP), DIMENSION(n,n) :: a
  REAL(DP), DIMENSION(n) :: b
  INTEGER, DIMENSION(n) :: indx

  ii=0
  !     First we solve equation 2.3.6 of numerical recipes 
  DO i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     IF (ii /= 0)THEN
        DO j=ii,i-1
           sum=sum-a(i,j)*b(j)
        ENDDO
     ELSEIF (sum /= 0.) THEN
        ii=i
     ENDIF
     b(i)=sum
  ENDDO
  !     then we solve equation 2.3.7
  DO i=n,1,-1
     sum=b(i)
     IF (i < n) THEN
        DO j=i+1,n
           sum=sum-a(i,j)*b(j)
        ENDDO
     ENDIF
     !     store a component of the solution x in the same place as b
     b(i)=sum/a(i,i)
  ENDDO

END SUBROUTINE lu_linear_equation
! 
! routine for calculating the principal square root ( and inverse ) 
! of a general matrix AA, i.e. solving the matrix equation: AA - A = 0. 
! The routine is built on the property : 
!
!            | 0  AA |      | 0    A | 
! sign(B) = (|       |) =  (|        |)
!            | 1   0 |      | A^-1 0 |
! 
! the sign of the matrix AA is calculated using Newtons iteration method,
! see Higham et. al. ref. Numerical Algorithms 15 (1997) 227-242 
! 
SUBROUTINE sqrtmat( aa, a, a_inv ,n ) 
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: aa
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT) :: a, a_inv
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: x, x0, x1, x2
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: i_mat, temp
  INTEGER :: i,j,k
  COMPLEX*16 :: d

  ALLOCATE( x(n+n,n+n), x0(n+n,n+n), x1(n+n,n+n),x2(n+n,n+n))
  ALLOCATE( i_mat(n,n),temp(n,n))
  ! setup real identity matrix only
  i_mat = (0.D0,0.D0)
  DO i = 1, n
     i_mat(i,i) = (1.d0, 0.d0)
  ENDDO
  DO i = 1, 2*n
     DO j = 1, 2*n
        x0(j,i) = (0.d0, 0.d0)
        x1(j,i) = (0.d0, 0.d0)
        x2(j,i) = (0.d0, 0.d0)
        x(j,i) = (0.d0, 0.d0)
     ENDDO
  ENDDO
  DO i = 1, n
     DO j = 1, n
        temp(j,i) = (0.d0, 0.d0)
     ENDDO
  ENDDO
  DO i = n+1, 2*n
     DO j = 1, n
        x0(j,i) = aa(j,i-n)
     ENDDO
  ENDDO
  DO i = 1, n
     DO j = n+1, 2*n
        x0(j,i) = i_mat(j-n,i)
     ENDDO
  ENDDO
  k = 0
  DO WHILE( MAXVAL(ABS(temp-aa)) > 1.D-14 .AND.  k < 1000 )
     x1 = x0 
     x2 = x0 
     CALL cmplxmatinv(x2,n+n,d)
     x = 0.5d0 * ( x1 + x2 ) 
     x0 = x 
     k = k + 1
     DO i = 1, n
        DO j =1, n
           a(i,j) = x(i,j+n)
           a_inv(i,j) = x(i+n,j) 
        ENDDO
     ENDDO
     temp = MATMUL( a,a )
  ENDDO
  DEALLOCATE(i_mat,temp); DEALLOCATE(x,x0,x1,x2)

END SUBROUTINE sqrtmat
!
!    F90 program library, adapted from Numerical Recipes
!    All functions have been translated to F90 from F77
!    This is the complex*16 version of the 
!    routines to do matrix inversion, from Numerical
!    Recipes, Teukolsky et al. Routines included
!    below are MATINV, LUDCMP and LUBKSB. See chap 2
!    of Numerical Recipes for further details
!
SUBROUTINE cmplxmatinv(a,n,d)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  INTEGER :: i, j
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT)  :: a
  COMPLEX*16, ALLOCATABLE :: y(:,:)
  COMPLEX*16 :: d
  INTEGER, ALLOCATABLE :: indx(:)

  ALLOCATE (y( n, n))  ; ALLOCATE ( indx (n))
  y=0.
  !     setup identity matrix
  DO i=1,n
     y(i,i)=(1.d0, 0.d0) 
  ENDDO
  !     LU decompose the matrix just once
  CALL  cmplxlu_decompose(a,n,indx,d)

  !     Find inverse by columns
  DO j=1,n
     CALL cmplxlu_linear_equation(a,n,indx,y(:,j))
  ENDDO
  !     The original matrix a was destroyed, now we equate it with the inverse y 
  a=y

  DEALLOCATE ( y ); DEALLOCATE ( indx )

END SUBROUTINE cmplxmatinv
!
!     Given an NxN matrix A(N,N), this routine replaces it by the LU 
!     decomposed one, where the matrix elements are stored in the same 
!     matrix A. The array indx is  an output vector which records the row
!     permutation effected by the partial pivoting. d is the determinant
!
SUBROUTINE cmplxlu_decompose(a,n,indx,d)
  USE constants
  IMPLICIT NONE
  INTEGER :: n, i, j, k, imax
  COMPLEX*16 :: sum, dum, tiny, d, aamax
  COMPLEX*16, DIMENSION(n,n) :: a
  INTEGER, DIMENSION(n) :: indx
  COMPLEX*16, ALLOCATABLE :: vv(:)

  tiny= ( 1.0D-20, 0.d0 ) 
  ALLOCATE ( vv(n) )
  D=1.
  DO i=1,n
     aamax=(0.D0,0.D0)
     DO j=1,n
        IF (ABS(a(i,j)) > ABS(aamax) ) aamax=ABS(a(i,j))
     ENDDO
     !     Zero is the largest element
     IF (aamax == (0.D0,0.D0)) STOP 'Singular matrix.'
     !     No nonzero largest element
     vv(i)=1./aamax
  ENDDO
  !     loop over columns
  DO j=1,n
     !     solves equation 2.3.12 except for i=j of Numerical Recipes
     IF (j > 1) THEN
        DO i=1,j-1
           sum=a(i,j)
           IF (i > 1)THEN
              DO k=1,i-1
                 sum=sum-a(i,k)*a(k,j)
              ENDDO
              a(i,j)=sum
           ENDIF
        ENDDO
     ENDIF
     !    start searching for largest pivot element
     aamax=(0.D0,0.D0)
     DO i=j,n
        sum=a(i,j)
        IF (j > 1)THEN
           DO k=1,j-1
              sum=sum-a(i,k)*a(k,j)
           ENDDO
           a(i,j)=sum
        ENDIF
        dum=vv(i)*ABS(sum)
        IF (ABS( dum ) >= ABS( aamax) ) THEN
           imax=i
           aamax=dum
        ENDIF
     ENDDO
     !    interchange of rows
     IF (j /= imax)THEN
        DO k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        ENDDO
        !    change of parity for determinant
        d=-d
        vv(imax)=vv(j)
     ENDIF
     indx(j)=imax
     IF(j /= n) THEN
        IF(a(j,j) == 0.) a(j,j)=tiny
        dum=1./a(j,j)
        DO i=j+1,n
           a(i,j)=a(i,j)*dum
        ENDDO
     ENDIF
     !    set up determinant
     d=d*a(j,j)
  ENDDO
  IF(a(n,n) == (0.d0,0.d0) )  a(n,n)=tiny
  DEALLOCATE ( vv)

END SUBROUTINE cmplxlu_decompose
!
!     Solves set of linear equations Ax=b, A is input as an LU decompomsed
!     matrix and indx keeps track of the permutations of the rows. b is input
!     as the right-hand side vector b and returns the solution x. A, n and indx
!     are not modified by this routine. This function takes into that b can contain
!     many zeros and is therefore suitable for matrix inversion
!
SUBROUTINE cmplxlu_linear_equation(a,n,indx,b)
  USE constants
  IMPLICIT NONE
  INTEGER :: n, ii, ll, i, j
  COMPLEX*16 :: sum 
  COMPLEX*16, DIMENSION(n,n) :: a
  COMPLEX*16, DIMENSION(n) :: b
  INTEGER, DIMENSION(n) :: indx

  ii=0
  !     First we solve equation 2.3.6 of numerical recipes 
  DO i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     IF (ii /= 0)THEN
        DO j=ii,i-1
           sum=sum-a(i,j)*b(j)
        ENDDO
     ELSEIF (sum /= 0.) THEN
        ii=i
     ENDIF
     b(i)=sum
  ENDDO
  !     then we solve equation 2.3.7
  DO i=n,1,-1
     sum=b(i)
     IF (i < n) THEN
        DO j=i+1,n
           sum=sum-a(i,j)*b(j)
        ENDDO
     ENDIF
     !     store a component of the solution x in the same place as b
     b(i)=sum/a(i,i)
  ENDDO

END SUBROUTINE cmplxlu_linear_equation
!
!     Does the lagrangian interpolation                    
!     see abramowitz & stegun, p.878 & p.882            
!     x=pt. at which you evaluate                       
!     w= interpolation points                           
!     f= values of f at w's                             
!
SUBROUTINE interpolate(x,w,f,val)
  USE constants
  IMPLICIT NONE   
  REAL(DP) :: x, w, f, val, pip, p_i
  INTEGER :: i, k
  DIMENSION w(n_startenergy_g), f(n_startenergy_g), pip(n_startenergy_g)

  pip=1.
  p_i=1.0d0
  val=0.0d0
  DO i=1,n_startenergy_g
     IF (x == w(i)) val=f(i)
  ENDDO
  IF(val == 0.) THEN
     DO i=1,n_startenergy_g
        p_i=p_i*(x-w(i))
     ENDDO
     DO k=1,n_startenergy_g
        pip(k)=PRODUCT(w(k)-w, mask=(w(k) /= w) )
     ENDDO
     DO i=1,n_startenergy_g
        val = val + p_i*f(i)/((x-w(i))*pip(i))
     ENDDO
  ENDIF

END SUBROUTINE interpolate

