module chiral_constants_nnlo_opt
  use mconstants
  
  integer,parameter,private :: wp=selected_real_kind(15)
  !nuclear constants
  real(wp), parameter :: proton_mass  = 938.272_wp   !in [MeV]
  real(wp), parameter :: neutron_mass = 939.5653_wp  !in [MeV]
  real(wp), parameter :: pi_m_mass    = 139.5702_wp  !in [MeV]
  real(wp), parameter :: pi_p_mass    = 139.5702_wp  !in [MeV]
  real(wp), parameter :: pi_0_mass    = 134.9766_wp  !in [MeV]

  ! chiral order definition
  INTEGER, PARAMETER :: LO   = 0
  ! all contributions vanish at order 1
  ! due to parity and time-reversal  invariance
  INTEGER, PARAMETER :: NLO  = 2
  INTEGER, PARAMETER :: NNLO = 3
  INTEGER, PARAMETER :: N3LO = 4
  
  real*8, parameter :: LEC_c1 =  -0.9186395287347203d0
  real*8, parameter :: LEC_c2 =   0.d0 
  real*8, parameter :: LEC_c3 =  -3.888687492763241d0
  real*8, parameter :: LEC_c4 =   4.310327160829740d0
  real*8, parameter :: LEC_c1_3NF =  -0.9186395287347203d0
  real*8, parameter :: LEC_c2_3NF =   0.d0 
  real*8, parameter :: LEC_c3_3NF =  -3.888687492763241d0
  real*8, parameter :: LEC_c4_3NF =   4.310327160829740d0
  
  
  real*8, parameter :: gA = 1.29d0 
  real*8, parameter :: gA2 = gA*gA 
  real*8, parameter :: gA4 = gA2*gA2 
  real*8, parameter :: fpi = 92.4d0 
  real*8, parameter :: fpi2 = fpi*fpi
  real*8, parameter :: fpi4 = fpi2*fpi2

  
  real*8 :: const(100) 
  real*8 :: const_3N(20)
  
  real*8 :: mnuc(-1:1), mpi(-1:2)
  REAL*8 :: mnuc2(-1:1)   ! nucleon mass squared
  REAL*8 :: mpi2(-1:2)    ! pion mass squared
  REAL*8 :: mpi4(-1:2)    ! mpi2 squared
  REAL*8 :: mpi5(-1:2)    ! mpi^5
  REAL*8 :: twompi(-1:2)  ! two times pion mass
  REAL*8 :: fourmpi2(-1:2)! four times pion mass squared
  COMPLEX*16 :: sigma_x(2,2), sigma_y(2,2), sigma_z(2,2)
  REAL*8 :: c1,c2,c3,c4
  REAL*8 :: C1_3NF, C2_3NF, C3_3NF, C4_3NF
  REAL*8, PARAMETER :: lambda = 500.D0
  REAL*8, PARAMETER :: lambda3N = 500.D0
  
  REAL*8 :: ilambda3,ilambda3NF_3
  REAL*8, PARAMETER :: sfr = 700.d0 
  REAL*8, PARAMETER :: sfr2 = sfr*sfr
  CHARACTER(LEN=2), PARAMETER :: chp_itope = 'EM'
  REAL*8 :: sfr_heavyside(-1:2) ! THETA(sfr-twompi)
  REAL*8 :: CS(-1:1), CT(-1:1), cnlo(1:7) 
  REAL*8 :: CE, CD, Econst, Dconst 
  
  !-- added by tn
  real*8 :: mom_scale,vol,a_lat,escale,iescale
  real*8 :: loop_s
  
CONTAINS 
      function delta (a,b) result(d)
        integer :: a,b
        real*8 ::d
        d=0.0d0
        if(a == b) d=1.0d0

    end function delta
    
  subroutine init_chp_constants_nnlo_opt(rho,npart,e_scale)
    use ang_mom_functions, only : commons_to_angmom
    implicit none   
    real*8,intent(in)  :: rho
    integer,intent(in) :: npart
    real*8,intent(out)  :: e_scale    
    
    double precision :: c1s0(-1:1), c3s1(-1:1), cnlo_pw(1:7)
    real*8 :: lambdachi
    integer :: i 
    
    vol = real(npart)/rho
    a_lat=vol**(1.d0/3.d0)
    mom_scale = (twopi/a_lat)
!    write(6,*) "kf = ",2.d0*mom_scale,neutron_mass,pi,a_lat
! #ifdef ISOSPIN_ON
    escale = hbarc*(hbarc/(0.5d0*(proton_mass+neutron_mass)))*mom_scale*mom_scale
! #else
!     escale = hbarc*(hbarc/neutron_mass)*mom_scale*mom_scale
! #endif

    e_scale = escale
    
    iescale=hbarc**3/vol/escale
    
    mnuc(-1) = proton_mass
    mnuc(0)  = 0.5d0*(proton_mass+neutron_mass)
    mnuc(+1) = neutron_mass
    
    mpi(-1)  = pi_m_mass
    mpi(0)   = pi_0_mass
    mpi(+1)  = pi_p_mass
    mpi(+2)  = ( mpi(-1)+mpi(0)+mpi(+1) ) /3.d0  
    
    
    mnuc2(:)   = mnuc*mnuc 
    mpi2(:)    = mpi*mpi   
    mpi4(:)    = mpi2*mpi2
    mpi5(:)    = mpi4*mpi 
    twompi(:)  = 2.0D0*mpi 
    fourmpi2(:)= 4.0D0*mpi2
    
    
    ilambda3= lambda**(-6.d0)
    
    ilambda3NF_3= lambda3N**(-6.d0)
    
    loop_s = chp_NLO_sfr_two_pion_exchange_loop_s(2)

    c1 = LEC_c1*1.0D-3
    c2 = LEC_c2*1.0D-3
    c3 = LEC_c3*1.0D-3
    c4 = LEC_c4*1.0D-3

    c1_3NF = LEC_c1_3NF*1.0D-3
    c2_3NF = LEC_c2_3NF*1.0D-3
    c3_3NF = LEC_c3_3NF*1.0D-3
    c4_3NF = LEC_c4_3NF*1.0D-3
  
    const = 0.0D0
    !1: 1/(2pi)^3
    const(1) = 1.0D0/(twopi**3)
    !2: gA^2/(4*fpi^2)
    const(2) = gA2/(4.0D0*fpi2)
    !3: 384pi^2*fpi^4
    const(3) = 384.0D0*pi2*fpi4
    !4: 5gA^4-4gA^2-1
    const(4) = 5.0D0*gA4-4.0D0*gA2-1.0D0
    !5: 23gA^4-10gA^2-1
    const(5) = 23.0D0*gA4-10.0D0*gA2 - 1.0D0
    !6: 48gA^4
    const(6) = 48.0D0*gA4
    !7: 3gA^4
    const(7) = 3.0D0*gA4
    !8: 64pi^2fpi^4
    const(8) = 64.0D0*pi2*fpi4
    !9: 3gA^2/16pifpi^4
    const(9) = 3.0D0*gA2/(16.0D0*pi*fpi4)
    !10: ga^2/16
    const(10) = ga2/16.0D0
    !11: 2.0D0*(2c1-c3)
    const(11) = 2.0D0*(2.0D0*c1-c3)
    !12 : const(7)/256pifpi^4
    const(12) = const(7)/(256.0D0*pi*fpi4)
    !13: gA^2/(128pifpi^4)
    const(13) = gA2/(128.0D0*pi*fpi4)
    !14: gA^4/(128pifpi^4)
    const(14) = gA4/(128.0D0*pi*fpi4)
    !15: 3gA^4/(512pifpi^4)
    const(15) = 3.0D0*gA4/(512.0D0*pi*fpi4)
    !16: gA2/(32pifpi^4)
    const(16) = gA2/(32.0D0*pi*fpi4)
    !17: gA2/8
    const(17) = gA2/8.0D0
    !18: gA4/(256pifpi^4)
    const(18) = gA4/(256.0D0*pi*fpi4)
    !19: 3gA4/(32pifpi^4)
    const(19) = 3.0D0*gA4/(32.0D0*pi*fpi4)
    !20: const(16)*(1-gA2)
    const(20) = const(16)*(1.0D0-gA2)
    !21: gA/(8*fpi^2)
    const(21) = gA/(8.0D0*fpi2)
    
    !1: gA^2/(8*fpi^2)
    const_3N(1) = gA2/(8.0D0*fpi2)
    !2: -4c1_3NF*mpi2/fpi^2
    const_3N(2) = -4.d0*c1_3NF*mpi2(2)/fpi2
    !3: 2c3_3NF/fpi^2
    const_3N(3) = 2.d0*c3_3NF/fpi2
    !4: c4/fpi2
    const_3N(4) = c4_3NF/fpi2
    
    
    sigma_x = 0.d0
    sigma_y = 0.d0
    sigma_z = 0.d0 
    
    sigma_x(1,2) = 1.d0 
    sigma_x(2,1) = 1.d0 
    
    sigma_z(1,1) = 1.d0 
    sigma_z(2,2) = -1.d0 
    
    sigma_y(1,2) = dcmplx(0.d0,-1.d0) 
    sigma_y(2,1) = dcmplx(0.d0, 1.d0) 
    
    DO i=-1,2
       IF ( (sfr-twompi(i)) <  0.0D0 ) sfr_heavyside(i) = 0.0D0
       IF ( (sfr-twompi(i)) >= 0.0D0 ) sfr_heavyside(i) = 1.0D0
    END DO
    

    !
    ! leading order contacts in PW
    !
    c1s0(-1) = -0.1513660372031080D+00
    c1s0(0)  = -0.1521410882366787D+00
    c1s0(1)  = -0.1517647459006913D+00
    
    c3s1(-1) = -0.1584341766228121D+00
    c3s1(0)  = -0.1584341766228121D+00
    c3s1(1)  = -0.1584341766228121D+00
    
    ! LO-CONTACTS IF THIS IS MAXIMUM ORDER
    ! THEN NO CIB
    !c1s0 =  -0.1521410882366787D+00
    !c3s1 =  -0.1584341766228121D+00

    ! the LO CIBcontacts are input in units of 10^4/GeV^2
    ! this = 10^-2/MeV^2 
    c1s0 = c1s0*0.01d0 
    c3s1 = c3s1*0.01d0 
    
    CS = (c1s0 + 3.d0*c3s1) /16.d0/pi
    CT = (c3s1 - c1s0) /16.d0/pi
    
!     write(6,*) "CS: ",CS
!     write(6,*) "CT: ",CT
    
    !
    ! next-to-leading order contacts in PW
    !
    cnlo_pw(1) =  0.2404021944134705D+01
    cnlo_pw(2) =  0.1263390763475578D+01
    cnlo_pw(3) =  0.4170455420556486D+00
    cnlo_pw(4) = -0.7826584999752046D+00
    cnlo_pw(5) =  0.9283846626623043D+00
    cnlo_pw(6) =  0.6181414190474576D+00
    cnlo_pw(7) = -0.6778085114063558D+00
    
    ! NLO
    ! the NLO contacts are input in units of 10^4/GeV^4
    ! this = 10^-8/MeV^4
    DO i=1,7
       CNLO_PW(i) = CNLO_PW(i) * 1.D-08
    END DO
    
    cnlo(1) = (-5.d0*cnlo_pw(7)+6.d0*cnlo_pw(5)-3.d0*cnlo_pw(4)-3.d0*cnlo_pw(3)-cnlo_pw(2)+2.d0*cnlo_pw(1))/(64.d0*pi)
    cnlo(2) = ( 5.d0*cnlo_pw(7)+6.d0*cnlo_pw(5)+3.d0*cnlo_pw(4)+3.d0*cnlo_pw(3)+cnlo_pw(2)+2.d0*cnlo_pw(1))/(16.d0*pi)
    cnlo(3) = -( 2.d0*cnlo_pw(7)-2.d0*dsqrt(2.d0)*cnlo_pw(6)-2.d0*cnlo_pw(5)-3.d0*cnlo_pw(3)+cnlo_pw(2)+2.d0*cnlo_pw(1))/(64.d0*pi)
    cnlo(4) = -(-2.d0*cnlo_pw(7)-2.d0*dsqrt(2.d0)*cnlo_pw(6)-2.d0*cnlo_pw(5)+3.d0*cnlo_pw(3)-cnlo_pw(2)+2.d0*cnlo_pw(1))/(16.d0*pi)
    cnlo(5) = -(-5.d0*cnlo_pw(7)+3.d0*cnlo_pw(4)+2.d0*cnlo_pw(2))/(16.d0*pi)
    cnlo(6) = ( cnlo_pw(7)-6.d0*dsqrt(2.d0)*cnlo_pw(6)-3.d0*cnlo_pw(4)+2.d0*cnlo_pw(2))/(64.d0*pi)
    cnlo(7) = -(cnlo_pw(7)+6.d0*dsqrt(2.d0)*cnlo_pw(6)-3.d0*cnlo_pw(4)+2.d0*cnlo_pw(2))/(16.d0*pi) 
    
    
    !WRITE(6,"(A12,F30.16)") 'C1', cnlo(1)* 1.D8
    !WRITE(6,"(A12,F30.16)") 'C2', cnlo(2)* 1.D8
    !WRITE(6,"(A12,F30.16)") 'C3', cnlo(3)* 1.D8
    !WRITE(6,"(A12,F30.16)") 'C4', cnlo(4)* 1.D8
    !WRITE(6,"(A12,F30.16)") 'C5', cnlo(5)* 1.D8
    !WRITE(6,"(A12,F30.16)") 'C6', cnlo(6)* 1.D8
    !WRITE(6,"(A12,F30.16)") 'C7', cnlo(7)* 1.D8
    
    
    ! NNLO 3NF constants
    !
    !
    !read(5,*);read(5,*) cE, cD 
    

    !cE = -0.398
    !cD = -0.39
    lambdaChi = 700d0 ! MeV
    Econst = cE/(fpi4*lambdaChi)
    Dconst = cD/(fpi2*lambdaChi)
    
    call commons_to_angmom()
    

  end subroutine init_chp_constants_nnlo_opt
 
  subroutine init_cEcD
    
    real*8 :: lambdachi
    lambdaChi = 700 ! MeV
    Econst = cE/(fpi4*lambdaChi)
    Dconst = cD/(fpi2*lambdaChi)


  end subroutine init_cEcD

  ! static one pion exchange, [1] eq 4.5
  ! without isospin structure
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_one_pion_exchange(q2, impi) RESULT(res)
  
    IMPLICIT NONE 
    REAL*8 , INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = -1.0D0 * const(2)/ (q2 + mpi2(impi))
    
        
  END FUNCTION chp_one_pion_exchange
  

  ! NLO function Eq 4.9
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_two_pion_exchange_Wc(q2, L, w, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: L
    REAL*8 , INTENT(IN) :: w
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = -L *(fourmpi2(impi)*const(4) + q2*const(5) + const(6)*mpi4(impi)/(w*w))/const(3) 
!    write(6,*) L,fourmpi2(impi),const(4)
    
  END FUNCTION chp_NLO_two_pion_exchange_Wc
  
  ! NLO function Eq 4.10
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_two_pion_exchange_Vs(q2,L,impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = +const(7)*L*q2/const(8)
    
  END FUNCTION chp_NLO_two_pion_exchange_Vs

  ! NLO Vt function Eq 4.10
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_two_pion_exchange_VT(L,impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0 
    res = -const(7)*L/const(8)
    
  END FUNCTION chp_NLO_two_pion_exchange_VT

  ! NLO loop function w [1] Eq 4.12 (DR and SFR)
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_two_pion_exchange_loop_w(q2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = DSQRT(fourmpi2(impi) + q2)
    
  END FUNCTION chp_NLO_two_pion_exchange_loop_w

  ! NLO SFR loop function s [2] Eq 2.16
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_sfr_two_pion_exchange_loop_s(impi) RESULT(res)
    
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = DSQRT(sfr2 - fourmpi2(impi))

  END FUNCTION chp_NLO_sfr_two_pion_exchange_loop_s
  
  ! NLO SFR loop function L [2] Eq 2.16
  ! q   : momentum transfer
  ! q2  : momentum transfer squared
  ! w,s : SFR loop functions
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_sfr_two_pion_exchange_loop_L(q, q2, w, s, impi) RESULT(res)
    
    REAL*8 , INTENT(INOUT) :: q
    REAL*8 , INTENT(INOUT) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: s
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res, eps 
    
    eps = 1.0D-8
    res = 0.0D0
    IF (sfr_heavyside(impi) == 0.0D0) return
    
    if ( q == 0.d0 ) then  
       q = q + eps 
       q2 = q*q 
!        res=s/sfr
!     else
      endif
       res = w * dlog( (sfr2*w*w+q2*s*s+2.0D0*sfr*q*w*s)/( fourmpi2(impi)*(sfr2+q2) ) )/(2.0D0*q)
!     endif
    !write(6,'(3g)') q, q2, res !(sfr-twompi(i)), sfr, sfr2, q,q2 ! (sfr2*w*w+q2*s*s+2.0D0*sfr*q*w*s), ( fourmpi2(impi)*(sfr2+q2) )
    !write(6,*) 0.5d0*( 2.d0*s/sfr/dsqrt(fourmpi2(impi)) - sfr2*fourmpi2(impi)) * w 
    
  END FUNCTION chp_NLO_sfr_two_pion_exchange_loop_L
  
   ! NNLO loop function wtilde SQUARED [1] Eq 4.20 (DR)
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NNLO_two_pion_exchange_loop_wtilde2(q2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = 2.0D0*mpi2(impi) + q2
    
  END FUNCTION chp_NNLO_two_pion_exchange_loop_wtilde2

  ! NNLO loop function wtilde [2] Eq 2.17 (SFR)
  ! q   : momentum transfer
  ! q2  : momentum transfer squared
  ! impi: determines which mpi to use, 
  FUNCTION chp_NNLO_sfr_two_pion_exchange_loop_A(q, q2, impi) RESULT(res)
    
    REAL*8 , INTENT(INOUT) :: q
    REAL*8 , INTENT(INOUT) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    REAL*8 :: eps
    
    res = 0.0D0
    eps = 1.0D-8
    IF (sfr_heavyside(impi) == 0.0D0) return
    
    if ( q == 0.d0 ) then  
       q = q + eps 
       q2 = q*q 
    end if
    
!    res = datan( q*(sfr-twompi(impi) )/(q2 + sfr*twompi(impi) ) )/(2.0D0*q)
    res = datan2( q*(sfr-twompi(impi) ),(q2 + sfr*twompi(impi) ) )/(2.0D0*q)
    
  END FUNCTION chp_NNLO_sfr_two_pion_exchange_loop_A

  ! NNLO function Eq [1] 4.13 and 4.21
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_Vc(q2, w, A, wt2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    REAL*8 , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = const(9)*(const(10)*mpi5(impi)/(mnuc(0)*w*w) - &
         (mpi2(impi)*const(11) - q2*c3 - q2*3.0D0*const(10)/mnuc(0) ) *wt2*A)
    
    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope == 'EM') THEN
       res = res - const(12)*(mpi(impi)*w*w+wt2*wt2*A)/mnuc(0)
    END IF

  END FUNCTION chp_NNLO_two_pion_exchange_Vc

  ! NNLO function Eq [1] 4.14 and 4.22
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_Wc(q2, w, A, wt2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    REAL*8 , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = const(13) * (3.0D0*gA2*mpi5(impi)/(w*w) - & 
         (fourmpi2(impi) + 2.0D0*q2 - gA2*(fourmpi2(impi)+3.0D0*q2))*wt2*A)/mnuc(0)
             
    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope == 'EM') THEN
       res = res + const(14)*(mpi(impi)*w*w + wt2*wt2*A)/mnuc(0)
    END IF
    
  END FUNCTION chp_NNLO_two_pion_exchange_Wc
  
  ! NNLO function Eq [1] 4.15 and 4.23
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_VT(w, A, wt2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    REAL*8 , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = 3.0D0*const(15)*wt2*A/mnuc(0)
             
    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope == 'EM') THEN
       res = res + const(15)*(mpi(impi) + w*w*A )/mnuc(0)
    END IF
    
  END FUNCTION chp_NNLO_two_pion_exchange_VT

  ! NNLO function Eq [1] 4.15 and 4.23
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_Vs(q2, w, A, wt2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    REAL*8 , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = -3.0D0*q2*const(15)*wt2*A/mnuc(0)

    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope == 'EM') THEN
       res = res - 1.0D0*const(15)*q2*(mpi(impi) + w*w*A)/mnuc(0)
    END IF
    
    !WRITE(*,*)q2
    !WRITE(*,*)const(15)
    !WRITE(*,*)wt2
    !WRITE(*,*)A
    !WRITE(*,*)mnuc(0)
    !WRITE(*,*)const(15)
    !WRITE(*,*)mpi(impi)
    !WRITE(*,*)w*w
    !WRITE(*,*)res
    
  END FUNCTION chp_NNLO_two_pion_exchange_Vs
  
  ! NNLO function Eq [1] 4.16 and 4.24
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_WT(q2, w, A, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = -1.0D0*const(16)*A*( (c4 + 1.0D0/(4.0D0*mnuc(0)))*w*w - &
         const(17)*(10.0D0*mpi2(impi) + 3.0D0*q2)/mnuc(0))
             
    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope == 'EM') THEN
       res = res - const(18)*(mpi(impi) + w*w*A)/mnuc(0)
    END IF
    
  END FUNCTION chp_NNLO_two_pion_exchange_WT 

  ! NNLO function Eq [1] 4.16 and 4.24
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_Ws(q2, w, A, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = q2*const(16)*A*( (c4+1.0D0/(4.0D0*mnuc(0)))*w*w - & 
         const(17)*(10.0D0*mpi2(impi) +3.0D0*q2)/mnuc(0))
             
    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope == 'EM') THEN
       res = res + const(18)*q2*(mpi(impi) + w*w*A)/mnuc(0)
    END IF
    
  END FUNCTION chp_NNLO_two_pion_exchange_Ws

  ! NNLO function Eq [1] 4.17
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  FUNCTION chp_NNLO_two_pion_exchange_VLS(A, wt2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: A
    REAL*8 , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = const(19) * wt2*A/mnuc(0)
             
  END FUNCTION chp_NNLO_two_pion_exchange_VLS

  ! NNLO function Eq [1] 4.18
  ! w   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! A     : NNLO loop function A Eq. 4.19
  ! impi  : determines which mpi2 to use
  FUNCTION chp_NNLO_two_pion_exchange_WLS(w, A, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = const(20)*w*w*A/mnuc(0)
             
  END FUNCTION chp_NNLO_two_pion_exchange_WLS

  ! regulator, eq 4.63 
  ! pfinal : final momentum
  ! pinit  : initial momentum
  ! n      : cutoff order 
  ! LAMBDA is accessed from the chp constant chp_lambda
  FUNCTION freg(pfinal, pinit, n) RESULT(res)
    
    REAL*8 , INTENT(IN) :: pfinal(3), pinit(3)
    REAL*8 , INTENT(IN) :: n
    REAL*8 :: res,exponent, p2, pp2 
    
    res = 0.0D0
    p2 = sum( pfinal*pfinal) 
    pp2 = sum( pinit*pinit )
    
    exponent = (p2**n/lambda**(2.d0*n) + &
         pp2**n/lambda**(2.0D0*n) )
    
    res = dexp(-exponent)
    
  END FUNCTION freg
  
  !in this case n=3d0
  FUNCTION fregnnlo(p2, pp2) RESULT(res)
    
    REAL*8 , INTENT(IN) :: p2, pp2 
    REAL*8 :: res,exponent
    
    exponent = (p2**3 + pp2**3)*ilambda3
    res = dexp(-exponent)
    
  END FUNCTION fregnnlo

  !in this case n=3d0
  FUNCTION freg3NF_local(p1, p2, p3) RESULT(res)
    
    REAL*8 , INTENT(IN) :: p2, p1,p3 
    REAL*8 :: res,exponent
    
    exponent = (p1**3 + p2**3+p3**3)*ilambda3NF_3
    res = dexp(-exponent)
    
  END FUNCTION freg3NF_local
  
  !in this case n=3d0
  FUNCTION freg3NF_nnlocal(p2, q2) RESULT(res)
    
    REAL*8 , INTENT(IN) :: p2, q2 
    REAL*8 :: res,exponent
    
    exponent = (((4.d0*p2 + 3.d0*q2)*0.25d0)**3)*ilambda3NF_3
    res = dexp(-exponent)
    
  END FUNCTION freg3NF_nnlocal
  
end module chiral_constants_nnlo_opt



module chiral_potentials_nnlo_opt
!     use mtypes,   only : wp
!     use non_int , only : sparray,bc_bvec
    use spin_isospin_operators

  implicit none 
  
  TYPE, PUBLIC :: chp_int_type
     INTEGER          :: val
     CHARACTER(LEN=12):: name
     LOGICAL          :: set
  END TYPE chp_int_type

  TYPE, PUBLIC :: chp_real_type
     REAL*8           :: val
     CHARACTER(LEN=12):: name
     LOGICAL          :: set
  END TYPE chp_real_type

  TYPE(chp_int_type), PRIVATE :: chp_chiral_order ! LO, NLO, NNLO, N3LO
  TYPE(chp_real_type), PRIVATE :: chp_regcut, chp_regcut_nnlo

! type for single-particle states
    type :: t_vector
       integer :: s,t
       real*8 :: k(3)
    end type
    
    public :: t_vector
  
!   integer,parameter,private :: wp=selected_real_kind(15)
  
contains 
    
  complex*16 function chiral_pot_nnlo_opt_np(stateP,stateQ,stateR,stateS)
    use chiral_constants_nnlo_opt
    use ang_mom_functions, only : tjs 
    
    implicit none 
    type(t_vector), intent(in) :: stateP,stateQ,stateR,stateS
!     INTEGER, INTENT(IN) :: p,q,r,s 
    INTEGER :: m1,m2,m3,m4, t1,t2,t3,t4, Tiso
    REAL*8 :: k1(3), k2(3), k3(3), k4(3), kmean(3) 
    REAL*8 :: qtrans(3), prel(3), pprel(3), qxk(3)
    REAL*8 :: q2, p2, kmean2, qabs, pp2, cg1, cg2 
    REAL*8 :: nucleon_mass, relativity_factor
    
    COMPLEX*16 :: vlo, vnlo, vnnlo, vn3lo, vdir
    COMPLEX*16 :: cont_lo, cont_nlo, cont_n3lo!, v1pe_tiso(0:1)
    ! LIST OF OPERATOR STRUCTURES
    COMPLEX*16 :: Vc
    COMPLEX*16 :: Wc
    COMPLEX*16 :: Vs
    COMPLEX*16 :: Ws
    COMPLEX*16 :: VLS
    COMPLEX*16 :: WLS
    COMPLEX*16 :: VsigL
    COMPLEX*16 :: WsigL
    COMPLEX*16 :: VT
    COMPLEX*16 :: WT
    COMPLEX*16 :: Vsigk
    COMPLEX*16 :: Wsigk
    COMPLEX*16 :: term_CS
    COMPLEX*16 :: term_CT
    COMPLEX*16 :: term_C(1:7)
    REAL*8 :: loop_w
    REAL*8 :: loop_wtilde2
    REAL*8 :: loop_L
    REAL*8 :: loop_A
    
    REAL*8 :: sds,dspins,dispins,isds
    COMPLEX*16 :: s1qds2q,sdqxk
    
    REAL*8 :: freg_nnlo
    LOGICAL :: lsds,ldspins,ldispins,lisds

    chp_chiral_order%val = NNLO 
    
    Vc      = 0.0D0
    Wc      = 0.0D0
    Vs      = 0.0D0
    Ws      = 0.0D0
    VLS     = 0.0D0
    WLS     = 0.0D0
    VsigL   = 0.0D0
    WsigL   = 0.0D0
    VT      = 0.0D0
    WT      = 0.0D0
    Vsigk   = 0.0D0
    Wsigk   = 0.0D0
  
    term_C = 0.0D0

    vlo = 0.0D0 ; vnlo = 0.0D0 ; vnnlo = 0.0D0 ; vn3lo = 0.0D0
    
    cont_lo = 0.0D0 ; cont_nlo = 0.0D0 ; cont_n3lo = 0.0D0

    sds=0.d0; isds=0.d0; sdqxk=0.d0; s1qds2q=0.d0
    lsds=.false.; lisds=.false.; ldspins=.false.; ldispins=.false.
    
    ! momenta in MeV
    k1(:) = mom_scale*hbarc*stateP%k(:)
    k2(:) = mom_scale*hbarc*stateQ%k(:)
    k3(:) = mom_scale*hbarc*stateR%k(:)
    k4(:) = mom_scale*hbarc*stateS%k(:)

    m1 = stateP%s 
    m2 = stateQ%s
    m3 = stateR%s
    m4 = stateS%s
  
    t1 = stateP%t
    t2 = stateQ%t
    t3 = stateR%t
    t4 = stateS%t
    ! 
    ! RELATIVE MOMENTA <prel |v| pprel > 
    ! 
    prel  = 0.5d0*(k1-k2)
    pprel = 0.5d0*(k3-k4)
    p2  = sum(prel*prel) 
    pp2 = sum(pprel*pprel)

    freg_nnlo=fregnnlo(p2, pp2)
!WARNING maybe dangerous... 1.d-10 should be safe enough(~4*lambda)
    if(freg_nnlo < 1.d-10) then
       vdir=0d0
       return
    endif
    
    !
    ! AVERAGE MOMENTA kav = 1/2(prel+pprel)
    !
    kmean = 0.5d0*( prel + pprel ) 
    kmean2 = sum(kmean*kmean)
    
    nucleon_mass = mnuc((t1+t2)/2)
    relativity_factor = nucleon_mass/ & 
         sqrt( sqrt(( nucleon_mass**2 + p2) * (nucleon_mass**2 + pp2 )) )
    
    !
    ! momentum transfer 
    !
    qtrans = prel - pprel
    q2 = sum(qtrans*qtrans) 
    qabs = dsqrt(q2)
    
    !
    !  cross product between momentum transfer and average momenta q X k 
    !
    qxk(1) = qtrans(2)*kmean(3)-qtrans(3)*kmean(2) 
    qxk(2) = qtrans(3)*kmean(1)-qtrans(1)*kmean(3) 
    qxk(3) = qtrans(1)*kmean(2)-qtrans(2)*kmean(1) 
    
    
    IF (chp_chiral_order%val == LO ) THEN
       
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !
       ! LEADING ORDER 
       ! 
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       
       chp_regcut%val = 2.0D0
       
       !
       ! leading order one-pion exchange 
       ! 
       WT = sigma1_dot_q_sigma2_dot_q(m1,m2,m3,m4,qtrans)*&
            chp_one_pion_exchange(q2, 2)*tau_dot_tau(t1,t2,t3,t4)
       
       vlo = WT
       
       
       !
       ! leading order contacts 
       !
       term_CS = CS(0)*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4)
       term_CT = CT(0)*sigma_dot_sigma(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) 
       cont_lo = term_CS + term_CT
       
    END IF
    
    IF (chp_chiral_order%val >= NLO) THEN
       
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !
       ! NEXT TO LEADING ORDER 
       ! 
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       
       chp_regcut%val = 2.0D0
       
       !
       ! NLO pion-exchange 
       ! 
       ! irreducible, or non-polynomial, NLO two-pion exchanges
       
       loop_w = chp_NLO_two_pion_exchange_loop_w(q2,2)
       loop_L = chp_NLO_sfr_two_pion_exchange_loop_L(qabs, q2, loop_w, loop_s, 2)
       
       lsds=.false.
       sds = sigma_dot_sigma(m1,m2,m3,m4)
       if(int(sds) /= 0) lsds=.true.
       lisds=.false.
       isds = tau_dot_tau(t1,t2,t3,t4)
       if(int(isds) /= 0) lisds=.true.
       ldspins=.false.
       dspins = delta(m1,m3)*delta(m2,m4)
       if(int(dspins)/=0) ldspins = .true.
       ldispins=.false.
       dispins = delta(t1,t3)*delta(t2,t4)
       if(int(dispins)/=0) ldispins = .true.
       s1qds2q = sigma1_dot_q_sigma2_dot_q(m1,m2,m3,m4,qtrans)
       sdqxk = spin_dot_qxk(m1,m2,m3,m4,qxk)
       
       ! LO CIB one-pion exchanges
       ! [1] Eq. 4.77-4.79
       
       !
       ! Note: the LO contribution is modified at NLO
       ! to include CIB effects. Therefore, should you go beyond
       ! LO, make sure to not include the vlo contribution 
       ! in the sum defining vdir.
       !
       
       
       ! 
       ! Below tjs is the 3-j symbol and is given in the module "ang_mom_functions" 
       ! You need to call "commons_to_angmom" to initialize the 3-j, 6-j symboles etc. 
       ! 
       WT = 0.d0 
       IF (t1 + t2 == 0) THEN 
          WT = 0.d0 
          DO Tiso = 0, 2, 2
             cg1 = tjs(1,1,Tiso,t1,t2,-(t1+t2))*(1.d0-2.d0*mod((t1+t2)/2,2))*sqrt(Tiso+1.d0)
             cg2 = tjs(1,1,Tiso,t3,t4,-(t3+t4))*(1.d0-2.d0*mod((t3+t4)/2,2))*sqrt(Tiso+1.d0)
             WT = WT + cg1*cg2*(-chp_one_pion_exchange(q2, 0) + (1.d0-2.d0*mod(Tiso/2+1,2))*2.0D0*chp_one_pion_exchange(q2,1)) 
          end DO
          

       END IF
       IF (t1 + t2 /= 0) THEN
          WT = chp_one_pion_exchange(q2, 0) 
       END IF

       vlo = WT*s1qds2q
      
       !
       ! leading order CIB contacts 
       !
       term_CS = 0.d0
       term_CT = 0.d0
       if(ldispins) then
           if(ldspins) term_CS = CS((t1+t2)/2)
           if(lsds) term_CT = CT((t1+t2)/2) * sds
       endif
       cont_lo = term_CS + term_CT
       
       Wc=0d0
       if(ldspins) Wc = isds*chp_NLO_two_pion_exchange_Wc(q2, loop_L, loop_w, 2)

       Vs=0d0
       VT=0d0
       if(ldispins) then
           if(lsds) Vs = chp_NLO_two_pion_exchange_Vs(q2, loop_L, 2)*sds 
           VT = chp_NLO_two_pion_exchange_VT(loop_L, 2)*s1qds2q
       endif
       
       vnlo = ( WT + WC + VS + VT ) !############ WARNING (check if WT should be here) #################
       
       !
       ! next-to-leading order contacts 
       !
       
       cont_nlo = 0.d0
       if(ldispins) then
         term_C(1) = cnlo(1)*q2*dspins
         term_C(2) = cnlo(2)*kmean2*dspins
         term_C(3) = cnlo(3)*q2*sds
         term_C(4) = cnlo(4)*kmean2*sds
         term_C(5) = cnlo(5)*sdqxk
         term_C(6) = cnlo(6)*s1qds2q
         term_C(7) = cnlo(7)*sigma1_dot_q_sigma2_dot_q(m1,m2,m3,m4,kmean)
         cont_nlo = SUM(term_C)
       endif      

    END IF
    
    IF (chp_chiral_order%val >= NNLO) THEN
       
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !
       ! NEXT TO NEXT TO LEADING ORDER 
       ! 
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       
       chp_regcut_nnlo%val = 3.0D0
       
       ! NNLO loop functions
       
       loop_wtilde2 = chp_NNLO_two_pion_exchange_loop_wtilde2(q2,2)
       loop_A = chp_NNLO_sfr_two_pion_exchange_loop_A(qabs, q2, 2)
       
       ! chp_itope type is set in the header of this function
       VLS = 0d0
       VT  = 0d0
       Vs = 0d0
       Vc = 0d0
       if(ldispins) then 
           if(ldspins) Vc  = chp_NNLO_two_pion_exchange_Vc (q2, loop_w, loop_A, loop_wtilde2, 2)
           VLS = chp_NNLO_two_pion_exchange_VLS(loop_A, loop_wtilde2, 2)*&
            sdqxk
           VT  = chp_NNLO_two_pion_exchange_VT (loop_w, loop_A, loop_wtilde2, 2)*&
             s1qds2q
           if(lsds) Vs  = chp_NNLO_two_pion_exchange_Vs (q2, loop_w, loop_A, loop_wtilde2, 2)*&
             sds
       endif
       
       Wc=0d0
       WLS=0d0
       WT=0d0
       Ws=0d0
       if(lisds) then
         if(ldspins) Wc  = chp_NNLO_two_pion_exchange_Wc (q2, loop_w, loop_A, loop_wtilde2, 2)*isds
         WLS = chp_NNLO_two_pion_exchange_WLS(loop_w, loop_A, 2)*sdqxk*isds
         WT  = chp_NNLO_two_pion_exchange_WT (q2, loop_w, loop_A, 2)*s1qds2q*isds
         if(lsds) Ws  = chp_NNLO_two_pion_exchange_Ws (q2, loop_w, loop_A, 2)*sds*isds
       endif
       vnnlo = (Vc + Wc + VLS + WLS + VT + WT + Vs + Ws)
       
    END IF
    
    IF (chp_chiral_order%val >= N3LO) THEN
       
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !
       ! NEXT TO NEXT TO NEXT TO LEADING ORDER 
       ! 
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

       chp_regcut_nnlo%val = 3.0D0
       
       vn3lo = 0.0D0
       cont_n3lo = 0.0D0

    END IF
    
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !
    ! DONE WITH CHIRAL ORDER CONTRIBUTIONS
    !
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    ! sum up all orders 
    ! regulator and relativity factor 
    IF (chp_chiral_order%val >= NNLO) THEN
       vdir = (vlo + cont_lo + vnlo + cont_nlo + vnnlo + cont_n3lo + vn3lo) * &
           relativity_factor * freg_nnlo
    ELSE
       vdir = (vlo + cont_lo + vnlo + cont_nlo + vnnlo + cont_n3lo + vn3lo) * &
             relativity_factor * freg(prel, pprel, chp_regcut%val)
    ENDIF

    chiral_pot_nnlo_opt_np =  vdir*iescale

  end function chiral_pot_nnlo_opt_np
  
end module chiral_potentials_nnlo_opt

