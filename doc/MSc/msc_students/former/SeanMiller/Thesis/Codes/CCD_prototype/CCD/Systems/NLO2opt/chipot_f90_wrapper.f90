subroutine chipot_f90_wrapper(matel_real,matel_im,npart,rho,ps,pt,ppx,ppy,ppz,qs,qt,qpx,qpy,qpz,rs,rt,rpx,rpy,rpz,ss,st,spx,spy,spz)
! NNLO_opt chiral interaction
use chiral_constants_nnlo_opt , only : init_chp_constants_nnlo_opt
use chiral_potentials_nnlo_opt

implicit none
integer, INTENT(IN) :: ps,pt,qs,qt,rs,rt,ss,st
integer, INTENT(IN) :: ppx,ppy,ppz
integer, INTENT(IN) :: qpx,qpy,qpz
integer, INTENT(IN) :: rpx,rpy,rpz
integer, INTENT(IN) :: spx,spy,spz

integer,parameter :: wp=selected_real_kind(15)
! energy scale
real(wp) :: e_scale
real(wp), INTENT(OUT) :: matel_real
real(wp), INTENT(OUT) :: matel_im
! rho and npart are only needed to get the volume which goes in the definition of e_scale
integer, INTENT(IN)  :: npart
real(wp), INTENT(IN) :: rho
type(t_vector) :: stateP,stateQ,stateR,stateS
! matrix element
complex*16 :: matel
! these are needed only to set the energy scale
!npart =4
!rho=0.08
!write(6,*) "rho: ",rho
!write(6,*) "npart: ",npart
! do I need to make momentum 0._wp? 15 zeros?
! right now they are integers.
! initialize potential
call init_chp_constants_nnlo_opt(rho,npart,e_scale)
! assign spin/isospin to states
stateP%s=ps;   stateP%t=pt
stateQ%s=qs;   stateQ%t=qt
stateR%s=rs;   stateR%t=rt
stateS%s=ss;   stateS%t=st
! assign momenta to states
!stateP%k(:)=0._wp
stateP%k(1) = ppx
stateP%k(2) = ppy
stateP%k(3) = ppz
!stateQ%k(:)=0._wp; stateQ%k(1)=1._wp
stateQ%k(1) = qpx
stateQ%k(2) = qpy
stateQ%k(3) = qpz
!stateR%k(:)=0._wp; stateR%k(1)=1._wp
stateR%k(1) = rpx
stateR%k(2) = rpy
stateR%k(3) = rpz
!stateS%k(:)=0._wp
stateS%k(1) = spx
stateS%k(2) = spy
stateS%k(3) = spz
! get matrix element
! -NOTE: we need to check beforehand momentum and spin conservation
matel=chiral_pot_nnlo_opt_np(stateP,stateQ,stateR,stateS)
!write(*,*) *, ppx, qpx, rpx, spx
!write(*,*) ppy, qpy, rpy, spy
!write(*,*) ppz, qpz, rpz, spz
!write(*,*) matel
! restore energy units
!matel=matel*e_scale
matel_real = RealPart(matel)
matel_im = ImagPart(matel)
!write(6,*) "matel: ",matel
!write(6,*) "matel_real: ",matel_real
!write(6,*) "matel_im: ",matel_im

end
