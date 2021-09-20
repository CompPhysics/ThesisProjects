program chipot_test
! NNLO_opt chiral interaction
use chiral_constants_nnlo_opt , only : init_chp_constants_nnlo_opt
use chiral_potentials_nnlo_opt

implicit none
integer,parameter :: wp=selected_real_kind(15)
! energy scale
real(wp) :: e_scale
! rho and npart are only needed to get the volume which goes in the definition of e_scale
integer:: npart
real(wp) :: rho
type(t_vector) :: stateP,stateQ,stateR,stateS
! matrix element
complex*16 :: matel
! these are needed only to set the energy scale
npart =4
rho=0.08
! initialize potential
call init_chp_constants_nnlo_opt(rho,npart,e_scale)
! assign spin/isospin to states
stateP%s=1;   stateP%t=1
stateQ%s=-1; stateQ%t=1
stateR%s=1;   stateR%t=1
stateS%s=-1;  stateS%t=1
! assign momenta to states
stateP%k(:)=0._wp
stateQ%k(:)=0._wp; stateQ%k(1)=1._wp
stateR%k(:)=0._wp; stateR%k(1)=1._wp
stateS%k(:)=0._wp
! get matrix element
! -NOTE: we need to check beforehand momentum and spin conservation
matel=chiral_pot_nnlo_opt_np(stateP,stateQ,stateR,stateS)
! restore energy units
matel=matel*e_scale

write(6,*) "matel: ",matel

end program
