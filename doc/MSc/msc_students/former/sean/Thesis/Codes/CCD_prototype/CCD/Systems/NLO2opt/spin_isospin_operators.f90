!
! This module contains spin isospin operators that appear in usual nuclear potentials
!
module spin_isospin_operators

implicit none
public :: sigma_dot_sigma,tau_dot_tau,tau_dot_tau_vec_tau,sigma_dot_q
public :: sigma1_dot_q_sigma2_dot_q,sigma1_dot_q1_sigma2_dot_q2,spin_dot_qxk
private

  COMPLEX*16,parameter :: sigma_x(2,2)=reshape((/0.d0,1.d0,0.d0,1.d0/),(/2,2/))
  COMPLEX*16,parameter :: sigma_y(2,2)=reshape((/dcmplx(0.d0,0.d0),dcmplx(0.d0,-1.d0),dcmplx(0.d0,0.d0),dcmplx(0.d0,1.d0)/),(/2,2/))
  COMPLEX*16,parameter :: sigma_z(2,2)=reshape((/1.d0,0.d0,0.d0,-1.d0/),(/2,2/))

contains
  FUNCTION sigma_dot_sigma(ms1,ms2,ms3,ms4) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8 :: res
    
    res = 0.d0
    if((ms1 /= ms3).and.(ms2 /= ms4)) then
      res = dcmplx(1.d0,real(ms3))*dcmplx(1.d0,real(ms4))
    else if((ms1 == ms3).and.(ms2 == ms4)) then
      res = real(ms3*ms4)
    endif

  END FUNCTION sigma_dot_sigma

  FUNCTION tau_dot_tau(ms1,ms2,ms3,ms4) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8 :: res
    
    res = 0.0D0
    if((ms1 /= ms3).and.(ms2 /= ms4)) then
      res = dcmplx(1.d0,real(ms3))*dcmplx(1.d0,real(ms4))
    else if((ms1 == ms3).and.(ms2 == ms4)) then
      res = real(ms3*ms4)
    endif
    
  END FUNCTION tau_dot_tau

  !this is t_i .dot. (t_j x t_k) 
  FUNCTION tau_dot_tau_vec_tau(ms1,ms2,ms3,ms4,ms5,ms6) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4,ms5,ms6
    REAL*8 :: res
    complex*16 :: tjxtk(3)
    
    
    !first do (t_j x t_k)
    res = 0.0D0
    if(ms2==ms5) then
      if(ms3==ms6) then !in this case we have 0
       return
      else
        tjxtk(1)=dcmplx(0.d0,real(ms2)*real(ms3))
        tjxtk(2)=real(ms2)
        tjxtk(3)=0.d0
      endif
    else
      if(ms3==ms6) then !in this case we have 0
        tjxtk(1)=-dcmplx(0.d0,real(ms2)*real(ms3))
        tjxtk(2)=-real(ms3)
        tjxtk(3)=dcmplx(0.d0,real(ms2)*real(ms3))
      else
        tjxtk(1)=0.d0
        tjxtk(2)=0.d0
        tjxtk(3)=-dcmplx(0.d0,real(ms3))
      endif
    endif
    
    if(ms1==ms4) then !only z component
      res=real(ms1)*tjxtk(3)
    else !only x and y components
      res=tjxtk(1)-dcmplx(0.d0,real(ms1))*tjxtk(2)
    endif
    
  END FUNCTION tau_dot_tau_vec_tau
  
  FUNCTION sigma_dot_q(ms1,ms3,q) RESULT(res)

    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1, ms3
    REAL*8, INTENT(IN) :: q(3) 
    COMPLEX*16 :: res
    COMPLEX*16 :: res1
    COMPLEX*16 :: chi1(2), chi3(2), mat(2,2) 
    INTEGER :: i1 
    
!     chi1 = 0.d0; chi3 = 0.d0
!     i1 = nint(1.5-0.5*ms1) 
!     chi1(i1) = 1.d0 
!     i1 = nint(1.5-0.5*ms3) 
!     chi3(i1) = 1.d0 
    
!     mat = q(1)*sigma_x+q(2)*sigma_y+q(3)*sigma_z 
!     res1 = dot_product(chi1, matmul( mat, chi3))
!     res = res1 

    res=0.d0
    if(ms1 /= ms3) then
      res1 = dcmplx(q(1),q(2)*real(ms3))
    else if(ms1 == ms3) then
      res1 = real(ms3)*q(3)
    endif
    
  end FUNCTION sigma_dot_q

  FUNCTION sigma1_dot_q_sigma2_dot_q(ms1,ms2,ms3,ms4,q) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8, INTENT(IN) :: q(3) 
    COMPLEX*16 :: res
    COMPLEX*16 :: res1, res2 

    res1 = 0.d0
    if(ms1 /= ms3) then
      res1 = dcmplx(q(1),q(2)*real(ms3))
    else if(ms1 == ms3) then
      res1 = real(ms3)*q(3)
    endif
    
    res2 = 0.d0
    if(ms2 /= ms4) then
      res2 = dcmplx(q(1),q(2)*real(ms4))
    else if(ms2 == ms4) then
      res2 = real(ms4)*q(3)
    endif
    
    res = res1 * res2  

  END FUNCTION sigma1_dot_q_sigma2_dot_q

  FUNCTION sigma1_dot_q1_sigma2_dot_q2(ms1,ms2,q1,ms3,ms4,q2) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8, INTENT(IN) :: q1(3),q2(3) 
    COMPLEX*16 :: res
    COMPLEX*16 :: res1, res2 

    res1 = 0.d0
    if(ms1 /= ms3) then
      res1 = dcmplx(q1(1),q1(2)*real(ms3))
    else if(ms1 == ms3) then
      res1 = real(ms3)*q1(3)
    endif
    
    res2 = 0.d0
    if(ms2 /= ms4) then
      res2 = dcmplx(q2(1),q2(2)*real(ms4))
    else if(ms2 == ms4) then
      res2 = real(ms4)*q2(3)
    endif
    
    res = res1 * res2  

  END FUNCTION sigma1_dot_q1_sigma2_dot_q2  

  
    FUNCTION spin_dot_qxk(ms1,ms2,ms3,ms4,qxk) RESULT(res)

    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8, INTENT(IN) :: qxk(3)
    COMPLEX*16 :: res
    COMPLEX*16 :: res1, res2 
    
    res1 = 0.d0
    if(ms2 == ms4) then
      if(ms1 /= ms3) then
        res1 = dcmplx(qxk(1),qxk(2)*real(ms3))
      else if(ms1 == ms3) then
        res1 = real(ms3)*qxk(3)
      endif
    endif
    
    res2 = 0.d0
    if(ms1 == ms3) then
      if(ms2 /= ms4) then
        res2 = dcmplx(qxk(1),qxk(2)*real(ms4))
      else if(ms2 == ms4) then
        res2 = real(ms4)*qxk(3)
      endif
    endif
    
    res = -dcmplx(0.d0,0.5D0) *(res1+res2)
    
  END FUNCTION spin_dot_qxk
end module spin_isospin_operators
