SUBROUTINE GENERAL_MESH( nk1, nk2, angle, b_trans, contour )
  USE CONSTANTS
  USE SETUP_CMPLX_MESH
  use wave_functions
  
  implicit none
  integer, intent(in) :: nk1, nk2
  double precision, intent(in) :: angle, b_trans
  character(len=100), intent(in) :: contour
  double precision :: klim_max, const, c_trans, xlim1, xlim2,xlim3, xlim4, step
  double precision, dimension(nk1) :: k1, wk1, k3, wk3, x1, wx1
  double precision, dimension(nk2) :: k2, wk2, x2, wx2
  integer :: i, nk
  

  ! total number of meshpoints
  if ( contour /= 'triangle' ) nk = nk1 + nk2
  if ( contour == 'triangle' ) nk = 2*nk1 + nk2
  
  !if ( nk /= n_lab ) write(6,*) 'wrong number of meshpoints in lab system'
  

  if ( contour == 'rotation_translation' ) then
     
     ! translation into the lower half k-plane by -c
     c_trans = sin(angle)*b_trans
     ! weighting of meshpoints
     !const = 1.5d0
     !const = 5.d0
     const = 2.d0
     
     klim_max = b_trans*cos(angle)
         
     ! x1 = (4./PI)*atan(klim_max/const) - 1.d0
     ! x2 = 1.d0
    
     xlim1 = klim_max 
     xlim2 = 4.d0 
     
     !x1 = -1.d0
     !x2 = 1.d0
  
     
     ! setup meshpoints along the line L1
     CALL gauss_legendre( 0.d0, b_trans, k1, wk1, nk1 )
     ! mesh for L2
     CALL gauss_legendre( xlim1, xlim2, x2, wx2, nk2 )
     
     do i= 1, nk2
        
        !k2(i) =   const * TAN( ( PI*0.25d0 ) * (1.d0 + x(i) ) )
        !wk2(i) =  const * (PI*0.25d0)*wx(i)* &
        !     ( 1.d0 + TAN( ( PI * 0.25d0 ) * ( 1.d0 + x(i) ) )**2 )
        k2(i) = x2(i)
        wk2(i) = wx2(i)
        
     end do
     
     !setup complex mesh along line L1 and line L2, and momenta for the t-matrix
     !phase = exp(-dcmplx(0.,1.)*angle)
     do i = 1, nk
        ! complex mesh along L1 & L2
        if (i <= nk1 ) then
           
           k_cmplx(i) = k1(i) * phase
           wk_cmplx(i) = wk1(i) * phase
        elseif ( i > nk1 .and. i <= nk2+nk1 ) then
           
           k_cmplx(i) =  dcmplx(k2(i-nk1),-c_trans)
           wk_cmplx(i) = dcmplx(wk2(i-nk1),0.)
        end if
        
     end do

  elseif ( contour == 'rotation' ) then


     const = 2.d0
     xlim1 = 0.d0 !-1.d0
     xlim2 = 2.d0 ! 1.d0
     xlim3 = 2.d0
     xlim4 = 10.d0
     CALL gauss_legendre( xlim1, xlim2, x1, wx1, nk1 )
     CALL gauss_legendre( xlim3, xlim4, x2, wx2, nk2 )
     
     do i= 1, nk1
        
     !   k2(i) =   const * TAN( ( PI*0.25d0 ) * (1.d0 + x(i) ) )
     !   wk2(i) =  const * (PI*0.25d0)*wx(i)* &
     !        ( 1.d0 + TAN( ( PI * 0.25d0 ) * ( 1.d0 + x(i) ) )**2 )
        k1(i) = x1(i)
        wk1(i) = wx1(i)
        
     end do
     
     do i= 1, nk2
        
        !   k2(i) =   const * TAN( ( PI*0.25d0 ) * (1.d0 + x(i) ) )
        !   wk2(i) =  const * (PI*0.25d0)*wx(i)* &
        !        ( 1.d0 + TAN( ( PI * 0.25d0 ) * ( 1.d0 + x(i) ) )**2 )
        k2(i) = x2(i)
        wk2(i) = wx2(i)
     end do
     
     
     do i = 1, nk
        ! complex mesh along L1 & L2
        if (i <= nk1 ) then
           
           k_cmplx(i) = k1(i) * phase
           wk_cmplx(i) = wk1(i) * phase
        elseif ( i > nk1 .and. i <= nk2+nk1 ) then
           
           k_cmplx(i) =  k2(i-nk1)*phase
           wk_cmplx(i) = wk2(i-nk1)*phase
        end if
        !write(6,*) k_cmplx(i)
     end do

  elseif ( contour == 'triangle' ) then
     
     ! translation into the lower half k-plane by -c
     c_trans = sin(angle)*b_trans
     ! weighting of meshpoints
     !const = 1.5d0
     !const = 5.d0
     const = 2.d0
     
     klim_max = b_trans*cos(angle)
    
     ! x1 = (4./PI)*atan(klim_max/const) - 1.d0
     ! x2 = 1.d0
     
     xlim1 = 2.d0*klim_max 
     xlim2 = 4.d0 
     
     !x1 = -1.d0
     !x2 = 1.d0
     !k_limit = xlim1*hbarc

     ! setup meshpoints along the line L1 and L2
     CALL gauss_legendre( 0.d0, b_trans, k1, wk1, nk1 )
     ! mesh for L3
     CALL gauss_legendre( xlim1, xlim2, x2, wx2, nk2 )
     
     do i= 1, nk2
        
        !k2(i) =   const * TAN( ( PI*0.25d0 ) * (1.d0 + x(i) ) )
        !wk2(i) =  const * (PI*0.25d0)*wx(i)* &
        !     ( 1.d0 + TAN( ( PI * 0.25d0 ) * ( 1.d0 + x(i) ) )**2 )
        k2(i) = x2(i)
        wk2(i) = wx2(i)
        
     end do
     
     !setup complex mesh along line L1 and line L2, and momenta for the t-matrix
     phase = exp(-dcmplx(0.,1.)*angle)
     do i = 1, nk
        ! complex mesh along L1 & L2
        if (i <= nk1 ) then
           
           k_cmplx(i) = k1(i) * phase
           wk_cmplx(i) = wk1(i) * phase
        elseif ( i > nk1 .and. i <= 2*nk1 ) then
           k_cmplx(i) = ( k1(i-nk1)) * conjg( phase ) + dcmplx(klim_max,-c_trans )
           wk_cmplx(i) = wk1(i-nk1) * conjg( phase )
        elseif ( i > 2*nk1 .and. i <= nk2+2*nk1 ) then
           
           k_cmplx(i) =  dcmplx(k2(i-2*nk1),0.d0)
           wk_cmplx(i) = dcmplx(wk2(i-2*nk1),0.d0)
        end if
                
     end do
     k_cmplx = k_cmplx*hbarc
     wk_cmplx = wk_cmplx*hbarc
     
  end if


end SUBROUTINE GENERAL_MESH
