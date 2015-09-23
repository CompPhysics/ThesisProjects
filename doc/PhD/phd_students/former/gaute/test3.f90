module arrays
  integer :: n_rel
  double precision, dimension(100) :: rgkk
  double precision, dimension(200) :: ra, w_ra
end module arrays

program interpolate
  use arrays
  implicit none
  double precision, dimension(100,100) :: g_bhf
  double precision :: bessi,x
  integer :: i, j, ncoup
!  external bessi
  
!  OPEN(9, FILE='eff_pot.dat', STATUS='old' )

  OPEN(7, FILE='out1.dat')
  open(6,FILE= 'out2.dat')
  
  ! read in original potential
!  read(9,*)  rgkk, g_bhf
!  do i = 1, 100
!     write(6,*) rgkk(i), g_bhf(i,i)   
!  end do
 
  
  ! number of datapoints
  n_rel = 100
  rgkk(1) = 0.0001
  do i = 2, n_rel
     rgkk(i) =  (i-1)*0.4d0
  end do
  
  call v_pot_yukawa(g_bhf)
  do i = 1, 100
     write(7,*) rgkk(i), g_bhf(i,i)   
  end do

  ! the interpolation points ra
  
  CALL GAULEG( 0.D0, 10.d0, ra, w_ra, n_rel )
  
  
  ncoup = 1
  call g_interpolate(ncoup,g_bhf)


  do i = 1, 100
     write(6,*) ra(i), g_bhf(i,i)   
  end do
  
end program interpolate



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!funksjonen interpolate kaller 2-dim lagrange interpol
!arrayen ra inneholder de nye punktene.
!i funksjonen lagrange-2dim inngår de gamle punktene 
!i vektoren rgkk. begge disse arrayene er i mitt tilfelle
!definert i modulen  wave_functions, så jeg definerer de ikke 
!i g_interpolate og lagrange-2dim
!n_rel er antall gitterpunkter.
!før kalling må du altså ha satt opp ra og rgkk sjølsagt!


SUBROUTINE g_interpolate(ncoup,g_bhf)
  use arrays
  IMPLICIT NONE
  INTEGER :: i, ncoup, i1,ev,cm, j, ndim, j1, lim1, lim2, lim3, lim4
  DOUBLE PRECISION, DIMENSION(ncoup*n_rel,ncoup*n_rel), INTENT(INOUT) ::  g_bhf
  DOUBLE PRECISION, ALLOCATABLE :: a(:,:), b(:,:)
  DOUBLE PRECISION :: mtx_element

  ndim=ncoup*n_rel

  ALLOCATE(a (ncoup*n_rel,ncoup*n_rel))
  ALLOCATE(b (ncoup*n_rel,ncoup*n_rel))
  a = 0.
  b = g_bhf
  DO i=1,ndim
     lim1 = 1; lim2 = n_rel 
     i1=i
     IF(i > n_rel) THEN 
       i1=i-n_rel
       lim1 = n_rel+1
       lim2 = n_rel+n_rel
     ENDIF
     DO j = 1, ndim
        lim3 = 1; lim4 = n_rel
        j1 = j
        IF(j > n_rel ) THEN 
           j1 = j -n_rel
           lim3 = n_rel+1
           lim4 = n_rel+n_rel
        ENDIF
        CALL lagrange_2dim(ra(j1),ra(i1),b(lim3:lim4,lim1:lim2),mtx_element)
        a(j, i) = mtx_element
     ENDDO
  ENDDO
  g_bhf = a

  DEALLOCATE(a, b )

END SUBROUTINE g_interpolate




!             *************************************************
!             *   Two-dimensional lagrange interpol           *
!             *************************************************

SUBROUTINE lagrange_2dim(x,y,a,fxy)
  USE arrays
      IMPLICIT NONE
      INTEGER :: nx0, ny0, i0, j0, k, l, m
      DOUBLE PRECISION, INTENT(INOUT) :: fxy
      DOUBLE PRECISION, INTENT(IN) :: a(n_rel,n_rel)
      DOUBLE PRECISION, INTENT(IN) :: x, y
      DOUBLE PRECISION, DIMENSION(n_rel) :: q
      DOUBLE PRECISION :: fxk, p
      fxy=0.d0
      nx0=n_rel
      ny0=n_rel
      q = rgkk
      IF(X > 0.5*(q(nx0-2)+q(nx0-1))) THEN
        i0=nx0-3
      ELSE 
        i0=-1
 400    i0=i0+1
        IF (x > 0.5*(q(i0+2)+q(i0+3))) GOTO 400
      ENDIF  
      IF(y > 0.5*(q(ny0-2)+q(ny0-1))) THEN
        j0=ny0-3
      ELSE 
        j0=-1
 410    j0=j0+1
        IF (y > 0.5*(q(j0+2)+q(j0+3))) GOTO 410
      ENDIF  
      DO k=1,3
         fxk=0.
         DO l=1,3
            p=1.d0
            DO m=1,3
               IF (m /= l) p=p*(x-q(i0+m))/(q(i0+l)-q(i0+m))
            ENDDO
            fxk=fxk+p*a(i0+l,j0+k)
         ENDDO
         p=1.d0
         DO l=1,3
            IF (k /= l) p=p*(y-q(j0+l))/(q(j0+k)-q(j0+l))
         ENDDO
         fxy=fxy+fxk*p
      ENDDO

      END SUBROUTINE lagrange_2dim


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
  SUBROUTINE gauleg(x1,x2,x,w,n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j, m
    DOUBLE PRECISION, INTENT(IN) :: x1, x2
    DOUBLE PRECISION, INTENT(INOUT) :: x, w
    DOUBLE PRECISION :: eps
    DIMENSION :: x(n), w(n)
    PARAMETER (eps=3.D-14)
    DOUBLE PRECISION :: p1,p2,p3,pp,xl,xm,z,z1

    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    DO i=1,m
       z1=0.
       z=COS(3.141592654d0*(i-.25d0)/(n+.5d0))
       DO WHILE ( ABS(z-z1) > EPS)
          p1=1.
          p2=0.
          DO j=1,n
             p3=p2
             p2=p1
             p1=((2.*j-1.)*z*p2-(j-1.)*p3)/j
          ENDDO
          pp=n*(z*p1-p2)/(z*z-1.)
          z1=z
          z=z-p1/pp
       ENDDO
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.*xl/((1.-z*z)*pp*pp)
       w(n+1-i)=w(i)
    ENDDO

  END SUBROUTINE gauleg


      SUBROUTINE v_pot_yukawa(vkk)
      USE arrays
      IMPLICIT NONE
      INTEGER i,j
      DOUBLE PRECISION  vkk, mu1, mu2, mu3, v_1, v_2, v_3, a, b, fac
      PARAMETER(mu1=0.49d0,v_1=-0.252d0)
      PARAMETER(mu2=7.84d0,v_2=-39.802d0)
      PARAMETER(mu3=24.01d0,v_3=156.359d0)
      DIMENSION vkk(n_rel,n_rel)

      DO i=1,n_rel     ! set up the free potential
         DO j=1,i
            a=(rgkk(j)+rgkk(i))**2
            b=(rgkk(j)-rgkk(i))**2
            fac=1./(4.d0*rgkk(i)*rgkk(j))
            vkk(j,i)=v_1*fac*DLOG((a+mu1)/(b+mu1))+ &
                    v_2*fac*DLOG((a+mu2)/(b+mu2))+&
                    v_3*fac*DLOG((a+mu3)/(b+mu3))
            vkk(i,j)=vkk(j,i)
         ENDDO
      ENDDO

      END  SUBROUTINE v_pot_yukawa






