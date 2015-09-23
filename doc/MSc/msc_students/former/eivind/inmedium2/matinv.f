
!
!            Routines to do mtx inversion, from Numerical
!            Recepies, Teukolsky et al. Routines included
!            below are MATINV, LUDCMP and LUBKSB. See chap 2
!            of Numerical Recipes for further details
!            Recoded in FORTRAN 90 by M. Hjorth-Jensen
!
SUBROUTINE matinv(a,n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  INTEGER :: i, j
  DOUBLE PRECISION, DIMENSION(n,n), INTENT(INOUT)  :: a
  DOUBLE PRECISION, ALLOCATABLE :: y(:,:)
  DOUBLE PRECISION :: d
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
  IMPLICIT NONE
  INTEGER :: n, i, j, k, imax
  DOUBLE PRECISION :: sum , tiny, aamax, dum, d
  DOUBLE PRECISION, DIMENSION(n,n) :: a
  INTEGER, DIMENSION(n) :: indx
  DOUBLE PRECISION, ALLOCATABLE :: vv(:)

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
  IMPLICIT NONE
  INTEGER :: n, ii, ll, i, j
  DOUBLE PRECISION :: sum 
  DOUBLE PRECISION, DIMENSION(n,n) :: a
  DOUBLE PRECISION, DIMENSION(n) :: b
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

!     determine eigenvalues and eigenvectors of a real symmetric
!     tri-diagonal matrix, or a real, symmetric matrix previously
!     reduced by function tred2 to tri-diagonal form. On input,
!     d[] contains the diagonal element and e[] the sub-diagonal
!     of the tri-diagonal matrix. On output d[] contains the
!     eigenvalues and  e[] is destroyed. If eigenvectors are
!     desired z[][] on input contains the identity matrix. If
!     eigenvectors of a matrix reduced by tred2() are required,
!     then z[][] on input is the matrix output from tred2().
!     On output, the k'th column returns the normalized eigenvector
!     corresponding to d[k]. 
!     The function is modified from the version in Numerical recipe.

SUBROUTINE tqli(d,e,n,z)
  IMPLICIT NONE
  INTEGER :: n 
  DOUBLE PRECISION  :: d(n),e(n),z(n,n)
  INTEGER :: i,iter,k,l,m
  DOUBLE PRECISION  :: b,c,dd,f,g,p,r,s,pythag,one

  DO i=2,n
     e(i-1)=e(i)
  ENDDO
  one=1.
  DO l=1,n
     iter=0
     ITERATE : DO
        DO m=l,n-1
           dd=ABS(d(m))+ABS(d(m+1))
           IF (ABS(e(m))+dd == dd) EXIT
        ENDDO
        IF(m == l) EXIT ITERATE
        IF(iter == 30) STOP 'too many iterations in tqli'
        iter=iter+1
        g=(d(l+1)-d(l))/(2.*e(l))
        r=pythag(g,one)
        g=d(m)-d(l)+e(l)/(g+sign(r,g))
        s=1.
        c=1.
        p=0.
        DO i=m-1,l,-1
           f=s*e(i)
           b=c*e(i)
           r=pythag(f,g)
           e(i+1)=r
           IF(r == 0.) THEN
              d(i+1)=d(i+1)-p
              e(m)=0.
              CYCLE ITERATE
           ENDIF
           s=f/r
           c=g/r
           g=d(i+1)-p
           r=(d(i)-g)*s+2.*c*b
           p=s*r
           d(i+1)=g+p
           g=c*r-b
           !     Omit lines from here ...
           DO k=1,n
              f=z(k,i+1)
              z(k,i+1)=s*z(k,i)+c*f
              z(k,i)=c*z(k,i)-s*f
           ENDDO
           !     ... to here when finding only eigenvalues.
        ENDDO
        d(l)=d(l)-p
        e(l)=g
        e(m)=0.
     ENDDO ITERATE
  ENDDO

END SUBROUTINE tqli


DOUBLE PRECISION FUNCTION pythag(a,b)
  DOUBLE PRECISION  :: a,b
  DOUBLE PRECISION  :: absa,absb
  absa=ABS(a)
  absb=ABS(b)
  IF(absa > absb) THEN
     pythag=absa*sqrt(1.+(absb/absa)**2)
  ELSE
     IF(absb == 0.) THEN
        pythag=0.
     ELSE
        pythag=absb*sqrt(1.+(absa/absb)**2)
     ENDIF
  ENDIF

END FUNCTION pythag

!    perform a Housholder reduction of a real symmetric matrix
!    a[][]. On output a[][] is replaced by the orthogonal matrix 
!    effecting the transformation. d[] returns the diagonal elements
!    of the tri-diagonal matrix, and e[] the off-diagonal elements, 
!    with e[0] = 0.
!    The function is modified from the version in Numerical recipes.

SUBROUTINE tred2(a,n,d,e,eigen_vectors)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: n
  DOUBLE PRECISION, INTENT(INOUT) :: a(n,n)
  DOUBLE PRECISION, INTENT(OUT) :: d(n),e(n)
  LOGICAL, INTENT(IN) :: eigen_vectors
  INTEGER :: i,j,k,l
  DOUBLE PRECISION ::  f,g,h,hh,scale

  DO i=n,2,-1
     l=i-1
     h=0.
     IF (l > 1) THEN
        scale=SUM(ABS(a(i,1:l)))
        IF (scale == 0.) THEN
           e(i)=a(i,l)
        ELSE
           a(i,1:l)=a(i,1:l)/scale
           h=sum(a(i,1:l)**2)
           f=a(i,l)
           g=-SIGN(SQRT(h),f)
           e(i)=scale*g
           h=h-f*g
           a(i,l)=f-g
           f=0.
           IF ( eigen_vectors) a(1:l,i)=a(i,1:l)/h
           DO j=1,l
              e(j)=(DOT_PRODUCT(a(j,1:j),a(i,1:j)) &
                   +DOT_PRODUCT(a(j+1:l,j),a(i,j+1:l)))/h 
           ENDDO
           f=DOT_PRODUCT(e(1:l),a(i,1:l))
           hh=f/(h+h)
           e(1:l)=e(1:l)-hh*a(i,1:l)
           DO j=1,l
              a(j,1:j)=a(j,1:j)-a(i,j)*e(1:j)-e(j)*a(i,1:j)
           ENDDO
        ENDIF
     ELSE
        e(i)=a(i,l)
     ENDIF
     d(i)=h
  ENDDO
  !     Omit following line if finding only eigenvalues.
  IF ( eigen_vectors) d(1)=0.
  e(1)=0.
  DO i=1,n
     IF ( eigen_vectors) THEN
        l=i-1
        IF (d(i) /= 0.) THEN
           DO j=1,l
              g=0.
              DO k=1,l
                 g=g+a(i,k)*a(k,j)
              ENDDO
              DO k=1,l
                 a(k,j)=a(k,j)-g*a(k,i)
              ENDDO
           ENDDO
        ENDIF
        a(i,i)=1.
        DO j=1,l
           a(i,j)=0.
           a(j,i)=0.
        ENDDO
     ELSE
        d(i)=a(i,i)
     ENDIF
  ENDDO

END SUBROUTINE tred2




