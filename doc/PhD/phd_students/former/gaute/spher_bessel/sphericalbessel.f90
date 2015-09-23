PROGRAM bessel
IMPLICIT NONE 
INTEGER :: n 
DOUBLE PRECISION :: k,  x, y, sj, a, b

WRITE(*,*) 'read in order of sph. besselfunction'
READ(*,*) n 
!WRITE(*,*) 'read in real momentum k'
!READ(*,*) k
!WRITE(*,*) 'read in interval for sphbes to be evaluated'
!READ(*,*) a, b 

! create output file

!OPEN(6,FILE='output.dat')
!DO y=0.1,2,0.1
x=2

!call subroutine sphbes, with input n, x and output sj, sy, sjp and syp

CALL sphbes(n,x,sj)

! write output to the  output.dat file

WRITE(*,*) x, sj
 
!END DO
END PROGRAM bessel

