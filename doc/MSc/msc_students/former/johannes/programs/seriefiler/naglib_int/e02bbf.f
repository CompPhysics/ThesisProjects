      SUBROUTINE E02BBF(NCAP7, K, C, X,	S, IFAIL)
C     NAG LIBRARY SUBROUTINE  E02BBF
C
C     E02BBF  EVALUATES	A CUBIC	SPLINE FROM ITS
C     B-SPLINE REPRESENTATION.
C
C     DE BOOR*S	METHOD OF CONVEX COMBINATIONS.
C
C     USES NAG LIBRARY ROUTINE	P01AAF.
C
C     STARTED -	1973.
C     COMPLETED	- 1976.
C     AUTHOR - MGC AND JGH.
C
C     NAG COPYRIGHT 1975
C     MARK 5 RELEASE
C     MARK 7 REVISED IER-141 (DEC 1978)
      INTEGER NCAP7, IFAIL, P01AAF, IERROR, J, J1, L
C$P 1
      DOUBLE PRECISION SRNAME
      DOUBLE PRECISION K(NCAP7), C(NCAP7), X, S, C1, C2, C3, E2, E3, E4,
     * E5,K1, K2, K3, K4, K5, K6
      DATA SRNAME /8H E02BBF /
C
      IERROR = 0
      IF (NCAP7.GE.8) GO TO 10
      IERROR = 2
      GO TO 100
   10 IF (X.GE.K(4) .AND. X.LE.K(NCAP7-3)) GO TO 20
      IERROR = 1
      S	= 0.0D0
      GO TO 100
C
C     DETERMINE	 J  SUCH THAT  K(J + 3)	.LE. X .LE. K(J	+ 4).
C
   20 J1 = 0
      J	= NCAP7	- 7
   40 L	= (J1+J)/2
      IF (J-J1.LE.1) GO	TO 80
      IF (X.GE.K(L+4)) GO TO 60
      J	= L
      GO TO 40
   60 J1 = L
      GO TO 40
C
C     USE THE METHOD OF	CONVEX COMBINATIONS TO COMPUTE	S(X).
C
   80 K1 = K(J+1)
      K2 = K(J+2)
      K3 = K(J+3)
      K4 = K(J+4)
      K5 = K(J+5)
      K6 = K(J+6)
      E2 = X - K2
      E3 = X - K3
      E4 = K4 -	X
      E5 = K5 -	X
      C2 = C(J+1)
      C3 = C(J+2)
      C1 = ((X-K1)*C2+E4*C(J))/(K4-K1)
      C2 = (E2*C3+E5*C2)/(K5-K2)
      C3 = (E3*C(J+3)+(K6-X)*C3)/(K6-K3)
      C1 = (E2*C2+E4*C1)/(K4-K2)
      C2 = (E3*C3+E5*C2)/(K5-K3)
      S	= (E3*C2+E4*C1)/(K4-K3)
  100 IF (IERROR) 120, 140, 120
  120 IFAIL = P01AAF(IFAIL,IERROR,SRNAME)
      RETURN
  140 IFAIL = 0
      RETURN
      END
