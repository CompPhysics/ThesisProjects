      SUBROUTINE E01ACY(NP1, X,	Y, W, D, IG)
C     MARK 1 RELEASE.  NAG COPYRIGHT 1971
C     MARK 3 REVISED.
C     MARK 4 REVISED. ER AM3-62
C     MARK 4.5 REVISED
      INTEGER NM, N, NP1, N1, I, I1, I2, IG
      DOUBLE PRECISION A, B, C,	E, DI1,	X(NP1),	Y(IG), W(IG), D(IG)
      N	= NP1 -	1
      NM = N - 1
      N1 = N + 1
      DO 20 I=1,NM
	 D(I) =	0.0D0
   20 CONTINUE
      I	= N
   40 CONTINUE
      I	= I - 1
      I1 = I + 1
      I2 = I + 2
      A	= X(I2)	- X(I1)
      B	= X(I1)	- X(I)
      C	= Y(I2)	- Y(I1)
      E	= Y(I1)	- Y(I)
      IF (I-N+1) 80, 60, 80
   60 D(I) = (X(I2)-X(I))/3.0D0
      W(I1) = C/A - E/B
      GO TO 100
   80 CONTINUE
      DI1 = D(I1)
      D(I) = (12.0D0*DI1*(X(I2)-X(I))-A*A)/(36.0D0*DI1)
      W(I1) = C/A - E/B	- A*W(I2)/(6.0D0*DI1)
  100 CONTINUE
      IF (I-1) 120, 120, 40
  120 CONTINUE
      W(1) = 0.0D0
      W(N1) = 0.0D0
      DO 180 I=1,NM
	 I1 = I	+ 1
	 IF (I-1) 160, 140, 160
  140	 CONTINUE
	 W(I1) = W(I1)/D(I)
	 GO TO 180
  160	 CONTINUE
	 W(I1) = (6.0D0*W(I1)-(X(I1)-X(I))*W(I))/(6.0D0*D(I))
  180 CONTINUE
      RETURN
      END
