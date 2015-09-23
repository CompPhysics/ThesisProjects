      SUBROUTINE E01ACF(A, B, X, Y, F, VAL, VALL, IFAIL, XX, WORK,W, D,
     * IG1, M1,	N1)
C     MARK 1 RELEASE.  NAG COPYRIGHT 1971
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
      INTEGER P01AAF, ISAVE, IFAIL, N1,	M1, M, N, J, I,	IG1
C$P 1
      DOUBLE PRECISION SRNAME
      DOUBLE PRECISION A, B, VAL, VALL,	E01ACZ,	X(N1), Y(M1), F(N1,M1),
     *XX(IG1), WORK(IG1), W(IG1), D(IG1)
      DATA SRNAME /8H E01ACF /
      ISAVE = IFAIL
      IF ((A.LT.X(1)) .OR. (A.GT.X(N1))) GO TO 100
      IF ((B.LT.Y(1)) .OR. (B.GT.Y(M1))) GO TO 100
      IFAIL = 0
      M	= M1 - 1
      N	= N1 - 1
      DO 40 J=1,M1
	 DO 20 I=1,N1
	    WORK(I) = F(I,J)
   20	 CONTINUE
	 CALL E01ACY(N1, X, WORK, W, D,	IG1)
	 XX(J) = E01ACZ(N1,X,WORK,W,A,IG1)
   40 CONTINUE
      CALL E01ACY(M1, Y, XX, W,	D, IG1)
      VAL = E01ACZ(M1,Y,XX,W,B,IG1)
      DO 80 I=1,N1
	 DO 60 J=1,M1
	    WORK(J) = F(I,J)
   60	 CONTINUE
	 CALL E01ACY(M1, Y, WORK, W, D,	IG1)
	 XX(I) = E01ACZ(M1,Y,WORK,W,B,IG1)
   80 CONTINUE
      CALL E01ACY(N1, X, XX, W,	D, IG1)
      VALL = E01ACZ(N1,X,XX,W,A,IG1)
      GO TO 120
  100 IFAIL = P01AAF(ISAVE,1,SRNAME)
  120 RETURN
      END


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


      DOUBLE PRECISION FUNCTION	E01ACZ(N, X, Y,	W, T, IG)
C     MARK 1 RELEASE.  NAG COPYRIGHT 1971
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
      INTEGER K, N2, N,	K1, IG
      DOUBLE PRECISION T, A, B,	C, V, X(N), Y(IG), W(IG)
      K	= 2
      N2 = N
   20 CONTINUE
      IF ((T.LE.X(K)) .OR. (K.EQ.N2)) GO TO 40
      K	= K + 1
      GO TO 20
   40 K1 = K - 1
      A	= X(K) - X(K1)
      B	= X(K) - T
      C	= T - X(K1)
      V	= (W(K1)*B*B*B+W(K)*C*C*C+(6.0D0*Y(K)-W(K)*A*A)*C+(6.0D0*Y(K1)-
     *W(K1)*A*A)*B)/(6.0D0*A)
      E01ACZ = V
      RETURN
      END
