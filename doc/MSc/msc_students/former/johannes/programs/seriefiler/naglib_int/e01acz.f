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
