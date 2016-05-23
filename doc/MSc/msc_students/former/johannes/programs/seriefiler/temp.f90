 !
  !     H.O. functions using Kummers function   
  !
  REAL(DP) FUNCTION rnl(n,l,z)
    IMPLICIT NONE
    INTEGER :: lll, nn
    INTEGER, INTENT(IN) :: l, n
    REAL(DP) :: y, dl, gamfaa, dfll, gamfab, dfnn
    REAL(DP), INTENT(IN) :: z

    rnl=0. ; y=0.5_dp*z*z
    IF(y > 60.0_dp) RETURN
    dl = l
    IF((ABS(z) < 1.0d-6) .AND. (l == 0)) rnl = 1.0_dp
    IF( ABS(z) > 1.0d-6) rnl = (z**l) * EXP(-y) * hypkum(n,dl+1.5_dp,z*z)
    gamfaa = 0.5_dp * SQRT(pi)
    IF(l /= 0) THEN
       DO lll = 1, l
          dfll = lll - 1
          gamfaa = gamfaa * (dfll + 1.5_dp)
       ENDDO
    ENDIF
    gamfab = gamfaa
    IF(n /= 0) THEN
       dfll = dl + 0.5_dp
       DO nn = 1, n
          dfnn = nn
          gamfab = gamfab * ((dfnn + dfll) / dfnn)
       ENDDO
    ENDIF
    rnl = rnl * (SQRT(2.0_dp * gamfab) / gamfaa)

  END FUNCTION rnl

   !
  !     Kummers function, Abramowitz & Stegun   
  !     exp. 13.1.2. a(there) equals (-n)       
  !  
  REAL(DP) FUNCTION hypkum(n,b,z)
    IMPLICIT NONE
    INTEGER :: nmax, nf
    INTEGER, INTENT(IN)  :: n
    REAL(DP) :: af, bf, zf, term, dfnf, xadd, sum
    REAL(DP), INTENT(IN) :: b, z

    IF(n < 0) WRITE (6,*)' error exit in hypkum ',  n,b,z
    hypkum = 1.0_dp
    IF(n == 0) RETURN
    nmax = n ; af = - n ; bf = b ; zf = z ; sum = 1.0 ; term = 1.0_dp
    DO nf = 1, nmax
       dfnf = nf
       xadd = dfnf - 1.0_dp
       term = term * ((af + xadd) / (bf + xadd)) * (zf / dfnf)
       IF(ABS(term) <  1.0d-12) EXIT
       sum = sum + term
    ENDDO
    hypkum = sum

  END FUNCTION hypkum
