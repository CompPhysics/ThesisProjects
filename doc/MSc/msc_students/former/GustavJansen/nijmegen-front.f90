MODULE myconst
    INTEGER , PUBLIC, PARAMETER :: nitot = 0
    INTEGER , PUBLIC, PARAMETER :: nbase = 2
    REAL(kind(1.0D0)), PUBLIC, PARAMETER :: pi = 3.141592741012573

    INTEGER , PUBLIC :: s, j, nqchn, nymod, ns, ncoup

END MODULE myconst

SUBROUTINE nijmegen_nn(q,nq,vkk,nv, trans, nt)
    USE myconst
    INTEGER :: i, i1, i2, k, k1, k2, j1
    INTEGER, INTENT(IN) :: nq, nv, nt
    REAL(KIND=8), INTENT(IN) :: q(nq)
    INTEGER, INTENT(IN) :: trans(nt)
    REAL(KIND=8), INTENT(OUT):: vkk(nv,nv)
    REAL(KIND=8) :: vnnpot(10,8)
    REAL(KIND=8) :: fac

    common/potsm0/vnnpot

    fac = (2.*pi)**3./(4.*pi)

    j1 = j+1
        
    do i = 1,nq
        i2 = i + nq
        do k = 1,nq
            CALL pwnnpot(nymod, j, q(i), q(k), nbase, nitot, nqchn, ns)
            k2 = k + nq
            IF (s == 0) THEN
                vkk(k,i) = vnnpot(j1,1)/fac
            ELSEIF (s == 1) THEN
                IF (j == 0) THEN
                    vkk(k2,i2) = vnnpot(j1,5)/fac
                ELSEIF(j > 0) THEN
                    IF (ncoup == 1) THEN
                        vkk(k,i) = vnnpot(j1,2)/fac
                    ELSE
                        vkk(k,i) = vnnpot(j1,3)/fac
                        vkk(k,i2) = vnnpot(j1,4)/fac
                        vkk(k2,i) = vnnpot(j1,6)/fac
                        vkk(k2,i2) = vnnpot(j1,5)/fac
                    ENDIF
                ENDIF
            ENDIF
        ENDDO
    ENDDO
END SUBROUTINE nijmegen_nn

SUBROUTINE nijmegen_yn(q,nq,vkk,nv,trans,nt)
    USE myconst
    INTEGER :: i, i1, i2, j1, k, k1, k2, a, b, s1, s2, nl
    INTEGER, INTENT(IN) :: nq, nv, nt
    REAL(KIND=8), INTENT(IN) :: q(nq)
    INTEGER, INTENT(IN) :: trans(nt)
    REAL(KIND=8), INTENT(OUT):: vkk(nv,nv)
    REAL(KIND=8) :: vynpot(10,8,3,3)

    common/potsm1/vynpot
    fac = (2.*pi)**3./(4.*pi)

    j1 = j+1
    nl = nq*nt
        
    do i = 1,nq
        do k = 1,nq
            CALL pwhnpot(nymod, j, q(i), q(k), nbase, nitot, nqchn, ns)
            do s1 = 1,nt
                a = trans(s1) +1
                do s2 = 1,nt
                    b = trans(s2) +1
                    i1 = i+nq*(s1-1)
                    i2 = i1+nl
                    k1 = k+nq*(s2-1)
                    k2 = k1+nl
                    IF (s == 0) THEN
                        vkk(k1,i1) = vynpot(j1,1,a,b)/fac
                    ELSEIF (s == 1) THEN
                        IF (j == 0) THEN
                            vkk(k2,i2) = vynpot(j1,5,a,b)/fac
                        ELSEIF (j > 0) THEN
                            IF (ncoup == 1) THEN
                                vkk(k1,i1) = vynpot(j1,2,a,b)/fac
                            ELSE
                                vkk(k1,i1) = vynpot(j1,3,a,b)/fac
                                vkk(k1,i2) = vynpot(j1,4,a,b)/fac
                                vkk(k2,i1) = vynpot(j1,6,a,b)/fac
                                vkk(k2,i2) = vynpot(j1,5,a,b)/fac
                            ENDIF
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ENDDO
END SUBROUTINE nijmegen_yn

SUBROUTINE nijmegen_yy(q,nq,vkk, nv, trans, nt)
    USE myconst
    INTEGER :: i, i1, i2, j1, k, k1, k2, a, b, s1, s2
    INTEGER, INTENT(IN) :: nq, nv, nt
    REAL(KIND=8), INTENT(IN) :: q(nq)
    INTEGER, INTENT(IN) :: trans(nt)
    REAL(KIND=8), INTENT(OUT):: vkk(nv,nv)
    REAL(KIND=8) :: vcnpot(10,8,6,6)

    common/potsm2/vcnpot
    fac = (2.*pi)**3./(4.*pi)

    j1 = j+1
    nl = nq*nt
        
    do i = 1,nq
        do k = 1,nq
            CALL pwxnpot(nymod, j, q(i), q(k), nbase, nitot, nqchn, ns)
            do s1 = 1,nt
                a = trans(s1)+1
                do s2 = 1,nt
                    b = trans(s2)+1
                    i1 = i+nq*(s1-1)
                    i2 = i1+nl
                    k1 = k+nq*(s2-1)
                    k2 = k1+nl
                    IF (s == 0) THEN
                        vkk(k1,i1) = vcnpot(j1,1,a,b)/fac
                    ELSEIF (s == 1) THEN
                        IF (j == 0) THEN
                            vkk(k2,i2) = vcnpot(j1,5,a,b)/fac
                        ELSEIF (j > 0) THEN
                            IF (ncoup == 1) THEN
                                vkk(k1,i1) = vcnpot(j1,2,a,b)/fac
                            ELSE
                                vkk(k1,i1) = vcnpot(j1,3,a,b)/fac
                                vkk(k1,i2) = vcnpot(j1,4,a,b)/fac
                                vkk(k2,i1) = vcnpot(j1,6,a,b)/fac
                                vkk(k2,i2) = vcnpot(j1,5,a,b)/fac
                            ENDIF
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ENDDO
        
END SUBROUTINE nijmegen_yy

