C	MACMILLAN'S VARIATIONAL MONTECARLO IN LIQUID HELIUM-4 WITH
C	LENNARD-JONES POTENTIAL
	DIMENSION X(3,1000),XOLD(3),XNEW(3)
	REAL L,LSQU,LSQDIV4
	COMMON X,NP,L,LSQU,HL,LSQDIV4 !Arguments for subroutines
	WRITE(*,*) ' INTEGER FOR RANDOMNESS, NUMBER PARTICLES AND CORRELATION'
	READ(*,*)NR,NP,C	!C=B/SIGMA, Dimensionless
	DENS=0.3648*0.9		!SCHIFF-VERLET in SIGMA units
	L=(NP/DENS)**0.3333333	!Simulation cube
	LSQ=L*L
	LSQDIV4=LSQU/4.
	HL=L/2.
	SCAL=L/4.		!Scale of Metropolis moves
	C5=C*C*C*C*C
	WRITE(*,*)' ENTER NUMBER BLOCKS and NUMBER OF MOVES'
	READ(*,*)NBLK,NMOV
	DO 10 I=1,NP		!Starting postition
	   DO 10 J=1,3		!At random in simulation box
	      X(J,I)=L*RAN1(NR)
 10	   CONTINUE
C	METROPOLIS MOVES
	   ETOT=0		!Block average
	   ETOT2=0
	   DO 100 IBLK=0,NBLK
	      EAV=0
	      DO 200 IMOV=1,NMOV
		 DO 300 IPAR=1,NP !Move particles sequentially
		    DO 1000 I=1,3
		       XOLD(I)=X(I,IPAR) !Old position
		       XNEW(I)=XOLD(I)+SCAL*(RAN1(NR)-0.5) !Attempted move
		       IF(XNEW(I).LT.0.) XNEW(I)=XNEW(I)+L !Move particle to the
		       IF(XNEW(I).GT.0.) XNEW(I)=XNEW(I)-L !Fundamental cell
 1000		    CONTINUE
		    SUM=0.	!Old and New W.F.
		    DO 320 JPAR=1,NP
		       IF (JPAR.EQ.IPAR) GOTO 320 !Skip particle IPAR
		       CALL CLOSE (XOLD,JPAR,R2) !Old neighbour: R2=squared
		       IF (R2.GT.LSQDIV4) R2=LSQDIV4 !Cutoff at L/2
		       CALL CLOSE (XNEW,JPAR,S2) !New neighbour: S2=squared
		       IF (S2.GT.LSQDIV4) S2=LSQDIV4 !and cutoff
		       SUM=SUM+1./(R2*R2*SQRT(R2))-1./(S2*S2*SQRT(S2))
 320		    CONTINUE
C	QUOTIENT OF PROB: (W.F.NEW/W.F.OLD)**2
		    EXPR=C5*SUM	!Exponent
		    IF(EXPR.GT.0.) GOTO 330 !Accepted
		    PROBQ=EXP(C5*SUM) !Probability quotient
		    IF(PROBQ.LT.RAN1(NR)) GOTO 300 !Test Acceptance-Rejection
 330		    DO 1010 I=1,3 !Accepted move
 1010		       X(I,IPAR)=XNEW(I) !Update position
 300		    CONTINUE
		    CALL ENER(E,C) !Sample energy
		    EAV=EAV+E
 200		 CONTINUE
		 EAV=EAV/NMOV
		 WRITE(*,*)IBLK,EAV !Print iterations
		 IF (IBLK.EQ.0) GOTO 100 !Neglect first block
		 ETOT=ETOT+EAV
		 ETOT2=ETOT2+EAV**2
 100	      CONTINUE
	      ETOT=ETOT/NBLK/NP	!Block average
	      ETOT2=ETOT2/NBLK/NP**2
	      ERR=SQRT((ETOT2-ETOT**2)/(NBLK-1))
	      PRINT*,' E =',ETOT, ' ERR =',ERR
	      END
C	=================================================================
	SUBROUTINE CLOSE (XX,JPAR,R2)
C	DETERMINES THE IMAGE OF PARTICLE JPAR CLOSEST TO IPAR
C	THE RESULT R2 IS THE SQUARE OF THE DISTANCE
	DIMENSION X(3,1000),XX(3),XCLOSE(3)
	REAL L,LSQU,LSQDIV4
	COMMON X,NP,L,LSQU,HL,LSQDIV4
	COMMON/COORD/XCLOSE	!Components of (Ri-Rj)
	R2=0.
	DO 1000 I=1,3
	   XCLOSE(I)=ABS(XX(I)-X(I,JPAR)) !Always positive
	   IF(XCLOSE(I).GT.HL) XCLOSE(I)=L-XCLOSE(I)
 1000	   R2=R2+XCLOSE(I)**2	!Squared distance
	   RETURN
	   END
C	=================================================================
	SUBROUTINE ENER(E,C)
	DIMENSION X(3,1000),XCLOSE(3)
	REAL L,LSQU,LSQDIV4
	COMMON X,NP,L,LSQU,HL,LSQDIV4
	COMMON/COORD/XCLOSE
	!DATA A,B/9.27647,40.88	!A is 5*(HBAR**2/M)/SIGMA**2
	A=9.27647
	B=40.88
	E=0			!B is 4*EPSILON
	CKIN=A*C*C*C*C*C
	DO 10 IPAR=1,NP-1	!Count only PAIRS
	   DO 10 JPAR=IPAR+1,NP
	      CALL CLOSE (X(1,IPAR),JPAR,2.)
	      IF(R2.GT.LSQDIV4) GOTO 100 !Cutoff correlation
	      E=E+CKIN/(R2*R2*R2*SQRT(r2)) !Kinetic energy, only
C                                          !nearest neighbour
 100	      DO 20 K1=-1,0	!Check of images, X axis
		 Z1=(XCLOSE(1)+K1*L)**2
		 DO 20 K2=-1,0
		    Z2=(XCLOSE(2)+K2*L)**2 ! - Y axis and
		    DO 20 K3=-1,0
		       Z3=(XCLOSE(3)+K3*L)**2 ! - Z axis
		       R2=Z1+Z2+Z3
		       R6=R2**3	!This is R**6
		       IF (R2.GT.LSQU) GOTO 20 !Cutoff of potential at L
		       E=E+B*(1./R6-1.)/R6
 20		    CONTINUE
 10		 CONTINUE
		 RETURN
		 END
	FUNCTION RAN1(idum)
	INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
	REAL ran1,AM,EPS,RNMX
	PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32)
	PARAMETER (NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
	INTEGER j,k,iv(NTAB),iy
	SAVE iv,iy
	DATA iv /NTAB*0/, iy /0/
	if (idum.le.0.or.iy.eq.0) then
	   idum=max(-idum,1)
	   do 11 j=NTAB+8,1,-1
	      k=idum/IQ
	      idum=IA*(idum-k*IQ)-IR*k
	      if (idum.lt.0) idum=idum+IM
	      if (j.le.NTAB) iv(j)=idum
 11	   continue
	   iy=iv(1)
	endif
	k=idum/IQ
	idum=IA*(idum-k*IQ)-IR*k
	if (idum.lt.0) idum=idum+IM
	j=1+iy/NDIV
	iy=iv(j)
	iv(j)=idum
	ran1=min(AM*iy,RNMX)
	return
	END
	
