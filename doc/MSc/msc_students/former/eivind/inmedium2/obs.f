C        READPH                                                         00000010
c************** this is the correct version ************
C                                                                       00000020
C     ROUTINE READS PHASESHIFTS FROM DATASET AND CALCULATES             00000030
C            SCATTERING OBSERVABLES                                     00000040
c            the phase shifts have to be in the Stapp convention
C                                                                       00000050
C     ROUTINE NEEDS THE FOLLOWING INPUT PARAMETERS:                     00000060
C       1) NUMBER OF J  (AT MOST 90)                                    00000070
C       2) MASSES OF PARTICLES                                          00000080
C       3) NUMBER OF ELABS  (AT MOST 50)                                00000090
C       4) IRAD = 0       : PHASE-SHIFTS ON DATASET ARE IN DEGREES      00000100
C          IRAD NON-EQ. 0 : PHASE-SHIFTS ON DATASET ARE IN RADIANS      00000110
C       5) ICSB = 0       : NEGLECT CHARGE SYMMETRY BREAKING            00000120
C          ICSB NON-EQ. 0 : INCLUDE CHARGE SYMMETRY BREAKING            00000130
C       6) ISYM = 0       : ALL PARTIAL WAVES                           00000140
C          ISYM NON-EQ. 0 : ONLY ANTISYMMETRIC LSJ-STATES               00000150
C       7) INELA = 0       : ELASTIC PHASE-SHIFTS                       00000160
C          INELA NON-EQ. 0 : READ INELASTIC PHASE-SHIFTS                00000170
C                                                                       00000180
C********  ROUTINE CAN NOT READ THE MIXING PARAMETER GAMMO              00000190
C          IN THE PRESENT STATE, BECAUSE THE DATA-FILE                  00000200
C          NEEDS MORE THAN LRECL =80                                    00000210
C          ONE HAS TO CHANGE THE WRITE STATEMENT IN PROGRAM             00000220
C********  PHASES                                                       00000230
C                                                                       00000240
C                                                                       00000250
      IMPLICIT REAL*8 (A-H,O-Z)                                         00000260
C                                                                       00000270
      COMMON /CRDWRT/ KREAD,KWRITE,KPUNCH,KDA(9)                        00000280
      COMMON /COBSER/ DELO(5),JOBS,JEND,ELAB(51),KOBS,MELA,Q0OBS(51),   00000290
     1                WMASS2,GAMMO,ETAJ(5),CSB,ISOSYM                   00000300
C                                                                       00000310
      LOGICAL INPUT/.FALSE./                                            00000320
      logical firstl/.true./
      LOGICAL CSB,ISOSYM,INELAS                                         00000330
C                                                                       00000340
      CHARACTER*4 NNAME(17)                                             00000350
C                                                                       00000360
      DATA PI/3.141592653589793D0/                                      00000540
C                                                                       00000370
C                                                                       00000380
C                                                                       00000390
C     DEFINE INPUT PARAMETERS ONLY ONCE                                 00000400
C                                                                       00000410
      IF(INPUT) GO TO 100                                               00000420
      INPUT=.TRUE.                                                      00000430
C                                                                       00000440
C                                                                       00000450
C******************* INPUT PARAMETERS ***************************       00000460
C                                                                       00000470
      JB = 0                                                            00000480
      JE = 20                                                           00000490
      WMASS1 = 938.926d0!938.919D0                                                00000500
      WMASS2 = 938.926d0!938.919D0                                                00000510
C                                                                       00000520
C     NUMBER OF ELABS                                                   00000530
C                                                                       00000540
      MELA = 17                                                          00000550
      IRAD = 0                                                          00000560
      ICSB = 0                                                          00000570
c isym=1: pp ; isym=0: np
      ISYM = 0                                                          00000580
      INELA =1                                                          00000590
C                                                                       00000600
      KREAD=5                                                           00000610
      KWRITE=6                                                          00000620
      open (unit=5,file='TobarePaul.d')
      open (unit=6,file='obsdat.d')
C                                                                       00000630
C****************************************************************       00000640
      WRITE (KWRITE,10000) JB,JE,WMASS1,WMASS2,MELA,IRAD                00000650
C                                                                       00000660
C                                                                       00000670
      CSB=.FALSE.                                                       00000680
      ISOSYM=.FALSE.                                                    00000690
      INELAS=.FALSE.                                                    00000700
      IF (ICSB.NE.0) CSB=.TRUE.                                         00000710
      IF (ISYM.NE.0) ISOSYM=.TRUE.                                      00000720
      IF (INELA.NE.0) INELAS=.TRUE.                                     00000730
C                                                                       00000740
C                                                                       00000750
      WRITE (KWRITE,10010) ICSB,ISYM,INELA                              00000760
      WRITE (KWRITE,10004)                                              00000770
      read (kread,10030) nname
      WRITE (KWRITE,10006) nname                                        00000780
      IF (IRAD.NE.0) GO TO 80                                           00000790
      WRITE (KWRITE,10001)                                              00000800
      GO TO 90                                                          00000810
   80 WRITE (KWRITE,10005)                                              00000820
   90 JEND=JE                                                           00000830
      JB1=JB+1                                                          00000840
      JE1=JE+1                                                          00000850
      RD=PI/180.D0                                                      00000860
C                                                                       00000870
      WW=2.D0*WMASS1*WMASS1*WMASS2                                      00000880
      WV=WMASS1+WMASS2                                                  00000890
      WV=WV*WV                                                          00000900
      WVW=WW/WV                                                         00000910
C                                                                       00000920
      GAMMO=0.D0                                                        00000930
      DO 1 I=1,4                                                        00000940
    1 ETAJ(I)=1.D0                                                      00000950
      ETAJ(5)=0.D0                                                      00000960
  100 CONTINUE                                                          00000970
C                                                                       00000980
C                                                                       00000990
10000 FORMAT(' INPUT PARAMETERS FOR READPH'/1H ,28(1H-)//5X,            00001000
     1       'JB =',I3,'  JE =',I3//5X,'W1 =',F9.3,'  W2 =',F9.3//2X,   00001010
     2       '   NUMBER OF ELABS  ',I3//5X,'IRAD =',I3)                 00001020
10001 FORMAT(//' P H A S E - S H I F T S (DEGREES)'/1H ,33(1H-))        00001030
10005 FORMAT(//' P H A S E - S H I F T S (RADIANS)'/1H ,33(1H-))        00001040
10002 FORMAT(1X,F8.2,5F14.6)                                            00001050
10004 FORMAT(//8(1H+),' OBSERVABLES ARE CALCULATED WITH PHASE-SHIFTS GIV
     1EN BELOW ',8(1H+))                                                00001080
10006 FORMAT(//' NUCLEON-NUCLEON DATASET :   ',17a4)                    00001090
10003 FORMAT(/' J =',I3/' ELAB(MEV)',5X,'SINGLET',7X,'TRIPLET',         00001100
     1       8X,'L=J-1',9X,'L=J+1',9X,'EPS')                            00001110
10007 FORMAT(/' J =',I3/'0ELAB(MEV)',5X,'SINGLET',7X,'TRIPLET',         00001120
     1       8X,'L=J-1',9X,'L=J+1',9X,'EPS',9X,'GA')                    00001130
10010 FORMAT(/5X,'ICSB =',I3//5X,'ISYM =',I3//5X,'INELA =',I3)          00001140
10020 FORMAT(1H ,F8.2,5F14.6/1H ,8X,F14.6)                              00001150
10021 FORMAT(1H ,F8.2,5D14.6)                                           00001160
10022 FORMAT(1H ,8X,5D14.6)                                             00001170
10030 FORMAT(17A4)                                                      00001180
C                                                                       00001190
C                                                                       00001200
C                                                                       00001210
C                                                                       00001220
C                                                                       00001230
C                                                                       00001240
C     READ PHASESHIFTS FROM DATASET AND                                 00001250
C          CALL OBS                                                     00001260
C                                                                       00001270
C     DO LOOP OF J                                                      00001280
C                                                                       00001290
      DO 200 J1=JB1,JE1                                                 00001300
C                                                                       00001310
      nkr=3
      if (firstl) then
      nkr=2
      firstl=.false.
      end if
c
      DO 101 KR=1,nkr                                                   00001320
  101 READ (KREAD,10030) NNAME                                          00001330
C                                                                       00001340
      JOBS=J1-1                                                         00001350
      IF (CSB) GO TO 110                                                00001360
      WRITE (KWRITE,10003) JOBS                                         00001370
      GO TO 120                                                         00001380
C                                                                       00001390
  110 WRITE (KWRITE,10007) JOBS                                         00001400
C                                                                       00001410
C     DO LOOP OF ELAB                                                   00001420
C                                                                       00001430
  120 DO 300 KOBS=1,MELA                                                00001440
C                                                                       00001450
      IF (CSB) GO TO 150                                                00001460
      IF (INELAS) GO TO 140
C                                                                       00001480
      READ (KREAD,10002) ELAB(KOBS),(DELO(I),I=1,5)                     00001490
      WRITE (KWRITE,10002) ELAB(KOBS),(DELO(I),I=1,5)                   00001500
      GO TO 190                                                         00001510
C                                                                       00001520
C                                                                       00001530
  140 READ (KREAD,10021) ELAB(KOBS),(DELO(I),I=1,5)                     00001540
      WRITE (KWRITE,10021) ELAB(KOBS),(DELO(I),I=1,5)                   00001550
      READ (KREAD,10022) (ETAJ(I),I=1,5)                                00001560
      WRITE (KWRITE,10022) (ETAJ(I),I=1,5)                              00001570
      GO TO 190                                                         00001580
C                                                                       00001590
C                                                                       00001600
  150 READ (KREAD,10020) ELAB(KOBS),(DELO(I),I=1,5),GAMMO               00001610
      WRITE (KWRITE,10020) ELAB(KOBS),(DELO(I),I=1,5),GAMMO             00001620
C                                                                       00001630
  190 IF (IRAD.NE.0) GO TO 202                                          00001640
C                                                                       00001650
C                                                                       00001660
C     PHASESHIFTS IN RAD                                                00001670
C                                                                       00001680
      DO 201 K=1,5                                                      00001690
      DELO(K)=RD*DELO(K)                                                00001700
  201 CONTINUE                                                          00001710
      IF (.NOT.CSB) GO TO 202                                           00001720
      GAMMO=GAMMO*RD                                                    00001730
  202 CONTINUE                                                          00001740
C                                                                       00001750
C     CALCULATE Q0OBS                                                   00001760
C                                                                       00001770
      Q0OBS(KOBS)=DSQRT(WVW*ELAB(KOBS))                                 00001780
C                                                                       00001790
C                                                                       00001800
      CALL OBS                                                          00001810
C                                                                       00001820
  300 CONTINUE                                                          00001830
  200 CONTINUE                                                          00001840
      END                                                               00001850
      SUBROUTINE OBS                                                    00000010
C                                                                       00000020
C                                                                       00000030
C         ROUTINE COMPUTES SCATTERING OBSERVABLES                       00000040
C                                                                       00000050
C                                                                       00000060
C                                                                       00000070
C         INFORMATIONS ABOUT THE VARIABLES USED:                        00000080
C                                                                       00000090
C           D-FUNCTIONS:   D(LL,3,2) = D1,0                             00000100
C                          D(LL,3,3) = D1,1                             00000110
C                          D(LL,1,3) = D-1,1                            00000120
C                                                                       00000130
C           COMPLEX R-MATRIX ELEMENTS IN LSJ-REPRESENTATION:            00000140
C                          RS  --- R0  --- RJ                           00000150
C                          RT0 --- R1  --- RJJ                          00000160
C                          RTM --- R+  --- RJ+1,J                       00000170
C                          RTP --- R-  --- RJ-1,J                       00000180
C                          RT  --- R+-                                  00000190
C                          RST --- RST                                  00000200
C                                                                       00000210
C                                                                       00000220
C           DELO(1) ... DELO(4) :   PHASE-SHIFTS                        00000230
C           ETAJ(1) ... ETAJ(4) :   INELASTICITY COEFFICIENTS           00000240
C                       DELO(5) :   TENSOR MIXING-PARAMETER (REAL PART) 00000250
C                       ETAJ(5) :   TENSOR MIXING-PARAMETER (IMAGINARY  00000260
C                                   PART)                               00000270
C                       GAMMO   :   SINGLET-TRIPLET MIXING-PARAMETER    00000280
C                                                                       00000290
C                                                                       00000300
C                                                                       00000310
C                                                                       00000320
C                                                                       00000330
      IMPLICIT COMPLEX*16 (Z), REAL*8 (A-H, O-Q, S-Y), INTEGER (I-N)    00000340
C                                                                       00000350
      COMMON  /CRDWRT/ KREAD,KWRITE,KPUNCH,KDA(9)                       00000360
      COMMON  /COBSER/ DELO(5),JOBS,JEND,ELAB(51),KOBS,MELA,Q0OBS(51),  00000370
     1                 WMASS2,GAMMO,ETAJ(5),CSB,ISOSYM                  00000380
C                                                                       00000390
      COMPLEX*16 RS,RTP,RT0,RTM,RT,RST,MS,MT0,A2                        00000400
      COMPLEX*16 ZPHI1(92,51), ZPHI2(92,51), ZPHI3(92,51), ZPHI4(92,51),00000410
     1           ZPHI5(92,51),ZPHI7(92,51)                              00000420
      COMPLEX*16 MTM(92,51),MTP(92,51),MT(92,51)                        00000430
      REAL*8 THETA(181),P(92),D(92,3,3),OBSERV(15,181,51),SI0TOT(51),   00000440
     1       SI1TOT(51),SI2TOT(51),DELLON(51),DELTRA(51),SITOAV(51),    00000450
     2       PLAB(51)                                                   00000460
      REAL*8 AMSQ(92,51),AMT0Q(92,51),DELAST(51),DELASL(51),SIGEAV(51), 00000470
     1       DINELT(51),DINELL(51),DINELA(51)                           00000480
C                                                                       00000490
      LOGICAL CSB,ISOSYM                                                00000500
      LOGICAL INDTHE /.FALSE./                                          00000510
C                                                                       00000520
      DATA UF/197.3270D0/                                               00000530
      DATA PI/3.141592653589793D0/                                      00000540
C                                                                       00000550
C                                                                       00000560
C                                                                       00000570
10001 FORMAT(//15X,15(1H-),'SCATTERING OBSERVABLES',15(1H-))            00000580
10002 FORMAT(//' LABOR ENERGY =',F9.3,'MEV')                            00000590
10003 FORMAT(/' THETA(DEGREES)',5X,'DSG IO (1)',2X,                     00000600
     1       'POLAR. PN0 (2)',3X,'DEPOL. DNN (3)',5X,                   00000610
     2       'IOP(13)')                                                 00000620
10004 FORMAT(5f15.6) 
10005 FORMAT(/' THETA(DEGREES)',8X,'R (4)',10X,'A (5)',10X,             00000640
     1       'RS (6)',9X,'AS (7)')                                      00000650
10006 FORMAT(/' THETA(DEGREES)',8X,'CKP (8)',8X,'CNN (9)',8X,           00000660
     1       'KNN (10)',7X,'RT (11)')                                   00000670
10007 FORMAT(5f15.6) 
10008 FORMAT(//20X,'TOTAL CROSS SECTIONS IN MBARN',/20X,29(1H=))        00000690
10009 FORMAT(//' SIGTOT = SIG0TOT + SIG1TOT*(P1.P2) + SIG2TOT*(P1.K 1)(
     1P2.K)')                                                           00000710
10010 FORMAT(//' LAB ENERGY(MEV)',2X,'LAB MOMENTUM(MEV)'                00000720
     1      ,3X,'SIG0TOT',8X,'SIG1TOT',8X,'SIG2TOT')                    00000730
C0011 FORMAT(9X,F7.3,15X,F7.3,15X,D12.4,9X,D12.4,9X,D12.4)              00000740
10011 FORMAT(5f15.6) 
10012 FORMAT(//' DELSIGL: DIFFERENCE BETWEEN TOTAL CROSS SECTIONS'/     00000760
     1 10x,'FOR ANTIPARALLEL AND PARALLEL SPIN STATES (LONGITUDINAL)'   00000770
     2          /' DELSIGT: DIFFERENCE BETWEEN TOTAL CROSS SECTIONS'/   00000780
     3 10x,'FOR ANTIPARALLEL AND PARALLEL SPIN STATES (TRANSVERSE)'     00000790
     4          /' SIGAVER: SPIN AVERAGED TOTAL CROSS SECTION '/        00000800
     5//15X,'( E L A S T I C  P L U S  I N E L A S T I C )')            00000810
10013 FORMAT(/' LAB ENERGY(MEV)',2X,'LAB MOMENTUM(MEV)',2X,'DELSIGL',
     1 8X,'DELSIGT',8X,'SIGAVER')                                       00000830
10014 FORMAT(//25X,'( E L A S T I C )')                                 00000840
10015 FORMAT(//23X,'( I N E L A S T I C )')                             00000850
10020 FORMAT(//' THETA(DEGREES)',8X,'AT (12)',7X,'P0N (14)',3X,         00000860
     1      'ANAL.-POW. DIFF. (PN0-P0N) (15)')                          00000870
10021 FORMAT(5f15.6) 
C                                                                       00000890
C                                                                       00000900
C                                                                       00000910
C                                                                       00000920
C     DEFINE THETA ONLY ONCE  (MAX. 181)                                00000930
C      (IN DEGREES)                                                     00000940
C                                                                       00000950
      IF (INDTHE) GO TO 10                                              00000960
      INDTHE=.TRUE.                                                     00000970
C                                                                       00000980
C                                                                       00000990
      ITMAX=11                                                          00001000
C                                                                       00001010
      THETA(1)=10.D0
      THETA(2)=20.D0
      THETA(3)=30.D0                                                    00001020
      THETA(4)=35.D0                                                    00001020
      THETA(5)=40.D0                                                    00001030
      THETA(6)=45.D0                                                    00001030
      THETA(7)=50.D0                                                    00001040
      THETA(8)=60.D0                                                    00001050
      THETA(9)=70.D0                                                    00001060
      THETA(10)=80.D0                                                   00001070
      THETA(11)=90.D0                                                   00001080
C                                                                       00001270
C                                                                       00001280
C                                                                       00001290
   10 CONTINUE                                                          00001300
C                                                                       00001310
C                                                                       00001320
C                                                                       00001330
C                                                                       00001340
C                                                                       00001350
      ALAMB=WMASS2                                                      00001360
      ALAMBQ=ALAMB*ALAMB                                                00001370
      PRO=1.D0                                                          00001380
      IF (ISOSYM) PRO=2.D0                                              00001390
      JEND1=JEND+1                                                      00001400
      JEND2=JEND+2                                                      00001410
C                                                                       00001420
C     CALCULATE MATRIX ELEMENTS FROM PHASESHIFTS                        00001430
C                                                                       00001440
C                                                                       00001450
      Z2=(0.D0,2.D0)*DELO(1)                                            00001460
      Z3=(0.D0,2.D0)*DELO(3)                                            00001470
      Z4=(0.D0,2.D0)*DELO(2)                                            00001480
      Z5=(0.D0,2.D0)*DELO(4)                                            00001490
      Z6=(0.D0,1.D0)*(DELO(3)+DELO(4))                                  00001500
      ZEP=DCMPLX(DELO(5),ETAJ(5))                                       00001510
      A2=2.D0*ZEP                                                       00001520
C                                                                       00001530
      IF (.NOT.CSB) GO TO 20                                            00001540
      Z7=(0.D0,1.D0)*(DELO(1)+DELO(2))                                  00001550
      A3=2.D0*GAMMO                                                     00001560
      RS=DCOS(A3)*CDEXP(Z2)-1.D0                                        00001570
      RT0=DCOS(A3)*CDEXP(Z4)-1.D0                                       00001580
      RST=-(0.D0,1.D0)*DSIN(A3)*CDEXP(Z7)                               00001590
      GO TO 30                                                          00001600
C                                                                       00001610
   20 RS=ETAJ(1)*CDEXP(Z2)-1.D0                                         00001620
      RT0=ETAJ(2)*CDEXP(Z4)-1.D0                                        00001630
C                                                                       00001640
   30 RTM=ETAJ(4)*CDCOS(A2)*CDEXP(Z5)-1.D0                              00001650
      RTP=ETAJ(3)*CDCOS(A2)*CDEXP(Z3)-1.D0                              00001660
      RT=-(0.D0,1.D0)*DSQRT(ETAJ(3)*ETAJ(4))*CDSIN(A2)*CDEXP(Z6)        00001670
C                                                                       00001680
C                                                                       00001690
      IF (.NOT.ISOSYM) GO TO 35                                         00001700
      IF (MOD(JOBS,2).EQ.0) GO TO 31                                    00001710
      ASY1=0.D0                                                         00001720
      ASY2=1.D0                                                         00001730
      GO TO 32                                                          00001740
   31 ASY1=1.D0                                                         00001750
      ASY2=0.D0                                                         00001760
   32 RS=ASY1*RS                                                        00001770
      RT0=ASY2*RT0                                                      00001780
      RTM=ASY1*RTM                                                      00001790
      RTP=ASY1*RTP                                                      00001800
      RT=ASY1*RT                                                        00001810
C                                                                       00001820
C                                                                       00001830
C        COMPUTE AND STORE M-MATRIX ELEMENTS ACCORDING TO W.M.KLOET     00001840
C        ET AL.,NUCL.PHYS. A364, 346 (1981) FOR CALCULATION OF          00001850
C        ELASTIC CROSS-SECTIONS                                         00001860
C                                                                       00001870
   35 MS=-(0.D0,0.5D0)*RS                                               00001880
      MT0=-(0.D0,0.5D0)*RT0                                             00001890
C                                                                       00001900
      J1OBS=JOBS+1                                                      00001910
      MTM(J1OBS,KOBS)=-(0.D0,0.5D0)*RTM                                 00001920
      MTP(J1OBS,KOBS)=-(0.D0,0.5D0)*RTP                                 00001930
      MT(J1OBS,KOBS)=-(0.D0,0.5D0)*RT                                   00001940
C                                                                       00001950
      AMSQ(J1OBS,KOBS)=DCONJG(MS)*MS                                    00001960
      AMT0Q(J1OBS,KOBS)=DCONJG(MT0)*MT0                                 00001970
C                                                                       00001980
C                                                                       00001990
C                                                                       00002000
C     TRANSFORM R-MATRIX ELEMENTS INTO HELICITY BASIS                   00002010
C                                                                       00002020
      AJ=DFLOAT(JOBS)                                                   00002030
      AJ1=AJ+1.D0                                                       00002040
      A2J1=2.D0*AJ+1.D0                                                 00002050
      AJJ=DSQRT(AJ*AJ1)                                                 00002060
C                                                                       00002070
      ZPHI12=(AJ1*RTM+AJ*RTP-2.D0*AJJ*RT)/A2J1                          00002080
      ZPHI34=(AJ*RTM+AJ1*RTP+2.D0*AJJ*RT)/A2J1                          00002090
      ZPHI57=(AJJ*(RTP-RTM)-RT)/A2J1                                    00002100
C                                                                       00002110
C     STORE MATRIX ELEMENTS                                             00002120
C                                                                       00002130
      ZPHI1(J1OBS,KOBS)=0.5D0*(RS+ZPHI12)                               00002140
      ZPHI2(J1OBS,KOBS)=0.5D0*(-RS+ZPHI12)                              00002150
      ZPHI3(J1OBS,KOBS)=0.5D0*(RT0+ZPHI34)                              00002160
      ZPHI4(J1OBS,KOBS)=0.5D0*(-RT0+ZPHI34)                             00002170
C                                                                       00002180
      IF (.NOT.CSB) GO TO 40                                            00002190
      ZPHI5(J1OBS,KOBS)=0.5D0*(-RST+ZPHI57)                             00002200
      ZPHI7(J1OBS,KOBS)=0.5D0*(RST+ZPHI57)                              00002210
      GO TO 50                                                          00002220
C                                                                       00002230
   40 ZPHI5(J1OBS,KOBS)=0.5D0*ZPHI57                                    00002240
C                                                                       00002250
C                                                                       00002260
C     CONTINUE ONLY IF ALL PARTIAL-WAVES ARE PRESENT                    00002270
C                                                                       00002280
C                                                                       00002290
   50 IF (JOBS.EQ.JEND.AND.KOBS.EQ.MELA) GO TO 100                      00002300
      GO TO 999                                                         00002310
C                                                                       00002320
C                                                                       00002330
  100 CONTINUE                                                          00002340
C                                                                       00002350
C                                                                       00002360
C     LOOP OF ELABS                                                     00002370
      ISTOT=0                                                           00002380
  400 DO 500 KELA=1,MELA                                                00002390
C                                                                       00002450
C                                                                       00002460
C     LOOP OF ANGLES                                                    00002470
C                                                                       00002480
      DO 600 ITE=1,ITMAX                                                00002490
C                                                                       00002400
C     IF ISTOT=1 CALCULATE TOTAL CROSS SECTIONS                         00002410
C     (BASED ON THE OPTICAL THEOREM)                                    00002420
C                                                                       00002430
      IF (ISTOT.EQ.1) GO TO 750                                         00002440
C                                                                       00002500
C     ANGLES IN RADIANDS                                                00002510
C                                                                       00002520
      THET=THETA(ITE)*PI/180.D0                                         00002530
      U=DCOS(THET)                                                      00002540
c
      V=DSIN(THET)                                                      00002550
      THETH=0.5D0*THET                                                  00002560
      UU=DCOS(THETH)                                                    00002570
      VV=DSIN(THETH)                                                    00002580
C                                                                       00002590
C                                                                       00002600
C     CALCULATE LEGENDRE POLYNOMIALS AND D-FUNCTIONS                    00002610
C                                                                       00002620
      P(1)=0.D0                                                         00002630
      P(2)=1.D0                                                         00002640
C                                                                       00002650
      DO 701 LL=2,JEND2                                                 00002660
      AL=LL-2                                                           00002670
      P(LL+1)=((2.D0*AL+1.D0)*U*P(LL)-AL*P(LL-1))/(AL+1.D0)             00002680
c
  701 CONTINUE                                                          00002690
C                                                                       00002700
C                                                                       00002710
      DO 702 LL=2,JEND2                                                 00002720
      AL=LL-2                                                           00002730
      A1=AL*(AL+1.D0)                                                   00002740
      A2=DSQRT(A1)                                                      00002750
      D(LL,3,2)=A2*(P(LL+1)-P(LL-1))/((2.D0*AL+1.D0)*V)                 00002760
      D(LL,3,3)=(P(LL)+(AL*P(LL+1)+(AL+1.D0)*P(LL-1))/(2.D0*AL+1.D0))/  00002770
     1(1.D0+U)                                                          00002780
      D(LL,1,3)=(-P(LL)+(AL*P(LL+1)+(AL+1.D0)*P(LL-1))/(2.D0*AL+1.D0))/ 00002790
     1(1.D0-U)                                                          00002800
  702 CONTINUE                                                          00002810
      GO TO 760                                                         00002820
C                                                                       00002830
C                                                                       00002840
C                                                                       00002850
  750 U=1.D0                                                            00002860
      V=0.D0                                                            00002870
C                                                                       00002880
      DO 751 LL=2,JEND2                                                 00002890
      P(LL)=1.D0                                                        00002900
      D(LL,3,3)=1.D0                                                    00002910
      D(LL,1,3)=0.D0                                                    00002920
  751 D(LL,3,2)=0.D0                                                    00002930
C                                                                       00002940
C                                                                       00002950
C     CALCULATE HELICITY AMPLITUDES                                     00002960
C                                                                       00002970
  760 ZZPHI1=(0.D0,0.D0)                                                00002980
      ZZPHI2=(0.D0,0.D0)                                                00002990
      ZZPHI3=(0.D0,0.D0)                                                00003000
      ZZPHI4=(0.D0,0.D0)                                                00003010
      ZZPHI5=(0.D0,0.D0)                                                00003020
      IF (ISTOT.EQ.0) GO TO 765                                         00003030
      XSIN=0.D0                                                         00003040
      XTRI=0.D0                                                         00003050
      XMIX=0.D0                                                         00003060
  765 IF (.NOT.CSB) GO TO 790                                           00003070
      ZZPHI7=(0.D0,0.D0)                                                00003080
C                                                                       00003090
  790 DO 800 JJ=1,JEND1                                                 00003100
      JJ1=JJ+1                                                          00003110
      AJX=DFLOAT(JJ-1)                                                  00003120
      AJX1=AJX+1.D0                                                     00003130
      A2JX1=2.D0*AJX+1.D0                                               00003140
      A2JX1H=A2JX1/2.D0                                                 00003150
      AJXJ1=DSQRT(AJX*AJX1)                                             00003160
C                                                                       00003170
      ZZPHI1=ZZPHI1+A2JX1H*ZPHI1(JJ,KELA)*P(JJ1)                        00003180
      ZZPHI2=ZZPHI2+A2JX1H*ZPHI2(JJ,KELA)*P(JJ1)                        00003190
      ZZPHI3=ZZPHI3+A2JX1H*ZPHI3(JJ,KELA)*D(JJ1,3,3)                    00003200
      ZZPHI4=ZZPHI4+A2JX1H*ZPHI4(JJ,KELA)*D(JJ1,1,3)                    00003210
      ZZPHI5=ZZPHI5+A2JX1H*ZPHI5(JJ,KELA)*D(JJ1,3,2)                    00003220
C                                                                       00003230
      IF (ISTOT.EQ.0) GO TO 795                                         00003240
      AMTMQ=DCONJG(MTM(JJ,KELA))*MTM(JJ,KELA)                           00003250
      AMTPQ=DCONJG(MTP(JJ,KELA))*MTP(JJ,KELA)                           00003260
      AMTQ=DCONJG(MT(JJ,KELA))*MT(JJ,KELA)                              00003270
C                                                                       00003280
      XSIN=XSIN+A2JX1*AMSQ(JJ,KELA)                                     00003290
      XTRI=XTRI+A2JX1*(AMT0Q(JJ,KELA)+AMTPQ+AMTMQ+2.D0*AMTQ)            00003300
      XMIX=XMIX+AJX*AMTPQ+AJX1*AMTMQ+A2JX1*AMTQ-2.D0*AJXJ1*             00003310
     1     DREAL((MTP(JJ,KELA)+MTM(JJ,KELA))*DCONJG(MT(JJ,KELA)))       00003320
C                                                                       00003330
  795 IF (.NOT.CSB) GO TO 800                                           00003340
      ZZPHI7=ZZPHI7-A2JX1H*ZPHI7(JJ,KELA)*D(JJ1,3,2)                    00003350
  800 CONTINUE                                                          00003360
C                                                                       00003370
      ZZPHI1=-(0.D0,1.D0)*ZZPHI1/Q0OBS(KELA)                            00003380
      ZZPHI2=-(0.D0,1.D0)*ZZPHI2/Q0OBS(KELA)                            00003390
      ZZPHI3=-(0.D0,1.D0)*ZZPHI3/Q0OBS(KELA)                            00003400
      ZZPHI4=-(0.D0,1.D0)*ZZPHI4/Q0OBS(KELA)                            00003410
      ZZPHI5=-(0.D0,1.D0)*ZZPHI5/Q0OBS(KELA)                            00003420
      IF (.NOT.CSB) GO TO 810                                           00003430
      ZZPHI7=-(0.D0,1.D0)*ZZPHI7/Q0OBS(KELA)                            00003440
C                                                                       00003450
  810 ZFAC=4.D0*ZZPHI5                                                  00003460
      IF (CSB) ZFAC=2.D0*(ZZPHI5-ZZPHI7)                                00003470
C                                                                       00003480
C                                                                       00003490
C     CALCULATE SCATTERING OBSERVABLES                                  00003500
C                                                                       00003510
      ZA=.5D0*((ZZPHI1+ZZPHI2+ZZPHI3-ZZPHI4)*U-ZFAC*V)                  00003520
      ZB=.5D0*(ZZPHI1-ZZPHI2+ZZPHI3+ZZPHI4)                             00003530
      ZC=.5D0*(-ZZPHI1+ZZPHI2+ZZPHI3+ZZPHI4)                            00003540
      ZD=.5D0*(ZZPHI1+ZZPHI2-ZZPHI3+ZZPHI4)                             00003550
      IF (ISTOT.EQ.1) GO TO 850                                         00003560
      ZE=-(0.D0,.5D0)*((ZZPHI1+ZZPHI2+ZZPHI3-ZZPHI4)*V                  00003570
     1+ZFAC*U)                                                          00003580
C                                                                       00003590
      ZF=(0.D0,0.D0)                                                    00003600
      IF (CSB) ZF=(0.D0,1.D0)*(ZZPHI5+ZZPHI7)                           00003610
C                                                                       00003620
      AA=CDABS(ZA)**2                                                   00003630
      BB=CDABS(ZB)**2                                                   00003640
      CC=CDABS(ZC)**2                                                   00003650
      DD=CDABS(ZD)**2                                                   00003660
      EE=CDABS(ZE)**2                                                   00003670
      FF=CDABS(ZF)**2                                                   00003680
C                                                                       00003690
      AI0=.5D0*(AA+BB+CC+DD+EE+FF)                                      00003700
C                                                                       00003710
      Z1=DCONJG(ZA)*ZB+DCONJG(ZC)*ZD-DCONJG(ZE)*ZF                      00003720
      DKK=DREAL(Z1)/AI0                                                 00003730
      Z2=DCONJG(ZB)*ZE+DCONJG(ZA)*ZF                                    00003740
      DPK=-DIMAG(Z2)/AI0                                                00003750
      Z3=DCONJG(ZA)*ZB-DCONJG(ZC)*ZD-DCONJG(ZE)*ZF                      00003760
      DPP=DREAL(Z3)/AI0                                                 00003770
      AR=DKK*UU+DPK*VV                                                  00003780
      A=-DKK*VV+DPK*UU                                                  00003790
      Z4=DCONJG(ZD)*ZE-DCONJG(ZC)*ZF                                    00003800
      CKP=DIMAG(Z4)/AI0                                                 00003810
      CNN=.5D0*(AA-BB-CC+DD+EE-FF)/AI0                                  00003820
      DNN=.5D0*(AA+BB-CC-DD+EE+FF)/AI0                                  00003830
      Z5=DCONJG(ZA)*ZE+DCONJG(ZB)*ZF                                    00003840
      PN0=DREAL(Z5)/AI0                                                 00003850
      Z9=DCONJG(ZA)*ZE-DCONJG(ZB)*ZF                                    00003860
      P0N=DREAL(Z9)/AI0                                                 00003870
      Z6=DCONJG(ZB)*ZE+DCONJG(ZA)*ZF                                    00003880
      DKP=DIMAG(Z6)/AI0                                                 00003890
      ARS=DKP*UU+DPP*VV                                                 00003900
      AS=-DKP*VV+DPP*UU                                                 00003910
      AKNN=.5D0*(AA-BB+CC-DD+EE-FF)/AI0                                 00003920
      Z7=DCONJG(ZC)*ZE-DCONJG(ZD)*ZF                                    00003930
      AKKP=DIMAG(Z7)/AI0                                                00003940
      Z8=DCONJG(ZA)*ZC-DCONJG(ZB)*ZD                                    00003950
      AKPP=DREAL(Z8)/AI0                                                00003960
      ART=AKKP*UU+AKPP*VV                                               00003970
      AT=-AKKP*VV+AKPP*UU                                               00003980
C                                                                       00003990
C     DIFF. CROSS SECTION IN MB                                         00004000
C                                                                       00004010
c****  factor of 4 for pp diff. cross sections
      AI0=10.D0*UF*UF*AI0*pro*pro                                       00004020
C                                                                       00004030
      OBSERV(1,ITE,KELA)=AI0                                            00004040
      OBSERV(2,ITE,KELA)=PN0                                            00004050
      OBSERV(3,ITE,KELA)=DNN                                            00004060
      OBSERV(4,ITE,KELA)=AR                                             00004070
      OBSERV(5,ITE,KELA)=A                                              00004080
      OBSERV(6,ITE,KELA)=ARS                                            00004090
      OBSERV(7,ITE,KELA)=AS                                             00004100
      OBSERV(8,ITE,KELA)=CKP                                            00004110
      OBSERV(9,ITE,KELA)=CNN                                            00004120
      OBSERV(10,ITE,KELA)=AKNN                                          00004130
      OBSERV(11,ITE,KELA)=ART                                           00004140
      OBSERV(12,ITE,KELA)=AT                                            00004150
      OBSERV(13,ITE,KELA)=AI0*PN0                                       00004160
      OBSERV(14,ITE,KELA)=P0N                                           00004170
      OBSERV(15,ITE,KELA)=PN0-P0N                                       00004180
C                                                                       00004190
  600 CONTINUE                                                          00004200
      GO TO 500                                                         00004210
C                                                                       00004220
C                                                                       00004230
C     COMPUTE TOTAL CROSS SECTIONS                                      00004240
C                                                                       00004250
C                                                                       00004260
  850 ZAB=ZA+ZB                                                         00004270
      ZCD=ZC+ZD                                                         00004280
      AZAB=DIMAG(ZAB)                                                   00004290
      AZCD=DIMAG(ZCD)                                                   00004300
      AZD =DIMAG(ZD)                                                    00004310
      AZAB=AZAB/Q0OBS(KELA)                                             00004320
      AZCD=AZCD/Q0OBS(KELA)                                             00004330
      AZD =AZD/Q0OBS(KELA)                                              00004340
      AZAB=2.D0*PI*AZAB                                                 00004350
      AZCD=2.D0*PI*AZCD                                                 00004360
      AZD =-4.D0*PI*AZD                                                 00004370
C                                                                       00004380
C     TOTAL CROSS SECTIONS IN MB                                        00004390
C                                                                       00004400
C          (TOTAL)                                                      00004410
C                                                                       00004420
      AZAB=AZAB*UF*UF*10.D0                                             00004430
      AZCD=AZCD*UF*UF*10.D0                                             00004440
      AZD =AZD*UF*UF*10.D0                                              00004450
C                                                                       00004460
C     FACTOR 2 IN CASE OF TWO PROTONS                                   00004470
C                                                                       00004480
      SI0TOT(KELA)=AZAB*PRO                                             00004490
      SI1TOT(KELA)=AZCD*PRO                                             00004500
      SI2TOT(KELA)=AZD*PRO                                              00004510
C                                                                       00004520
      DELLON(KELA)=-2.D0*(AZCD+AZD)*PRO                                 00004530
      DELTRA(KELA)=-2.D0*AZCD*PRO                                       00004540
      SITOAV(KELA)=AZAB*PRO                                             00004550
C                                                                       00004560
C                                                                       00004570
C         (ELASTIC AND INELASTIC)                                       00004580
C                                                                       00004590
      Q0OBSQ=Q0OBS(KELA)*Q0OBS(KELA)                                    00004600
      DELAST(KELA)=2.D0*PI*PRO*(XSIN-XMIX)/Q0OBSQ                       00004610
      DELASL(KELA)=2.D0*PI*PRO*(XSIN-XTRI+2.D0*XMIX)/Q0OBSQ             00004620
      SIGEAV(KELA)=PI*PRO*(XSIN+XTRI)/Q0OBSQ                            00004630
C                                                                       00004640
      DELAST(KELA)=DELAST(KELA)*UF*UF*10.D0                             00004650
      DELASL(KELA)=DELASL(KELA)*UF*UF*10.D0                             00004660
      SIGEAV(KELA)=SIGEAV(KELA)*UF*UF*10.D0                             00004670
C                                                                       00004680
      DINELT(KELA)=DELTRA(KELA)-DELAST(KELA)                            00004690
      DINELL(KELA)=DELLON(KELA)-DELASL(KELA)                            00004700
      DINELA(KELA)=SITOAV(KELA)-SIGEAV(KELA)                            00004710
C                                                                       00004720
  500 CONTINUE                                                          00004730
      ISTOT=ISTOT+1                                                     00004740
      IF (ISTOT.EQ.1) GO TO 400                                         00004750
C                                                                       00004760
C                                                                       00004770
C                                                                       00004780
C     WRITE RESULTS                                                     00004790
C                                                                       00004800
C                                                                       00004810
      WRITE (KWRITE,10001)                                              00004820
C                                                                       00004830
      DO 901 K=1,MELA                                                   00004840
      WRITE (KWRITE,10002) ELAB(K)                                      00004850
C                                                                       00004860
      II=0                                                              00004870
  900 DO 902 IT=1,ITMAX                                                 00004880
C                                                                       00004890
      IF (IT.EQ.1.AND.II.EQ.0) WRITE (KWRITE,10003)                     00004900
      IF (IT.EQ.1.AND.II.EQ.1) WRITE (KWRITE,10005)                     00004910
      IF (IT.EQ.1.AND.II.EQ.2) WRITE (KWRITE,10006)                     00004920
      IF (IT.EQ.1.AND.II.EQ.3) WRITE (KWRITE,10020)                     00004930
C                                                                       00004940
      IF (II.EQ.1) GO TO 910                                            00004950
      IF (II.EQ.2) GO TO 920                                            00004960
      IF (II.EQ.3) GO TO 930                                            00004970
C                                                                       00004980
      WRITE (KWRITE,10004) THETA(IT),OBSERV(1,IT,K),OBSERV(2,IT,K),     00004990
     1                         OBSERV(3,IT,K),OBSERV(13,IT,K)           00005000
      GO TO 902                                                         00005010
C                                                                       00005020
  910 WRITE (KWRITE,10004) THETA(IT),OBSERV(4,IT,K),OBSERV(5,IT,K),     00005030
     1                         OBSERV(6,IT,K),OBSERV(7,IT,K)            00005040
      GO TO 902                                                         00005050
C                                                                       00005060
  920 WRITE (KWRITE,10007) THETA(IT),OBSERV(8,IT,K),OBSERV(9,IT,K),     00005070
     1                         OBSERV(10,IT,K),OBSERV(11,IT,K)          00005080
      GO TO 902                                                         00005090
C                                                                       00005100
  930 WRITE (KWRITE,10021) THETA(IT),OBSERV(12,IT,K),OBSERV(14,IT,K),   00005110
     1                         OBSERV(15,IT,K)                          00005120
C                                                                       00005130
  902 CONTINUE                                                          00005140
      II=II+1                                                           00005150
      IF (II.GT.3) GO TO 901                                            00005160
      GO TO 900                                                         00005170
  901 CONTINUE                                                          00005180
C                                                                       00005190
C     WRITE TOTAL CROSS SECTIONS                                        00005200
      WRITE (KWRITE,10008)                                              00005210
      WRITE (KWRITE,10009)                                              00005220
      WRITE (KWRITE,10010)                                              00005230
      DO 950 I=1,MELA                                                   00005240
C                                                                       00005250
C     LAB MOMENTUM                                                      00005260
C                                                                       00005270
      AAA=ALAMB+ELAB(I)                                                 00005280
      AAAQ=AAA*AAA                                                      00005290
      ASD=AAAQ-ALAMBQ                                                   00005300
      PLAB(I)=DSQRT(ASD)                                                00005310
C                                                                       00005320
      WRITE (KWRITE,10011) ELAB(I),PLAB(I),SI0TOT(I),SI1TOT(I),SI2TOT(I)00005330
  950 CONTINUE                                                          00005340
C                                                                       00005350
C                                                                       00005360
      WRITE (KWRITE,10012)                                              00005370
      WRITE (KWRITE,10013)                                              00005380
C                                                                       00005390
      DO 960 I=1,MELA                                                   00005400
  960 WRITE (KWRITE,10011) ELAB(I),PLAB(I),DELLON(I),DELTRA(I),SITOAV(I)00005410
C                                                                       00005420
C                                                                       00005430
      WRITE (KWRITE,10014)                                              00005440
      WRITE (KWRITE,10013)                                              00005450
C                                                                       00005460
      DO 970 I=1,MELA                                                   00005470
  970 WRITE (KWRITE,10011) ELAB(I),PLAB(I),DELASL(I),DELAST(I),SIGEAV(I)00005480
C                                                                       00005490
C                                                                       00005500
      WRITE (KWRITE,10015)                                              00005510
      WRITE (KWRITE,10013)                                              00005520
C                                                                       00005530
      DO 980 I=1,MELA                                                   00005540
  980 WRITE (KWRITE,10011) ELAB(I),PLAB(I),DINELL(I),DINELT(I),DINELA(I)00005550
C                                                                       00005560
  999 RETURN                                                            00005570
      END  
