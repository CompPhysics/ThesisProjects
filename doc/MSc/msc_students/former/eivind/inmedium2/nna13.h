      subroutine nna13                                                   00000010
c
c       this is the code for a NN potential model that describes
c       NN scattering up to 1000 MeV quantitatively.
c
c       besides an OBE part that is similar to CD-Bonn,
c       this model includes N-Delta and Delta-Delta box
c       diagrams. the latter diagrams
c       are important to generate the inelasticity above
c       pion-production theshold (about 300 MeV) and to create
c       certain "resonance-like" structures in particular
c       partial waves (like, 1D2, 3F3, ...).
c       due to those Delta box diagrams, the model is
c       energy-dependent and the application of this
c       code is somewhat more involved than the application
c       of a primitive NN potential code. please ask the author
c       for technical help if needed.
c
c       this package is selfcontained and includes all
c       subroutines needed for the NN interaction.
c
c       the present version is a neutron-proton potential.
c
c
c       author: r. machleidt
c               dept. of physics
c               univ. of idaho
c               moscow, id 83844
c               u. s. a.
c               machleid@uidaho.edu
c
c       February 29, 2000
c
c
c
c                                                                       00000020
      implicit real*8(a,b,d-h,o-z), complex*16  (c)                     00000030
c                                                                       00000040
      common /crdwrt/ kread,kwrite,kpunch,kda(9)                        00000050
      common /cpot/   v(6),xmev,ymev                                    00000060
      common /cpotc/  cv(6),cxmev,cymev                                 00000070
      common /cstate/ j,heform,sing,trip,coup,endep,label               00000080
      common /cpoted/ q0qmev,qfmev,pmev,uanspp,wsnspp,ucnspp,udnspp,    00000090
     1                znrl,zrel,smev,noced                              00000100
c                                                                       00000110
      dimension cvv(6)                                                  00000120
      logical heform,sing,trip,coup,endep                               00000130
      logical noced                                                     00000140
      logical nnoced                                                    00000150
c
c
c**** statement functions **********
      real(cc)=dble(cc)
      cmplx(x,y)=dcmplx(x,y)
c
c                                                                       00000170
      do 95 iv=1,6                                                      00000180
   95 cvv(iv)=(0.d0,0.d0)                                               00000190
      nnoced=noced                                                      00000200
c                                                                       00000210
c                                                                       00000220
      xmev= real(cxmev)                                                 00000230
      ymev= real(cymev)                                                 00000240
c                                                                       00000340
c                                                                       00000350
c                                                                       00000360
c                                                                       00000370
c        call two-boson-exchange                                        00000380
c                                                                       00000390
c                                                                       00000400
      if (j.ge.30) go to 900
c                                                                       00000420
c                                                                       00000430
      noced=nnoced                                                      00000440
c                                                                       00000450
      call tbibdc                                                       00000460
c                                                                       00000470
      do 200 iv=1,6                                                     00000480
  200 cvv(iv)=cvv(iv)+cv(iv)                                            00000490
c                                                                       00000500
c                                                                       00000510
  900 continue
c
c
c        call one-boson-exchange                                        00000270
c                                                                       00000280
c                                                                       00000290
      noced=nnoced                                                      00000440
c
      call cdbond                                                       00000300
c                                                                       00000310
      do 100 iv=1,6                                                     00000320
  100 cvv(iv)=cvv(iv)+ cmplx(v(iv),0.d0)                                00000330
c                                                                       00000520
      endep=.true.
c                                                                       00000530
      do 1005 iv=1,6                                                    00000540
 1005 cv(iv)=cvv(iv)                                                    00000550
c
c
      return                                                            00000580
      end                                                               00000590
      subroutine cdbond
c
c        this is a special code to be used in conjunction
c        with a model for the NN interaction that includes
c        intermediate delta isobars. This code provides
c        the OBE part of that model.
c
c        the model describes NN scattering up to 1000 MeV lab energy
c        quantitatively.
c
c        the present version of this model ignores charge-dependence
c        and is strictly speaking a model for the
c
c               neutron-proton interaction.
c
c******************************************************************
c
c        this package is self-contained and includes
c        all subroutines needed.
c        only cdbond needs to be called by the user.
c        all codes are consistently in double precision.
c        when working on an UNIX system, it is crucial
c        to compile this code with the  -static  option.
c
c
c************************************************************************
c
c
c        author:      Ruprecht Machleidt
c                     department of physics
c                     university of idaho
c                     moscow, idaho 83844
c                     u. s. a.
c                     e-mail: machleid@uidaho.edu
c
c                     formerly:
c                     institut fuer theoretische kernphysik der
c                     universitaet bonn
c                     nussallee 14-16
c                     d - 5300  bonn, w. germany
c
c*********************************************************************
c
c        this version of the code uses the Legendre functions
c        ----------------------------------------------------
c        of the second kind for the partial wave decomposition
c        -----------------------------------------------------
c        and includes the meson-parameters in data statements.
c        -----------------------------------------------------
c
c
c        an earlier version of this code has been published in
c        "computational nuclear physics 2 -- nuclear reactions",
c        K. Langanke, J.A. Maruhn, and S.E. Koonin, eds.
c        (Springer, New York, 1993), Chapter 1, pp. 1-29.
c
c        This code is a slight modification of the earlier, published
c        version. However, the mathematical formulae, as well as the
c        general organization of this code is the same as described in
c        the above-referenced book-chapter.
c        in this version of the code, the integrals, Eqs. (1.68) of the
c        above reference, are solved analytically by means of the
c        Legendre functions of the second kind, see Eqs. (E.44) of
c        R. Machleidt et al., Phys. Rep. 149, 1 (1987).
c
c        Still, the above-referenced article may serve as a good
c        introduction into this code.
c
c*********************************************************************
c
c
      implicit real*8 (a-h,o-z)
c
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
c        arguments and values of this subroutine:
c
      common /cpot/   v(6),xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep,label
      common /cnn/ inn
c
c
c        this has been the end of the common-blocks containing
c        the arguments and values of this subroutine
c
c        specifications for these two common blocks
c
      logical heform,sing,trip,coup,endep
c
c        THE ABOVE FOUR COMMON BLOCKS IS ALL THE USER NEEDS
c        TO BE FAMILIAR WITH.
c
c
c        xmev and ymev are the final and initial relative momenta,
c        respectively, in units of mev/c.
c        v is the potential in units of mev**(-2).
c        concerning units and factor of pi etc.,
c        cf. with the partial-wave Lippmann-Schwinger equation, Eq. (1.32),
c        and with the phase shift relation, Eq. (1.41) of
c        R. Machleidt, in: Computational Nuclear Physics 2
c        -- Nuclear Reactions, Langanke et al., eds.
c        (Springer, New York, 1993), Chapter 1, pp. 1-29.
c
c        the partial-wave Lippmann-Schwinger equation for the
c        K-matrix reads:
c
c        K(q',q) = V(q',q) + M P \int dk k^2 V(q',k) K(k,q)/(q^2-k^2)
c
c        with M the nucleon mass in MeV and P denoting the principal value;
c        V(q',q) as provided by this code in common block /cpot/;
c        all momenta in MeV.
c
c        the phase-shift relation is:
c
c        tan \delta_L = -(pi/2) M q K_L(q,q)
c
c        with M and q in units of MeV, K_L in MeV**(-2) like V.
c
c
c        if heform=.true., v contains the 6 matrix elements
c        associated with one j in the helicity formalism
c        in the following order:
c        0v, 1v, 12v, 34v, 55v, 66v (for notation see above article).
c        if heform=.false., v contains the 6 matrix elements
c        associated with one j in the lsj formalism
c        in the following order:
c        0v(singlet), 1v(uncoupled triplet), v++, v--, v+-, v-+ (coupled)
c        (see above article for notation).
c        j is the total angular momentum. there is essentially no upper
c        limit for j.
c        sing, trip, and coup should in general be .true..
c        endep and label can be ignored.
c        it is customary, to set kread=5 and kwrite=6;
c        ignore kpunch and kda(9).
c
c
c        THIS IS ESSENTIALLY ALL THE USER NEEDS TO KNOW.
c
c**********************************************************************
c
c
c        common block for all ob-subroutines
c
      common /cobd/   vj(32,50),c(20,50),fff,ff,f(52),aa(96),ai(19,15),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                ez1,ez2,ct(96),wt(96),
     3                ic(20,50),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(12,3),mggo(12,3),ima(15,12,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(10,50),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c
c        further specifications
c
      data pi/3.141592653589793d0/
      character*4 mesong(12)
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     1            '1-  ','1-t ','1-tt','1-st','1-ss',
     2            '1+  ','2+  '/
      logical index
      data index/.false./
      data jj/-1/
      data innn/-1/
      dimension vl(4),adminv(4,4),ldminv(4),mdminv(4)
c
c
c
c
      inter=1
c
c       in its present form, this code is restricted to np;
c       therefore, no matter what the user specified,
c       we always set:
c
      inn=2
c
c
      if (inn.lt.1.or.inn.gt.3) then
      write (kwrite,19001) inn
19001 format (////' error in cdbonn: potential type  inn =',i3,
     1'  unknown.'/' execution terminated.'////)
      stop
      endif
      if (j.lt.0) then
      write (kwrite,19002)
19002 format (////' error in cdbonn: total angular momentum j',
     1' is negative.'/' execution terminated.'////)
      stop
      endif
c
c
c
c
c        call obpard whenever j or inn has changed
c
c
      if (j.eq.jj.and.inn.eq.innn) go to 50
      innn=inn
c
c
c
      call obpard
c     -----------
c     -----------
c
c
      dwn=1.d0/wnn(inter)
c
c        prepare constant over-all factor
c
      fac=1.d0/(2.d0*pi)*dwn*dwn
c     --------------------------
c
c
      if (index) go to 30
      index=.true.
c
      iftgo=ift(inter)+1
      iman=imaa(inter)
      imen=imea(inter)
c
      imanm1=iman-1
c
      iman1=imanm1+1
      iman2=imanm1+2
      iman3=imanm1+3
      iman4=imanm1+4
      iman5=imanm1+5
      iman6=imanm1+6
      iman7=imanm1+7
      iman8=imanm1+8
      iman9=imanm1+9
      iman10=imanm1+10
      iman11=imanm1+11
      iman12=imanm1+12
      iman13=imanm1+13
      iman14=imanm1+14
      iman15=imanm1+15
      iman16=imanm1+16
      iman17=imanm1+17
      iman18=imanm1+18
c
c
c
c
   30 if (j.eq.jj) go to 50
      jj=j
      if (j.eq.0) go to 50
      aj=dfloat(j)
      aj1=dfloat(j+1)
      a2j1=dfloat(2*j+1)
      aaj6=dsqrt(aj*aj1)
c
c        coefficient matrix for the translations into lsj formalism
c
      adminv(1,1)=aj1
      adminv(1,2)=aj
      adminv(1,3)=-aaj6
      adminv(1,4)=-aaj6
      adminv(2,1)=aj
      adminv(2,2)=aj1
      adminv(2,3)=aaj6
      adminv(2,4)=aaj6
      adminv(3,1)=aaj6
      adminv(3,2)=-aaj6
      adminv(3,3)=aj1
      adminv(3,4)=-aj
      adminv(4,1)=aaj6
      adminv(4,2)=-aaj6
      adminv(4,3)=-aj
      adminv(4,4)=aj1
c
c       inversion
c
      call dminv (adminv,4,deter,ldminv,mdminv)
c
c
c
c
c
c
c
c        prepare expressions depending on x and y
c        ----------------------------------------
c        ----------------------------------------
c
c
c
c
   50 x=xmev*dwn
      y=ymev*dwn
      indxy=.false.
c
c
      if (xmev.lt.0.d0) then
      write (kwrite,19003)
19003 format (////' error in cdbonn: momentum xmev',
     1' is negative.'/' execution terminated.'////)
      stop
      endif
      if (ymev.lt.0.d0) then
      write (kwrite,19004)
19004 format (////' error in cdbonn: momentum ymev',
     1' is negative.'/' execution terminated.'////)
      stop
      endif
c
c
      xx=x*x
      yy=y*y
      xy2=x*y*2.d0
      xxpyy=xx+yy
      ex=dsqrt(1.d0+xx)
      ey=dsqrt(1.d0+yy)
      eem12=(ex*ey-1.d0)*2.d0
c
c
c
c
      xy=xy2*0.5d0
      ee=ex*ey
      ree=dsqrt(ee)
      eem1=ee-1.d0
      eme=ex-ey
      emeh=eme*0.5d0
      emehq=emeh*emeh
      eep1=ee+1.d0
       epe=ex+ey
      xxyy=xx*yy
c
c
c
c
c        prepare over-all factor
c
c
      go to (70,71,72),iftgo
c
c        no additional factor
c
   70 fff=fac
      go to 90
c
c        minimal relativity
c
   71 fff=fac/ree
      go to 90
c
c        factor m/e*m/e
c
   72 fff=fac/ee
c
c
c
c
c
c
   90 do 93 iv=1,6
   93 v(iv)=0.d0
      do 95 il=1,50
      do 95 iv=1,6
   95 vj(iv,il)=0.d0
c
c
c
c
c        contributions of mesons
c        -----------------------
c        -----------------------
c
c
c
c
      do 1995 img=1,mge
      mg=mggo(img,inter)
      if (mg.eq.0) go to 2000
      if (mg.gt.7) go to 9000
      me=mgg(mg,inter)
      go to (100,9000,9000,400,9000,9000,700),mg
c
c
c
c        0-  , pseudo-scalar coupling
c        ----------------------------
c
c
c
c
  100 mc=1
c
      ff=1.d0
      f(1)=eem1
      f(2)=-xy
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-eme
      f(8)=-f(7)
c
      call obstrd(1,1,me)
      go to 1995
c
c
c
c
c        0+  , scalar coupling
c        ---------------------
c
c
c
c
  400 mc=1
c
      ff=1.d0
      f(1)=-eep1
      f(2)=xy
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=epe
      f(8)=f(7)
c
      call obstrd(1,1,me)
      go to 1995
c
c
c
c
c        1-t , vector mesons
c        -------------------
c
c
c
c
c        vector-vector coupling
c
c
c
c
  700 mc=1
c
      ff=2.d0
      f(1)=eem1+ee
      f(2)=0.d0
      f(3)=ee
      f(4)=xy
      f(5)=xy2
      f(6)=1.d0
      f(7)=-ey
      f(8)=-ex
c
      call obstrd(1,1,me)
c
c
c
c
c        tensor-tensor coupling
c
c
c
c
      mc=2
c
      ff=0.25d0
      f(1)=(3.d0*ee+1.d0)*xxpyy
      f(2)=-(6.d0*ee+2.d0-xxpyy)*xy
      f(3)=eem1*xxpyy+4.d0*xxyy
      f(4)=-(4.d0*ee+xxpyy)*xy
      f(5)=(4.d0-3.d0*xxpyy)*xy
      f(6)=6.d0*xxyy-(ee+3.d0)*xxpyy
      f(7)=(ex+3.d0*ey)*xx+eme*yy
      f(8)=(ey+3.d0*ex)*yy-eme*xx
c        factors for additional terms
      f(9)=-2.d0*xxyy
      f(10)=eep1*xy2
      f(11)=-epe*xy2
c
      call obstrd(2,1,me)
c
c
c
c
c        vector-tensor coupling
c
c
c
c
      mc=3
c
      ff=1.d0
      f(1)=xxpyy
      f(2)=-xy2
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=6.d0*xy
      f(6)=3.d0*f(3)
      f(7)=(ex*yy+3.d0*ey*xx)
      f(8)=(ey*xx+3.d0*ex*yy)
c
      call obstrd(1,1,me)
      go to 1995
c
c
c
c
c        this has been the end of the contributions of mesons
c        ----------------------------------------------------
c
c
c
c
c        error exit
c        ----------
c
c
c
c
 9000 write (kwrite,19000) mesong(mg)
19000 format(////' error in cdbonn:   meson-group   ',a4,'  does not exi
     1st in this program.'/' execution terminated.'
     2////)
      stop
c
c
c
c
 1995 continue
c
c
c
c
c        add up contributions of mesons
c        ------------------------------
c
 2000 continue
      do 2005 iv=1,6
 2005 v(iv)=vj(iv,iman1)+vj(iv,iman2)
c
c
      if (j.eq.0) then
      v(1)=vj(1,iman1)+vj(1,iman3)
      v(3)=vj(3,iman1)+vj(3,iman3)
      end if
c
c
      if (j.eq.1) then
      v(1)=vj(1,iman1)+vj(1,iman4)
      v(2)=vj(2,iman1)+vj(2,iman4)
      do  iv=3,6
      v(iv)=vj(iv,iman1)+vj(iv,iman3)
      enddo
      end if
c
c
      if (j.eq.2) then
      do 2007 iv=3,6
 2007 v(iv)=vj(iv,iman1)+vj(iv,iman4)
      end if
c
c
      if (mod(j,2).eq.1) go to 2020
c
c
c        j even
c        ------
c
      v(1)=v(1)+vj(1,iman5)+vj(1,iman6)
      v(2)=v(2)+vj(2,iman9)+vj(2,iman10)
c
      if (j.eq.2) then
c
c        the pions for 3P2-3F2
      do 2014 iv=3,6
 2014 v(iv)=v(iv)+vj(iv,iman7)+vj(iv,iman8)
      else
c
c        the pions in all other T=1 coupled cases
      do 2015 iv=3,6
 2015 v(iv)=v(iv)+vj(iv,iman5)+vj(iv,iman6)
      end if
      go to 2030
c
c
c        j odd
c        -----
c
 2020 v(1)=v(1)+vj(1,iman9)+vj(1,iman10)
      v(2)=v(2)+vj(2,iman5)+vj(2,iman6)
c
      do 2025 iv=3,6
 2025 v(iv)=v(iv)+vj(iv,iman9)+vj(iv,iman10)
c
c
c        for all j
c        _________
c
 2030 v(1)=v(1)+vj(1,iman11)+vj(1,iman12)
      v(2)=v(2)+vj(2,iman13)+vj(2,iman14)
      v(3)=v(3)+vj(3,iman15)+vj(3,iman16)
      v(4)=v(4)+vj(4,iman17)+vj(4,iman18)
      do 2035 iv=5,6
 2035 v(iv)=v(iv)+vj(iv,iman17)+vj(iv,iman18)
c
c
c
c
      if (j.eq.0.or..not.heform) go to 4000
c
c
c         translation into (combinations of) helicity states
c
c
      do 3005 i=1,4
 3005 vl(i)=v(i+2)
c
      do 3020 ii=1,4
      iii=ii+2
      v(iii)=0.d0
c
      do 3015 i=1,4
 3015 v(iii)=v(iii)+adminv(ii,i)*vl(i)
 3020 v(iii)=v(iii)*a2j1
c
c
c
c
 4000 return
      end
      subroutine obpard
c
c        obpard provides the parameters for the
c        charge-dependent Bonn potential.
c
c
      implicit real*8 (a-h,o-z)
c
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
      character*4 label
c
c
c        common block for all ob-subroutines
c
      common /cobd/   vj(32,50),c(20,50),fff,ff,f(52),aa(96),ai(19,15),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                ez1,ez2,ct(96),wt(96),
     3                ic(20,50),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(12,3),mggo(12,3),ima(15,12,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(10,50),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c
      common /cnn/ inn
c
c
c        further specifications
c
      dimension cc(5)
      integer name(4),nname(15)
      dimension wscale(3)
      integer imga(3)
      character*4 nucnuc(3)
      data nucnuc/'BNpp','BNnp','BNnn'/
      integer cut,end
      data cut/'cut '/,end/'end '/
      integer two
      data two/'2   '/
      integer mesong(12)
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     1            '1-  ','1-t ','1-tt','1-st','1-ss',
     2            '1+  ','2+  '/
      logical index
      data index/.false./
      data innn/-1/
c
c
c
c
c        parameter tables
c        ----------------
c        ----------------
c
c
c        identification labels
c        ---------------------
      integer ntab1(4,10)
      data ntab1/
     1 '1-t ',' ','rho','    ',
     2 'cut ',' ','   ','    ',
     3 '1-t ',' ','rho','    ',
     4 'cut ',' ','   ','    ',
     5 '1-t ',' ','ome','ga  ',
     6 'cut ',' ','   ','    ',
     7 '1-t ',' ','ome','ga  ',
     8 '0-  ','2','pio','ns  ',
     9 '0-  ','2','pio','ns  ',
     * '0-  ','2','pio','ns  '/
c
      integer ntab2(4)
      data ntab2/
     1 '0+  ','2','sig','mas '/
c
c
c        global parameters
c        -----------------
c        -----------------
      dimension tab1(5,10,3)
      data tab1/
c
c proton-proton potential
c -----------------------
     1        0.84    , 6.1       , 769.9     , 1.        , 0.,
     2        2.      , 0.        , 2.        , 1310.     , 0.1,
     3       20.0     , 0.        , 781.94    , 0.        , 0.,
     4        2.      , 0.        , 2.        , 1500.     , 0.,
     5       20.0     , 0.        , 781.94    , 0.        , 0.,
     6        2.      , 0.        , 2.        , 2000.     , 0.,
     7       20.0     , 0.        , 781.94    , 0.        , 0.,
c t=1:
     8      -13.6     , 134.9764  , 27.2      , 139.56995 , 1720.,
     9      -13.6     , 134.9764  , 27.2      , 139.56995 , 3000.,
c t=0:
     *      -13.6     , 134.9764  , -27.2     , 139.56995 , 1720.,
c
c neutron-proton potential
c ------------------------
     1        0.84    , 6.1       , 769.9     , 1.        , 0.,
     2        2.      , 0.        , 2.        , 1310.     , 0.1,
     3       20.0     , 0.        , 781.94    , 0.        , 0.,
     4        2.      , 0.        , 2.        , 1500.     , 0.,
     5       20.0     , 0.        , 781.94    , 0.        , 0.,
     6        2.      , 0.        , 2.        , 2000.     , 0.,
     7       20.0     , 0.        , 781.94    , 0.        , 0.,
c t=1:
     8      -13.6     , 134.9764  , 27.2      , 139.56995 , 1720.,
     9      -13.6     , 134.9764  , 27.2      , 139.56995 , 3000.,
c t=0:
     *      -13.6     , 134.9764  , -27.2     , 139.56995 , 1720.,
c
c neutron-neutron potential
c -------------------------
     1        0.84    , 6.1       , 769.9     , 1.        , 0.,
     2        2.      , 0.        , 2.        , 1310.     , 0.1,
     3       20.0     , 0.        , 781.94    , 0.        , 0.,
     4        2.      , 0.        , 2.        , 1500.     , 0.,
     5       20.0     , 0.        , 781.94    , 0.        , 0.,
     6        2.      , 0.        , 2.        , 2000.     , 0.,
     7       20.0     , 0.        , 781.94    , 0.        , 0.,
c t=1:
     8      -13.6     , 134.9764  , 27.2      , 139.56995 , 1720.,
     9      -13.6     , 134.9764  , 27.2      , 139.56995 , 3000.,
c t=0:
     *      -13.6     , 134.9764  , -27.2     , 139.56995 , 1720./
c
c
c         partial-wave dependent parameters
c         ---------------------------------
c         ---------------------------------
      dimension tab2(5,4,7,3)
      data tab2/
c
c proton-proton potential
c -----------------------
c j=0:
     1        8.5067,   650.,        0.0 ,      1225.,      1500.,
     2        0.0   ,   350.,        0.0 ,      1225.,      1500.,
     3        3.5   ,   500.,        0.0 ,      1225.,      1500.,
     4        0.0   ,   350.,        0.0 ,      1225.,      1500.,
c j=1:
     1        0.74  ,   350.,       40.  ,      1225.,      2500.,
     2        0.58  ,   350.,        0.0 ,      1225.,      1500.,
     3        0.4   ,   350.,        0.0 ,      1225.,      1500.,
     4        14.08875, 700.,        0.0 ,      1225.,      1500.,
c j=2:
     1        0.78  ,   350.,       75.  ,      1225.,      2500.,
     2        0.4   ,   350.,       130. ,      1225.,      2500.,
     3        0.    ,   350.,       74.  ,       793.,      2500.,
     4        1.2   ,   400.,       11.08,       793.,      2500.,
c j=3:
     1        0.60  ,   350.,       16.0 ,       793.,      2500.,
     2        31.0  ,   700.,        0.0 ,       793.,      1500.,
     3        0.0   ,   350.,       37.0 ,       793.,      2500.,
     4        0.9   ,   350.,       56.0 ,      1225.,      2500.,
c j=4:
     1        0.40  ,   350.,       40.0 ,       793.,      2500.,
     2        5.0   ,   540.,        0.0 ,       793.,      1500.,
     3        0.35  ,   350.,        0.0 ,       793.,      1500.,
     4        0.70  ,   350.,       16.1 ,       793.,      2500.,
c j=5:
     1        0.65  ,   350.,       0.0  ,       793.,      1500.,
     2        12.   ,   540.,       0.0  ,       793.,      1500.,
     3        0.1   ,   350.,       0.0  ,       793.,      1500.,
     4        0.83  ,   350.,       0.0  ,       793.,      1500.,
c j=6:
     1        0.5   ,   350.,       0.0  ,       793.,      1500.,
     2        0.5   ,   350.,       0.0  ,       793.,      1500.,
     3        0.5   ,   350.,       0.0  ,       793.,      1500.,
     4        0.5   ,   350.,       0.0  ,       793.,      1500.,
c
c neutron-proton potential
c ------------------------
c j=0:
     1        8.5067,   650.,        0.0 ,      1225.,      1500.,
     2        0.0   ,   350.,        0.0 ,      1225.,      1500.,
     3        3.5   ,   500.,        0.0 ,      1225.,      1500.,
     4        0.0   ,   350.,        0.0 ,      1225.,      1500.,
c j=1:
     1        0.74  ,   350.,       40.  ,      1225.,      2500.,
     2        0.58  ,   350.,        0.0 ,      1225.,      1500.,
     3        0.4   ,   350.,        0.0 ,      1225.,      1500.,
     4        14.08875, 700.,        0.0 ,      1225.,      1500.,
c j=2:
     1        0.78  ,   350.,       75.  ,      1225.,      2500.,
     2        0.4   ,   350.,       130. ,      1225.,      2500.,
     3        0.    ,   350.,       74.  ,       793.,      2500.,
     4        1.2   ,   400.,       11.08,       793.,      2500.,
c j=3:
     1        0.60  ,   350.,       16.0 ,       793.,      2500.,
     2        31.0  ,   700.,        0.0 ,       793.,      1500.,
     3        0.0   ,   350.,       37.0 ,       793.,      2500.,
     4        0.9   ,   350.,       56.0 ,      1225.,      2500.,
c j=4:
     1        0.40  ,   350.,       40.0 ,       793.,      2500.,
     2        5.0   ,   540.,        0.0 ,       793.,      1500.,
     3        0.35  ,   350.,        0.0 ,       793.,      1500.,
     4        0.70  ,   350.,       16.1 ,       793.,      2500.,
c j=5:
     1        0.65  ,   350.,       0.0  ,       793.,      1500.,
     2        12.   ,   540.,       0.0  ,       793.,      1500.,
     3        0.1   ,   350.,       0.0  ,       793.,      1500.,
     4        0.83  ,   350.,       0.0  ,       793.,      1500.,
c j=6:
     1        0.5   ,   350.,       0.0  ,       793.,      1500.,
     2        0.5   ,   350.,       0.0  ,       793.,      1500.,
     3        0.5   ,   350.,       0.0  ,       793.,      1500.,
     4        0.5   ,   350.,       0.0  ,       793.,      1500.,
c
c neutron-neutron potential
c -------------------------
c j=0:
     1        8.5067,   650.,        0.0 ,      1225.,      1500.,
     2        0.0   ,   350.,        0.0 ,      1225.,      1500.,
     3        3.5   ,   500.,        0.0 ,      1225.,      1500.,
     4        0.0   ,   350.,        0.0 ,      1225.,      1500.,
c j=1:
     1        0.74  ,   350.,       40.  ,      1225.,      2500.,
     2        0.58  ,   350.,        0.0 ,      1225.,      1500.,
     3        0.4   ,   350.,        0.0 ,      1225.,      1500.,
     4        14.08875, 700.,        0.0 ,      1225.,      1500.,
c j=2:
     1        0.78  ,   350.,       75.  ,      1225.,      2500.,
     2        0.4   ,   350.,       130. ,      1225.,      2500.,
     3        0.    ,   350.,       74.  ,       793.,      2500.,
     4        1.2   ,   400.,       11.08,       793.,      2500.,
c j=3:
     1        0.60  ,   350.,       16.0 ,       793.,      2500.,
     2        31.0  ,   700.,        0.0 ,       793.,      1500.,
     3        0.0   ,   350.,       37.0 ,       793.,      2500.,
     4        0.9   ,   350.,       56.0 ,      1225.,      2500.,
c j=4:
     1        0.40  ,   350.,       40.0 ,       793.,      2500.,
     2        5.0   ,   540.,        0.0 ,       793.,      1500.,
     3        0.35  ,   350.,        0.0 ,       793.,      1500.,
     4        0.70  ,   350.,       16.1 ,       793.,      2500.,
c j=5:
     1        0.65  ,   350.,       0.0  ,       793.,      1500.,
     2        12.   ,   540.,       0.0  ,       793.,      1500.,
     3        0.1   ,   350.,       0.0  ,       793.,      1500.,
     4        0.83  ,   350.,       0.0  ,       793.,      1500.,
c j=6:
     1        0.5   ,   350.,       0.0  ,       793.,      1500.,
     2        0.5   ,   350.,       0.0  ,       793.,      1500.,
     3        0.5   ,   350.,       0.0  ,       793.,      1500.,
     4        0.5   ,   350.,       0.0  ,       793.,      1500./
c
c        this has been the end of all tables
c        -----------------------------------
c
c
c
c
10002 format (' type/meson       m e s o n    p a r a m e t e r s')
10004 format (1h ,a4,a1,a3,a4,4f12.5,f12.2)
10008 format (1h ,61(1h-))
10011 format (///' CDBONN: The Charge-Dependent Bonn NN Potential ("CDBo
     1nn2000")')
10017 format (///)
10018 format (' Potential Type:  ',a4)
c
c
c
c
      if (index) go to 50
      index=.true.
c
      x=-1.d0
      y=-1.d0
c
c
c
c
c        maxima of certain indices related to the dimension as follows:
c        dimension c(mme,imee),ic(mice,imee),indc(mindce,imee),
c                  mgg(mge,3),mggo(mge,3),mesong(mge),vj(32,imee),
c                  ima(mee,mge,3)
c
      mge=12
      mee=15
      mme=20
      mice=20
      mindce=10
      imee=50
c        mme always ge mice, mindce
      imb=1
      endep=.false.
c
c
c        set all meson-parameters and indices to zero or .false.
c
      do 1 int=1,3
      imga(int)=0
      indpar(int)=.false.
      do 1 mgx=1,mge
      mgg(mgx,int)=0
    1 mggo(mgx,int)=0
c
c
      do 2 il=1,imee
      do 2 mm=1,mme
      if (mm.le.mindce) indc(mm,il)=.false.
      if (mm.le.mice) ic(mm,il)=0
    2 c(mm,il)=0.d0
c
c
c        write headline
c
c**** write (kwrite,10011)
c**** write (kwrite,10008)
c**** write (kwrite,10008)
c**** write (kwrite,10017)
c
c
      ift(inter)=1
c
c        scaling mass
c
      wscale(inter)=938.27231d0
c
c
c
c
   50 if (inn.eq.innn) go to 55
      innn=inn
c
c
      imga(inter)=0
      do 11 mgx=1,mge
      mgg(mgx,inter)=0
   11 mggo(mgx,inter)=0
c
      ime=0
      line=0
      jj=-1
c
c
c
c
c**** write (kwrite,10018) nucnuc(inn)
      label=nucnuc(inn)
c
c
      if (inn.eq.1) then
      wn1=938.27231d0
      wn2=938.27231d0
      else
      if (inn.eq.2) then
      wn1=939.56563d0
      wn2=938.27231d0
      else
      wn1=939.56563d0
      wn2=939.56563d0
      end if
      end if
c
c
      wnn(inter)=dsqrt(wn1*wn2)
c
c
      wn=wnn(inter)
      wnq=wn*wn
      dwn=1.d0/wn
      dwnq=dwn*dwn
c
c
c        write headline for meson parameters
c
c**** write (kwrite,10008)
c**** write (kwrite,10008)
c**** write (kwrite,10002)
c**** write (kwrite,10008)
c**** write (kwrite,10008)
      go to 61
c
c
c
c
   55 ime=10
      mgx=mggo(imga(inter),inter)
      mgg(mgx,inter)=0
      mggo(imga(inter),inter)=0
      imga(inter)=imga(inter)-1
c
c        write headline for meson parameters
c
c**** write (kwrite,10008)
c**** write (kwrite,10008)
c**** write (kwrite,10002)
c**** write (kwrite,10008)
c
c
c
c
c        read, write, and store meson parameters
c        ---------------------------------------
c        ---------------------------------------
c
c
c
   61 if (ime.eq.18) go to 2000
c
c
c        instead of reading, get meson parameters from data tables
c        ---------------------------------------------------------
c
c
      if (ime.lt.10) then
      line=line+1
      do 63 i=1,5
      if (i.le.4) then
      name(i)=ntab1(i,line)
      end if
   63 cc(i)=tab1(i,line,inn)
c
c
      else
      if (j.eq.jj) then
      line=line+1
      else
      line=1
      jj=j
      j1=j+1
      if (j.gt.6) j1=7
      end if
c
      do 65 i=1,5
      if (i.le.4) then
      name(i)=ntab2(i)
      end if
   65 cc(i)=tab2(i,line,j1,inn)
      end if
c
c
c        check if record just read contains cut-off parameters
c
      if (name(1).eq.cut) go to 80
c
c
c
c
c        write meson-parameters
c        ----------------------
c
c
c
c
c**** write (kwrite,10004) name,cc
c
c        find out number of meson-group mg
c
      do 73 mg=1,mge
      if (name(1).eq.mesong(mg)) go to 74
   73 continue
      go to 9000
c
c
   74 if (name(2).eq.two) go to 1000
c
c
c
c
c        store meson parameters, which are no cut-off parameters
c        -------------------------------------------------------
c
c
c
c
      ime=ime+1
      if (ime.gt.imee) go to 9011
      mgg(mg,inter)=mgg(mg,inter)+1
      m=mgg(mg,inter)
      if (m.gt.mee) go to 9001
      ima(m,mg,inter)=ime
      if (m.ne.1) go to 76
      imga(inter)=imga(inter)+1
      mggo(imga(inter),inter)=mg
   76 continue
c
c        store coupling constant g**2/4pi
      c(1,ime)=cc(1)
c        store coupling constant f*g/4pi
      c(3,ime)=cc(1)*cc(2)*wn/wscale(inter)
c        store coupling constant f**2/4pi
      c(2,ime)=cc(2)*c(3,ime)*wn/wscale(inter)
c        store meson mass squared in units of nucleon mass squared
      c(4,ime)=cc(3)*cc(3)*dwnq
c
c        get iso-spin
      icc=cc(4)
      if (icc.ne.0.and.icc.ne.1) go to 9004
c         store isospin as logical constant
      if (icc.eq.1) indc(1,ime)=.true.
c        store parameter for meson propagator (iprop)
      ic(1,ime)=cc(5)
      if (ic(1,ime).ne.0) go to 9005
c
c        index values for further storing
      mi=4
      mm=5
      go to 61
c
c
c
c
c        write cut-off parameters
c        ------------------------
c
c
c
c
   80 continue
c**** write (kwrite,10004) name,cc
c
c
      if (ime.eq.1) eps=cc(5)
c
c
c
c        if cutoff type = 0, ignore cutoff
      if (cc(1).eq.0.d0) go to 61
c
c
c
c
c        store cut-off parameters
c        ------------------------
c
c
c
c
c        store type of cut-off
      ic(mi,ime)=cc(1)
      ityp=ic(mi,ime)
      if (ityp.ne.2) go to 9002
c        store and test type of denominator of cut-off
      ic(mi+1,ime)=cc(2)
      if (ic(mi+1,ime).ne.0) go to 9006
c
c
c        cut-off of monopole/dipole type
c        *******************************
c
c
c        store and test exponent of cut-off
      ic(mi+2,ime)=cc(3)
      if (ic(mi+2,ime).lt.0) go to 9009
      if (ic(mi+2,ime).gt.0) go to 101
c        exponent is zero, omit cut-off
      ic(mi,ime)=0
      ic(mi+1,ime)=0
      go to 999
  101 if (ic(mi+2,ime).ne.2) go to 9012
c        store first cut-off mass
      c(mm,ime)=(cc(4)+eps)**2*dwnq
c        store second cut-off mass
      c(mm+1,ime)=(cc(4)-eps)**2*dwnq
      mi=mi+3
      mm=mm+2
c
c
c
c
c        end cut-offs
c        ************
c
c        test dimensions
  999 if (mi.gt.mice.or.mm.gt.mme) go to 9010
c
c
      go to 61
c
c
c
c
c        two mesons on one input line
c        ----------------------------
c
c
c        store input parameters and set defaults
c
c
 1000 do 1995 ii=1,2
      ime=ime+1
      if (ime.gt.imee) go to 9011
      mgg(mg,inter)=mgg(mg,inter)+1
      m=mgg(mg,inter)
      if (m.gt.mee) go to 9001
      ima(m,mg,inter)=ime
      if (m.ne.1) go to 1076
      imga(inter)=imga(inter)+1
      mggo(imga(inter),inter)=mg
 1076 continue
c
c        store coupling constant g**2/4pi
      if (ii.eq.1) then
      c(1,ime)=cc(1)
      else
      c(1,ime)=cc(3)
      end if
c
c        scale the pi-NN coupling constant
      if (ime.ge.5.and.ime.le.10) then
      c(1,ime)=c(1,ime)*(wn/wscale(inter))**2
      end if
c
c        set coupling constant f*g/4pi
      c(3,ime)=0.d0
c        set coupling constant f**2/4pi
      c(2,ime)=0.d0
c        store meson mass squared in units of nucleon mass squared
      if (ii.eq.1) then
      c(4,ime)=cc(2)*cc(2)*dwnq
      else
      c(4,ime)=cc(4)*cc(4)*dwnq
      end if
c
c         set isospin-0 as logical constant
      indc(1,ime)=.false.
c        set parameter for meson propagator (iprop=0)
      ic(1,ime)=0
c
c        index values for further storing
      mi=4
      mm=5
c
c
c        store and set cut-off parameters
c
c        set type of cut-off
      ic(mi,ime)=2
c        set type of denominator of cut-off
      ic(mi+1,ime)=0
c
c
c        cut-off of monopole/dipole type
c        *******************************
c
c
c        set exponent of cut-off
      ic(mi+2,ime)=2
c        store first cut-off mass
      c(mm,ime)=(cc(5)+eps)**2*dwnq
c        store second cut-off mass
      c(mm+1,ime)=(cc(5)-eps)**2*dwnq
      mi=mi+3
      mm=mm+2
c
c
c
c
c        end cut-offs
c        ************
c
c        test dimensions
      if (mi.gt.mice.or.mm.gt.mme) go to 9010
c
 1995 continue
c
      go to 61
c
c
c
c
c        end of mesons for one j
c        -----------------------
c
c
 2000 imaa(inter)=imb
      imea(inter)=ime
c**** write (kwrite,10008)
c**** write (kwrite,10008)
c
      return
c
c
c
c
c        errors
c        ------
c        ------
c
c
c
c
 9000 write (kwrite,19000) name(1)
19000 format (/////' error in cdbonn:  meson-group   ',a4,'   does not
     1exist in this program.'/' execution terminated.'////)
      go to 9999
c
c
 9001 write (kwrite,19001)
19001 format (/////' error in cdbonn:   too many mesons within a meson-g
     1roup with respect to'/' the given dimensions.    execution termina
     2ted.'////)
      go to 9999
c
c
 9002 write (kwrite,19002) cc(1)
19002 format (/////' error in cdbonn:  cut-off type',f10.4,'  does not e
     1xist in this program.'/' execution terminated.'////)
      go to 9999
c
c
 9004 write (kwrite,19004) cc(4)
19004 format (/////' error in cdbonn:  isospin has the non-permissible
     1value',f10.4,'  .'/' execution terminated.'////)
      go to 9999
c
c
 9005 write (kwrite,19005) cc(5)
19005 format (/////' error in cdbonn:    iprop has the non-permissible
     1value',f10.4,'  .'/' execution terminated.'////)
      go to 9999
c
c
 9006 write (kwrite,19006) cc(2)
19006 format (/////' error in cdbonn: the index for the denominator of
     1 the cut-off has the'/' non-permissible value',f10.4,' . execution
     2 terminated.'////)
      go to 9999
c
c
 9009 write (kwrite,19009)
19009 format (/////' error in cdbonn:   the exponent of the cut-off is
     1less than zero.'/' execution terminated.'////)
      go to 9999
c
c
 9010 write (kwrite,19010)
19010 format (/////' error in cdbonn: too many cut-off parameters with
     1 respect to the given'/' dimensions. execution terminated.'////)
      go to 9999
c
c
 9011 write (kwrite,19011)
19011 format (/////' error in cdbonn:  too many mesons with respect to
     1 the dimensions given'/' in this program. execution terminated.'
     2////)
      go to 9999
c
c
 9012 write (kwrite,19012)
19012 format (/////' error in cdbonn:   the exponent of the cut-off is
     1not two.'/' execution terminated.'////)
      go to 9999
c
c
 9999 stop
      end
      subroutine obstrd (icase,max,mex)
c
c        obstrd computes the structure of one-boson-exchanges
c
c
      implicit real*8 (a-h,o-z)
c
c
c        common blocks
c
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
c
c
c        common block for all ob-subroutines
c
      common /cobd/   vj(32,50),c(20,50),fff,ff,f(52),aa(96),ai(19,15),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                ez1,ez2,ct(96),wt(96),
     3                ic(20,50),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(12,3),mggo(12,3),ima(15,12,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(10,50),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c
      common /cnn/ inn
c
c
c     further specifications
c
      dimension vv(32)
      dimension tt(2,3)
      data jj/-1/
      logical index
      data index/.false./
      logical indiso
c
c
c
c
      if (index) go to 50
      index=.true.
c
c
      tt(1,1)=1.d0
      tt(2,1)=-3.d0
c
      do 1 ii=2,3
      do 1 i=1,2
    1 tt(i,ii)=1.d0
c
c
c
c
c
   50 do 1095 m=max,mex
      im=ima(m,mg,inter)
c
      if (mg.le.5.and.c(mc,im).eq.0.d0) go to 1095
c
c
      if (mc.ne.1) go to 60
c
c
c
c
c        call integrals
c        --------------
c
c
c
c
      call obaid
c
c
c
c
   60 continue
c
      if (c(mc,im).eq.0.d0) go to 1095
c
c
c
c
c        nn-nn helicity amplitudes
c        -------------------------
c
c
c        vv(1), ..., vv(6) contain in the following order:
c        0v, 1v, 12v, 34v, 55v, 66v.
c
c
c        basic structure
c
c
  100 ive=6
c
      vv(1)=f(1)*ai(1,m)+f(2)*ai(2,m)
      vv(2)=f(3)*ai(1,m)+f(4)*ai(3,m)
      vv(3)=f(5)*ai(1,m)+f(6)*ai(2,m)
      vv(4)=f(4)*ai(1,m)+f(3)*ai(3,m)
      vv(5)=f(7)*ai(4,m)
      vv(6)=f(8)*ai(4,m)
c
c
      go to (1000,120),icase
c
c
c        additional terms for the case of tensor-tensor coupling
c
c
  120 vv(1)=vv(1)+f(9)*ai(5,m)
      vv(2)=vv(2)+f(10)*ai(2,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(10)*ai(5,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(10)*ai(6,m)
         e1=f(11)*ai(7,m)
      vv(5)=vv(5)+e1
      vv(6)=vv(6)+e1
      go to 1000
c
c
c
c
c        set certain cases to zero
c
 1000 if (j.ne.0) go to 1021
      vv(2)=0.d0
      vv(4)=0.d0
      vv(5)=0.d0
      vv(6)=0.d0
c
 1021 mmod=mod(j,2)
      if (.not.sing.or.(mmod.eq.1.and.inn.ne.2)) vv(1)=0.d0
      if (.not.trip.or.(mmod.eq.0.and.inn.ne.2)) vv(2)=0.d0
      if (coup.and.(mmod.eq.0.or.inn.eq.2)) go to 1030
      do 1025 iv=3,6
 1025 vv(iv)=0.d0
c
 1030 continue
c
c
c        transformation into lsj-formalism
c
      if (j.eq.jj) go to 1035
      jj=j
      aj=dfloat(j)
      aj1=dfloat(j+1)
      d2j1=1.d0/dfloat(2*j+1)
      arjj1=dsqrt(aj*aj1)
c
 1035 v3=vv(3)
      v4=vv(4)
      v5=vv(5)
      v6=vv(6)
      v34=arjj1*(v3-v4)
      v56=arjj1*(v5+v6)
      vv(3)=d2j1*(aj1*v3+aj*v4-v56)
      vv(4)=d2j1*(aj*v3+aj1*v4+v56)
      vv(5)=d2j1*(v34+aj1*v5-aj*v6)
      vv(6)=d2j1*(v34-aj*v5+aj1*v6)
c
c        after transformation into lsj formalism,
c        vv(3), ..., vv(6) contain:
c        v++, v--, v+-, v-+.
c
c
c
c
c        multiply with factors
c        ---------------------
c
c
c
c
      is=mod(j,2)+1
      it=mod(is,2)+1
      indiso=indc(1,im)
c        get coupling constant
      cmc=c(mc,im)
      fc=fff*ff*cmc
      do 1045 iv=1,ive
c
c        multiply with coupling-constant and factors fff and ff
c
      vv(iv)=vv(iv)*fc
c
c        multiply with isospin factor
c
      if (.not.indiso) go to 1045
      if (iv.eq.2) go to 1043
      vv(iv)=vv(iv)*tt(is,inter)
      go to 1045
 1043 vv(iv)=vv(iv)*tt(it,inter)
c
c
c        add up in case of several couplings for one meson and store
 1045 vj(iv,im)=vj(iv,im)+vv(iv)
c
c
 1095 continue
c
c
      return
      end
      subroutine obaid
c
c        obaid performs the integration over angle theta
c        (necessary for the partial wave decomposition)
c        in analytic form by using the Legendre functions of the
c        second kind.
c
c
      implicit real*8 (a-h,o-z)
c
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
c
c
c        common block for all ob-subroutines
c
      common /cobd/   vj(32,50),c(20,50),fff,ff,f(52),aa(96),ai(19,15),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                ez1,ez2,ct(96),wt(96),
     3                ic(20,50),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(12,3),mggo(12,3),ima(15,12,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(10,50),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c
      dimension gi(5,7)
      data jj/-1/
      logical index
      data index/.false./
c
c
c
c
      if (index) go to 50
      index=.true.
c
      sqr2=dsqrt(2.d0)
c
c
c
c
c
   50 if (j.eq.jj) go to 70
      jj=j
c
c
      aj=dfloat(j)
      aj1=dfloat(j+1)
      dj1=1.d0/aj1
      ajdj1=aj*dj1
      aaj=dsqrt(ajdj1)
c
c
      if (j.eq.0) then
      delj0=1.d0
      else
      delj0=0.d0
      end if
c
      if (j.eq.1) then
      delj1=1.d0
      else
      delj1=0.d0
      end if
c
c
   70 continue
c
c
c
c
      mi=4
      mm=3
      ityp=ic(mi,im)
      nexp=ic(mi+2,im)
c
      nterms=nexp+1
      if (ityp.eq.0) nterms=1
c
c
      do 555 i=1,nterms
      mmi=mm+i
c
c        calculate the argument for the legendre function
c
      if (x.eq.y) then
      zstamm=1.d0
      zdelta=c(mmi,im)/xy2
      else
      zstamm=(xxpyy+c(mmi,im))/xy2
      zdelta=0.d0
      end if
c
      z=zstamm+zdelta
c
c
c        call legendre functions of the second kind
c
      if (j.eq.0) then
c
      call legen2 (qj,qjp1,zzq1m,1,zstamm,zdelta)
      qjm1=0.d0
c
      else
c
      call legen2 (qjm1,qj,zzq1m,j,zstamm,zdelta)
c
      end if
c
c
      gi(i,1)=qj
c
      if (j.eq.0) then
      gi(i,2)=qjp1
      else
      gi(i,2)=z*qj-delj0
      end if
c
      gi(i,3)=ajdj1*z*qj+dj1*qjm1
      gi(i,4)=aaj*(z*qj-qjm1)
c
      if (j.eq.1) then
      gi(i,5)=zzq1m
      gi(i,6)=0.5d0*(zzq1m+qj)
      gi(i,7)=0.5d0*sqr2*(zzq1m-qj)
      else
      gi(i,5)=z*gi(i,2)-delj1/3.d0
      gi(i,6)=z*gi(i,3)-2.d0*delj1/3.d0
      gi(i,7)=z*gi(i,4)+sqr2*delj1/3.d0
      end if
c
c
      if (i.eq.1) then
      fact=1.d0
      else
      ix=1
      if (i.eq.3) ix=-1
      fact=(c(mmi+ix,im)-c(4,im))/(c(mmi,im)-c(mmi+ix,im))
      end if
c
c
      do 545 ii=1,7
      gi(i,ii)=fact*gi(i,ii)
  545 continue
c
c
  555 continue
c
c
      do 725 ii=1,7
      ai(ii,m)=0.d0
      do 715 i=1,nterms
  715 ai(ii,m)=ai(ii,m)+gi(i,ii)
  725 continue
c
c
      dxy=2.d0/xy2
      do 2015 ii=1,7
 2015 ai(ii,m)=ai(ii,m)*dxy
c
c
c
c
      return
      end
      subroutine legen2 (qjm1,qj,zzq1m,j,zstamm,zdelta)
c
c**** legendre-funktionen der zweiten art ****
c
c**** notation: qjm1 = Q_{J-1}, qj = Q_J, z=zstamm+zdelta,
c**** in the case of j=1, this program also provides
c****                     zzq1m =  z*z*Q_1-1./3.
c
c
c     author:
c              R. Machleidt
c              Institut fuer Theoretische Kernphysik
c              der Universitaet Bonn
c              Nussallee 14-16
c              D-5300 Bonn
c              W. Germany
c
c              original version: April 1972
c              last revision: April 1995
c
c
c**** genauigkeit:
c**** j kleiner gleich 10   15 stellen
c**** j gleich  11 bis 30   mindestens 13 stellen
c**** j gleich 31 bis 100   mindestens 12 stellen
c**** eine dimension der koeffizienten von  me=40000  c(2,40001)
c**** ist gut bis j=50;
c**** fuer j=100 wird  me=150000  c(2,150001) benoetigt.
c
c
      implicit real*8 (a-h,o-z)
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
      data tolr/1.d-16/
      dimension c(2,40001)
      data me/40000/
      data jj/-1/
c
c
c**** berechnung des arguments ****
      z=zstamm+zdelta
c
      qjm1=0.d0
      qj=0.d0
      zzq1m=0.d0
      if (j.lt.0) go to 123
      if (z.le.1.d0) go to 113
c
c
c**** fallunterscheidung ****
      if (j.ne.0) go to 2
      if (z-10.d0) 10,10,11
    2 if (j.ne.1) go to 3
      if (z-1.5d0) 10,10,11
    3 if (j.ne.2) go to 4
      if (z-1.2d0) 10,10,11
    4 zcut=1.d0+dfloat(j)**(-2.d0)
      if (z-zcut) 10,10,11
c
c**** rekursive berechnung mit dem logarithmus ****
   10 zdel=zstamm-1.d0
      zdel  =zdel+zdelta
      zz=2.d0/zdel  +1.d0
      qjm1=0.5d0*dlog(zz)
      if (j.eq.0) then
      qj=qjm1
      qjm1=0.d0
      return
      end if
      qj=z*qjm1-1.d0
      if (j.eq.1) then
      zzq1m=z*z*qj-1.d0/3.d0
      return
      end if
      do  7 i=2,j
      qq=dfloat(2*i-1)/dfloat(i)*z*qj-dfloat(i-1)/dfloat(i)*qjm1
      qjm1=qj
    7 qj=qq
      return
c
c**** berechnung mit reihe ****
c**** der laufende index m ist immer mue plus eins ****
   11 zqinv=z**(-2)
      zzqinv=1.d0
      qjm1=0.d0
      qj=0.d0
      zzq1m=0.d0
      if (j.eq.jj) go to 12
      jj=j
      ma=1
      go to 14
   12 do 13 m=1,mme
      cz1=c(1,m)*zzqinv
      cz= c(2,m)*zzqinv
      qjm1=cz1+qjm1
      qj=  cz+qj
      if (j.eq.1) then
      if (m.eq.1) then
      zzq1m=0.d0
      else
      zzq1m=zzq1m+cz
      end if
      if (cz.lt.tolr*zzq1m) go to 62
      go to 13
      end if
      if (cz.lt.tolr*qj) go to 62
   13 zzqinv=zzqinv*zqinv
      ma=mme+1
c
c**** verteiler ****
   14 if (j.le.1) go to 20
      if (j.eq.2) go to 30
      if (mod(j,2)) 50,40,50
c
c**** die faelle j gleich null und j gleich eins ****
   20 if (ma.ne.1) go to 22
      ma=2
      c(1,1)=1.d0
      c(2,1)=1.d0/3.d0
      qjm1=c(1,1)
      qj=c(2,1)
      zzq1m=0.d0
      zzqinv=zqinv
   22 do 21 m=ma,me
      c(1,m)=c(2,m-1)
      c(2,m)=1.d0/dfloat(2*m+1)
      cz1=c(1,m)*zzqinv
      cz= c(2,m)*zzqinv
      qjm1=cz1+qjm1
      qj=  cz+qj
      if (j.eq.1) then
      zzq1m=zzq1m+cz
      if (cz.lt.tolr*zzq1m) go to 61
      go to 21
      end if
      if (cz.lt.tolr*qj) go to 61
   21 zzqinv=zzqinv*zqinv
      go to 60
c
c**** fall j gleich zwei ****
   30 do 31 m=ma,me
      m2=2*m
      c(1,m)=1.d0/dfloat(m2+1)
      c(2,m)=c(1,m)*dfloat(m2)/dfloat(m2+3)
      qjm1=   c(1,m)*zzqinv+qjm1
      cz=     c(2,m)*zzqinv
      qj=cz+qj
      if (cz.lt.tolr*qj) go to 61
   31 zzqinv=zzqinv*zqinv
      go to 60
c
c**** fall j ist gerade ****
   40 do 41 m=ma,me
      m2=2*m
c**** zaehler ****
      aehler=1.d0
      ka=m2
      kez=m2+j-4
      do 42 k=ka,kez,2
   42 aehler=aehler*dfloat(k)
c**** nenner ****
      aenner=1.d0
      ka=m2+j-1
      ken=m2+2*j-3
      do 43 k=ka,ken,2
   43 aenner=aenner*dfloat(k)
      c(1,m)=aehler/aenner
      c(2,m)=c(1,m)*dfloat(kez+2)/dfloat(ken+2)
      qjm1=   c(1,m)*zzqinv+qjm1
      cz=     c(2,m)*zzqinv
      qj=cz+qj
      if (cz.lt.tolr*qj) go to 61
   41 zzqinv=zzqinv*zqinv
      go to 60
c
c**** fall j ist ungerade ****
   50 do 51 m=ma,me
      m2=2*m
c**** zaehler ****
      aehler=1.d0
      ka=m2
      ke=m2+j-3
      do 52 k=ka,ke,2
   52 aehler=aehler*dfloat(k)
      if (m.ne.1) go to 55
      m2=0
      go to 54
   56 m2=2
   55 c(1,m)=aehler/aenner
c**** nenner ****
   54 aenner=1.d0
      ka=m2+j
      ke=m2+2*j-1
      do 53 k=ka,ke,2
   53 aenner=aenner*dfloat(k)
      if (m2) 57,56,57
   57 c(2,m)=aehler/aenner
      qjm1=   c(1,m)*zzqinv+qjm1
      cz=     c(2,m)*zzqinv
      qj=cz+qj
      if (cz.lt.tolr*qj) go to 61
   51 zzqinv=zzqinv*zqinv
c
c
   60 mme=me
      write (kwrite,1131)
 1131 format (/////' warning in legen2. the dimension for the'/
     1' coefficients is too small. the Legendre function of the'/
     2' second kind may be inaccurate.'/////)
      go to 62
c
   61 mme=m
c
c**** schlussrechnung ****
   62 zmj1=z**(-j-1)
      if (j.eq.0) go to 68
      qj=qj*zmj1
      qjm1=qjm1*zmj1*z
      return
   68 qj=qjm1*zmj1
      qjm1=0.d0
      return
c
c**** fehlermeldung ****
  113 write (kwrite,1130)
 1130 format (/////' error in legen2. the argument of the'/
     1' Legendre function of the second kind is smaller or'/
     2' equal one. the function is set to zero.'/
     3' results may be wrong.'/////)
      return
  123 write (kwrite,1230)
 1230 format (/////' error in legen2. the parameter j of the'/
     1' Legendre function of the second kind is smaller zero.'/
     2' the function is set to zero.'/
     3' results may be wrong.'/////)
      return
      end
c***************** this program consistently in double precision ****
c name:    dminv
c        programmbibliothek rhrz bonn        28/11/78       dminv
c                                            fortran iv     ibm 370/168
c
c purpose:
c
c invert a matrix
c
c usage:   call dminv (a,n,d,l,m)
c
c parameters:
c
c a:       input matrix, destroyed in computation and replaced by
c          resultant inverse.
c          double precision required.
c
c n:       order of matrix a
c
c d:       resultant determinant
c          double precision required.
c
c l:       work vector of length n
c
c m:       work vector of length n
c
c remarks: matrix a must be a general matrix
c
c method:
c
c the standard gauss-jordan method is used. the determinant
c is also calculated. a determinant of zero indicates that
c the matrix is singular.
c
c programs required:
c          none
c
c author:  ibm, ssp iii
c
c**********************************************************************
      subroutine dminv (a,n,d,l,m)
      implicit real*8 (a-h,o-z)
      dimension a(1),l(1),m(1)
c
c
c        search for largest element
c
      d=1.d0
      nk=-n
      do 80 k=1,n
      nk=nk+n
      l(k)=k
      m(k)=k
      kk=nk+k
      biga=a(kk)
      do 20 j=k,n
      iz=n*(j-1)
      do 20 i=k,n
      ij=iz+i
   10 if (dabs(biga)-dabs(a(ij)))  15,20,20
   15 biga=a(ij)
      l(k)=i
      m(k)=j
   20 continue
c
c        interchange rows
c
      j=l(k)
      if(j-k) 35,35,25
   25 ki=k-n
      do 30 i=1,n
      ki=ki+n
      hold=-a(ki)
      ji=ki-k+j
      a(ki)=a(ji)
   30 a(ji) =hold
c
c        interchange columns
c
   35 i=m(k)
      if(i-k) 45,45,38
   38 jp=n*(i-1)
      do 40 j=1,n
      jk=nk+j
      ji=jp+j
      hold=-a(jk)
      a(jk)=a(ji)
   40 a(ji) =hold
c
c        divide column by minus pivot (value of pivot element is
c        contained in biga)
c
   45 if(biga) 48,46,48
   46 d=0.d0
      return
   48 do 55 i=1,n
      if(i-k) 50,55,50
   50 ik=nk+i
      a(ik)=a(ik)/(-biga)
   55 continue
c
c        reduce matrix
c
      do 65 i=1,n
      ik=nk+i
      hold=a(ik)
      ij=i-n
      do 65 j=1,n
      ij=ij+n
      if(i-k) 60,65,60
   60 if(j-k) 62,65,62
   62 kj=ij-i+k
      a(ij)=hold*a(kj)+a(ij)
   65 continue
c
c        divide row by pivot
c
      kj=k-n
      do 75 j=1,n
      kj=kj+n
      if(j-k) 70,75,70
   70 a(kj)=a(kj)/biga
   75 continue
c
c        product of pivots
c
      d=d*biga
c
c        replace pivot by reciprocal
c
      a(kk)=1.d0/biga
   80 continue
c
c        final row and column interchange
c
      k=n
  100 k=(k-1)
      if(k) 150,150,105
  105 i=l(k)
      if(i-k) 120,120,108
  108 jq=n*(k-1)
      jr=n*(i-1)
      do 110 j=1,n
      jk=jq+j
      hold=a(jk)
      ji=jr+j
      a(jk)=-a(ji)
  110 a(ji) =hold
  120 j=m(k)
      if(j-k) 100,100,125
  125 ki=k-n
      do 130 i=1,n
      ki=ki+n
      hold=a(ki)
      ji=ki-k+j
      a(ki)=-a(ji)
  130 a(ji) =hold
      go to 100
  150 return
      end
      subroutine obpar
c
c        obpar reads writes and stores the parameter for all ob-subrout.
c
c
      implicit real*8 (a-h,o-z)
c
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
c
c
c        common block for all ob-subroutines
c
      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                ez1,ez2,ct(96),wt(96),
     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c
c        further specifications
c
      dimension cc(5)
      integer name(3),nname(15)
      integer imga(3)
      integer cut,cutg,end
      data cut/'cut '/,cutg/'cutg'/,end/'end '/
      integer mesong(12)
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     1                   '1-  ','1-t ','1-tt','1-st','1-ss',
     2                    '1+  ','2+  '/
      logical index
      data index/.false./
      logical zerocp,indcut
      data zerocp/.true./,indcut/.false./
      data uf/197.327d0/
c
c
c
c
c
10000 format (2a4,a2,15a4)
10001 format (1h1)
10002 format (1h0/' jp  name      g**2      f/g       mass    iso-spin
     1iprop/spe'/9x,'cut typ     c u t - o f f   p a r a m e t e r s')
10003 format (2a4,a2,5f10.4)
10004 format (1h0,2a4,a2,3f10.4,2(f7.1,3x))
10005 format (1h ,2a4,a2,f3.1,f11.1,f9.4,f14.4,f10.4)
10006 format (2a4,a2,3i3)
10007 format (1h0,2a4,a2,3i3)
10008 format (1h ,61(1h-))
10011 format (1h1  //' obnn:  one-boson-exchange nn-nn interaction (nume
     1r. integra.)')
10012 format (1h1  //' obnd:  one-boson-exchange nn-nd interaction (nume
     1r. integra.)')
10013 format (1h1  //' obdd:  one-boson-exchange nn-dd interaction (nume
     1r. integra.)')
10015 format ('0input-parameter-set:'/1h ,20(1h-))
10016 format (1h0,2a4,a2,15a4)
c
c
10020 format (2a4,a2,f10.7,4f10.4)
10021 format (1h0,2a4,a2,f10.7,f10.4,f9.2,1x,2(f7.1,3x))
c
c
      if (index) go to 50
      index=.true.
c
      x=-1.d0
      y=-1.d0
c
c
c
c
c        maxima of certain indices related to the dimension as follows:
c        dimension c(mme,imee),ic(mice,imee),indc(mindce,imee),
c                  mgg(mge,3),mggo(mge,3),mesong(mge),vj(32,imee),
c                  ima(mee,mge,3)
c
      mge=12
      mee=5
      mme=10
      mice=10
      mindce=2
      imb=1
      ime=0
      imee=15
      imec=0
c        mme always ge mice, mindce
c
c        set all meson-parameters and indices to zero or .false.
c
      do 1 int=1,3
      imga(int)=0
      indpar(int)=.false.
      do 1 mgx=1,mge
      mgg(mgx,int)=0
    1 mggo(mgx,int)=0
c
c
      do 2 il=1,imee
      do 2 mm=1,mme
      if (mm.le.mindce) indc(mm,il)=.false.
      if (mm.le.mice) ic(mm,il)=0
    2 c(mm,il)=0.d0
      endep=.false.
c
c
c        indc(2,il) is reserved for information concerning the eikonal
c        form-factor
c
c
c
c
c
c
c        reading and writing of first 4,5 cards
c        --------------------------------------
c        --------------------------------------
c
c
c
c        write headline and read and write name of parameter set
c
   50 go to (51,52,53),inter
   51 write (kwrite,10011)
      go to 55
   52 write (kwrite,10012)
      go to 55
   53 write (kwrite,10013)
   55 write (kwrite,10008)
      write (kwrite,10015)
      read  (kread, 10000) name,nname
      write (kwrite,10016) name,nname
      if (inter.eq.1) label=name(1)
      indpar(inter)=.true.
c
c        read and write index-parameter concerning the factor of the
c        potential
c
      read  (kread, 10006) name,ift(inter)
      write (kwrite,10007) name,ift(inter)
      iftyp=ift(inter)
c**** if (iftyp.lt.0.or.iftyp.gt.4.or.inter.eq.1.and.iftyp.gt.2)goto9003
c
c        read and write parameters for numerical integration
c
      read  (kread, 10006) name,mint(inter),maxt(inter)
      write (kwrite,10007) name,mint(inter),maxt(inter)
c
c        read and write mass of nucleon
c
      read  (kread, 10003) name,wn
      write (kwrite,10004) name,wn
      wnq=wn*wn
      dwn=1.d0/wn
      dwnq=dwn*dwn
      wnn(inter)=wn
      if (inter.lt.2) go to 60
c
c        read and write mass of n-star
c
      read  (kread, 10003) name,wdd(inter)
      write (kwrite,10004) name,wdd(inter)
c
c        write headline for meson parameters
c
   60 write (kwrite,10002)
      write (kwrite,10008)
c
c
c
c
c        read, write and store meson parameters
c        --------------------------------------
c        --------------------------------------
c
c
c
   61 read  (kread, 10003) name,cc
c
c        check if data-card just read contains cut-off parameters
c
      if (name(1).eq.cut.or.name(1).eq.cutg) go to 70
c
c        check if end of mesons
c
      if (name(1).eq.end) go to 2000
c
c
c
c
c        write meson-parameters, which are no cut-off parameters
c        -------------------------------------------------------
c
c
c
c
      indcut=.false.
c
      write (kwrite,10004) name,cc
c
c        check if coupling constants are zero
c
      if (cc(1).ne.0.d0) go to 62
      zerocp=.true.
      go to 61
c
   62 zerocp=.false.
c
c        find out number of meson-group mg
c
      do 63 mg=1,mge
      if (name(1).eq.mesong(mg)) go to 64
   63 continue
      go to 9000
c
c
c
c
c        store meson parameters, which are no cut-off parameters
c        -------------------------------------------------------
c
c
c
c
   64 ime=ime+1
      if (ime.gt.imee) go to 9011
      mgg(mg,inter)=mgg(mg,inter)+1
      m=mgg(mg,inter)
      if (m.gt.mee) go to 9001
      ima(m,mg,inter)=ime
      if (m.ne.1) go to 65
      imga(inter)=imga(inter)+1
      mggo(imga(inter),inter)=mg
   65 continue
c
c        store coupling constant g**2
      c(1,ime)=cc(1)
c        store coupling constant f*g
      c(3,ime)=cc(1)*cc(2)
      if (inter.eq.2.and.mg.eq.10) c(1,ime)=c(1,ime)+c(3,ime)
c        store coupling constant f**2
      c(2,ime)=cc(2)*c(3,ime)
      if (inter.eq.1.and.mg.eq.10)
     1  c(1,ime)=c(1,ime)+c(3,ime)*2.d0+c(2,ime)
c        store meson mass sqare in units of nucleon mass square
      c(4,ime)=cc(3)*cc(3)*dwnq
c
c        test iso-spin
      icc=cc(4)
      if (icc.ne.0.and.icc.ne.1) go to 9004
c         store isospin as logical constant
      if (icc.eq.1) indc(1,ime)=.true.
c        store and test iprsp for meson/delta/nucleon
      icc=cc(5)
      iccc=mod(icc,100)
c        iprsp for nucleons
      ic(1,ime)=mod(iccc,10)
c        ispd for deltas
      ic(2,ime)=iabs(iccc/10)
c        ispm for mesons
      ic(3,ime)=iabs(icc/100)
      if (iabs(ic(1,ime)).gt.8) go to 9005
      if (iabs(ic(1,ime)).ge.2.and.iabs(ic(1,ime)).le.7) endep=.true.
c
c        index values for further storing
      mi=4
      mm=5
      go to 61
c
c
c
c
c        write cut-off parameters
c        ------------------------
c
c
c
c
   70 write (kwrite,10005) name,cc
c
c
c        check if individuel cut or general cut
c
      if (name(1).eq.cut) go to 73
c        case of general cut-off
      if (indcut) go to 90
      if (imec.ge.ime) go to 61
      imac=imec+1
      imec=ime
      if (imac.lt.imb) imac=imb
      go to 90
c        case of individuel cut-off
   73 imac=ime
      imec=ime
      if (zerocp) go to 61
c
c        save present values of indices
c
   90 indcut=.true.
      if (cc(1).eq.0.d0) go to 61
      mix=mi
      mmx=mm
c
c        start loop of mesons, which present cut-off refers to
c
      do 1095 im=imac,imec
      mi=mix
      mm=mmx
c
c
c
c
c        store cut-off parameters
c        ------------------------
c
c
c
c
c        store typ of cut-off
      ic(mi,im)=cc(1)
      ityp=ic(mi,im)
      if (ityp.lt.1.or.ityp.gt.9) go to 9002
c        store and test typ of propagator of cut-off
      ic(mi+1,im)=cc(2)
      if (ic(mi+1,im).lt.0.or.ic(mi+1,im).gt.8) go to 9006
      if (ic(mi+1,im).ge.2.and.ic(mi+1,im).le.7) endep=.true.
      go to (100,100,300,400,400,400,700,800,900),ityp
c
c
c        cut-off of dipole type
c        **********************
c
c
c        store and test exponent of cut-off
  100 ic(mi+2,im)=cc(3)
      if (ic(mi+2,im).lt.0) go to 9009
      if (ic(mi+2,im).gt.0) go to 101
c        exponent is zero, omit cut-off
      ic(mi,im)=0
      ic(mi+1,im)=0
      go to 1000
c        store cut-off mass for denominator
  101 c(mm+1,im)=cc(4)*cc(4)*dwnq
c        store numerator of cut-off
      c(mm,im)=c(mm+1,im)
      if (ityp.eq.2)     c(mm,im)=c(mm,im)-c(4,im)
      mi=mi+3
      mm=mm+2
      go to 1000
c
c
c        cut-off of regge type /schierholz/
c        **********************************
c
c
  300 c(mm,im)=2.d0/dsqrt(4.d0-cc(3)*cc(3)*dwnq)
      c(mm+1,im)=cc(4)-1.d0
      c(mm+2,im)=cc(5)*wnq*1.d-6
      mi=mi+2
      mm=mm+3
      go to 1000
c
c
c        eikonal form factor
c        *******************
c
c
c        store gamma as -4.*gamma
  400 c(mm,im)=-cc(3)*4.d0
      ieik=ityp-3
c
c        compute and store normalization factor of t-form
      d=-c(4,im)
      d1=dsqrt(-d)
      d2=dsqrt(4.d0+d)
      c(mm+1,im)=dexp(4.d0*cc(3)*(2.d0+d)/(d1*d2)*datan(d1/d2))
c
      go to (490,412,412),ieik
c
c        compute and store normalization factor of u-form
c        and t- * u-form
  412 ieik2=cc(4)
      if (ieik2.eq.0) ieik2=1
      ic(mi+2,im)=ieik2
      go to (421,422,423,421,422,423),ieik2
c
  421 ess=4.d0*(1.d0+cc(5)*dwn)**2
      go to 450
c
  422 endep=.true.
c
  423 go to 490
c
  450 d=-(4.d0-ess)
      if (ieik2.le.3) d=d+c(4,im)
      if (d.eq.0.d0) go to 460
      d2=dsqrt(4.d0+d)
      if (d.lt.0.d0) go to 470
      d1=dsqrt(d)
      c(mm+2,im)=dexp(4.d0*cc(3)*(2.d0+d)/(d1*d2)*dlog(0.5d0*(d1+d2)))
      go to 480
c
  460 c(mm+2,im)=dexp(2.d0*cc(3))
      go to 480
c
  470 d1=dsqrt(-d)
      c(mm+2,im)=dexp(4.d0*cc(3)*(2.d0+d)/(d1*d2)*datan(d1/d2))
c
  480 if (ieik.eq.3) c(mm+2,im)=c(mm+2,im)*c(mm+1,im)
c
  490 mi=mi+3
      mm=mm+3
      go to 1000
c
c
c        exponential form factor
c        ***********************
c
c
c        check exponent
  700 if (cc(3).lt.0.d0) go to 9009
      if (cc(3).gt.0.d0) go to 701
c        exponent is zero, omit cutoff
      ic (mi,im)=0
      ic (mi+1,im)=0
      go to 1000
c        compute constant factor for argument of exponential function
  701 c(mm,im)=cc(3)*wnq/(cc(4)*cc(4))
      mi=mi+2
      mm=mm+1
      go to 1000
c
c
c        cloudy bag form factor
c        ***********************
c
c
c        check exponent
  800 if (cc(3).lt.0.d0) go to 9009
      if (cc(3).gt.0.d0) go to 801
c        exponent is zero, omit cutoff
      ic (mi,im)=0
      ic (mi+1,im)=0
      go to 1000
c        store exponent
  801 ic(mi+2,im)=cc(3)
c        store cutoff radius
      c(mm,im)=cc(4)*wn/uf
      mi=mi+3
      mm=mm+1
      go to 1000
c
c
c        propagator of mass-distributed meson
c        ************************************
c
c
c        spin of meson
  900 icc=cc(3)
      ic(mi+2,im)=icc
c        full width
      c(mm,im)=cc(4)*dwn
c        2*pion-mass squared
      c(mm+1,im)=cc(5)*cc(5)*dwnq
c        recalculate width
      d=c(4,im)-c(mm+1,im)
      c(mm,im)=c(mm,im)*dsqrt(c(4,im))/dsqrt(d)
      if (icc.lt.1) go to 910
      do 905 i=1,icc
  905 c(mm,im)=c(mm,im)/d
  910 mi=mi+3
      mm=mm+2
      go to 1000
c
c
c
c
c        end cut-offs
c        ************
c
c        test dimensions
 1000 if (mi.gt.mice.or.mm-1.gt.mme) go to 9010
c
c
 1095 continue
      go to 61
c
c
c
c
c        last card
c        ---------
c        ---------
c
c
c
c
c        write end mesons
 2000 imaa(inter)=imb
      imea(inter)=ime
      imb=ime+1
      write (kwrite,10004) name
      write (kwrite,10008)
      write (kwrite,10008)
c
c
c
c
      return
c
c
c
c        errors
c        ------
c        ------
c
c
c
c
 9000 write (kwrite,19000) name(1)
19000 format (1h0/////'0error in obpar:  meson-group   ',a4,'   does not
     1 exist in this program.'/'0execution terminated.'////)
      go to 9999
c
c
 9001 write (kwrite,19001)
19001 format (1h0/////'0error in obpar: too many mesons within a meson-g
     1roup with respect to'/'0the given dimensions. execution terminatee
     2.'////)
      go to 9999
c
c
 9002 write (kwrite,19002) cc(1)
19002 format (1h0/////'0error in obpar: cut-off typ',f10.4,'  does not e
     1xist in this program.'/'0execution terminated.'////)
      go to 9999
c
c
 9003 write (kwrite,19003) iftyp
19003 format (1h0/////'0error in obpar: factor typ has the non-permissib
     1le value',i4,' .'/'0execution terminated.'////)
      go to 9999
c
c
 9004 write (kwrite,19004) cc(4)
19004 format (1h0/////'0error in obpar: isospin has the non-permissible
     1value',f10.4,'  .'/'0execution terminated.'////)
      go to 9999
c
c
 9005 write (kwrite,19005) cc(5)
19005 format (1h0/////'0error in obpar: iprop/spe has the non-permissibl
     1e value',f10.4,'  .'/'0execution terminated.'////)
      go to 9999
c
c
 9006 write (kwrite,19006) cc(2)
19006 format (1h0/////'0error in obpar: the index for the propagator of
     1the cut-off has the'/'0non-permissible value',f10.4,'  . execution
     2 terminated.'////)
      go to 9999
c
c
 9009 write (kwrite,19009)
19009 format (1h0/////'0error in obpar: the exponent of the cut-off is l
     1ess than zero.'/'0execution terminated.'////)
      go to 9999
c
c
 9010 write (kwrite,19010)
19010 format (1h0/////'0error in obpar: too many cut-off parameters with
     1 respect to the given'/'0dimensions. execution terminated.'////)
      go to 9999
c
c
 9011 write (kwrite,19011)
19011 format (1h0/////'0error in obpar:  too many mesons with respect to
     1 the dimensions given'/'0to this program. execution terminated.'
     2////)
      go to 9999
c
c
 9999 stop
      end
      subroutine obnn
c
c
c        one-boson-exchange nn-nn interaction;
c        version which uses numerical integration
c
c        author:      r. machleidt
c                     institut fuer theoretische kernphysik bonn
c                     nussallee 14-16
c                     d - 5300  bonn, w. germany
c
c
      implicit real*8 (a-h,o-z)
c
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
c        arguments and values of this subroutine
c
      common /cpot/   v(6),xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep,label
c
c
c        this has been the end of the common-blocks containing
c        the arguments and values of this subroutine in the case of
c        no energy-dependence of the potential;
c        in case of energy-dependence look for the common-block /cped/
c        in obai and obaa.
c
c        specifications for these two common blocks
c
      logical heform,sing,trip,coup,endep
c
c
c        common block for all ob-subroutines
c
      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                ez1,ez2,ct(96),wt(96),
     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c
c        further specifications
c
      data pi/3.141592653589793d0/
      data  d3/0.3333333333333333d0/
      data td3/0.6666666666666667d0/
      data fd3/1.3333333333333333d0/
      character*4 mesong(12)
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     1                   '1-  ','1-t ','1-tt','1-st','1-ss',
     2                    '1+  ','2+  '/
      logical index
      data index/.false./
      logical indmg(12)
      data indmg/12*.false./
c
c
c
c
      inter=1
c
c
c
c
c        call subroutine obpar once and only once
c
c
      if (index) go to 50
      index=.true.
      if (indpar(inter)) go to 40
c
c
      call obpar
c
c
   40 wdd(inter)=0.d0
      iftgo=ift(inter)+1
      dwn=1.d0/wnn(inter)
      iman=imaa(inter)
      imen=imea(inter)
c
c
c        prepare constant over-all factor
c
      fac=1.d0/(2.d0*pi)*dwn*dwn
c     --------------------------
c
c
c
c
c
c
c
c        prepare expressions depending on x and y
c        ----------------------------------------
c        ----------------------------------------
c
c
c
c
   50 xa=xmev*dwn
      ya=ymev*dwn
      indxy=.false.
      x=xa
      xx=x*x
      y=ya
      yy=y*y
      xy2=x*y*2.d0
      xxpyy=xx+yy
      ex=dsqrt(1.d0+xx)
      ey=dsqrt(1.d0+yy)
      eem12=(ex*ey-1.d0)*2.d0
c
c
c
c
   55 xy=xy2*0.5d0
c     dxy=1.d0/xy
      ee=ex*ey
      ree=dsqrt(ee)
      eem1=ee-1.d0
      eme=ex-ey
      emeh=eme*0.5d0
      emehq=emeh*emeh
      eep1=ee+1.d0
       epe=ex+ey
      xxyy=xx*yy
c
c
c
c
c        prepare over-all factor
c
c
      go to (70,71,72),iftgo
c
c        no additional factor
c
   70 fff=fac
      go to 90
c
c        minimal relativity
c
   71 fff=fac/ree
      go to 90
c
c        factor m/e*m/e
c
   72 fff=fac/ee
c
c
c
c
c
c
   90 do 93 iv=1,6
   93 v(iv)=0.d0
      do 95 il=iman,imen
      do 95 iv=1,32
   95 vj(iv,il)=0.d0
c
c
c
c
c        contributions of mesons
c        -----------------------
c        -----------------------
c
c
c
c
      do 1995 img=1,mge
      mg=mggo(img,inter)
      if (mg.eq.0) go to 2000
      me=mgg(mg,inter)
      go to (100,200,300,400,500,600,600,600,900,1000,1100,1200),mg
c
c
c
c        0-  , pseudo-scalar mesons
c        --------------------------
c
c
c
c
  100 mc=1
c
      ff=1.d0
      f(1)=eem1
      f(2)=-xy
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-eme
      f(8)=-f(7)
c
      call obstr(1,1,me)
      go to 1995
c
c
c
c        0-t , pseudo-vector mesons
c        --------------------------
c
c
c
c
  200 mc=1
c
      ff=1.d0
      f(1)=eem1+emehq*(ee+3.d0)
      f(2)=-xy+emehq*xy
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-eme-eme*(emehq+eem1)
      f(8)=-f(7)
c
      call obstr(1,1,me)
      go to 1995
c
c
c
c        0-st, pseudo-scalar mesons in static limit
c        ------------------------------------------
c
c
c
c
  300 mc=1
c
      ff=1.d0
      f(1)=xxpyy*0.5d0
      f(2)=-xy
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-(xx-yy)*0.5d0
      f(8)=-f(7)
c
      call obstr(1,1,me)
      go to 1995
c
c
c
c
c        0+  , scalar mesons
c        -------------------
c
c
c
c
  400 mc=1
c
      ff=1.d0
      f(1)=-eep1
      f(2)=xy
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=epe
      f(8)=f(7)
c
      call obstr(1,1,me)
      go to 1995
c
c
c
c
c        0+st, scalar mesons in static limit
c        -----------------------------------
c
c
c
c
  500 mc=1
c
c
c
c        central term  '1'  only
c
c
c
      ff=1.d0
      f(1)=-2.d0
      f(2)=0.d0
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-f(1)
      f(8)=f(7)
      f(9)=f(2)
      f(10)=f(2)
      f(11)=f(2)
c
c
c
c        central term  ' k**2 + p**2 '
c
c
c
      f(2)=f(2)+xy
      f(10)=f(10)+xy
      f(11)=f(11)-xy
c
c
c
c        spin orbit
c
c
c
      f(4)=f(4)+xy
      f(5)=f(5)+xy
      f(10)=f(10)-xy
      f(11)=f(11)+xy
c
      call obstr(2,1,me)
c
c
c
c        case of additional sigma-l
c
c
c
c     mc=-1
c     xxyy8=xxyy/8.d0
c     f(1)=-xxyy8
c     f(2)=0.d0
c     f(3)=f(1)
c     f(4)=f(2)
c     f(5)=f(2)
c     f(6)=xxyy8
c     f(7)=f(1)
c     f(8)=f(7)
c     f(9)=f(6)*2.d0
c
c     call obstr(4,1,me)
      go to 1995
c
c
c
c
c        1-  , vector mesons
c        -------------------
c
c
c
c
c        vector coupling
c
c
c
c
  600 mc=1
c
      ff=2.d0
      f(1)=eem1+ee
      f(2)=0.d0
      f(3)=ee
      f(4)=xy
      f(5)=xy2
      f(6)=1.d0
      f(7)=-ey
      f(8)=-ex
c
      call obstr(1,1,me)
      if (mg.eq.7) go to 720
      if (mg.eq.8) go to 810
c
c
c
c
c        tensor coupling
c
c
c
c
      mc=2
c
c
      ff=0.5d0
        e1=-xxpyy+xxyy
      f(1)=eem1*(xxpyy+6.d0)+e1
      f(2)=-(xxpyy+4.d0)*xy
      f(3)=eem1*(xxpyy+2.d0)+e1
        e2=xxpyy+ee
      f(4)=-(1.d0+e2)*xy
      f(5)= (3.d0-e2)*xy
      f(6)=eem1*(xxpyy-2.d0)-xxpyy
        e3=2.d0*eme
        e4=yy*(ey-2.d0*ex)+xx*(ex-2.d0*ey)
      f(7)= e3-e4
      f(8)=-e3-e4
c        factors for additional terms
      f(9)=-xxyy
      f(10)=(1.d0+ee)*xy
      f(11)=-epe*xy
c
      call obstr(2,1,me)
c
c
c
c
c        vector-tensor coupling
c
c
c
c
      mc=3
c
      ff=2.d0
      f(1)=3.d0*eem1-xxpyy
      f(2)=-xy
      f(3)=eem1-xxpyy
      f(4)=-f(2)
      f(5)=3.d0*xy
      f(6)=-(eem1+xxpyy)
        e1=yy*ex+xx*ey
      f(7)= eme+e1
      f(8)=-eme+e1
c
      call obstr(1,1,me)
      go to 1995
c
c
c
c
c        1-t , vector mesons a la gtg
c        ----------------------------
c
c
c
c
c        tensor coupling
c
c
c
c
  720 mc=2
c
      ff=0.25d0
      f(1)=(3.d0*ee+1.d0)*xxpyy
      f(2)=-(6.d0*ee+2.d0-xxpyy)*xy
      f(3)=eem1*xxpyy+4.d0*xxyy
      f(4)=-(4.d0*ee+xxpyy)*xy
      f(5)=(4.d0-3.d0*xxpyy)*xy
      f(6)=6.d0*xxyy-(ee+3.d0)*xxpyy
      f(7)=(ex+3.d0*ey)*xx+eme*yy
      f(8)=(ey+3.d0*ex)*yy-eme*xx
c        factors for additional terms
      f(9)=-2.d0*xxyy
      f(10)=eep1*xy2
      f(11)=-epe*xy2
c
      call obstr(2,1,me)
c
c
c
c
c        vector-tensor coupling
c
c
c
c
      mc=3
c
      ff=1.d0
      f(1)=xxpyy
      f(2)=-xy2
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=6.d0*xy
      f(6)=3.d0*f(3)
      f(7)=(ex*yy+3.d0*ey*xx)
      f(8)=(ey*xx+3.d0*ex*yy)
c
      call obstr(1,1,me)
      go to 1995
c
c
c
c
c        1-tt, vector mesons with consequent ignoration of retardation
c        -------------------------------------------------------------
c
c
c
c
c        vector coupling (additional term)
c
c
c
c
  810 mc=0
c
      f(1)=eep1
      f(2)=xy
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-epe
      f(8)=f(7)
        e1=eme*eme
c
      do 815 mx=1,me
      ff=e1/c(4,ima(mx,mg,inter))
  815 call obstr(1,mx,mx)
c
c
c
c
c
c
c        tensor coupling and vector-tensor coupling
c
c
c
c
      go to 720
c
c
c
c
c        1-st, vector mesons in static limit
c        -----------------------------------
c
c
c
c        vector coupling
c
c
c
c
  900 mc=1
c
c
c
c        central term  '1'  only
c
c
c
      ff=1.d0
      f(1)=2.d0
      f(2)=0.d0
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-f(1)
      f(8)=f(7)
      f(9)=f(2)
      f(10)=f(2)
      f(11)=f(2)
      xxpyyh=xxpyy*0.5d0
c
c
c**** goto 990
c
c
c        central term  ' k**2 + p**2 '
c
c
c
      f(1)=f(1)+xxpyyh
      f(2)=f(2)+xy2
      f(3)=f(3)+xxpyyh
      f(6)=f(6)+xxpyyh
      f(7)=f(7)-xxpyyh
      f(8)=f(8)-xxpyyh
      f(10)=f(10)+xy2
      f(11)=f(11)-xy2
c
c
c
c        spin - spin
c
c
c
      f(1)=f(1)+xxpyyh*3.d0
      f(2)=f(2)-xy*3.d0
      f(3)=f(3)-xxpyyh
      f(6)=f(6)-xxpyyh
      f(7)=f(7)+xxpyyh
      f(8)=f(8)+xxpyyh
      f(10)=f(10)+xy
      f(11)=f(11)-xy
c
c
c
c        tensor  with  spin-spin-contribution
c
c
c
      f(1)=f(1)-xxpyyh
      f(2)=f(2)+xy
      f(3)=f(3)+xxpyyh
      f(4)=f(4)-xy
      f(5)=f(5)+xy
      f(6)=f(6)-xxpyyh
      f(7)=f(7)+(xx-yy)*0.5d0
      f(8)=f(8)-(xx-yy)*0.5d0
c
c
c
c        spin - orbit
c
c
c
      xy3=xy*3.d0
      f(4)=f(4)+xy3
      f(5)=f(5)+xy3
      f(10)=f(10)-xy3
      f(11)=f(11)+xy3
c
  990 continue
c
      call obstr(2,1,me)
c
c**** goto 1995
c
c
c        case of additional sigma-l
c
c
c
c     mc=-1
c     xxyy8=xxyy/8.d0
c     f(1)=xxyy8
c     f(2)=0.d0
c     f(3)=f(1)
c     f(4)=f(2)
c     f(5)=f(2)
c     f(6)=-xxyy8
c     f(7)=f(1)
c     f(8)=f(7)
c     f(9)=f(6)*2.d0
c
c     call obstr(4,1,me)
c
c
c
c
c        tensor coupling
c
c      ( 1-ss , rho-meson in static limit, tensor coupling only )
c
c
 1000 mc=2
c
      if (mg.eq.10) mc=1
c
c
c        spin - spin
c
c
c
      f(1)=xxpyyh*3.d0
      f(2)=-xy*3.d0
      f(3)=-xxpyyh
      f(4)=0.d0
      f(5)=0.d0
      f(6)=-xxpyyh
      f(7)=xxpyyh
      f(8)=xxpyyh
      f(9)=0.d0
      f(10)=xy
      f(11)=-xy
c
c
c
c        tensor  with  spin-spin-contribution
c
c
c
      f(1)=f(1)-xxpyyh
      f(2)=f(2)+xy
      f(3)=f(3)+xxpyyh
      f(4)=f(4)-xy
      f(5)=f(5)+xy
      f(6)=f(6)-xxpyyh
      f(7)=f(7)+(xx-yy)*0.5d0
      f(8)=f(8)-(xx-yy)*0.5d0
c
      call obstr(2,1,me)
      if (mg.eq.10) go to 1995
c
c
c
c        case of additional sigma-l
c
c
c
c     f(1)=xxyy
c     f(2)=0.d0
c     f(3)=f(1)
c     f(4)=f(2)
c     f(5)=f(2)
c     f(6)=-xxyy
c     f(7)=f(1)
c     f(8)=f(7)
c     f(9)=f(6)*2.d0
c
c     call obstr(4,1,me)
c
c
c
c
c        vector-tensor coupling
c
c
c
c
      mc=3
c
c
c
c        central term  ' k**2 + p**2 '
c
c
c
      f(1)=-xxpyy
      f(2)=xy2
      f(3)=f(1)
      f(4)=0.d0
      f(5)=0.d0
      f(6)=f(1)
      f(7)=-f(1)
      f(8)=f(7)
      f(9)=0.d0
      f(10)=f(2)
      f(11)=-f(2)
c
c
c
c        spin - spin
c
c
c
      f(1)=f(1)+xxpyy*3.d0
      f(2)=f(2)-xy2*3.d0
      f(3)=f(3)-xxpyy
      f(6)=f(6)-xxpyy
      f(7)=f(7)+xxpyy
      f(8)=f(8)+xxpyy
      f(10)=f(10)+xy2
      f(11)=f(11)-xy2
c
c
c
c        tensor  with  spin-spin-contribution
c
c
c
      f(1)=f(1)-xxpyy
      f(2)=f(2)+xy2
      f(3)=f(3)+xxpyy
      f(4)=f(4)-xy2
      f(5)=f(5)+xy2
      f(6)=f(6)-xxpyy
      f(7)=f(7)+(xx-yy)
      f(8)=f(8)-(xx-yy)
c
c
c
c        spin - orbit
c
c
c
      xy4=xy*4.d0
      f(4)=f(4)+xy4
      f(5)=f(5)+xy4
      f(10)=f(10)-xy4
      f(11)=f(11)+xy4
c
      call obstr(2,1,me)
c
c
c
c        case of additional sigma-l
c
c
c
c     f(1)=xxyy
c     f(2)=0.d0
c     f(3)=f(1)
c     f(4)=f(2)
c     f(5)=f(2)
c     f(6)=-xxyy
c     f(7)=f(1)
c     f(8)=f(7)
c     f(9)=f(6)*2.d0
c
c     call obstr(4,1,me)
      go to 1995
c
c
c
c
c
c        1+  , a1-meson
c        --------------
c
c
c
c
 1100 mc=1
c
      ff=2.d0
      f(1)=eep1+ee
      f(2)=0.d0
      f(3)=-ee
      f(4)=-xy
      f(5)=xy2
      f(6)=-1.d0
      f(7)=ey
      f(8)=ex
c
      call obstr (1,1,me)
      go to 1995
c
c
c
c
c
c        2+  , f-meson
c        -------------
c
c
c
c
c        g1**2 coupling
c
c
c
c
 1200 mc=1
c
      ff=-8.d0
      ee4=ee*4.d0
      xxyy4=xxyy*4.d0
      eep143=eep1*fd3
        e1=2.d0*(xxpyy+td3)
      f(1)= eep143 +(3.d0*ee+2.d0)*xxpyy+xxyy4
      f(2)=(ee4+td3+xxpyy)*xy
      f(3)=3.d0*xxyy+eep1*e1
      f(4)=(2.d0*xxpyy+3.d0*eep1-d3)*xy
      f(5)=(ee4+11.d0*d3+3.d0*xxpyy)*xy
      f(6)= eep143 +(ee+2.d0)*xxpyy+xxyy4
        e2=-epe*e1
      f(7)=e2+xx*ex
      f(8)=e2+yy*ey
c        factors for additional terms
      f(9)=xxyy
      f(10)=ee*xy
      f(11)=xy
      f(12)=-ey*xy
      f(13)=-ex*xy
c
      call obstr(3,1,me)
      go to 1995
c
c
c
c
c        this has been the end of the contributions of mesons
c        ----------------------------------------------------
c
c
c
c
c        errors and warnings
c        -------------------
c
c
c
c
 9000 if (indmg(mg)) go to 1995
      write (kwrite,19000) mesong(mg)
19000 format(1h0////'0warning in obnn: meson-group  ',a4,'  does not exi
     1st in this program.'/'0contribution ignored. execution continued.'
     2////)
      indmg(mg)=.true.
c
c
c
c
 1995 continue
c
c
c
c
c        add up contributions of mesons
c        ------------------------------
c
c
c
c
 2000 do 2005 il=iman,imen
      do 2005 iv=1,6
 2005 v(iv)=v(iv)+vj(iv,il)
c
c
c
c
      return
      end
      subroutine obstr (icase,max,mex)
c
c        obstr computes the structure of one-boson-exchanges
c
c
      implicit real*8 (a-h,o-z)
c
c
c        common blocks
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
c
c
c        common block for all ob-subroutines
c
      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                ez1,ez2,ct(96),wt(96),
     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c     further specifications
c
      dimension vv(32)
      dimension tt(2,3)
      data jj/-1/
      logical index
      data index/.false./
      logical indiso
c
c
c
c
      if (index) go to 50
      index=.true.
c
c
      tt(1,1)=1.d0
      tt(2,1)=-3.d0
c
      do 1 ii=2,3
      do 1 i=1,2
    1 tt(i,ii)=1.d0
c
c
c
c
c
   50 do 1095 m=max,mex
      im=ima(m,mg,inter)
c
c
      if (mc.ne.1) go to 60
c
c
c
c
c        call integrals
c        --------------
c
c
c
c
      call obai
c
c
c
c
   60 if (mc.lt.1) mc=1
c
      if (c(mc,im).eq.0.d0) go to 1095
c
c
c
c
      go to (100,200,300),inter
c
c
c
c
c        nn-nn helicity amplitudes /combinations/
c        ----------------------------------------
c
c
c
c
c        ground structure (a factor of 2 is included in v5 and v6)
c
c
  100 ive=6
c
      vv(1)=f(1)*ai(1,m)+f(2)*ai(2,m)
      vv(2)=f(3)*ai(1,m)+f(4)*ai(3,m)
      vv(3)=f(5)*ai(1,m)+f(6)*ai(2,m)
      vv(4)=f(4)*ai(1,m)+f(3)*ai(3,m)
      vv(5)=f(7)*ai(4,m)
      vv(6)=f(8)*ai(4,m)
c
c
      go to (1000,120,130,140),icase
c
c
c        additional terms  in case of tensor coupling
c
c
  120 vv(1)=vv(1)+f(9)*ai(5,m)
      vv(2)=vv(2)+f(10)*ai(2,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(10)*ai(5,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(10)*ai(6,m)
         e1=f(11)*ai(7,m)
      vv(5)=vv(5)+e1
      vv(6)=vv(6)+e1
      go to 1000
c
c
c        additional terms in case of 2+ mesons
c
c
  130 vv(2)=vv(2)+f(10)*ai(2,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(11)*ai(5,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(10)*ai(6,m)
      vv(5)=vv(5)+f(12)*ai(7,m)
      vv(6)=vv(6)+f(13)*ai(7,m)
      go to 1000
c
c
c        additional terms in case of sigma-l in static limit
c
c
  140 vv(1)=vv(1)+f(6)*ai(5,m)
      vv(2)=vv(2)+f(1)*ai(5,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(1)*ai(11,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(1)*ai(12,m)
      vv(5)=vv(5)+f(6)*ai(13,m)
      vv(6)=vv(6)+f(6)*ai(13,m)
      go to 1000
c
c
c
c
c        nn-nd helicity amplitudes
c        -------------------------
c
c
c
c
  200 ive=16
c
      ai3m1=ai(3,m)-ai(1,m)
      ai6m2=ai(6,m)-ai(2,m)
      ai3p1=ai(3,m)+ai(1,m)
      ai6p2=ai(6,m)+ai(2,m)
c
c
      vv( 1)= f( 1)* ai( 4,m) + f( 2)* ai( 7,m)
      vv( 2)= f( 3)* ai( 1,m) + f( 4)* ai( 2,m) + f( 5)* ai( 5,m)
      vv( 3)= f( 6)* ai( 4,m) + f( 7)* ai( 7,m)
      vv( 4)= f( 8)* ai( 8,m)
      vv( 5)= f( 9)* ai( 8,m)
      vv( 6)=-f(10)* ai( 4,m) + f(11)* ai( 7,m)
      vv( 7)= f(12)* ai( 1,m) - f(13)* ai( 2,m) + f(14)* ai( 5,m)
      vv( 8)=-f(15)* ai( 4,m) + f(16)* ai( 7,m)
c
      vv( 9)= f(17)* ai3m1    + f(18)* ai6m2
      vv(10)= f(19)* ai( 4,m) + f(20)* ai( 7,m)
      vv(11)= f(21)* ai3p1    + f(22)* ai6p2
      vv(12)= f(23)* ai( 9,m)
      vv(13)=-f(24)* ai(10,m)
      vv(14)=-f(25)* ai3m1    + f(26)* ai6m2
      vv(15)=-f(27)* ai( 4,m) + f(28)* ai( 7,m)
      vv(16)=-f(29)* ai3p1    + f(30)* ai6p2
      go to 1000
c
c
c
c        nn-dd helicity amplitudes
c        -------------------------
c
c
c
c
  300 ive=32
c
      ai31p=ai( 3,m)+ai( 1,m)
      ai31m=ai( 3,m)-ai( 1,m)
      ai62p=ai( 6,m)+ai( 2,m)
      ai62m=ai( 6,m)-ai( 2,m)
      aic5p=ai(12,m)+ai( 5,m)
      aic5m=ai(12,m)-ai( 5,m)
c
c
      vv( 1)= f( 1)*(ai(1,m)+ai(2,m)) + f( 2)*(ai(2,m)+ai(5,m)) +
     1        f( 3)*(ai(5,m)+ai(11,m))
c
      vv( 2)= f( 4)* ai( 4,m) + f( 5)* ai( 7,m) + f( 6)* ai(13,m)
      vv( 3)= f( 7)* ai( 8,m) + f( 8)* ai(14,m)
      vv( 4)= f( 9)* ai(17,m)
      vv( 5)=vv( 2)
      vv( 6)= f(10)* ai( 1,m) + f(11)* ai( 2,m) + f(12)* ai( 5,m) +
     1        f(13)* ai(11,m)
      vv( 7)= f(14)* ai( 4,m) + f(15)* ai( 7,m) + f(16)* ai(13,m)
      vv( 8)=-f(17)* ai( 8,m) + f(18)* ai(14,m)
      vv( 9)=vv( 3)
      vv(10)=vv( 7)
      vv(11)=-f(19)* ai( 1,m) + f(20)* ai( 2,m) - f(21)* ai( 5,m) +
     1        f(22)* ai(11,m)
      vv(12)= f(23)* ai( 4,m) - f(24)* ai( 7,m) + f(25)* ai(13,m)
      vv(13)=vv( 4)
      vv(14)=vv( 8)
      vv(15)=vv(12)
c
      vv(16)=-f(26)*(ai(1,m)-ai(2,m)) + f(27)*(ai(2,m)-ai(5,m)) -
     1        f(28)*(ai(5,m)-ai(11,m))
c
c
      vv(17)= f(29)* ai( 4,m) + f(30)* ai( 7,m) + f(31)* ai(13,m)
      vv(18)= f(32)* ai31p    + f(33)* ai62p    + f(34)* aic5p
      vv(19)= f(35)* ai( 9,m) + f(36)* ai(15,m)
      vv(20)= f(37)* ai(18,m)
      vv(21)= f(44)* ai31m    - f(45)* ai62m    + f(46)* aic5m
      vv(22)= f(38)* ai( 4,m) + f(39)* ai( 7,m) + f(40)* ai(13,m)
      vv(23)= f(41)* ai31p    + f(42)* ai62p    + f(43)* aic5p
      vv(24)=vv(19)
      vv(25)= f(47)* ai(10,m) - f(48)* ai(16,m)
      vv(26)= f(49)* ai31m    - f(50)* ai62m    + f(51)* aic5m
      vv(27)=vv(22)
      vv(28)=vv(18)
      vv(29)= f(52)* ai(19,m)
      vv(30)=vv(25)
      vv(31)=vv(21)
      vv(32)=vv(17)
c
c
c
c
 1000 if (inter.ge.2) go to 1040
c
c
c
c
c        set certain cases to zero in case of inter=1
c
      if (j.ne.0) go to 1021
      vv(2)=0.d0
      vv(4)=0.d0
      vv(5)=0.d0
      vv(6)=0.d0
c
 1021 if (.not.sing) vv(1)=0.d0
      if (.not.trip) vv(2)=0.d0
      if (coup) go to 1030
      do 1025 iv=3,6
 1025 vv(iv)=0.d0
c
 1030 if (heform) go to 1040
c
c
c        transformation into lsj-formalism in case of inter=1
c        (if requested)
      if (j.eq.jj) go to 1035
      jj=j
      aj=dfloat(j)
      aj1=dfloat(j+1)
      d2j1=1.d0/dfloat(2*j+1)
      arjj1=dsqrt(aj*aj1)
c
 1035 v3=vv(3)
      v4=vv(4)
      v5=vv(5)
      v6=vv(6)
      v34=-arjj1*(v3-v4)
      v56=arjj1*(v5+v6)
      vv(3)=d2j1*(aj1*v3+aj*v4-v56)
      vv(4)=d2j1*(aj*v3+aj1*v4+v56)
      vv(5)=d2j1*(v34-aj1*v5+aj*v6)
      vv(6)=d2j1*(v34+aj*v5-aj1*v6)
c
c
c        possible different sign depending on the convention used
      vv(5)=-vv(5)
      vv(6)=-vv(6)
c
c
c
c
c        multiply with factors
c        ---------------------
c
c
c
c
 1040 is=mod(j,2)+1
      it=mod(is,2)+1
      indiso=indc(1,im)
      cmc=c(mc,im)
      fc=fff*ff*cmc
      do 1045 iv=1,ive
c
c        multiply with coupling-constant and factors fff and ff
c
      vv(iv)=vv(iv)*fc
      if (inter.ge.2) go to 1045
c
c        multiply with isospin factor
c
      if (.not.indiso) go to 1045
      if (iv.eq.2) go to 1043
      vv(iv)=vv(iv)*tt(is,inter)
      go to 1045
 1043 vv(iv)=vv(iv)*tt(it,inter)
c
c     add up in case of several couplings for one meson-exchange
c     and store
 1045 vj(iv,im)=vj(iv,im)+vv(iv)
c
c
c        if single contributions to one meson are to be printed:
c     write (kwrite,10001) mesong(mg),mc
c     write (kwrite,10002) (vv(iv),iv=1,ive)
10001 format (' contribution from ',a4,' mc =',i2)
10002 format (33x,6d16.6)
c
c
 1095 continue
c
c
      return
      end
      subroutine obai
c
c        obai integrates over theta
c
c
      implicit real*8 (a-h,o-z)
c
      common /cpot/   v(6),xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep,label
      common /cpoted/ q0qmev,qfmev,pmev,uanspp,wsnspp,ucnspp,udnspp,
     1                znrl,zrel,smev,noced
      logical heform,sing,trip,coup,endep
      logical noced
c
c
c        common block for all ob-subroutines
c
      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                ez1,ez2,ct(96),wt(96),
     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c
c        further specifications
      dimension gi(7)
      data ige/7/
c
      dimension pj(7,96)
      real*4 aez,axy2,aomq,am1,am2,am
      data nnt/-1/,iinter/-1/,jj/-1/
      logical indj,indint,indee,indepe,indz,indepz
c
c
c
c
      if (inter.eq.iinter) go to 60
      iinter=inter
      indint=.false.
      min=mint(inter)
      max=maxt(inter)
c
      go to (51,51,53),inter
   51 igeint=5
      go to 55
   53 igeint=7
   55 continue
c
      wn=wnn(inter)
      dwn=1.d0/wn
      wnq=wn*wn
c        wd is the mass of the delta
c        wdd(..) is the array of the masses of the delta
      wd=wdd(inter)
      wdq=wd*wd
c
c
c
c
   60 if (j.eq.jj) go to 70
      jj=j
      indj=.false.
c
c
      aj=dfloat(j)
      aj1=dfloat(j+1)
      dj1=1.d0/aj1
      ajdj1=aj*dj1
      aaj=dsqrt(ajdj1)
c
c
      aj2=dfloat(j+2)
      aj3=dfloat(j+3)
      ajm1=dfloat(j-1)
      ajm2=dfloat(j-2)
      ajm3=dfloat(j-3)
c
c
      ajj1=aj*aj1
      ajj2=ajm1*aj2
      ajj3=ajm2*aj3
      ajja=aj*ajm3
      ajjb=aj*ajm1
c
      aajj=0.d0
      if (j.gt.1)
     1aajj=aj/dsqrt(ajj1*ajj2)
c
      aaj1=aajj*ajm1
      aaj2=aajj*aj1
      aaj3=aajj*2.d0
c
      if (j.gt.1) go to 62
      aajj=0.d0
      go to 63
   62 aajj=1.d0/(aj1*dsqrt(ajj2))
c
   63 aaj4=aajj*ajjb
      aaj5=aajj*aj1*2.d0
      aaj6=aajj*(ajj1+2.d0)
      aaj7=aajj*ajj2
c
      if (j.gt.2) go to 64
      aajj=0.d0
      go to 65
   64 aajj=-aj/dsqrt(ajj1*ajj2*ajj3)
c
   65 aaj8=aajj*(ajj1+6.d0)
      aaj9=aajj*ajj2
      aaj10=aajj*(ajja+2.d0)
      aaj11=aajj*(ajja-6.d0)
c
      if (j.gt.2) go to 66
      aajj=0.d0
      go to 67
   66 aajj=-1.d0/(aj1*dsqrt(ajj2*ajj3))
c
   67 aaj12=aajj*ajjb*ajm2
      aaj13=aajj*(aj*ajjb+4.d0*aj+12.d0)
      aaj14=aajj*(5.d0*ajj1+6.d0)
      aaj15=aajj*3.d0*ajj2
      aaj16=aajj*ajj3*ajm1
      aaj17=aajj*aj1*ajj3
      aaj18=aajj*2.d0*ajj3
c
c
c
c
c
c        find out appropriate number of gauss-points
c        -------------------------------------------
c
c
c
c
   70 c4=c(4,im)
      iprsp=ic(1,im)
c
c
c        prepare starting energy
      if (noced) go to 73
      noced=.true.
      indz=.false.
      indepz=.false.
      if (iprsp.lt.2) go to 74
      if (iprsp.ge.4) go to 72
   71 indz=.true.
      z=2.d0*dsqrt(wnq+q0qmev)
      go to 74
   72 if (iprsp.eq.8) go to 74
      indepz=.true.
      eppq=pmev*pmev
      epz=zrel
      go to 74
c
c
   73 if (indxy.and.indint) go to 85
   74 indint=.true.
      indee=.false.
      indepe=.false.
c
      ispm=ic(3,im)
      ez1=0.d0
      if (ispm.eq.1.and.iprsp.lt.2) go to 90
c
c         delta-nucleon mass difference in propagator
c
      ez1=(wd-wn)*dwn
c
c
      if (iprsp.lt.2) go to 90
c
c        prepare for propagators of non-covariant perturbation theory
   75 if (iprsp.ge.4) go to 76
      if (.not.indz) go to 71
      pmevq=0.d0
      ez=z
      go to 77
   76 if (iprsp.eq.8) go to 1800
      if (.not.indepz) go to 72
      pmevq=eppq
      ez=epz
   77 qfmevq=qfmev*qfmev
      qqx=xmev*xmev+pmevq
      ymevq=ymev*ymev
      qqy=ymevq+pmevq
      go to (1100,1200,1300),inter
c
c        ez1 for the nn case
c
 1100 ez1=spe(qqx,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,2,iprsp,0)
     1  + spe(qqy,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,2,iprsp,0)
     2  - ez
      ez1=ez1*dwn
      go to 80
c
c        ez1 for the nd case
c
 1200 pq=4.d0*pmevq/(wn+wd)**2
      ispd=ic(2,im)
      ed=spe(ymevq+pq*wdq,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wd,2,
     1 ispd,0)
      en=spe(ymevq+pq*wnq,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,2,
     1 iprsp,0)
      eza=spe(qqx,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,2,iprsp,0)
     1 - ez
      ez1=eza+en
      ez2=eza+ed
      ez1=ez1*dwn
      ez2=ez2*dwn
      go to 80
c
c        ez1 for the dd case
c
 1300 ispdd=ic(2,im)
      ez1=spe(qqx,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,2,iprsp,0)
     1  + spe(qqy,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wd,2,ispdd,0)
     2  - ez
      ez1=ez1*dwn
      go to 80
c
c
c        case iprsp=8
c
 1800 ez1=ex-ey
c
c
c        store ez1 and ez2
c
c
   80 if (iprsp.ge.3) go to 81
      indee=.true.
      ee1=ez1
      ee2=ez2
      go to 89
   81 iiprsp=iprsp
      indepe=.true.
      epe1=ez1
      epe2=ez2
      go to 89
c
c
c        get stored ez1 and ez2
c
c
   85 if (iprsp.lt.2) go to 90
      if (iprsp.ge.3) go to 86
      if (.not.indee) go to 75
      ez1=ee1
      ez2=ee2
      go to 89
   86 if (iprsp.ne.iiprsp) go to 75
      if (.not.indepe) go to 75
      ez1=epe1
      ez2=epe2
c
c
   89 aez=ez1
c
c
c        compute am
c
c
   90 axy2=xy2
      if (iprsp.ne.1) go to 91
      aomq=eem12+c4
      go to 92
   91 aomq=xxpyy+c4
c
   92 am=axy2/aomq
c
      if (iprsp.lt.2) go to 93
      am1=am
      am2=axy2/(aomq+aez*abs(aez))
      if (am2.lt.0.) go to 94
      am=amax1(am1,am2)
c
c
c        compute number of gausspoints (nt)
c
c
   93 if (am.gt.0.999) go to 94
c
c
      if (am.gt.0.85) am=am**(-log(1.-am)-0.9)
c
c
      nt=float(min)/(1.-am)+0.9
c
c
      if (nt.gt.max) nt=max
      go to 95
c
c
   94 nt=max
c
c
   95 nt=nt+j
c
c        compute nt, which is suitable for gset
c
      if (nt.le.16) go to 98
      if (nt.gt.24) go to 96
      nt=4*(nt/4)
      go to 98
   96 if (nt.gt.48) go to 97
      nt=8*(nt/8)
      go to 98
   97 nt=16*(nt/16)
      if (nt.gt.96) nt=96
c
   98 if (nt.eq.nnt.and.indj) go to 100
c
c
c
c
c        call gauss-points
c        -----------------
c
c
c
c
      call gset (-1.d0,1.d0,nt,ct,wt)
      nnt=nt
c
c
c
c
c        call legendre-polynoms if necessary
c        -----------------------------------
c
c
c
c
      indxy=.false.
      indj=.true.
      do 99 i=1,nt
      t=ct(i)
      call legp (pj(1,i),pj(3,i),t,j)
      pj(2,i)=pj(1,i)*t
      pj(4,i)=pj(2,i)*t
      pj(6,i)=pj(4,i)*t
      pj(5,i)=pj(3,i)*t
   99 pj(7,i)=pj(5,i)*t
c
c
c
c
c        call integrand
c        --------------
c
c
c
c
  100 call obaa
c
c
c
c
c        prepare for integration
c
c
c
c
      do 2001 ig=1,igeint
 2001 gi(ig)=0.d0
c
c
c
c
c        integration-loop of theta
c        -------------------------
c
c
c
c
      do 2005 i=1,nt
      do 2005 ig=1,igeint
 2005 gi(ig)=gi(ig)+pj(ig,i)*aa(i)
c
c
c
      if (j.ne.0) go to 2010
      gi(3)=0.d0
      gi(5)=0.d0
      gi(7)=0.d0
c
c
c
c
c        combinations of integrals
c        -------------------------
c
c
c
c
 2010 ai(1,m)=gi(1)
c
      ai(2,m)=gi(2)
      ai(3,m)= ajdj1*gi(2)+dj1*gi(3)
      gi23m  =gi(2)-gi(3)
      ai(4,m)=aaj*gi23m
c
c
      ai(5,m)=gi(4)
      ai(6,m)= ajdj1*gi(4)+dj1*gi(5)
      gi45m  =gi(4)-gi(5)
      ai(7,m)=aaj*gi45m
c
c
      if (inter.eq.1) go to 3000
c
c
      ai( 8,m)= aaj1*gi(4)-aaj2*gi(1)+aaj3*gi(5)
      aai1    = aaj4*gi(4)+aaj5*gi(1)-aaj6*gi(5)
      aai2    = aaj7*gi23m
      ai( 9,m)= aai2+aai1
      ai(10,m)= aai2-aai1
c
c
      if (inter.ne.3) go to 3000
c
c
      ai(11,m)=gi(6)
      ai(12,m)=ajdj1*gi(6)+dj1*gi(7)
      ai(13,m)=aaj*(gi(6)-gi(7))
c
      ai(14,m)= aaj1*gi(6)-aaj2*gi(2)+aaj3*gi(7)
      aai1    = aaj4*gi(6)+aaj5*gi(2)-aaj6*gi(7)
      aai2    = aaj7*gi45m
      ai(15,m)= aai2+aai1
      ai(16,m)= aai2-aai1
c
      ai(17,m)= aaj8*gi(7)-aaj9*gi(3)-aaj10*gi(6)+aaj11*gi(2)
      aai1    =-aaj12*gi(6)+aaj13*gi(2)-aaj14*gi(7)+aaj15*gi(3)
      aai2    = aaj16*gi(4)-aaj17*gi(1)+aaj18*gi(5)
      ai(18,m)= aai1-aai2
      ai(19,m)= aai1+aai2
c
c
c
 3000 return
      end
      subroutine obaa
c
c        obaa computes the propagators and the cutoffs of ob-exchanges
c
c
      implicit real*8 (a-h,o-z)
c
      common /cstate/ j,heform,sing,trip,coup,endep,label
      common /cpoted/ q0qmev,qfmev,pmev,uanspp,wsnspp,ucnspp,udnspp,
     1                znrl,zrel,smev,noced
      logical heform,sing,trip,coup,endep
      logical noced
c
c
c        common block for all ob-subroutines
c
      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                ez1,ez2,ct(96),wt(96),
     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c
c
c        further specifications
      dimension deltaq(96,3),cut(96)
      data iinter/-1/
      logical indp2
      data ssmev/-1.d0/
c
c
c
c
      if (inter.eq.iinter) go to 60
      iinter=inter
c
c        dwn is needed for the eikonal cutoff
c
      dwn=1.d0/wnn(inter)
c
c
c
c
c        delta square
c        ------------
c
c
c
c
   60 if (indxy) go to 1000
      indxy=.true.
      do 65 i=1,nt
      xy2t=xy2*ct(i)
c
c        retardation ignored
c
      deltaq(i,1)=xy2t-xxpyy
c     ----------------------
c
c
c        retardation incorporated
c
      deltaq(i,2)=xy2t-eem12
c     ----------------------
c
c
c        for on shell cutoff  ( ferchlaender )
c
   65 deltaq(i,3)=-xy2t-xxpyy
c     -----------------------
c
c
c
c        propagator
c        ----------
c        ----------
c
c
c
c
 1000 c4=c(4,im)
      iprsp=ic(1,im)
      if (iprsp.lt.0) go to 1400
      if (iprsp.ge.2) go to 1050
      iret=iprsp+1
c
      go to (1010,1020,1030), inter
c
c         propagator for the nn case
 1010 do 1011 i=1,nt
 1011 aa(i)=wt(i)/(c4-deltaq(i,iret))
      go to 1500
c         propagator for the nd case
 1020 do 1021 i=1,nt
      omq=c4-deltaq(i,iret)
      om=dsqrt(omq)
 1021 aa(i)=wt(i)*(1.d0/omq+1.d0/(om*(om+ez1)))*0.5d0
C**** APPROXIMATION USED BY CANTON:
C1021 aa(i)=wt(i)*(1.d0/omq+1.d0/(om*om+dsqrt(c4)*ez1))*0.5d0
      go to 1500
c         propagator for the dd case
 1030 do 1031 i=1,nt
      omq=c4-deltaq(i,iret)
      om=dsqrt(omq)
 1031 aa(i)=wt(i)/(om*(om+ez1))
C**** APPROXIMATION USED BY CANTON:
C1031 aa(i)=wt(i)/(om*om+dsqrt(c4)*ez1)
      go to 1500
c
c
c        starting energy dependent propagator
c
c
 1050 ispm=ic(3,im)
      go to (1100,1200,1300),inter
c
c        the propagator for the nn case
c
 1100 do 1105 i=1,nt
      omq=c4-deltaq(i,1)
      om=dsqrt(omq)
      oms=om
      if (ispm.le.2) go to 1105
      if (dabs(ez1).lt.1.d-12) go to 1105
      oms=oms+smp(-deltaq(i,1),qfmev,ispm)
 1105 aa(i)=wt(i)/(om*(oms+ez1))
      go to 1500
c
c        the propagator for the nd case
c
 1200 do 1205 i=1,nt
      omq=c4-deltaq(i,1)
      om=dsqrt(omq)
      oms=om
      if (ispm.le.2) go to 1205
      oms=oms+smp(-deltaq(i,1),qfmev,ispm)
 1205 aa(i)=wt(i)*(1.d0/(om*(oms+ez1))+1.d0/(om*(oms+ez2)))*0.5d0
      go to 1500
c
c        the propagator for the dd case
c
 1300 do 1305 i=1,nt
      omq=c4-deltaq(i,1)
      om=dsqrt(omq)
      oms=om
      if (ispm.le.2) go to 1305
      oms=oms+smp(-deltaq(i,1),qfmev,ispm)
 1305 aa(i)=wt(i)/(om*(oms+ez1))
      go to 1500
c
c
c        "no propagator"
c
 1400 do 1405 i=1,nt
 1405 aa(i)=wt(i)
c
c
 1500 continue
c
c
c
c
c
c        cut-offs
c        --------
c        --------
c
c
c
c
      mi=4
      mm=5
c
c
  999 ityp=ic(mi,im)
      if (ityp.eq.0) go to 2000
      iprspc=ic(mi+1,im)
      iret=iprspc+1
      if (iprspc.eq.3) iret=1
      go to (100,100,300,400,400,400,700,800,900),ityp
c
c
c
c
c        cut-off of dipole type
c        **********************
c
c
c
c
  100 c5=c(mm,im)
      c6=c(mm+1,im)
      nexp=ic(mi+2,im)
c
      do 105 i=1,nt
c
      aaa=c5/(c6-deltaq(i,iret))
c     -------------------------
c
      do 105 ii=1,nexp
  105 aa(i)=aa(i)*aaa
c
      if (iprspc.le.2) go to 120
c
      do 110 i=1,nt
c
      aaa=c5/(c6-deltaq(i,iprspc))
c     ----------------------------
c
      do 110 ii=1,nexp
  110 aa(i)=aa(i)*aaa
c
c
  120 mi=mi+3
      mm=mm+2
      go to 999
c
c
c
c
c        cut-off of regge type /schierholz/
c        **********************************
c
c
c
c
  300 ax=ex*c(mm,im)
      ay=ey*c(mm,im)
      expo=dlog((ax+dsqrt(ax*ax-1.d0))*(ay+dsqrt(ay*ay-1.d0)))
      expo1=c(mm+1,im)*expo
      expo2=c(mm+2,im)*expo
      do 305 i=1,nt
      expon=expo1+expo2*deltaq(i,iret)
      if (expon.lt.-50.d0) go to 302
c
      aa(i)=aa(i)*dexp(expon)
c     ---------------------
c
c
c     ---------------------
c
      go to 305
  302 aa(i)=0.d0
  305 continue
      mi=mi+2
      mm=mm+3
      go to 999
c
c
c
c
c        eikonal form factor
c        *******************
c
c
c
c
  400 ieik=ityp-3
c
      eikc5=c(mm,im)
      do 407 i=1,nt
      expon=0.d0
      ieik3=0
      go to (401,402,403),ieik
c
c        t-form
  401 d=-deltaq(i,iret)
      go to 404
c
c        u-form
  402 d=2.d0*xy2*ct(i)-deltaq(i,iret)
      go to 404
c
c        t-form * u-form
  403 ieik3=ieik3+1
      go to (401,402,405),ieik3
c
  404 d1=dsqrt(d)
      d2=dsqrt(4.d0+d)
c
      expon=eikc5*(2.d0+d)/(d1*d2)*dlog(0.5d0*(d1+d2))+expon
c     ------------------------------------------------
c
      go to (405,405,403),ieik
  405 if (expon.lt.-50.d0) go to 406
c
      cut(i)=dexp(expon)
c     ------------------
c
      go to 407
  406 cut(i)=0.d0
  407 continue
c
c        get or calculate normalization factor
c
      go to (411,412,412),ieik
c
c        normalization of t-form
  411 c6=c(mm+1,im)
      go to 490
c
c        normalization of u-form and t- * u-form
  412 ieik2=ic(mi+2,im)
      go to (421,422,423,421,422,423),ieik2
c
  421 c6=c(mm+2,im)
      go to 490
c
  422 if (smev.eq.ssmev) go to 442
      ssmev=smev
      do 441 il=1,ime
  441 indc(2,il)=.false.
  442 if (indc(2,im)) go to 421
      indc(2,im)=.true.
      ess=smev*dwn
      ess=ess*ess
      go to 450
c
  423 ess=4.d0*ex*ey
c
  450 d=-(4.d0-ess)
      if (ieik2.le.3) d=d+c(4,im)
      if (d.eq.0.d0) go to 460
      d2=dsqrt(4.d0+d)
      if (d.lt.0.d0) go to 470
      d1=dsqrt(d)
      c6=dexp(-eikc5*(2.d0+d)/(d1*d2)*dlog(0.5d0*(d1+d2)))
      go to 480
c
  460 c6=dexp(-eikc5*0.5d0)
      go to 480
c
  470 d1=dsqrt(-d)
      c6=dexp(-eikc5*(2.d0+d)/(d1*d2)*datan(d1/d2))
c
  480 if (ieik.eq.3) c6=c6*c(mm+1,im)
      if (ieik2.eq.3.or.ieik2.eq.6) go to 490
      c(mm+2,im)=c6
c
c        compute form factor
c
  490 do 495 i=1,nt
  495 aa(i)=aa(i)*cut(i)*c6
c     -------------------
c
      mi=mi+3
      mm=mm+3
      go to 999
c
c
c        exponential form factor
c        ***********************
c
c
  700 c5=c(mm,im)
      do 705 i=1,nt
c
c**** expo=(deltaq(i,iret)-c4)*c5
      expo=(deltaq(i,iret))*c5
c     ----------------------
c
c**** if (expo.lt.-50.d0) go to 704
      if (expo.lt.-50.d0) expo=-50.d0
c
      aa(i)=aa(i)*dexp(expo)
c     ----------------------
c
      go to 705
  704 aa(i)=0.d0
  705 continue
      mi=mi+2
      mm=mm+1
      go to 999
c
c
c        cloudy bag form factor
c        ***********************
c
c
  800 c5=c(mm,im)
      nexp=ic(mi+2,im)
      do 805 i=1,nt
c
      arg=dsqrt(-deltaq(i,iret))*c5
      argc=arg*arg*arg
      aaa=3.d0*(dsin(arg)-arg*dcos(arg))/argc
c
      do 805 ii=1,nexp
  805 aa(i)=aa(i)*aaa
c
      mi=mi+3
      mm=mm+1
      go to 999
c
c
c        propagator of mass-distributed meson
c        ************************************
c
c
  900 c5=c(mm,im)
      c6=c(mm+1,im)
      nspin=ic(mi+2,im)
      indp2=.false.
      if (iprspc.le.1) go to 901
      indp2=.true.
      iret=1
c
  901 do 915 i=1,nt
      d=c6-deltaq(i,iret)
      if (nspin.eq.0) go to 903
      d1=-d*deltaq(i,iret)/c4
  903 d=dsqrt(d)
      if (nspin.eq.0) go to 907
      do 905 ii=1,nspin
  905 d=d*d1
c
  907 omq=c4-deltaq(i,iret)+c5*d
      if (indp2) go to 910
      aa(i)=aa(i)/omq
      go to 915
c
  910 om=dsqrt(omq)
      aa(i)=aa(i)/(om*(om+ez1))
c
  915 continue
c
      mi=mi+3
      mm=mm+2
      go to 999
c
c
c
c
 2000 return
      end
      subroutine legp (pj,pjm1,x,j)
c
c
c        subroutine legp   computes the legendre polynoms
c
      real*8 pj,pjm1,x,a,b
c
c
c
c        compute legendre polynom for j equals zero
c
c
c
      pjm1=1.d0
      if (j.gt.0) go to 1
      pj=1.d0
      return
c
c
c
c        compute legendre polynoms for j equals one
c
c
c
    1 pj=x
      if (j.eq.1) return
c
c
c
c        compute legendre polynom for j greater or equal two
c
c
c
      do 2 i=2,j
      a=x*pj
      b=a-pjm1
      pjm1=pj
    2 pj=-b/dfloat(i)+b+a
c
c
      return
      end
      subroutine gset(ax,bx,n,z,w)
      implicit real*8 (a-h,o-z)
c
c     n-point gauss zeros and weights for the interval (ax,bx) are
c           stored in  arrays z and w respectively.
c
      dimension     a(273),x(273),ktab(96)
      dimension z(2),w(2)
c
c-----table of initial subscripts for n=2(1)16(4)96
      data ktab(2)/1/
      data ktab(3)/2/
      data ktab(4)/4/
      data ktab(5)/6/
      data ktab(6)/9/
      data ktab(7)/12/
      data ktab(8)/16/
      data ktab(9)/20/
      data ktab(10)/25/
      data ktab(11)/30/
      data ktab(12)/36/
      data ktab(13)/42/
      data ktab(14)/49/
      data ktab(15)/56/
      data ktab(16)/64/
      data ktab(20)/72/
      data ktab(24)/82/
      data ktab(28)/82/
      data ktab(32)/94/
      data ktab(36)/94/
      data ktab(40)/110/
      data ktab(44)/110/
      data ktab(48)/130/
      data ktab(52)/130/
      data ktab(56)/130/
      data ktab(60)/130/
      data ktab(64)/154/
      data ktab(68)/154/
      data ktab(72)/154/
      data ktab(76)/154/
      data ktab(80)/186/
      data ktab(84)/186/
      data ktab(88)/186/
      data ktab(92)/186/
      data ktab(96)/226/
c
c-----table of abscissae (x) and weights (a) for interval (-1,+1).
c
c**** n=2
      data x(1)/0.577350269189626  d0/, a(1)/1.000000000000000  d0/
c**** n=3
      data x(2)/0.774596669241483  d0/, a(2)/0.555555555555556  d0/
      data x(3)/0.000000000000000  d0/, a(3)/0.888888888888889  d0/
c**** n=4
      data x(4)/0.861136311594053  d0/, a(4)/0.347854845137454  d0/
      data x(5)/0.339981043584856  d0/, a(5)/0.652145154862546  d0/
c**** n=5
      data x(6)/0.906179845938664  d0/, a(6)/0.236926885056189  d0/
      data x(7)/0.538469310105683  d0/, a(7)/0.478628670499366  d0/
      data x(8)/0.000000000000000  d0/, a(8)/0.568888888888889  d0/
c**** n=6
      data x(9)/0.932469514203152  d0/, a(9)/0.171324492379170  d0/
      data x(10)/0.661209386466265 d0/, a(10)/0.360761573048139 d0/
      data x(11)/0.238619186083197 d0/, a(11)/0.467913934572691 d0/
c**** n=7
      data x(12)/0.949107912342759 d0/, a(12)/0.129484966168870 d0/
      data x(13)/0.741531185599394 d0/, a(13)/0.279705391489277 d0/
      data x(14)/0.405845151377397 d0/, a(14)/0.381830050505119 d0/
      data x(15)/0.000000000000000 d0/, a(15)/0.417959183673469 d0/
c**** n=8
      data x(16)/0.960289856497536 d0/, a(16)/0.101228536290376 d0/
      data x(17)/0.796666477413627 d0/, a(17)/0.222381034453374 d0/
      data x(18)/0.525532409916329 d0/, a(18)/0.313706645877887 d0/
      data x(19)/0.183434642495650 d0/, a(19)/0.362683783378362 d0/
c**** n=9
      data x(20)/0.968160239507626 d0/, a(20)/0.081274388361574 d0/
      data x(21)/0.836031107326636 d0/, a(21)/0.180648160694857 d0/
      data x(22)/0.613371432700590 d0/, a(22)/0.260610696402935 d0/
      data x(23)/0.324253423403809 d0/, a(23)/0.312347077040003 d0/
      data x(24)/0.000000000000000 d0/, a(24)/0.330239355001260 d0/
c**** n=10
      data x(25)/0.973906528517172 d0/, a(25)/0.066671344308688 d0/
      data x(26)/0.865063366688985 d0/, a(26)/0.149451349150581 d0/
      data x(27)/0.679409568299024 d0/, a(27)/0.219086362515982 d0/
      data x(28)/0.433395394129247 d0/, a(28)/0.269266719309996 d0/
      data x(29)/0.148874338981631 d0/, a(29)/0.295524224714753 d0/
c**** n=11
      data x(30)/0.978228658146057 d0/, a(30)/0.055668567116174 d0/
      data x(31)/0.887062599768095 d0/, a(31)/0.125580369464905 d0/
      data x(32)/0.730152005574049 d0/, a(32)/0.186290210927734 d0/
      data x(33)/0.519096129206812 d0/, a(33)/0.233193764591990 d0/
      data x(34)/0.269543155952345 d0/, a(34)/0.262804544510247 d0/
      data x(35)/0.000000000000000 d0/, a(35)/0.272925086777901 d0/
c**** n=12
      data x(36)/0.981560634246719 d0/, a(36)/0.047175336386512 d0/
      data x(37)/0.904117256370475 d0/, a(37)/0.106939325995318 d0/
      data x(38)/0.769902674194305 d0/, a(38)/0.160078328543346 d0/
      data x(39)/0.587317954286617 d0/, a(39)/0.203167426723066 d0/
      data x(40)/0.367831498998180 d0/, a(40)/0.233492536538355 d0/
      data x(41)/0.125233408511469 d0/, a(41)/0.249147045813403 d0/
c**** n=13
      data x(42)/0.984183054718588 d0/, a(42)/0.040484004765316 d0/
      data x(43)/0.917598399222978 d0/, a(43)/0.092121499837728 d0/
      data x(44)/0.801578090733310 d0/, a(44)/0.138873510219787 d0/
      data x(45)/0.642349339440340 d0/, a(45)/0.178145980761946 d0/
      data x(46)/0.448492751036447 d0/, a(46)/0.207816047536889 d0/
      data x(47)/0.230458315955135 d0/, a(47)/0.226283180262897 d0/
      data x(48)/0.000000000000000 d0/, a(48)/0.232551553230874 d0/
c**** n=14
      data x(49)/0.986283808696812 d0/, a(49)/0.035119460331752 d0/
      data x(50)/0.928434883663574 d0/, a(50)/0.080158087159760 d0/
      data x(51)/0.827201315069765 d0/, a(51)/0.121518570687903 d0/
      data x(52)/0.687292904811685 d0/, a(52)/0.157203167158194 d0/
      data x(53)/0.515248636358154 d0/, a(53)/0.185538397477938 d0/
      data x(54)/0.319112368927890 d0/, a(54)/0.205198463721296 d0/
      data x(55)/0.108054948707344 d0/, a(55)/0.215263853463158 d0/
c**** n=15
      data x(56)/0.987992518020485 d0/, a(56)/0.030753241996117 d0/
      data x(57)/0.937273392400706 d0/, a(57)/0.070366047488108 d0/
      data x(58)/0.848206583410427 d0/, a(58)/0.107159220467172 d0/
      data x(59)/0.724417731360170 d0/, a(59)/0.139570677926154 d0/
      data x(60)/0.570972172608539 d0/, a(60)/0.166269205816994 d0/
      data x(61)/0.394151347077563 d0/, a(61)/0.186161000015562 d0/
      data x(62)/0.201194093997435 d0/, a(62)/0.198431485327111 d0/
      data x(63)/0.000000000000000 d0/, a(63)/0.202578241925561 d0/
c**** n=16
      data x(64)/0.989400934991650 d0/, a(64)/0.027152459411754 d0/
      data x(65)/0.944575023073233 d0/, a(65)/0.062253523938648 d0/
      data x(66)/0.865631202387832 d0/, a(66)/0.095158511682493 d0/
      data x(67)/0.755404408355003 d0/, a(67)/0.124628971255534 d0/
      data x(68)/0.617876244402644 d0/, a(68)/0.149595988816577 d0/
      data x(69)/0.458016777657227 d0/, a(69)/0.169156519395003 d0/
      data x(70)/0.281603550779259 d0/, a(70)/0.182603415044924 d0/
      data x(71)/0.095012509837637 d0/, a(71)/0.189450610455069 d0/
c**** n=20
      data x(72)/0.993128599185094 d0/, a(72)/0.017614007139152 d0/
      data x(73)/0.963971927277913 d0/, a(73)/0.040601429800386 d0/
      data x(74)/0.912234428251325 d0/, a(74)/0.062672048334109 d0/
      data x(75)/0.839116971822218 d0/, a(75)/0.083276741576704 d0/
      data x(76)/0.746331906460150 d0/, a(76)/0.101930119817240 d0/
      data x(77)/0.636053680726515 d0/, a(77)/0.118194531961518 d0/
      data x(78)/0.510867001950827 d0/, a(78)/0.131688638449176 d0/
      data x(79)/0.373706088715419 d0/, a(79)/0.142096109318382 d0/
      data x(80)/0.227785851141645 d0/, a(80)/0.149172986472603 d0/
      data x(81)/0.076526521133497 d0/, a(81)/0.152753387130725 d0/
c**** n=24
      data x(82)/0.995187219997021 d0/, a(82)/0.012341229799987 d0/
      data x(83)/0.974728555971309 d0/, a(83)/0.028531388628933 d0/
      data x(84)/0.938274552002732 d0/, a(84)/0.044277438817419 d0/
      data x(85)/0.886415527004401 d0/, a(85)/0.059298584915436 d0/
      data x(86)/0.820001985973902 d0/, a(86)/0.073346481411080 d0/
      data x(87)/0.740124191578554 d0/, a(87)/0.086190161531953 d0/
      data x(88)/0.648093651936975 d0/, a(88)/0.097618652104113 d0/
      data x(89)/0.545421471388839 d0/, a(89)/0.107444270115965 d0/
      data x(90)/0.433793507626045 d0/, a(90)/0.115505668053725 d0/
      data x(91)/0.315042679696163 d0/, a(91)/0.121670472927803 d0/
      data x(92)/0.191118867473616 d0/, a(92)/0.125837456346828 d0/
      data x(93)/0.064056892862605 d0/, a(93)/0.127938195346752 d0/
c**** n=32
      data x(94)/0.997263861849481 d0/, a(94)/0.007018610009470 d0/
      data x(95)/0.985611511545268 d0/, a(95)/0.016274394730905 d0/
      data x(96)/0.964762255587506 d0/, a(96)/0.025392065309262 d0/
      data x(97)/0.934906075937739 d0/, a(97)/0.034273862913021 d0/
      data x(98)/0.896321155766052 d0/, a(98)/0.042835898022226 d0/
      data x(99)/0.849367613732569 d0/, a(99)/0.050998059262376 d0/
      data x(100)/0.794483795967942d0/, a(100)/0.058684093478535d0/
      data x(101)/0.732182118740289d0/, a(101)/0.065822222776361d0/
      data x(102)/0.663044266930215d0/, a(102)/0.072345794108848d0/
      data x(103)/0.587715757240762d0/, a(103)/0.078193895787070d0/
      data x(104)/0.506899908932229d0/, a(104)/0.083311924226946d0/
      data x(105)/0.421351276130635d0/, a(105)/0.087652093004403d0/
      data x(106)/0.331868602282127d0/, a(106)/0.091173878695763d0/
      data x(107)/0.239287362252137d0/, a(107)/0.093844399080804d0/
      data x(108)/0.144471961582796d0/, a(108)/0.095638720079274d0/
      data x(109)/0.048307665687738d0/, a(109)/0.096540088514727d0/
c**** n=40
      data x(110)/0.998237709710559d0/, a(110)/0.004521277098533d0/
      data x(111)/0.990726238699457d0/, a(111)/0.010498284531152d0/
      data x(112)/0.977259949983774d0/, a(112)/0.016421058381907d0/
      data x(113)/0.957916819213791d0/, a(113)/0.022245849194166d0/
      data x(114)/0.932812808278676d0/, a(114)/0.027937006980023d0/
      data x(115)/0.902098806968874d0/, a(115)/0.033460195282547d0/
      data x(116)/0.865959503212259d0/, a(116)/0.038782167974472d0/
      data x(117)/0.824612230833311d0/, a(117)/0.043870908185673d0/
      data x(118)/0.778305651426519d0/, a(118)/0.048695807635072d0/
      data x(119)/0.727318255189927d0/, a(119)/0.053227846983936d0/
      data x(120)/0.671956684614179d0/, a(120)/0.057439769099391d0/
      data x(121)/0.612553889667980d0/, a(121)/0.061306242492928d0/
      data x(122)/0.549467125095128d0/, a(122)/0.064804013456601d0/
      data x(123)/0.483075801686178d0/, a(123)/0.067912045815233d0/
      data x(124)/0.413779204371605d0/, a(124)/0.070611647391286d0/
      data x(125)/0.341994090825758d0/, a(125)/0.072886582395804d0/
      data x(126)/0.268152185007253d0/, a(126)/0.074723169057968d0/
      data x(127)/0.192697580701371d0/, a(127)/0.076110361900626d0/
      data x(128)/0.116084070675255d0/, a(128)/0.077039818164247d0/
      data x(129)/0.038772417506050d0/, a(129)/0.077505947978424d0/
c**** n=48
      data x(130)/0.998771007252426d0/, a(130)/0.003153346052305d0/
      data x(131)/0.993530172266350d0/, a(131)/0.007327553901276d0/
      data x(132)/0.984124583722826d0/, a(132)/0.011477234579234d0/
      data x(133)/0.970591592546247d0/, a(133)/0.015579315722943d0/
      data x(134)/0.952987703160430d0/, a(134)/0.019616160457355d0/
      data x(135)/0.931386690706554d0/, a(135)/0.023570760839324d0/
      data x(136)/0.905879136715569d0/, a(136)/0.027426509708356d0/
      data x(137)/0.876572020274247d0/, a(137)/0.031167227832798d0/
      data x(138)/0.843588261624393d0/, a(138)/0.034777222564770d0/
      data x(139)/0.807066204029442d0/, a(139)/0.038241351065830d0/
      data x(140)/0.767159032515740d0/, a(140)/0.041545082943464d0/
      data x(141)/0.724034130923814d0/, a(141)/0.044674560856694d0/
      data x(142)/0.677872379632663d0/, a(142)/0.047616658492490d0/
      data x(143)/0.628867396776513d0/, a(143)/0.050359035553854d0/
      data x(144)/0.577224726083972d0/, a(144)/0.052890189485193d0/
      data x(145)/0.523160974722233d0/, a(145)/0.055199503699984d0/
      data x(146)/0.466902904750958d0/, a(146)/0.057277292100403d0/
      data x(147)/0.408686481990716d0/, a(147)/0.059114839698395d0/
      data x(148)/0.348755886292160d0/, a(148)/0.060704439165893d0/
      data x(149)/0.287362487355455d0/, a(149)/0.062039423159892d0/
      data x(150)/0.224763790394689d0/, a(150)/0.063114192286254d0/
      data x(151)/0.161222356068891d0/, a(151)/0.063924238584648d0/
      data x(152)/0.097004699209462d0/, a(152)/0.064466164435950d0/
      data x(153)/0.032380170962869d0/, a(153)/0.064737696812683d0/
c**** n=64
      data x(154)/0.999305041735772d0/, a(154)/0.001783280721696d0/
      data x(155)/0.996340116771955d0/, a(155)/0.004147033260562d0/
      data x(156)/0.991013371476744d0/, a(156)/0.006504457968978d0/
      data x(157)/0.983336253884625d0/, a(157)/0.008846759826363d0/
      data x(158)/0.973326827789910d0/, a(158)/0.011168139460131d0/
      data x(159)/0.961008799652053d0/, a(159)/0.013463047896718d0/
      data x(160)/0.946411374858402d0/, a(160)/0.015726030476024d0/
      data x(161)/0.929569172131939d0/, a(161)/0.017951715775697d0/
      data x(162)/0.910522137078502d0/, a(162)/0.020134823153530d0/
      data x(163)/0.889315445995114d0/, a(163)/0.022270173808383d0/
      data x(164)/0.865999398154092d0/, a(164)/0.024352702568710d0/
      data x(165)/0.840629296252580d0/, a(165)/0.026377469715054d0/
      data x(166)/0.813265315122797d0/, a(166)/0.028339672614259d0/
      data x(167)/0.783972358943341d0/, a(167)/0.030234657072402d0/
      data x(168)/0.752819907260531d0/, a(168)/0.032057928354851d0/
      data x(169)/0.719881850171610d0/, a(169)/0.033805161837141d0/
      data x(170)/0.685236313054233d0/, a(170)/0.035472213256882d0/
      data x(171)/0.648965471254657d0/, a(171)/0.037055128540240d0/
      data x(172)/0.611155355172393d0/, a(172)/0.038550153178615d0/
      data x(173)/0.571895646202634d0/, a(173)/0.039953741132720d0/
      data x(174)/0.531279464019894d0/, a(174)/0.041262563242623d0/
      data x(175)/0.489403145707052d0/, a(175)/0.042473515123653d0/
      data x(176)/0.446366017253464d0/, a(176)/0.043583724529323d0/
      data x(177)/0.402270157963991d0/, a(177)/0.044590558163756d0/
      data x(178)/0.357220158337668d0/, a(178)/0.045491627927418d0/
      data x(179)/0.311322871990210d0/, a(179)/0.046284796581314d0/
      data x(180)/0.264687162208767d0/, a(180)/0.046968182816210d0/
      data x(181)/0.217423643740007d0/, a(181)/0.047540165714830d0/
      data x(182)/0.169644420423992d0/, a(182)/0.047999388596458d0/
      data x(183)/0.121462819296120d0/, a(183)/0.048344762234802d0/
      data x(184)/0.072993121787799d0/, a(184)/0.048575467441503d0/
      data x(185)/0.024350292663424d0/, a(185)/0.048690957009139d0/
c**** n=80
      data x(186)/0.999553822651630d0/, a(186)/0.001144950003186d0/
      data x(187)/0.997649864398237d0/, a(187)/0.002663533589512d0/
      data x(188)/0.994227540965688d0/, a(188)/0.004180313124694d0/
      data x(189)/0.989291302499755d0/, a(189)/0.005690922451403d0/
      data x(190)/0.982848572738629d0/, a(190)/0.007192904768117d0/
      data x(191)/0.974909140585727d0/, a(191)/0.008683945269260d0/
      data x(192)/0.965485089043799d0/, a(192)/0.010161766041103d0/
      data x(193)/0.954590766343634d0/, a(193)/0.011624114120797d0/
      data x(194)/0.942242761309872d0/, a(194)/0.013068761592401d0/
      data x(195)/0.928459877172445d0/, a(195)/0.014493508040509d0/
      data x(196)/0.913263102571757d0/, a(196)/0.015896183583725d0/
      data x(197)/0.896675579438770d0/, a(197)/0.017274652056269d0/
      data x(198)/0.878722567678213d0/, a(198)/0.018626814208299d0/
      data x(199)/0.859431406663111d0/, a(199)/0.019950610878141d0/
      data x(200)/0.838831473580255d0/, a(200)/0.021244026115782d0/
      data x(201)/0.816954138681463d0/, a(201)/0.022505090246332d0/
      data x(202)/0.793832717504605d0/, a(202)/0.023731882865930d0/
      data x(203)/0.769502420135041d0/, a(203)/0.024922535764115d0/
      data x(204)/0.744000297583597d0/, a(204)/0.026075235767565d0/
      data x(205)/0.717365185362099d0/, a(205)/0.027188227500486d0/
      data x(206)/0.689637644342027d0/, a(206)/0.028259816057276d0/
      data x(207)/0.660859898986119d0/, a(207)/0.029288369583267d0/
      data x(208)/0.631075773046871d0/, a(208)/0.030272321759557d0/
      data x(209)/0.600330622829751d0/, a(209)/0.031210174188114d0/
      data x(210)/0.568671268122709d0/, a(210)/0.032100498673487d0/
      data x(211)/0.536145920897131d0/, a(211)/0.032941939397645d0/
      data x(212)/0.502804111888784d0/, a(212)/0.033733214984611d0/
      data x(213)/0.468696615170544d0/, a(213)/0.034473120451753d0/
      data x(214)/0.433875370831756d0/, a(214)/0.035160529044747d0/
      data x(215)/0.398393405881969d0/, a(215)/0.035794393953416d0/
      data x(216)/0.362304753499487d0/, a(216)/0.036373749905835d0/
      data x(217)/0.325664370747701d0/, a(217)/0.036897714638276d0/
      data x(218)/0.288528054884511d0/, a(218)/0.037365490238730d0/
      data x(219)/0.250952358392272d0/, a(219)/0.037776364362001d0/
      data x(220)/0.212994502857666d0/, a(220)/0.038129711314477d0/
      data x(221)/0.174712291832646d0/, a(221)/0.038424993006959d0/
      data x(222)/0.136164022809143d0/, a(222)/0.038661759774076d0/
      data x(223)/0.097408398441584d0/, a(223)/0.038839651059051d0/
      data x(224)/0.058504437152420d0/, a(224)/0.038958395962769d0/
      data x(225)/0.019511383256793d0/, a(225)/0.039017813656306d0/
c**** n=96
      data x(226)/0.999689503883230d0/, a(226)/0.000796792065552d0/
      data x(227)/0.998364375863181d0/, a(227)/0.001853960788946d0/
      data x(228)/0.995981842987209d0/, a(228)/0.002910731817934d0/
      data x(229)/0.992543900323762d0/, a(229)/0.003964554338444d0/
      data x(230)/0.988054126329623d0/, a(230)/0.005014202742927d0/
      data x(231)/0.982517263563014d0/, a(231)/0.006058545504235d0/
      data x(232)/0.975939174585136d0/, a(232)/0.007096470791153d0/
      data x(233)/0.968326828463264d0/, a(233)/0.008126876925698d0/
      data x(234)/0.959688291448742d0/, a(234)/0.009148671230783d0/
      data x(235)/0.950032717784437d0/, a(235)/0.010160770535008d0/
      data x(236)/0.939370339752755d0/, a(236)/0.011162102099838d0/
      data x(237)/0.927712456722308d0/, a(237)/0.012151604671088d0/
      data x(238)/0.915071423120898d0/, a(238)/0.013128229566961d0/
      data x(239)/0.901460635315852d0/, a(239)/0.014090941772314d0/
      data x(240)/0.886894517402420d0/, a(240)/0.015038721026994d0/
      data x(241)/0.871388505909296d0/, a(241)/0.015970562902562d0/
      data x(242)/0.854959033434601d0/, a(242)/0.016885479864245d0/
      data x(243)/0.837623511228187d0/, a(243)/0.017782502316045d0/
      data x(244)/0.819400310737931d0/, a(244)/0.018660679627411d0/
      data x(245)/0.800308744139140d0/, a(245)/0.019519081140145d0/
      data x(246)/0.780369043867433d0/, a(246)/0.020356797154333d0/
      data x(247)/0.759602341176647d0/, a(247)/0.021172939892191d0/
      data x(248)/0.738030643744400d0/, a(248)/0.021966644438744d0/
      data x(249)/0.715676812348967d0/, a(249)/0.022737069658329d0/
      data x(250)/0.692564536642171d0/, a(250)/0.023483399085926d0/
      data x(251)/0.668718310043916d0/, a(251)/0.024204841792364d0/
      data x(252)/0.644163403784967d0/, a(252)/0.024900633222483d0/
      data x(253)/0.618925840125468d0/, a(253)/0.025570036005349d0/
      data x(254)/0.593032364777572d0/, a(254)/0.026212340735672d0/
      data x(255)/0.566510418561397d0/, a(255)/0.026826866725591d0/
      data x(256)/0.539388108324357d0/, a(256)/0.027412962726029d0/
      data x(257)/0.511694177154667d0/, a(257)/0.027970007616848d0/
      data x(258)/0.483457973920596d0/, a(258)/0.028497411065085d0/
      data x(259)/0.454709422167743d0/, a(259)/0.028994614150555d0/
      data x(260)/0.425478988407300d0/, a(260)/0.029461089958167d0/
      data x(261)/0.395797649828908d0/, a(261)/0.029896344136328d0/
      data x(262)/0.365696861472313d0/, a(262)/0.030299915420827d0/
      data x(263)/0.335208522892625d0/, a(263)/0.030671376123669d0/
      data x(264)/0.304364944354496d0/, a(264)/0.031010332586313d0/
      data x(265)/0.273198812591049d0/, a(265)/0.031316425596861d0/
      data x(266)/0.241743156163840d0/, a(266)/0.031589330770727d0/
      data x(267)/0.210031310460567d0/, a(267)/0.031828758894411d0/
      data x(268)/0.178096882367618d0/, a(268)/0.032034456231992d0/
      data x(269)/0.145973714654896d0/, a(269)/0.032206204794030d0/
      data x(270)/0.113695850110665d0/, a(270)/0.032343822568575d0/
      data x(271)/0.081297495464425d0/, a(271)/0.032447163714064d0/
      data x(272)/0.048812985136049d0/, a(272)/0.032516118713868d0/
      data x(273)/0.016276744849602d0/, a(273)/0.032550614492363d0/
c
c
c-----test n
      alpha=0.5d0*(ax+bx)
      beta=0.5d0*(bx-ax)
      if( n.lt.1 .or. n.gt.96 ) go to 100
      if(n.ne.1) go to 1
      z(1)=alpha
      w(1)=bx-ax
      return
c
    1 if (n.le.16) go to 3
      if (n.gt.24) go to 4
      n=4*(n/4)
      go to 3
    4 if (n.gt.48) go to 5
      n=8*(n/8)
      go to 3
    5 n=16*(n/16)
c
c----- set k equal to initial subscript and store results
    3 k=ktab(n)
      m=n/2
      do 2 j=1,m
      jtab=k-1+j
      wtemp=beta*a(jtab)
      delta=beta*x(jtab)
      z(j)=alpha-delta
      w(j)=wtemp
      jp=n+1-j
      z(jp)=alpha+delta
      w(jp)=wtemp
    2 continue
      if((n-m-m).eq.0) return
      z(m+1)=alpha
      jmid=k+m
      w(m+1)=beta*a(jmid)
      return
c
  100 zn=n
      write(6,200) zn
  200 format(1h0/////'0error in gset. n has the non-permissible value',
     1e11.3/'0execution terminated.')
      stop
      end




























      
      function spe (qq,qfq,ua,ws,uc,ud,wn,iprop,ispex,ipaho)
c
c        single particle energy of one particle below or above
c        the fermi surface of a medium
c
c        spe is called by:
c                          matbnd,
c                          matusf,
c                          matgmt,
c                          obaa,
c                          tbibnn,
c                          tbibd,
c                          tbaann,
c                          tbainn.
c
c        ispe = 0 : same as ispe=2
c        ispe = 1 : same as ispe=2
c        ispe = 2 : continous spectrum with free energies below and
c                   above fermi surface.
c        ispe = 3 : same as ispe=2 plus constant shift of particle
c                   potential above fermi surface by uc.
c        ispe = 4 : conventional choice for sp energy:
c                   bound energy below the fermi surface and
c                   free energy above the fermi surface (i.e. gap).
c        ispe = 5 : continous spectrum with bound energy below the ferm
c                   surface and continous continuation above the fermi
c                   surface until u=0., free energies after
c        ispe = 6 : continuous for ever
c        ispe = 7 : continuous choice with k-dependent parameters
c                   according to fua in case of iprop=2;
c                   in case of iprop=1 same as ispe=6.
c
c
      implicit real*8 (a-h,o-z)
c
      data wwn/938.919d0/
c
c
c
c
      ispe=ispex
      if (ispe.lt.2) ispe=2
c
c
c
c
      go to (1000,2000),iprop
c
c
c
c
 1000 if (ipaho.eq.0) go to 1001
      go to (1100,1200),ipaho
 1001 if (qq.gt.qfq) go to 1200
c
c
 1100 go to (1110,1110,1110,1140,1140,1140,1170),ispe
 1110 spe=0.5d0*qq/wn+wn-wwn
      go to 9000

 1140 spe=0.5d0*qq/ws+ua
      if (wn.ne.wwn) spe=spe+0.5d0*qq*(1.d0/wn-1.d0/wwn)+wn-wwn
      go to 9000
 1170 wsk=wn+fua(qq,qfq,wn)
      spe=0.5d0*qq/wn+wsk+fub(qq,qfq,wn)
      go to 9000
c
c
 1200 go to (1210,1210,1210,1210,1250,1140,1170),ispe
 1210 spe=0.5d0*qq/wn+wn-wwn
      go to 8000
 1250 spe1=0.5d0*qq/ws+ua
      if (wn.ne.wwn) spe1=spe1+0.5d0*qq*(1.d0/wn-1.d0/wwn)+wn-wwn
      spe2=0.5d0*qq/wn+wn-wwn
      spe=dmin1(spe1,spe2)
      go to 9000
c
c
c
c
 2000 if (ipaho.eq.0) go to 2001
      go to (2100,2200),ipaho
 2001 if (qq.gt.qfq) go to 2200
c
c
 2100 go to (2110,2110,2110,2140,2140,2140,2270),ispe
 2110 spe=dsqrt(qq+wn*wn)
      go to 9000
 2140 spe=dsqrt(qq+ws*ws)+wwn-ws+ua
      if (wn.ne.wwn) spe=spe+dsqrt(qq+wn*wn)-dsqrt(qq+wwn*wwn)
      go to 9000
c
c
 2200 go to (2210,2210,2210,2210,2250,2140,2270),ispe
 2210 spe=dsqrt(qq+wn*wn)
      go to 8000
 2250 spe1=dsqrt(qq+ws*ws)+wwn-ws+ua
      if (wn.ne.wwn) spe1=spe1+dsqrt(qq+wn*wn)-dsqrt(qq+wwn*wwn)
      spe2=dsqrt(qq+wn*wn)
      spe=dmin1(spe1,spe2)
      go to 9000
 2270 wsk=wn+fua(qq,qfq,wn)
      spe=dsqrt(qq+wsk*wsk)+fub(qq,qfq,wn)
c****      if (qq.lt.25.*qfq) go to 9000
c****      spe=dsqrt(qq+wn*wn)
      go to 9000
c
c
c
c
 8000 if (mod(ispe,2).eq.1) spe=spe+uc
c
c
c
c
 9000 return
      end
      function derspe (q,qq,qfq,ws,wn,iprop,ispex,ipaho)
c
c        derivative of single particle energy of one particle below or a
c        the fermi surface of a medium
c
c        note: qq is not just the square of q in the argument list,
c              qq=pmevq+q*q
c
c        derspe is needed for  the computation of the principal value
c        integral in the brueckner equation if the single particle poten
c        is also evaluated above the fermi surface.
c
c
c
c
c        ispe = 0 : same as ispe=2
c        ispe = 1 : same as ispe=2
c        ispe = 2 : continous spectrum with free energies below and
c                   above fermi surface.
c        ispe = 3 : same as ispe=2 plus constant shift of particle
c                   potential above fermi surface by uc.
c        ispe = 4 : conventional choice for sp energy:
c                   bound energy below the fermi surface and
c                   free energy above the fermi surface (i.e. gap).
c        ispe = 5 : continous spectrum with bound energy below the ferm
c                   surface and continous continuation above the fermi
c                   surface until u=0., free energies after
c        ispe = 6 : continuous for ever
c        ispe = 7 : continuous choice with k-dependent parameters
c                   according to fua in case of iprop=2;
c                   in case of iprop=1 same as ispe=6.
c
c
      implicit real*8 (a-h,o-z)
c
      data wwn/938.926d0/
c
c
c
c
      ispe=ispex
      if (ispe.lt.2) ispe=2
c
c
c
c
      go to (1000,2000),iprop
c
c
c
c
 1000 if (ipaho.eq.0) go to 1001
      go to (1100,1200),ipaho
 1001 if (qq.gt.qfq) go to 1200
c
c
 1100 go to (1110,1110,1110,1140,1140,1140,1140),ispe
 1110 derspe=q/wn
      go to 9000
 1140 derspe=q/ws
      if (wn.ne.wwn) derspe=derspe+q*(1.d0/wn-1.d0/wwn)
      go to 9000
c
c
 1200 go to (1210,1210,1210,1210,1250,1140,1140),ispe
 1210 derspe=q/wn
      go to 8000
 1250 spe1=q/ws
      if (wn.ne.wwn) spe1=spe1+q*(1.d0/wn-1.d0/wwn)
      spe2=q/wn
      derspe=dmin1(spe1,spe2)
      go to 9000
c
c
c
c
 2000 if (ipaho.eq.0) go to 2001
      go to (2100,2200),ipaho
 2001 if (qq.gt.qfq) go to 2200
c
c
 2100 go to (2110,2110,2110,2140,2140,2140,2270),ispe
 2110 derspe=q/dsqrt(qq+wn*wn)
      go to 9000
 2140 derspe=q/dsqrt(qq+ws*ws)
      if (wn.ne.wwn) derspe=derspe+q/dsqrt(qq+wn*wn)-q/dsqrt(qq+wwn*wwn)
      go to 9000
c
c
 2200 go to (2210,2210,2210,2210,2250,2140,2270),ispe
 2210 derspe=q/dsqrt(qq+wn*wn)
      go to 8000
 2250 spe1=q/dsqrt(qq+ws*ws)
      if (wn.ne.wwn) spe1=spe1+q/dsqrt(qq+wn*wn)-q/dsqrt(qq+wwn*wwn)
      spe2=q/dsqrt(qq+wn*wn)
      derspe=dmin1(spe1,spe2)
      go to 9000
 2270 wsk=wn+fua(qq,qfq,wn)
      derspe=(q+wsk*derfua(q,qq,qfq,wn))/dsqrt(qq+wsk*wsk)
     1 +derfub(q,qq,qfq,wn)
      go to 9000
c
c
c
c
 8000 continue
c
c
c
c
 9000 return
      end
      function fua (qq,qfq,wn)
c
c        functional form of a;
c        to be applied in case ispe = 7 only.
c
c
      implicit real*8 (a-h,o-z)
c
c
      common /cfuab/ a(5),b(5),c(5)
c
c
c
c
      qqa=qq/qfq
      aa2=(a(2)/a(1))
      fua=a(1)+a(2)*qqa
c
      if (aa2.gt.0.d0) go to 200
c
      if (fua.gt.0.d0) fua=0.d0
      go to 1000
c
c
  200 aa3=aa2/81.d0
      fua=fua*dexp(-aa3*qqa*qqa)
c
c
c****      ak=a(1)+a(2)*qqa
c****      ck=c(1)+c(2)*qqa
c
c****      fua=(ak-ck)/(wn+ck)*wn
c
c
c
c
c
 1000 return
      end
      function fub (qq,qfq,wn)
c
c        functional form of b;
c        to be applied in case ispe = 7 only.
c
c
      implicit real*8 (a-h,o-z)
c
      common /cfuab/ a(5),b(5),c(5)
c
c
c
c
      qqa=qq/qfq
      bb2=(b(2)/b(1))
      fub=b(1)+b(2)*qqa
c
      if (bb2.gt.0.d0) go to 200
c
      if (fub.lt.0.d0) fub=0.d0
      go to 1000
c
c
  200 bb3=bb2/81.d0
      fub=fub*dexp(-bb3*qqa*qqa)
c
c
c****      ak=a(1)+a(2)*qqa
c****      bk=b(1)+b(2)*qqa
c****      ck=c(1)+c(2)*qqa
c
c****      cpk=1.d0+ck/wn
c****      wsk=wn+ak
c
c****      fub=(bk*wn+(dsqrt(qq*cpk*cpk+wsk*wsk)+bk)*ck)/(wn+ck)
c
c
c
c
c
 1000 return
      end
      function derfua (q,qq,qfq,wn)
c
c        derivative of functional form of a;
c        to be applied in case ispe = 7 only.
c
c
      implicit real*8 (a-h,o-z)
c
      common /cfuab/ a(5),b(5),c(5)
c
c
c
c
      qa=q/qfq
      qqa=qq/qfq
      aa2=(a(2)/a(1))
      fua=a(1)+a(2)*qqa
      derfua=2.d0*a(2)*qa
c
      if (aa2.gt.0.d0) go to 200
c
      if (fua.gt.0.d0) derfua=0.d0
      go to 1000
c
c
  200 aa3=aa2/81.d0
      dexpo=dexp(-aa3*qqa*qqa)
      derfua=derfua*dexpo+fua*dexpo*(-aa3*qa*qqa*4.d0)
c
c
c****      ak=a(1)+a(2)*qqa
c****      ck=c(1)+c(2)*qqa
c
c****      wnpck=wn+ck
c
c****      derak=2.d0*a(2)*qa
c****      derck=2.d0*c(2)*qa
c
c****      derfua=wn*((derak-derck)*wnpck-(ak-ck)*derck)/
c****     1 (wnpck*wnpck)
c
c
c
c
 1000 return
      end
      function derfub (q,qq,qfq,wn)
c
c        derivative of functional form of b;
c        to be applied in case ispe = 7 only.
c
c
      implicit real*8 (a-h,o-z)
c
      common /cfuab/ a(5),b(5),c(5)
c
c
c
c
      qa=q/qfq
      qqa=qq/qfq
      bb2=(b(2)/b(1))
      fub=b(1)+b(2)*qqa
      derfub=2.d0*b(2)*qa
c
      if (bb2.gt.0.d0) go to 200
c
c
      if (fub.lt.0.d0) derfub=0.d0
      go to 1000
c
c
  200 bb3=bb2/81.d0
      dexpo=dexp(-bb3*qqa*qqa)
      derfub=derfub*dexpo+fub*dexpo*(-bb3*qa*qqa*4.d0)
c
c
c****      ak=a(1)+a(2)*qqa
c****      bk=b(1)+b(2)*qqa
c****      ck=c(1)+c(2)*qqa
c
c****      wnpck=wn+ck
c
c****      derak=2.d0*a(2)*qa
c****      derbk=2.d0*b(2)*qa
c****      derck=2.d0*c(2)*qa
c
c****      cpk=1.d0+ck/wn
c****      wsk=wn+ak
c
c****      root=dsqrt(qq*cpk*cpk+wsk*wsk)
c****      bra=root+bk
c****      zah=wn*bk+bra*ck
c
c****      derfub=((wn*derbk+derck*bra+((q*cpk*cpk+qq*cpk*derck/wn
c****     1 +wsk*derak)/root+derbk)*ck)*wnpck-zah*derck)/(wnpck*wnpck)
c
c
c
c
 1000 return
      end
      function smp (qq,qfmev,ismp)
      implicit real*8 (a-h,o-z)
 1000 smp=0.d0
 2000 return
      end
      subroutine obnd                                                   00000010
c                                                                       00000020
c                                                                       00000030
c        one-boson-exchange nn-nd interaction;                          00000040
c        version which uses numerical integration                       00000050
c                                                                       00000060
c        author:      r. machleidt                                      00000070
c                     institut fuer theoretische kernphysik bonn        00000080
c                     nussallee 14-16                                   00000090
c                     d - 5300  bonn, w. germany                        00000100
c                                                                       00000110
c                                                                       00000120
      implicit real*8 (a-h,o-z)                                         00000130
c                                                                       00000140
c                                                                       00000150
      common /crdwrt/ kread,kwrite,kpunch,kda(9)                        00000160
c                                                                       00000170
c                                                                       00000180
c        arguments and values of this subroutine                        00000190
c                                                                       00000200
      common /cpot/   v(6),xmev,ymev                                    00000210
      common /cstate/ j,heform,sing,trip,coup,endep,label               00000220
c                                                                       00000230
c                                                                       00000240
c        this has been the end of the common-blocks containing          00000250
c        the arguments and values of this subroutine in the case of     00000260
c        no energy-dependence of the potential;                         00000270
c        in case of energy-dependence look for the common-block /cped/  00000280
c        in obai and obaa.                                              00000290
c                                                                       00000300
c        specifications for these two common blocks                     00000310
c                                                                       00000320
      logical heform,sing,trip,coup,endep                               00000330
c                                                                       00000340
c                                                                       00000350
c        common block for all ob-subroutines                            00000360
c                                                                       00000370
      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),  00000380
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,    00000390
     2                ez1,ez2,ct(96),wt(96),                            00000400
     3                ic(10,15),ift(3),mint(3),maxt(3),nt,              00000410
     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),             00000420
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,    00000430
     6                indc(2,15),indpar(3),indxy                        00000440
c                                                                       00000450
c         specifications for this common block                          00000460
c                                                                       00000470
      logical indc,indxy,indpar                                         00000480
c                                                                       00000490
c                                                                       00000500
c        further specifications                                         00000510
c                                                                       00000520
      logical indmg(12),index,indfa                                     00000530
      character*4 mesong(12)
      dimension fa(15)                                                  00000550
      data pi/3.141592653589793d0/                                      00000560
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',                   00000570
     1                   '1-  ','1-t ','1-tt','1-st','1-ss',            00000580
     2                    '1+  ','2+  '/                                00000590
      data indmg/12*.false./,index/.false./                             00000600
c
c
c                                                                       00000630
c                                                                       00000640
      inter=2                                                           00000650
c                                                                       00000660
c                                                                       00000670
c                                                                       00000680
c                                                                       00000690
c        call subroutine obpar once and only once                       00000700
c                                                                       00000710
c                                                                       00000720
      if (index) go to 50                                               00000730
      index=.true.                                                      00000740
      if (indpar(inter)) go to 40                                       00000750
c                                                                       00000760
c                                                                       00000770
      call obpar                                                        00000780
c                                                                       00000790
c                                                                       00000800
   40 iftgo=ift(inter)+1                                                00000810
      wn=wnn(inter)                                                     00000820
      dwn=1.d0/wn                                                       00000830
c        in this program ws is the mass of the delta;                   00000840
c        wsn is ws divided by wn                                        00000850
      ws=wdd(inter)                                                     00000860
      dws=1.d0/ws                                                       00000870
      wsn=ws*dwn                                                        00000880
      dwsn=1.d0/wsn                                                     00000890
      dr6=1.d0/dsqrt(6.d0)                                              00000900
      dr12=1.d0/dsqrt(12.d0)                                            00000910
      r3=dsqrt(3.d0)                                                    00000920
      imad=imaa(inter)                                                  00000930
      imed=imea(inter)                                                  00000940
c                                                                       00000950
c                                                                       00000960
c        prepare constant over-all factor                               00000970
c                                                                       00000980
      fac=1.d0/(2.d0*pi)*dwn*dwn                                        00000990
c     --------------------------                                        00001000
c                                                                       00001010
c                                                                       00001020
c                                                                       00001030
c                                                                       00001040
c                                                                       00001050
c                                                                       00001060
c                                                                       00001070
c        prepare expressions depending on x and y                       00001080
c        ----------------------------------------                       00001090
c        ----------------------------------------                       00001100
c                                                                       00001110
c                                                                       00001120
c                                                                       00001130
c                                                                       00001140
   50 xa=xmev*dwn                                                       00001150
      ya=ymev*dwn                                                       00001160
      indxy=.false.                                                     00001180
      x=xa                                                              00001190
      xx=x*x                                                            00001200
      y=ya                                                              00001210
      yy=y*y                                                            00001220
      xy2=x*y*2.d0                                                      00001230
      xxpyy=xx+yy                                                       00001240
      ex=dsqrt(1.d0+xx)                                                 00001250
      ey=dsqrt(1.d0+yy)                                                 00001260
      eem12=(ex*ey-1.d0)*2.d0                                           00001270
c                                                                       00001280
c                                                                       00001290
c                                                                       00001300
c                                                                       00001310
   55 xy=xy2*0.5d0                                                      00001320
c     dxy=1.d0/xy                                                       00001330
      ee=ex*ey                                                          00001340
      eem1=ee-1.d0                                                      00001350
      eme=ex-ey                                                         00001360
      eep1=ee+1.d0                                                      00001370
       epe=ex+ey                                                        00001380
      xxyy=xx*yy                                                        00001390
c                                                                       00001400
c                                                                       00001410
      xh=x*0.5d0                                                        00001420
      xpy=x+y                                                           00001430
      xmy=x-y                                                           00001440
      xhpy=xh+y                                                         00001450
      xhmy=xh-y                                                         00001460
      y2=y*2.d0                                                         00001470
c                                                                       00001480
      wx=ex+1.d0                                                        00001490
      xw=x/wx                                                           00001500
      wy=ey+1.d0                                                        00001510
      yw=y/wy                                                           00001520
      ww=wx*wy                                                          00001530
      dww=1.d0/ww                                                       00001540
      xyww=xy*dww                                                       00001550
c                                                                       00001560
c                                                                       00001570
      ys=ymev*dws                                                       00001580
      ys2=ys*2.d0                                                       00001590
      ysys=ys*ys                                                        00001600
      eyss=dsqrt(1.d0+ysys)                                             00001610
      wyss=eyss+1.d0                                                    00001620
      eys=eyss*wsn                                                      00001630
      wys=eys+wsn                                                       00001640
      yws=y/wys                                                         00001650
      eysx=eyss*x                                                       00001660
      eysy=eyss*y                                                       00001670
      wws=wx*wys                                                        00001680
      dwws=1.d0/wws                                                     00001690
      xywws=xy*dwws                                                     00001700
      eerel=ex*dsqrt(ey*eyss)                                           00001710
c                                                                       00001720
c        the fundamental expressions e1, e2, e3, and e4                 00001730
c         and e5, e6, e7 and e8                                         00001740
c                                                                       00001750
      e1=1.d0-xywws                                                     00001760
      e2=1.d0+xywws                                                     00001770
      e3=xw-yw                                                          00001780
      e4=xw+yw                                                          00001790
      e5=1.d0-xyww                                                      00001800
      e6=1.d0+xyww                                                      00001810
      e7=xw-yws                                                         00001820
      e8=xw+yws                                                         00001830
c                                                                       00001840
c                                                                       00001850
c                                                                       00001860
c                                                                       00001870
c        prepare over-all factor                                        00001880
c                                                                       00001890
c                                                                       00001900
      go to (70,71,72,73,74),iftgo                                      00001910
c                                                                       00001920
c        no additional factor                                           00001930
c                                                                       00001940
   70 fff=fac                                                           00001950
      go to 90                                                          00001960
c                                                                       00001970
c        minimal relativity                                             00001980
c                                                                       00001990
   71 fff=fac/dsqrt(ee)                                                 00002000
      go to 90                                                          00002010
c                                                                       00002020
c        factor m/e*m/e                                                 00002030
c                                                                       00002040
   72 fff=fac/ee                                                        00002050
      go to 90                                                          00002060
c                                                                       00002070
c        minimal relativity for nd                                      00002080
c                                                                       00002090
   73 fff=fac/dsqrt(eerel)                                              00002100
      go to 90                                                          00002110
c                                                                       00002120
c        factor m/e*m/e for nd                                          00002130
c                                                                       00002140
   74 fff=fac/eerel                                                     00002150
c                                                                       00002160
c                                                                       00002170
c                                                                       00002180
c                                                                       00002190
c                                                                       00002200
c                                                                       00002210
   90 do 93 iv=1,6                                                      00002220
   93 v(iv)=0.d0                                                        00002230
      do 95 il=imad,imed                                                00002240
      do 95 iv=1,32                                                     00002250
   95 vj(iv,il)=0.d0                                                    00002260
c                                                                       00002270
c                                                                       00002280
c                                                                       00002290
c                                                                       00002300
c        contributions of mesons                                        00002310
c        -----------------------                                        00002320
c        -----------------------                                        00002330
c                                                                       00002340
c                                                                       00002350
c                                                                       00002360
c                                                                       00002370
      do 1995 img=1,mge                                                 00002380
      mg=mggo(img,inter)                                                00002390
      if (mg.eq.0) go to 2000                                           00002400
      me=mgg(mg,inter)                                                  00002410
      go to (100,9000,300,9000,9000,600,600,800,900,900,9000,9000),mg   00002420
c                                                                       00002430
c                                                                       00002440
c                                                                       00002450
c                                                                       00002460
c        0-  , pseudo-scalar mesons                                     00002470
c        --------------------------                                     00002480
c                                                                       00002490
c                                                                       00002500
c                                                                       00002510
c                                                                       00002520
  100 mc=1                                                              00002530
c                                                                       00002540
c                                                                       00002550
      indfa=.false.                                                     00002560
c                                                                       00002570
c        abbreviations                                                  00002580
c                                                                       00002590
  111 e13=e1*e3                                                         00002600
      e23=e2*e3                                                         00002610
      e13xh=e13*xh                                                      00002620
      e23xh=e23*xh                                                      00002630
      e13sy=e13*eysy                                                    00002640
      e23sy=e23*eysy                                                    00002650
      e13sx=e13*eysx                                                    00002660
      e23sx=e23*eysx                                                    00002670
c                                                                       00002680
c        groundstructure of factors                                     00002690
c                                                                       00002700
      fa(1)=-e13xh*r3                                                   00002710
      fa(2)=fa(1)                                                       00002720
      fa(3)=-e13sy-e23xh                                                00002730
      fa(4)=-e13sy+e13sx                                                00002740
      fa(5)= e23xh+e13sx                                                00002750
      fa(6)=-e13xh+e23sy                                                00002760
      fa(7)=-e23sx-e13xh                                                00002770
      fa(8)= e23xh*r3                                                   00002780
c                                                                       00002790
c        insert (-y) for y                                              00002800
      a =e1                                                             00002810
      e1=e2                                                             00002820
      e2=a                                                              00002830
      a =e3                                                             00002840
      e3=e4                                                             00002850
      e4=a                                                              00002860
      eysy=-eysy                                                        00002870
c                                                                       00002880
      if (indfa) go to 113                                              00002890
c                                                                       00002900
c        factors f from 1 to 8                                          00002910
      do 112 i=1,8                                                      00002920
  112 f(i)=fa(i)                                                        00002930
      indfa=.true.                                                      00002940
      go to 111                                                         00002950
c                                                                       00002960
c        factors f from 9 to 16                                         00002970
  113 f( 9)=fa(8)                                                       00002980
      f(10)=fa(6)                                                       00002990
      f(11)=fa(7)                                                       00003000
      f(12)=fa(3)                                                       00003010
      f(13)=fa(4)                                                       00003020
      f(14)=fa(5)                                                       00003030
      f(15)=fa(1)                                                       00003040
      f(16)=fa(2)                                                       00003050
c                                                                       00003060
c        factors f from 17 to 30                                        00003070
      f(17)=f( 8)                                                       00003080
      f(18)=f( 8)                                                       00003090
      f(19)=f( 6)                                                       00003100
      f(20)=f( 7)                                                       00003110
      f(21)=f( 3)                                                       00003120
      f(22)=f( 5)                                                       00003130
      f(23)=f( 1)                                                       00003140
      f(24)=f(15)                                                       00003150
      f(25)=f(12)                                                       00003160
      f(26)=f(14)                                                       00003170
      f(27)=f(10)                                                       00003180
      f(28)=f(11)                                                       00003190
      f(29)=f( 9)                                                       00003200
      f(30)=f( 9)                                                       00003210
c                                                                       00003220
c                                                                       00003230
      ffc=-wx*dsqrt(wy*wyss)*dr12                                       00003240
      if (mg.eq.3) ffc=-4.d0*dr12                                       00003250
c                                                                       00003260
      do 117 mx=1,me                                                    00003270
      ff=ffc/dsqrt(c(4,ima(mx,mg,inter)))                               00003280
  117 call obstr (1,mx,mx)                                              00003290
      go to 1995                                                        00003300
c                                                                       00003310
c                                                                       00003320
c                                                                       00003330
c                                                                       00003340
c        0-st, pseudo-scalar mesons in static limit                     00003350
c        ------------------------------------------                     00003360
c                                                                       00003370
c                                                                       00003380
c                                                                       00003390
c                                                                       00003400
  300 mc=1                                                              00003410
c                                                                       00003420
      indfa=.false.                                                     00003430
c                                                                       00003440
c        abbreviations                                                  00003450
c                                                                       00003460
  311 e3=(x-y)*0.5d0                                                    00003470
      e4=(x+y)*0.5d0                                                    00003480
      e3xh=e3*xh                                                        00003490
      e3x=e3*x                                                          00003500
      e3y=e3*y                                                          00003510
c                                                                       00003520
c        groundstructure of factors                                     00003530
c                                                                       00003540
      fa(1)=-e3xh*r3                                                    00003550
      fa(2)=fa(1)                                                       00003560
      fa(3)=-e3y-e3xh                                                   00003570
      fa(4)=-e3y+e3x                                                    00003580
      fa(5)= e3xh+e3x                                                   00003590
      fa(6)=-e3xh+e3y                                                   00003600
      fa(7)=-e3x-e3xh                                                   00003610
      fa(8)=-fa(1)                                                      00003620
c                                                                       00003630
c        insert (-y) for y                                              00003640
      a =e3                                                             00003650
      e3=e4                                                             00003660
      e4=a                                                              00003670
      y =-y                                                             00003680
c                                                                       00003690
      if(indfa) go to 113                                               00003700
c                                                                       00003710
c        factors f from 1 to 8                                          00003720
c                                                                       00003730
      do 312 i=1,8                                                      00003740
  312 f(i)=fa(i)                                                        00003750
      indfa=.true.                                                      00003760
      go to 311                                                         00003770
c                                                                       00003780
c                                                                       00003790
c                                                                       00003800
c                                                                       00003810
c        1-  , vector mesons                                            00003820
c        -------------------                                            00003830
c                                                                       00003840
c                                                                       00003850
c                                                                       00003860
c                                                                       00003870
  600 emea=eme                                                          00003880
      eme=0.d0                                                          00003890
c                                                                       00003900
c                                                                       00003910
c                                                                       00003920
c                                                                       00003930
c        vector coupling                                                00003940
c                                                                       00003950
c                                                                       00003960
c                                                                       00003970
c                                                                       00003980
  610 mc=1                                                              00003990
c                                                                       00004000
c                                                                       00004010
      indfa=.false.                                                     00004020
c                                                                       00004030
c        abbreviations                                                  00004040
c                                                                       00004050
  611 e14=e1*e4                                                         00004060
      e24=e2*e4                                                         00004070
      e47=e4*e7                                                         00004080
      e48=e4*e8                                                         00004090
      e67=e6*e7                                                         00004100
      e68=e6*e8                                                         00004110
c                                                                       00004120
      e14y=e14*y                                                        00004130
      e24y=e24*y                                                        00004140
      e47a=e47*x*ys2                                                    00004150
      e47b=e47*y*ys2                                                    00004160
      e48a=e48*xpy*ys2                                                  00004170
c                                                                       00004180
      e67y=e67*y                                                        00004190
      e68y=e68*y                                                        00004200
c                                                                       00004210
      e68p24=e68+e24                                                    00004220
      e68m24=e68-e24                                                    00004230
      e67p14=e67+e14                                                    00004240
      e67m14=e67-e14                                                    00004250
c                                                                       00004260
      e68py =e68p24*y                                                   00004270
      e68pxh=e68p24*xh                                                  00004280
      e68psx=e68p24*eysx                                                00004290
      e68mxh=e68m24*xh                                                  00004300
      e68msx=e68m24*eysx                                                00004310
c                                                                       00004320
      e67py =e67p14*y                                                   00004330
      e67px =e67p14*x                                                   00004340
      e67pxh=e67p14*xh                                                  00004350
      e67psx=e67p14*eysx                                                00004360
      e67my =e67m14*y                                                   00004370
      e67mxh=e67m14*xh                                                  00004380
      e67msx=e67m14*eysx                                                00004390
      if (mg.eq.7) go to 710                                            00004400
c                                                                       00004410
c        ground structure of factors                                    00004420
c                                                                       00004430
      fa( 1)=-(e68pxh+e24y)*r3                                          00004440
      fa( 2)=- e68mxh*r3                                                00004450
      fa( 3)=  e67mxh-e67py-e68msx-e48a                                 00004460
      fa( 4)=  e67px-e67my+e48a                                         00004470
      fa( 5)=  e67mxh+e68msx                                            00004480
      fa( 6)=  e68pxh+e68y-e67psx-e47a                                  00004490
      fa( 7)=- e68mxh-e67msx                                            00004500
      fa( 8)=  e67mxh*r3                                                00004510
      fa( 9)=-(e67mxh-e14y)*r3                                          00004520
      fa(10)=- e67pxh*r3                                                00004530
      fa(11)=- e68pxh-e68py+e67msx+e47b                                 00004540
      fa(12)=  e68pxh+e67psx                                            00004550
      fa(13)=- e67mxh+e67y+e68psx                                       00004560
      fa(14)=- e67pxh-e68psx                                            00004570
      fa(15)=  e68pxh*r3                                                00004580
c                                                                       00004590
c                                                                       00004600
c        insert (-y) for y                                              00004610
  613 a =e2                                                             00004620
      e2=e1                                                             00004630
      e1=a                                                              00004640
      a =e4                                                             00004650
      e4=e3                                                             00004660
      e3=a                                                              00004670
      a =e6                                                             00004680
      e6=e5                                                             00004690
      e5=a                                                              00004700
      a =e7                                                             00004710
      e7=e8                                                             00004720
      e8=a                                                              00004730
      y=-y                                                              00004740
      ys=-ys                                                            00004750
      ys2=-ys2                                                          00004760
      eysy=-eysy                                                        00004770
      a=xpy                                                             00004780
      xpy=xmy                                                           00004790
      xmy=a                                                             00004800
c                                                                       00004810
      if (indfa) go to 615                                              00004820
c                                                                       00004830
c        factors f from 1 to 8                                          00004840
c         and     from 17 to 23                                         00004850
      do 614 i=1,8                                                      00004860
      if (i.ge.8) go to 614                                             00004870
      f(16+i)=fa(8+i)                                                   00004880
  614 f(i)=fa(i)                                                        00004890
      indfa=.true.                                                      00004900
      go to 611                                                         00004910
c                                                                       00004920
c        factors f from 9 to 16                                         00004930
  615 f( 9)=fa(8)                                                       00004940
      f(10)=fa(6)                                                       00004950
      f(11)=fa(7)                                                       00004960
      f(12)=fa(3)                                                       00004970
      f(13)=fa(4)                                                       00004980
      f(14)=fa(5)                                                       00004990
      f(15)=fa(1)                                                       00005000
      f(16)=fa(2)                                                       00005010
c                                                                       00005020
c        factors f from 24 to 30                                        00005030
      f(24)=fa(15)                                                      00005040
      f(25)=fa(13)                                                      00005050
      f(26)=fa(14)                                                      00005060
      f(27)=fa(11)                                                      00005070
      f(28)=fa(12)                                                      00005080
      f(29)=fa( 9)                                                      00005090
      f(30)=fa(10)                                                      00005100
c                                                                       00005110
c                                                                       00005120
      ffc=-wx*dsqrt(wy*wyss)*dr12                                       00005130
      if (mg.eq.9.or.mg.eq.10) go to 940                                00005140
c                                                                       00005150
      do 617 mx=1,me                                                    00005160
      ff=ffc/dsqrt(c(4,ima(mx,mg,inter)))                               00005170
  617 call obstr (1,mx,mx)                                              00005180
c                                                                       00005190
c                                                                       00005200
c                                                                       00005210
c                                                                       00005220
c        tensor coupling                                                00005230
c                                                                       00005240
c                                                                       00005250
c                                                                       00005260
c                                                                       00005270
      mc=3                                                              00005280
c                                                                       00005290
c                                                                       00005300
      indfa=.false.                                                     00005310
c                                                                       00005320
c        abbreviations                                                  00005330
c                                                                       00005340
  621 e13=e1*e3*0.5d0                                                   00005350
      e15=e1*e5*0.5d0                                                   00005360
      e23=e2*e3*0.5d0                                                   00005370
      e25=e2*e5*0.5d0                                                   00005380
      e37=e3*e7*0.5d0                                                   00005390
      e38=e3*e8*0.5d0                                                   00005400
      e57=e5*e7*0.5d0                                                   00005410
      e58=e5*e8*0.5d0                                                   00005420
c                                                                       00005430
      a1=(epe*eyss+ys*y2)*x                                             00005440
c                                                                       00005450
c                                                                       00005460
      e13xh =e13*xh                                                     00005470
      e13xhp=e13*xhpy                                                   00005480
      e13sx =e13*eysx                                                   00005490
      e13eme=e13*eme                                                    00005500
      e13ea =e13eme*xh                                                  00005510
      e13eb =e13eme*xhpy                                                00005520
      e13ec =e13eme*xpy                                                 00005530
c                                                                       00005540
      e15xy =e15*x*y                                                    00005550
c                                                                       00005560
      e23xh =e23*xh                                                     00005570
      e23xhp=e23*xhpy                                                   00005580
      e23sx =e23*eysx                                                   00005590
      e23eme=e23*eme                                                    00005600
      e23ea =e23eme*xh                                                  00005610
      e23eb =e23eme*xhpy                                                00005620
c                                                                       00005630
      e25xy =e25*x*y                                                    00005640
c                                                                       00005650
      e37xys=e37*x*ys2                                                  00005660
      e37yys=e37*y*ys2                                                  00005670
c                                                                       00005680
      e38xya=e38*xpy*ys2                                                00005690
      e38xye=e38xya*eme                                                 00005700
c                                                                       00005710
      e57epe=e57*epe                                                    00005720
      e57ea =e57epe*xh                                                  00005730
      e57eb =e57epe*xhmy                                                00005740
      e57ec =e57epe*xmy                                                 00005750
      e57a1 =e57*a1                                                     00005760
c                                                                       00005770
      e58epe=e58*epe                                                    00005780
      e58ea =e58epe*xh                                                  00005790
      e58eb =e58epe*xhpy                                                00005800
      e58a1 =e58*a1                                                     00005810
c                                                                       00005820
      ea=e58ea+e25xy                                                    00005830
      eb=e15xy +e58a1                                                   00005840
      ec=e25xy -e57a1                                                   00005850
      ed=e15xy -e57ea                                                   00005860
      ef=e58eb+ec                                                       00005870
      eg=-e58ea+ec                                                      00005880
      eh=(e23xh +e13sx )*eme                                            00005890
      ei=(e13xh +e23sx )*eme                                            00005900
      ej=eb-ei                                                          00005910
      if (mg.eq.7) go to 720                                            00005920
c                                                                       00005930
c        groundstructure for factors                                    00005940
c                                                                       00005950
      fa( 1)= (ea-e23eb)*r3                                             00005960
      fa( 2)= (ea+e23ea)*r3                                             00005970
      fa( 3)= -e57eb+eb-(e13xhp+e38xya-e23sx )*eme                      00005980
      fa( 4)= -e57ec+e13ec+e38xye                                       00005990
      fa( 5)= -e57ea-eb-ei                                              00006000
      fa( 6)= -ef+(e23xh -e13sx -e37xys)*eme                            00006010
      fa( 7)= -eg+eh                                                    00006020
      fa( 8)= (ed-e13ea)*r3                                             00006030
      fa( 9)=-(ed-e13eb)*r3                                             00006040
      fa(10)=-(ed+e13ea)*r3                                             00006050
      fa(11)=  ef-(e23xhp+e13sx -e37yys)*eme                            00006060
      fa(12)=  eg+eh                                                    00006070
      fa(13)=  e57eb-ej                                                 00006080
      fa(14)=  e57ea+ej                                                 00006090
      fa(15)=-(ea-e23ea)*r3                                             00006100
c                                                                       00006110
c                                                                       00006120
c        insert (-y) for y                                              00006130
  623 a =e1                                                             00006140
      e1=e2                                                             00006150
      e2=a                                                              00006160
      a =e3                                                             00006170
      e3=e4                                                             00006180
      e4=a                                                              00006190
      a =e5                                                             00006200
      e5=e6                                                             00006210
      e6=a                                                              00006220
      a =e7                                                             00006230
      e7=e8                                                             00006240
      e8=a                                                              00006250
      y=-y                                                              00006260
      y2=-y2                                                            00006270
      ys=-ys                                                            00006280
      ys2=-ys2                                                          00006290
      eysy=-eysy                                                        00006300
      a=xpy                                                             00006310
      xpy=xmy                                                           00006320
      xmy=a                                                             00006330
      a=xhpy                                                            00006340
      xhpy=xhmy                                                         00006350
      xhmy=a                                                            00006360
c                                                                       00006370
c                                                                       00006380
      if (indfa) go to 625                                              00006390
c                                                                       00006400
c        factors f from  1 to  8                                        00006410
c        and   from     17 to 23                                        00006420
      do 624 i=1,8                                                      00006430
      if (i.ge.8) go to 624                                             00006440
      ii=16+i                                                           00006450
      f(ii)=f(ii)+fa(8+i)                                               00006460
  624 f(i)=f(i)+fa(i)                                                   00006470
      indfa=.true.                                                      00006480
      go to 621                                                         00006490
c                                                                       00006500
c        factors f from 9 to 16                                         00006510
  625 f( 9)=f( 9)+fa(8)                                                 00006520
      f(10)=f(10)+fa(6)                                                 00006530
      f(11)=f(11)+fa(7)                                                 00006540
      f(12)=f(12)+fa(3)                                                 00006550
      f(13)=f(13)+fa(4)                                                 00006560
      f(14)=f(14)+fa(5)                                                 00006570
      f(15)=f(15)+fa(1)                                                 00006580
      f(16)=f(16)+fa(2)                                                 00006590
c                                                                       00006600
c        factors f from 24 to 30                                        00006610
      f(24)=f(24)+fa(15)                                                00006620
      f(25)=f(25)+fa(13)                                                00006630
      f(26)=f(26)+fa(14)                                                00006640
      f(27)=f(27)+fa(11)                                                00006650
      f(28)=f(28)+fa(12)                                                00006660
      f(29)=f(29)+fa( 9)                                                00006670
      f(30)=f(30)+fa(10)                                                00006680
c                                                                       00006690
c                                                                       00006700
      ffc=-wx*dsqrt(wy*wyss)*dr12                                       00006710
      if (mg.eq.9.or.mg.eq.10) ffc=-4.d0*dr12                           00006720
c                                                                       00006730
      do 627 mx=1,me                                                    00006740
      ff=ffc/dsqrt(c(4,ima(mx,mg,inter)))                               00006750
  627 call obstr (1,mx,mx)                                              00006760
      eme=emea                                                          00006770
      go to 1995                                                        00006780
c                                                                       00006790
c                                                                       00006800
c                                                                       00006810
c                                                                       00006820
c        1-t ,vector mesons with lagrangian from brown et al.           00006830
c        ----------------------------------------------------           00006840
c                                                                       00006850
c                                                                       00006860
c                                                                       00006870
c        vector coupling                                                00006880
c                                                                       00006890
c                                                                       00006900
c                                                                       00006910
c                                                                       00006920
c        additional abbreviations                                       00006930
c                                                                       00006940
  710 e26=e2*e6                                                         00006950
      e16=e1*e6                                                         00006960
      e14xh=e14*xh                                                      00006970
      e24xh=e24*xh                                                      00006980
      e26a=e26*xmy*ys                                                   00006990
      e16a=e16*xpy*ys                                                   00007000
      e67xh=e67*xh                                                      00007010
      e68xh=e68*xh                                                      00007020
      e24sx=e24*eysx                                                    00007030
      e24sy=e24*eysy                                                    00007040
      e68sx=e68*eysx                                                    00007050
      e68sy=e68*eysy                                                    00007060
      e67sy=e67*eysy                                                    00007070
      e14sx=e14*eysx                                                    00007080
      e14sy=e14*eysy                                                    00007090
c                                                                       00007100
c        groundstructure of factors                                     00007110
c                                                                       00007120
      fa( 1)=-(e68pxh+e24y)*r3                                          00007130
      fa( 2)=- e68mxh*r3                                                00007140
      fa( 3)=- e26a-e24sx-2.d0*e24sy-e68sy-e67xh+3.d0*e14xh+e14y        00007150
      fa( 4)=- e26a+e68sx-e68sy+2.d0*(e24sx+e24sy-e14xh)-e14y           00007160
      fa( 5)=  e67mxh+e68msx                                            00007170
      fa( 6)=- e68xh-3.d0*e24xh+e67sy+e14sx+e16a                        00007180
      fa( 7)=- e68mxh-e67msx                                            00007190
      fa( 8)=  e67mxh*r3                                                00007200
      fa( 9)=-(e67mxh-e14y)*r3                                          00007210
      fa(10)=- e67pxh*r3                                                00007220
      fa(11)=  e68xh-e24xh+e24y-e16a-e67sy-e14sx-2.d0*e14sy             00007230
      fa(12)=  e68pxh+e67psx                                            00007240
      fa(13)=  e67xh+e14xh+e26a+e68sy+e24sx                             00007250
      fa(14)=- e67pxh-e68psx                                            00007260
      fa(15)=  e68pxh*r3                                                00007270
      go to 613                                                         00007280
c                                                                       00007290
c                                                                       00007300
c                                                                       00007310
c                                                                       00007320
c        tensor coupling                                                00007330
c                                                                       00007340
c                                                                       00007350
c                                                                       00007360
c        additional abbreviations                                       00007370
c                                                                       00007380
  720 e13ed=e13eme*y                                                    00007390
      e15ea=e15*epe*xpy*ys                                              00007400
      e15ey=e15*epe*y*ys                                                00007410
      e15sxy=2.d0*e15*eysx*y                                            00007420
      e23ec=e23eme*xhmy                                                 00007430
      e23ex=e23eme*eysx                                                 00007440
      e23ey=e23eme*eysy                                                 00007450
      e25ea=e25*epe*xmy*ys                                              00007460
      e25sxy=2.d0*e25*eysx*y                                            00007470
      e57ey=e57*epe*eysy                                                00007480
      e57ex=e57*epe*eysx                                                00007490
      e58ey=e58*epe*eysy                                                00007500
      e58ex=e58*epe*eysx                                                00007510
c                                                                       00007520
c        ground structure for factor                                    00007530
c                                                                       00007540
      fa( 1)= (ea-e23eb)*r3                                             00007550
      fa( 2)= (ea+e23ea)*r3                                             00007560
      fa( 3)=  e57ea+e58ey-e15xy+e25ea+e25sxy+3.d0*e13ea+e13ed          00007570
     1         -2.d0*e23eb*eyss                                         00007580
      fa( 4)=  e25ea-e58ex+e58ey-e13ec+2.d0*(e23ex+e23ey)               00007590
      fa( 5)=- e57ea-e58ex+e15xy-e25sxy-ei                              00007600
      fa( 6)=  e58ea-e57ey+e25xy-e15ea-e15sxy-3.d0*e23ea+e13sx*eme      00007610
      fa( 7)=  e58ea+e57ex+e25xy-e15sxy+e23ea+e13sx*eme                 00007620
      fa( 8)= (ed-e13ea)*r3                                             00007630
      fa( 9)=-(ed-e13eb)*r3                                             00007640
      fa(10)=-(ed+e13ea)*r3                                             00007650
      fa(11)=- e58ea+e57ey+e15sxy+e15ea-e25xy-e23ec-e13eb*eyss*2.d0     00007660
      fa(12)=- e58ea-e57ex-e25xy+e15sxy+e23ea+e13ea*eyss*2.d0           00007670
      fa(13)=- e57ea-e58ey+e15xy-e25sxy-e25ea+e13ea+e23ex               00007680
      fa(14)=  e57ea+e58ex-e15xy+e25sxy-e13ea-e23ex                     00007690
      fa(15)=-(ea-e23ea)*r3                                             00007700
      go to 623                                                         00007710
c                                                                       00007720
c                                                                       00007730
c                                                                       00007740
c                                                                       00007750
c        1-tt, vector mesons with eme-term                              00007760
c        ---------------------------------                              00007770
c                                                                       00007780
c                                                                       00007790
c                                                                       00007800
c                                                                       00007810
  800 emea=eme                                                          00007820
      go to 610                                                         00007830
c                                                                       00007840
c                                                                       00007850
c                                                                       00007860
c                                                                       00007870
c        1-st, vector mesons in static limit                            00007880
c        -----------------------------------                            00007890
c                                                                       00007900
c                                                                       00007910
c        vector coupling                                                00007920
c                                                                       00007930
c                                                                       00007940
  900 mc=1                                                              00007950
c                                                                       00007960
      indfa=.false.                                                     00007970
c                                                                       00007980
c        abbreviations                                                  00007990
c                                                                       00008000
  911 e3=(x-y)*0.5d0                                                    00008010
      e4=(x+y)*0.5d0                                                    00008020
      exx=xx                                                            00008030
      eyy=y*y                                                           00008040
      e3y=e3*y                                                          00008050
      e4x=e4*x                                                          00008060
      e4y=e4*y                                                          00008070
      exy=x*y                                                           00008080
      exyh=exy*0.5d0                                                    00008090
      emea=eme                                                          00008100
c                                                                       00008110
c        groundstructure of factors                                     00008120
c                                                                       00008130
      fa( 1)=-r3*(e4x+e4y)                                              00008140
      fa( 2)= 0.d0                                                      00008150
      fa( 3)=-exyh-exy                                                  00008160
      fa( 4)= exx+eyy                                                   00008170
      fa( 5)=-exyh                                                      00008180
      fa( 6)= e4x+e4y-exx                                               00008190
      fa( 7)= exy                                                       00008200
      fa( 8)=-r3*exyh                                                   00008210
      fa( 9)= r3*(exyh+e4y)                                             00008220
      fa(10)=-r3*exx*0.5d0                                              00008230
      fa(11)=-e4x-2.d0*e4y-exy                                          00008240
      fa(12)= e4x+exx                                                   00008250
      fa(13)= exyh+e3y+e4x*2.d0                                         00008260
      fa(14)=-exx*0.5d0-e4x*2.d0                                        00008270
      fa(15)= r3*e4x                                                    00008280
c                                                                       00008290
c        insert (-y) for y                                              00008300
c                                                                       00008310
      a =e3                                                             00008320
      e3=e4                                                             00008330
      e4=a                                                              00008340
      a =e4                                                             00008350
      y=-y                                                              00008360
c                                                                       00008370
      if (indfa) go to 615                                              00008380
c                                                                       00008390
c        factors f from 1 to 8                                          00008400
c         and     from 17 to 23                                         00008410
      do 914 i=1,8                                                      00008420
      if (i.ge.8) go to 914                                             00008430
      f(16+i)=fa(8+i)                                                   00008440
  914 f(i)=fa(i)                                                        00008450
      indfa=.true.                                                      00008460
      go to 911                                                         00008470
c                                                                       00008480
  940 if(mg.eq.10) go to 950                                            00008490
c                                                                       00008500
      ffc=-4.d0*dr12                                                    00008510
      do 941 mx=1,me                                                    00008520
      ff=ffc/dsqrt(c(4,ima(mx,mg,inter)))                               00008530
  941 call obstr(1,mx,mx)                                               00008540
c                                                                       00008550
c                                                                       00008560
c        tensor coupling                                                00008570
c                                                                       00008580
c        ( 1-ss ,vector mesons in static limit )                        00008590
c                                                                       00008600
c                                                                       00008610
  950 if (mg.eq.9) mc=3                                                 00008620
c                                                                       00008630
      indfa=.false.                                                     00008640
c                                                                       00008650
c        abbreviations                                                  00008660
c                                                                       00008670
  951 e3=x-y                                                            00008680
      e4=x+y                                                            00008690
      e3h=e3*0.5d0                                                      00008700
      e4h=e4*0.5d0                                                      00008710
      e3x=e3h*x                                                         00008720
      e4x=e4h*x                                                         00008730
      e3xh=e3x*0.5d0                                                    00008740
      e4xh=e4x*0.5d0                                                    00008750
      exyh=x*y*0.5d0                                                    00008760
      e3xhy=xh-y                                                        00008770
      e4xhy=xh+y                                                        00008780
      e3xe3=e3h*e3xhy                                                   00008790
      e4xe4=e4h*e4xhy                                                   00008800
c                                                                       00008810
c        groundstructure of factors                                     00008820
c                                                                       00008830
  955 fa( 1)= r3*(e4xh+exyh)                                            00008840
      fa( 2)= fa(1)                                                     00008850
      fa( 3)=-e3xe3+e4x+exyh                                            00008860
      fa( 4)=-e3*e3h                                                    00008870
      fa( 5)=-e3xh-e4x-exyh                                             00008880
      fa( 6)=-e4xe4-exyh+e3x                                            00008890
      fa( 7)= e4xh-exyh+e3x                                             00008900
      fa( 8)= r3*(exyh-e3xh)                                            00008910
      fa( 9)=-fa(8)                                                     00008920
      fa(10)= fa(9)                                                     00008930
      fa(11)=-fa(6)                                                     00008940
      fa(12)=-fa(7)                                                     00008950
      fa(13)=-fa(3)                                                     00008960
      fa(14)=-fa(5)                                                     00008970
      fa(15)=-fa(1)                                                     00008980
c                                                                       00008990
c        insert (-y) for y                                              00009000
c                                                                       00009010
      a =e3                                                             00009020
      e3=e4                                                             00009030
      e4=a                                                              00009040
      y=-y                                                              00009050
c                                                                       00009060
c                                                                       00009070
      if (indfa) go to 625                                              00009080
c                                                                       00009090
c        factors f from  1 to  8                                        00009100
c        and   from     17 to 23                                        00009110
      do 956 i=1,8                                                      00009120
      if (i.ge.8) go to 956                                             00009130
      ii=16+i                                                           00009140
      f(ii)=f(ii)+fa(8+i)                                               00009150
  956 f(i)=f(i)+fa(i)                                                   00009160
      indfa=.true.                                                      00009170
      go to 951                                                         00009180
c                                                                       00009190
c                                                                       00009200
c                                                                       00009210
c                                                                       00009220
c        this has been the end of the contributions of mesons           00009230
c                                                                       00009240
c                                                                       00009250
c                                                                       00009260
c                                                                       00009270
c        errors and warnings                                            00009280
c        -------------------                                            00009290
c                                                                       00009300
c                                                                       00009310
c                                                                       00009320
c                                                                       00009330
 9000 if (indmg(mg)) go to 1995                                         00009340
      write (kwrite,19000) mesong(mg)                                   00009350
19000 format(1h0////'0warning in obnd: meson-group  ',a4,'  does not exi00009360
     1st in this program.'/'0contribution ignored. execution continued.'00009370
     2////)                                                             00009380
      indmg(mg)=.true.                                                  00009390
c                                                                       00009400
c                                                                       00009410
c                                                                       00009420
c                                                                       00009430
 1995 continue                                                          00009440
c                                                                       00009450
c                                                                       00009460
c                                                                       00009470
c                                                                       00009480
 2000 continue                                                          00009490
c
c        antisymmetrize the nn-nd transition potential
c        it is always assumed that t=1
c
      do 2009 il=imad,imed
      if (mod(j,2).ne.0) then
c        j odd
      do 2003 iv=1,8
 2003 vj(iv,il)=0.d0
      end if
c
      do 2005 iv=9,12
      if (mod(j,2).eq.0) then
c        j even
      vj(iv,il)=(vj(iv,il)-vj(25-iv,il))*0.5d0
      vj(25-iv,il)=-vj(iv,il)
      else
c        j odd
      vj(iv,il)=(vj(iv,il)+vj(25-iv,il))*0.5d0
      vj(25-iv,il)=vj(iv,il)
      end if
 2005 continue
 2009 continue
c
c
      return                                                            00009490
      end                                                               00009500
      subroutine obdd                                                   00000010
c                                                                       00000020
c                                                                       00000030
c        one-boson-exchange nn-dd interaction;                          00000040
c        version which uses numerical integration                       00000050
c                                                                       00000060
c        author:      r. machleidt                                      00000070
c                     institut fuer theoretische kernphysik bonn        00000080
c                     nussallee 14-16                                   00000090
c                     d - 5300  bonn, w. germany                        00000100
c                                                                       00000110
c                                                                       00000120
      implicit real*8 (a-h,o-z)                                         00000130
c                                                                       00000140
c                                                                       00000150
      common /crdwrt/ kread,kwrite,kpunch,kda(9)                        00000160
c                                                                       00000170
c                                                                       00000180
c        arguments and values of this subroutine                        00000190
c                                                                       00000200
      common /cpot/   v(6),xmev,ymev                                    00000210
      common /cstate/ j,heform,sing,trip,coup,endep,label               00000220
c                                                                       00000230
c                                                                       00000240
c        this has been the end of the common-blocks containing          00000250
c        the arguments and values of this subroutine in the case of     00000260
c        no energy-dependence of the potential;                         00000270
c        in case of energy-dependence look for the common-block /cped/  00000280
c        in obai and obaa.                                              00000290
c                                                                       00000300
c        specifications for these two common blocks                     00000310
c                                                                       00000320
      logical heform,sing,trip,coup,endep                               00000330
c                                                                       00000340
c                                                                       00000350
c        common block for all ob-subroutines                            00000360
c                                                                       00000370
      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),  00000380
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,    00000390
     2                ez1,ez2,ct(96),wt(96),                            00000400
     3                ic(10,15),ift(3),mint(3),maxt(3),nt,              00000410
     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),             00000420
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,    00000430
     6                indc(2,15),indpar(3),indxy                        00000440
c                                                                       00000450
c         specifications for this common block                          00000460
c                                                                       00000470
      logical indc,indxy,indpar                                         00000480
c                                                                       00000490
c                                                                       00000500
c        further specifications                                         00000510
c                                                                       00000520
      logical indmg(12),index,indfa                                     00000530
      character*4 mesong(12)
      dimension fa(31)                                                  00000550
      data pi/3.141592653589793d0/                                      00000560
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',                   00000570
     1                   '1-  ','1-t ','1-tt','1-st','1-ss',            00000580
     2                    '1+  ','2+  '/                                00000590
      data indmg/12*.false./,index/.false./                             00000600
c
c
c                                                                       00000630
c                                                                       00000640
      inter=3                                                           00000650
c                                                                       00000660
c                                                                       00000670
c                                                                       00000680
c                                                                       00000690
c        call subroutine obpar once and only once                       00000700
c                                                                       00000710
c                                                                       00000720
      if (index) go to 50                                               00000730
      index=.true.                                                      00000740
      if (indpar(inter)) go to 40                                       00000750
c                                                                       00000760
c                                                                       00000770
      call obpar                                                        00000780
c                                                                       00000790
c                                                                       00000800
   40 iftgo=ift(inter)+1                                                00000810
      wn=wnn(inter)                                                     00000820
      dwn=1.d0/wn                                                       00000830
c        in this program ws is the mass of the delta                    00000840
c        wsn is ws devided by wn                                        00000850
      ws=wdd(inter)                                                     00000860
      dws=1.d0/ws                                                       00000870
      wsn=ws*dwn                                                        00000880
      dwsn=1.d0/wsn                                                     00000890
      r3=dsqrt(3.d0)                                                    00000900
      r33=r3/3.d0                                                       00000910
      r332=r33*2.d0                                                     00000920
      dr34=1.d0/(4.d0*r3)                                               00000930
      imadd=imaa(inter)                                                 00000940
      imedd=imea(inter)                                                 00000950
c                                                                       00000960
c                                                                       00000970
c        prepare constant over-all factor                               00000980
c                                                                       00000990
      fac=1.d0/(2.d0*pi)*dwn*dwn                                        00001000
c     --------------------------                                        00001010
c                                                                       00001020
c                                                                       00001030
c                                                                       00001040
c                                                                       00001050
c                                                                       00001060
c                                                                       00001070
c                                                                       00001080
c        prepare expressions depending on x and y                       00001090
c        ----------------------------------------                       00001100
c        ----------------------------------------                       00001110
c                                                                       00001120
c                                                                       00001130
c                                                                       00001140
c                                                                       00001150
   50 xa=xmev*dwn                                                       00001160
      ya=ymev*dwn                                                       00001170
      indxy=.false.                                                     00001190
      x=xa                                                              00001200
      xx=x*x                                                            00001210
      y=ya                                                              00001220
      yy=y*y                                                            00001230
      xy2=x*y*2.d0                                                      00001240
      xxpyy=xx+yy                                                       00001250
      ex=dsqrt(1.d0+xx)                                                 00001260
      ey=dsqrt(1.d0+yy)                                                 00001270
      eem12=(ex*ey-1.d0)*2.d0                                           00001280
c                                                                       00001290
c                                                                       00001300
c                                                                       00001310
c                                                                       00001320
   55 xy=xy2*0.5d0                                                      00001330
c     dxy=1.d0/xy                                                       00001340
      ee=ex*ey                                                          00001350
      eem1=ee-1.d0                                                      00001360
      eme=ex-ey                                                         00001370
      eep1=ee+1.d0                                                      00001380
       epe=ex+ey                                                        00001390
      xxyy=xx*yy                                                        00001400
c                                                                       00001410
c                                                                       00001420
      wx=ex+1.d0                                                        00001430
      xw=x/wx                                                           00001440
c                                                                       00001450
      ys=ymev*dws                                                       00001460
      ys2=ys*2.d0                                                       00001470
      ysys=ys*ys                                                        00001480
      ysys8=ysys*8.d0                                                   00001490
      eyss=dsqrt(1.d0+ysys)                                             00001500
      wyss=eyss+1.d0                                                    00001510
      eyss2=eyss*2.d0                                                   00001520
      eys=eyss*wsn                                                      00001530
      wys=eys+wsn                                                       00001540
      yws=y/wys                                                         00001550
      wws=wx*wys                                                        00001560
      dwws=1.d0/wws                                                     00001570
      xywws=xy*dwws                                                     00001580
      eerel=ex*eyss                                                     00001590
c                                                                       00001600
c                                                                       00001610
c        the fundamental expressions e1, e2, e7 and e8                  00001620
c                                                                       00001630
      e1=1.d0-xywws                                                     00001640
      e2=1.d0+xywws                                                     00001650
      e7=xw-yws                                                         00001660
      e8=xw+yws                                                         00001670
c                                                                       00001680
c                                                                       00001690
c                                                                       00001700
c                                                                       00001710
c        prepare over-all factor                                        00001720
c                                                                       00001730
c                                                                       00001740
      go to (70,71,72,73,74),iftgo                                      00001750
c                                                                       00001760
c        no additional factor                                           00001770
c                                                                       00001780
   70 fff=fac                                                           00001790
      go to 90                                                          00001800
c                                                                       00001810
c        minimal relativity                                             00001820
c                                                                       00001830
   71 fff=fac/dsqrt(ee)                                                 00001840
      go to 90                                                          00001850
c                                                                       00001860
c        factor m/e*m/e                                                 00001870
c                                                                       00001880
   72 fff=fac/ee                                                        00001890
      go to 90                                                          00001900
c                                                                       00001910
c        minimal relativity for dd                                      00001920
c                                                                       00001930
   73 fff=fac/dsqrt(eerel)                                              00001940
      go to 90                                                          00001950
c                                                                       00001960
c        factor m/e*m/e for dd                                          00001970
c                                                                       00001980
   74 fff=fac/eerel                                                     00001990
c                                                                       00002000
c                                                                       00002010
c                                                                       00002020
c                                                                       00002030
c                                                                       00002040
c                                                                       00002050
   90 do 93 iv=1,6                                                      00002060
   93 v(iv)=0.d0                                                        00002070
      do 95 il=imadd,imedd                                              00002080
      do 95 iv=1,32                                                     00002090
   95 vj(iv,il)=0.d0                                                    00002100
c                                                                       00002110
c                                                                       00002120
c                                                                       00002130
c                                                                       00002140
c        contributions of mesons                                        00002150
c        -----------------------                                        00002160
c        -----------------------                                        00002170
c                                                                       00002180
c                                                                       00002190
c                                                                       00002200
c                                                                       00002210
      do 1995 img=1,mge                                                 00002220
      mg=mggo(img,inter)                                                00002230
      if (mg.eq.0) go to 2000                                           00002240
      me=mgg(mg,inter)                                                  00002250
      go to (100,9000,300,9000,9000,600,600,9000,900,9000,9000,9000),mg 00002260
c                                                                       00002270
c                                                                       00002280
c                                                                       00002290
c                                                                       00002300
c        0-  , pseudo-scalar mesons                                     00002310
c        --------------------------                                     00002320
c                                                                       00002330
c                                                                       00002340
c                                                                       00002350
c                                                                       00002360
c                                                                       00002370
  100 mc=1                                                              00002380
c                                                                       00002390
c                                                                       00002400
      indfa=.false.                                                     00002410
c                                                                       00002420
c        abbreviations                                                  00002430
c                                                                       00002440
  111 e11=e1*e1                                                         00002450
      e12=e1*e2                                                         00002460
      e22=e2*e2                                                         00002470
c                                                                       00002480
c                                                                       00002490
      e11xx=e11*xx                                                      00002500
      e12xx=e12*xx                                                      00002510
      e22xx=e22*xx                                                      00002520
c                                                                       00002530
c                                                                       00002540
      e11s=e11*eyss2                                                    00002550
      e12s=e12*eyss2                                                    00002560
      e22s=e22*eyss2                                                    00002570
c                                                                       00002580
c                                                                       00002590
      e11z=e11s*eyss2                                                   00002600
      e12z=e12s*eyss2                                                   00002610
      e22z=e22s*eyss2                                                   00002620
c                                                                       00002630
c                                                                       00002640
      e11sxx=e11s*xx                                                    00002650
      e12sxx=e12s*xx                                                    00002660
      e22sxx=e22s*xx                                                    00002670
c                                                                       00002680
      e11sxy=e11s*xy                                                    00002690
      e12sxy=e12s*xy                                                    00002700
      e22sxy=e22s*xy                                                    00002710
c                                                                       00002720
      e11zxx=e11z*xx                                                    00002730
      e12zxx=e12z*xx                                                    00002740
      e22zxx=e22z*xx                                                    00002750
c                                                                       00002760
      e11zyy=e11z*yy                                                    00002770
      e12zyy=e12z*yy                                                    00002780
      e22zyy=e22z*yy                                                    00002790
c                                                                       00002800
      e11zxy=e11z*xy                                                    00002810
      e12zxy=e12z*xy                                                    00002820
      e22zxy=e22z*xy                                                    00002830
c                                                                       00002840
c                                                                       00002850
c        groundstructure of factors                                     00002860
c                                                                       00002870
      fa(1)= e11xx*r3                                                   00002880
      fa(2)=0.d0                                                        00002890
      fa(3)=-fa(1)                                                      00002900
      fa(4)=-e11sxy-e12xx                                               00002910
      fa(5)=-e11sxy+e11sxx                                              00002920
      fa(6)= e12xx+e11sxx                                               00002930
      fa(7)=-e11xx+e12sxy                                               00002940
      fa(8)=-e11xx-e12sxx                                               00002950
      if (.not.indfa)                                                   00002960
     1fa(9)= e12xx*r3                                                   00002970
      fa(10)=-e22xx-e11zyy-e12sxy*2.d0                                  00002980
      fa(11) =-e11zyy+e11zxy*2.d0+e22xx+e12sxx*2.d0                     00002990
      fa(12)= e22xx+e11zxy*2.d0-e11zxx+e12sxy*2.d0                      00003000
      fa(13)=-e11zxx-e22xx-e12sxx*2.d0                                  00003010
      if (indfa) go to 112                                              00003020
      fa(14)=-e12xx+e12zyy-e11sxy+e22sxy                                00003030
      fa(15)=-e12zxy*2.d0-e22sxx-e22sxy+e11sxx-e11sxy                   00003040
      fa(16)= e12xx+e12zxx+e11sxx+e22sxx                                00003050
  112 continue                                                          00003060
c                                                                       00003070
c                                                                       00003080
c                                                                       00003090
      do 113 i=10,16                                                    00003100
  113 fa(i)=fa(i)*r33                                                   00003110
c                                                                       00003120
c        insert (-y) for y                                              00003130
      a=e1                                                              00003140
      e1=e2                                                             00003150
      e2=a                                                              00003160
      xy=-xy                                                            00003170
c                                                                       00003180
      if (indfa) go to 115                                              00003190
c                                                                       00003200
      do 114 i=1,16                                                     00003210
  114 f(i)=fa(i)                                                        00003220
      indfa=.true.                                                      00003230
      go to 111                                                         00003240
c                                                                       00003250
c        factors f                                                      00003260
c                                                                       00003270
  115 f(17)=fa( 7)                                                      00003280
      f(18)=fa( 8)                                                      00003290
      f(19)=fa(10)                                                      00003300
      f(20)=fa(11)                                                      00003310
      f(21)=fa(12)                                                      00003320
      f(22)=fa(13)                                                      00003330
      f(23)=fa( 4)                                                      00003340
      f(24)=fa( 5)                                                      00003350
      f(25)=fa( 6)                                                      00003360
      f(26)=fa( 1)                                                      00003370
      f(27)=fa( 2)                                                      00003380
      f(28)=fa( 3)                                                      00003390
c                                                                       00003400
c                                                                       00003410
      f(29)=-f(9)                                                       00003420
      f(30)=0.d0                                                        00003430
      f(31)=f( 9)                                                       00003440
      f(32)=-f(7)                                                       00003450
      f(33)=f( 7)-f( 8)                                                 00003460
      f(34)=f( 8)                                                       00003470
      f(35)=f( 4)                                                       00003480
      f(36)=f(6)                                                        00003490
      f(37)=-f(1)                                                       00003500
      f(38)=f(14)                                                       00003510
      f(39)=f(15)                                                       00003520
      f(40)=f(16)                                                       00003530
      f(41)=f(10)                                                       00003540
      f(42)=f(11)-f(10)                                                 00003550
      f(43)=f(13)                                                       00003560
c                                                                       00003570
      f(44)=-f(17)                                                      00003580
      f(45)=f(17)-f(18)                                                 00003590
      f(46)=f(18)                                                       00003600
      f(47)=f(23)                                                       00003610
      f(48)=f(25)                                                       00003620
      f(49)=f(19)                                                       00003630
      f(50)=f(20)-f(19)                                                 00003640
      f(51)=f(22)                                                       00003650
      f(52)=-f(26)                                                      00003660
c                                                                       00003670
c                                                                       00003680
      ffc=-dr34*wx*wyss*0.5d0                                           00003690
      if (mg.eq.3) ffc=-dr34*2.d0                                       00003700
c                                                                       00003710
      do 117 mx=1,me                                                    00003720
      ff=ffc/c(4,ima(mx,mg,inter))                                      00003730
  117 call obstr (1,mx,mx)                                              00003740
      go to 1995                                                        00003750
c                                                                       00003760
c                                                                       00003770
c                                                                       00003780
c                                                                       00003790
c        0-st, pseudo-scalar mesons in static limit                     00003800
c        ------------------------------------------                     00003810
c                                                                       00003820
c                                                                       00003830
c                                                                       00003840
  300 mc=1                                                              00003850
c                                                                       00003860
      indfa=.false.                                                     00003870
c                                                                       00003880
c        abbreviations                                                  00003890
c                                                                       00003900
  311 exx=xx                                                            00003910
      e2xx=2.d0*xx                                                      00003920
      e2xy=2.d0*xy                                                      00003930
      e3xx=3.d0*xx                                                      00003940
      e4xx=2.d0*e2xx                                                    00003950
      e4xy=2.d0*e2xy                                                    00003960
      e4yy=4.d0*yy                                                      00003970
c                                                                       00003980
c                                                                       00003990
c        groundstructure of factors                                     00004000
c                                                                       00004010
      fa(1) =exx*r3                                                     00004020
      fa(2) =0.d0                                                       00004030
      fa(3) =-fa(1)                                                     00004040
      fa(4) =-e2xy-exx                                                  00004050
      fa(5) =e2xx-e2xy                                                  00004060
      fa(6) =e3xx                                                       00004070
      fa(7) =-exx+e2xy                                                  00004080
      fa(8) =-e3xx                                                      00004090
      if (.not.indfa)                                                   00004100
     1fa(9) = exx*r3                                                    00004110
      fa(10)=-exx-e4yy-e4xy                                             00004120
      fa(11)=e3xx+e2xx+2.d0*e4xy-e4yy                                   00004130
      fa(12)=3.d0*e4xy-e3xx                                             00004140
      fa(13)=-3.d0*e3xx                                                 00004150
c                                                                       00004160
      if(indfa) go to 312                                               00004170
      fa(14)=-exx+e4yy                                                  00004180
      fa(15)=-3.d0*e4xy                                                 00004190
      fa(16)= 3.d0*e3xx                                                 00004200
  312 continue                                                          00004210
c                                                                       00004220
c                                                                       00004230
      do 313 i=10,16                                                    00004240
  313 fa(i)=fa(i)*r33                                                   00004250
c                                                                       00004260
c        insert (-y) for y                                              00004270
c                                                                       00004280
      xy=-xy                                                            00004290
c                                                                       00004300
c        factors f                                                      00004310
c                                                                       00004320
      if (indfa) go to 115                                              00004330
c                                                                       00004340
      do 314 i=1,16                                                     00004350
  314 f(i)=fa(i)                                                        00004360
      indfa=.true.                                                      00004370
      go to 311                                                         00004380
c                                                                       00004390
c                                                                       00004400
c                                                                       00004410
c                                                                       00004420
c        1-  , vector mesons                                            00004430
c        -------------------                                            00004440
c                                                                       00004450
c                                                                       00004460
c                                                                       00004470
c                                                                       00004480
  600 mc=1                                                              00004490
c                                                                       00004500
c                                                                       00004510
      indfa=.false.                                                     00004520
c                                                                       00004530
c        abbreviations                                                  00004540
c                                                                       00004550
  611 e11=e1*e1                                                         00004560
      e12=e1*e2                                                         00004570
      e17=e1*e7                                                         00004580
      e18=e1*e8                                                         00004590
      e22=e2*e2                                                         00004600
      e27=e2*e7                                                         00004610
      e28=e2*e8                                                         00004620
      e77=e7*e7                                                         00004630
      e78=e7*e8                                                         00004640
      e88=e8*e8                                                         00004650
c                                                                       00004660
      e11s=e11*eyss2                                                    00004670
      e12s=e12*eyss2                                                    00004680
      e17s=e17*eyss2                                                    00004690
      e18s=e18*eyss2                                                    00004700
      e22s=e22*eyss2                                                    00004710
      e27s=e27*eyss2                                                    00004720
      e28s=e28*eyss2                                                    00004730
      e77s=e77*eyss2                                                    00004740
      e78s=e78*eyss2                                                    00004750
      e88s=e88*eyss2                                                    00004760
c                                                                       00004770
      e11z=e11s*eyss2                                                   00004780
      e12z=e12s*eyss2                                                   00004790
      e17z=e17s*eyss2                                                   00004800
      e18z=e18s*eyss2                                                   00004810
      e22z=e22s*eyss2                                                   00004820
      e27z=e27s*eyss2                                                   00004830
      e28z=e28s*eyss2                                                   00004840
      e77z=e77s*eyss2                                                   00004850
      e78z=e78s*eyss2                                                   00004860
      e88z=e88s*eyss2                                                   00004870
c                                                                       00004880
      e11t=e11*ys2                                                      00004890
      e12t=e12*ys2                                                      00004900
      e17t=e17*ys2                                                      00004910
      e18t=e18*ys2                                                      00004920
      e22t=e22*ys2                                                      00004930
      e27t=e27*ys2                                                      00004940
      e28t=e28*ys2                                                      00004950
      e77t=e77*ys2                                                      00004960
      e78t=e78*ys2                                                      00004970
      e88t=e88*ys2                                                      00004980
c                                                                       00004990
      e11u=e11t*eyss2                                                   00005000
      e12u=e12t*eyss2                                                   00005010
      e17u=e17t*eyss2                                                   00005020
      e18u=e18t*eyss2                                                   00005030
      e22u=e22t*eyss2                                                   00005040
      e27u=e27t*eyss2                                                   00005050
      e28u=e28t*eyss2                                                   00005060
      e77u=e77t*eyss2                                                   00005070
      e78u=e78t*eyss2                                                   00005080
      e88u=e88t*eyss2                                                   00005090
c                                                                       00005100
      e11xx=e11*xx                                                      00005110
      e12xx=e12*xx                                                      00005120
      e22xx=e22*xx                                                      00005130
      e77xx=e77*xx                                                      00005140
      e78xx=e78*xx                                                      00005150
      e88xx=e88*xx                                                      00005160
c                                                                       00005170
      e11yy=e11*yy                                                      00005180
      e12yy=e12*yy                                                      00005190
      e22yy=e22*yy                                                      00005200
      e77yy=e77*yy                                                      00005210
      e78yy=e78*yy                                                      00005220
      e88yy=e88*yy                                                      00005230
c                                                                       00005240
      e11xy=e11*xy                                                      00005250
      e12xy=e12*xy                                                      00005260
      e22xy=e22*xy                                                      00005270
      e77xy=e77*xy                                                      00005280
      e78xy=e78*xy                                                      00005290
      e88xy=e88*xy                                                      00005300
c                                                                       00005310
      e11sxx=e11s*xx                                                    00005320
      e12sxx=e12s*xx                                                    00005330
      e22sxx=e22s*xx                                                    00005340
      e77sxx=e77s*xx                                                    00005350
      e78sxx=e78s*xx                                                    00005360
      e88sxx=e88s*xx                                                    00005370
c                                                                       00005380
      e12zxx=e12z*xx                                                    00005390
      e22zxx=e22z*xx                                                    00005400
      e78zxx=e78z*xx                                                    00005410
      e88zxx=e88z*xx                                                    00005420
c                                                                       00005430
      e17txx=e17t*xx                                                    00005440
      e18txx=e18t*xx                                                    00005450
      e27txx=e27t*xx                                                    00005460
      e28txx=e28t*xx                                                    00005470
c                                                                       00005480
      e18uxx=e18u*xx                                                    00005490
      e27uxx=e27u*xx                                                    00005500
      e28uxx=e28u*xx                                                    00005510
c                                                                       00005520
      e78zyy=e78z*yy                                                    00005530
c                                                                       00005540
      e17tyy=e17t*yy                                                    00005550
      e18tyy=e18t*yy                                                    00005560
      e27tyy=e27t*yy                                                    00005570
      e28tyy=e28t*yy                                                    00005580
c                                                                       00005590
      e18uyy=e18u*yy                                                    00005600
      e27uyy=e27u*yy                                                    00005610
c                                                                       00005620
      e11sxy=e11s*xy                                                    00005630
      e12sxy=e12s*xy                                                    00005640
      e22sxy=e22s*xy                                                    00005650
      e77sxy=e77s*xy                                                    00005660
      e78sxy=e78s*xy                                                    00005670
      e88sxy=e88s*xy                                                    00005680
c                                                                       00005690
      e17txy=e17t*xy                                                    00005700
      e18txy=e18t*xy                                                    00005710
      e27txy=e27t*xy                                                    00005720
      e28txy=e28t*xy                                                    00005730
c                                                                       00005740
      e18uxy=e18u*xy                                                    00005750
      e27uxy=e27u*xy                                                    00005760
      e28uxy=e28u*xy                                                    00005770
      if (mg.eq.7) go to 710                                            00005780
c                                                                       00005790
c                                                                       00005800
c        groundstructure of factors                                     00005810
c                                                                       00005820
      fa( 1)= (e22xx+e22yy*2.d0+e88xx)*r3                               00005830
      fa( 2)=- e22xy*4.d0*r3                                            00005840
      fa( 3)= (e22xx-e88xx)*r3                                          00005850
      fa( 4)=- e12xx-e12yy*2.d0+e22sxy-e28txx-e28tyy*2.d0-e28txy        00005860
     1        +e78xx-e78xy*2.d0-e88sxx                                  00005870
      fa( 5)=  e12xy*4.d0-e22sxx+e22sxy+e28txx+e28txy*3.d0              00005880
     1        +(e78xx-e78xy)*2.d0                                       00005890
      fa( 6)=- e12xx+e78xx-e22sxx+e88sxx                                00005900
      fa( 7)=  e22xx+e88xx+e88xy*2.d0-e12sxy-e78sxx-e27txx-e27txy       00005910
      fa( 8)=  e22xx-e88xx+e12sxx-e78sxx                                00005920
      if (.not.indfa)                                                   00005930
     1fa( 9)=-(e12xx-e78xx)*r3                                          00005940
      fa(10)= 2.d0*e12sxy-e11xx-e11yy*2.d0-e22zxx+(e28uxx+e28uxy        00005950
     1        - e18txx-e18txy-e18tyy*2.d0+e78sxx)*2.d0-e78sxy*4.d0      00005960
     2        -e77xx+(e77xy-e77yy)*4.d0   -e88zxx-(e88xx+e88xy+e88yy)   00005970
     3        *ysys8                                                    00005980
      fa(11)=-2.d0*e12sxx   +e11xx+e11xy*4.d0+e11yy*2.d0-e22zxx+(e18txx 00005990
     1        +e18txy*2.d0+e18tyy)*4.d0+2.d0*e78sxx-e77xx*3.d0+(e77xy   00006000
     2        *2.d0-e77yy)*4.d0+e88zxx+(e88xx+e88xy*2.d0+e88yy)*ysys8   00006010
      fa(12)=-2.d0*e12sxy-e11xx-e11xy*4.d0+e22zxx-(e28uxx+e28uxy+e18txx 00006020
     1        +e18txy*3.d0+e78sxx)*2.d0+4.d0*e78sxy-e77xx*3.d0          00006030
     2        +e77xy*4.d0+e88zxx-e88xy*ysys8                            00006040
      fa(13)=  e11xx-e77xx+e22zxx-e88zxx+(e12sxx-e78sxx)*2.d0           00006050
      if (indfa) go to 612                                              00006060
      fa(14)=  e12xx+e12zxx+e22sxy-e11sxy+e28txx-e28txy+e27uxx-e27uxy   00006070
     1        +e77sxx-e77sxy*2.d0-e18uxx-e18uxy-e17txx-e17txy           00006080
     2        -e78zxx-e78xx*ysys8-e78xx+e78yy*4.d0                      00006090
     3        +e88sxx+e88sxy*2.d0                                       00006100
      fa(15)=- e22sxx+e22sxy+e11sxx+e11sxy+e28txx-e28txy+e27uxx-e27uxy  00006110
     1        +e18uxx+e18uxy+e17txx+e17txy-e78xy*(4.d0-ysys8)           00006120
     2        -(e88sxx+e88sxy)*2.d0+(e77sxx-e77sxy)*2.d0                00006130
      fa(16)=- e12xx-e12zxx+e78xx+e78zxx-e11sxx+e77sxx-e22sxx+e88sxx    00006140
c                                                                       00006150
      fa(17)= (e78xx-e12xx-e12yy*2.d0)*r3                               00006160
      fa(18)=  e12xy*4.d0*r3                                            00006170
      fa(19)=-(e12xx+e78xx)*r3                                          00006180
  612 continue                                                          00006190
      fa(20)=  e12sxy+e22xx+e22yy*2.d0-e27txx+e27txy-e27tyy*2.d0-e78sxx 00006200
     1        +e88xx+e88xy*2.d0                                         00006210
      fa(21)= -e12sxx-e12sxy-e22xy*4.d0-e27txx+e27txy*3.d0-(e88xx+e88xy)00006220
     1        *2.d0                                                     00006230
      fa(22)=  e22xx+e88xx+e12sxx+e78sxx                                00006240
      fa(23)=  e12xx+e22sxy+e28txx-e28txy-e78xx+e78xy*2.d0+e88sxx       00006250
      fa(24)= -e12xx-e78xx-e22sxx-e88sxx                                00006260
      fa(25)= (e22xx+e88xx)*r3                                          00006270
      if (indfa) go to 613                                              00006280
      fa(26)= -e12xx-e12yy*2.d0+e12zxx+e22sxy-e11sxy-e28txx-e28txy      00006290
     1        -e28tyy*2.d0+e27uxx-e27uxy-e18uxx-e18uxy+e17txx-e17txy    00006300
     2        +e17tyy*2.d0+e78zxx*1.5d0-e78zyy*0.5d0-(e78xx-e78yy*5.d0) 00006310
     3        *ysys*2.d0-e78xx-e78yy*2.d0-e88sxx-e88sxy*2.d0-e77sxx     00006320
     4        +e77sxy*2.d0                                              00006330
      fa(27)=  e12xy*4.d0-e22sxx+e22sxy+e11sxx+e11sxy+e28txx+e28txy*3.d000006340
     1        +e27uxx-e27uxy+e18uxx+e18uxy+e17txx-e17txy*3.d0+e78xy*    00006350
     2        (4.d0-ysys8)+(e88sxx+e88sxy-e77sxx+e77sxy)*2.d0           00006360
      fa(28)= -e12xx-e12zxx-e78xx-e78zxx-e11sxx-e77sxx-e22sxx-e88sxx    00006370
  613 continue                                                          00006380
      fa(29)=  e12sxy*2.d0-e22zxx+e11xx+(e28uxx+e28uxy+e18txx-e18txy    00006390
     1        -e78sxx+e78sxy*2.d0)*2.d0+e88zxx-e88xy*ysys8+e77xx        00006400
     2        -(e77xy-e77yy)*4.d0                                       00006410
      fa(30)=-(e12sxx+e12sxy+e11xx+e28uxx+e28uxy+e18txx-e18txy+e78sxy   00006420
     1        *2.d0+e88zxx-e88xy*ysys*4.d0-e77xx+e77xy*2.d0)*2.d0       00006430
      fa(31)=  e11xx+e77xx+e22zxx+e88zxx+(e12sxx+e78sxx)*2.d0           00006440
c                                                                       00006450
c                                                                       00006460
c                                                                       00006470
c                                                                       00006480
c                                                                       00006490
c                                                                       00006500
  614 do 615 i=10,16                                                    00006510
  615 fa(i)=fa(i)*r33                                                   00006520
      do 616 i=26,31                                                    00006530
  616 fa(i)=fa(i)*r33                                                   00006540
c                                                                       00006550
c        insert (-y) for y                                              00006560
      a=e1                                                              00006570
      e1=e2                                                             00006580
      e2=a                                                              00006590
      a=e7                                                              00006600
      e7=e8                                                             00006610
      e8=a                                                              00006620
      ys2=-ys2                                                          00006630
      xy=-xy                                                            00006640
c                                                                       00006650
      if (indfa) go to 618                                              00006660
c                                                                       00006670
c        factors f                                                      00006680
c                                                                       00006690
      do 617 i=1,16                                                     00006700
      if (i.ge.16) go to 617                                            00006710
      f(i+28)=fa(i+16)                                                  00006720
  617 f(i)=fa(i)                                                        00006730
      indfa=.true.                                                      00006740
      go to 611                                                         00006750
c                                                                       00006760
c                                                                       00006770
  618 f(17)=fa( 7)                                                      00006780
      f(18)=fa( 8)                                                      00006790
      f(19)=fa(10)                                                      00006800
      f(20)=fa(11)                                                      00006810
      f(21)=fa(12)                                                      00006820
      f(22)=fa(13)                                                      00006830
      f(23)=fa( 4)                                                      00006840
      f(24)=fa( 5)                                                      00006850
      f(25)=fa( 6)                                                      00006860
      f(26)=fa( 1)                                                      00006870
      f(27)=fa( 2)                                                      00006880
      f(28)=fa( 3)                                                      00006890
c                                                                       00006900
      f(44)=fa(20)                                                      00006910
      f(45)=fa(21)                                                      00006920
      f(46)=fa(22)                                                      00006930
      f(47)=fa(23)                                                      00006940
      f(48)=fa(24)                                                      00006950
      f(49)=fa(29)                                                      00006960
      f(50)=fa(30)                                                      00006970
      f(51)=fa(31)                                                      00006980
      f(52)=fa(25)                                                      00006990
c                                                                       00007000
c                                                                       00007010
      ffc=-dr34*wx*wyss*0.5d0                                           00007020
      if (mg.eq.9) ffc=-dr34*2.d0                                       00007030
c                                                                       00007040
      do 619 mx=1,me                                                    00007050
      ff=ffc/c(4,ima(mx,mg,inter))                                      00007060
  619 call obstr (1,mx,mx)                                              00007070
      go to 1995                                                        00007080
c                                                                       00007090
c                                                                       00007100
c                                                                       00007110
c                                                                       00007120
c        1-t ,vector mesons with lagrangian from brown et al.           00007130
c        ----------------------------------------------------           00007140
c                                                                       00007150
c                                                                       00007160
c                                                                       00007170
c                                                                       00007180
c        additional abbreviations                                       00007190
c                                                                       00007200
  710 e11syy=e11s*yy                                                    00007210
      e12syy=e12s*yy                                                    00007220
      e12txx=e12t*xx                                                    00007230
      e12tyy=e12t*yy                                                    00007240
      e12zxy=e12z*xy                                                    00007250
      e12zyy=e12z*yy                                                    00007260
      e22syy=e22s*yy                                                    00007270
      e22txx=e22t*xx                                                    00007280
      e22txy=e22t*xy                                                    00007290
      e22tyy=e22t*yy                                                    00007300
      e22zxy=e22z*xy                                                    00007310
      e22zyy=e22z*yy                                                    00007320
      e28uyy=e28u*yy                                                    00007330
      e78zxy=e78z*xy                                                    00007340
      e88zxy=e88z*xy                                                    00007350
      e88zyy=e88z*yy                                                    00007360
c                                                                       00007370
c        ground structure for factors                                   00007380
c                                                                       00007390
      fa( 1)= (e22xx+e22yy*2.d0+e88xx)*r3                               00007400
      fa( 2)=- e22xy*4.d0*r3                                            00007410
      fa( 3)= (e22xx-e88xx)*r3                                          00007420
      fa( 4)=- e22sxx-2.d0*e22syy+e12xx+2.d0*(e12xy+e12yy)-e28txx+e28txy00007430
     1         -e78xx-e88sxy                                            00007440
      fa( 5)=  4.d0*e22sxy-2.d0*(e12xx+e12xy)-e28txx+e28txy+e88sxx      00007450
     1         -e88sxy                                                  00007460
      fa( 6)=- e12xx+e78xx-e22sxx+e88sxx                                00007470
      fa( 7)=- e22xx-2.d0*e22xy-e88xx+e12sxx+e78sxy+e18txx+e18txy       00007480
      fa( 8)=  e22xx-e88xx+e12sxx-e78sxx                                00007490
      if (.not.indfa)                                                   00007500
     1fa( 9)=-(e12xx-e78xx)*r3                                          00007510
      fa(10)=- e22zxx-2.d0*e22zyy-(e22txx-2.d0*e22txy+e22tyy)*ys2       00007520
     1         +2.d0*(e12sxx+2.d0*(e12sxy+e12syy))-5.d0*e11xx-4.d0*e11xy00007530
     2         -2.d0*(e11yy+e28uxy-e28uyy+e27txx-e27txy)-e88zyy         00007540
     3         -2.d0*e78sxy-e77xx                                       00007550
      fa(11)=  e22zxx+4.d0*e22zxy+2.d0*e22zyy-(e22txx-2.d0*e22txy       00007560
     1         +e22tyy)*ys2-(6.d0*e12sxx+8.d0*e12sxy+4.d0*e12syy)       00007570
     2         +e11xx+4.d0*e11xy+2.d0*e11yy+2.d0*(e28uxx-2.d0*e28uxy    00007580
     3         +e28uyy)+2.d0*e88zxy-e88zyy+2.d0*e78sxx+e77xx            00007590
      fa(12)=- e22zxx-4.d0*e22zxy+2.d0*e12sxx+4.d0*e12sxy+3.d0*e11xx    00007600
     1         +2.d0*(e28uxx-e28uxy+e27txx-e27txy)-e88zxx+2.d0*(e88zxy  00007610
     2         +e78sxy)+e77xx                                           00007620
      fa(13)=  e11xx-e77xx+e22zxx-e88zxx+(e12sxx-e78sxx)*2.d0           00007630
      if (indfa) go to 712                                              00007640
      fa(14)=- e22sxx-2.d0*e22sxy+(5.d0*e12xx+e12zxx+(e12txx-e12tyy)    00007650
     1         *ys2)-e11sxx+2.d0*e11sxy-e28txx+e28txy+e27uxy-e27uyy     00007660
     2         +e18uxy+e18uyy+e17txx+e17txy-e88sxy-e78xx+e78zyy+e77sxy  00007670
      fa(15)=  2.d0*(e22sxx+e22sxy)-4.d0*e12xy-2.d0*(e11sxx-e11sxy)     00007680
     1         -e28txx+e28txy-e27uxx+e27uxy-e18uxx-e18uxy-e17txx-e17txy 00007690
     2         +e88sxx-e88sxy-2.d0*e78zxy-e77sxx-e77sxy                 00007700
      fa(16)=- e12xx-e12zxx+e78xx+e78zxx-e11sxx+e77sxx-e22sxx+e88sxx    00007710
c                                                                       00007720
      fa(17)= (e78xx-e12xx-e12yy*2.d0)*r3                               00007730
      fa(18)=  e12xy*4.d0*r3                                            00007740
      fa(19)=-(e12xx+e78xx)*r3                                          00007750
  712 continue                                                          00007760
      fa(20)=- e22xx+2.d0*(e22xy-e22yy)+e12sxx+2.d0*e12syy+e18txx       00007770
     1         +e18txy-e88xx+e78sxy                                     00007780
      fa(21)=- 2.d0*(e22xx-e22xy)-4.d0*e12sxy-e18txx-e18txy-e78sxx      00007790
     1         -e78sxy                                                  00007800
      fa(22)=  e22xx+e88xx+e12sxx+e78sxx                                00007810
      fa(23)=  e22sxx-e12xx+2.d0*e12xy+e28txx-e28txy+e88sxy+e78xx       00007820
      fa(24)= -e12xx-e78xx-e22sxx-e88sxx                                00007830
      fa(25)= (e22xx+e88xx)*r3                                          00007840
      if (indfa) go to 713                                              00007850
      fa(26)=  e22sxx-2.d0*(e22sxy-e22syy)+3.d0*e12xx-2.d0*e12yy-e12zxx 00007860
     1         -2.d0*e12zyy-(e12txx-e12tyy)*ys2+e11sxx+2.d0*(e11sxy     00007870
     2         +e11syy)+e28txx-e28txy-e27uxy+e27uyy-e18uxy-e18uyy-e17txx00007880
     3         -e17txy+e88sxy+e78xx-e78zyy-e77sxy                       00007890
      fa(27)=  2.d0*(e22sxx-e22sxy)+4.d0*e12zxy-2.d0*(e11sxx+e11sxy)    00007900
     1         +e28txx-e28txy+e27uxx-e27uxy+e18uxx+e18uxy+e17txx+e17txy 00007910
     2         -e88sxx+e88sxy+2.d0*e78zxy+e77sxx+e77sxy                 00007920
      fa(28)= -e12xx-e12zxx-e78xx-e78zxx-e11sxx-e77sxx-e22sxx-e88sxx    00007930
  713 continue                                                          00007940
      fa(29)=  e22zxx+(e22txx-2.d0*e22txy+e22tyy)*ys2-2.d0*e12sxx       00007950
     1         +4.d0*e12sxy-3.d0*e11xx-4.d0*e11xy+2.d0*(e28uxy-e28uyy   00007960
     2         +e27txx-e27txy)+e88zyy+2.d0*e78sxy+e77xx                 00007970
      fa(30)=- 2.d0*(e22zxx+2.d0*e12sxy-e11xx-2.d0*e11xy+e28uxx-e28uxy  00007980
     1         +e27txx-e27txy+e88zxy+e78sxx+e78sxy+e77xx)               00007990
      fa(31)=  e11xx+e77xx+e22zxx+e88zxx+(e12sxx+e78sxx)*2.d0           00008000
      go to 614                                                         00008010
c                                                                       00008020
c                                                                       00008030
c                                                                       00008040
c                                                                       00008050
c        1-st, vector mesons in static limit                            00008060
c        -----------------------------------                            00008070
c                                                                       00008080
c                                                                       00008090
c                                                                       00008100
c                                                                       00008110
  900 mc=1                                                              00008120
c                                                                       00008130
      indfa=.false.                                                     00008140
c                                                                       00008150
c        abbreviations                                                  00008160
c                                                                       00008170
  911 exx=xx                                                            00008180
      e2xx=2.d0*xx                                                      00008190
      e3xx=3.d0*xx                                                      00008200
      e4xx=2.d0*e2xx                                                    00008210
      e2xy=2.d0*xy                                                      00008220
      e4xy=2.d0*e2xy                                                    00008230
      e2yy=2.d0*yy                                                      00008240
      e4yy=2.d0*e2yy                                                    00008250
c                                                                       00008260
c        groundstructure of factors                                     00008270
c                                                                       00008280
      fa( 1)= (exx+e2yy)*r3                                             00008290
      fa( 2)=-e4xy*r3                                                   00008300
      fa( 3)= exx*r3                                                    00008310
      fa( 4)=-exx-e2yy+e2xy                                             00008320
      fa( 5)= e4xy+e2xy-e2xx                                            00008330
      fa( 6)=-e3xx                                                      00008340
      fa( 7)= exx-e2xy                                                  00008350
      fa( 8)= e3xx                                                      00008360
      if (.not.indfa)                                                   00008370
     1fa( 9)=-exx*r3                                                    00008380
      fa(10)=-e4xx-exx-e2yy+e4xy                                        00008390
      fa(11)= e2yy+e4xy-e4xx-e3xx                                       00008400
      fa(12)= e3xx-2.d0*e4xy                                            00008410
      fa(13)= 3.d0*e3xx                                                 00008420
      if (indfa) go to 912                                              00008430
      fa(14)= e3xx+e2xx                                                 00008440
      fa(15)= e4xy                                                      00008450
      fa(16)=-3.d0*e3xx                                                 00008460
      fa(17)=-r3*(exx+e2yy)                                             00008470
      fa(18)= r3*e4xy                                                   00008480
      fa(19)=-r3*exx                                                    00008490
  912 continue                                                          00008500
      fa(20)= exx+e2xy+e2yy                                             00008510
      fa(21)=-e2xx-e4xy-e2xy                                            00008520
      fa(22)= e3xx                                                      00008530
      fa(23)= exx+e2xy                                                  00008540
      fa(24)=-e3xx                                                      00008550
      fa(25)= exx*r3                                                    00008560
      if (indfa) go to 913                                              00008570
      fa(26)= e3xx-e2yy                                                 00008580
      fa(27)= e4xy+e4xy                                                 00008590
      fa(28)=-3.d0*e3xx                                                 00008600
  913 continue                                                          00008610
      fa(29)= e4xy-e3xx                                                 00008620
      fa(30)=-e4xx-e2xx-e4xy                                            00008630
      fa(31)= 3.d0*e3xx                                                 00008640
c                                                                       00008650
c                                                                       00008660
c                                                                       00008670
      do 915 i=10,16                                                    00008680
  915 fa(i)=fa(i)*r33                                                   00008690
      do 916 i=26,31                                                    00008700
  916 fa(i)=fa(i)*r33                                                   00008710
c                                                                       00008720
c        insert (-y) for y                                              00008730
c                                                                       00008740
      xy=-xy                                                            00008750
c                                                                       00008760
      if (indfa) go to 618                                              00008770
c                                                                       00008780
c        factors f                                                      00008790
c                                                                       00008800
      do 917 i=1,16                                                     00008810
      if (i.ge.16) go to 917                                            00008820
      f(i+28)=fa(i+16)                                                  00008830
  917 f(i)=fa(i)                                                        00008840
      indfa=.true.                                                      00008850
      go to 911                                                         00008860
c                                                                       00008870
c                                                                       00008880
c        this has been the end of the contributions of mesons           00008890
c                                                                       00008900
c                                                                       00008910
c                                                                       00008920
c                                                                       00008930
c        errors and warnings                                            00008940
c        -------------------                                            00008950
c                                                                       00008960
c                                                                       00008970
c                                                                       00008980
c                                                                       00008990
 9000 if (indmg(mg)) go to 1995                                         00009000
      write (kwrite,19000) mesong(mg)                                   00009010
19000 format(1h0////'0warning in obdd: meson-group  ',a4,'  does not exi00009020
     1st in this program.'/'0contribution ignored. execution continued.'00009030
     2////)                                                             00009040
      indmg(mg)=.true.                                                  00009050
c                                                                       00009060
c                                                                       00009070
c                                                                       00009080
c                                                                       00009090
 1995 continue                                                          00009100
c                                                                       00009110
c                                                                       00009120
c                                                                       00009130
c                                                                       00009140
 2000 return                                                            00009150
      end                                                               00009160
c**** this package is consistently in double precision *********
c**** July 22, 1993
c
c**** the compiled version of this is in libtm.a
c
c**** this package contains all subroutines needed
c**** namely: tbibdc, delim, obas, cspe;
c**** except the transition potentials which are
c**** contained in the package OBALLTM.F.
c
c
c**** tbibdc.f requires the antisymmetrized obnd
c**** which is contained in oballtm.f
c**** (do not use oball.f)                 5/24/95
c
c
      subroutine tbibdc                                                 00000010
c                                                                       00000020
c                                                                       00000030
c        iterative complex box-diagrams with                            00000040
c          nd and dd intermediate states                                00000050
c                                                                       00000060
c        including delta self energy contributions                      00000070
c          calculated in delta c.m. system                              00000080
c                                                                       00000090
c                                                                       00000100
c                                                                       00000110
c        author:  r. machleidt and ch. elster                           00000120
c                 institut fuer theoretische kernphysik bonn            00000130
c                 nussallee 14-16                                       00000140
c                 d-5300 bonn, w. germany                               00000150
c                 node  : dbnrhrz1                                      00000160
c                 userid: unq306                                        00000170
c                                                                       00000180
c                                                                       00000190
c                                                                       00000200
c        before calling this program for the first time a number of nx1 00000210
c        points have to be stored in the array qx(..) of common/cpts/,  00000220
c        if transition potentials are real or in the array cqx(..) of   00000230
c        common/cptsc/ if transition potentials are complex;            00000240
c        these points - except the point qx(nx1), which gets extra      00000250
c        attention by this program - have to be all arguments xmev,     00000260
c        ymev or cxmev, cymev for which tbibdc will be called for one   00000270
c        quantum number j                                               00000280
c                                                                       00000290
c                                                                       00000300
c                                                                       00000310
      implicit real*8 (a,b,d-h,o-z), complex*16 (c)                     00000320
c
c
      parameter (inx1=97,iimex=41,inr=96,inr1=inr+1)
c
      parameter (lvvnd=16*inx1*inr1,lvvdd=2*lvvnd)
      parameter (irsdd=2*inr1*iimex)
c
c         inx1 : number of momenta from phases/deuter/matter
c         iimex: number of elabs
c         inr  : momenta in loop integration
c                                                                       00000330
      common /crdwrt/ kread,kwrite,kpunch,kda(9)                        00000340
      common /cptsc/  cqx(97)                                           00000350
      common /cpts/   qx(97),dx,nx1,ix,iy                               00000360
c                                                                       00000370
c        common blocks which contain the arguments and values of the    00000380
c        ob-subroutines                                                 00000390
c                                                                       00000400
      common /cpot/   v(6),xmev,ymev                                    00000410
      common /cpotc/  cv(6),cxmev,cymev                                 00000420
      common /cstate/ j,heform,sing,trip,coup,endep,label               00000430
      common /cpoted/ q0qmev,qfmev,pmev,uanspp,wsnspp,ucnspp,udnspp,    00000440
     1                znrl,zrel,smev,noced                              00000450
c                                                                       00000460
c        specifications for these common blocks                         00000470
      logical heform,sing,trip,coup,endep                               00000480
      logical noced                                                     00000490
c                                                                       00000500
c                                                                       00000510
c        common block for all ob-subroutines                            00000520
c                                                                       00000530
      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),  00000540
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,    00000550
     2                ez1,ez2,ct(96),wt(96),                            00000560
     3                ic(10,15),ift(3),mint(3),maxt(3),nt,              00000570
     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),             00000580
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,    00000590
     6                indc(2,15),indpar(3),indxy                        00000600
c                                                                       00000610
c         specifications for this common block                          00000620
c                                                                       00000630
      logical indc,indxy,indpar                                         00000640
      real*8 c,ct                                                       00000650
c                                                                       00000660
c        common block for all obc-subroutines                           00000670
c                                                                       00000680
      common /cobc/   cvj(32,25),cf(52),caa(96),cai(19,5),cfff,cff,     00000690
     1                cx,cxx,cy,cyy,cxy2,cxxpyy,cex,cey,ceem12,         00000700
     2                cez1,cez2,wnni(3),wddi(3)                         00000710
c                                                                       00000720
c                                                                       00000730
c                                                                       00000740
c         common block which contains information of tbibdc             00000750
c                                                                       00000760
c        dimension drsnd(2*(n+1)*(iime-1),2), drsdd(2*(n+1)*(iime-1),3) 00000770
c        dimension of drsnd and drsdd for n=24 and iime=13 (600)        00000780
c                                                                       00000790
      common /ctbsc/  drsnd(irsdd,2),drsdd(irsdd,3),wimnd(irsdd,2),     00000800
     1                wimdd(irsdd,3),cq(inr1),q(inr1),                  00000810
     2                zmevm,zzmevm,zpol,zzpol,wn,wd,wpi,n,issd,         00000820
     3                wdim,wddim,indrnd,indrdd                          00000830
c                                                                       00000840
c         specification for this common block                           00000850
c                                                                       00000860
      logical indrnd,indrdd                                             00000870
      complex*16 drsnd,drsdd,wimnd,wimdd                                00000880
c                                                                       00000890
c                                                                       00000900
c        further specifications                                         00000910
c                                                                       00000920
      dimension cvv(6)                                                  00000930
      dimension cs(inr),cqq(inr1)                                       00000940
      complex*16 vil(32),vill(32)                                       00000950
      dimension qq(inr1),wq(inr),qc(inr),wqc(inr)                       00000960
      dimension tt(2,3,2),itt(8)                                        00000970
      complex*16 u(inr1),uu(inr1)                                       00000980
      real*8 cqi                                                        00000990
c                                                                       00001000
c        dimension vvnd(16*n1*nx1), vvdd(32*n1*nx1)                     00001010
c        dimension of vvnd and vvdd for n=32 and nx1=25                 00001020
      complex*16 vvnd(lvvnd,2),vvdd(lvvdd,2)                            00001030
c                                                                       00001040
c                                                                       00001050
c        dimension of vvnd and vvdd for n=24 and nx1=25                 00001060
c     complex*16 vvnd(10000,2),vvdd(20000,2)                            00001070
c                                                                       00001080
c        dimension of vvnd and vvdd for n=32 and nx1=33                 00001090
c     complex*16 vvnd(17424,2),vvdd(34848,2)                            00001100
c                                                                       00001110
c                                                                       00001120
c                                                                       00001130
      integer nj(20)                                                    00001140
      data nj/20*0/
      data jj/-1/
      integer nname(17),name(3)                                         00001150
      logical index
      logical swed,swd,swdd
      logical indnd,inddd
      logical realnd,realdd
      logical indsnd,indsdd
      logical indwnd,indwdd
      logical indmnd,indmdd
c                                                                       00001240
      logical indvd
      logical indbet
      logical indcxp
      logical indsel
      logical edint
      logical indpnd,indpdd                                             00001300
c                                                                       00001310
      logical indtj
c
c
      data index/.false./
      data swed/.false./,swd/.false./,swdd/.false./                     00001180
      data indnd/.false./,inddd/.false./                                00001190
      data realnd/.true./,realdd/.true./                                00001200
      data indsnd/.false./,indsdd/.false./                              00001210
      data indwnd/.false./,indwdd/.false./                              00001220
      data indmnd/.false./,indmdd/.false./                              00001230
c                                                                       00001240
      data indvd/.false./                                               00001250
      data indbet/.false./                                              00001260
      data indcxp/.false./                                              00001270
      data indsel/.true./                                               00001280
      data edint/.false./                                               00001290
c                                                                       00001310
      data indtj/.false./                                               00001320
c                                                                       00001330
      data q0mev/-1.d0/                                                 00001340
      data pih/1.570796326794897d0/                                     00001350
      data pi/3.141592653589793d0/                                      00001360
      data rr/0.81d0/,gamma/71.d0/                                      00001370
c
c
c**** statement functions ***********
      cmplx(x,y)=dcmplx(x,y)
      float(i)=dfloat(i)
      abs(a)=dabs(a)
      cabs(ca)=zabs(ca)
      sqrt(a)=dsqrt(a)
      cos(a)=dcos(a)
      tan(a)=dtan(a)
      log(a)=dlog(a)
c                                                                       00001380
c                                                                       00001390
c                                                                       00001400
c                                                                       00001410
10000 format (17a4)                                                     00001420
10001 format (1h0,17a4)                                                 00001430
10002 format (2a4,a2,20i3)                                              00001440
10003 format (1h0,2a4,a2,20i3)                                          00001450
10004 format (2a4,a2,6f10.4)                                            00001460
10005 format (1h0,2a4,a2,6f10.4)                                        00001470
10010 format(1h1//' tbibdc: two-boson-exchange - iterative box diagrams'00001480
     1/1x,51(1h-)/18x,'with nd and dd intermediate states'/             00001490
     2 18x,34(1h-)/'0input-parameter-set:'/1h ,20(1h-)/1h0,15a4)        00001500
c                                                                       00001510
c                                                                       00001520
c                                                                       00001530
c                                                                       00001540
      if (index) go to 50                                               00001550
      index=.true.                                                      00001560
c                                                                       00001570
      write (kwrite,10010)                                              00001580
      read  (kread ,10000) nname                                        00001590
      write (kwrite,10001) nname                                        00001600
c                                                                       00001610
c        iprnd=0 no nd box,                                             00001620
c        iprnd=1 non-rel. propagator in nd box-diagram,                 00001630
c        iprnd=2 relativistic propagator in nd box-diagram;             00001640
c        ispnd: parameter for single particle energy of nucleons in nd; 00001650
c        ispd:  parameter for single particle energy of deltas in nd;   00001660
c                                                                       00001670
      read  (kread ,10002) name,iprnd,ispnd,ispd                        00001680
      write (kwrite,10003) name,iprnd,ispnd,ispd                        00001690
c                                                                       00001700
c        iprdd, ispndd, ispdd for the dd box-diagram analogous to nd    00001710
c                                                                       00001720
      read  (kread ,10002) name,iprdd,ispndd,ispdd                      00001730
      write (kwrite,10003) name,iprdd,ispndd,ispdd                      00001740
c                                                                       00001750
c        isnd=1 : self energy corrections in nd box diagram             00001760
c        ipmnd=1: propagator modification                               00001770
c                 for static obe-propagators ipmnd=0                    00001780
c        ipwnd=1: delta width                                           00001790
c                                                                       00001800
      read  (kread ,10002) name,isnd,ipmnd,ipwnd                        00001810
      write (kwrite,10003) name,isnd,ipmnd,ipwnd                        00001820
c                                                                       00001830
c        isdd, ipmdd, ipwdd for the dd box-diagram analogous to nd      00001840
c                                                                       00001850
      read  (kread ,10002) name,isdd,ipmdd,ipwdd                        00001860
      write (kwrite,10003) name,isdd,ipmdd,ipwdd                        00001870
c                                                                       00001880
c        itt=0: iso-spin factor for boxes only                          00001890
c        itt=1: iso-scalar iso-spin factor which includes               00001900
c                   crossed boxes in case of pi-rho                     00001910
c                                                                       00001920
      read  (kread ,10002) name,ida,ide,(itt(i),i=1,4)                  00001930
      write (kwrite,10003) name,ida,ide,(itt(i),i=1,4)                  00001940
      ide1=ide+1                                                        00001950
      ide4=ide+4                                                        00001960
      read  (kread ,10002) name,idda,idde,(itt(i),i=ide1,ide4)          00001970
      write (kwrite,10003) name,idda,idde,(itt(i),i=ide1,ide4)          00001980
      read  (kread ,10002) name,(nj(j1),j1=1,20)                        00001990
      write (kwrite,10003) name,(nj(j1),j1=1,20)                        00002000
      read  (kread ,10004) name,cqi                                     00002010
      write (kwrite,10005) name,cqi                                     00002020
      read  (kread ,10004) name,wn                                      00002030
      write (kwrite,10005) name,wn                                      00002040
      read  (kread ,10004) name,wd                                      00002050
      write (kwrite,10005) name,wd                                      00002060
      read  (kread ,10004) name,wpi                                     00002070
      write (kwrite,10005) name,wpi                                     00002080
c                                                                       00002090
c        iprnd/dd=1 :  icmplx=0,  arguments  xmev, ymev  for obnd/dd    00002100
c        iprnd/dd=2 :  icmplx=1,  arguments cxmev, cymev for obndc/ddc  00002110
c                                                                       00002120
      read  (kread ,10002) name,icmplx                                  00002130
      write (kwrite,10003) name,icmplx                                  00002140
c                                                                       00002150
c        parameter for complex contour deformation in                   00002160
c        case of iprnd/dd=2                                             00002170
c                                                                       00002180
      read  (kread ,10004) name,beta                                    00002190
      write (kwrite,10005) name,beta                                    00002200
c                                                                       00002210
c         complex delta mass:                                           00002220
c        itj=0: no complex delta mass                                   00002230
c         parametrization of bransden-moorhouse                         00002240
c        itj=1: complex delta mass in box propagator only               00002250
c        itj=2: complex delta mass in box- and potential propagators    00002260
c         imaginary part of delta self energy diagram                   00002270
c        itj=4: thresholds  nd: 2.*wn+wpi,  dd: 2.*wn+2.*wpi            00002280
c                                                                       00002290
      read  (kread ,10002) name,itj                                     00002300
      write (kwrite,10003) name,itj                                     00002310
c                                                                       00002320
c        wd  is the mass of the delta                                   00002330
c        wpi is the mass of the pion                                    00002340
c                                                                       00002350
      if (iprnd.ne.0) indnd=.true.                                      00002360
      if (iprdd.ne.0) inddd=.true.                                      00002370
      if (.not.indnd) ide=0                                             00002400
      if (.not.inddd) idde=0                                            00002410
      if (beta.ne.0.d0) indbet=.true.                                   00002420
      if (indbet) realnd=.false.
      if (indbet) realdd=.false.
      if (icmplx.ne.0) indcxp=.true.                                    00002430
      if (itj.ne.0) indtj=.true.                                        00002440
      idaa=ide+idda                                                     00002450
      idee=ide+idde                                                     00002460
      do 10 i=1,8                                                       00002470
   10 itt(i)=itt(i)+1                                                   00002480
c                                                                       00002490
      if (isnd.ne.0) indsnd=.true.                                      00002500
      if (isdd.ne.0) indsdd=.true.                                      00002510
      if (indsnd.or.indsdd) indsel=.false.                              00002520
      if (ipwnd.ne.0) indwnd=.true.                                     00002530
      if (ipwdd.ne.0) indwdd=.true.                                     00002540
      if (ipmnd.ne.0) indmnd=.true.                                     00002550
      if (ipmdd.ne.0) indmdd=.true.                                     00002560
      indrnd=.false.                                                    00002570
      indrdd=.false.                                                    00002580
c                                                                       00002590
c                                                                       00002600
c        get meson parameters                                           00002610
c                                                                       00002620
      inter=1                                                           00002630
      if (indpar(inter)) go to 15                                       00002640
c                                                                       00002650
      call obpar                                                        00002660
c                                                                       00002670
   15 inter=2                                                           00002680
      if (indpar(inter)) go to 20                                       00002690
c                                                                       00002700
      call obpar                                                        00002710
c                                                                       00002720
   20 inter=3                                                           00002730
      if (indpar(inter)) go to 25                                       00002740
c                                                                       00002750
      call obpar                                                        00002760
c                                                                       00002770
c                                                                       00002780
c                                                                       00002790
   25 dwn=1.d0/wn                                                       00002800
      wnq=wn*wn                                                         00002810
      wd2=2.d0*wd                                                       00002820
      wdq=wd*wd                                                         00002830
      wdmn=wd-wn                                                        00002840
      wdpn=wd+wn                                                        00002850
      wdtn=wd*wn                                                        00002860
      wdtnq=wdq*wnq                                                     00002870
      wdmn2=2.d0*wdmn                                                   00002880
c                                                                       00002890
      dwpi=1.d0/wpi                                                     00002900
      wpimn=wpi-wn                                                      00002910
      wpipn=wpi+wn                                                      00002920
      wpimnq=wpimn*wpimn                                                00002930
      wpipnq=wpipn*wpipn                                                00002940
c                                                                       00002950
      wnim=0.d0                                                         00002960
      wdim=0.d0                                                         00002970
      wddim=0.d0                                                        00002980
      wddi(2)=0.d0                                                      00002990
      wddi(3)=0.d0                                                      00003000
c                                                                       00003010
      thresh=wn+wn+wpi                                                  00003020
c
c
c
c
      nii=0
c
      do 35 ii=1,3
      do 35 i=1,irsdd
      wimdd(i,ii)=0.d0
      if (ii.gt.2) go to 35
      wimnd(i,ii)=0.d0
   35 continue
c                                                                       00003030
c                                                                       00003040
c                                                                       00003050
c        iso-spin factors                                               00003060
c                                                                       00003070
c                                                                       00003080
c        iso-spin factors for boxes only                                00003090
c         nd-boxes                                                      00003100
c         t=1                                                           00003110
      tt(1,2,1)=8.d0/3.d0                                               00003120
c         t=0                                                           00003130
      tt(2,2,1)=0.d0                                                    00003140
c         dd-boxes                                                      00003150
c         t=1                                                           00003160
      tt(1,3,1)=10.d0/9.d0                                              00003170
c         t=0                                                           00003180
      tt(2,3,1)=2.d0                                                    00003190
c                                                                       00003200
c        iso-scalar iso-spin factors which include crossed boxes        00003210
c                                                  in case of pi-rho    00003220
c                                                                       00003230
c         nd-boxes                                                      00003240
      tt(1,2,2)=4.d0                                                    00003250
      tt(2,2,2)=4.d0                                                    00003260
c         dd-boxes                                                      00003270
      tt(1,3,2)=8.d0/3.d0                                               00003280
      tt(2,3,2)=8.d0/3.d0                                               00003290
c                                                                       00003300
c                                                                       00003310
c                                                                       00003320
c                                                                       00003330
c                                                                       00003340
c                                                                       00003350
   50 if (j.eq.jj) go to 55                                             00003360
      jj=j                                                              00003370
      issd=-1                                                           00003380
      indvd=.false.                                                     00003390
      j1=j+1                                                            00003400
      aj= float(j)                                                      00003410
      aj1= float(j+1)                                                   00003420
      d2j1=1.d0/ float(2*j+1)                                           00003430
      arjj1= sqrt(aj*aj1)                                               00003440
c
      n=96
      if (j1.gt.20) go to 51
      if (nj(j1).eq.0) nj(j1)=nj(j)                                     00003450
      n=nj(j1)                                                          00003460
   51 continue
      n1=n+1                                                            00003470
c                                                                       00003480
c                                                                       00003490
   55 if (q0mev.eq.qx(nx1)) go to 200                                   00003500
c                                                                       00003510
c                                                                       00003520
c        get gauss points and weights                                   00003530
c                                                                       00003540
c                                                                       00003550
      call gset (0.d0,1.d0,n,q,wq)                                      00003560
c                                                                       00003570
c        transform gauss points and weights                             00003580
c                                                                       00003590
      do 60 i=1,n                                                       00003600
      qaa=pih*q(i)                                                      00003610
      q(i)= tan(qaa)*cqi                                                00003620
      qq(i)=q(i)*q(i)                                                   00003630
      qaa=1.d0/ cos(qaa)                                                00003640
   60 wq(i)=pih*cqi*qaa*qaa*wq(i)                                       00003650
      q(n1)=0.d0                                                        00003660
      cq(n1)=(0.d0,0.d0)                                                00003670
c                                                                       00003680
      q0mev=qx(nx1)                                                     00003690
c                                                                       00003700
      if (zrel.lt.thresh) go to 105                                     00003710
      if (realnd.and.realdd) go to 105                                  00003720
c                                                                       00003730
c                                                                       00003740
c        prepare integration in complex plane                           00003750
c                                                                       00003760
c                                                                       00003770
      naa=n/2                                                           00003780
      nab=n-naa                                                         00003790
c                                                                       00003800
      call gset (0.d0,q0mev,naa,q,wq)                                   00003810
      call gset (0.d0,1.d0,nab,qc,wqc)                                  00003820
c                                                                       00003830
      do 100 i=1,nab                                                    00003840
      xx=pih*qc(i)                                                      00003850
      dcc=1.d0/ cos(xx)                                                 00003860
      q(i+naa)=cqi* tan(xx)+q0mev                                       00003870
  100 wq(i+naa)=pih*cqi*wqc(i)*dcc*dcc                                  00003880
  105 continue                                                          00003890
c                                                                       00003900
c                                                                       00003910
      do 150 i=1,n                                                      00003920
      if (zrel.ge.thresh.and.indbet.and.q(i).lt.q0mev) go to 140        00003930
      qre=q(i)                                                          00003940
      qim=0.d0                                                          00003950
      sre=wq(i)                                                         00003960
      sim=0.d0                                                          00003970
      go to 145                                                         00003980
c                                                                       00003990
  140 qdq0=q(i)/q0mev                                                   00004000
      qre=q(i)*(2.d0-qdq0)                                              00004010
      qim=beta*qdq0*(q(i)-q0mev)                                        00004020
      sre=2.d0*(1.d0-qdq0)*wq(i)                                        00004030
      sim=beta*(2.d0*qdq0-1.d0)*wq(i)                                   00004040
c                                                                       00004050
  145 cq(i)= cmplx(qre,qim)                                             00004060
      cs(i)= cmplx(sre,sim)                                             00004070
  150 cqq(i)=cq(i)*cq(i)                                                00004080
c                                                                       00004090
c                                                                       00004100
      if (.not.indtj) go to 205                                         00004110
      if (itj.ge.3) go to 180                                           00004120
c                                                                       00004130
c                                                                       00004140
c         complex delta mass: parametrization of bransden-moorhouse     00004150
c                                                                       00004160
c                                                                       00004170
      spndr=zrel-wn                                                     00004180
c                                                                       00004190
c        nd                                                             00004200
c                                                                       00004210
      spnd=spndr*spndr                                                  00004220
      sqndq=(spnd-wpimnq)*(spnd-wpipnq)*0.25d0/spnd                     00004230
c                                                                       00004240
      if (sqndq.gt.0.d0) go to 170                                      00004250
      wdim=0.d0                                                         00004260
      wddim=0.d0                                                        00004270
      go to 175                                                         00004280
c                                                                       00004290
  170 sqnd=dwpi*rr* sqrt(sqndq)                                         00004300
      sqndq=sqnd*sqnd                                                   00004310
      wdim=-gamma*sqndq*sqnd/(1.d0+sqndq)                               00004320
c                                                                       00004330
c        dd                                                             00004340
c                                                                       00004350
      spndr=spndr-wpi                                                   00004360
      spnd=spndr*spndr                                                  00004370
      sqddq=(spnd-wpimnq)*(spnd-wpipnq)*0.25d0/spnd                     00004380
c                                                                       00004390
      if (sqddq.le.0.d0) go to 175                                      00004400
      sqdd=dwpi*rr* sqrt(sqddq)                                         00004410
      sqddq=sqdd*sqdd                                                   00004420
      wddim=-gamma*sqddq*sqdd/(1.d0+sqddq)                              00004430
c                                                                       00004440
  175 if (itj.ne.2) go to 205                                           00004450
      wddi(2)=wdim                                                      00004460
      wddi(3)=wddim                                                     00004470
      go to 205                                                         00004480
c                                                                       00004490
c                                                                       00004500
c        complex delta mass: imaginary part of delta self energy        00004510
c                            diagram  ( delta c.m.system, i.e. q=0 )    00004520
c                                                                       00004530
c                                                                       00004540
  180 call delim (itj)                                                  00004550
c
c                                                                       00004560
      go to 205                                                         00004570
c                                                                       00004580
c                                                                       00004590
c                                                                       00004600
c                                                                       00004610
c                                                                       00004620
c        compute propagators                                            00004630
c        -------------------                                            00004640
c                                                                       00004650
c                                                                       00004660
c                                                                       00004670
c                                                                       00004680
  200 if (noced) go to 1000                                             00004690
      if (edint) go to 203                                              00004700
      noced=.true.                                                      00004710
      go to 205                                                         00004720
  203 indvd=.false.                                                     00004730
  205 continue                                                          00004740
c**** write (6,25000) wdim,wddim
25000 format (' wdim, wddim',2d20.10)
      qfmevq=qfmev*qfmev                                                00004750
c                                                                       00004760
c                                                                       00004770
c                                                                       00004780
c        prepare starting energies for nd                               00004790
c                                                                       00004800
c                                                                       00004810
      if (.not.indnd) go to 300                                         00004820
      pmevq=0.d0                                                        00004830
      pq=0.d0                                                           00004840
      go to (221,222),iprnd                                             00004850
  221 zmevm=q0qmev*dwn                                                  00004860
      go to 250                                                         00004870
  222 zmevm=2.d0* sqrt(wnq+q0qmev)                                      00004880
c                                                                       00004890
c                                                                       00004900
c        prepare pole of propagator for nd                              00004910
c                                                                       00004920
  250 zpls=0.d0                                                         00004930
      czpln=(0.d0,0.d0)                                                 00004940
      zpol=0.d0                                                         00004950
      indpnd=.true.                                                     00004960
      if (indtj.or.indwnd) go to 300                                    00004970
c                                                                       00004980
      go to (260,270),iprnd                                             00004990
c                                                                       00005000
  260 if (zmevm.lt.wdmn) go to 300                                      00005010
c                                                                       00005020
      zp1=zmevm-wdmn                                                    00005030
      zplq=2.d0*zp1*wdtn/wdpn                                           00005040
      zpol= sqrt(zplq)                                                  00005050
      enzp=spe(pq*wnq+zplq,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,iprnd, 00005060
     1 ispnd,0)                                                         00005070
      edzp=spe(pq*wdq+zplq,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wd,iprnd, 00005080
     1 ispd,0)                                                          00005090
      zpls=zpol*(zmevm+enzp+edzp)                                       00005100
      zpln1=zpol*wdtn/wdpn*pi                                           00005110
      czpln= cmplx(0.d0,-zpln1)                                         00005120
      indpnd=.false.                                                    00005130
      go to 300                                                         00005140
c                                                                       00005150
c                                                                       00005160
  270 if (zmevm.lt.wdpn) go to 300                                      00005170
c                                                                       00005180
      dzmev=1.d0/zmevm                                                  00005190
      zp1=zmevm*zmevm-wdq-wnq                                           00005200
      zplq=(zp1*zp1-4.d0*wdtnq)*dzmev*dzmev*0.25d0                      00005210
      zpol= sqrt(zplq)                                                  00005220
      if (zrel.ge.thresh.and.indbet.and.zpol.lt.q0mev) go to 300        00005230
      enzp=spe(pq*wnq+zplq,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,iprnd, 00005240
     1 ispnd,0)                                                         00005250
      edzp=spe(pq*wdq+zplq,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wd,iprnd, 00005260
     1 ispd,0)                                                          00005270
      zpls=zpol*enzp*edzp                                               00005280
      zpln1= abs(zmevm/wdpn-1.d0)                                       00005290
      czpln=zpls*dzmev* cmplx( log(zpln1),-pi)                          00005300
      indpnd=.false.                                                    00005310
c                                                                       00005320
c                                                                       00005330
c        prepare starting energies for dd                               00005340
c                                                                       00005350
  300 if (.not.inddd) go to 400                                         00005360
      ppmevq=0.d0                                                       00005370
      go to (321,322),iprdd                                             00005380
  321 zzmevm=q0qmev*dwn                                                 00005390
      go to 350                                                         00005400
  322 zzmevm=2.d0* sqrt(wnq+q0qmev)                                     00005410
c                                                                       00005420
c                                                                       00005430
c        prepare pole of propagator for dd                              00005440
c                                                                       00005450
  350 zzpls=0.d0                                                        00005460
      czzpln=(0.d0,0.d0)                                                00005470
      zzpol=0.d0                                                        00005480
      indpdd=.true.                                                     00005490
      if (indtj.or.indwdd) go to 400                                    00005500
c                                                                       00005510
      go to (360,370),iprdd                                             00005520
c                                                                       00005530
  360 if (zzmevm.lt.wdmn2) go to 400                                    00005540
c                                                                       00005550
      zzplq=wd*(zzmevm-wdmn2)                                           00005560
      zzpol= sqrt(zzplq)                                                00005570
      edzzp=spe(ppmevq+zzplq,qfmevq,uanssp,wsnspp,ucnspp,udnspp,wd,iprdd00005580
     1,ispdd,0)                                                         00005590
      zzpls=zzpol*(zzmevm+2.d0*edzzp)                                   00005600
      zzpln=pi*zzpol*0.5d0*wd                                           00005610
      czzpln= cmplx(0.d0,-zzpln)                                        00005620
      indpdd=.false.                                                    00005630
      go to 400                                                         00005640
c                                                                       00005650
c                                                                       00005660
  370 if (zzmevm.lt.wd2) go to 400                                      00005670
c                                                                       00005680
      zzplq=0.25d0*zzmevm*zzmevm-wdq                                    00005690
      zzpol= sqrt(zzplq)                                                00005700
      if (zrel.ge.thresh.and.indbet.and.zzpol.lt.q0mev) go to 400       00005710
      edzzp=spe(ppmevq+zzplq,qfmevq,uanssp,wsnspp,ucnspp,udnspp,wd,iprdd00005720
     1,ispdd,0)                                                         00005730
       zzpls=zzplq*(2.d0*edzzp+zzmevm)                                  00005740
       czzpln= cmplx(0.d0,-pi*zzplq)                                    00005750
       indpdd=.false.                                                   00005760
c                                                                       00005770
c                                                                       00005780
c                                                                       00005790
  400 cuq0=(0.d0,0.d0)                                                  00005800
      cuuq0=(0.d0,0.d0)                                                 00005810
      u(n1)=(0.d0,0.d0)                                                 00005820
      uu(n1)=(0.d0,0.d0)                                                00005830
      cwnd=(0.d0,0.d0)                                                  00005840
      cwdd=(0.d0,0.d0)                                                  00005850
c                                                                       00005860
c                                                                       00005870
c                                                                       00005880
c                                                                       00005890
c                                                                       00005900
c        compute self energy corrections for nd and dd                  00005910
c        ---------------------------------------------                  00005920
c                                                                       00005930
c                                                                       00005940
c                                                                       00005950
c                                                                       00005960
      if (indsel) go to 500                                             00005970
c                                                                       00005980
      issd=issd+1                                                       00005990
      nii=issd*n1                                                       00006000
c                                                                       00006010
      if (indrnd.or.indrdd) go to 500                                   00006020
c                                                                       00006030
c**** this call is disabled ************
c**** call deldsd                                                       00006040
      write (6,20000)
20000 format (' Watch it: you are using a part of TBIBDC that is
     1 disabled.')
c     -----------                                                       00006050
c                                                                       00006060
c                                                                       00006070
c                                                                       00006080
c                                                                       00006090
  500 do 595 i=1,n                                                      00006100
c                                                                       00006110
      iiid=nii+i                                                        00006120
c                                                                       00006130
      if (.not.indnd) go to 550                                         00006140
c                                                                       00006150
c                                                                       00006160
c                                                                       00006170
c        nd                                                             00006180
c                                                                       00006190
      cen=cspe(cqq(i),wn,wnim,iprnd,ispnd)                              00006200
      ced=cspe(cqq(i),wd,wdim,iprnd,ispd)                               00006210
      csx=cen+ced-zmevm                                                 00006220
      csxx=cs(i)*cqq(i)                                                 00006230
      if (.not.indwnd) go to 510                                        00006240
      cwnd=wimnd(iiid,1)*(0.d0,1.d0)                                    00006250
  510 u(i)=csxx/(csx-cwnd)                                              00006260
      if (indwnd) go to 550                                             00006270
      csdq=cs(i)/csx                                                    00006280
      cseq=1.d0/(cen*ced)                                               00006290
      if (iprnd.eq.1) cseq=1.d0/(cen+ced+zmevm)                         00006300
      cuq0=cuq0+cq(i)*csdq*cseq                                         00006310
c                                                                       00006320
c                                                                       00006330
c                                                                       00006340
c        dd                                                             00006350
c                                                                       00006360
  550 if (.not.inddd) go to 595                                         00006370
      ced=cspe(cqq(i),wd,wddim,iprdd,ispdd)                             00006380
      csx=2.d0*ced-zzmevm                                               00006390
      csxx=cs(i)*cqq(i)                                                 00006400
      if (.not.indwdd) go to 560                                        00006410
      cwdd=wimdd(iiid,1)*(0.d0,1.d0)                                    00006420
  560 uu(i)=csxx/(csx-cwdd)                                             00006430
      if (indwdd) go to 595                                             00006440
      csdq=cs(i)/csx                                                    00006450
      cseq=1.d0/(2.d0*ced+zzmevm)                                       00006460
      if (iprdd.eq.1) cseq=cseq*cq(i)                                   00006470
      cuuq0=cuuq0+csdq*cseq                                             00006480
c                                                                       00006490
c                                                                       00006500
c                                                                       00006510
  595 continue                                                          00006520
c                                                                       00006530
      if (indpnd) go to 610                                             00006540
      cuq0=-cuq0*zpls                                                   00006550
      u(n1)=cuq0-czpln                                                  00006560
      indvd=.false.                                                     00006570
c                                                                       00006580
  610 if (indpdd) go to 800                                             00006590
      cuuq0=-cuuq0*zzpls                                                00006600
      uu(n1)=cuuq0-czzpln                                               00006610
      indvd=.false.                                                     00006620
c                                                                       00006630
c                                                                       00006640
c        this has been the end of the computation of the energy denom.  00006650
c                                                                       00006660
c                                                                       00006670
c                                                                       00006680
c                                                                       00006690
c                                                                       00006700
  800 if (indsel) go to 1000                                            00006710
c                                                                       00006720
c                                                                       00006730
c                                                                       00006740
c        multiply propagator with dressing factor                       00006750
c        ----------------------------------------                       00006760
c                                                                       00006770
c                                                                       00006780
c                                                                       00006790
c        case of nd                                                     00006800
c                                                                       00006810
c                                                                       00006820
      if (.not.indsnd) go to 900                                        00006830
      if (.not.indmnd) go to 900                                        00006840
c                                                                       00006850
c                                                                       00006860
      do 850 i=1,n1                                                     00006870
      iiid=nii+i                                                        00006880
  850 u(i)=u(i)*drsnd(iiid,1)                                           00006890
c                                                                       00006900
c                                                                       00006910
c                                                                       00006920
c        case of dd                                                     00006930
c                                                                       00006940
c                                                                       00006950
  900 if (.not.indsdd) go to 1000                                       00006960
      if (.not.indmdd) go to 1000                                       00006970
c                                                                       00006980
c                                                                       00006990
      do 950 i=1,n1                                                     00007000
      iiid=nii+i                                                        00007010
  950 uu(i)=uu(i)*drsdd(iiid,1)                                         00007020
c                                                                       00007030
c                                                                       00007040
c                                                                       00007050
c                                                                       00007060
c                                                                       00007070
c                                                                       00007080
c        call nd and dd - transition-potentials and store them          00007090
c        -----------------------------------------------------          00007100
c                                                                       00007110
c                                                                       00007120
c                                                                       00007130
c                                                                       00007140
 1000 if (indvd) go to 1005                                             00007150
      nx1a=1                                                            00007160
      go to 1010                                                        00007170
c                                                                       00007180
 1005 if (ix.ne.nx1.and.iy.ne.nx1) go to 2000                           00007190
      if (q0meva.eq.qx(nx1)) go to 2000                                 00007200
      nx1a=nx1                                                          00007210
c                                                                       00007220
c                                                                       00007230
 1010 if (indcxp) go to 1500                                            00007240
c                                                                       00007250
c                                                                       00007260
c                                                                       00007270
c                                                                       00007280
c        transition potentials real                                     00007290
c        --------------------------                                     00007300
c                                                                       00007310
c                                                                       00007320
c                                                                       00007330
c        store arguments                                                00007340
c                                                                       00007350
      xm=xmev                                                           00007360
      ym=ymev                                                           00007370
      q0meva=qx(nx1)                                                    00007380
c                                                                       00007390
      do 1395 k=nx1a,nx1                                                00007400
      ik=n1*(k-1)                                                       00007410
      xmev=qx(k)                                                        00007420
c                                                                       00007430
c                                                                       00007440
      do 1395 i=1,n1                                                    00007450
      im1=i-1                                                           00007460
      ymev=q(i)                                                         00007470
c                                                                       00007480
c                                                                       00007490
c        obnd                                                           00007500
c                                                                       00007510
      if (.not.indnd) go to 1300                                        00007520
      if (i.eq.n1) ymev=zpol                                            00007530
c                                                                       00007540
      call obnd                                                         00007550
c                                                                       00007560
      ik1=16*(im1+ik)                                                   00007570
      if (swd) go to 1210                                               00007580
      swd=.true.                                                        00007590
      imad=imaa(2)                                                      00007600
      imed=imea(2)                                                      00007610
      ld=imad-1                                                         00007620
      imed=imed-imad+1                                                  00007630
 1210 do 1225 il=1,imed                                                 00007640
      iil=il+ld                                                         00007650
      do 1225 ii=1,16                                                   00007660
      if ( abs(vj(ii,iil)).lt.1.d-30) vj(ii,iil)=0.d0                   00007670
 1225 vvnd(ii+ik1,il)= cmplx(vj(ii,iil),0.d0)                           00007680
c                                                                       00007690
c                                                                       00007700
c        obdd                                                           00007710
c                                                                       00007720
 1300 if (.not.inddd) go to 1395                                        00007730
      if (i.eq.n1) ymev=zzpol                                           00007740
c                                                                       00007750
      call obdd                                                         00007760
c                                                                       00007770
      ik1=32*(im1+ik)                                                   00007780
      if (swdd) go to 1310                                              00007790
      swdd=.true.                                                       00007800
      imadd=imaa(3)                                                     00007810
      imedd=imea(3)                                                     00007820
      ldd=imadd-1                                                       00007830
      imedd=imedd-imadd+1                                               00007840
 1310 do 1325 il=1,imedd                                                00007850
      iil=il+ldd                                                        00007860
      do 1325 ii=1,32                                                   00007870
      if ( abs(vj(ii,iil)).lt.1.d-30) vj(ii,iil)=0.d0                   00007880
 1325 vvdd(ii+ik1,il)= cmplx(vj(ii,iil),0.d0)                           00007890
c                                                                       00007900
 1395 continue                                                          00007910
      indvd=.true.                                                      00007920
c                                                                       00007930
c        get information about endep of transition potentials           00007940
c                                                                       00007950
      if (swed) go to 1405                                              00007960
      edint=endep                                                       00007970
      swed=.true.                                                       00007980
c                                                                       00007990
c        set endep according to the fact that tbibd is endep            00008000
      endep=.true.                                                      00008010
 1405 continue                                                          00008020
c                                                                       00008030
c        restore arguments                                              00008040
c                                                                       00008050
      xmev=xm                                                           00008060
      ymev=ym                                                           00008070
      go to 2000                                                        00008080
c                                                                       00008090
c                                                                       00008100
c                                                                       00008110
c                                                                       00008120
c        transition potentials complex                                  00008130
c        -----------------------------                                  00008140
c                                                                       00008150
c                                                                       00008160
c                                                                       00008170
c        store arguments                                                00008180
c                                                                       00008190
 1500 cxm=cxmev                                                         00008200
      cym=cymev                                                         00008210
      q0meva=qx(nx1)                                                    00008220
c                                                                       00008230
      do 1795 k=nx1a,nx1                                                00008240
      ik=n1*(k-1)                                                       00008250
      cxmev=cqx(k)                                                      00008260
c                                                                       00008270
c                                                                       00008280
      do 1795 i=1,n1                                                    00008290
      im1=i-1                                                           00008300
      cymev=cq(i)                                                       00008310
c                                                                       00008320
c                                                                       00008330
c        obndc                                                          00008340
c                                                                       00008350
      if (.not.indnd) go to 1700                                        00008360
      if (i.eq.n1) cymev= cmplx(zpol,0.d0)                              00008370
c                                                                       00008380
c**** this call is diabled ********
c**** call obndc                                                        00008390
      write (6,20000)
c                                                                       00008400
      ik1=16*(im1+ik)                                                   00008410
      if (swd) go to 1610                                               00008420
      swd=.true.                                                        00008430
      imad=imaa(2)                                                      00008440
      imed=imea(2)                                                      00008450
      ld=imad-1                                                         00008460
      imed=imed-imad+1                                                  00008470
 1610 do 1625 il=1,imed                                                 00008480
      iil=il+ld                                                         00008490
      do 1625 ii=1,16                                                   00008500
c**** if-statement disabled **********
c**** if ( cabs(cvj(ii,iil)).lt.1.d-20) cvj(ii,iil)=(0.d0,0.d0)         00008510
 1625 vvnd(ii+ik1,il)=cvj(ii,iil)                                       00008520
c                                                                       00008530
c                                                                       00008540
c        obddc                                                          00008550
c                                                                       00008560
 1700 if (.not.inddd) go to 1795                                        00008570
      if (i.eq.n1) cymev= cmplx(zzpol,0.d0)                             00008580
c                                                                       00008590
c**** this call is disabled ************
c**** call obddc                                                        00008600
      write (6,20000)
c                                                                       00008610
      ik1=32*(im1+ik)                                                   00008620
      if (swdd) go to 1710                                              00008630
      swdd=.true.                                                       00008640
      imadd=imaa(3)                                                     00008650
      imedd=imea(3)                                                     00008660
      ldd=imadd-1                                                       00008670
      imedd=imedd-imadd+1                                               00008680
 1710 do 1725 il=1,imedd                                                00008690
      iil=il+ldd                                                        00008700
      do 1725 ii=1,32                                                   00008710
c**** if-statement below is disabled ***********
c**** if ( cabs(cvj(ii,iil)).lt.1.d-20) cvj(ii,iil)=(0.d0,0.d0)         00008720
 1725 vvdd(ii+ik1,il)=cvj(ii,iil)                                       00008730
c                                                                       00008740
 1795 continue                                                          00008750
      indvd=.true.                                                      00008760
c                                                                       00008770
c        get information about endep of transition potentials           00008780
c                                                                       00008790
      if (swed) go to 1805                                              00008800
      edint=endep                                                       00008810
      swed=.true.                                                       00008820
c                                                                       00008830
c        set endep according to the fact that tbibd is endep            00008840
      endep=.true.                                                      00008850
 1805 continue                                                          00008860
c                                                                       00008870
c        restore arguments                                              00008880
c                                                                       00008890
      cxmev=cxm                                                         00008900
      cymev=cym                                                         00008910
c                                                                       00008920
c                                                                       00008930
c                                                                       00008940
c                                                                       00008950
c                                                                       00008960
c                                                                       00008970
c        compute box-diagrams                                           00008980
c        --------------------                                           00008990
c                                                                       00009000
c                                                                       00009010
c                                                                       00009020
c                                                                       00009030
 2000 do 2003 iv=1,6                                                    00009040
 2003 cv(iv)=(0.d0,0.d0)                                                00009050
      do 2005 id=1,idee                                                 00009060
      do 2005 iv=1,32                                                   00009070
 2005 cvj(iv,id)=(0.d0,0.d0)                                            00009080
c                                                                       00009090
c                                                                       00009100
c                                                                       00009110
c                                                                       00009120
c        loop of integration                                            00009130
c                                                                       00009140
c                                                                       00009150
c                                                                       00009160
c                                                                       00009170
      ixx=n1*(ix-1)                                                     00009180
      iyy=n1*(iy-1)                                                     00009190
      do 2395 i=1,n1                                                    00009200
      im1=i-1                                                           00009210
c                                                                       00009220
c                                                                       00009230
c                                                                       00009240
c                                                                       00009250
c        nd - box-diagram                                               00009260
c        ----------------                                               00009270
c                                                                       00009280
c                                                                       00009290
c                                                                       00009300
c                                                                       00009310
      if (.not.indnd) go to 2300                                        00009320
c                                                                       00009330
c                                                                       00009340
      ixx1=16*(im1+ixx)                                                 00009350
      iyy1=16*(im1+iyy)                                                 00009360
      id=0                                                              00009370
      do 2245 il=1,imed                                                 00009380
      do 2245 ill=1,imed                                                00009390
      id=id+1                                                           00009400
      if (id.lt.ida) go to 2245                                         00009410
      if (id.gt.ide) go to 2300                                         00009420
c                                                                       00009430
c                                                                       00009440
      do 2205 iv=1,6                                                    00009450
 2205 cvv(iv)=(0.d0,0.d0)                                               00009460
c                                                                       00009470
c                                                                       00009480
c                                                                       00009510
c                                                                       00009520
      do 2211 ii = 1,16                                                 00009530
      vil(ii)=vvnd(ii+ixx1,il)                                          00009540
 2211 vill(ii)=vvnd(ii+iyy1,ill)                                        00009550
c                                                                       00009560
c                                                                       00009570
c                                                                       00009690
c                                                                       00009700
c        sum over discrete intermediate states                          00009710
c                                                                       00009730
      do 2215 ii=1,4                                                    00009740
      ii1=4+ii                                                          00009750
      ii2=8+ii                                                          00009760
      ii3=12+ii                                                         00009770
      ii4=9-ii                                                          00009780
      ii5=17-ii                                                         00009790
      cvv(1)=cvv(1)+(vil(ii)+vil(ii4))*(vill(ii)+vill(ii4))             00009800
      cvv(2)=cvv(2)+(vil(ii2)+vil(ii5))*(vill(ii2)+vill(ii5))           00010080
      cvv(3)=cvv(3)+(vil(ii)-vil(ii4))*(vill(ii)-vill(ii4))             00009810
      cvv(4)=cvv(4)+(vil(ii2)-vil(ii5))*(vill(ii2)-vill(ii5))           00009820
      cvv(5)=cvv(5)+vil(ii)*vill(ii2)                                   00009830
      cvv(5)=cvv(5)+vil(ii1)*vill(ii3)                                  00009840
      cvv(6)=cvv(6)+vil(ii2)*vill(ii)                                   00009850
 2215 cvv(6)=cvv(6)+vil(ii3)*vill(ii1)                                  00009860
c                                                                       00010000
c                                                                       00010010
c                                                                       00010090
c                                                                       00010100
c        one integration step over intermediate momentum                00010110
c                                                                       00010120
c                                                                       00010130
      do 2235 iv=1,6                                                    00010140
 2235 cvj(iv,id)=cvj(iv,id)-cvv(iv)*u(i)                                00010150
c                                                                       00010160
c                                                                       00010170
 2245 continue                                                          00010180
c                                                                       00010190
c                                                                       00010200
c                                                                       00010210
c                                                                       00010220
c        dd - box-diagram                                               00010230
c        ----------------                                               00010240
c                                                                       00010250
c                                                                       00010260
c                                                                       00010270
c                                                                       00010280
 2300 if (.not.inddd) go to 2395                                        00010290
c                                                                       00010300
c                                                                       00010310
      ixx1=32*(im1+ixx)                                                 00010320
      iyy1=32*(im1+iyy)                                                 00010330
      id=ide                                                            00010340
      do 2345 il=1,imedd                                                00010350
      do 2345 ill=1,imedd                                               00010360
      id=id+1                                                           00010370
      if (id.lt.idaa) go to 2345                                        00010380
      if (id.gt.idee) go to 2395                                        00010390
c                                                                       00010400
c                                                                       00010410
      do 2315 ii=1,32                                                   00010420
      vil(ii)=vvdd(ii+ixx1,il)                                          00010430
 2315 vill(ii)=vvdd(ii+iyy1,ill)                                        00010440
c                                                                       00010450
c                                                                       00010460
      do 2321 iv=1,6                                                    00010470
 2321 cvv(iv)=(0.d0,0.d0)                                               00010480
c                                                                       00010490
      do 2325 ii=1,8                                                    00010500
      ii1=8+ii                                                          00010510
      ii2=16+ii                                                         00010520
      ii3=24+ii                                                         00010530
      ii4=17-ii                                                         00010540
      ii5=33-ii                                                         00010550
      cvv(1)=cvv(1)+(vil(ii)-vil(ii4))*(vill(ii)-vill(ii4))             00010560
      cvv(3)=cvv(3)+(vil(ii)+vil(ii4))*(vill(ii)+vill(ii4))             00010570
      cvv(2)=cvv(2)+(vil(ii2)-vil(ii5))*(vill(ii2)-vill(ii5))           00010580
      cvv(4)=cvv(4)+(vil(ii2)+vil(ii5))*(vill(ii2)+vill(ii5))           00010590
      cvv(5)=cvv(5)+vil(ii)*vill(ii2)                                   00010600
      cvv(5)=cvv(5)+vil(ii1)*vill(ii3)                                  00010610
      cvv(6)=cvv(6)+vil(ii2)*vill(ii)                                   00010620
 2325 cvv(6)=cvv(6)+vil(ii3)*vill(ii1)                                  00010630
c                                                                       00010640
c                                                                       00010650
c        integrate                                                      00010660
c                                                                       00010670
      do 2335 iv=1,6                                                    00010680
 2335 cvj(iv,id)=cvj(iv,id)-cvv(iv)*uu(i)                               00010690
c                                                                       00010700
c                                                                       00010710
 2345 continue                                                          00010720
c                                                                       00010730
c                                                                       00010740
c                                                                       00010750
c                                                                       00010760
 2395 continue                                                          00010770
c        this has been the end of the loop of integration               00010780
c                                                                       00010790
c                                                                       00010800
c                                                                       00010810
c                                                                       00010820
c                                                                       00010830
c                                                                       00010840
      int=2                                                             00010850
      do 3055 id=1,idee                                                 00010860
c                                                                       00010870
c        special factors for v5 and v6                                  00010880
c                                                                       00010890
      cvj(5,id)=2.d0*cvj(5,id)                                          00010900
      cvj(6,id)=2.d0*cvj(6,id)                                          00010910
c                                                                       00010920
c                                                                       00010930
      if (j.ne.0) go to 3021                                            00010940
      cvj(2,id)=(0.d0,0.d0)                                             00010950
      cvj(4,id)=(0.d0,0.d0)                                             00010960
      cvj(5,id)=(0.d0,0.d0)                                             00010970
      cvj(6,id)=(0.d0,0.d0)                                             00010980
c                                                                       00010990
c                                                                       00011000
 3021 if (.not.sing) cvj(1,id)=(0.d0,0.d0)                              00011010
      if (.not.trip) cvj(2,id)=(0.d0,0.d0)                              00011020
      if (coup) go to 3030                                              00011030
      do 3025 iv=3,6                                                    00011040
 3025 cvj(iv,id)=(0.d0,0.d0)                                            00011050
c                                                                       00011060
 3030 if (heform) go to 3040                                            00011070
c                                                                       00011080
c                                                                       00011090
c        transformation into lsj-formalism                              00011100
c        ---------------------------------                              00011110
c                                                                       00011120
      cv3=cvj(3,id)                                                     00011130
      cv4=cvj(4,id)                                                     00011140
      cv5=cvj(5,id)                                                     00011150
      cv6=cvj(6,id)                                                     00011160
      cv34=-arjj1*(cv3-cv4)                                             00011170
      cv56=arjj1*(cv5+cv6)                                              00011180
      cvj(3,id)=d2j1*(aj1*cv3+aj*cv4-cv56)                              00011190
      cvj(4,id)=d2j1*(aj*cv3+aj1*cv4+cv56)                              00011200
      cvj(5,id)=d2j1*(cv34-aj1*cv5+aj*cv6)                              00011210
      cvj(6,id)=d2j1*(cv34+aj*cv5-aj1*cv6)                              00011220
c                                                                       00011230
c                                                                       00011240
c        possible different sign depending on the convention used       00011250
      cvj(5,id)=-cvj(5,id)                                              00011260
      cvj(6,id)=-cvj(6,id)                                              00011270
c                                                                       00011280
c                                                                       00011290
c                                                                       00011300
c        multiply with iso-spin factor                                  00011310
c        -----------------------------                                  00011320
c                                                                       00011330
 3040 is=mod(j,2)+1                                                     00011340
      it=mod(is,2)+1                                                    00011350
      if (id.gt.ide) int=3                                              00011360
c                                                                       00011370
c                                                                       00011380
      do 3045 iv=1,6                                                    00011390
      if (iv.eq.2) go to 3043                                           00011400
c                                                                       00011410
c        neglect 3s1-3d1 in case of nd-pi-rho (it.+nonit.) (bagnoud)    00011420
c                                                                       00011430
      if(id.le.ide.and.(id.le.2.or.id.eq.3).and.j.eq.1.and.iv.ge.3)     00011440
     1  go to 3045                                                      00011450
c                                                                       00011460
      cvj(iv,id)=cvj(iv,id)*tt(is,int,itt(id))                          00011470
      go to 3044                                                        00011480
 3043 cvj(iv,id)=cvj(iv,id)*tt(it,int,itt(id))                          00011490
c                                                                       00011500
 3044 cv(iv)=cv(iv)+cvj(iv,id)                                          00011510
 3045 continue                                                          00011520
c                                                                       00011530
 3055 continue                                                          00011540
c                                                                       00011550
c                                                                       00011560
c                                                                       00011570
      return                                                            00011580
      end                                                               00011590
      subroutine delim (itj)                                            00000010
c                                                                       00000020
c                                                                       00000030
c        imaginary part of delta self energy diagram                    00000040
c        -------------------------------------------                    00000050
c           ( q=0 , i.e. delta c.m. system )                            00000060
c                                                                       00000070
c        thresholds:                                                    00000080
c        -----------                                                    00000090
c        itj=3:   nd: 2.*wn+wpi,  dd: wn+wd+wpi                         00000100
c        itj=4:   nd: 2.*wn+wpi,  dd: 2.*wn+2.*wpi                      00000110
c                                                                       00000120
c                                                                       00000130
c                                                                       00000140
      implicit real*8 (a-h,o-z)                                         00000150
c
c
      parameter (inr=96,inr1=inr+1,iimex=41)
      parameter (irsdd=2*inr1*iimex)
c
c                                                                       00000160
      common /crdwrt/ kread,kwrite,kpunch,kda(9)                        00000170
      common /cpoted/ q0qmev,qfmev,pmev,uanspp,wsnspp,ucnspp,udnspp,    00000180
     1                znrl,zrel,smev,noced                              00000190
c                                                                       00000200
c                                                                       00000210
c        common block for all ob-subroutines                            00000220
c                                                                       00000230
      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),  00000240
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,    00000250
     2                ez1,ez2,ct(96),wt(96),                            00000260
     3                ic(10,15),ift(3),mint(3),maxt(3),nt,              00000270
     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),             00000280
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,    00000290
     6                indc(2,15),indpar(3),indxy                        00000300
c                                                                       00000310
c         specifications for this common block                          00000320
c                                                                       00000330
      logical indc,indxy,indpar                                         00000340
c                                                                       00000350
c                                                                       00000360
c         common for self-energy corrections                            00000370
c                                                                       00000380
      common /csthet/ a1(96,2),a2(96,2),a3(96,2),hthi(3,3),xq,xqn,      00000390
     1                eeqa,eeqb,qk,eek,wpi,fs,epz,mgd,mcd,nw,nwt,       00000400
     2                dbcutq,dbcpl,xqvq,
     3                indwt,inddbc,inddbl,indqav
c                                                                       00000410
c         specifications for this common block                          00000420
c                                                                       00000430
      logical indwt,inddbc,inddbl,indqav                                00000440
c                                                                       00000450
c         common block which contains information of tbibdc             00000460
c                                                                       00000470
c        dimension of drsnd and drsdd for n=32 and iime=13              00000490
c                                                                       00000500
c                                                                       00000510
      common /ctbsc/  drsnd(irsdd,2),drsdd(irsdd,3),wimnd(irsdd,2),     00000520
     1                wimdd(irsdd,3),cqx(inr1),qx(inr1),                00000530
     2                zmevm,zzmevm,zpol,zzpol,wn,wd,wpx,nx,issd,        00000540
     3                wdim,wddim,indrnd,indrdd                          00000550
c                                                                       00000560
c         specification for this common block                           00000570
c                                                                       00000580
      logical indrnd,indrdd                                             00000590
      complex*16 drsnd,drsdd,wimnd,wimdd,cqx                            00000600
c                                                                       00000610
c                                                                       00000620
c        further specifications                                         00000630
c                                                                       00000640
c                                                                       00000650
      dimension tdind(2),tdidd(2)                                       00000660
      logical index
c                                                                       00000670
      data index/.false./
      data d23/0.666666666666667d0/                                     00000680
      data xq0/0.d0/                                                    00000690
c
c
c                                                                       00000700
c                                                                       00000710
c        prepare constants                                              00000720
c                                                                       00000730
c
      inter=3
c
c
c                                                                       00000740
      if (index) go to 10
      index=.true.
c
c
c                                                                       00000770
      mgd=mggo(1,inter)                                                 00000780
      med=1                                                             00000790
      mcd=1                                                             00000800
c                                                                       00000810
      wn=wnn(inter)                                                     00000820
      wd=wdd(inter)                                                     00000830
      wnq=wn*wn                                                         00000840
      wdq=wd*wd                                                         00000850
c                                                                       00000860
      dwn=1.d0/wn                                                       00000870
      dwnq=dwn*dwn                                                      00000880
      dwd=1.d0/wd                                                       00000890
      dwdq=dwd*dwd                                                      00000900
      wds=wd*dwn                                                        00000910
c                                                                       00000920
      wpi=wpx                                                           00000930
      wpiq=wpi*wpi                                                      00000940
      dwpiq=1.d0/wpiq                                                   00000950
c                                                                       00000960
      nw=0                                                              00000970
      nwt=-1                                                            00000980
      indwt=.false.                                                     00000990
      indqav=.false.
      inddbc=.false.
      inddbl=.false.
c                                                                       00001000
c                                                                       00001010
c        prepare q dependent expressions                                00001020
c                                                                       00001030
c                                                                       00001040
   10 xq=xq0                                                            00001050
      xqn=xq
      xqq=xq*xq                                                         00001060
c                                                                       00001070
      een=dsqrt(xqq+wnq)                                                00001080
      eed=dsqrt(xqq+wdq)                                                00001090
c                                                                       00001100
      eeqa=eed                                                          00001110
      eeqb=eed                                                          00001120
      deed=1.d0/eed                                                     00001130
      eef1=eed*eed*dwdq*d23                                             00001140
      eef2=eef1*3.d0-0.5d0                                              00001150
      omqq=wpiq+xqq                                                     00001160
      omq=dsqrt(omqq)                                                   00001170
c                                                                       00001180
      thresi=wn+omq                                                     00001190
c                                                                       00001200
      dnd=zrel-een                                                      00001210
c                                                                       00001220
      if (itj.eq.4) go to 1010                                          00001230
      ddd=zrel-eed                                                      00001240
      go to 1040                                                        00001250
c                                                                       00001260
 1010 ddd=zrel-een-omq                                                  00001270
c                                                                       00001280
 1040 do 1050 ii=1,2                                                    00001290
      tdind(ii)=0.d0                                                    00001300
 1050 tdidd(ii)=0.d0                                                    00001310
c                                                                       00001320
      wdim=0.d0                                                         00001330
      wddim=0.d0                                                        00001340
c                                                                       00001350
c                                                                       00001360
      if (dnd.lt.thresi) go to 2000                                     00001370
c                                                                       00001380
c                                                                       00001390
c        delta function for nd                                          00001400
c        ---------------------                                          00001410
c                                                                       00001420
c                                                                       00001430
      dndq=dnd*dnd                                                      00001440
      and=dndq+wnq-omqq                                                 00001450
      qkqd=and*and*0.25d0/dndq-wnq                                      00001460
c                                                                       00001470
      epz=zrel                                                          00001480
      qk=dsqrt(qkqd)                                                    00001490
      eek=dsqrt(qkqd+wnq)                                               00001500
c                                                                       00001510
      qkde=qkqd/eek                                                     00001520
      fs=deed*0.25d0                                                    00001530
c                                                                       00001540
      htx=eeqa*eek*dwnq                                                 00001550
c                                                                       00001560
      hthi(1,inter)=htx+wds                                             00001570
      hthi(2,inter)=eef1                                                00001580
      hthi(3,inter)=eef2                                                00001590
c                                                                       00001600
      call obas (1,med)                                                 00001610
c                                                                       00001620
      do 1200 ii=1,2                                                    00001630
      do 1100 iw=1,nwt                                                  00001640
      tdind(ii)=tdind(ii)+a3(iw,ii)                                     00001650
      a1(iw,ii)=0.d0                                                    00001660
      a2(iw,ii)=0.d0                                                    00001670
 1100 a3(iw,ii)=0.d0                                                    00001680
c                                                                       00001690
 1200 tdind(ii)=tdind(ii)*qkde                                          00001700
c                                                                       00001710
      wdim=-tdind(1)                                                    00001720
c                                                                       00001730
c                                                                       00001740
      if (ddd.lt.thresi) go to 2000                                     00001750
c                                                                       00001760
c                                                                       00001770
c        delta function for dd                                          00001780
c        ---------------------                                          00001790
c                                                                       00001800
c                                                                       00001810
      dddq=ddd*ddd                                                      00001820
      and=dddq+wnq-omqq                                                 00001830
      qkqd=and*and*0.25d0/dddq-wnq                                      00001840
c                                                                       00001850
      qk=dsqrt(qkqd)                                                    00001860
      eek=dsqrt(qkqd+wnq)                                               00001870
c                                                                       00001880
      qkde=2.d0*qkqd/eek                                                00001890
c                                                                       00001900
      htx=eeqa*eek*dwnq                                                 00001910
c                                                                       00001920
      hthi(1,inter)=htx+wds                                             00001930
c                                                                       00001940
      call obas (1,med)                                                 00001950
c                                                                       00001960
      do 1400 ii=1,2                                                    00001970
      do 1300 iw=1,nwt                                                  00001980
      tdidd(ii)=tdidd(ii)+a3(iw,ii)                                     00001990
      a1(iw,ii)=0.d0                                                    00002000
      a2(iw,ii)=0.d0                                                    00002010
 1300 a3(iw,ii)=0.d0                                                    00002020
c                                                                       00002030
 1400 tdidd(ii)=tdidd(ii)*qkde                                          00002040
c                                                                       00002050
      wddim=-tdidd(1)                                                   00002060
c                                                                       00002070
c                                                                       00002080
c                                                                       00002090
 2000 return                                                            00002100
      end                                                               00002110
      subroutine obas (mea,mee)                                         00000010
c                                                                       00000020
c          integration over theta for self-energy diagrams              00000030
c                                                                       00000040
c                                                                       00000050
      implicit real*8 (a-h,o-z)                                         00000060
c                                                                       00000070
      common /crdwrt/ kread,kwrite,kpunch,kda(9)                        00000080
c                                                                       00000090
c        common blocks which contain the arguments and values of the    00000100
c        ob-subroutines                                                 00000110
c                                                                       00000120
      common /cpoted/ q0qmev,qfmev,pmev,uanspp,wsnspp,ucnspp,udnspp,    00000130
     1                znrl,zrel,smev,noced                              00000140
c                                                                       00000150
c        specifications for these common blocks                         00000160
      logical noced                                                     00000170
c                                                                       00000180
c                                                                       00000190
c        common block for all ob-subroutines                            00000200
c                                                                       00000210
      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),  00000220
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,    00000230
     2                ez1,ez2,ct(96),wt(96),                            00000240
     3                ic(10,15),ift(3),mint(3),maxt(3),nt,              00000250
     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),             00000260
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,idee,idde,   00000270
     6                indc(2,15),indpar(3),indxy                        00000280
c                                                                       00000290
c         specifications for this common block                          00000300
c                                                                       00000310
      logical indc,indxy,indpar                                         00000320
c                                                                       00000330
c         common for self-energy corrections                            00000340
c                                                                       00000350
      common /csthet/ a1(96,2),a2(96,2),a3(96,2),hthi(3,3),xq,xqn,      00000360
     1                eeqa,eeqb,qk,eek,wpi,fs,epz,mgd,mcd,nw,nwt,       00000370
     2                dbcutq,dbcpl,xqvq,                                00000380
     3                indwt,inddbc,inddbl,indqav                        00000390
c                                                                       00000400
c                                                                       00000410
c         specification for this common block                           00000420
c                                                                       00000430
      logical indwt,inddbc,inddbl,indqav                                00000440
c                                                                       00000450
c        further specifications                                         00000460
c                                                                       00000470
      dimension aa1(96),aa2(96),aa3(96),deltak(96,3)                    00000480
      dimension cwt(96),swt(96),xnum(3)                                 00000490
      data nnt/-1/,iinter/-1/
      real*4 axy2,aomq,am1,am2,am,aez                                   00000510
      logical indret(3)
      data indret/3*.false./
      logical indrt                                                     00000530
      logical indnn                                                     00000540
      logical indiso                                                    00000550
      data d23/0.666666666666d0/                                        00000560
      data xq0/0.d0/                                                    00000570
c
c
c
c                                                                       00000580
      if (mcd.ne.1) go to 1000                                          00000590
c                                                                       00000600
c                                                                       00000610
      if (inter.eq.iinter) go to 50                                     00000620
      iinter=inter                                                      00000630
c                                                                       00000640
      wn=wnn(inter)                                                     00000650
      dwn=1.d0/wn                                                       00000660
      dwnq=dwn*dwn                                                      00000670
c                                                                       00000680
      indnn=indret(inter)                                               00000690
      indrt=indret(inter)                                               00000700
      if (inter.eq.1) indnn=.true.                                      00000710
      if (xqn.eq.0.d0.and.inter.eq.3) indrt=.true.                      00000720
c                                                                       00000730
      min=mint(inter)                                                   00000740
      max=maxt(inter)                                                   00000750
c                                                                       00000760
   50 x1=xq*dwn                                                         00000770
      y1=qk*dwn                                                         00000780
      x1n=xqn*dwn                                                       00000790
      xqq=x1*x1                                                         00000800
      xqqn=x1n*x1n                                                      00000810
      xkk=y1*y1                                                         00000820
      xqk=x1*y1                                                         00000830
      xqkn=x1n*y1                                                       00000840
      xqk2=xqk*2.d0                                                     00000850
c                                                                       00000860
      xkkh=xkk*0.5d0                                                    00000870
      xkkd=xkkh*d23                                                     00000880
      qqpkk=xqq+xkk                                                     00000890
c                                                                       00000900
      eqmek=(eeqa-eek)*dwn                                              00000910
      zmeqk=(epz-eeqb-eek)*dwn                                          00000920
c                                                                       00000930
      if (indwt) go to 99                                               00000940
c                                                                       00000950
c                                                                       00000960
c         compute am                                                    00000970
c                                                                       00000980
      axy2=xqk2                                                         00000990
      aomq=qqpkk+wpi*wpi*dwnq                                           00001000
      am1=axy2/aomq                                                     00001010
      aez=-zmeqk                                                        00001020
      am2=axy2/(aomq+aez*abs(aez))                                      00001030
      if (am2.lt.0) go to 94                                            00001040
      am=amax1(am1,am2)                                                 00001050
c                                                                       00001060
c         compute number of gausspoints nwt                             00001070
c                                                                       00001080
      if (am.gt.0.999) go to 94                                         00001090
      if (am.gt.0.85) am=am**(-log(1.-am)-0.9)                          00001100
      nwt=float(min)/(1.-am)+0.9                                        00001110
      if (nwt.gt.max) nwt=max                                           00001120
      go to 95                                                          00001130
   94 nwt=max                                                           00001140
c                                                                       00001150
c         compute nwt suitable for gset                                 00001160
c                                                                       00001170
   95 if (nwt.le.16) go to 100                                          00001180
      if (nwt.gt.24) go to 96                                           00001190
      nwt=4*(nwt/4)                                                     00001200
      go to 100                                                         00001210
   96 if (nwt.gt.48) go to 97                                           00001220
      nwt=8*(nwt/8)                                                     00001230
      go to 100                                                         00001240
   97 nwt=16*(nwt/16)                                                   00001250
      if (nwt.gt.96) nwt=96                                             00001260
      go to 100                                                         00001270
c                                                                       00001280
   99 nwt=nw                                                            00001290
c                                                                       00001300
  100 if (nwt.eq.nnt) go to 105                                         00001310
c                                                                       00001320
c         call gauss points and weihts                                  00001330
c                                                                       00001340
      call gset(-1.d0,1.d0,nwt,cwt,swt)                                 00001350
c                                                                       00001360
      nnt=nwt                                                           00001370
  105 continue                                                          00001380
c                                                                       00001390
      do 110 i=1,nwt                                                    00001400
      xk2t=xqk2*cwt(i)                                                  00001410
      deltak(i,1)= xk2t-qqpkk                                           00001420
      deltak(i,2)= -xkk                                                 00001430
  110 deltak(i,3)=-xk2t-qqpkk                                           00001440
c                                                                       00001450
c                                                                       00001460
c                                                                       00001470
 1000 do 5000 ms=mea,mee                                                00001480
      ims=ima(ms,mgd,inter)                                             00001490
c                                                                       00001500
      if (mcd.ne.1) go to 2000                                          00001510
c                                                                       00001520
c                                                                       00001530
c         compute propagators                                           00001540
c                                                                       00001550
c                                                                       00001560
      c4=c(4,ims)                                                       00001570
      iprsp=ic(1,ims)                                                   00001580
      iret=iprsp+1                                                      00001590
c                                                                       00001600
      if (inter.ne.1) go to 1200                                        00001610
c                                                                       00001620
c                                                                       00001630
c        nucleon                                                        00001640
c        -------                                                        00001650
c                                                                       00001660
c                                                                       00001670
      do 1110 i=1,nwt                                                   00001680
      omvq=c4-deltak(i,2)+xqvq                                          00001690
      omv=dsqrt(omvq)                                                   00001700
      omq=c4-deltak(i,1)                                                00001710
      om=dsqrt(omq)                                                     00001720
      omep=eqmek-om                                                     00001730
      omepq=omep*omep                                                   00001740
      if (indqav) go to 1109                                            00001750
      omv=om                                                            00001760
 1109 omez=(zmeqk-omv)*omepq*om                                         00001770
      aa3(i)=swt(i)/(om*omep)                                           00001780
      aa2(i)=swt(i)/(om*omepq)                                          00001790
 1110 aa1(i)=swt(i)/omez                                                00001800
      go to 1250                                                        00001810
c                                                                       00001820
c                                                                       00001830
c         delta                                                         00001840
c         -----                                                         00001850
c                                                                       00001860
c                                                                       00001870
 1200 iret=1                                                            00001880
      if (indrt) iret=2                                                 00001890
      do 1210 i=1,nwt                                                   00001900
      omq=c4-deltak(i,1)                                                00001910
      om=dsqrt(omq)                                                     00001920
      omkq=c4-deltak(i,iret)                                            00001930
      omk=dsqrt(omkq)                                                   00001940
      omvq=omkq+xqvq                                                    00001950
      omv=dsqrt(omvq)                                                   00001960
      if (indqav) go to 1209                                            00001970
      omv=om                                                            00001980
 1209 omez=zmeqk-omv                                                    00001990
      som=swt(i)/omk                                                    00002000
      aa3(i)=som                                                        00002010
      aa2(i)=som                                                        00002020
 1210 aa1(i)=som/omez                                                   00002030
c                                                                       00002040
c                                                                       00002050
c                                                                       00002060
c                                                                       00002070
c         cut-off                                                       00002080
c         -------                                                       00002090
c                                                                       00002100
c                                                                       00002110
 1250 mi=4                                                              00002120
      mm=5                                                              00002130
c                                                                       00002140
 1300 ityp=ic(mi,ims)                                                   00002150
      if (ityp.eq.0) go to 1999                                         00002160
      iprspc=ic(mi+1,ims)                                               00002170
      iret=iprspc+1                                                     00002180
      if (indrt) iret=2                                                 00002190
      if (iprspc.eq.3) iret=1                                           00002200
      go to (1400,1400,1999,1999,1999,1999,1700,1800,1999), ityp        00002210
c                                                                       00002220
c                                                                       00002230
c                                                                       00002240
c         cutoff of dipole type                                         00002250
c         *********************                                         00002260
c                                                                       00002270
c                                                                       00002280
c                                                                       00002290
 1400 c5=c(mm,ims)                                                      00002300
      c6=c(mm+1,ims)                                                    00002310
      nexp=ic(mi+2,ims)                                                 00002320
c                                                                       00002330
      if (indnn) go to 1405                                             00002340
      if (.not.inddbc) go to 1405                                       00002350
      c6=dbcutq                                                         00002360
      c5=c6-c4                                                          00002370
c                                                                       00002380
 1405 do 1410 i=1,nwt                                                   00002390
c                                                                       00002400
      aaa=c5/(c6-deltak(i,iret))                                        00002410
c     --------------------------                                        00002420
c                                                                       00002430
      do 1410 ii=1,nexp                                                 00002440
      aa1(i)=aa1(i)*aaa                                                 00002450
      aa2(i)=aa2(i)*aaa                                                 00002460
 1410 aa3(i)=aa3(i)*aaa                                                 00002470
c                                                                       00002480
      if (iprspc.le.2) go to 1430                                       00002490
c                                                                       00002500
      do 1420 i=1,nwt                                                   00002510
c                                                                       00002520
      aaa=c5/(c6-deltak(i,iprspc))                                      00002530
c                                                                       00002540
      do 1420 ii=1,nexp                                                 00002550
      aa1(i)=aa1(i)*aaa                                                 00002560
      aa2(i)=aa2(i)*aaa                                                 00002570
 1420 aa3(i)=aa3(i)*aaa                                                 00002580
c                                                                       00002590
c                                                                       00002600
 1430 mi=mi+3                                                           00002610
      mm=mm+2                                                           00002620
      go to 1300                                                        00002630
c                                                                       00002640
c                                                                       00002650
c                                                                       00002660
c        exponential form factor                                        00002670
c        ***********************                                        00002680
c                                                                       00002690
c                                                                       00002700
c                                                                       00002710
 1700 c5=c(mm,ims)                                                      00002720
      do 1710 i=1,nwt                                                   00002730
c                                                                       00002740
      expo=deltak(i,iret)                                               00002750
c     -------------------                                               00002760
c                                                                       00002770
      if (expo.lt.-50.d0) go to 1704                                    00002780
c                                                                       00002790
      aaa=dexp(expo)                                                    00002800
c     --------------                                                    00002810
c                                                                       00002820
      aa1(i)=aa1(i)*aaa                                                 00002830
      aa2(i)=aa2(i)*aaa                                                 00002840
      aa3(i)=aa3(i)*aaa                                                 00002850
      go to 1710                                                        00002860
c                                                                       00002870
 1704 aa1(i)=0.d0                                                       00002880
      aa2(i)=0.d0                                                       00002890
      aa3(i)=0.d0                                                       00002900
c                                                                       00002910
c                                                                       00002920
 1710 continue                                                          00002930
      mi=mi+2                                                           00002940
      mm=mm+1                                                           00002950
      go to 1300                                                        00002960
c                                                                       00002970
c                                                                       00002980
c                                                                       00002990
c        cloudy bag form factor                                         00003000
c        **********************                                         00003010
c                                                                       00003020
c                                                                       00003030
c                                                                       00003040
 1800 c5=c(mm,ims)                                                      00003050
      nexp=ic(mi+2,ims)                                                 00003060
c                                                                       00003070
      do 1805 i=1,nwt                                                   00003080
      arg=dsqrt(-deltak(i,iret))*c5                                     00003090
      argc=arg*arg*arg                                                  00003100
c                                                                       00003110
      aaa=3.d0*(dsin(arg)-arg*dcos(arg))/argc                           00003120
c     ---------------------------------------                           00003130
c                                                                       00003140
      do 1805 ii=1,nexp                                                 00003150
      aa1(i)=aa1(i)*aaa                                                 00003160
      aa2(i)=aa2(i)*aaa                                                 00003170
 1805 aa3(i)=aa3(i)*aaa                                                 00003180
c                                                                       00003190
      mi=mi+3                                                           00003200
      mm=mm+1                                                           00003210
      go to 1300                                                        00003220
c                                                                       00003230
c                                                                       00003240
 1999 continue                                                          00003250
c                                                                       00003260
c                                                                       00003270
c                                                                       00003280
c         numerators and factors                                        00003290
c         ----------------------                                        00003300
c                                                                       00003310
c                                                                       00003320
c                                                                       00003330
 2000 go to (3000,4000,4000), inter                                     00003340
c                                                                       00003350
c                                                                       00003360
c                                                                       00003370
 3000 cmc=c(mcd,ims)                                                    00003380
      indiso=indc(1,ims)                                                00003390
      fc=fs*cmc                                                         00003400
      if (indiso) fc=fc*3.d0                                            00003410
      hthi1=hthi(1,inter)                                               00003420
      hthi2=hthi(2,inter)                                               00003430
c                                                                       00003440
c                                                                       00003450
      do 3500 i=1,nwt                                                   00003460
c                                                                       00003470
      xkt=xqk*cwt(i)                                                    00003480
c                                                                       00003490
      go to (3110,3200,3200,3140,3200,3160,3200,3200,3200,3200,3200,    00003500
     1       3200),mgd                                                  00003510
c                                                                       00003520
 3110 xnum(1)=hthi1-xkt                                                 00003530
      go to 3200                                                        00003540
c                                                                       00003550
 3140 xnum(1)=hthi1-xkt                                                 00003560
      go to 3200                                                        00003570
c                                                                       00003580
 3160 xnum(1)=hthi1-xkt                                                 00003590
      xnn=hthi2-xkt                                                     00003600
      xnum(2)=hthi1-4.d0*xkt-0.5d0*xnn*xnn                              00003610
      xnum(3)=hthi1-xkt                                                 00003620
c                                                                       00003630
 3200 continue                                                          00003640
c                                                                       00003650
c                                                                       00003660
      a1(i,1)=a1(i,1)+aa1(i)*xnum(mcd)*fc*dwnq                          00003670
      a2(i,1)=a2(i,1)+aa2(i)*xnum(mcd)*fc*dwn                           00003680
 3500 a3(i,1)=a3(i,1)+aa3(i)*xnum(mcd)*fc                               00003690
c                                                                       00003700
      go to 5000                                                        00003710
c                                                                       00003720
c                                                                       00003730
c                                                                       00003740
c                                                                       00003750
 4000 cmc=c(mcd,ims)                                                    00003760
      if (inddbl) cmc=dbcpl                                             00003770
      fc=fs*cmc/c4                                                      00003780
c                                                                       00003790
      hthi1=hthi(1,inter)                                               00003800
      hthi2=hthi(2,inter)                                               00003810
      hthi3=hthi(3,inter)                                               00003820
c                                                                       00003830
      do 4200 i=1,nwt                                                   00003840
c                                                                       00003850
      xkt=xqkn*cwt(i)                                                   00003860
      xkt2=xkt*2.d0                                                     00003870
      cwq=cwt(i)*cwt(i)                                                 00003880
      xnn=hthi1-xkt                                                     00003890
      xnum(1)=(1.d0-cwq)*xkkh                                           00003900
      xnum(2)=(0.5d0+hthi3*cwq)*xkkd+(xqqn-xkt2)*hthi2                  00003910
c                                                                       00003920
      ab1=aa1(i)*xnn*fc                                                 00003930
      ab2=aa2(i)*xnn*fc*wn                                              00003940
      ab3=aa3(i)*xnn*fc*wn                                              00003950
c                                                                       00003960
      do 4200 in=1,2                                                    00003970
c                                                                       00003980
      a1(i,in)=a1(i,in)+ab1*xnum(in)                                    00003990
      a2(i,in)=a2(i,in)+ab2*xnum(in)                                    00004000
 4200 a3(i,in)=a3(i,in)+ab3*xnum(in)                                    00004010
c                                                                       00004020
c                                                                       00004030
c                                                                       00004040
 5000 continue                                                          00004050
c                                                                       00004060
      return                                                            00004070
      end                                                               00004080
      function cspe (cqq,wn,wnim,iprop,ispex)                           00000010
c                                                                       00000020
c                                                                       00000030
c         single particle energy with complex particle mass             00000040
c                                                                       00000050
c         cspe is called by:  tbibdc                                    00000060
c                             tbibsc                                    00000070
c                                                                       00000080
c                                                                       00000090
      implicit real*8 (a,b,d-h,o-z), complex*16 (c)                     00000100
c                                                                       00000110
      data wwn/938.926d0/                                               00000120
c                                                                       00000130
c                                                                       00000140
c**** statement functions ************
      cdsqrt(cc)=zsqrt(cc)
c
c
c
c                                                                       00000150
      ispe=ispex                                                        00000160
      cwn=dcmplx(wn,wnim)                                               00000170
c                                                                       00000180
      go to (1000,2000), iprop                                          00000190
c                                                                       00000200
c                                                                       00000210
 1000 cwmn=cwn-dcmplx(wwn,0.d0)                                         00000220
      cspe=0.5d0*cqq/cwn+cwmn                                           00000230
      go to 9000                                                        00000240
c                                                                       00000250
c                                                                       00000260
c                                                                       00000270
 2000 cwnq=cwn*cwn                                                      00000280
      cspe=cdsqrt(cqq+cwnq)                                             00000290
c                                                                       00000300
c                                                                       00000310
c                                                                       00000320
 9000 return                                                            00000330
      end                                                               00000340
