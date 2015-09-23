
      subroutine cmplx_cdbonn
c     
c     alias  "CDBonn2000"
c     
c******************************************************************
c     
c     FINAL VERSION: JANUARY 1, 2000
c     (FORTRAN improved by W. Schadow 7/23/00)
c     
c     this potential is published in:
c     R. Machleidt, PRC 63, 024001 (2001).
c     
c******************************************************************
c     
c     This code computes the
c     Charge-Dependent Bonn NN-Potential (`CDBonn2000'),
c     ------------------------------------------------
c     in momentum space.
c     
c     This CD-Bonn potential includes CSB and CIB effects derived
c     from irreducible 2pi and pi-rho exchanges with nucleon
c     and delta(1232) intermediate states. CIB also includes
c     the effects of irreducible pi-gamma exchange.
c     Besides this, the usual CIB due to OPE is included.
c     Note that the latter is contained in any modern high-precision
c     potential, while the former is not.
c     
c*******************************************************************
c     
c     
c     this package is self-contained and includes
c     all subroutines needed.
c     only cdbonn needs to be called by the user.
c     all codes are consistently in double precision.
c     when working on an UNIX system, it is crucial
c     to compile this code with the  -static  option.
c     
c     NOTE: as compared to the earlier version of cdbonn, there is
c     a minor practical change in the code: this code does not
c     read-in anything. the type of potential to be used
c     (pp, np, or nn) is now controlled by a new common block,
c     
c     common /cnn/ inn  ,
c     where
c     inn=1  means pp potential,
c     inn=2  means np potential, and
c     inn=3  means nn potential.
c     
c     the user needs to include this common block in his/her code,
c     and specify which potential he/she wants. inn can be
c     changed at any time, i.e., all three potentials can be
c     used in one run, which is the advantage of this new regulation.
c     
c     
c************************************************************************
c     
c     
c     author:      Ruprecht Machleidt
c     department of physics
c     university of idaho
c     moscow, idaho 83844
c     u. s. a.
c     e-mail: machleid@uidaho.edu
c     
c     formerly:
c     institut fuer theoretische kernphysik der
c     universitaet bonn
c     nussallee 14-16
c     d - 5300  bonn, w. germany
c     
c*********************************************************************
c     
c     this version of the code uses the Legendre functions
c     ----------------------------------------------------
c     of the second kind for the partial wave decomposition
c     -----------------------------------------------------
c     and includes the meson-parameters in data statements.
c     -----------------------------------------------------
c     
c     
c     an earlier version of this code has been published in
c     "computational nuclear physics 2 -- nuclear reactions",
c     K. Langanke, J.A. Maruhn, and S.E. Koonin, eds.
c     (Springer, New York, 1993), Chapter 1, pp. 1-29.
c     
c     This code is a slight modification of the earlier, published
c     version. However, the mathematical formulae, as well as the
c     general organization of this code is the same as described in
c     the above-referenced book-chapter.
c     in this version of the code, the integrals, Eqs. (1.68) of the
c     above reference, are solved analytically by means of the
c     Legendre functions of the second kind, see Eqs. (E.44) of
c     R. Machleidt et al., Phys. Rep. 149, 1 (1987).
c     
c     Still, the above-referenced article may serve as a good
c     introduction into this code.
c     
c*********************************************************************
c     
c     
      implicit complex*16 (a-h,o-z)
      real*8 dwn, wnn, fac , aj, aj1, a2j1, aaj6, deter
      real*8 vl,adminv,ldminv,mdminv
c     
c     
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c     
c     arguments and values of this subroutine:
c     
      common /cpot/   v(6),xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep,label
      common /cnn/ inn
c     
c     
c     this has been the end of the common-blocks containing
c     the arguments and values of this subroutine
c     
c     specifications for these two common blocks
c     
      logical heform,sing,trip,coup,endep
c     
c     THE ABOVE FOUR COMMON BLOCKS IS ALL THE USER NEEDS
c     TO BE FAMILIAR WITH.
c     
c     
c     xmev and ymev are the final and initial relative momenta,
c     respectively, in units of mev/c.
c     v is the potential in units of mev**(-2).
c     concerning units and factor of pi etc.,
c     cf. with the partial-wave Lippmann-Schwinger equation, Eq. (1.32),
c     and with the phase shift relation, Eq. (1.41) of
c     R. Machleidt, in: Computational Nuclear Physics 2
c     -- Nuclear Reactions, Langanke et al., eds.
c     (Springer, New York, 1993), Chapter 1, pp. 1-29.
c     
c     the partial-wave Lippmann-Schwinger equation for the
c     K-matrix reads:
c     
c     K(q',q) = V(q',q) + M P \int dk k^2 V(q',k) K(k,q)/(q^2-k^2)
c     
c     with M the nucleon mass in MeV and P denoting the principal value;
c     V(q',q) as provided by this code in common block /cpot/;
c     all momenta in MeV.
c     
c     the phase-shift relation is:
c     
c     tan \delta_L = -(pi/2) M q K_L(q,q)
c     
c     with M and q in units of MeV, K_L in MeV**(-2) like V.
c     
c     
c     if heform=.true., v contains the 6 matrix elements
c     associated with one j in the helicity formalism
c     in the following order:
c     0v, 1v, 12v, 34v, 55v, 66v (for notation see above article).
c     if heform=.false., v contains the 6 matrix elements
c     associated with one j in the lsj formalism
c     in the following order:
c     0v(singlet), 1v(uncoupled triplet), v++, v--, v+-, v-+ (coupled)
c     (see above article for notation).
c     j is the total angular momentum. there is essentially no upper
c     limit for j.
c     sing, trip, and coup should in general be .true..
c     endep and label can be ignored.
c     it is customary, to set kread=5 and kwrite=6;
c     ignore kpunch and kda(9).
c     inn=1  means pp potential,
c     inn=2  means np potential, and
c     inn=3  means nn potential.
c     
c     
c     THIS IS ESSENTIALLABS(Y).ALL THE USER NEEDS TO KNOW.
c     
c**********************************************************************
c     
c     
c     common block for all ob-subroutines
c     
      common /cobq/   vj(32,50),c(20,50),fff,ff,f(52),aa(96),ai(19,15),
     1     wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2     ez1,ez2,ct(96),wt(96),
     3     ic(20,50),ift(3),mint(3),maxt(3),nt,
     4     mge,mgg(12,3),mggo(12,3),ima(15,12,3),
     5     imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6     indc(10,50),indpar(3),indxy
c     
c     specifications for this common block
c     
      logical indc,indxy,indpar
c     
c     
c     further specifications
c     


      character*4 mesong(12)
      logical index

      dimension vl(4),adminv(4,4),ldminv(4),mdminv(4)


      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     1     '1-  ','1-t ','1-tt','1-st','1-ss',
     2     '1+  ','2+  '/
      data pi/3.141592653589793d0/
      data index/.false./
      data jj/-1/
      data innn/-1/
      save
c     
c     
c     
c     
      inter=1


c     
c     
      if (inn.lt.1.or.inn.gt.3) then
         write (kwrite,19001) inn
19001    format (////' error in cdbonn: potential type  inn =',i3,
     1        '  unknown.'/' execution terminated.'////)
         stop
      endif
      if (j.lt.0) then
         write (kwrite,19002)
19002    format (////' error in cdbonn: total angular momentum j',
     1        ' is negative.'/' execution terminated.'////)
         stop
      endif
c     
c     
c     
c     
c     call obparq whenever j or inn has changed
c     
c     
      if (j.eq.jj.and.inn.eq.innn) go to 50
      innn=inn
c     
c     
c     
      call obparq
c     -----------
c     -----------
c     
c     
      dwn=dble(1.d0/wnn(inter))
c     
c     prepare constant over-all factor
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
 30   if (j.eq.jj) go to 50
      jj=j
      if (j.eq.0) go to 50
      aj=dble(j)
      aj1=dble(j+1)
      a2j1=dble(2*j+1)
      aaj6=sqrt(aj*aj1)
c     
c     coefficient matrix for the translations into lsj formalism
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
c     inversion
c     
      call dminv (adminv,4,deter,ldminv,mdminv)
!      call  cmplxmatinv(adminv,4)
c     
c     
c     
c     
c     
c     
c     
c     prepare expressions depending on x and y
c     ----------------------------------------
c     ----------------------------------------
c     
c     
c     
c     
 50   x=xmev*dwn                   
      y=ymev*dwn
      indxy=.false.
c     
c     
c     
      xx=x*x
      yy=y*y
      xy2=x*y*2.d0
      xxpyy=xx+yy
      ex=sqrt(1.d0+xx)
      ey=sqrt(1.d0+yy)
      eem12=(ex*ey-1.d0)*2.d0
c     
c     
c     
c     
      xy=xy2*0.5d0
      ee=ex*ey
      ree=sqrt(ee)
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
c     prepare over-all factor
c     
c     
      go to (70,71,72),iftgo
c     
c     no additional factor
c     
 70   fff=fac
      go to 90
c     
c     minimal relativity
c     
 71   fff=fac/ree
      go to 90
c     
c     factor m/e*m/e
c     
 72   fff=fac/ee
c     
c     
c     
c     
c     
c     
 90   do 93 iv=1,6
 93      v(iv)=0.d0
         do 95 il=1,50
            do 95 iv=1,6
 95            vj(iv,il)=0.d0
c     
c     
c     
c     
c     contributions of mesons
c     -----------------------
c     -----------------------
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
c     0-  , pseudo-scalar coupling
c     ----------------------------
c     
c     
c     
c     
 100              mc=1
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
                  call obstrq(1,1,me)
                  go to 1995
c     
c     
c     
c     
c     0+  , scalar coupling
c     ---------------------
c     
c     
c     
c     
 400              mc=1
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
                  call obstrq(1,1,me)
                  go to 1995
c     
c     
c     
c     
c     1-t , vector mesons
c     -------------------
c     
c     
c     
c     
c     vector-vector coupling
c     
c     
c     
c     
 700              mc=1
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
                  call obstrq(1,1,me)
c     
c     
c     
c     
c     tensor-tensor coupling
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
c     factors for additional terms
                  f(9)=-2.d0*xxyy
                  f(10)=eep1*xy2
                  f(11)=-epe*xy2
c     
                  call obstrq(2,1,me)
c     
c     
c     
c     
c     vector-tensor coupling
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
                  call obstrq(1,1,me)
                  go to 1995
c     
c     
c     
c     
c     this has been the end of the contributions of mesons
c     ----------------------------------------------------
c     
c     
c     
c     
c     error exit
c     ----------
c     
c     
c     
c     
 9000             write (kwrite,19000) mesong(mg)
19000             format(////' error in cdbonn:   meson-group   ',a4,'  does not exi
     1                 st in this program.'/' execution terminated.'
     2                 ////)
                  stop
c     
c     
c     
c     
 1995          continue
c     
c     
c     
c     
c     add up contributions of mesons
c     ------------------------------
c     
c     
c     
c     
 2000          continue

               do 2005 iv=1,6
 2005             v(iv)=vj(iv,iman1)+vj(iv,iman3)
c     
c     
                  if (j.eq.1) then
                     v(1)=vj(1,iman1)+vj(1,iman4)
                  end if

c     
c     
                  if (j.eq.2) then
                     do 2007 iv=3,6
 2007                   v(iv)=vj(iv,iman2)+vj(iv,iman3)
                     end if
c     
c     
                     if (mod(j,2).eq.1) go to 2020
c     
c     
c     j even
c     ------
c     
                     v(1)=v(1)+vj(1,iman5)+vj(1,iman6)
                     v(2)=v(2)+vj(2,iman9)+vj(2,iman10)
c     
                     if (j.eq.2) then
c     
c     the pions for 3P2-3F2
                        do 2014 iv=3,6
 2014                      v(iv)=v(iv)+vj(iv,iman7)+vj(iv,iman8)
                        else
c     
c     the pions in all other T=1 coupled cases
                           do 2015 iv=3,6
 2015                         v(iv)=v(iv)+vj(iv,iman5)+vj(iv,iman6)
                           end if
                           go to 2030
c     
c     
c     j odd
c     -----
c     
 2020                      v(1)=v(1)+vj(1,iman9)+vj(1,iman10)
                           v(2)=v(2)+vj(2,iman5)+vj(2,iman6)
c     
                           do 2025 iv=3,6
 2025                         v(iv)=v(iv)+vj(iv,iman9)+vj(iv,iman10)
c     
c     
c     for all j
c     _________
c     
 2030                         v(1)=v(1)+vj(1,iman11)+vj(1,iman12)
                              v(2)=v(2)+vj(2,iman13)+vj(2,iman14)
                              v(3)=v(3)+vj(3,iman15)+vj(3,iman16)
                              v(4)=v(4)+vj(4,iman17)+vj(4,iman18)

                              do 2035 iv=5,6
 2035                            v(iv)=v(iv)+vj(iv,iman17)+vj(iv,iman18)
c     
c     
c     
c     
                                 if (j.eq.0.or..not.heform) go to 4000
c     
c     
c     translation into (combinations of) helicity states
c     
c     
                                 do 3005 i=1,4
 3005                               vl(i)=v(i+2)
c     
                                    do 3020 ii=1,4
                                       iii=ii+2
                                       v(iii)=0.d0
c     
                                       do 3015 i=1,4
 3015               v(iii)=v(iii)+adminv(ii,i)*vl(i)
 3020              v(iii)=v(iii)*a2j1
c     
c     
c     
c     
 4000                                     return
                                          end
      subroutine obparq
c
c        obparq provides the parameters for the
c        charge-dependent Bonn potential.
c
c
      implicit complex*16 (a-h,o-z)
      real*8 dwn, wnn, fac

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
      common /cobq/   vj(32,50),c(20,50),fff,ff,f(52),aa(96),ai(19,15),
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
      character*4 name(4),nname(15)
      dimension wscale(3)
      integer imga(3)

      character*4 nucnuc(3)
c      integer nucnuc(3)

      character*4 cut,end
c      integer cut,end

      character*4 two
c      integer two

      logical index

      character*4 mesong(12)
c      integer mesong(12)

      character*4 ntab1(4,10)
c      integer ntab1(4,10)

      character*4 ntab2(4)
c      integer ntab2(4)

      dimension tab1(5,10,3)
      dimension tab2(5,4,7,3)

      data cut/'cut '/,end/'end '/
      data two/'2   '/
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     1            '1-  ','1-t ','1-tt','1-st','1-ss',
     2            '1+  ','2+  '/
      data nucnuc/'BNpp','BNnp','BNnn'/
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

      data ntab2/
     1 '0+  ','2','sig','mas '/
c
c
c        global parameters
c        -----------------
c        -----------------

      data tab1/
c
c proton-proton potential
c -----------------------
     1        0.84    , 6.1       , 769.9     , 1.        , 0.,
     2        2.      , 0.        , 2.        , 1310.     , 0.1,
     3        0.84    , 6.1       , 769.9     , 1.        , 0.,
     4        2.      , 0.        , 2.        , 1310.     , 0.,
     5       20.0     , 0.        , 781.94    , 0.        , 0.,
     6        2.      , 0.        , 2.        , 1500.     , 0.,
     7       20.0     , 0.        , 781.94    , 0.        , 0.,
c t=1:
     8       13.6     , 134.9764  , 0.0       , 139.56995 , 1720.,
     9       13.6     , 134.9764  , 0.0       , 139.56995 , 3000.,
c t=0:
     *        0.0     , 134.9764  , 0.0       , 139.56995 , 1720.,
c
c neutron-proton potential
c ------------------------
     1        0.84    , 6.1       , 769.9     , 1.        , 0.,
     2        2.      , 0.        , 2.        , 1310.     , 0.1,
     3        0.862   , 6.1       , 769.9     , 1.        , 0.,
     4        2.      , 0.        , 2.        , 1310.     , 0.,
     5       20.0     , 0.        , 781.94    , 0.        , 0.,
     6        2.      , 0.        , 2.        , 1500.     , 0.,
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
     3        0.844   , 6.1       , 769.9     , 1.        , 0.,
     4        2.      , 0.        , 2.        , 1310.     , 0.,
     5       20.0     , 0.        , 781.94    , 0.        , 0.,
     6        2.      , 0.        , 2.        , 1500.     , 0.,
     7       20.0     , 0.        , 781.94    , 0.        , 0.,
c t=1:
     8       13.6     , 134.9764  , 0.0       , 139.56995 , 1720.,
     9       13.6     , 134.9764  , 0.0       , 139.56995 , 3000.,
c t=0:
     *        0.0     , 134.9764  , 0.0       , 139.56995 , 1720./
c
c
c         partial-wave dependent parameters
c         ---------------------------------
c         ---------------------------------

      data tab2/
c
c proton-proton potential
c -----------------------
c j=0:
     1        4.24591,  452.,       17.61,      1225.,      2500.,
     2        0.0    ,  350.,       0.0  ,       793.,      2500.,
     3        7.866  ,  560.,       17.61,      1225.,      2500.,
     4        0.0    ,  452.,       0.0  ,      1225.,      2500.,
c j=1:
     1        0.0    ,  350.,       0.0  ,      1225.,      2500.,
     2        2.303  ,  424.,       17.61,      1225.,      2500.,
     3        0.0    ,  350.,       0.0  ,       793.,      2500.,
     4        0.0    ,  350.,       0.0  ,       793.,      2500.,
c j=2:
     1        2.225  ,  400.,       190.7,      1225.,      2500.,
     2        0.0    ,  350.,       0.0  ,      1225.,      2500.,
     3        1.5    ,  452.,       56.21,       793.,      2500.,
     4        4.166  ,  470.,       24.80,      1225.,      2500.,
c j=3:
     1        0.0    ,  350.,       0.0  ,       793.,      2500.,
     2        1.5    ,  452.,       74.44,       793.,      2500.,
     3        0 0    ,  350.,       0.0  ,       793.,      2500.,
     4        0.0    ,  452.,       0.0  ,       793.,      2500.,
c j=4:
     1        4.24591,  452.,       0.0  ,      1225.,      2500.,
     2        0.0    ,  350.,       0.0  ,       793.,      2500.,
     3        3.8    ,  452.,       17.61,      1225.,      2500.,
     4        3.8    ,  452.,       17.61,      1225.,      2500.,
c j=5:
     1        0.0    ,  350.,       0.0  ,       793.,      2500.,
     2        4.24591,  452.,       0.0  ,      1225.,      2500.,
     3        0.0    ,  350.,       0.0  ,       793.,      2500.,
     4        0.0    ,  350.,       0.0  ,       793.,      2500.,
c j=6:
     1        2.3    ,  452.,       0.0  ,      1225.,      2500.,
     2        2.3    ,  452.,       0.0  ,      1225.,      2500.,
     3        2.3    ,  452.,       0.0  ,      1225.,      2500.,
     4        2.3    ,  452.,       0.0  ,      1225.,      2500.,
c
c neutron-proton potential
c ------------------------
c j=0:
     1        3.96451,  452.,       22.50007,   1225.,      2500.,
     2        0.0    ,  350.,       0.0     ,    793.,      2500.,
     3        7.866  ,  560.,       5.8     ,   1225.,      2500.,
     4        0.0    ,  452.,       0.0     ,   1225.,      2500.,
c j=1:
     1        0.81   ,  350.,       71.5    ,   1225.,      2500.,
     2        2.346  ,  424.,       19.22   ,   1225.,      2500.,
     3        0.575  ,  350.,       0.0     ,    793.,      2500.,
     4        0.51673,  350.,       14.01164,    793.,      2500.,
c j=2:
     1        2.236  ,  400.,       189.7   ,   1225.,      2500.,
     2        0.53   ,  350.,       154.5   ,   1225.,      2500.,
     3        1.573  ,  452.,       56.21   ,    793.,      2500.,
     4        4.194  ,  470.,       24.562  ,   1225.,      2500.,
c j=3:
     1        0.73   ,  350.,       0.0     ,    793.,      2500.,
     2        1.53   ,  452.,       74.85   ,    793.,      2500.,
     3        0.29   ,  350.,       0.0     ,    793.,      2500.,
     4        3.4    ,  452.,       0.0     ,    793.,      2500.,
c j=4:
     1        4.27591,  452.,       0.0     ,   1225.,      2500.,
     2        0.62   ,  350.,       0.0     ,    793.,      2500.,
     3        3.85   ,  452.,       17.61   ,   1225.,      2500.,
     4        3.8115 ,  452.,       17.61   ,   1225.,      2500.,
c j=5:
     1        0.51673,  350.,       0.0     ,    793.,      2500.,
     2        4.24591,  452.,       0.0     ,   1225.,      2500.,
     3        0.96   ,  350.,       0.0     ,    793.,      2500.,
     4        0.96   ,  350.,       0.0     ,    793.,      2500.,
c j=6:
     1        2.3    ,  452.,       0.0     ,   1225.,      2500.,
     2        2.3    ,  452.,       0.0     ,   1225.,      2500.,
     3        2.3    ,  452.,       0.0     ,   1225.,      2500.,
     4        2.3    ,  452.,       0.0     ,   1225.,      2500.,
c
c neutron-neutron potential
c -------------------------
c j=0:
     1        4.26338,  452.,       17.540,     1225.,      2500.,
     2        0.0    ,  350.,       0.0   ,      793.,      2500.,
     3        7.892  ,  560.,       16.747,     1225.,      2500.,
     4        0.0    ,  452.,       0.0   ,     1225.,      2500.,
c j=1:
     1        0.0    ,  350.,       0.0   ,     1225.,      2500.,
     2        2.326  ,  424.,       17.61 ,     1225.,      2500.,
     3        0.0    ,  350.,       0.0   ,      793.,      2500.,
     4        0.0    ,  350.,       0.0   ,      793.,      2500.,
c j=2:
     1        2.241  ,  400.,       190.7 ,     1225.,      2500.,
     2        0.0    ,  350.,       0.0   ,     1225.,      2500.,
     3        1.522  ,  452.,       56.28 ,      793.,      2500.,
     4        4.180  ,  470.,       24.737,     1225.,      2500.,
c j=3:
     1        0.0    ,  350.,       0.0   ,      793.,      2500.,
     2        1.53   ,  452.,       74.44 ,      793.,      2500.,
     3        0 0    ,  350.,       0.0   ,      793.,      2500.,
     4        0.0    ,  452.,       0.0   ,      793.,      2500.,
c j=4:
     1        4.284  ,  452.,       0.0   ,     1225.,      2500.,
     2        0.0    ,  350.,       0.0   ,      793.,      2500.,
     3        3.83   ,  452.,       17.61 ,     1225.,      2500.,
     4        3.81   ,  452.,       17.61 ,     1225.,      2500.,
c j=5:
     1        0.0    ,  350.,       0.0   ,      793.,      2500.,
     2        4.24591,  452.,       0.0   ,     1225.,      2500.,
     3        0.0    ,  350.,       0.0   ,      793.,      2500.,
     4        0.0    ,  350.,       0.0   ,      793.,      2500.,
c j=6:
     1        2.3    ,  452.,       0.0   ,     1225.,      2500.,
     2        2.3    ,  452.,       0.0   ,     1225.,      2500.,
     3        2.3    ,  452.,       0.0   ,     1225.,      2500.,
     4        2.3    ,  452.,       0.0   ,     1225.,      2500./
c
c        this has been the end of all tables
c        -----------------------------------
c
c
c
c

      save

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
      write (kwrite,10011)
      write (kwrite,10008)
      write (kwrite,10008)
      write (kwrite,10017)
c hit går det bra!!     
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
c      label=nucnuc(inn)
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
      wnn(inter)=sqrt(wn1*wn2)
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
      subroutine obstrq (icase,max,mex)
c
c        obstrq computes the structure of one-boson-exchanges
c
c
      implicit complex*16 (a-h,o-z)
      real*8 dwn, wnn, fac
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
      common /cobq/   vj(32,50),c(20,50),fff,ff,f(52),aa(96),ai(19,15),
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
      logical index
      logical indiso
      data jj/-1/
      data index/.false./
      save
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
      call obaiq
c
c
c
c
   60 continue
c
      if (dble(c(mc,im)).eq.0.d0) go to 1095
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
      aj=dble(j)
      aj1=dble(j+1)
      d2j1=1.d0/dble(2*j+1)
      arjj1=sqrt(aj*aj1)
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
      subroutine obaiq
c
c        obaiq performs the integration over angle theta
c        (necessary for the partial wave decomposition)
c        in analytic form by using the Legendre functions of the
c        second kind.
c
c
      implicit complex*16 (a-h,o-z)
      real*8 dwn, wnn, fac
c
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
c
c
c        common block for all ob-subroutines
c
      common /cobq/   vj(32,50),c(20,50),fff,ff,f(52),aa(96),ai(19,15),
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
      logical index
      data jj/-1/
      data index/.false./
      save
c
c
c
c
      if (index) go to 50
      index=.true.
c
      sqr2=sqrt(2.d0)
c
c
c
c
c
   50 if (j.eq.jj) go to 70
      jj=j
c
c
      aj=dble(j)
      aj1=dble(j+1)
      dj1=1.d0/aj1
      ajdj1=aj*dj1
      aaj=sqrt(ajdj1)
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
      zstamm= dcmplx(1.,0.)
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
c      call legen22 (qj,qjp1,zzq1m,1,zstamm,zdelta)
     
c      qj = 0.5d0*log( (z+1.)/(z-1.) )
c      qjp1 =  0.5d0*z*log( (z+1.)/(z-1.) ) - 1.d0
c      zzq1m = z*z*qjp1 - 1.D0/3.D0
      call legen2(qj,qjp1,zzq1m,1,zstamm,zdelta)
         
      qjm1= 0.d0
c
      else
c
c      write(*,*) zstamm+zdelta
      call legen2 (qjm1,qj,zzq1m,j,zstamm,zdelta)
c      write(*,*) zstamm+zdelta
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
      implicit complex*16 (a-h,o-z)
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
      dimension c(2,40001)

      data tolr/1.d-16/
      data me/40000/
      data jj/-1/
      save
c
c
c**** berechnung des arguments ****
      z=zstamm+zdelta
c
      qjm1=0.d0
      qj=0.d0
      zzq1m=0.d0
      if (j.lt.0) go to 123

c
c
c**** fallunterscheidung ****
      if (j.ne.0) go to 2
      if (ABS(z)-10.d0) 10,10,11
    2 if (j.ne.1) go to 3
      if (ABS(z)-1.5d0) 10,10,11
    3 if (j.ne.2) go to 4
      if (ABS(z)-1.2d0) 10,10,11
    4 zcut=1.d0+dble(j)**(-2.d0)
      if (ABS(z)-ABS(zcut)) 10,10,11
c
c**** rekursive berechnung mit dem logarithmus ****
   10 zdel=zstamm-1.d0
      zdel  =zdel+zdelta
      zz=2.d0/zdel  +1.d0
      qjm1=0.5d0*log(zz)
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

      qq=dble(2*i-1)/dble(i)*z*qj-dble(i-1)/dble(i)*qjm1
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
      if (ABS(cz).lt.ABS(tolr*zzq1m)) go to 62
      go to 13
      end if
      if (ABS(cz).lt.ABS(tolr*qj)) go to 62
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
      c(2,m)=1.d0/dble(2*m+1)
      cz1=c(1,m)*zzqinv
      cz= c(2,m)*zzqinv
      qjm1=cz1+qjm1
      qj=  cz+qj
      if (j.eq.1) then
      zzq1m=zzq1m+cz
      if (ABS(cz).lt.ABS(tolr*zzq1m)) go to 61
      go to 21
      end if
      if (ABS(cz).lt.ABS(tolr*qj)) go to 61
   21 zzqinv=zzqinv*zqinv
      go to 60
c
c**** fall j gleich zwei ****
   30 do 31 m=ma,me
      m2=2*m
      c(1,m)=1.d0/dble(m2+1)
      c(2,m)=c(1,m)*dble(m2)/dble(m2+3)
      qjm1=   c(1,m)*zzqinv+qjm1
      cz=     c(2,m)*zzqinv
      qj=cz+qj
      if (ABS(cz).lt.ABS(tolr*qj)) go to 61
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
   42 aehler=aehler*dble(k)
c**** nenner ****
      aenner=1.d0
      ka=m2+j-1
      ken=m2+2*j-3
      do 43 k=ka,ken,2
   43 aenner=aenner*dble(k)
      c(1,m)=aehler/aenner
      c(2,m)=c(1,m)*dble(kez+2)/dble(ken+2)
      qjm1=   c(1,m)*zzqinv+qjm1
      cz=     c(2,m)*zzqinv
      qj=cz+qj
      if (ABS(cz).lt.ABS(tolr*qj)) go to 61
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
   52 aehler=aehler*dble(k)
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
   53 aenner=aenner*dble(k)
      if (m2) 57,56,57
   57 c(2,m)=aehler/aenner
      qjm1=   c(1,m)*zzqinv+qjm1
      cz=     c(2,m)*zzqinv
      qj=cz+qj
      if (ABS(cz).lt.ABS(tolr*qj)) go to 61
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

  123 write (kwrite,1230)
 1230 format (/////' error in legen2. the parameter j of the'/
     1' Legendre function of the second kind is smaller zero.'/
     2' the function is set to zero.'/
     3' results may be wrong.'/////)
      return
      end
c********************************************************************
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
   10 if (abs(biga)-abs(a(ij)))  15,20,20
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
   45 if(ABS(biga)) 48,46,48
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
c**************** this is the end of the program cdbonn **********************




