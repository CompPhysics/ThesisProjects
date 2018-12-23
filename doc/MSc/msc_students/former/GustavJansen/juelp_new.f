*    there are two variables not defined (eerel,pmev)
*    i don't think it is relevant to our problem but it should be
*    checked every time
*    calling arguments same as nijmegen,  set ipole=0 
*
      subroutine juelp(jcalc,q,n,ipole)
      implicit real*8(a-h,o-z)
      common/info/ivvsm,iaasm,iap,iasp,ibp,ibsp,iwap,iwaap,
     1 isp,iwkp,iwkkp
      dimension s(96),u(96),eq1(96),eq2(96),eq3(96)
      dimension q(1),elab(21),elabs(21)
      common/print/iprint
       common /hpot/     vv11(70,70,8),vv12(70,70,8),
     1           vv21(70,70,8),vv22(70,70,8),
     2           vv(70,70,8)
      common /crdwrt/ kread,kwrite
      common /ctrans/ v11(8),v12(8),v21(8),v22(8),q0mev,q0smev,itra
      common /cpot2/  v2(8)
      common /cpot/   v(8),xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep
      common /cpoted/ zrel,noced
      common /cwnpwy/ wnpwl,wnpws
      common /ctyp/   ibary,indhal(18,90),indone(18,90)
c
      common /cobtb/  iob11,i11ds,i11ny,
     1                iob12,i12ds,i12ny,
     2                iob21,i21ds,i21ny,
     3                iob22,i22ds,i22ny,i22dl
c
c
c        specifications for these arguments
      logical heform,sing,trip,coup,csb,endep
      logical noced
      logical indhal,indone
c
c
c        further specifications
c
      real*8     uf/197.3271d0/
      real*8     pih/1.57079632679490d0/
      real*8     pi/3.14159265358979d0/
c        the following dimension for at most 20 elabs
      integer    memax/21/
      real*8     q0(20),q0q(20),q0s(20),q0sq(20),eq0(20),eq01(20),
     1           eq02(20),eq03(20),eq0s1(20)
c        the following dimension for at most j up to 19
      integer    nj(20)
c:
      character*10  name
      character*70  nname
      character*14  mname
      integer    nalt/-1/
c
      character*1 multi(4)/'1',3*'3'/
c
      logical indwrt/.false./
      logical ihead/.false./
      logical indj
      save
c
cccccccccccccccccccccccccccccccccccccc
c     suppress underflows
c
c     call xuflow(0)
ccccccccccccccccccccccccccccccccccccccc
c
10000 format (a10,20i3)
10001 format (1h0,a10,20i3)
10006 format (a14,20i3)
10007 format (1h0,a14,20i3)
10002 format (a10,6f10.4)
10003 format (1h0,a10,6f10.4)
c
10008 format (a70)
10009 format (1h //1h ,a70/1h ,16(1h-))
10014 format (1h ,a10,6f10.4)
10020 format (1h ,'number of gauss points nj = ',i3,' too high for given
     1 array dimensions'/1h ,'execution terminated')
10058 format (1h ,'order of particle channels:'//1h ,20x,
     1            '11 : nucleon lambda ---> nucleon lambda '/1h ,20x,
     2            '12 : nucleon sigma  ---> nucleon lambda '/1h ,20x,
     3            '21 : nucleon lambda ---> nucleon sigma  '/1h ,20x,
     4            '22a: nucleon sigma  ---> nucleon sigma  (i=1/2)'/1h ,
     5        20x,'22b: nucleon sigma  ---> nucleon sigma  (i=3/2)')
11015 format (/10x,'q'' = q0',5x,'q = ',f14.6,' mev')
11016 format (/10x,'q'' =',f14.6,' mev',5x,'q = q0')
11017 format (/10x,'q'' = q0',5x,'q = q0')
11021 format (/10x,'q'' =',f14.6,' mev',5x,'q =',f14.6,' mev')
11022 format (1h ,'v',a3,'(q'',q) =',8e14.6)
11020 format (//10x,'off-shell potential (1/mev**2) :'
     1             /t20,'0',t34,'1',t48,'+',
     2              t62,'-',t76,'+-',t90,'-+',t104,'st',t118,'ts'/)
11025 format (//10x,'half-off-shell potential (1/mev**2) :'
     1             /t20,'0',t34,'1',t48,'+',
     2              t62,'-',t76,'+-',t90,'-+',t104,'st',t118,'ts'/)
11026 format (//10x ,'on-shell potential (1/mev**2) :'
     1             /t20,'0',t34,'1',t48,'+',
     2              t62,'-',t76,'+-',t90,'-+',t104,'st',t118,'ts'/)
11030 format (///1h ,'j =',i2,t9,'lambda-elab = ',f6.2,' mev',3x,
     1 'lambda on-shell momentum = ',f12.4,' mev')
11031 format (1h ,t9,' sigma-elab = ',f6.2,' mev',3x,
     1 ' sigma on-shell momentum = ',f12.4,' mev')
11040 format (1h ,'je exceeds limit (19)')
12000 format (3x,'elab',i2,'  below sigma threshold')
12001 format (3x,'sigma-elab(mev)',i2,'  =',f10.4)
c
c        read and write input parameters for this program
c
c
c
      print *,' inside JUELP'
      read  (kread, 10008) nname
      iprint=0


      if(iprint.eq.1) write (kwrite,10009) nname
c
      jb=jcalc
      je=jcalc
      n1=n+1
      jb1=jb+1
      je1=je+1
      jee1=je1
      jend=je
      if (jee1.gt.20) jee1=20
      do 20 j1=1,jee1
      if (nj(j1).ge.ivvsm) then
      if(iprint.eq.1) write(kwrite,10020) nj(j1)
      stop
      endif
   20 continue
      if(iprint.eq.1) write (kwrite,10001) name,(nj(j1),j1=1,jee1)
      read  (kread ,10000) name,ising,itrip,icoup
      if(iprint.eq.1) write (kwrite,10001) name,ising,itrip,icoup
      read  (kread ,10000) name,icsb
      if(iprint.eq.1) write (kwrite,10001) name,icsb
c
      read  (kread ,10002) name,c
      if(iprint.eq.1) write (kwrite,10003) name,c
      read  (kread ,10000) name,ibary
      if(iprint.eq.1) write (kwrite,10001) name,ibary
      read  (kread ,10002) name,wn,w1,w2
      if(iprint.eq.1) write (kwrite,10003) name,wn,w1,w2
      if(iprint.eq.1) write (kwrite,10003)
      do 1 k=1,memax
      read  (kread ,10002) name,elab(k)
      if(iprint.eq.1) write (kwrite,10014) name,elab(k)
      if (elab(k).eq.0.d0) go to 2
    1 continue
    2 melab=k-1
      mela=melab
c
      read  (kread ,10000) name,iwrite
      if(iprint.eq.1) write (kwrite,10001) name,iwrite
c
c
c     iob.., i..ds, i..ny : potential for channel .. is switched
c            i..dl          on (iob.. etc. non-eq.0) or
c                           off (iob.. etc. eq.0)
c     iob..  : one-boson-exchange
c     i..ds  : delta-sigma box-diagram
c     i..ny  : nucleon-y-star box-diagram
c     i..dl  : delta-lambda box-diagram
      read  (kread ,10006) mname,iob11,i11ds,i11ny
      if(iprint.eq.1) write (kwrite,10007) mname,iob11,i11ds,i11ny
      read  (kread ,10006) mname,iob12,i12ds,i12ny
      if(iprint.eq.1) write (kwrite,10007) mname,iob12,i12ds,i12ny
      read  (kread ,10006) mname,iob21,i21ds,i21ny
      if(iprint.eq.1) write (kwrite,10007) mname,iob21,i21ds,i21ny
      read  (kread ,10006) mname,iob22,i22ds,i22ny,i22dl
      if(iprint.eq.1) write (kwrite,10007) mname,iob22,i22ds,i22ny,i22dl
c
c
c
c
c
c************************************
c
c        preparation of constants
c
c
      sing=.false.
      trip=.false.
      coup=.false.
      csb =.false.
      if (ising.ne.0) sing=.true.
      if (itrip.ne.0) trip=.true.
      if (icoup.ne.0) coup=.true.
      if (icsb.ne.0)  csb=.true.
c
      if (iwrite.ne.0) indwrt=.true.
c
      heform=.false.
      endep=.false.
c*** endep is true if potential propagators depend on start energy ***
c
c**************************************************************
c
c
      wnh=wn*0.5d0
      wnq=wn*wn
      wp=wn*pih
      wmre=w1*w2/(w1+w2)
      wred=wmre
      wnpwl=wn+w1
      wnpws=wn+w2
      wred1=w1*wn/(w1+wn)
      wred2=w2*wn/(w2+wn)
      wdiff=w1-w2
c
c
c        prepare energies
c
      zf1=2.d0*w1
      zf2=(wn+w1)**2
      wns2=4.d0*w2*w2*wn*wn
      zf3=w2*w2
      zf4=wnq+zf3
      do 10 k=1,melab
      q0q(k)=wnq*elab(k)*(elab(k)+zf1)/(zf2+2.d0*wn*elab(k))
      q0(k)=dsqrt(q0q(k))
      eq01(k)=dsqrt(q0q(k)+wn**2)
      eq02(k)=dsqrt(q0q(k)+w1**2)
      ens=eq01(k)+eq02(k)
      enss=ens*ens-wn*wn-w2*w2
      q0sq(k)=(enss*enss-wns2)/ens/ens/4.d0
      if (q0sq(k).lt.0.d0) q0sq(k)=0.d0
      q0s(k)=dsqrt(q0sq(k))
      eq0s1(k)=dsqrt(q0sq(k)+wn**2)
      eq03(k)=dsqrt(q0sq(k)+w2**2)
      elabs(k)=q0sq(k)/wn-w2+dsqrt(q0sq(k)**2/wnq+zf3+q0sq(k)*zf4/wnq)
      if (q0sq(k).eq.0.d0) then
      if(iprint.eq.1) write(6,12000) k
      elabs(k)=0.d0
      else
      if(iprint.eq.1) write(6,12001) k,elabs(k)
      endif
   10 continue
c
c
c
c        loop of total angular momentum j
c        --------------------------------
c        --------------------------------
c
      do 2000 j1=jb1,je1
c
c
      indj=.false.
      j2=j1+1
      j=j1-1
      ja=j
      aj=dfloat(j)
      aj1=dfloat(j+1)
      a2j1=dfloat(2*j+1)
      d2j1=1.d0/a2j1
      arjj1=dsqrt(aj*aj1)
      aaj=arjj1
c
c
      if (j1.gt.20) then
      if(iprint.eq.1) write(kwrite,11040)
      stop
      endif
  
  
  
  
  
  
c
c        loop of elabs
c        -------------
c        -------------
c
  300 do 1000 k=1,melab
c
c
c        define starting energy
c
      noced=.false.
c... we evaluate zrel in diswaveln --> NO! l'he activat (01.03.06) 
      zrel=eq01(k)+eq02(k)
      q0qmev=q0q(k)
      q0mev=q0(k)
      q0smev=q0s(k)
c here he anadido una c
c     q(n1)=q0(k)
c
c
      if (ipole.ne.0) go to 3148
c
c        compute potential matrix
c        ------------------------
c
c
      itra=1
c
      do 5401 ix=1,n
c
      xmev=q(ix)
c
      do 5401 iy=ix,n
c
      ymev=q(iy)
c
      call transi
c
      do 5401 ii=1,8
      vv11(ix,iy,ii)=v11(ii)
      vv12(ix,iy,ii)=v12(ii)
      vv21(ix,iy,ii)=v21(ii)
      vv22(ix,iy,ii)=v22(ii)
      vv  (ix,iy,ii)=v2 (ii)
 5401 continue

      do 5700 ix=2,n
      do 5700 iy=1,ix-1
      do 5700 iia=1,2
      do 5700 iib=1,4
      ii=(iia-1)*4+iib
      iii=ii+(iia-1)*(2*mod(ii,2)-1)
      vv11(ix,iy,ii)=vv11(iy,ix,iii)
      vv12(ix,iy,ii)=vv21(iy,ix,iii)
      vv21(ix,iy,ii)=vv12(iy,ix,iii)
      vv22(ix,iy,ii)=vv22(iy,ix,iii)
      vv  (ix,iy,ii)=vv  (iy,ix,iii)
 5700 continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     up to here:    < n y'| v | n y >
c
c     now for convenience: transormation to < y' n | v | y n >
c
c     leading to a phase factor:
c      (-1) for a sing-trip transition
c                 times
c      (-1) for a lamb-sig transition
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do 5800 ix=1,n
      do 5800 iy=1,n
c
      do 5810 ii=7,8
      vv11(ix,iy,ii)=-vv11(ix,iy,ii)
      vv22(ix,iy,ii)=-vv22(ix,iy,ii)
      vv  (ix,iy,ii)=-vv  (ix,iy,ii)
 5810 continue
c
      do 5820 ii=1,6
      vv12(ix,iy,ii)=-vv12(ix,iy,ii)
      vv21(ix,iy,ii)=-vv21(ix,iy,ii)
 5820 continue
c
 5800 continue
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     o u t p u t
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (.not.indwrt) go to 1000
      if (ihead) go to 5922
      ihead=.true.
c
      if(iprint.eq.1) write (kwrite,10058)
c
c
c
c
 5922 continue
      if(iprint.eq.1) write(6,11030) j,elab(k),q0(k)
      if(iprint.eq.1) write(6,11031)   elabs(k),q0s(k)
      if(iprint.eq.1) write (kwrite,11020)
c
      do 5930 ix=1,n
      do 5930 iy=1,n
  
      if(iprint.eq.1) write(kwrite,11021) q(ix),q(iy)
c
      if(iprint.eq.1) write(kwrite,11022) '11 ',(vv11(ix,iy,ii),ii=1,8)
      if(iprint.eq.1) write(kwrite,11022) '12 ',(vv12(ix,iy,ii),ii=1,8)
      if(iprint.eq.1) write(kwrite,11022) '21 ',(vv21(ix,iy,ii),ii=1,8)
      if(iprint.eq.1) write(kwrite,11022) '22a',(vv22(ix,iy,ii),ii=1,8)
      if(iprint.eq.1) write(kwrite,11022) '22b',(vv  (ix,iy,ii),ii=1,8)
c
 5930 continue
      goto 1000 
 3148 continue
c        compute potential vector
c        ------------------------
c
c
 5500 itra=2
      ix=n1
c
c
      do 5501 iy=1,n
c
c
      ymev=q(iy)
c
c
      call transi
c
c
c
      do 5501 ii=1,8
      vv11(ix,iy,ii)=v11(ii)
      vv12(ix,iy,ii)=v12(ii)
      vv21(ix,iy,ii)=v21(ii)
      vv22(ix,iy,ii)=v22(ii)
      vv  (ix,iy,ii)=v2 (ii)
 5501 continue
c
c
c        compute potential element
c        -------------------------
c
c
 5510 itra=3
      ix=n1
      iy=n1
c
c
      call transi
c
c
      do 5601 ii=1,8
      vv11(ix,iy,ii)=v11(ii)
      vv12(ix,iy,ii)=v12(ii)
      vv21(ix,iy,ii)=v21(ii)
      vv22(ix,iy,ii)=v22(ii)
      vv  (ix,iy,ii)=v2 (ii)
 5601 continue
c
c
c
      iy=n1
      do 5710 ix=1,n
      do 5710 iia=1,2
      do 5710 iib=1,4
      ii=(iia-1)*4+iib
      iii=ii+(iia-1)*(2*mod(ii,2)-1)
      vv11(ix,iy,ii)=vv11(iy,ix,iii)
      vv12(ix,iy,ii)=vv21(iy,ix,iii)
      vv21(ix,iy,ii)=vv12(iy,ix,iii)
      vv22(ix,iy,ii)=vv22(iy,ix,iii)
      vv  (ix,iy,ii)=vv  (iy,ix,iii)
 5710 continue
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     up to here:    < n y'| v | n y >
c
c     now for convenience: transormation to < y' n | v | y n >
c
c     leading to a phase factor:
c      (-1) for a sing-trip transition
c                 times
c      (-1) for a lamb-sig transition
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      iy=n1
      do 8900 ix=1,n1
c
      do 8910 ii=7,8
      vv11(ix,iy,ii)=-vv11(ix,iy,ii)
      vv22(ix,iy,ii)=-vv22(ix,iy,ii)
      vv  (ix,iy,ii)=-vv  (ix,iy,ii)
 8910 continue
c
      do 8920 ii=1,6
      vv12(ix,iy,ii)=-vv12(ix,iy,ii)
      vv21(ix,iy,ii)=-vv21(ix,iy,ii)
 8920 continue
c
 8900 continue

      ix=n1
      do 8930 iy=1,n1
c
      do 8940 ii=7,8
      vv11(ix,iy,ii)=-vv11(ix,iy,ii)
      vv22(ix,iy,ii)=-vv22(ix,iy,ii)
      vv  (ix,iy,ii)=-vv  (ix,iy,ii)
 8940 continue
c
      do 8950 ii=1,6
      vv12(ix,iy,ii)=-vv12(ix,iy,ii)
      vv21(ix,iy,ii)=-vv21(ix,iy,ii)
 8950 continue
c
 8930 continue
c
c
c
      if(iprint.eq.1) write (kwrite,11025)
  
      ix=n1
      do 5931 iy=1,n
  
      if(iprint.eq.1) write(kwrite,11015) q(iy)
c
      if(iprint.eq.1) write(kwrite,11022) '11 ',(vv11(ix,iy,ii),ii=1,8)
      if(iprint.eq.1) write(kwrite,11022) '12 ',(vv12(ix,iy,ii),ii=1,8)
      if(iprint.eq.1) write(kwrite,11022) '21 ',(vv21(ix,iy,ii),ii=1,8)
      if(iprint.eq.1) write(kwrite,11022) '22a',(vv22(ix,iy,ii),ii=1,8)
      if(iprint.eq.1) write(kwrite,11022) '22b',(vv  (ix,iy,ii),ii=1,8)
c
 5931 continue
c
      iy=n1
      do 5932 ix=1,n
  
      if(iprint.eq.1) write(kwrite,11016) q(ix)
c
      if(iprint.eq.1) write(kwrite,11022) '11 ',(vv11(ix,iy,ii),ii=1,8)
      if(iprint.eq.1) write(kwrite,11022) '12 ',(vv12(ix,iy,ii),ii=1,8)
      if(iprint.eq.1) write(kwrite,11022) '21 ',(vv21(ix,iy,ii),ii=1,8)
      if(iprint.eq.1) write(kwrite,11022) '22a',(vv22(ix,iy,ii),ii=1,8)
      if(iprint.eq.1) write(kwrite,11022) '22b',(vv  (ix,iy,ii),ii=1,8)
c
 5932 continue
c
c
      if(iprint.eq.1) write (kwrite,11026)
c
  
      if(iprint.eq.1) write(kwrite,11017)
c
      if(iprint.eq.1) write(kwrite,11022) '11 ',(vv11(n1,n1,ii),ii=1,8)
      if(iprint.eq.1) write(kwrite,11022) '12 ',(vv12(n1,n1,ii),ii=1,8)
      if(iprint.eq.1) write(kwrite,11022) '21 ',(vv21(n1,n1,ii),ii=1,8)
      if(iprint.eq.1) write(kwrite,11022) '22a',(vv22(n1,n1,ii),ii=1,8)
      if(iprint.eq.1) write(kwrite,11022) '22b',(vv  (n1,n1,ii),ii=1,8)
c
c
c
c
 1000 continue
c        this has been the end of the elab loop
c
 2000 continue
c        this has been the end of the j loop
c
c
  999 return
      end
c***********************************************************************
c
c***********************************************************************
      subroutine transi
c
c
c        routine calls transition-potentials
c
c
      implicit real*8 (a-h,o-z)
c
c
      common/print/iprint
      common /ctrans/ v11(8),v12(8),v21(8),v22(8),q0mev,q0smev,itra
      common /cpot2/  v2(8)
      common /cpot/   v(8),xmev,ymev
      common /cpoted/ zrel,noced
      common /cobtb/  iob11,i11ds,i11ny,
     1                iob12,i12ds,i12ny,
     2                iob21,i21ds,i21ny,
     3                iob22,i22ds,i22ny,i22dl
c
      logical noced,nnoced
      save
      nnoced=noced
c
      do 5 ivv=1,8
      v11(ivv)=0.d0
      v12(ivv)=0.d0
      v21(ivv)=0.d0
      v22(ivv)=0.d0
    5 v2(ivv)=0.d0
c
c
c
c
      if (iob11.eq.0.and.i11ds.eq.0.and.i11ny.eq.0) go to 100
      if (itra.eq.2) then
      xmev=q0mev
      else if (itra.eq.3) then
      xmev=q0mev
      ymev=q0mev
      end if
c
      call pot11
      noced=nnoced
      do 11 i=1,8
      v11(i)=v(i)
   11 continue
  100 continue
c
      if (itra.eq.3.and.q0smev.eq.0.d0) then
      do 101 iv=1,8
      v12(iv)=0.d0
      v21(iv)=0.d0
      v22(iv)=0.d0
  101 v2(iv)=0.d0
      go to 999
      end if
c
c
      if (iob12.eq.0.and.i12ds.eq.0.and.i12ny.eq.0) go to 200
      if (itra.eq.2) then
      xmev=q0mev
      else if (itra.eq.3) then
      xmev=q0mev
      ymev=q0smev
      end if
c
      call pot12
      noced=nnoced
      do 12 i=1,8
      v12(i)=v(i)
   12 continue
c
c
c
  200 if (iob21.eq.0.and.i21ds.eq.0.and.i21ny.eq.0) go to 300
      if (itra.eq.2.and.q0smev.gt.0.d0) then
      xmev=q0smev
      else if (itra.eq.3) then
      xmev=q0smev
      ymev=q0mev
      end if
c
      if (itra.eq.2.and.q0smev.eq.0.d0) then
      do 210 iv=1,8
  210 v21(iv)=0.d0
      go to 300
      end if
c
      call pot21
      noced=nnoced
      do 21 i=1,8
      v21(i)=v(i)
   21 continue
c
c
  300 if (iob22.eq.0.and.i22ds.eq.0.and.i22ny.eq.0
     1                             .and.i22dl.eq.0) go to 999
      if (itra.ne.1.and.q0smev.eq.0.d0) then
      do 310 iv=1,8
      v22(iv)=0.d0
  310 v2(iv)=0.d0
      go to 999
      end if
c
      if (itra.eq.2) then
      xmev=q0smev
      else if (itra.eq.3) then
      xmev=q0smev
      ymev=q0smev
      end if
c
      call pot22
      noced=nnoced
      do 22 i=1,8
      v22(i)=v(i)
   22 continue
c
  999 return
      end
c***********************************************************************
c
c***********************************************************************
      subroutine pot11
c
c
c        routine calls one and two boson exchange diagrams
c        for channel 11
c
      implicit real*8 (a-h,o-z)
c
      common/print/iprint
      common /cpot/   v(8),xmev,ymev
      common /cpot2/  v2(8)
      common /cstate/ j,heform,sing,trip,coup,endep
      common /cwnpwy/ wnpwl,wnpws
      common /cwingo/ wingo
      common /cobtb/  iob11,i11ds,i11ny,
     1                iob12,i12ds,i12ny,
     2                iob21,i21ds,i21ny,
     3                iob22,i22ds,i22ny,i22dl
c
      dimension vv(8)
c
      logical heform,sing,trip,coup,endep
      save
c
      wingo=wnpwl
c
      do 5 iv=1,8
    5 vv(iv)=0.d0
c
      if (iob11.eq.0) go to 15
c
c     nucleon lambda   (inter = 1)
c     ----------------------------
      call obep(1)
      do 10 i=1,8
      vv(i)=vv(i)+v(i)
   10 continue
c
c
c
   15 if (i11ds.eq.0) go to 25
c
c     delta sigma box  (inter = 5)
c     ----------------------------
      call tbep(5)
      do 20 i=1,8
      vv(i)=vv(i)+v(i)
   20 continue
c
c
c
   25 if (i11ny.eq.0) go to 35
c
c     nucleon y-star box  (inter = 7)
c     -------------------------------
      call tbep(7)
      do 30 i=1,8
      vv(i)=vv(i)+v(i)
   30 continue
   35 continue
c
c
      do 50 iv=1,8
      v2(iv)=0.d0
   50 v(iv)=vv(iv)
c
ccccc
c
c
c
c
  999 continue
      return
      end
c***********************************************************************
c
c***********************************************************************
      subroutine pot12
c
c
c        routine calls one and two boson exchange diagrams
c        for channel 12
c
      implicit real*8 (a-h,o-z)
c
c
      common/print/iprint
      common /cpot/   v(8),xmev,ymev
      common /cpot2/  v2(8)
      common /cstate/ j,heform,sing,trip,coup,endep
      common /cwnpwy/ wnpwl,wnpws
      common /cwingo/ wingo
      common /cobtb/  iob11,i11ds,i11ny,
     1                iob12,i12ds,i12ny,
     2                iob21,i21ds,i21ny,
     3                iob22,i22ds,i22ny,i22dl
c
      dimension vv(8)
c
      logical heform,sing,trip,coup,endep
      save
c
      wingo=(wnpws+wnpwl)/2.d0
c
      do 5 iv=1,8
    5 vv(iv)=0.d0
c
      if (iob12.eq.0) go to 15
c
c        (inter = 2)
c     ----------------------------
      call obep(2)
      do 10 i=1,8
      vv(i)=vv(i)+v(i)
   10 continue
c
c
c
   15 if (i12ds.eq.0) go to 25
c
c     delta sigma box  (inter = 11)
c     ----------------------------
      call tbep(11)
      do 20 i=1,8
      vv(i)=vv(i)+v(i)
   20 continue
c
c
c
   25 if (i12ny.eq.0) go to 35
c
c     nucleon y-star box  (inter = 13)
c     -------------------------------
      call tbep(13)
      do 30 i=1,8
      vv(i)=vv(i)+v(i)
   30 continue
   35 continue
c
c
      do 50 iv=1,8
      v2(iv)=0.d0
   50 v(iv)=vv(iv)
c
ccccc
c
c
c
c
  999 continue
      return
      end
c***********************************************************************
c
c***********************************************************************
      subroutine pot21
c
c
c        routine calls one and two boson exchange diagrams
c        for channel 21
c
      implicit real*8 (a-h,o-z)
c
c
      common/print/iprint
      common /cpot/   v(8),xmev,ymev
      common /cpot2/  v2(8)
      common /cstate/ j,heform,sing,trip,coup,endep
      common /cwnpwy/ wnpwl,wnpws
      common /cwingo/ wingo
      common /cobtb/  iob11,i11ds,i11ny,
     1                iob12,i12ds,i12ny,
     2                iob21,i21ds,i21ny,
     3                iob22,i22ds,i22ny,i22dl
c
      dimension vv(8)
c
      logical heform,sing,trip,coup,endep
      save
c
      wingo=(wnpwl+wnpws)/2.d0
c
      do 5 iv=1,8
    5 vv(iv)=0.d0
c
      if (iob21.eq.0) go to 15
c
c        (inter = 3)
c     ----------------------------
      call obep(3)
      do 10 i=1,8
      vv(i)=vv(i)+v(i)
   10 continue
c
c
c
   15 if (i21ds.eq.0) go to 25
c
c     delta sigma box  (inter = 15)
c     ----------------------------
      call tbep(15)
      do 20 i=1,8
      vv(i)=vv(i)+v(i)
   20 continue
c
c
c
   25 if (i21ny.eq.0) go to 35
c
c     nucleon y-star box  (inter = 17)
c     -------------------------------
      call tbep(17)
      do 30 i=1,8
      vv(i)=vv(i)+v(i)
   30 continue
   35 continue
c
c
      do 50 iv=1,8
      v2(iv)=0.d0
   50 v(iv)=vv(iv)
c
ccccc
c
c
c
c
  999 continue
      return
      end
c***********************************************************************
c
c***********************************************************************
      subroutine pot22
c
c
c        routine calls one and two boson exchange diagrams
c        for channel 22
c
      implicit real*8 (a-h,o-z)
c
c
      common/print/iprint
      common /cpot2/  v2(8)
      common /cpot/   v(8),xmev,ymev
      common /cwnpwy/ wnpwl,wnpws
      common /cwingo/ wingo
      common /cstate/ j,heform,sing,trip,coup,endep
      common /cobtb/  iob11,i11ds,i11ny,
     1                iob12,i12ds,i12ny,
     2                iob21,i21ds,i21ny,
     3                iob22,i22ds,i22ny,i22dl
c
      dimension vv(8),vv2(8)
c
      logical heform,sing,trip,coup,endep
      save
c
      wingo=wnpws
c
      do 5 iv=1,8
      vv(iv)=0.d0
    5 vv2(iv)=0.d0
c
      if (iob22.eq.0) go to 15
c
c     nucleon sigma  (i=1/2)   (inter = 4)
c     ------------------------------------
      call obep(4)
      do 10 i=1,8
      vv(i)=vv(i)+v(i)
   10 continue
c
c
c
c     nucleon sigma  (i=3/2)  (inter = 9)
c     ------------------------------------
      call obep(9)
      do 11 i=1,8
      vv2(i)=vv2(i)+v2(i)
   11 continue
c
c
c
c
c
c
   15 if (i22ds.eq.0) go to 25
c
c     delta sigma box  (inter = 6)
c     ----------------------------
      call tbep(6)
      do 20 i=1,8
      vv(i)=vv(i)+v(i)
      vv2(i)=vv2(i)+v2(i)
   20 continue
c
c
c
   25 if (i22ny.eq.0) go to 35
c
c     nucleon y-star box  (inter = 8)
c     -------------------------------
      call tbep(8)
      do 30 i=1,8
      vv(i)=vv(i)+v(i)
      vv2(i)=vv2(i)+v2(i)
   30 continue
c
c
c
c
   35 if (i22dl.eq.0) go to 45
c
c     delta lambda (i=3/2)   (inter = 10)
c     -----------------------------------
      call tbep(10)
      do 40 i=1,8
      vv2(i)=vv2(i)+v2(i)
   40 continue
   45 continue
c
c
c
      do 50 iv=1,8
      v(iv)=vv(iv)
   50 v2(iv)=vv2(iv)
c
ccccc
c
c
c
c
  999 continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE OBEP(INTERA)
C
C
C        ONE-BOSON-EXCHANGE NY-NY INTERACTION;
C        VERSION WHICH USES NUMERICAL INTEGRATION
C
C
C
C
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
c     older version contains the following comented line
c     common/print/iprint 
      COMMON /CRDWRT/ KREAD,KWRITE
C===================================================================
C        COMMON BLOCKS WHICH CONTAIN THE ARGUMENTS AND VALUES OF THE
C        OB-SUBROUTINES
C
      COMMON /CPOT2/   V2(8)
      COMMON /CPOT/   V(8),XMEV,YMEV
      COMMON /CSTATE/ J,HEFORM,SING,TRIP,COUP,ENDEP
      COMMON /CTYP/   IBARY,INDHAL(18,90),INDONE(18,90)
C
C        IN CASE OF ENERGY-DEPENDENCE LOOK FOR THE COMMON-BLOCK /CPOTED/
C        IN OBAI AND OBAA.
C
C
C
C        SPECIFICATIONS FOR THESE COMMON BLOCKS
C
      LOGICAL HEFORM,SING,TRIP,COUP,ENDEP
C
C=======================================================================
C        COMMON BLOCK FOR ALL OB-SUBROUTINES
C
      COMMON /COB2/   VJ2(32,90)
      COMMON /COB/    VJ(32,90),C(10,90),FFF,FF,F(52),AA(96),AI(19,5),
     1                WNN(18),X,XX,Y,YY,XY2,XXPYY,EX,EY,EEM12,
     2                EZ1,EZ2,CT(96),WT(96),WPPS(4,18),WEX,
     3                IC(10,90),IFT(18),MINT(18),MAXT(18),NT,
     4                MGE,MGG(11,18),MGGO(11,18),IMA(5,11,18),
     5                IMAA(18),IMEA(18),IME,IM,MC,M,MG,INTER,
     6                INDC(2,90),INDXY
      COMMON /CCOUL/  VCOUL(8),RCUT,JCOUL
C
C
C         SPECIFICATIONS FOR THIS COMMON BLOCK
C
      LOGICAL INDC,INDXY,INDHAL,INDONE
      DIMENSION WPS(4)
C
C
C        FURTHER SPECIFICATIONS
C
      CHARACTER*4 MESONG(11)/'0-  ','0-E ','0+E ','0+  ','1-AT',
     1                   '1-  ','1-T ','1-E ','1-TE',
     2                    '1-C ','2+  '/
      LOGICAL INDEX(18)
      LOGICAL INDMG(11)
      save
C
      DATA PI/3.14159265358979D0/
      DATA  D3/0.333333333333333D0/
      DATA TD3/0.666666666666667D0/
      DATA FD3/1.333333333333333D0/
      DATA INDEX/18*.FALSE./
      DATA INDMG/11*.FALSE./
C
C
C
C
      INTER=INTERA
C
C
C
C
C        CALL SUBROUTINE OBPAR ONCE AND ONLY ONCE
C
C
      IF (INDEX(INTERA)) GO TO 50
      INDEX(INTERA)=.TRUE.
C
C
      CALL OBPAR
C
C
   50 CONTINUE
      IFTGO=IFT(INTER)+1
      DWN=1.D0/WNN(INTER)
C
      DO 10 IN=1,4
   10 WPS(IN)=WPPS(IN,INTER)
C
      WY1=WPS(1)*WPS(1)
      WY2=WPS(2)*WPS(2)
      WY3=WPS(3)*WPS(3)
      WY4=WPS(4)*WPS(4)
      IMAN=IMAA(INTER)
      IMEN=IMEA(INTER)
C
C
C        PREPARE CONSTANT OVER-ALL FACTOR
C
      FAC=1.D0/(2.D0*PI)*DWN*DWN
C     --------------------------
C
C
C
C
C
C
C
C        PREPARE EXPRESSIONS DEPENDING ON X AND Y
C        ----------------------------------------
C        ----------------------------------------
C
C
C
C
      XA=XMEV*DWN
      YA=YMEV*DWN
      IF (XA.EQ.X.AND.YA.EQ.Y) GO TO 55
      INDXY=.FALSE.
      X=XA
      XX=X*X
      Y=YA
      YY=Y*Y
      XY2=X*Y*2.D0
      XXPYY=XX+YY
      EX=DSQRT(1.D0+XX)
      EY=DSQRT(1.D0+YY)
      EEM12=(EX*EY-1.D0)*2.D0
C
C
C
C
   55 XY=XY2*0.5D0
      DXY=1.D0/XY
      EX1=DSQRT(WY3+XX)
      EX2=DSQRT(WY4+XX)
      EY1=DSQRT(WY1+YY)
      EY2=DSQRT(WY2+YY)
c     in the older version the next 4 line were commented
CCCCCCCCCC  APPROXIMATION :  E = M  CCCCCCCCCCCC
c      EX1 = WPS(3)
c      EX2 = WPS(4)
c      EY1 = WPS(1)
c      EY2 = WPS(2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      EX=EX1+EX2
      EY=EY1+EY2
      W12=WPS(1)+WPS(2)
      W34=WPS(3)+WPS(4)
      EE=EY/W12*EX/W34
      REE=DSQRT(EE)
      EEM1=EE-1.D0
      EME=EX-EY
      EMEH=EME*0.5D0
      EMEHQ=EMEH*EMEH
      EEP1=EE+1.D0
       EPE=EX+EY
      XXYY=XX*YY
C
C
      EPX1=EX1+WPS(3)
      EPX2=EX2+WPS(4)
      EPY1=EY1+WPS(1)
      EPY2=EY2+WPS(2)
C
      OVF1=1.D0+XX/EPX1/EPX2
      OVF2=1.D0+YY/EPY1/EPY2
      OVF3=1.D0-XX/EPX1/EPX2
      OVF4=1.D0-YY/EPY1/EPY2
      OVF5=XX/EPX1/EPX2+YY/EPY1/EPY2
      OVF6=1.D0/EPX1+1.D0/EPX2
      OVF7=1.D0/EPY1+1.D0/EPY2
      OVF8=1.D0/EPX1-1.D0/EPX2
      OVF9=1.D0/EPY1-1.D0/EPY2
      OVF10=1.D0/EPX1/EPY1+1.D0/EPX2/EPY2
      OVF11=1.D0/EPX1/EPY2+1.D0/EPX2/EPY1
C
C        PREPARE OVER-ALL FACTOR
C
C
      GO TO (70,71,72),IFTGO
C
C        NO ADDITIONAL FACTOR
C
   70 FFF=FAC
      GO TO 90
C
C        MINIMAL RELATIVITY
C
   71 FFF=FAC/REE
      GO TO 90
C
C        FACTOR (M1/E1*M2/E2*M3/E3*M4/E4)**1/2
C
C  72 FFF = FAC/EE
C
C
   72 EES1 = DSQRT(WPS(1)/EY1)
      EES2 = DSQRT(WPS(2)/EY2)
      EES3 = DSQRT(WPS(3)/EX1)
      EES4 = DSQRT(WPS(4)/EX2)
      EESS = EES1*EES2*EES3*EES4
      FFF = FAC*EESS
C
C
C
C
C
C
   90 DO 93 IV=1,8
      VCOUL(IV)=0.D0
      V2(IV)=0.D0
   93 V(IV)=0.D0
      DO 95 IL=IMAN,IMEN
      DO 95 IV=1,32
      VJ2(IV,IL)=0.D0
   95 VJ(IV,IL)=0.D0
C
C
C
C
C        CONTRIBUTIONS OF MESONS
C        -----------------------
C        -----------------------
C
C
C
C
      DO 1995 IMG=1,MGE
      MG=MGGO(IMG,INTER)
      IF (MG.EQ.0) GO TO 2000
      ME=MGG(MG,INTER)
C
C
      GO TO (100,200,300,400,500,600,720,601,721,1000,1100),MG
      GOTO 9000
C
C
C        0-  , PSEUDO-SCALAR MESONS (DIRECT DIAGRAM)
C        -------------------------------------------
C
C
C
C
  100 MC=1
C
      FF=EPX1*EPX2*EPY1*EPY2/WPS(1)/WPS(2)/WPS(3)/WPS(4)
      FF=DSQRT(FF)/2.D0
      PSM1=XX/EPX1/EPX2-YY/EPY1/EPY2
      PSM2=1.D0/EPX1/EPY2-1.D0/EPX2/EPY1
C
      F(1)=OVF5
      F(2)=-XY*OVF11
      F(3)=-F(1)
      F(4)=-F(2)
      F(5)=F(2)
      F(6)=F(1)
      F(7)=-PSM1
      F(8)=-F(7)
      F(12)=-XY*PSM2
      F(13)=-F(12)
C
      CALL OBSTR(1,1,ME)
      GO TO 1995
C
C
C
C        0-E , PSEUDO-SCALAR MESONS (EXCHANGE DIAGRAM)
C        ---------------------------------------------
C
C
C
C
  200 MC=1
C
      FF=EPX1*EPX2*EPY1*EPY2/WPS(1)/WPS(2)/WPS(3)/WPS(4)
      FF=DSQRT(FF)/2.D0
      PSM1=XX/EPX1/EPX2-YY/EPY1/EPY2
      PSM2=1.D0/EPX2/EPY2-1.D0/EPX1/EPY1
C
      F(1)=OVF5
      F(2)=-XY*OVF10
      F(3)=-F(1)
      F(4)=-F(2)
      F(5)=F(2)
      F(6)=F(1)
      F(7)=-PSM1
      F(8)=-F(7)
      F(12)=-XY*PSM2
      F(13)=-F(12)
C
      CALL OBSTR(3,1,ME)
      GO TO 1995
C
C
C
C        0-ST, PSEUDO-SCALAR MESONS IN STATIC LIMIT
C        ------------------------------------------
C
C     Range 300 corresponds to next lines in older version
C
C
c 300  MC=1
C
c      FF=1.D0
c      F(1)=XXPYY*0.5D0
c      F(2)=-XY
c      F(3)=-F(1)
c      F(4)=-F(2)
c      F(5)=F(2)
c      F(6)=F(1)
c      F(7)=-(XX-YY)*0.5D0
c      F(8)=-F(7)
C
c      CALL OBSTR(1,1,ME)
c      GO TO 1995
C
C
C        0+E , SCALAR MESONS  EXCHANGE-DIAGRAM 
C        -------------------
C
C
C
C
  300 MC=1
C
      FF=EPX1*EPX2*EPY1*EPY2/WPS(1)/WPS(2)/WPS(3)/WPS(4)
      FF=DSQRT(FF)/2.D0
      SM1=1.D0+XXYY/EPX1/EPX2/EPY1/EPY2
      SM2=1.D0-XXYY/EPX1/EPX2/EPY1/EPY2
      SM3=1.D0/EPX2/EPY1-1.D0/EPX1/EPY2
C
      F(1)=-SM1
      F(2)=XY*OVF11
      F(3)=F(1)
      F(4)=F(2)
      F(5)=F(2)
      F(6)=F(1)
      F(7)=SM2
      F(8)=F(7)
      F(12)=-XY*SM3
      F(13)=F(12)
C
      CALL OBSTR(3,1,ME)
      GO TO 1995
C
C
C
C
C        0+  , SCALAR MESONS
C        -------------------
C
C
C
C
  400 MC=1
C
      FF=EPX1*EPX2*EPY1*EPY2/WPS(1)/WPS(2)/WPS(3)/WPS(4)
      FF=DSQRT(FF)/2.D0
      SM1=1.D0+XXYY/EPX1/EPX2/EPY1/EPY2
      SM2=1.D0-XXYY/EPX1/EPX2/EPY1/EPY2
      SM3=1.D0/EPX1/EPY1-1.D0/EPX2/EPY2
C
      F(1)=-SM1
      F(2)=XY*OVF10
      F(3)=F(1)
      F(4)=F(2)
      F(5)=F(2)
      F(6)=F(1)
      F(7)=SM2
      F(8)=F(7)
      F(12)=-XY*SM3
      F(13)=F(12)
C
      CALL OBSTR(1,1,ME)
      GO TO 1995
C
C
C
C
C        1-AT, ADDITIONAL TERMS FOR VECTOR MESONS
C              IN CASE OF DIFFERENT MASSES
C        ----------------------------------------
C
C
C
  500 MC=1
C
C
C
      FF=EPX1*EPX2*EPY1*EPY2/WPS(1)/WPS(2)/WPS(3)/WPS(4)
      FF=DSQRT(FF)/2.D0
CCCCC FF=-(WPS(4)-WPS(1))*(WPS(3)-WPS(2))*FF
CCCCC FF=FF/WEX
C
      SM1=1.D0+XXYY/EPX1/EPX2/EPY1/EPY2
      SM2=1.D0-XXYY/EPX1/EPX2/EPY1/EPY2
      SM3=1.D0/EPX2/EPY1-1.D0/EPX1/EPY2
C
      F(1)=-SM1
      F(2)=XY*OVF11
      F(3)=F(1)
      F(4)=F(2)
      F(5)=F(2)
      F(6)=F(1)
      F(7)=SM2
      F(8)=F(7)
      F(12)=-XY*SM3
      F(13)=F(12)
C
      CALL OBSTR(3,1,ME)
      GO TO 1995
C
C
C
C
C        1-,  VECTOR MESONS  (DIRECT DIAGRAM)
C        ------------------------------------
C
C
C
C
C        VECTOR COUPLING (G1*G2)
C
C
C
C
  600 MC=1
C
      FF=EPX1*EPX2*EPY1*EPY2/WPS(1)/WPS(2)/WPS(3)/WPS(4)
      FF=DSQRT(FF)/2.D0
C     FF=0.D0
      VC1=OVF1*OVF2
C
      F(1)=VC1+2.D0*OVF5
      F(2)=XY*OVF8*OVF9
      F(3)=VC1
      F(4)=XY*OVF6*OVF7
      F(5)=XY*(OVF10+3.D0*OVF11)
      F(6)=OVF3*OVF4
      F(7)=-OVF3*OVF2
      F(8)=-OVF1*OVF4
      F(12)=-XY*OVF8*OVF7
      F(13)=-XY*OVF6*OVF9
C
C
      CALL OBSTR(1,1,ME)
C
C
C
C
C        TENSOR COUPLING (F1*F2)
C
C
C
C
      MC=2
C
C
      FF=EPX1*EPX2*EPY1*EPY2/WPS(1)/WPS(2)/WPS(3)/WPS(4)
      FF=DSQRT(FF)/8.D0
C     FF=0.D0
      OWF1=WPS(3)+WPS(1)
      OWF2=WPS(4)+WPS(2)
      OWF3=1.D0+XXYY/EPX1/EPX2/EPY1/EPY2
      OWF4=EX1+EX2+EY1+EY2-OWF1
      OWF5=EX1+EX2+EY1+EY2-OWF2
      OWF6=1.D0/EPX1/EPY1-1.D0/EPX2/EPY2
      OWF7=(EX1-WPS(3))*(EY1-WPS(1))
      OWF8=(EX2-WPS(4))*(EY2-WPS(2))
      OWF9=(EX2+EY2)/EPX1/EPY1+1.D0/EPX1+1.D0/EPY1
      OWF10=(EX1+EY1)/EPX2/EPY2+1.D0/EPX2+1.D0/EPY2
      OWF11=OWF9-OWF4/EPX2/EPY2
      OWF12=OWF10-OWF5/EPX1/EPY1
      OWF13=XX/EPX1/EPX2-YY/EPY1/EPY2
      OUF1=(EX1+EY1)*(EX2+EY2)+XXPYY
      OUF2=1.D0-XXYY/EPX1/EPX2/EPY1/EPY2
C
      TC1=OWF1*OWF2
      TC2=OWF1*(OWF4-OWF8*OWF9)
      TC3=OWF2*(OWF5-OWF7*OWF10)
      TC4=OWF1*OWF11+OWF2*OWF12
      TC5=OUF1*OVF10
      TC6=TC1*OVF1*OVF2
      TC7=OUF1*OWF3-TC2-TC3
      TC8=2.D0*(OWF7+OWF8)
      TC9=OWF1*((WPS(4)-WPS(2))*OWF13+(EX2+EY2)*OVF1*OVF2)+OWF2*((
     1WPS(3)-WPS(1))*OWF13+(EX1+EY1)*OVF1*OVF2)-OUF1*OUF2
      TC10=OWF1*OWF1*OWF6+OWF1*(EX2+EX1+EY2+EY1)*OVF10+OWF2*OWF2*OWF6
     1-OWF2*(EX2+EX1+EY1+EY2)*OVF10+OUF1*OWF6
      TC11=2.D0*(OWF7-OWF8)
C
      F(1)=OWF3*(OUF1+TC1)+3.D0*TC1*OVF5-TC2-TC3
      F(2)=XY*(TC1*OVF8*OVF9-TC4+2.D0*OWF3-TC5)
      F(3)=TC6+TC7
      F(4)=XY*(TC1*OVF6*OVF7-TC4-TC5)
      F(5)=XY*(TC1*(OVF10+3.D0*OVF11)-TC4-TC5)
      F(6)=TC1*OVF3*OVF4+TC7-TC8
      F(7)=TC9-TC1*OVF3*OVF2
      F(8)=TC9-TC1*OVF1*OVF4
      F(9)=-TC8
      F(10)=2.D0*XY*OWF3
      F(11)=-2.D0*XY*OUF2
      F(12)=XY*(TC10-TC1*OVF8*OVF7)
      F(13)=XY*(TC10-TC1*OVF6*OVF9)
      F(14)=TC11
C
      CALL OBSTR(2,1,ME)
      GO TO 1995
C
C
C
C
C
C
C        1-E, VECTOR MESONS  (EXCHANGE DIAGRAM)
C        -------------------------------------
C
C
C
C
C        VECTOR COUPLING (G1*G2)
C
C
C
C
  601 MC=1
C
      FF=EPX1*EPX2*EPY1*EPY2/WPS(1)/WPS(2)/WPS(3)/WPS(4)
      FF=DSQRT(FF)/2.D0
C     FF=0.D0
      VC1=OVF1*OVF2
C
      F(1)=VC1+2.D0*OVF5
      F(2)=-XY*OVF8*OVF9
      F(3)=VC1
      F(4)=XY*OVF6*OVF7
      F(5)=XY*(OVF11+3.D0*OVF10)
      F(6)=OVF3*OVF4
      F(7)=-OVF3*OVF2
      F(8)=-OVF1*OVF4
      F(12)=XY*OVF8*OVF7
      F(13)=-XY*OVF6*OVF9
C
C
      CALL OBSTR(3,1,ME)
C
C
C
C
C        TENSOR COUPLING (F1*F2)
C
C
C
C
      MC=2
C
C
      FF=EPX1*EPX2*EPY1*EPY2/WPS(1)/WPS(2)/WPS(3)/WPS(4)
      FF=DSQRT(FF)/8.D0
C     FF=0.D0
      OWF1=WPS(4)+WPS(1)
      OWF2=WPS(3)+WPS(2)
      OWF3=1.D0+XXYY/EPX1/EPX2/EPY1/EPY2
      OWF4=EX1+EX2+EY1+EY2-OWF1
      OWF5=EX1+EX2+EY1+EY2-OWF2
      OWF6=1.D0/EPX2/EPY1-1.D0/EPX1/EPY2
      OWF7=(EX2-WPS(4))*(EY1-WPS(1))
      OWF8=(EX1-WPS(3))*(EY2-WPS(2))
      OWF9=(EX1+EY2)/EPX2/EPY1+1.D0/EPX2+1.D0/EPY1
      OWF10=(EX2+EY1)/EPX1/EPY2+1.D0/EPX1+1.D0/EPY2
      OWF11=OWF9-OWF4/EPX1/EPY2
      OWF12=OWF10-OWF5/EPX2/EPY1
      OWF13=XX/EPX1/EPX2-YY/EPY1/EPY2
      OUF1=(EX2+EY1)*(EX1+EY2)+XXPYY
      OUF2=1.D0-XXYY/EPX1/EPX2/EPY1/EPY2
C
      TC1=OWF1*OWF2
      TC2=OWF1*(OWF4-OWF8*OWF9)
      TC3=OWF2*(OWF5-OWF7*OWF10)
      TC4=OWF1*OWF11+OWF2*OWF12
      TC5=OUF1*OVF11
      TC6=TC1*OVF1*OVF2
      TC7=OUF1*OWF3-TC2-TC3
      TC8=2.D0*(OWF7+OWF8)
      TC9=OWF1*((WPS(3)-WPS(2))*OWF13+(EX1+EY2)*OVF1*OVF2)+OWF2*((
     1WPS(4)-WPS(1))*OWF13+(EX2+EY1)*OVF1*OVF2)-OUF1*OUF2
      TC10=OWF1*OWF1*OWF6+OWF1*(EX1+EX2+EY2+EY1)*OVF11+OWF2*OWF2*OWF6
     1-OWF2*(EX1+EX2+EY1+EY2)*OVF11+OUF1*OWF6
      TC11=2.D0*(OWF7-OWF8)
C
      F(1)=OWF3*(OUF1+TC1)+3.D0*TC1*OVF5-TC2-TC3
      F(2)=XY*(-TC1*OVF8*OVF9-TC4+2.D0*OWF3-TC5)
      F(3)=TC6+TC7
      F(4)=XY*(TC1*OVF6*OVF7-TC4-TC5)
      F(5)=XY*(TC1*(OVF11+3.D0*OVF10)-TC4-TC5)
      F(6)=TC1*OVF3*OVF4+TC7-TC8
      F(7)=TC9-TC1*OVF3*OVF2
      F(8)=TC9-TC1*OVF1*OVF4
      F(9)=-TC8
      F(10)=2.D0*XY*OWF3
      F(11)=-2.D0*XY*OUF2
      F(12)=XY*(TC10+TC1*OVF8*OVF7)
      F(13)=XY*(TC10-TC1*OVF6*OVF9)
      F(14)=TC11
C
      CALL OBSTR(2,1,ME)
      GO TO 1995
C
C
C
C
C
C
C        1-T , VECTOR-TENSOR COUPLING FOR DIFFERENT MASSES
C        ----------------------------------------------------------
C                       (DIRECT DIAGRAM)
C
C
C        VECTOR-TENSOR COUPLING (G1*F2)
C
C
C
C
  720 MC=1
C
      FF=EPX1*EPX2*EPY1*EPY2/WPS(1)/WPS(2)/WPS(3)/WPS(4)
      FF=DSQRT(FF)/4.D0
C     FF=0.D0
      OWF1=WPS(3)+WPS(1)
      OWF2=WPS(4)+WPS(2)
      OWF3=1.D0+XXYY/EPX1/EPX2/EPY1/EPY2
      OWF4=EX1+EX2+EY1+EY2-OWF1
      OWF5=EX1+EX2+EY1+EY2-OWF2
      OWF6=1.D0/EPX1/EPY1-1.D0/EPX2/EPY2
      OWF7=(EX1-WPS(3))*(EY1-WPS(1))
      OWF8=(EX2-WPS(4))*(EY2-WPS(2))
      OWF9=(EX2+EY2)/EPX1/EPY1+1.D0/EPX1+1.D0/EPY1
      OWF10=(EX1+EY1)/EPX2/EPY2+1.D0/EPX2+1.D0/EPY2
      OWF11=OWF9-OWF4/EPX2/EPY2
      OWF12=OWF10-OWF5/EPX1/EPY1
      OWF13=XX/EPX1/EPX2-YY/EPY1/EPY2
      VTC1=OWF8*OWF9-OWF4
      VTC2=(WPS(4)-WPS(2))*OWF13+(EX2+EY2)*OVF1*OVF2
      VTC3=OWF2*OVF8
      VTC4=OWF2*OVF1
      VTC5=OWF2*OVF6
      VTC6=OWF2*OVF3
      VTC7=OWF1*OWF6+(EX2+EY2+EX1+EY1)*OVF10
C
      F(1)=OWF2*(OWF3+3.D0*OVF5)+VTC1
      F(2)=XY*(VTC3*OVF9-OWF11)
      F(3)=VTC4*OVF2+VTC1
      F(4)=XY*(VTC5*OVF7-OWF11)
      F(5)=XY*(OWF2*(OVF10+3.D0*OVF11)-OWF11)
      F(6)=VTC6*OVF4+VTC1
      F(7)=VTC2-VTC6*OVF2
      F(8)=VTC2-VTC4*OVF4
      F(12)=XY*(VTC7-VTC3*OVF7)
      F(13)=XY*(VTC7-VTC5*OVF9)
C
C
      CALL OBSTR(1,1,ME)
C
C
C
C
C     VECTOR-TENSOR COUPLING (G2*F1)
C
C
C
C
      MC=3
C
      FF=EPX1*EPX2*EPY1*EPY2/WPS(1)/WPS(2)/WPS(3)/WPS(4)
      FF=DSQRT(FF)/4.D0
C     FF=0.D0
      VTD1=OWF7*OWF10-OWF5
      VTD2=OWF1*OVF8
      VTD3=OWF1*OVF1
      VTD4=OWF1*OVF6
      VTD5=OWF1*OVF3
      VTD6=(WPS(3)-WPS(1))*OWF13+(EX1+EY1)*OVF1*OVF2
      VTD7=OWF2*OWF6-(EX2+EY2+EX1+EY1)*OVF10
C
      F(1)=OWF1*(OWF3+3.D0*OVF5)+VTD1
      F(2)=XY*(VTD2*OVF9-OWF12)
      F(3)=VTD3*OVF2+VTD1
      F(4)=XY*(VTD4*OVF7-OWF12)
      F(5)=XY*(OWF1*(OVF10+3.D0*OVF11)-OWF12)
      F(6)=VTD5*OVF4+VTD1
      F(7)=VTD6-VTD5*OVF2
      F(8)=VTD6-VTD3*OVF4
      F(12)=XY*(VTD7-VTD2*OVF7)
      F(13)=XY*(VTD7-VTD4*OVF9)
C
C
      CALL OBSTR(1,1,ME)
      GO TO 1995
C
C
C
C
C
C
C
C        1-TE , VECTOR-TENSOR COUPLING FOR DIFFERENT MASSES
C        ----------------------------------------------------------
C                       (EXCHANGE DIAGRAM)
C
C
C        VECTOR-TENSOR COUPLING (G1*F2)
C
C
C
C
  721 MC=1
C
      FF=EPX1*EPX2*EPY1*EPY2/WPS(1)/WPS(2)/WPS(3)/WPS(4)
      FF=DSQRT(FF)/4.D0
C     FF=0.D0
      OWF1=WPS(4)+WPS(1)
      OWF2=WPS(3)+WPS(2)
      OWF3=1.D0+XXYY/EPX1/EPX2/EPY1/EPY2
      OWF4=EX1+EX2+EY1+EY2-OWF1
      OWF5=EX1+EX2+EY1+EY2-OWF2
      OWF6=1.D0/EPX2/EPY1-1.D0/EPX1/EPY2
      OWF7=(EX2-WPS(4))*(EY1-WPS(1))
      OWF8=(EX1-WPS(3))*(EY2-WPS(2))
      OWF9=(EX1+EY2)/EPX2/EPY1+1.D0/EPX2+1.D0/EPY1
      OWF10=(EX2+EY1)/EPX1/EPY2+1.D0/EPX1+1.D0/EPY2
      OWF11=OWF9-OWF4/EPX1/EPY2
      OWF12=OWF10-OWF5/EPX2/EPY1
      OWF13=XX/EPX2/EPX1-YY/EPY1/EPY2
      VTC1=OWF8*OWF9-OWF4
      VTC2=(WPS(3)-WPS(2))*OWF13+(EX1+EY2)*OVF1*OVF2
      VTC3=-OWF2*OVF8
      VTC4=OWF2*OVF1
      VTC5=OWF2*OVF6
      VTC6=OWF2*OVF3
      VTC7=OWF1*OWF6+(EX2+EY2+EX1+EY1)*OVF11
C
      F(1)=OWF2*(OWF3+3.D0*OVF5)+VTC1
      F(2)=XY*(VTC3*OVF9-OWF11)
      F(3)=VTC4*OVF2+VTC1
      F(4)=XY*(VTC5*OVF7-OWF11)
      F(5)=XY*(OWF2*(OVF11+3.D0*OVF10)-OWF11)
      F(6)=VTC6*OVF4+VTC1
      F(7)=VTC2-VTC6*OVF2
      F(8)=VTC2-VTC4*OVF4
      F(12)=XY*(VTC7-VTC3*OVF7)
      F(13)=XY*(VTC7-VTC5*OVF9)
C
C
      CALL OBSTR(3,1,ME)
C
C
C
C
C     VECTOR-TENSOR COUPLING (G2*F1)
C
C
C
C
      MC=3
C
      FF=EPX1*EPX2*EPY1*EPY2/WPS(1)/WPS(2)/WPS(3)/WPS(4)
      FF=DSQRT(FF)/4.D0
C     FF=0.D0
      VTD1=OWF7*OWF10-OWF5
      VTD2=-OWF1*OVF8
      VTD3=OWF1*OVF1
      VTD4=OWF1*OVF6
      VTD5=OWF1*OVF3
      VTD6=(WPS(4)-WPS(1))*OWF13+(EX2+EY1)*OVF1*OVF2
      VTD7=OWF2*OWF6-(EX2+EY2+EX1+EY1)*OVF11
C
      F(1)=OWF1*(OWF3+3.D0*OVF5)+VTD1
      F(2)=XY*(VTD2*OVF9-OWF12)
      F(3)=VTD3*OVF2+VTD1
      F(4)=XY*(VTD4*OVF7-OWF12)
      F(5)=XY*(OWF1*(OVF11+3.D0*OVF10)-OWF12)
      F(6)=VTD5*OVF4+VTD1
      F(7)=VTD6-VTD5*OVF2
      F(8)=VTD6-VTD3*OVF4
      F(12)=XY*(VTD7-VTD2*OVF7)
      F(13)=XY*(VTD7-VTD4*OVF9)
C
C
      CALL OBSTR(3,1,ME)
      GO TO 1995
C
C
C
C
C        1-C , COULOMB POT.  (PHOTON)
C        ----------------------------
C
C
C
C
 1000 MC=1
C
      FF=2.D0/EESS
      F(1)=1.D0
      F(2)=0.D0
      F(3)=1.D0
      F(4)=0.D0
      F(5)=0.D0
      F(6)=1.D0
      F(7)=-1.D0
      F(8)=-1.D0
      F(12)=0.D0
      F(13)=0.D0
C
      CALL OBSTR (1,1,ME)
      GO TO 1995
C
C
C
C
C
C        2+  , F-MESON
C        -------------
C      MESON GROUP 2+ IS NOT GENERALIZED TO YN-SCATTERING
C      **************************************************
C
C
C        G1**2 COUPLING
C
C
C
C
 1100 MC=1
C
C     FF=-8.D0
C     EE4=EE*4.D0
C     XXYY4=XXYY*4.D0
C     EEP143=EEP1*FD3
C       E1=2.D0*(XXPYY+TD3)
C     F(1)= EEP143 +(3.D0*EE+2.D0)*XXPYY+XXYY4
C     F(2)=(EE4+TD3+XXPYY)*XY
C     F(3)=3.D0*XXYY+EEP1*E1
C     F(4)=(2.D0*XXPYY+3.D0*EEP1-D3)*XY
C     F(5)=(EE4+11.D0*D3+3.D0*XXPYY)*XY
C     F(6)= EEP143 +(EE+2.D0)*XXPYY+XXYY4
C       E2=-EPE*E1
C     F(7)=E2+XX*EX
C     F(8)=E2+YY*EY
C        FACTORS FOR ADDITIONAL TERMS
C     F(9)=XXYY
C     F(10)=EE*XY
C     F(11)=XY
C     F(12)=-EY*XY
C     F(13)=-EX*XY
C
C     CALL OBSTR(3,1,ME)
C     GO TO 1995
C
C
C
C
C        THIS HAS BEEN THE END OF THE CONTRIBUTIONS OF MESONS
C        ----------------------------------------------------
C
C
C
C
C        ERRORS AND WARNINGS
C        -------------------
C
C
C
C
 9000 IF (INDMG(MG)) GO TO 1995
      WRITE (KWRITE,19000) MESONG(MG)
19000 FORMAT(1H0////'0WARNING IN OBNN: MESON-GROUP  ',A4,'  DOES NOT EXI
     1ST IN THIS PROGRAM.'/'0CONTRIBUTION IGNORED. EXECUTION CONTINUED.'
     2////)
      INDMG(MG)=.TRUE.
C
C
C
C
 1995 CONTINUE
C
C
C
C
C        ADD UP CONTRIBUTIONS OF MESONS
C        ------------------------------
C
C
C
C
 2000 DO 2005 IL=IMAN,IMEN
      DO 2005 IV=1,8
      V(IV)=V(IV)+VJ(IV,IL)
      V2(IV)=V2(IV)+VJ2(IV,IL)
 2005 CONTINUE
C
C
C
C
C
      RETURN
      END
C***********************************************************************
C
C***********************************************************************
      SUBROUTINE OBPAR
C
C        OBPAR READS WRITES AND STORES THE PARAMETER FOR ALL OB-SUBROUT.
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     next common wansn't in this new version. maybe it is not needed
      common/print/iprint
      COMMON /CRDWRT/ KREAD,KWRITE
C
      COMMON /CSTATE/ J,HEFORM,SING,TRIP,COUP,ENDEP
      LOGICAL HEFORM,SING,TRIP,COUP,ENDEP
C
C
C        COMMON BLOCK FOR ALL OB-SUBROUTINES
C
      COMMON /COB2/   VJ2(32,90)
      COMMON /COB/    VJ(32,90),C(10,90),FFF,FF,F(52),AA(96),AI(19,5),
     1                WNN(18),X,XX,Y,YY,XY2,XXPYY,EX,EY,EEM12,
     2                EZ1,EZ2,CT(96),WT(96),WPPS(4,18),WEX,
     3                IC(10,90),IFT(18),MINT(18),MAXT(18),NT,
     4                MGE,MGG(11,18),MGGO(11,18),IMA(5,11,18),
     5                IMAA(18),IMEA(18),IME,IM,MC,M,MG,INTER,
     6                INDC(2,90),INDXY
      COMMON /CTYP/   IBARY,INDHAL(18,90),INDONE(18,90)
      COMMON /CCOUL/  VCOUL(8),RCUT,JCOUL
C
C         SPECIFICATIONS FOR THIS COMMON BLOCK
C
      LOGICAL INDC,INDXY
      LOGICAL INDHAL,INDONE
C
C
C        FURTHER SPECIFICATIONS
C
      DIMENSION CC(5)
      CHARACTER*10 NAME
      CHARACTER*70 NNAME
      CHARACTER*27 SNAME
      INTEGER IMGA(18)
      CHARACTER*4 CUT/'CUT '/,CUTG/'CUTG'/,END/'END '/
      CHARACTER*4 MESONG(11)/'0-  ','0-E ','0+E ','0+  ','1-AT',
     1                   '1-  ','1-T ','1-E ','1-TE',
     2                    '1-C ','2+  '/
      LOGICAL INDEX/.FALSE./
      LOGICAL ZEROCP/.TRUE./,INDCUT/.FALSE./
      DATA UF/197.3271D0/
      save
C
C
C
C
10000 FORMAT (A70)
10001 FORMAT (1H1)
10002 FORMAT (1H0/' JP  NAME        A         B       MASS       ISOSPIN
     1    IPROP'/9X,'CUT TYP     C U T - O F F   P A R A M E T E R S')
10003 FORMAT (A10,5F10.4)
10004 FORMAT (1H0,A10,5F10.3)
10005 FORMAT (1H ,A10,F3.1,F11.1,F9.1,F13.1,1X,F8.1)
10006 FORMAT (A10,3I3)
10007 FORMAT (1H0,A10,3I3)
10008 FORMAT (1H ,61(1H-))
10009 FORMAT (1H0,A10,2F10.4,F10.3,2F10.1)
11011 FORMAT (1H1  //' OB11:   ONE-BOSON-EXCHANGE NL-NL INTERACTION')
11012 FORMAT (1H1  //' OB12:   ONE-BOSON-EXCHANGE NL-NS INTERACTION')
11013 FORMAT (1H1  //' OB21:   ONE-BOSON-EXCHANGE NS-NL INTERACTION')
11014 FORMAT (1H1  //' OB22A:  ONE-BOSON-EXCHANGE NS-NS INTERACTION')
11015 FORMAT (1H1  //' OB11DS: ONE-BOSON-EXCHANGE NL-DS INTERACTION')
11016 FORMAT (1H1  //' OB22DS: ONE-BOSON-EXCHANGE NS-DS INTERACTION')
11017 FORMAT (1H1  //' OB11NY*:ONE-BOSON-EXCHANGE NL-NY* INTERACTION')
11018 FORMAT (1H1  //' OB22NY*:ONE-BOSON-EXCHANGE NS-NY* INTERACTION')
11019 FORMAT (1H1  //' OB22B:  ONE-BOSON-EXCHANGE NS-NS INTERACTION')
11020 FORMAT (1H1  //' OB22DL: ONE-BOSON-EXCHANGE NS-DL INTERACTION')
11021 FORMAT (1H1  //' OB12DS: ONE-BOSON-EXCHANGE NL-DS INTERACTION')
11022 FORMAT (1H1  //' OB12DS: ONE-BOSON-EXCHANGE NS-DS INTERACTION')
11023 FORMAT (1H1  //' OB12NY*:ONE-BOSON-EXCHANGE NL-NY* INTERACTION')
11024 FORMAT (1H1  //' OB12NY*:ONE-BOSON-EXCHANGE NS-NY* INTERACTION')
11025 FORMAT (1H1  //' OB21DS: ONE-BOSON-EXCHANGE NS-DS INTERACTION')
11026 FORMAT (1H1  //' OB21DS: ONE-BOSON-EXCHANGE NL-DS INTERACTION')
11027 FORMAT (1H1  //' OB21NY*:ONE-BOSON-EXCHANGE NS-NY* INTERACTION')
11028 FORMAT (1H1  //' OB21NY*:ONE-BOSON-EXCHANGE NL-NY* INTERACTION')
10015 FORMAT ('0INPUT-PARAMETER-SET:'/1H ,20(1H-))
10016 FORMAT (1H0,A70)
10018 FORMAT (A27,I2)
10019 FORMAT (1H0,A27,I2)
C
C
C
C
      IF (INDEX) GO TO 50
      INDEX=.TRUE.
C
      X=-1.D0
      Y=-1.D0
C
C
C
C
C        MAXIMA OF CERTAIN INDICES RELATED TO THE DIMENSION AS FOLLOWS:
C        DIMENSION C(MME,IMEE),IC(MICE,IMEE),INDC(MINDCE,IMEE),
C                  MGG(MGE,18),MGGO(MGE,18),MESONG(MGE),VJ(32,IMEE),
C                  IMA(MEE,MGE,18)
C
      MGE=11
      MEE=5
      MME=10
      MICE=10
      MINDCE=2
      IMB=1
      IME=0
      IMEE=90
      IMEC=0
C        MME ALWAYS GE MICE, MINDCE
C
C        SET ALL MESON-PARAMETERS AND INDICES TO ZERO OR .FALSE.
C
      DO 1 INT=1,18
      IMGA(INT)=0
      DO 1 MGX=1,MGE
      MGG(MGX,INT)=0
    1 MGGO(MGX,INT)=0
C
C
      DO 2 IL=1,IMEE
      DO 3 INT=1,18
      INDHAL(INT,IL)=.FALSE.
    3 INDONE(INT,IL)=.FALSE.
      DO 2 MM=1,MME
      IF (MM.LE.MINDCE) INDC(MM,IL)=.FALSE.
      IF (MM.LE.MICE) IC(MM,IL)=0
    2 C(MM,IL)=0.D0
C
C
C
C        INDC(2,IL) IS RESERVED FOR INFORMATION CONCERNING THE EIKONAL
C        FORM-FACTOR
C
C
C
C
C
C
C        READING AND WRITING OF FIRST 4,5 CARDS
C        --------------------------------------
C        --------------------------------------
C
C
C
C        WRITE HEADLINE AND READ AND WRITE NAME OF PARAMETER SET
C
   50 GO TO (21,22,23,24,25,26,27,28,29,30,31,32,33,
     1       34,35,36,37,38),INTER
   21 WRITE (KWRITE,11011)
      GO TO 55
   22 WRITE (KWRITE,11012)
      GO TO 55
   23 WRITE (KWRITE,11013)
      GO TO 55
   24 WRITE (KWRITE,11014)
      GO TO 55
   25 WRITE (KWRITE,11015)
      GO TO 55
   26 WRITE (KWRITE,11016)
      GO TO 55
   27 WRITE (KWRITE,11017)
      GO TO 55
   28 WRITE (KWRITE,11018)
      GO TO 55
   29 WRITE (KWRITE,11019)
      GO TO 55
   30 WRITE (KWRITE,11020)
      GO TO 55
   31 WRITE (KWRITE,11021)
      GO TO 55
   32 WRITE (KWRITE,11022)
      GO TO 55
   33 WRITE (KWRITE,11023)
      GO TO 55
   34 WRITE (KWRITE,11024)
      GO TO 55
   35 WRITE (KWRITE,11025)
      GO TO 55
   36 WRITE (KWRITE,11026)
      GO TO 55
   37 WRITE (KWRITE,11027)
      GO TO 55
   38 WRITE (KWRITE,11028)
   55 WRITE (KWRITE,10008)
      WRITE (KWRITE,10015)
      READ  (KREAD, 10000) NNAME
      WRITE (KWRITE,10016) NNAME
C
C       READ AND WRITE TYP OF BARYON-BARYON INTERACTION
C             IBARY  = 0   :   NUC.-NUC.
C             IBARY  = 1   :   COUPLED CHANNEL
C
      READ  (KREAD,10018) SNAME,IBARY
      WRITE (KWRITE,10019) SNAME,IBARY
C
C        READ AND WRITE INDEX-PARAMETER CONCERNING THE FACTOR OF THE
C        POTENTIAL
C
      READ  (KREAD, 10006) NAME,IFT(INTER)
      WRITE (KWRITE,10007) NAME,IFT(INTER)
      IFTYP=IFT(INTER)
      IF (IFTYP.LT.0.OR.IFTYP.GT.2)GOTO9003
C
C        READ AND WRITE PARAMETERS FOR NUMERICAL INTEGRATION
C
      READ  (KREAD, 10006) NAME,MINT(INTER),MAXT(INTER)
      WRITE (KWRITE,10007) NAME,MINT(INTER),MAXT(INTER)
C
C        READ AND WRITE MASS OF PARTICLES
C
      READ  (KREAD, 10003) NAME,WN,WPPS(1,INTER),WPPS(2,INTER),
     1                             WPPS(3,INTER),WPPS(4,INTER)
      WRITE (KWRITE,10004) NAME,WN,WPPS(1,INTER),WPPS(2,INTER),
     1                             WPPS(3,INTER),WPPS(4,INTER)
      WNQ=WN*WN
      DWN=1.D0/WN
      DWNQ=DWN*DWN
      WNN(INTER)=WN
      WPPS(1,INTER)=WPPS(1,INTER)*DWN
      WPPS(2,INTER)=WPPS(2,INTER)*DWN
      WPPS(3,INTER)=WPPS(3,INTER)*DWN
      WPPS(4,INTER)=WPPS(4,INTER)*DWN
C
C
C
C        WRITE HEADLINE FOR MESON PARAMETERS
      WRITE (KWRITE,10002)
      WRITE (KWRITE,10008)
C
C
C
C
C        READ, WRITE AND STORE MESON PARAMETERS
C        --------------------------------------
C        --------------------------------------
C
C
C
   61 READ  (KREAD, 10003) NAME,CC
C
C        CHECK IF DATA-CARD JUST READ CONTAINS CUT-OFF PARAMETERS
C
      IF (NAME(1:4).EQ.CUT.OR.NAME(1:4).EQ.CUTG) GO TO 70
C
C        CHECK IF END OF MESONS
C
      IF (NAME(1:4).EQ.END) GO TO 2000
C
C
C
C
C        WRITE MESON-PARAMETERS, WHICH ARE NO CUT-OFF PARAMETERS
C        -------------------------------------------------------
C
C
C
C
      INDCUT=.FALSE.
C
      WRITE (KWRITE,10009) NAME,CC
C
C        CHECK IF COUPLING CONSTANTS ARE ZERO
C
      IF (CC(1).NE.0.D0.OR.CC(2).NE.0.D0) GO TO 62
      ZEROCP=.TRUE.
      GO TO 61
C
   62 ZEROCP=.FALSE.
C
C        FIND OUT NUMBER OF MESON-GROUP MG
C
      DO 63 MG=1,MGE
      IF (NAME(1:4).EQ.MESONG(MG)) GO TO 64
   63 CONTINUE
      GO TO 9000
C
C
C
C
C        STORE MESON PARAMETERS, WHICH ARE NO CUT-OFF PARAMETERS
C        -------------------------------------------------------
C
C
C
C
   64 IME=IME+1
      IF (IME.GT.IMEE) GO TO 9011
      MGG(MG,INTER)=MGG(MG,INTER)+1
      M=MGG(MG,INTER)
      IF (M.GT.MEE) GO TO 9001
      IMA(M,MG,INTER)=IME
      IF (M.NE.1) GO TO 65
      IMGA(INTER)=IMGA(INTER)+1
      MGGO(IMGA(INTER),INTER)=MG
   65 CONTINUE
C        STORE COUPLING CONSTANT G1*G2   OR G1*F2
      C(1,IME)=CC(1)
C        STORE COUPLING CONSTANT G2*F1
      C(3,IME)=CC(2)
C        STORE COUPLING CONSTANT F1*F2
      C(2,IME)=CC(2)
C        STORE MESON MASS SQARE IN UNITS OF NUCLEON MASS SQUARE
      C(4,IME)=CC(3)*CC(3)*DWNQ
      WEX=C(4,IME)
C        TEST ISO-SPIN
      IF (IBARY.NE.0) GO TO 66
      ICC=CC(4)
      IF (ICC.NE.0.AND.ICC.NE.1) GO TO 9004
C         STORE ISOSPIN AS LOGICAL CONSTANT
      IF (ICC.EQ.1) INDC(1,IME)=.TRUE.
      GO TO 69
C
   66 IF (CC(4).NE.0.0.AND.CC(4).NE.0.5.AND.CC(4).NE.1.0) GO TO 9004
      IF (CC(4).EQ.0.5) INDHAL(INTER,IME)=.TRUE.
      IF (CC(4).EQ.1.0) INDONE(INTER,IME)=.TRUE.
C
C        STORE AND TEST IPRSP FOR MESON/DELTA/NUCLEON
C***IC(1,IME)=0 <--> STATIC PROPAGATOR FOR MESONS
C***         =1 <--> START ENERGY Z APPR. BY MASS OF INCOMING PARTICLES
C***         =2 <--> COMPLETE TOPT MESON PROPAGATOR
C***IC(1,IME)=3 <--> NO PROPAGATOR FOR MESONS (CONTACT INTERACTION)
C
   69 ICC=CC(5)
      IC(1,IME)=ICC
      IC(2,IME)=0
      IC(3,IME)=0
      IF (IC(1,IME).LT.0 .OR. IC(1,IME).GT.3) GO TO 9005
      IF (IABS(IC(1,IME)).GE.2.AND.IABS(IC(1,IME)).LE.7) ENDEP=.TRUE.
      IF (IABS(IC(1,IME)).EQ.2) ENDEP=.TRUE.
C
C        INDEX VALUES FOR FURTHER STORING
      MI=4
      MM=5
      GO TO 61
C
C
C
C
C        WRITE CUT-OFF PARAMETERS
C        ------------------------
C
C
C
C
   70 WRITE (KWRITE,10005) NAME,CC
C
C
C     STORE RCUT FOR COULOMB POT.
C
      RCUT=CC(5)
      RCUT=RCUT*WN/UF
C
C        CHECK IF INDIVIDUAL CUT OR GENERAL CUT
C
      IF (NAME(1:4).EQ.CUT) GO TO 73
C        CASE OF GENERAL CUT-OFF
      IF (INDCUT) GO TO 90
      IF (IMEC.GE.IME) GO TO 61
      IMAC=IMEC+1
      IMEC=IME
      IF (IMAC.LT.IMB) IMAC=IMB
      GO TO 90
C        CASE OF INDIVIDUAL CUT-OFF
   73 IMAC=IME
      IMEC=IME
      IF (ZEROCP) GO TO 61
C
C        SAVE PRESENT VALUES OF INDICES
C
   90 INDCUT=.TRUE.
      IF (CC(1).EQ.0.D0) GO TO 61
      MIX=MI
      MMX=MM
C
C        START LOOP OF MESONS, WHICH PRESENT CUT-OFF REFERS TO
C
      DO 1095 IM=IMAC,IMEC
      MI=MIX
      MM=MMX
C
C
C
C
C        STORE CUT-OFF PARAMETERS
C        ------------------------
C
C
C
C
C        STORE TYP OF CUT-OFF
      IC(MI,IM)=CC(1)
      ITYP=IC(MI,IM)
      IF (ITYP.GT.2) GO TO 9002
C        STORE AND TEST TYP OF PROPAGATOR OF CUT-OFF
      IC(MI+1,IM)=CC(2)
      IF (IC(MI+1,IM).NE.0) GO TO 9006
c
c      in the older version, last line wasn there and instead:
c      IF (IC(MI+1,IM).LT.0.OR.IC(MI+1,IM).GT.8) GO TO 9006
c      IF (IC(MI+1,IM).GE.2.AND.IC(MI+1,IM).LE.7) ENDEP=.TRUE.
c
C
C
C        CUT-OFF OF DIPOLE TYPE
C        **********************
C
C
C        STORE AND TEST EXPONENT OF CUT-OFF
      IC(MI+2,IM)=CC(3)
      IF (IC(MI+2,IM).LT.0) GO TO 9009
      IF (IC(MI+2,IM).GT.0) GO TO 101
C        EXPONENT IS ZERO, OMIT CUT-OFF
      IC(MI,IM)=0
      IC(MI+1,IM)=0
      GO TO 1000
C        STORE CUT-OFF MASS FOR DENOMINATOR
  101 C(MM+1,IM)=CC(4)*CC(4)*DWNQ
C        STORE NUMERATOR OF CUT-OFF
      C(MM,IM)=C(MM+1,IM)
      IF (ITYP.EQ.2)     C(MM,IM)=C(MM,IM)-C(4,IM)
      MI=MI+3
      MM=MM+2
C
C
C
C
C        END CUT-OFFS
C        ************
C
C        TEST DIMENSIONS
 1000 IF (MI.GT.MICE.OR.MM-1.GT.MME) GO TO 9010
C
C
 1095 CONTINUE
      GO TO 61
C
C
C
C
C        LAST CARD
C        ---------
C        ---------
C
C
C
C
C        WRITE END MESONS
 2000 IMAA(INTER)=IMB
      IMEA(INTER)=IME
      IMB=IME+1
      WRITE (KWRITE,10004) NAME
      WRITE (KWRITE,10008)
      WRITE (KWRITE,10008)
C
      WRITE(KWRITE,50000)
50000 FORMAT('   IN CASE OF 0+,0-,1-,0-E,1-E,1-AT : A=G1*G2'/,T39,'B=F1*
     1F2'/'   IN CASE OF 1-T,1-TE              : A=G1*F2'/,T39,'B=G2*F1'
     2////)
C
C
C
      RETURN
C
C
C
C        ERRORS
C        ------
C        ------
C
C
C
C
 9000 WRITE (KWRITE,19000) NAME(1:4)
19000 FORMAT (1H0/////'0ERROR IN OBPAR:  MESON-GROUP   ',A4,'   DOES NOT
     1 EXIST IN THIS PROGRAM.'/'0EXECUTION TERMINATED.'////)
      GO TO 9999
C
C
 9001 WRITE (KWRITE,19001)
19001 FORMAT (1H0/////'0ERROR IN OBPAR: TOO MANY MESONS WITHIN A MESON-G
     1ROUP WITH RESPECT TO'/'0THE GIVEN DIMENSIONS. EXECUTION TERMINATED
     2.'////)
      GO TO 9999
C
C
 9002 WRITE (KWRITE,19002) CC(1)
19002 FORMAT (1H0/////'0ERROR IN OBPAR: CUT-OFF TYP',F10.4,'  DOES NOT E
     1XIST IN THIS PROGRAM.'/'0EXECUTION TERMINATED.'////)
      GO TO 9999
C
C
 9003 WRITE (KWRITE,19003) IFTYP
19003 FORMAT (1H0/////'0ERROR IN OBPAR: FACTOR TYP HAS THE NON-PERMISSIB
     1LE VALUE',I4,' .'/'0EXECUTION TERMINATED.'////)
      GO TO 9999
C
C
 9004 WRITE (KWRITE,19004) CC(4)
19004 FORMAT (1H0/////'0ERROR IN OBPAR: ISOSPIN HAS THE NON-PERMISSIBLE
     1VALUE',F10.4,'  .'/'0EXECUTION TERMINATED.'////)
      GO TO 9999
C
C
 9005 WRITE (KWRITE,19005) CC(5)
19005 FORMAT (1H0/////'0ERROR IN OBPAR: IPROP/SPE HAS THE NON-PERMISSIBL
     1E VALUE',F10.4,'  .'/'0EXECUTION TERMINATED.'////)
      GO TO 9999
C
C
 9006 WRITE (KWRITE,19006) CC(2)
19006 FORMAT (1H0/////'0ERROR IN OBPAR: THE INDEX FOR THE PROPAGATOR OF
     1THE CUT-OFF HAS THE'/'0NON-PERMISSIBLE VALUE',F10.4,'  . EXECUTION
     2 TERMINATED.'////)
      GO TO 9999
C
C
 9009 WRITE (KWRITE,19009)
19009 FORMAT (1H0/////'0ERROR IN OBPAR: THE EXPONENT OF THE CUT-OFF IS L
     1ESS THAN ZERO.'/'0EXECUTION TERMINATED.'////)
      GO TO 9999
C
C
 9010 WRITE (KWRITE,19010)
19010 FORMAT (1H0/////'0ERROR IN OBPAR: TOO MANY CUT-OFF PARAMETERS WITH
     1 RESPECT TO THE GIVEN'/'0DIMENSIONS. EXECUTION TERMINATED.'////)
      GO TO 9999
C
C
 9011 WRITE (KWRITE,19011)
19011 FORMAT (1H0/////'0ERROR IN OBPAR:  TOO MANY MESONS WITH RESPECT TO
     1 THE DIMENSIONS GIVEN'/'0TO THIS PROGRAM. EXECUTION TERMINATED.'
     2////)
      GO TO 9999
C
C
 9999 STOP
      END
C***********************************************************************
C
C***********************************************************************
      SUBROUTINE OBSTR (ICASE,MAX,MEX)
C
C        OBSTR COMPUTES THE STRUCTURE OF ONE-BOSON-EXCHANGES
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
C        COMMON BLOCKS
c     again, iprint common is not in the new version
c     common/print/iprint
      COMMON /CRDWRT/ KREAD,KWRITE
C
      COMMON /CSTATE/ J,HEFORM,SING,TRIP,COUP,ENDEP
      LOGICAL HEFORM,SING,TRIP,COUP,ENDEP
C
C
C        COMMON BLOCK FOR ALL OB-SUBROUTINES
C
      COMMON /COB2/   VJ2(32,90)
      COMMON /COB/    VJ(32,90),C(10,90),FFF,FF,F(52),AA(96),AI(19,5),
     1                WNN(18),X,XX,Y,YY,XY2,XXPYY,EX,EY,EEM12,
     2                EZ1,EZ2,CT(96),WT(96),WPPS(4,18),WEX,
     3                IC(10,90),IFT(18),MINT(18),MAXT(18),NT,
     4                MGE,MGG(11,18),MGGO(11,18),IMA(5,11,18),
     5                IMAA(18),IMEA(18),IME,IM,MC,M,MG,INTER,
     6                INDC(2,90),INDXY
      COMMON /CTYP/   IBARY,INDHAL(18,90),INDONE(18,90)
      COMMON /CCOUL/  VCOUL(8),RCUT,JCOUL
C
C         SPECIFICATIONS FOR THIS COMMON BLOCK
C
      LOGICAL INDC,INDXY
      LOGICAL INDHAL,INDONE
C
C     FURTHER SPECIFICATIONS
C
      DIMENSION VV(32),VV2(32)
      DIMENSION TT(2,3)
      DIMENSION TYNO(2,18),TYNH(2,18)
      INTEGER JJ/-1/
      LOGICAL INDEX/.FALSE./
      LOGICAL INDISO
      save
C
C
C      STORE SEVERAL ISOSPIN-FACTORS
C
      IF (INDEX) GO TO 50
      INDEX=.TRUE.
C
      IF (IBARY.NE.0) GO TO 30
C
      TT(1,1)=1.D0
      TT(2,1)=-3.D0
C
      DO 1 II=2,3
      DO 1 I=1,2
    1 TT(I,II)=1.D0
      GO TO 50
C
   30 R2=DSQRT(2.D0)
      R3=DSQRT(3.D0)
      R2D3=DSQRT(2.D0)/R3
      R5D3=DSQRT(5.D0)/R3
      TR2D3=2.D0*R2D3
C
C
C     ISOSPIN FACTORS IN COUPLED CHANNEL FORMALISM
C
C     EXPLANATION:
C     HYP. NUC. :   TYNO (ISOCHANNEL,INTER)    (ISOSPIN=1 - MESONS)
C                       ISOCHANNEL=1  : 1/2
C                       ISOCHANNEL=2  : 3/2
C
C     HYP. NUC. :   TYNH (ISOCHANNEL,INTER)    (ISOSPIN=1/2 - MESONS)
C                       ISOCHANNEL=1  : 1/2
C                       ISOCHANNEL=2  : 3/2
C
C
      TYNO(1,1)=0.D0
C
C FOR TESTING NN MODEL !!
C     TYNO(1,1)=-3.D0
C
      TYNO(2,1)=0.D0
      TYNO(1,2)=R3
      TYNO(2,2)=0.D0
      TYNO(1,3)=R3
      TYNO(2,3)=0.D0
      TYNO(1,4)=-2.D0
      TYNO(2,4)=0.D0
      TYNO(1,5)=-R2
      TYNO(2,5)=0.D0
      TYNO(1,6)=-R2D3
      TYNO(2,6)=-R5D3
      TYNO(1,7)=R3
      TYNO(2,7)=0.D0
      TYNO(1,8)=-2.D0
      TYNO(2,8)=1.D0
      TYNO(1,9)=0.D0
      TYNO(2,9)=1.D0
c*****
c next: new in this version
      TYNO(2,4)=1.D0
c*****
      TYNO(1,10)=0.D0
      TYNO(2,10)=1.D0
C
      TYNO(1,11)=-R2
      TYNO(1,12)=-R2D3
      TYNO(1,13)=R3
      TYNO(1,14)=-2.D0
      TYNO(1,15)=-R2D3
      TYNO(1,16)=-R2
      TYNO(1,17)=-2.D0
      TYNO(1,18)=R3
      DO 31 ITI=11,18
   31 TYNO(2,ITI)=0.D0
C
C
      TYNH(1,1)=1.D0
      TYNH(2,1)=0.D0
      TYNH(1,2)=R3
      TYNH(2,2)=0.D0
      TYNH(1,3)=R3
      TYNH(2,3)=0.D0
      TYNH(1,4)=-1.D0
      TYNH(2,4)=0.D0
      TYNH(1,5)=0.D0
      TYNH(2,5)=0.D0
      TYNH(1,6)=-TR2D3
      TYNH(2,6)=R5D3
      TYNH(1,7)=R3
      TYNH(2,7)=0.D0
      TYNH(1,8)=-1.D0
      TYNH(2,8)=2.D0
      TYNH(1,9)=0.D0
      TYNH(2,9)=2.D0
c****
c  new in this new version
      TYNH(2,4)=2.D0
c****
      TYNH(1,10)=0.D0
      TYNH(2,10)=1.D0
C
      TYNH(1,11)=0.D0
      TYNH(1,12)=-TR2D3
      TYNH(1,13)=R3
      TYNH(1,14)=-1.D0
      TYNH(1,15)=-TR2D3
      TYNH(1,16)=0.D0
      TYNH(1,17)=-1.D0
      TYNH(1,18)=R3
      DO 32 ITI=11,18
   32 TYNH(2,ITI)=0.D0
C
C
C
C
C
C
C
   50 DO 1095 M=MAX,MEX
      IM=IMA(M,MG,INTER)
C
C
      IF (MC.NE.1) GO TO 60
C
C
C
C
C        CALL INTEGRALS
C        --------------
C
C
C
C
      CALL OBAI
C
C
C
C
   60 IF (MC.LT.1) MC=1
C
      IF (C(MC,IM).EQ.0.D0) GO TO 1095
C
C
C
C
      GO TO (100,100,100,100,200,200,200,200,100,200,200,200,200,
     1       200,200,200,200,200),INTER
C
C
C
C
C        NY-NY HELICITY AMPLITUDES
C        -------------------------
C
C
C
C
C        GROUND STRUCTURE (A FACTOR OF 2 IS INCLUDED IN V5 AND V6)
C
C
  100 IVE=8
C
      VV(1)=F(1)*AI(1,M)+F(2)*AI(2,M)
      VV(2)=F(3)*AI(1,M)+F(4)*AI(3,M)
      VV(3)=F(5)*AI(1,M)+F(6)*AI(2,M)
      VV(4)=F(4)*AI(1,M)+F(3)*AI(3,M)
      VV(5)=F(7)*AI(4,M)
      VV(6)=F(8)*AI(4,M)
      VV(7)=F(12)*AI(4,M)
      VV(8)=F(13)*AI(4,M)
C
C
      GO TO (1000,120,130,140),ICASE
C
C
C        ADDITIONAL TERMS  IN CASE OF TENSOR COUPLING
C
C
  120 VV(1)=VV(1)+F(9)*AI(5,M)
      VV(2)=VV(2)+F(10)*AI(2,M)+F(9)*AI(6,M)
      VV(3)=VV(3)+F(10)*AI(5,M)
      VV(4)=VV(4)+F(9)*AI(2,M)+F(10)*AI(6,M)
         E1=F(11)*AI(7,M)
      VV(5)=VV(5)+E1
      VV(6)=VV(6)+E1
      VV(7)=VV(7)+F(14)*AI(7,M)
      VV(8)=VV(8)+F(14)*AI(7,M)
C
C        FOR EXCHANGE DIAGRAM MULTIPLY WITH PHASE FACTORS
C
      IF (MG.EQ.8) GO TO 130
      GO TO 1000
C
C
C        PHASE FACTORS FOR EXCHANGE DIAGRAM
C
C
  130 VV(1)=VV(1)*(-1.D0)**J
      VV(2)=-VV(2)*(-1.D0)**J
      VV(3)=VV(3)*(-1.D0)**J
      VV(4)=VV(4)*(-1.D0)**J
      VV(5)=VV(5)*(-1.D0)**J
      VV(6)=VV(6)*(-1.D0)**J
      VV(7)=VV(7)*(-1.D0)**J
      VV(8)=-VV(8)*(-1.D0)**J
      GO TO 1000
C
C
C        ADDITIONAL TERMS IN CASE OF SIGMA-L IN STATIC LIMIT
C
C
  140 VV(1)=VV(1)+F(6)*AI(5,M)
      VV(2)=VV(2)+F(1)*AI(5,M)+F(9)*AI(6,M)
      VV(3)=VV(3)+F(1)*AI(11,M)
      VV(4)=VV(4)+F(9)*AI(2,M)+F(1)*AI(12,M)
      VV(5)=VV(5)+F(6)*AI(13,M)
      VV(6)=VV(6)+F(6)*AI(13,M)
      GO TO 1000
C
C
C
C
C        NL-SD HELICITY AMPLITUDES
C        -------------------------
C
C
C
C
  200 IVE=16
C
      AI3M1=AI(3,M)-AI(1,M)
      AI6M2=AI(6,M)-AI(2,M)
      AI3P1=AI(3,M)+AI(1,M)
      AI6P2=AI(6,M)+AI(2,M)
C
C
      VV( 1)= F( 1)* AI( 4,M) + F( 2)* AI( 7,M)
      VV( 2)= F( 3)* AI( 1,M) + F( 4)* AI( 2,M) + F( 5)* AI( 5,M)
      VV( 3)= F( 6)* AI( 4,M) + F( 7)* AI( 7,M)
      VV( 4)= F( 8)* AI( 8,M)
      VV( 5)= F( 9)* AI( 8,M)
      VV( 6)= F(10)* AI( 4,M) + F(11)* AI( 7,M)
      VV( 7)= F(12)* AI( 1,M) + F(13)* AI( 2,M) + F(14)* AI( 5,M)
      VV( 8)= F(15)* AI( 4,M) + F(16)* AI( 7,M)
C
      VV( 9)= F(17)* AI3P1    + F(18)* AI6P2
      VV(10)= F(19)* AI( 4,M) + F(20)* AI( 7,M)
      VV(11)= F(21)* AI3M1    + F(22)* AI6M2
      VV(12)= F(23)* AI(10,M)
      VV(13)= F(24)* AI( 9,M)
      VV(14)= F(25)* AI3P1    + F(26)* AI6P2
      VV(15)= F(27)* AI( 4,M) + F(28)* AI( 7,M)
      VV(16)= F(29)* AI3M1    + F(30)* AI6M2
C
      IF (MG.EQ.2.OR.MG.EQ.8) GO TO 250
      GO TO 1000
C
C
C      PHASE FACTORS FOR EXCHANGE GRAPHS
C
  250 VV( 1) = VV( 1)*(-1.D0)**J
      VV( 2) = VV( 2)*(-1.D0)**J
      VV( 3) = VV( 3)*(-1.D0)**J
      VV( 4) = VV( 4)*(-1.D0)**J
      VV( 5) = VV( 5)*(-1.D0)**J
      VV( 6) = VV( 6)*(-1.D0)**J
      VV( 7) = VV( 7)*(-1.D0)**J
      VV( 8) = VV( 8)*(-1.D0)**J
C
      VV9=VV(9)
      VV10=VV(10)
      VV11=VV(11)
      VV12=VV(12)
      VV( 9) =-VV(16)*(-1.D0)**J
      VV(10) =-VV(15)*(-1.D0)**J
      VV(11) =-VV(14)*(-1.D0)**J
      VV(12) =-VV(13)*(-1.D0)**J
      VV(13) =-VV12*(-1.D0)**J
      VV(14) =-VV11*(-1.D0)**J
      VV(15) =-VV10*(-1.D0)**J
      VV(16) =-VV9*(-1.D0)**J
      GO TO 1000
C
C
C
C
C
C
C        NN-DD HELICITY AMPLITUDES
C        -------------------------
C
C
C
C
CCC  300 IVE=32
C
CCC      AI31P=AI( 3,M)+AI( 1,M)
CCC      AI31M=AI( 3,M)-AI( 1,M)
CCC      AI62P=AI( 6,M)+AI( 2,M)
CCC      AI62M=AI( 6,M)-AI( 2,M)
CCC      AIC5P=AI(12,M)+AI( 5,M)
CCC      AIC5M=AI(12,M)-AI( 5,M)
C
C
CCC      VV( 1)= F( 1)*(AI(1,M)+AI(2,M)) + F( 2)*(AI(2,M)+AI(5,M)) +
CCC     1        F( 3)*(AI(5,M)+AI(11,M))
C
CCC      VV( 2)= F( 4)* AI( 4,M) + F( 5)* AI( 7,M) + F( 6)* AI(13,M)
CCC      VV( 3)= F( 7)* AI( 8,M) + F( 8)* AI(14,M)
CCC      VV( 4)= F( 9)* AI(17,M)
CCC      VV( 5)=VV( 2)
CCC      VV( 6)= F(10)* AI( 1,M) + F(11)* AI( 2,M) + F(12)* AI( 5,M) +
CCC     1        F(13)* AI(11,M)
CCC      VV( 7)= F(14)* AI( 4,M) + F(15)* AI( 7,M) + F(16)* AI(13,M)
CCC      VV( 8)=-F(17)* AI( 8,M) + F(18)* AI(14,M)
CCC      VV( 9)=VV( 3)
CCC      VV(10)=VV( 7)
CCC      VV(11)=-F(19)* AI( 1,M) + F(20)* AI( 2,M) - F(21)* AI( 5,M) +
CCC     1        F(22)* AI(11,M)
CCC      VV(12)= F(23)* AI( 4,M) - F(24)* AI( 7,M) + F(25)* AI(13,M)
CCC      VV(13)=VV( 4)
CCC      VV(14)=VV( 8)
CCC      VV(15)=VV(12)
C
CCC      VV(16)=-F(26)*(AI(1,M)-AI(2,M)) + F(27)*(AI(2,M)-AI(5,M)) -
CCC     1        F(28)*(AI(5,M)-AI(11,M))
C
C
CCC      VV(17)= F(29)* AI( 4,M) + F(30)* AI( 7,M) + F(31)* AI(13,M)
CCC      VV(18)= F(32)* AI31P    + F(33)* AI62P    + F(34)* AIC5P
CCC      VV(19)= F(35)* AI( 9,M) + F(36)* AI(15,M)
CCC      VV(20)= F(37)* AI(18,M)
CCC      VV(21)= F(44)* AI31M    - F(45)* AI62M    + F(46)* AIC5M
CCC      VV(22)= F(38)* AI( 4,M) + F(39)* AI( 7,M) + F(40)* AI(13,M)
CCC      VV(23)= F(41)* AI31P    + F(42)* AI62P    + F(43)* AIC5P
CCC      VV(24)=VV(19)
CCC      VV(25)= F(47)* AI(10,M) - F(48)* AI(16,M)
CCC      VV(26)= F(49)* AI31M    - F(50)* AI62M    + F(51)* AIC5M
CCC      VV(27)=VV(22)
CCC      VV(28)=VV(18)
CCC      VV(29)= F(52)* AI(19,M)
CCC      VV(30)=VV(25)
CCC      VV(31)=VV(21)
CCC      VV(32)=VV(17)
C
C
C
C
 1000 IF (INTER.GE.5.AND.INTER.NE.9) GO TO 1040
C
C
C
C
C        SET CERTAIN CASES TO ZERO
C
      IF (J.NE.0) GO TO 1021
      VV(2)=0.D0
      VV(4)=0.D0
      VV(5)=0.D0
      VV(6)=0.D0
      VV(7)=0.D0
      VV(8)=0.D0
C
 1021 IF (.NOT.SING) VV(1)=0.D0
      IF (.NOT.TRIP) VV(2)=0.D0
      IF (COUP) GO TO 1030
      DO 1025 IV=3,6
 1025 VV(IV)=0.D0
C
 1030 IF (HEFORM) GO TO 1040
C
C
C        TRANSFORMATION INTO LSJ-FORMALISM IN CASE OF OBEP
C        (HEFORM = .FALSE. IN MAIN)
      IF (J.EQ.JJ) GO TO 1035
      JJ=J
      AJ=DFLOAT(J)
      AJ1=DFLOAT(J+1)
      D2J1=1.D0/DFLOAT(2*J+1)
      ARJJ1=DSQRT(AJ*AJ1)
C        M.E. ACCORDING TO HOLZENKAMP'S PAPERS
 1035 V3=VV(3)
      V4=VV(4)
      V5=VV(5)
      V6=VV(6)
      V34=-ARJJ1*(V3-V4)
      V56=ARJJ1*(V5+V6)
      VV(3)=D2J1*(AJ1*V3+AJ*V4-V56)
      VV(4)=D2J1*(AJ*V3+AJ1*V4+V56)
      VV(5)=D2J1*(V34-AJ1*V5+AJ*V6)
      VV(6)=D2J1*(V34+AJ*V5-AJ1*V6)
      VV(7)=-VV(7)
      VV(8)=-VV(8)
C
C        POSSIBLE DIFFERENT SIGN DEPENDING ON THE CONVENTION USED
C        (ACCORDING PHYS.REP.149 (1987) 1  EQU.<C.20>)
C     VV(5)=-VV(5)
C     VV(6)=-VV(6)
C
C
C
C
C        MULTIPLY WITH FACTORS
C        ---------------------
C
C
C
C
 1040 IS=MOD(J,2)+1
      IT=MOD(IS,2)+1
      INDISO=INDC(1,IM)
      CMC=C(MC,IM)
      FC=FFF*FF*CMC
      DO 1045 IV=1,IVE
C
C        MULTIPLY WITH COUPLING-CONSTANT AND FACTORS FFF AND FF
C
      VV(IV)=VV(IV)*FC
      VV2(IV)=VV(IV)
C
C        MULTIPLY WITH ISOSPIN FACTOR
C
      IF (IBARY.NE.0) GO TO 1044
      IF (.NOT.INDISO) GO TO 1049
      IF (IV.EQ.2) GO TO 1043
      VV(IV)=VV(IV)*TT(IS,INTER)
      GO TO 1049
 1043 VV(IV)=VV(IV)*TT(IT,INTER)
      GO TO 1049
 1044 IF (.NOT.INDHAL(INTER,IM).AND..NOT.INDONE(INTER,IM)) GO TO 1049
      IF (INDHAL(INTER,IM)) GO TO 1048
      VV(IV)=VV(IV)*TYNO(1,INTER)
      VV2(IV)=VV2(IV)*TYNO(2,INTER)
      GO TO 1049
 1048 VV(IV)=VV(IV)*TYNH(1,INTER)
      VV2(IV)=VV2(IV)*TYNH(2,INTER)
C
C     STORE COULOMB POT.
C
 1049 IF (MG.NE.10) GO TO 1050
      VCOUL(IV)=VV(IV)
      GO TO 1045
C
C     ADD UP IN CASE OF SEVERAL COUPLINGS FOR ONE MESON-EXCHANGE
C     AND STORE
 1050 VJ(IV,IM)=VJ(IV,IM)+VV(IV)
      VJ2(IV,IM)=VJ2(IV,IM)+VV2(IV)
 1045 CONTINUE
C
C
C
C
 1095 CONTINUE
C
C
      RETURN
      END
C***********************************************************************
C
C***********************************************************************
      SUBROUTINE OBAI
C
C        OBAI INTEGRATES OVER THETA
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /CPOT2/  V2(8)
      COMMON /CPOT/   V(8),XMEV,YMEV
      COMMON /CSTATE/ J,HEFORM,SING,TRIP,COUP,ENDEP
      COMMON /CPOTED/ ZREL,NOCED
      COMMON /CRDWRT/ KREAD,KWRITE
C***MASS OF PARTICLES IN STARTING CHANNEL (N+L ODER N+S)****
      COMMON /CWINGO/ WINGO
      LOGICAL HEFORM,SING,TRIP,COUP,ENDEP
      LOGICAL NOCED
C
C        COMMON BLOCK FOR ALL OB-SUBROUTINES
C
      COMMON /COB2/   VJ2(32,90)
      COMMON /COB/    VJ(32,90),C(10,90),FFF,FF,F(52),AA(96),AI(19,5),
     1                WNN(18),X,XX,Y,YY,XY2,XXPYY,EX,EY,EEM12,
     2                EZ1,EZ2,CT(96),WT(96),WPPS(4,18),WEX,
     3                IC(10,90),IFT(18),MINT(18),MAXT(18),NT,
     4                MGE,MGG(11,18),MGGO(11,18),IMA(5,11,18),
     5                IMAA(18),IMEA(18),IME,IM,MC,M,MG,INTER,
     6                INDC(2,90),INDXY
C
C         SPECIFICATIONS FOR THIS COMMON BLOCK
C
      LOGICAL INDC,INDXY
C
C
C        FURTHER SPECIFICATIONS
      DIMENSION GI(7)
      INTEGER IGE/7/
C
      DIMENSION PJ(7,96)
      REAL*4 AEZ,AXY2,AOMQ,AM1,AM2,AM
      INTEGER NNT/-1/,IINTER/-1/,JJ/-1/
      LOGICAL INDJ,INDZ,INDWIN,ITRAFO
      save
C
C
C
C
      IF (INTER.EQ.IINTER) GO TO 60
      IINTER=INTER
      MIN=MINT(INTER)
      MAX=MAXT(INTER)
C
      IGEINT=5
C
C
      WN=WNN(INTER)
      DWN=1.D0/WN
      WNQ=WN*WN
C
C
C
   60 IF (J.EQ.JJ) GO TO 70
      JJ=J
      INDJ=.FALSE.
C
C
      AJ=DFLOAT(J)
      AJ1=DFLOAT(J+1)
      DJ1=1.D0/AJ1
      AJDJ1=AJ*DJ1
      AAJ=DSQRT(AJDJ1)
C
C
      AJ2=DFLOAT(J+2)
      AJ3=DFLOAT(J+3)
      AJM1=DFLOAT(J-1)
      AJM2=DFLOAT(J-2)
      AJM3=DFLOAT(J-3)
C
C
      AJJ1=AJ*AJ1
      AJJ2=AJM1*AJ2
      AJJ3=AJM2*AJ3
      AJJA=AJ*AJM3
      AJJB=AJ*AJM1
C
      AAJJ=0.D0
      IF (J.GT.1)
     1AAJJ=AJ/DSQRT(AJJ1*AJJ2)
C
      AAJ1=AAJJ*AJM1
      AAJ2=AAJJ*AJ1
      AAJ3=AAJJ*2.D0
C
      IF (J.GT.1) GO TO 62
      AAJJ=0.D0
      GO TO 63
   62 AAJJ=1.D0/(AJ1*DSQRT(AJJ2))
C
   63 AAJ4=AAJJ*AJJB
      AAJ5=AAJJ*AJ1*2.D0
      AAJ6=AAJJ*(AJJ1+2.D0)
      AAJ7=AAJJ*AJJ2
C
      IF (J.GT.2) GO TO 64
      AAJJ=0.D0
      GO TO 65
   64 AAJJ=-AJ/DSQRT(AJJ1*AJJ2*AJJ3)
C
   65 AAJ8=AAJJ*(AJJ1+6.D0)
      AAJ9=AAJJ*AJJ2
      AAJ10=AAJJ*(AJJA+2.D0)
      AAJ11=AAJJ*(AJJA-6.D0)
C
      IF (J.GT.2) GO TO 66
      AAJJ=0.D0
      GO TO 67
   66 AAJJ=-1.D0/(AJ1*DSQRT(AJJ2*AJJ3))
C
   67 AAJ12=AAJJ*AJJB*AJM2
      AAJ13=AAJJ*(AJ*AJJB+4.D0*AJ+12.D0)
      AAJ14=AAJJ*(5.D0*AJJ1+6.D0)
      AAJ15=AAJJ*3.D0*AJJ2
      AAJ16=AAJJ*AJJ3*AJM1
      AAJ17=AAJJ*AJ1*AJJ3
      AAJ18=AAJJ*2.D0*AJJ3
C
C
C
   70 C4=C(4,IM)
c     next line new !
      WM=DSQRT(C4)
C
      IPRSP=IC(1,IM)
C
      IF (IPRSP.LT.0 .OR. IPRSP.GT.3) GO TO 9000
C
C        PREPARE STARTING ENERGY
      IF (NOCED) GO TO 74
      NOCED=.TRUE.
      INDZ=.FALSE.
      INDWIN=.FALSE.
      IF (IPRSP.EQ.0.OR.IPRSP.EQ.3) GO TO 74
      IF (IPRSP.EQ.1) GO TO 73
   71 INDZ=.TRUE.
      Z=ZREL
      GO TO 74
C
c  next lines disappeared in this new version
c    72 IF (IPRSP.EQ.8) GO TO 74
c      INDEPZ=.TRUE.
c      EPPQ=PMEV*PMEV
c      EPZ=ZREL
c      GO TO 74

C
   73 INDWIN=.TRUE.
      ZM=WINGO
C
C
c    the following few lines changed in the new version
c
   74 IF (IPRSP.EQ.0.OR.IPRSP.EQ.3) GO TO 90
C
C        PREPARE FOR PROPAGATORS OF NON-COVARIANT PERTURBATION THEORY
      IF (IPRSP.EQ.1) GO TO 75
      IF (.NOT.INDZ) GO TO 71
      EZ=Z
      GO TO 77
   75 IF (.NOT.INDWIN) GO TO 73
      EZ=ZM
   77 QQX=XMEV*XMEV
      QQY=YMEV*YMEV
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C***** INSERT CORRECT MASSES IN THE PROPAGATOR ********
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WP3=WPPS(3,INTER)
      WP4=WPPS(4,INTER)
      IF(MG.EQ.2.OR.MG.EQ.5.OR.MG.EQ.8.OR.MG.EQ.9) GO TO 78
      GO TO 79
   78 WWP3=WP3
      WP3=WP4
      WP4=WWP3
   79 EX1 = WP3*WP3*WNQ + QQX
      EX1 = DSQRT(EX1)
      EX2 = WP4*WP4*WNQ + QQX
      EX2 = DSQRT(EX2)
      EY1 = WPPS(1,INTER)*WPPS(1,INTER)*WNQ + QQY
      EY1 = DSQRT(EY1)
      EY2 = WPPS(2,INTER)*WPPS(2,INTER)*WNQ + QQY
      EY2 = DSQRT(EY2)
      EZ1 = EX1 + EY2 - EZ
      EZ1 = EZ1*DWN
      EZ2 = EX2 + EY1 - EZ
      EZ2 = EZ2*DWN
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  TEST WHETHER PROPAGATOR HAS SINGULARITY
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ET=DMIN1(EZ1,EZ2)
      ETEST=ET+DSQRT(C4+XXPYY-XY2)
      IF (ETEST.LE.0.) GO TO 9100
C
      AEZ=ET
C
C
 
   90 XKQMIN = XXPYY - XY2
      XKQMAX = XXPYY + XY2
C
      ITRAFO=.FALSE.
C
      IF(XKQMAX.GT.10.D0*XKQMIN .AND. XKQMAX.GT.10.D0*WM) GOTO 1000
C
C
C
C
C        FIND OUT APPROPRIATE NUMBER OF GAUSS-POINTS
C        -------------------------------------------
C
C        COMPUTE AM
C
      AXY2=XY2
C
      AOMQ=XXPYY+C4
C
      AM=AXY2/AOMQ
C
      IF (IPRSP.EQ.0.OR.IPRSP.EQ.3) GO TO 93
      AM1=AM
      AM2=AXY2/(AOMQ+AEZ*ABS(AEZ))
      AM=AMAX1(AM1,AM2)
C
C
C        COMPUTE NUMBER OF GAUSSPOINTS (NT)
C
   93 IF (AM.GT.0.99) GO TO 94
C
      IF (AM.GT.0.85) AM=AM**(-ALOG(1.-AM)-0.9)
C
C
      NT=FLOAT(MIN)/(1.-AM)+0.9
C
      IF (NT.GT.MAX) NT=MAX
      GO TO 95
C
C
   94 NT=MAX
C
C
   95 NT=NT+J
C
C        COMPUTE NT, WHICH IS SUITABLE FOR GSET
C
      IF (NT.GT.96) NT=96
C
      IF (NT.EQ.NNT.AND.INDJ) GO TO 100
C
C
C
C        CALL GAUSS-POINTS
C        -----------------
C
C
      CALL GSET (-1.D0,1.D0,NT,CT,WT)
      NNT=NT
C
      GO TO 1500
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IN ORDER TO GET STABLE RESULTS FOR DIAGONAL MATRIX ELEMENTS
C     WITH LARGE OFF-SHELL MOMENTA IN THIS CASE THE GAUSS POINTS
C     ARE NOT CHOSEN FOR COS(THETA) BUT FOR
C
C                      (        1         )
C                  X:= (------------------)   ** (0.5)
C                      ( 1  +  (K/LAM)**2 )
C    WITH K**2 = XMEV**2 + YMEV**2 - 2*XMEV*YMEV*COS(THETA)
C       LAM**2 = MESON MASS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
 1000 ITRAFO=.TRUE.
      NT = MAX+J
      IF (NT.GT.96) NT=96
      NNT = -1
C
      XLAMQ = WM
C
      XMIN = 1.D0/(1.D0+XKQMAX/XLAMQ)**.5D0
      XMAX = 1.D0/(1.D0+XKQMIN/XLAMQ)**.5D0
C
C        CALL GAUSS-POINTS
C        -----------------
      CALL GSET(XMIN,XMAX,NT,CT,WT)
C
C
C     TRANSFORM GAUSS POINTS TO COS(THETA)
C
      DO 1100, I=1,NT
      XKQ = XLAMQ * (1.D0/CT(I)**2-1.D0)
      CT(I) = ( XXPYY - XKQ ) / XY2
 1100 WT(I) = 2.D0*XLAMQ / XY2 * (1.D0+XKQ/XLAMQ)**1.5D0 * WT(I)
C
C
C
C        CALL LEGENDRE-POLYNOMS IF NECESSARY
C        -----------------------------------
C
C
 1500 INDXY=.FALSE.
      INDJ=.TRUE.
      DO 99 I=1,NT
      T=CT(I)
      CALL LEGP (PJ(1,I),PJ(3,I),T,J)
      PJ(2,I)=PJ(1,I)*T
      PJ(4,I)=PJ(2,I)*T
      PJ(6,I)=PJ(4,I)*T
      PJ(5,I)=PJ(3,I)*T
   99 PJ(7,I)=PJ(5,I)*T
C
C
C
C        CALL INTEGRAND
C        --------------
C
C
  100 CALL OBAA
C
C
C
C        PREPARE FOR INTEGRATION
C
C
      DO 2001 IG=1,IGEINT
 2001 GI(IG)=0.D0
C
C
C
C        INTEGRATION-LOOP OF THETA
C        -------------------------
C
      IF(ITRAFO) GO TO 2003
C
      FACT = AA(NT/2)/WT(NT/2)
      DO 2005 I=1,NT
      DO 2005 IG=1,IGEINT
 2005 GI(IG)=GI(IG) + PJ(IG,I)*AA(I) - PJ(IG,I)*FACT*WT(I)
C
C
C
      IF(J.GT.3) GO TO 2010
      IF(J.EQ.0) THEN
      GI(1)=GI(1) + 2.D0*FACT
      GI(3)=0.D0
      GI(4)=GI(4) + 2.D0/3.D0*FACT
      GI(5)=0.D0
      GI(7)=0.D0
      ELSEIF(J.EQ.1) THEN
      GI(2)=GI(2) + 2.D0/3.D0*FACT
      GI(3)=GI(3) + 2.D0*FACT
      GI(6)=GI(6) + 2.D0/5.D0*FACT
      GI(7)=GI(7) + 2.D0/3.D0*FACT
      ELSEIF(J.EQ.2) THEN
      GI(4)=GI(4) + 4.D0/15.D0*FACT
      GI(5)=GI(5) + 2.D0/3.D0*FACT
      ELSEIF(J.EQ.3) THEN
      GI(6)=GI(6) + 4.D0/35.D0*FACT
      GI(7)=GI(7) + 4.D0/15.D0*FACT
      ENDIF
      GO TO 2010
C
C
C
C
 2003 DO 2007 I=1,NT
      DO 2007 IG=1,IGEINT
 2007 GI(IG)=GI(IG)+PJ(IG,I)*AA(I)
C
C
C
      IF (J.NE.0) GO TO 2010
      GI(3)=0.D0
      GI(5)=0.D0
      GI(7)=0.D0
C
C
C
C        COMBINATIONS OF INTEGRALS
C        -------------------------
C
C
 2010 AI(1,M)=GI(1)
C
      AI(2,M)=GI(2)
      AI(3,M)= AJDJ1*GI(2)+DJ1*GI(3)
      GI23M  =GI(2)-GI(3)
      AI(4,M)=AAJ*GI23M
C
C
      AI(5,M)=GI(4)
      AI(6,M)= AJDJ1*GI(4)+DJ1*GI(5)
      GI45M  =GI(4)-GI(5)
      AI(7,M)=AAJ*GI45M
C
C
C*****IF (INTER.NE.5) GO TO 3000
C
C
      AI( 8,M)= AAJ1*GI(4)-AAJ2*GI(1)+AAJ3*GI(5)
      AAI1    = AAJ4*GI(4)+AAJ5*GI(1)-AAJ6*GI(5)
      AAI2    = AAJ7*GI23M
      AI( 9,M)= AAI2+AAI1
      AI(10,M)= AAI2-AAI1
C
C
C*****IF (INTER.NE.6) GO TO 3000
C*******FAKTOREN FUER NN-DD**************
C
C     AI(11,M)=GI(6)
C     AI(12,M)=AJDJ1*GI(6)+DJ1*GI(7)
C     AI(13,M)=AAJ*(GI(6)-GI(7))
C
C     AI(14,M)= AAJ1*GI(6)-AAJ2*GI(2)+AAJ3*GI(7)
C     AAI1    = AAJ4*GI(6)+AAJ5*GI(2)-AAJ6*GI(7)
C     AAI2    = AAJ7*GI45M
C     AI(15,M)= AAI2+AAI1
C     AI(16,M)= AAI2-AAI1
C
C     AI(17,M)= AAJ8*GI(7)-AAJ9*GI(3)-AAJ10*GI(6)+AAJ11*GI(2)
C     AAI1    =-AAJ12*GI(6)+AAJ13*GI(2)-AAJ14*GI(7)+AAJ15*GI(3)
C     AAI2    = AAJ16*GI(4)-AAJ17*GI(1)+AAJ18*GI(5)
C     AI(18,M)= AAI1-AAI2
C     AI(19,M)= AAI1+AAI2
C
C
C
 3000 RETURN
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     E R R O R S
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
 9000 WRITE (KWRITE,19000) IPRSP
19000 FORMAT (1H0/////'0ERROR IN OBAI: IPROP/SPE HAS THE NON-PERMISSIBLE
     1 VALUE',F10.4,'  .'/'0EXECUTION TERMINATED.'////)
      GO TO 9999
 9100 WRITE (6,19100)
19100 FORMAT(1H0//'0SINGULARITY IN THETA INTEGRATION (SUBROUTINE OBAI).'
     1,/'0EXECUTION TERMINATED.')
      GO TO 9999
C
 9999 STOP
      END
C***********************************************************************
C
C***********************************************************************
      SUBROUTINE OBAA
C
C        OBAA COMPUTES THE PROPAGATORS AND THE CUTOFFS OF OB-EXCHANGES
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
c      common/print/iprint
      COMMON /CRDWRT/ KREAD,KWRITE
      COMMON /CSTATE/ J,HEFORM,SING,TRIP,COUP,ENDEP
      COMMON /CPOTED/ ZREL,NOCED
      LOGICAL HEFORM,SING,TRIP,COUP,ENDEP
      LOGICAL NOCED
C
C
C        COMMON BLOCK FOR ALL OB-SUBROUTINES
C
      COMMON /COB2/   VJ2(32,90)
      COMMON /COB/    VJ(32,90),C(10,90),FFF,FF,F(52),AA(96),AI(19,5),
     1                WNN(18),X,XX,Y,YY,XY2,XXPYY,EX,EY,EEM12,
     2                EZ1,EZ2,CT(96),WT(96),WPPS(4,18),WEX,
     3                IC(10,90),IFT(18),MINT(18),MAXT(18),NT,
     4                MGE,MGG(11,18),MGGO(11,18),IMA(5,11,18),
     5                IMAA(18),IMEA(18),IME,IM,MC,M,MG,INTER,
     6                INDC(2,90),INDXY
      COMMON /CCOUL/  VCOUL(8),RCUT,JCOUL
C
C         SPECIFICATIONS FOR THIS COMMON BLOCK
C
      LOGICAL INDC,INDXY
C
C
C
C        FURTHER SPECIFICATIONS
      DIMENSION DELTAQ(96),CUT(96)
C
C
      save
      DATA SSMEV/-1.D0/
C
C
C
C
C
C
C
C
C
C
C
C        DELTA SQUARE
C        ------------
C
C
C
C
   60 IF (INDXY) GO TO 1000
      INDXY=.TRUE.
      DO 65 I=1,NT
      XY2T=XY2*CT(I)
C
C
C
   65 DELTAQ(I)=XY2T-XXPYY
C     ----------------------
C
C
C
C
C
C
C
C
C
C
C        PROPAGATOR
C        ----------
C        ----------
C
C
C
C
 1000 C4=C(4,IM)
      IPRSP=IC(1,IM)
      IF (IPRSP.GE.1.AND.IPRSP.LE.2) GO TO 1050
C
      DO 1005 I=1,NT
C
      AA(I)=WT(I)/(C4-DELTAQ(I))
      IF(IPRSP.EQ.3) AA(I)=WT(I)/C4
C     ----------------------------------
C
C     COULOMB FACTOR
C
      IF (RCUT.EQ.0.D0) GO TO 1005
      DELTA=DSQRT(-DELTAQ(I))
      COC=1.D0-DCOS(DELTA*RCUT)
      AA(I)=AA(I)*COC
 1005 CONTINUE
      GO TO 1500
C
C
C        STARTING ENERGY DEPENDENT PROPAGATOR
C
C
 1050 CONTINUE
      DO 1205 I=1,NT
      OMQ=C4-DELTAQ(I)
      OM=DSQRT(OMQ)
      OMS=OM
C
 1205 AA(I)=WT(I)*(1.D0/(OM*(OMS+EZ1))+1.D0/(OM*(OMS+EZ2)))*0.5D0
C
C
C
 1500 CONTINUE
C
C
C
C
C
C        CUT-OFF
C        --------
C        --------
C
C
C
C
      MI=4
      MM=5
C
C
  999 ITYP=IC(MI,IM)
      IF (ITYP.EQ.0 .AND. MI.GT.4) GO TO 2000
      IF (ITYP.GT.2) GO TO 9000
C
C
C
C
C        CUT-OFF OF DIPOLE TYPE
C        **********************
C
C
C
      C5=C(MM,IM)
      C6=C(MM+1,IM)
      NEXP=IC(MI+2,IM)
C
      DO 105 I=1,NT
C
      AAA=C5/(C6-DELTAQ(I))
C     -------------------------
C
      DO 105 II=1,NEXP
  105 AA(I)=AA(I)*AAA
C
C
      MI=MI+3
      MM=MM+2
      GO TO 999
C
C
C
 2000 RETURN
C
 9000 WRITE (KWRITE,19000) ITYP
19000 FORMAT (1H0/////'0ERROR IN OBAA: CUT-OFF TYP',I2,'  DOES NOT EXIST
     1 IN THIS PROGRAM.'/'0EXECUTION TERMINATED.'////)
 9999 STOP
C
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine tbep(intera)
c
c        iterative box-diagrams with
c        s=1/2  s=3/2 intermediate states
c
c
c
c
c
c
      implicit real*8 (a-h,o-z)
c
      common/print/iprint
      common /crdwrt/ kread,kwrite
c
c        common blocks which contain the arguments and values of the
c        ob-subroutines
c
      common /cpot2/  v2(8)
      common /cpot/   v(8),xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep
      common /cpoted/ zrel,noced
      common /cwingo/ wingo
c
c        specifications for these common blocks
      logical heform,sing,trip,coup,endep
      logical noced
c
c
c        common block for all ob-subroutines
c
      common /cob2/   vj2(32,90)
      common /cob/    vj(32,90),c(10,90),fff,ff,f(52),aa(96),ai(19,5),
     1                wnn(18),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                ez1,ez2,ct(96),wt(96),wpps(4,18),wex,
     3                ic(10,90),ift(18),mint(18),maxt(18),nt,
     4                mge,mgg(11,18),mggo(11,18),ima(5,11,18),
     5                imaa(18),imea(18),ime,im,mc,m,mg,inter,
     6                indc(2,90),indxy
c
c         specifications for this common block
c
      logical indc,indxy
c
      common /ctyp/   ibary,indhal(18,90),indone(18,90)
c
      logical indhal,indone
c
c        further specifications
c
      dimension vv(8),vv2(8)
      dimension vvj(16,90),vvjj(16,90),vvj2(16,90),vvjj2(16,90)
      dimension q(96),qq(96),wq(96),u(96),uu(96)
      dimension idno(10)
      dimension ipropa(18),idan(18),iden(18),idnoo(18,10),njj(18,10)
      dimension cqq(18),wwd(18),wws(18),wwn(18),wwl(18)
  
c
c
      integer nj(20)
      character*10 name
      character*68 nname
      integer jj/-1/,nn/-1/
      integer iitera/-1/
      logical index(18)
      save
      data index/18*.false./
      data pih/1.57079632679490d0/
      data nj/20*0/
c
c
c
c
10000 format (a68)
10001 format (1h0,a68)
10002 format (a10,20i3)
10003 format (1h0,a10,20i3)
10004 format (a10,6f10.4)
10005 format (1h0,a10,6f10.4)
10010 format(1h1//' tbibd:two-boson-exchange - iterative box diagrams'/
     11x,50(1h-)/17x,'with 3/2-1/2 intermediate states'/
     2 17x,34(1h-)/'0input-parameter-set:'/1h ,20(1h-)/1h0)
c
c
c
c
      if (index(intera)) go to 50
      index(intera)=.true.
c
      if(iprint.eq.1) write (kwrite,10010)
      read  (kread ,10000) nname
      if(iprint.eq.1) write (kwrite,10001) nname
c
c        iprop=0 no propagator
c        iprop=1 start energy z approx. by mass of incoming particles
c        iprop=2 relativistic propagator in box-diagram
c        iprop=3 nonrelativ. intermediate energies and
c                start energy z approx. by nonrelativ. energy in
c                lambda-nucleon channel
c
      read  (kread ,10002) name,ipropa(intera)
      if(iprint.eq.1) write (kwrite,10003) name,ipropa(intera)
c
c
      read  (kread ,10002) name,idan(intera),iden(intera)
      if(iprint.eq.1) 
     1    write (kwrite,10003) name,idan(intera),iden(intera)
      read  (kread ,10002) name,(idnoo(intera,ix),ix=1,10)
      if(iprint.eq.1) 
     1    write (kwrite,10003) name,(idnoo(intera,ix),ix=1,10)
      read  (kread ,10002) name,(njj(intera,j1),j1=1,10)
      if(iprint.eq.1) write (kwrite,10003) name,(njj(intera,j1),j1=1,10)
      read  (kread ,10004) name,cqq(intera)
      if(iprint.eq.1) write (kwrite,10005) name,cqq(intera)
      read  (kread ,10004) name,wwd(intera),wws(intera),
     1                          wwn(intera),wwl(intera)
      if(iprint.eq.1) write (kwrite,10005) name,wwd(intera),wws(intera),
     1                          wwn(intera),wwl(intera)
c
c
   50 continue
      iprop=ipropa(intera)
      if (iprop.ge.2) endep=.true.
      ida=idan(intera)
      ide=iden(intera)
      cq=cqq(intera)
      wd=wwd(intera)
      ws=wws(intera)
      wn=wwn(intera)
      wl=wwl(intera)
      do 51 i=1,10
      idno(i)=idnoo(intera,i)
   51 nj(i)=njj(intera,i)
c
c
      wnpl=wn+wl
      wspd=ws+wd
      renl=wn*wl/wnpl
      resd=ws*wd/wspd
      wnq=wn*wn
      wlq=wl*wl
      wsq=ws*ws
      wdq=wd*wd
c
c     cms momentum
      a1=zrel-wlq-wnq
      a2=4.d0*wlq*wnq
      a3=4.d0*zrel*zrel
      q0q=(a1*a1-a2)/a3
c
c
c
c
c
c
c
      if (j.eq.jj.and.iitera.eq.intera) go to 100
      jj=j
      j1=j+1
      aj=dfloat(j)
      aj1=dfloat(j+1)
      d2j1=1.d0/dfloat(2*j+1)
      arjj1=dsqrt(aj*aj1)
      if (nj(j1).eq.0) nj(j1)=nj(j)
      n=nj(j1)
c
c
      if (n.eq.nn.and.iitera.eq.intera) go to 100
      iitera=intera
c
c        get gauss points and weights
c
      call gset (0.d0,1.d0,n,q,wq)
      nn=n
c
c        transform gauss points and weights
c
      do 55 i=1,n
      qaa=pih*q(i)
      q(i)=dtan(qaa)*cq
      qq(i)=q(i)*q(i)
      qaa=1.d0/dcos(qaa)
   55 wq(i)=pih*cq*qaa*qaa*wq(i)
      go to 105
c
c
c
c        compute propagator
c        ------------------
c
c
  100 continue
  105 continue
c
c
c
c        prepare starting energies for sd
c
c
      if (iprop.eq.0) goto 250
      go to (211,212,213),iprop
c
  211 zmevm=wingo
      goto 215
c
  212 zmevm=zrel
      goto 215
c
  213 zmevm=q0q/(2.d0*renl)+wnpl
c
c
c
  215 do 220 i=1,n
      go to (221,221,222),iprop
c
  221 eis=dsqrt(wsq+qq(i))
      eid=dsqrt(wdq+qq(i))
      qmevm=eis+eid
      u(i)=wq(i)*qq(i)/(qmevm-zmevm)
      go to 220
c
  222 qmevm=qq(i)/(2.d0*resd)+wspd
      u(i)=wq(i)*qq(i)/(qmevm-zmevm)
c
  220 continue
      go to 1000
c
  250 do 230 i=1,n
  230 u(i)=wq(i)*qq(i)
c
c
c        this has been the end of the computation of the energy denom.
c
c
c
c
c        call - transition-potentials
c        --------------------------------------
c
c
c
c
 1000 do 1005 id=1,ide
      do 1005 iv=1,16
      vvjj2(iv,id)=0.d0
 1005 vvjj(iv,id)=0.d0
c
c
c        store arguments
c
      xm=xmev
      ym=ymev
c
c
c
      do 2395 i=1,n
      xmev=xm
      ymev=q(i)
c
c     upper transition potential
c     --------------------------
      call obd(intera)
c
      iman=imaa(inter)
      imen=imea(inter)
c
      do 2205 ill=iman,imen
      do 2205 iv=1,16
      vvj2(iv,ill)=vj2(iv,ill)
 2205 vvj(iv,ill)=vj(iv,ill)
c
c
      xmev=ym
      ymev=q(i)
c
      interu=intera
      if (intera.gt.10) interu=intera+1
c
c     lower transition potential
c     --------------------------
      call obd(interu)
c
      imanu=imaa(interu)
      imenu=imea(interu)
c
      id=0
      do 2245 ill=iman,imen
      do 2245 il=imanu,imenu
      id=id+1
      if (id.lt.ida) go to 2245
c
      do 2207 ix=1,10
      if (id.eq.idno(ix)) go to 2245
 2207 continue
c
      if (id.gt.ide) go to 2395
c
c
      do 2206 iv=1,8
      vv2(iv)=0.d0
 2206 vv(iv)=0.d0
c
c
c
      do 2215 ii=1,4
c**** ii1=4+ii
c**** ii3=12+ii
c**** ii6=13-ii
      ii2=8+ii
      ii4=9-ii
      ii5=17-ii
      vv(1)=vv(1)+(vvj(ii,ill)+vvj(ii4,ill))*(vj(ii,il)+vj(ii4,il))
      vv(2)=vv(2)+(vvj(ii2,ill)+vvj(ii5,ill))*(vj(ii2,il)+vj(ii5,il))
      vv(3)=vv(3)+(vvj(ii,ill)-vvj(ii4,ill))*(vj(ii,il)-vj(ii4,il))
      vv(4)=vv(4)+(vvj(ii2,ill)-vvj(ii5,ill))*(vj(ii2,il)-vj(ii5,il))
      vv(5)=vv(5)+(vvj(ii,ill)-vvj(ii4,ill))*(vj(ii2,il)-vj(ii5,il))
      vv(6)=vv(6)+(vvj(ii2,ill)-vvj(ii5,ill))*(vj(ii,il)-vj(ii4,il))
c
c        for comparison with nn-code the following formulars
c        for v5 and v6  ( no antisymmetrisation and no isospin-factors)
c        without the additional factor 2
c
c**** vv(5)=vv(5)+vvj(ii,ill)*vj(ii2,il)
c**** vv(5)=vv(5)+vvj(ii1,ill)*vj(ii3,il)
c**** vv(6)=vv(6)-vvj(ii,ill)*vj(ii5,il)
c**** vv(6)=vv(6)-vvj(ii1,ill)*vj(ii6,il)
      vv(7)=vv(7)+(vvj(ii,ill)+vvj(ii4,ill))*(vj(ii2,il)+vj(ii5,il))
      vv(8)=vv(8)+(vvj(ii2,ill)+vvj(ii5,ill))*(vj(ii,il)+vj(ii4,il))
c
c
      vv2(1)=vv2(1)+(vvj2(ii,ill)+vvj2(ii4,ill))*(vj2(ii,il)+vj2(ii4,il)
     1)
      vv2(2)=vv2(2)+(vvj2(ii2,ill)+vvj2(ii5,ill))*(vj2(ii2,il)+vj2(ii5,i
     1l))
      vv2(3)=vv2(3)+(vvj2(ii,ill)-vvj2(ii4,ill))*(vj2(ii,il)-vj2(ii4,il)
     1)
      vv2(4)=vv2(4)+(vvj2(ii2,ill)-vvj2(ii5,ill))*(vj2(ii2,il)-vj2(ii5,i
     1l))
      vv2(5)=vv2(5)+(vvj2(ii,ill)-vvj2(ii4,ill))*(vj2(ii2,il)-vj2(ii5,il
     1))
      vv2(6)=vv2(6)+(vvj2(ii2,ill)-vvj2(ii5,ill))*(vj2(ii,il)-vj2(ii4,il
     1))
      vv2(7)=vv2(7)+(vvj2(ii,ill)+vvj2(ii4,ill))*(vj2(ii2,il)+vj2(ii5,il
     1))
      vv2(8)=vv2(8)+(vvj2(ii2,ill)+vvj2(ii5,ill))*(vj2(ii,il)+vj2(ii4,il
     1))
c
 2215 continue
c
c
c
c         integrate
c
      do 2235 iv=1,8
      vvjj2(iv,id)=vvjj2(iv,id)-vv2(iv)*u(i)
 2235 vvjj(iv,id)=vvjj(iv,id)-vv(iv)*u(i)
c
c
 2245 continue
c
 2395 continue
c
c     restore arguments
c
      xmev=xm
      ymev=ym
c
c
c     final actions
c
      do 3005 iv=1,8
      v2(iv)=0.d0
 3005 v(iv)=0.d0
c
c
      do 3055 id=1,ide
c
c
      if (j.ne.0) go to 3021
      vvjj(2,id)=0.d0
      vvjj(4,id)=0.d0
      vvjj(5,id)=0.d0
      vvjj(6,id)=0.d0
      vvjj(7,id)=0.d0
      vvjj(8,id)=0.d0
c
c
      vvjj2(2,id)=0.d0
      vvjj2(4,id)=0.d0
      vvjj2(5,id)=0.d0
      vvjj2(6,id)=0.d0
      vvjj2(7,id)=0.d0
      vvjj2(8,id)=0.d0
c
c
 3021 if (.not.sing) vvjj(1,id)=0.d0
      if (.not.trip) vvjj(2,id)=0.d0
c
c
      if (.not.sing) vvjj2(1,id)=0.d0
      if (.not.trip) vvjj2(2,id)=0.d0
      if (coup) go to 3030
      do 3025 iv=3,6
      vvjj2(iv,id)=0.d0
 3025 vvjj(iv,id)=0.d0
c
 3030 if (heform) go to 3040
c
c
c        transformation into lsj-formalism
c        ---------------------------------
c
      v3=vvjj(3,id)
      v4=vvjj(4,id)
      v5=vvjj(5,id)
      v6=vvjj(6,id)
      v34=-arjj1*(v3-v4)
      v56=arjj1*(v5+v6)
      vvjj(3,id)=d2j1*(aj1*v3+aj*v4-v56)
      vvjj(4,id)=d2j1*(aj*v3+aj1*v4+v56)
      vvjj(5,id)=d2j1*(v34-aj1*v5+aj*v6)
      vvjj(6,id)=d2j1*(v34+aj*v5-aj1*v6)
      vvjj(7,id)=-vvjj(7,id)
      vvjj(8,id)=-vvjj(8,id)
c
c
c        possible different sign depending on the convention used
c     vvjj(5,id)=-vvjj(5,id)
c     vvjj(6,id)=-vvjj(6,id)
c
c
      v23=vvjj2(3,id)
      v24=vvjj2(4,id)
      v25=vvjj2(5,id)
      v26=vvjj2(6,id)
      v234=-arjj1*(v23-v24)
      v256=arjj1*(v25+v26)
      vvjj2(3,id)=d2j1*(aj1*v23+aj*v24-v256)
      vvjj2(4,id)=d2j1*(aj*v23+aj1*v24+v256)
      vvjj2(5,id)=d2j1*(v234-aj1*v25+aj*v26)
      vvjj2(6,id)=d2j1*(v234+aj*v25-aj1*v26)
      vvjj2(7,id)=-vvjj2(7,id)
      vvjj2(8,id)=-vvjj2(8,id)
c
c
c        possible different sign depending on the convention used
c     vvjj2(5,id)=-vvjj2(5,id)
c     vvjj2(6,id)=-vvjj2(6,id)
c
c
 3040 do 3045 iv=1,8
c
      vj2(iv,id)=vvjj2(iv,id)
      vj(iv,id)=vvjj(iv,id)
c
      v2(iv)=v2(iv)+vvjj2(iv,id)
 3045 v(iv)=v(iv)+vvjj(iv,id)
 3055 continue
c
c
c
      return
      end
c***********************************************************************
c
c***********************************************************************
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
      subroutine obd(intera)
c
c
c        one-boson-exchange yn-yd interaction;
c        (transition potentials)
c        version which uses numerical integration
c
c
c
      implicit real*8 (a-h,o-z)
c
c
      common/print/iprint
      common /crdwrt/ kread,kwrite
  
  
c        arguments and values of this subroutine
c
      common /cpot2/  v2(8)
      common /cpot/   v(8),xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep
c
c
c        this has been the end of the common-blocks containing
c        the arguments and values of this subroutine in the case of
c        no energy-dependence of the potential;
c        in case of energy-dependence look for the common-block /cpoted/
c        in obai and obaa.
c
c        specifications for these two common blocks
c
      logical heform,sing,trip,coup,endep
c
c
c        common block for all ob-subroutines
c
      common /cob2/   vj2(32,90)
      common /cob/    vj(32,90),c(10,90),fff,ff,f(52),aa(96),ai(19,5),
     1                wnn(18),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                ez1,ez2,ct(96),wt(96),wpps(4,18),wex,
     3                ic(10,90),ift(18),mint(18),maxt(18),nt,
     4                mge,mgg(11,18),mggo(11,18),ima(5,11,18),
     5                imaa(18),imea(18),ime,im,mc,m,mg,inter,
     6                indc(2,90),indxy
c
c         specifications for this common block
c
      logical indc,indxy
c
c
c        further specifications
c
      dimension fa(30)
      dimension wps(4)
      data pi/3.14159265358979d0/
      character*4 mesong(11)/'0-  ','0-e ','0-st','0+  ','1-at',
     1                   '1-  ','1-t ','1-e ','1-te',
     2                    '1+  ','2+  '/
      logical index(18)
      logical indmg(11)
      logical indfa
      save
      data indmg/11*.false./
      data index/18*.false./
c
c
c
c
      inter=intera
c
c
c
c
c        call subroutine obpar once and only once
c
c
      if (index(intera)) go to 50
      index(intera)=.true.
c
c
      call obpar
c
c
   50 continue
      iftgo=ift(inter)+1
      wn=wnn(inter)
      dwn=1.d0/wn
c
      do 10 in=1,4
   10 wps(in)=wpps(in,inter)
c
      dr6=1.d0/dsqrt(6.d0)
      dr12=1.d0/dsqrt(12.d0)
      r2=dsqrt(2.d0)
      dr2=1.d0/r2
      d2r2=1.d0/(2.d0*r2)
c
c
      wy1=wps(1)*wps(1)
      wy2=wps(2)*wps(2)
      wy3=wps(3)*wps(3)
      wy4=wps(4)*wps(4)
      wlps=wps(2)+wps(4)
c
c     w1 is delta mass
c
      w1=wps(1)
      dw1=1.d0/w1
c
      imad=imaa(inter)
      imed=imea(inter)
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
      xa=xmev*dwn
      ya=ymev*dwn
      if (xa.eq.x.and.ya.eq.y) go to 55
      indxy=.false.
      x=xa
      xx=x*x
      y=ya
      yy=y*y
      xy2=x*y*2.d0
      xxpyy=xx+yy
c
c
c
c
   55 xy=xy2*0.5d0
      dxy=1.d0/xy
      ex1=dsqrt(wy3+xx)
      ex2=dsqrt(wy4+xx)
      ey1=dsqrt(wy1+yy)
      ey2=dsqrt(wy2+yy)
cccccccccccc   e=m   cccccccccccccccc
c     ex1=wps(3)
c     ex2=wps(4)
c     ey1=wps(1)
c     ey2=wps(2)
ccccccccccccccccccccccccccccccccccccccccccc
c
c     e0 is delta energy
c
      e0=ey1
      ex=ex1+ex2
      ey=ey1+ey2
      eem12=(ex*ey-1.d0)*2.d0
c
      ee=ex*ey
      eem1=ee-1.d0
      eme=ex-ey
      eep1=ee+1.d0
       epe=ex+ey
      xxyy=xx*yy
c
      epx1=ex1+wps(3)
      epx2=ex2+wps(4)
      epy1=ey1+wps(1)
      epy2=ey2+wps(2)
c
      e2p2=ex2+ey2
      e2p2p1=e2p2+e0
c
c
      xh=x*0.5d0
      xpy=x+y
      xmy=x-y
      xw1=x/epx1
      yw1=y/epy1
      xw2=x/epx2
      yw2=y/epy2
      xyww1=xy/(epx1*epy1)
      xyww2=xy/(epx2*epy2)
      edm1=ey1/wps(1)
      edm1x=edm1*x
      edm1y=edm1*y
c
c
c        the fundamental expressions
c
c
      e1=1.d0-xyww1
      e2=1.d0+xyww1
      e3=yw2-xw2
      e4=yw2+xw2
c
c
c
c        prepare over-all factor
c
c
      go to (70,71,72,73,74),iftgo
c
c        no additional factor
c
   70 fff=fac
      go to 90
c
c        minimal relativity
c
   71 fff=fac/dsqrt(ee)
      go to 90
c
c        factor m/e*m/e
c
   72 ees1=dsqrt(wps(1)/ey1)
      ees2=dsqrt(wps(2)/ey2)
      ees3=dsqrt(wps(3)/ex1)
      ees4=dsqrt(wps(4)/ex2)
      eess=ees1*ees2*ees3*ees4
      fff=fac*eess
      go to 90
c
c        minimal relativity for nd
c
   73 fff=fac/dsqrt(eerel)
      go to 90
c
c        factor m/e*m/e for nd
c
   74 fff=fac/eerel
c
c
c
c
c
c
   90 do 93 iv=1,8
      v2(iv)=0.d0
   93 v(iv)=0.d0
      do 95 il=imad,imed
      do 95 iv=1,32
      vj2(iv,il)=0.d0
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
      if (mg.ne.6.and.mg.ne.8) go to 60
      if (mg.eq.8) go to 61
c
c        abbreviations for vector mesons
      a1p=e2
      a1m=e1
      b2p=e4
      b2m=-e3
      a2p=1.d0+xyww2
      a2m=1.d0-xyww2
      b1p=xw1+yw1
      b1m=xw1-yw1
      go to 62
c
c
c     abbreviations for strange vector mesons
c
   61 a1p=1.d0+xy/(epx2*epy1)
      a1m=1.d0-xy/(epx2*epy1)
      b2p=yw2+xw1
      b2m=-(yw2-xw1)
      a2p=1.d0+xy/(epx1*epy2)
      a2m=1.d0-xy/(epx1*epy2)
      b1p=xw2+yw1
      b1m=xw2-yw1
c
   62 a1pb2p=a1p*b2p
      a1pb2m=a1p*b2m
      a1mb2p=a1m*b2p
      a1mb2m=a1m*b2m
      a2pb1p=a2p*b1p
      a2pb1m=a2p*b1m
      a2mb1p=a2m*b1p
      a2mb1m=a2m*b1m
      a2ma1m=a2m*a1m
      a1pa2m=a1p*a2m
      a1pa2p=a1p*a2p
      a2pa1m=a2p*a1m
c
c
c
c     ground structur of factors for vector mesons
c
c
      fa(1 ) = dr2*(xh*a2pb1p+(xh+y)*a1pb2p)
      fa(2 ) = d2r2*x*(a2pb1p-a1pb2p)
      fa(3 ) = dr6*(edm1y*(a2pb1p+3.d0*a1pb2p)+xh*(a2pb1m-a1mb2p))
     1        +dr6*(xmy*dw1*a1p*(y*a2p+e0*b2p)-xpy*a1mb2p)
      fa(4 ) =-dr6*edm1*(xmy*a2pb1p+(3.d0*x+y)*a1pb2p)
     1        +dr6*(xmy*dw1*a1p*(y*a2p+e0*b2p)+xpy*a1mb2p)
      fa(5 ) = dr6*(edm1x*(a1pb2p-a2pb1p)-xh*(a2pb1m-a1mb2p))
      fa(6 ) =-dr6*(edm1y*(a2pb1m-a1mb2p)-xh*(a2pb1p+3.d0*a1pb2p))
     1        -dr6*xpy*dw1*a1m*(y*a2p+e0*b2p)
      fa(7 ) = dr6*(edm1x*(a2pb1m-a1mb2p)+xh*(a2pb1p-a1pb2p))
      fa(8 ) =-d2r2*x*(a2pb1m-a1mb2p)
      fa(9 ) = d2r2*x*(a1pb2m-a2mb1p)
      fa(10) = dr6*(edm1y*(a1pb2m-a2mb1p)-xh*(a2mb1m+3.d0*a1mb2m))
     1        +dr6*a1p*xmy*dw1*(e0*b2m-y*a2m)
      fa(11) = dr6*x*(edm1*(a2mb1p-a1pb2m)+0.5d0*(a2mb1m-a1mb2m))
      fa(12) =-dr6*(edm1y*(a2mb1m+3.d0*a1mb2m)-xh*(a2mb1p-a1pb2m))
     1        +dr6*(dw1*xpy*a1m*(e0*b2m-y*a2m)-a1pb2m*xmy)
      fa(13) = dr6*edm1*(xpy*a2mb1m+(3.d0*x-y)*a1mb2m)
     1        -dr6*(xpy*dw1*a1m*(e0*b2m-y*a2m)+a1pb2m*xmy)
      fa(14) =-dr6*(edm1x*(a2mb1m-a1mb2m)+xh*(a2mb1p-a1pb2m))
      fa(15) =-d2r2*x*(a2mb1m+3.d0*a1mb2m)
     1        +dr2*xpy*a1mb2m
      fa(16) = d2r2*x*(a2mb1m-a1mb2m)
c
c
      fa(17) = d2r2*x*(a2mb1p+a1pb2m)-dr2*a1pb2m*xmy
      fa(18) =-d2r2*x*(a2mb1p+a1pb2m)
      fa(19) =-dr6*(edm1y*(a2mb1p+a1pb2m)+xh*(a2mb1m+a1mb2m))
     1        -dr6*(xmy*dw1*a1p*(y*a2m-e0*b2m)-xpy*a1mb2m)
      fa(20) = dr6*x*(edm1*(a2mb1p+a1pb2m)+0.5d0*(a2mb1m+a1mb2m))
      fa(21) = dr6*(edm1y*(a2mb1m+a1mb2m)-xh*(a2mb1p+a1pb2m))
     1        +dr6*xpy*dw1*a1m*(y*a2m-e0*b2m)
      fa(22) =-dr6*x*(edm1*(a2mb1m+a1mb2m)+0.5d0*(a2mb1p+a1pb2m))
      fa(23) =-d2r2*x*(a2mb1m+a1mb2m)
      fa(24) = d2r2*x*(a2pb1p+a1pb2p)
      fa(25) = dr6*(edm1y*(a2pb1p+a1pb2p)+xh*(a2pb1m+a1mb2p))
     1        +dr6*xmy*dw1*a1p*(y*a2p+e0*b2p)
      fa(26) =-dr6*x*(edm1*(a2pb1p+a1pb2p)+0.5d0*(a2pb1m+a1mb2p))
      fa(27) =-dr6*(edm1y*(a2pb1m+a1mb2p)-xh*(a2pb1p+a1pb2p))
     1        -dr6*(xpy*dw1*a1m*(y*a2p+e0*b2p)+xmy*a1pb2p)
      fa(28) = dr6*x*(edm1*(a2pb1m+a1mb2p)+0.5d0*(a2pb1p+a1pb2p))
      fa(29) =-d2r2*x*(a2pb1m+a1mb2p)+dr2*xpy*a1mb2p
      fa(30) =-d2r2*x*(a2pb1m+a1mb2p)
c
c
c
   60 continue
      go to (100,200,9000,9000,9000,600,9000,800,9000,9000,9000),mg
c
c
c
c
c        0-  , pseudo-scalar mesons  (direct graph)
c        ------------------------------------------
c
c
c
c
  100 mc=1
c
c
c
c        abbreviations
      e13=e1*e3
      e23=e2*e3
      e14=e1*e4
      e24=e2*e4
c
c        factors for amplitudes 1 to 8
c
      f(1)=-d2r2*x*e13
      f(2)=f(1)
      f(3)=-dr6*e3*(edm1y*e1+xh*e2)
      f(4)=dr6*e13*edm1*xmy
      f(5)=dr6*x*e3*(edm1*e1+0.5d0*e2)
      f(6)=-dr6*e3*(xh*e1-edm1y*e2)
      f(7)=-dr6*e3*(xh*e1+edm1x*e2)
      f(8)=d2r2*x*e23
      f(9)=-d2r2*x*e14
      f(10)=-dr6*e4*(edm1y*e1+xh*e2)
      f(11)=dr6*e4*(edm1x*e1+xh*e2)
      f(12)=-dr6*e4*(edm1y*e2-xh*e1)
      f(13)=dr6*e24*edm1*xpy
      f(14)=-dr6*e4*(edm1x*e2+xh*e1)
      f(15)=-d2r2*x*e24
      f(16)=-f(15)
c
c
c     factors for amplitudes 9 to 16
c
      f(17)=f(9)
      f(18)=-f(9)
      f(19)=-f(10)
      f(20)=-f(11)
      f(21)=f(12)
      f(22)=-f(14)
      f(23)=-f(15)
      f(24)=-f(1)
      f(25)=-f(3)
      f(26)=-f(5)
      f(27)=-f(6)
      f(28)=-f(7)
      f(29)=-f(8)
      f(30)=-f(8)
c
c
c
c
      ffc=epx1*epx2*epy1*epy2/wps(1)/wps(2)/wps(3)/wps(4)
      ffc=dsqrt(ffc)/2.d0
c
c
c        for comparison with nn-code :
c
c**** ffc=ffc*dsqrt(2.d0)
c
c
c
      do 117 mx=1,me
      ff=ffc/dsqrt(c(4,ima(mx,mg,inter)))
  117 call obstr (1,mx,mx)
      if (mg.eq.2) go to 219
      go to 1995
c
c
c
c
c        0-e , pseudo-scalar mesons (exchange graph)
c        -------------------------------------------
c
c
c
c
c        store original values
  200 e1a=e1
      e2a=e2
      e3a=e3
      e4a=e4
c
c
c
c
c        change expressions
      e1=1.d0-xy/(epx2*epy1)
      e2=1.d0+xy/(epx2*epy1)
      e3=yw2-xw1
      e4=yw2+xw1
      go to 100
c
c        restore original values
  219 e1=e1a
      e2=e2a
      e3=e3a
      e4=e4a
      go to 1995
c
c
c
c        1-  , vector mesons   (direct graph)
c        ------------------------------------
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
c
c
      do 601 ifa=1,30
      f(ifa)= fa(ifa)
  601 continue
c
      ffc=epx1*epx2*epy1*epy2/wps(1)/wps(2)/wps(3)/wps(4)
      ffc=dsqrt(ffc)/2.d0
c
c
c
c        for comparison with nn-code :
c
c**** ffc=ffc*dsqrt(2.d0)
c
c
c
      do 617 mx=1,me
      ff=ffc/dsqrt(c(4,ima(mx,mg,inter)))
  617 call obstr (1,mx,mx)
c
c
c
c
c        tensor coupling
c
c
c
c
      mc=3
c
c
c
      f( 1) =-wlps*fa(1)+d2r2*x*a2m*(e2p2*b1p+xpy*a1p)
     1                  -d2r2*x*xmy*a1pa2m
      f( 2) =-wlps*fa(2)+d2r2*x*a2m*(e2p2*b1p+xpy*a1p)
     1                  -d2r2*x*xmy*a1pa2m
      f( 3) =-wlps*fa(3)+dr6*a2m*(edm1y*(e2p2*b1p+xpy*a1p)
     1                            +xh*(e2p2*b1m+xmy*a1m))
     2                  +dr6*(xmy*dw1*y*a1pa2m*e2p2p1-xh*xpy*a2ma1m)
      f( 4) =-wlps*fa(4)-dr6*edm1*xmy*a2m*(e2p2*b1p+xpy*a1p)
     1                  +dr6*xmy*dw1*a1pa2m*(y*e2p2+e0*xpy)
      f( 5) =-wlps*fa(5)-dr6*x*a2m*(edm1*(e2p2*b1p+xpy*a1p)
     1                              +0.5d0*(e2p2*b1m+xmy*a1m))
     2                  +dr6*x*a2m*(edm1*xmy*a1p+0.5d0*xpy*a1m)
      f( 6) =-wlps*fa(6)-dr6*a2m*(edm1y*(e2p2*b1m+xmy*a1m)
     1                            -xh*(e2p2*b1p+xpy*a1p))
     2                  -dr6*a2m*(xpy*y*dw1*a1m*e2p2p1+xh*xmy*a1p)
      f( 7) =-wlps*fa(7)+dr6*x*a2m*(edm1*(e2p2*b1m+xmy*a1m)
     1                              +0.5d0*(e2p2*b1p+xpy*a1p))
     2                  -dr6*x*a2m*(a1m*edm1*xpy+0.5d0*xmy)
      f( 8) =-wlps*fa(8)-d2r2*x*a2m*(e2p2*b1m+xmy*a1m)
     1                  +d2r2*x*xpy*a2ma1m
      f( 9) =-wlps*fa(9)-d2r2*x*a2p*(e2p2*b1p+xpy*a1p)
     1                  +d2r2*x*xmy*a1pa2p
      f(10) =-wlps*fa(10)-dr6*a2p*(edm1y*(e2p2*b1p+xpy*a1p)
     1                            +xh*(e2p2*b1m+xmy*a1m))
     2                  -dr6*a2p*(xmy*y*dw1*a1p*e2p2p1-xh*xpy*a1m)
      f(11) =-wlps*fa(11)+dr6*x*a2p*(edm1*(e2p2*b1p+xpy*a1p)
     1                               +0.5d0*(e2p2*b1m+xmy*a1m))
     2                   -dr6*x*a2p*(edm1*xmy*a1p+0.5d0*xpy*a1m)
      f(12) =-wlps*fa(12)-dr6*a2p*(edm1y*(e2p2*b1m+xmy*a1m)
     1                             -xh*(e2p2*b1p+xpy*a1p))
     2                   -dr6*a2p*(y*xpy*dw1*a1m*e2p2p1+xh*xmy*a1p)
      f(13) =-wlps*fa(13)+dr6*edm1*xpy*a2p*(e2p2*b1m+xmy*a1m)
     1                   +dr6*dw1*xpy*a2pa1m*(y*e2p2-e0*xmy)
      f(14) =-wlps*fa(14)-dr6*x*a2p*(edm1*(e2p2*b1m+xmy*a1m)
     1                               +0.5d0*(e2p2*b1p+xpy*a1p))
     2                   +dr6*x*a2p*(edm1*xpy*a1m+0.5d0*xmy*a1p)
      f(15) =-wlps*fa(15)-d2r2*x*a2p*(e2p2*b1m+xmy*a1m)
     1                   +d2r2*x*xpy*a2pa1m
      f(16) =-wlps*fa(16)+d2r2*x*a2p*(e2p2*b1m+xmy*a1m)
     1                   -d2r2*x*xpy*a2pa1m
c
c
      f(17) =-wlps*fa(17)+d2r2*x*a2p*(e2p2*b1p+xpy*a1p)
     1                   -d2r2*x*a1pa2p*xmy
      f(18) =-wlps*fa(18)-d2r2*x*a2p*(e2p2*b1p+xpy*a1p)
     1                   +d2r2*x*a1pa2p*xmy
      f(19) =-wlps*fa(19)-dr6*(edm1y*a2p*(e2p2*b1p+xpy*a1p)
     1                          +xh*a2p*(e2p2*b1m+xmy*a1m))
     2                    -dr6*a2p*(xmy*y*dw1*a1p*e2p2p1-xh*xpy*a1m)
      f(20) =-wlps*fa(20)+dr6*x*a2p*(edm1*(e2p2*b1p+xpy*a1p)
     1                               +0.5d0*(e2p2*b1m+xmy*a1m))
     2                   -dr6*x*a2p*(xmy*edm1*a1p+0.5d0*xpy*a1m)
      f(21) =-wlps*fa(21)+dr6*a2p*(edm1y*(e2p2*b1m+xmy*a1m)
     1                             -xh*(e2p2*b1p+xpy*a1p))
     2                   +dr6*a2p*(xpy*y*dw1*a1m*e2p2p1+xh*xmy*a1p)
      f(22) =-wlps*fa(22)-dr6*x*a2p*(edm1*(e2p2*b1m+xmy*a1m)
     1                               +0.5d0*(e2p2*b1p+xpy*a1p))
     2                   +dr6*x*a2p*(xpy*edm1*a1m+0.5d0*xmy*a1p)
      f(23) =-wlps*fa(23)-d2r2*x*a2p*(e2p2*b1m+xmy*a1m)
     1                   +d2r2*x*xpy*a2pa1m
      f(24) =-wlps*fa(24)+d2r2*x*a2m*(e2p2*b1p+xpy*a1p)
     1                   -d2r2*x*xmy*a1pa2m
      f(25) =-wlps*fa(25)+dr6*a2m*(edm1y*(e2p2*b1p+xpy*a1p)
     1                             +xh*(e2p2*b1m+xmy*a1m))
     2                   +dr6*a2m*(xmy*y*dw1*a1p*e2p2p1-xh*xpy*a1m)
      f(26) =-wlps*fa(26)-dr6*x*a2m*(edm1*(e2p2*b1p+xpy*a1p)
     1                               +0.5d0*(e2p2*b1m+xmy*a1m))
     2                   +dr6*x*a2m*(edm1*xmy*a1p+0.5d0*xpy*a1m)
      f(27) =-wlps*fa(27)-dr6*a2m*(edm1y*(e2p2*b1m+xmy*a1m)
     1                             -xh*(e2p2*b1p+xpy*a1p))
     2                   -dr6*a2m*(xpy*y*dw1*a1m*e2p2p1+xh*xmy*a1p)
      f(28) =-wlps*fa(28)+dr6*x*a2m*(edm1*(e2p2*b1m+xmy*a1m)
     1                               +0.5d0*(e2p2*b1p+xpy*a1p))
     2                   -dr6*x*a2m*(edm1*xpy*a1m+0.5d0*xmy*a1p)
      f(29) =-wlps*fa(29)-d2r2*x*a2m*(e2p2*b1m+xmy*a1m)
     1                   +d2r2*x*xpy*a2ma1m
      f(30) =-wlps*fa(30)-d2r2*x*a2m*(e2p2*b1m+xmy*a1m)
     1                   +d2r2*x*xpy*a2ma1m
c
      ffc=epx1*epx2*epy1*epy2/wps(1)/wps(2)/wps(3)/wps(4)
      ffc=-dsqrt(ffc)/4.d0
c
c
c
c        for comparison with nn-code :
c
c**** ffc=ffc*dsqrt(2.d0)
c
c
c
      do 627 mx=1,me
      ff=ffc/dsqrt(c(4,ima(mx,mg,inter)))
  627 call obstr (1,mx,mx)
      if (mg.eq.8) go to 819
      go to 1995
c
c
c
c
c        1-e ,vector mesons (exchange graph)
c        -----------------------------------
c
c
c     store original values
c
  800 wlpsa=wlps
      e2p2a=e2p2
      e2p2a1=e2p2p1
c
c
c     change expressions
c
      wlps=wps(2)+wps(3)
      e2p2=ex1+ey2
      e2p2p1=e2p2+e0
      go to 600
c
c     restore original values
c
  819 wlps=wlpsa
      e2p2=e2p2a
      e2p2p1=e2p2a1
      go to 1995
c
c        this has been the end of the contributions of mesons
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
      if(iprint.eq.1) write (kwrite,19000) mesong(mg)
19000 format(1h0////'0warning in obd2: meson-group  ',a4,'  does not exi
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
 2000 return
      end
      subroutine legp (pj,pjm1,x,j)
c
c        subroutine legp   computes the legendre polynoms
c
      real*8 pj,pjm1,x,a,b
      save
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
