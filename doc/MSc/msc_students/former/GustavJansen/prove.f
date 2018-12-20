      implicit real*8(a-h,o-z)
      dimension fkmev(50),help(50),au(50)
      dimension xk(50),xkwgt(50)

      common/info/ivvsm,iaasm,iap,iasp,ibp,ibsp,iwap,iwaap,
     1 isp,iwkp,iwkkp
      COMMON /CRDWRT/ KREAD,KWRITE
      common /hpot/VV11(70,70,8),VV12(70,70,8),
     1           VV21(70,70,8),VV22(70,70,8),
     2           VV(70,70,8)

      common/effmass/wneff,w1eff,w2eff,iter

      hbc=197.327d0
      hbc3=hbc**3.d0

c... info for the Juelich potential

      kread=32
      kwrite=16
      ivvsm=81
      iaasm=25920
      iap=162
      iasp=324
      ibp=162
      ibsp=324
      iwap=26568
      iwaap=105624
      isp=80
      iwkp=162
      iwkkp=324

      open(unit=kread,file='mass_2.data')

C .......

c ...  Logarithm Map

      xkmax=20.d0
      nk=24
      CALL GSET(0.D0,1.D0,nk,HELP,AU)
      FAC=(xkmax-0.D0)/DLOG((1.D0+HELP(NK))/(1.D0-HELP(NK)))
      DO 205 I=1,NK
      XKWGT(I)=AU(I)*FAC*2.D0/(1.D0-HELP(I)**2.D0)
 205  xk(I)=0.D0+FAC*DLOG((1.D0+HELP(I))/(1.D0-HELP(I)))

      do i=1,nk
      fkmev(i)=xk(i)*hbc ! We put all momenta in MeV
      enddo

C .....

      jc=0

      do iter=1,1

      if(iter.eq.1) then
      wneff=938.919d0
      w1eff=1115.6d0
      w2eff=1193.1d0
      elseif(iter.eq.2) then
      wneff=938.919d0/2.d0
      w1eff=1115.6d0/2.d0
      w2eff=1193.1d0/2.d0
      endif      

C ....

C ....

      call juelp(jc,fkmev,nk,0)
      rewind(kread)

      do ik=1,nk
      write(6,'(6d15.5)')xk(ik)
     &,vv11(ik,ik,1)*hbc3
     &,vv12(ik,ik,1)*hbc3
     &,vv21(ik,ik,1)*hbc3
     &,vv22(ik,ik,1)*hbc3
     &,vv(ik,ik,1)*hbc3
      enddo

      enddo

      stop 
      end
