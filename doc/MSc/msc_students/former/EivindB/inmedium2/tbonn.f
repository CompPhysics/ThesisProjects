C*************************************************************************      
C  Calculate the "phase shifts", "mixing factor" and "inelasticitys" for 
C  nn-nn scattering using BONN B potencial.By finding the "T" & "S"-matrix
C      
C  Uses lsj-formalism
C  J is choosen from main program (as a parameter)
C  "Blatt-Biedenharn" or "Stapp/bar"-formalsm is choosen from main prog     
C  Natural units 
C  Uses "constants.h" 
C  N is the number of meshpoints (from Gauss-Legendre)
C q is the momentum in MC-system
      
C DLINCG()="IMSL"-func: Compute the inverse of a complex general matrix      
C**************************************************************************

C==>  solved in momentum space      
C==>  c = chbar = 1      
      
C--   [phaseshift] = degrees
C-- Re[mixing factor] = degrees
C-- Im[mixing factor] = -
C--   [inelasticity] = -
C--   [mass] = MeV
C--   [chbar] = -
C--   [q] = MeV
C--   [V(q1,q2)] = 1/(MeV)**2      
C--   [qEx] = MeV
C--   [T-matix] = 1/(MeV)**2
C--   [S-matrix] = -

!      include "nna13in.h"
      include "nna13.h"
      include "bonn.h"
      
      program main
      implicit real*8 (a-h,o-z)
      include "constants.h" 
      
!      common /cpts/   qx(97) !"N+1 mesh points (used in nna13 potential)
!      common /cpoted/ q0qmev ! =q(N+1)**2 (used in the nna13 potential)
     
      parameter(NExper=17) !==> length of qEx-arrayen 
      parameter(Jmax=20)  !==> calculate all partiale waves with 0<=J<=Jmax 
      parameter(wp=938.27231d0,wn=939.56563d0) !from Phys.Rev C 63 024001 
      integer ipot
      real*8 q(N+1),w(N),qEx(NExper),qF
      complex*16 bepj  
      logical blatt,part1,part2,propag,kowalski,inmedium,degfixed
      logical vinmedium,scattering
      common/cinmedium/starM,qF,inmedium,vinmedium
      common/cmass/redM
      common/cmassp/w1,w2
      common/cprop/propag
      common/cipot/ipot
      common/clseq/kowalski
      common/boundfix/degfixed
      common/formalism/blatt     
      
      
      
C========== Calculate phase shifts at theese lab energies ===============

C-----Low energies (NExper=15)
!      data qEx/1.d-2,2.d-2,3.d-2,1.d-1,1.d0,5.d0,10.d0,25.d0,50.d0,
!     &         100.d0,150.d0,200.d0,250.d0,300.d0,350.d0/

C-----Hi energies (NExper=15)
!      data qEx/2.5d1,5.d1,1.d2,2.d2,3.d2, 4.d2, 5.d2, 6.d2, 
!     &         7.d2, 8.d2, 9.d2, 1.d3,1.1d3,1.2d3,1.4d3/

C-----In-medium energies (NExper=15)
!      data qEx/1.d-2, 1.d-1, 1.d0,1.d1,1.d2,1.5d2,
!     &          2.d2,3.d2,4.d2,5.d2,6.d2, 7.d2, 8.d2, 9.d2,1.d3/
      
C-----In-medium energies (NExper=9)
!      data qEx/2.d2,3.d2,4.d2,5.d2,6.d2,7.d2,8.d2,9.d2,1.d3/
      
C-----Hi-energy In medium (NExper=17)
      data qEx/2.d2, 2.5d2, 3.d2, 3.5d2, 4.d2, 4.5d2, 5.d2, 5.5d2,
     &   6.d2, 6.5d2, 7.d2, 7.5d2, 8.d2, 8.5d2, 9.d2, 9.5d2, 1.d3/
C========================================================================      
      
      
10000 format(f7.2, 8f12.6,2x, 2f10.6 )
10001 format(f7.2, 4f12.6,2x )       
10002 format(//9X,'******** 1S0 ********',4X,'******* 3P0 ********')
10003 format(//10X,'******* 1P1 ********',4X,
     & '******** 3P1 *******', 4X,'******* 3S1 ********',4X,
     & '******* 3D1 ********',4X,'******* e1 *******')
10004 format(//10X,'******* 1D2 ********',4X,
     & '******* 3D2 ********' 4X,'******** 3P2 *******',4X,
     & '******** 3F2 *******',4X,'******* e2 *******')
10005 format(//10X,'******* 1F3 ********',4X,
     & '******* 3F3 ********' 4X,'******** 3D3 *******',4X,
     & '******** 3G3 *******',4X,'******* e3 *******')
10006 format(//9X,' ******* 1',a,i2,' *******',3X,' ******* 3',a,i2,
     & ' *******'3X,' ******* 3',a,i2,' *******',3X,' ******* 3',a,i2,
     & ' *******',4X,'******* e',i2,' ******')  
10007 format('T-LAB',4x,'PHA SHIFT',4x,'INELASTI',4x,
     & 'PHA SHIFT',3x,'INELASTI')      
10008 format('T-LAB',5x,'PHA SHIFT',3x,'INELASTI',4x,
     & 'PHA SHIFT',3x,'INELASTI',4x,'PHA SHIFT',3x,'INELASTI',
     & 4x,'PHA SHIFT',3x,'INELASTI')
10009 format(//' '/,i3,' mesh points'/,a)
10010 format( '2 * redused mass = ',f12.6,' MeV')      
10011 format(f9.2,5d14.6,/9x,5d14.6) 
10012 format('J=',i2,/'ELAB(MEV)',8x,'SINGLET',8x,
     &        'TRIPLET',8x,'L=J-1',8x,'L=J+1',8x,'EPS')
10013 format(///'Negative values are poosible bound state energies:') 
      
C===================================================================



      
C-----If scattering=.true. then calculate scattering paremeres
C-----If scattering=.false. then calculate the binding energy      

      scattering=.true.

      

      
C-----if blatt=.true. then "Blatt-Biedenharn" -formalism
C-----if blatt=.false. then "Stapp/bar" -formalism
       
      blatt=.false.
      
      if(blatt)then
              write(*,10009)N,'"Blatt"-formalism'
      else
              write(*,10009)N,'"Stapp/bar"-formalism'
      end if 
      
      

C-----if "ipot=1" then bonn b potential
C-----if "ipot=2" then nna13  potential 
C-----if "ipot=3" then box  potential (only valid for the 1S0-state)
C-----if "ipot=4" then box  potential (only valid for the 1S0-state)     
C-----NOT IMPLEMENTED

      ipot=2
      
      if(ipot.eq.1) then
         write(*,*)'Uses "bonn" potential (good up til about 350 MeV.'
      else if(ipot.eq.2) then
         write(*,*)'Uses "nna13" potential (is good up til 1GeV.)'
         write(*,*)'(nna13 works only for np-np scattering)'
         if (N.gt.97) then           
           write(*,*)'Number of mesh points have to be less or equal 97'
           go to 99999
         end if
      else if(ipot.eq.3) then
         write(*,*)'Uses box potential (Only valid for 1S0-states)'
      else if(ipot.eq.4) then     
         write(*,*)'Uses Yukawa potentials (Only valid for 1S0-states)'
      else   
         write(*,*)'Wrong potential in main-program ("ipot"-value)'
      end if

      
C-----if .true. then proton
C-----if .false. then neutron     
      
      part1=.true.  !Particle moving in lab system
      part2=.false.  !Particle not moving in lab system
      
      if(part1 .eq. .true.) w1=wp
      if(part1 .eq. .false.)w1=wn
      if(part2 .eq. .true.) w2=wp
      if(part2 .eq. .false.)w2=wn

      if(part1.and.part2 .eq. .true.)then
              write(*,*)'pp-pp sactering' 
      else if(part1.or.part2 .eq. .true.)then
              write(*,*)'np-np sactering' 
      else
              write(*,*)'nn-nn sactering' 
      end if

Cdddd remove dddddddddd      Uses charge independence w1=w2=938.926MeV
      w1=938.926d0
      w2=938.926d0
Cdddddddddddddddddddddd
C-----redM=two times redused mass
      redM=2.d0*(w1*w2/(w1+w2))
      write(*,10010)redM
      

     
      
      
C-----if .true. then "minimal relativety"/ none relativistic propagator (BbS)
C-----if .false. then "Thompson" relativistic propagator (gTh)
      
      propag=.true.

      if(propag .eq. .true.) then 
              write(*,*)'"Minimal Relativety". (has special  dbonn.d)'
      else
              write(*,*)'"Thompson". (has its own special dbonn.d)'
      end if



C-----if .true. solve LS by K.L.Kowalski's metod(phys. rev.Lett. 15 1965)
C-----if .false. solve LS by normal metod (Brown and Jackton eq V(85.1)      
      
      kowalski=.false.
     
      if(kowalski .eq. .true.) then
             write(*,*)'Solving LS eq by using K.L.Kowalski method'
      else
             write(*,*)'Solving LS by Brown & Jackton eq V(85)'
      end if






      
C-----If .true. then calculate the NN phase shift in medium.
C-----If .false. then calculate the NN phase shift in free space 
      
      inmedium=.false.
      vinmedium=.false. !==> change nna13.h to nna13in.h
      
      if(inmedium .eq. .true.) then 
           starM=688.d0 !MeV 
           qF=1.35d0*chbar !MeV 
           write(*,*)'NN pasheshift in medium, uses Mstar=',starM,'MeV'
           write(*,*)'with Fermi mometum kF=',qF/chbar,'fm-1'
      else 
           starM=redM
           qF=0.d0 !MeV
           write(*,*)'Free NN pasheshifts' 
      end if



      
      
C-----makes special file later used to calculate cross sections
      write(*,*)'output for obs.f written in <phasedat.d>'
      open (unit=8,file='phasedat.d')
      kwrite=8
      write(kwrite,*)'potential nr in my code: ipot =',ipot
      write(kwrite,10012)J



      
C==================================================================     
C-----Calculates the bound energy for sing,trip and coubled states.
C-----The ii and JJ combination tells us which partial wave state 
C-----  we are calculating the bounding energy of.  
      if(scattering.eq..false.)then 
           ii=5 ! ii described in the potential
           J=1  !
           write(*,10013)
           call boundE(ii,J) 
           go to 99999  ! 99999="terminate prog"
      end if
C==================================================================
      
      
     
      
C-----finds meshpoits and weights. half the points will be in (0,const) 
      const=1000.d0
      
C-----------free NN-scattering (0-inf)      
            call gauleg(N,2,0.d0,1000.d0,q,w) !Gauss-Legendre method

C-----------In-medium scattering uses different integral (qF-inf)
            if(inmedium .eq. .true.) then
                    do 23,i=1,N
                    q(i)=q(i)+qF
   23               continue  
            end if


            
            
      if(meshBAD(q,qEx,NExper)) go to 99999 ! 99999="terminate prog"     
      
      
           
      write(*,10002)
      write(*,10007)
      
      
      
C--J=0      
      J=0

C--do-lop over diferent energies in lab system (Tlab=qEq(i) )      
      do 10, i=1,NExper
      
      q(N+1)=dsqrt( (2.d0*w2**2*w1*qEx(i)+w2**2*qEx(i)**2)
     &            /(w1**2+w2**2+2.d0*w2*(qEx(i)+w1)) )

C-----nna13 potential wants all the q in a different in a spesial way
      if(ipot.eq.2)  call nna13Spesial(q) 
      
      call GetPhaS(q,w,J,1,pha1,el1)  !in the bonn potesial==> "V(1)"       
      call GetPhaS(q,w,J,3,pha2,el2)  !in the bonn potesial==> "V(3)"                              
      
      write(*,10001)qEx(i),pha1,el1,pha2,el2
      write(kwrite,10011)qEx(i),pha1,0.d0,0.d0,pha2,0.d0,
     &                          el1,1.d0,1.d0,el2,0.d0 
      
   10 continue
      write(kwrite,*)' '
      
      
C--J=1,2,...      
      do 14, J=1,Jmax

C---------------------------------------------------------      
      if(J.eq.1)write(*,10003)
      if(J.eq.2)write(*,10004)
      if(J.eq.3)write(*,10005)
      if(J.ge.4) then
          write(*,10006) char(J+67),J,char(J+67),
     &                   J,char(J+66),J, char(J+68),J,J
      end if
      write(*,10008)
      
      
 
      write(kwrite,10012)J
C----------------------------------------------------------
            
      
C--do-lop over diferent energies in lab system (Tlab=qEq(i))      
      do 15, i=1,NExper             

      q(N+1)=dsqrt( (2.d0*w2**2*w1*qEx(i)+w2**2*qEx(i)**2)
     &            /(w1**2+w2**2+2.d0*w2*(qEx(i)+w1)) )
     
C-----nna13 potential wants all the q in a different in a spesial way      
      if(ipot.eq.2)  call nna13Spesial(q)      
                          
      call GetPhaS(q,w,J,1,pha1,el1)                          
      call GetPhaS(q,w,J,2,pha2,el2)   
      call GetPhaCC(q,w,J,pha3,el3,pha4,el4,bepj)
      
      
      write(*,10000)qEx(i),pha1,el1,pha2,el2,pha3,el3,pha4,el4,bepj 
      write(kwrite,10011)qEx(i),pha1,pha2,pha3,pha4,dble(bepj),
     &                          el1,el2,el3,el4,dimag(bepj)

   15 continue 
      write(kwrite,*)' '
   14 continue    

      
99999 stop
      end
C============= End of main-program==================================
C===================================================================




      
C*    ipot is defined in main program.                                **
C*    This subrutine finds  the potential choosen in the main program.**
C*    The program can be made faster by using the potential directly  ** 
C*     in the program, instead of using this subroutine.              **
C*    Mesh points are added qF if the phase shift is calculated       **
C*     in medium. (qF is the Fermi momentum).                         **    
C***********************************************************************    
      subroutine potential( q1,q2,JJ,ii,cv)
      real*8 q1, q2  
      complex*16 cv(6)
      common /cipot/ ipot 
 
C-----Chooses different potentials      
      if(ipot.eq.1) then 
               call bonnv(q1,q2,JJ,ii,cv) 
      else if(ipot.eq.2) then
               call nna13v(q1,q2,JJ,ii,cv)
      else if(ipot.eq.3) then 
               call cbox(q1,q2,cv) 
      else if(ipot.eq.4) then
               call cpara(q1,q2,cv)
      else
              write(*,*)'wrong potential in subroutine potential'
      end if
      
      
      return
      end


      

      
C*    Only used by the "nna13" potential                            **
C*    The nna13 potential wants the q written in a spesial way      **
C*    The input parameter differ if it is NN scattering in          **
C*     medium or free.                                              **      
C*********************************************************************      
      subroutine nna13Spesial(q)
      implicit real*8 (A-H,O-Z), complex*16 (c)
      include "constants.h"
      real*8 q(N+1),qx(97),q0qmev,qF,starM,m
      logical noced,vinmedium,inmedium
      common/cmass/redM
      common /cpts/   qx(97) !"N+1" mesh points 
      common/cinmedium/starM,qF,inmedium,vinmedium
      common /cpoted/ q0qmev,qfmev,pmev,uanspp,wsnspp,ucnspp,udnspp,
     &                znrl,zrel,smev,noced


C-----Input variables in the "nna13" potential are given their values.        

      if(vinmedium.eq..true.)then
              m=starM
              qfmev=qF
      else
              m=redM
              qfmev=0.d0
      end if
        
      do 1814,i=1,N+1
      qx(i)=q(i)  
1814  continue
         
      
      noced=.false. 
      q0qmev=qx(N+1)**2
      pmev=0.d0
      smev=0.d0
      znrl=0.d0
      uanspp=0.d0
      wsnspp=m
      ucnspp=0.d0
      udnspp=0.d0
      zrel=2.d0*dsqrt(q0qmev+m**2)
      return
      end
      


      

      
C*    NOT IN USE AT THE MOMENT                                      **
C*    Changes from "Stap/bar" too "Blatt-Biedenharn"- formalism     **
C*    Phase shift and mixing factor are given in degees             **
C*    pha3=delta-  and  pha4=delta+                                 **      
C*    This method does not work with complex pot (go back to def)   **
C*********************************************************************      
      subroutine useBlatt(pha3,pha4,epj)
      include "constants.h" 
      real*8 pha3,pha4,temp
      complex*16 epj,ctemp
      logical degfixed
      common/boundfix/degfixed

C-----because 3S1 has boundstate (for low enegies) something need to be fixed
      
 
      if(degfixed.eq..true.)then
!      if(pha3.gt.90.d0) then 
              epj=-dble(epj)*(1.d0,0.d0)+dimag(epj)*(0.d0,1.d0)
              pha3=pha3-180.d0
      end if
      
      
      pha3=pha3/180.d0*pi
      pha4=pha4/180.d0*pi
      
      epj=(1.d0,0.d0)*dble(epj)/180.d0*pi
!     &     +(0.d0,1.d0)*dimag(epj)
      
      ctemp=cdsin(2.d0*epj)/cdcos(2.d0*epj)/dsin(pha3-pha4)*(0.d0,1.d0)
      ctemp=1/4.d0*(0.d0,-1.d0)*cdlog((1.d0+ctemp)/(1.d0-ctemp))
      
      temp=0.5d0*(dasin(dble(cdsin(2.d0*epj)/cdsin(2.d0*ctemp)))
     &    +pha3+pha4)
           
      
      pha4=(pha4+pha3-temp)*180.d0/pi
      pha3=(temp)*180.d0/pi 
      epj=(1.d0,0.d0)*dble(ctemp)*180.d0/pi+dimag(ctemp)*(0.d0,1.d0) 

      if(degfixed.eq..true.) then
              pha3=pha3+180.d0
              epj=-dble(epj)*(1.d0,0.d0)+dimag(epj)*(0.d0,1.d0) 
      end if
      return
      end


      
C*    Finds u                                                      **
C*    In-medium is different than free scattering. Thesse          **
C*     diffrences are aproximated by Mstar (effective mass) whitch **
C*     is a Dispersian effect. And qF (Fermi momentum) due to      **
C*     the Pauli effect.                                           **
C********************************************************************      
      subroutine findu(q,w,u)
      include "constants.h"
      real*8 starM,qF
      logical propag,inmedium,vinmedium 
      real*8 q(N+1),w(N),sm,xint,yint
      complex*16 u(N+1) 
      common /cinmedium/ starM,qF,inmedium,vinmedium 
!      common/cmass/redM
      common/cprop/propag 

      u(N+1)=(0.d0,0.d0)
      
      sm=starM
      
C-----In medium or free NN scattering
      if(inmedium .eq. .true.) then
C----------Extra integral from the "add zero trick"           
           xint=-(1.d0,0.d0)*sm*q(N+1)/pi
     &          *dlog(dabs((qF+q(N+1))/(qF-q(N+1))))
           yint=-(1.d0,0.d0)*q(N+1)*dsqrt(sm**2+q(N+1)**2)/pi
     &          *dlog(dabs((qF+q(N+1))/(qF-q(N+1))))
      else 
           xint=0.d0
           yint=0.d0
      end if
     

      
C-----minimal relativety             
      if(propag .eq. .true.) then
          do 150, i=1,N        
          u(i)=(1.d0,0.d0)*2.d0/pi*w(i)*q(i)**2*sm/(q(N+1)**2-q(i)**2)
          u(N+1)=u(N+1)+(1.d0,0.d0)*w(i)/(q(N+1)**2-q(i)**2)
  150     continue
          u(N+1)=-u(N+1)*2.d0*q(N+1)**2*sm/pi-(0.d0,1.d0)*sm*q(N+1)+xint
      else
              
C-----Thompson
          do 151, i=1,N        
          u(i)=(1.d0,0.d0)*( w(i)/pi*q(i)**2*(dsqrt(sm**2+q(i)**2)
     &                +dsqrt(sm**2+q(N+1)**2) )/(q(N+1)**2-q(i)**2) ) 
          u(N+1)=u(N+1)+(1.d0,0.d0)*w(i)/(q(N+1)**2-q(i)**2)
  151     continue    
          u(N+1)=-u(N+1)*q(N+1)**2*dsqrt(sm**2+q(N+1)**2)*2.d0/pi
     &           -(0.d0,1.d0)*q(N+1)*dsqrt(sm**2+q(N+1)**2)+yint 
 
      end if
       
      return
      end

             
    


C*    Test if qEx=(any of the mesh points by +-0.05)              **
C*    Returns .false. if bad, and .true. if OK                    **
C*                                                                **             
C*    If a "mesh point"= qEx ==> propagator singulare and         **
C*     gauss-legendre method will not give a good answear         **      
C*******************************************************************
      logical function meshBAD(q,qEx,nqEx)
      real*8 w1,w2
      include "constants.h"
      integer nqEx
      real*8 q(N+1),qEx(nqEx)
      logical BAD
      common/cmassp/w1,w2

      
      BAD=.false.
      do 200,i=1,nqEx
      
      q(N+1)=dsqrt( (2.d0*w2**2*w1*qEx(i)+w2**2*qEx(i)**2)
     &            /(w1**2+w2**2+2.d0*w2*(qEx(i)+w1)) )
      
      
      do 210, j=1,N
      
      if(ANINT(1.d1*q(j)) .EQ. ANINT(1.d1*q(N+1)) ) then !d1=>+-0.05 
              write(*,*)'change N,mesh points ended up bad compared
     &                   with experimental data' 
              BAD=.true.
              go to 299     ! terminates program
      end if
  210 continue
  200 continue
  299 meshBAD=BAD
      return
      end  


      
C*    Calculates "T-Matrix"- elements: (simple matrix multipications) **
C*    T= A*V or AA*VV                                                 **
C*    Returns the interesting values of  T=A*V                        **       
C***********************************************************************     
      subroutine AmultV(A,V,ii,ndim,T)
      complex*16 T(ii,ii),  A(ndim,ndim),   V(ndim,ii)
      
      T(1,1)=0.d0
      if(ii .eq. 2) then 
          T(1,2)=0.d0 
          T(2,1)=0.d0
          T(2,2)=0.d0
      end if
          
      do 420,k=1,(ndim)
    
      if(ii .eq. 1) then
          T(1,1)=T(1,1)+ A(ndim,k)*V(k,1)
      else
          T(1,1)=T(1,1)+ A(ndim/2,k)*V(k,1)
          T(1,2)=T(1,2)+ A(ndim/2,k)*V(k,2)
          T(2,1)=T(2,1)+ A(ndim,k)*V(k,1)
          T(2,2)=T(2,2)+ A(ndim,k)*V(k,2)
      end if
  
  420 continue
      return
      end

      

C**   Calculates phase shifts(pha), inelasticity(el)                 **
C**   and complex mixing fator(bepj) for coupled channels            **
C**                                                                  ** 
C**   "propag" decides if calulate S-matrix from "minimal-           **
C**   relativety" or Thopson eq.                                     ** 
C**                                                                  **      
C**   Results in "STAPP/bar" or "Blatt-Biedenharn" formalism         **
C**                                                                  **
C**   Returns   pha3,el3,pha4,el4,bepj                               **
C**                                                                  **
C**   Phase shift and Re(mixing factor) are given in degees          **   
C**   Uses matric invertion fuction from IMSL: "DLINRG()"            ** 
C**   pha3=delta-  and  pha4=delta+                                  **
C**   el3=eta-  and  pha4=eta+                                       **
C**********************************************************************       

      
C1      
      subroutine GetPhaCC(q,w,JJ,pha3,el3,pha4,el4,bepj)
      logical propag,kowalski,degfixed,vinmedium
      real*8 starM
      include "constants.h"
      complex*16 A(2*N+2,2*N+2),Ainv(2*N+2,2*N+2),V(2*N+2,2)
      complex*16 T(2,2),S(2,2),x,bepj,y,z
      real*8 w(N),q(N+1),pha3,el3,pha4,el4,factor
      common /cinmedium/ starM,qF,inmedium,vinmedium
      common/clseq/kowalski
      common/cprop/propag
      common/boundfix/degfixed
      common/formalism/blatt
      
      n1=N+1
      n2=2*(N+1)
      

      
C-----Two different methods to solve LS. eq

      if (kowalski .eq. .true.) then
C1----Kowalski method
           call findCTkow(q,w,JJ,T)
           go to 2333

      else
C2----normal method
           
      call makeMM(A,V,q,w,JJ)  !makes A- and V-Matrixes
      
      call DLINCG(n2,A,n2,Ainv,n2) !returns A-invese
     
      call AmultV(Ainv,V,2,n2,T)  !Returns T= Ainv*V (="T-matix")
      
      end if
 
 2333 continue



      
       
      if(propag .eq. .true.) then
C-------------"minimal relativety"               
              factor=starM*2.d0*q(n1)
      else
              
C-------------Thompson 
              factor=dsqrt(starM**2+q(n1)**2)*2.d0*q(n1)
      end if
       


       
       
C------------ S-matrix 
       
      S(1,1)=(1.d0,0.d0)-(0.d0,1.d0)*T(1,1) *factor
      S(2,2)=(1.d0,0.d0)-(0.d0,1.d0)*T(2,2) *factor
      S(2,1)=           -(0.d0,1.d0)*T(2,1) *factor
      S(1,2)=           -(0.d0,1.d0)*T(1,2) *factor
      
     
          
      
      
C------------ Stapp formalism 
      if(blatt.eq..false.)then
      
        x=(S(2,1)+S(1,2))/(2.d0*cdsqrt(S(2,2)*S(1,1)))
        bepj=(0.d0,-1.d0)/4.d0*cdlog((1.d0+x)/(1.d0-x))
      
              
        pha3=0.5d0*(datan2(dimag(S(2,2)/cdcos(2.d0*bepj))
     &                   ,dble(S(2,2)/cdcos(2.d0*bepj))))*180.d0/pi

      
        pha4=0.5d0*(datan2(dimag(S(1,1)/cdcos(2.d0*bepj))
     &                   ,dble(S(1,1)/cdcos(2.d0*bepj))))*180.d0/pi
     
            

        el3=cdabs(S(2,2)/cdcos(2.d0*bepj))
        el4=cdabs(S(1,1)/cdcos(2.d0*bepj))
            
      else
              
C-----"Blatt-Biedenharn"- formalism
        S(2,2)=(S(2,2)-1.d0)/(0.d0,2.d0)
        S(1,1)=(S(1,1)-1.d0)/(0.d0,2.d0)
        S(1,2)=(S(1,2))/(0.d0,2.d0)
        S(2,1)=(S(2,1))/(0.d0,2.d0) 
        x=(0.d0,1.d0)*(S(1,2)+S(2,1))/(S(2,2)-S(1,1))
        bepj=(0.d0,-1.d0)/4.d0*cdlog((1.d0+x)/(1.d0-x))
        
        y=0.5d0*(S(1,1)+S(2,2))       
        z=cdsqrt(y*y-S(2,2)*S(1,1)+S(1,2)*S(1,2))
        
        pha3=0.5d0*datan2(dble(y+z),0.5d0-dimag(y+z))*180.d0/pi
        pha4=0.5d0*datan2(dble(y-z),0.5d0-dimag(y-z))*180.d0/pi
        el3=dsqrt(1.d0+4.d0*((dble(y+z))**2+(dimag(y+z))**2-dimag(y+z)))
        el4=dsqrt(1.d0+4.d0*(dble(y-z)**2+dimag(y-z)**2-dimag(y-z)))
      end if



      
      
C-----defenition of real and imag part of the mixing parameter
      bepj=(1.d0,0.d0)*dble(bepj)*180.d0/pi 
     &       +(0.d0,1.d0)*dimag(bepj) 
     



      
C-----fixes 3S1-state for low energies (because it's a bound state) 
        degfixed=.false.
        if(JJ.eq.1 .and.q(n1).le.6.d2) then !FIX BLATTTTTTTTTTTTTTTTTTTT
            if(dble(S(2,1)+S(1,2)).le.-1.d-10.and.
     $      dimag(cdsqrt((S(2,2)*S(1,1)))).le.-1.d-10) then
                bepj=-dble(bepj)*(1.d0,0.d0)+dimag(bepj)*(0.d0,1.d0)
                pha3= pha3+180.d0
                degfixed=.true.
            end if
        end if
      return
      end

      



      
C2      
C**   Calculates/Returns phase shifts(pha) and inelasticity(el)     **
C**                                                                 **      
C**   "propag" decides if calulate S-matrix from "minimal-          **
C**   relativety" or Thopson eq.                                    **
C**                                                                 **       
C**   Uncoubled channels                                            **
C**   Returns pha,el                                                **
C*********************************************************************      
      subroutine GetPhaS(q,w,JJ,ii,pha,el)
      logical propag,kowalski,vinmedium
      real*8 starM 
      include "constants.h"
      complex*16 A(N+1,N+1),Ainv(N+1,N+1),V(N+1),S,T
      real*8 q(N+1),w(N),pha,el,factor
      common /cinmedium/ starM,qF,inmedium,vinmedium
      common/cprop/propag
      common/clseq/kowalski
      
      n1=N+1



C-----Two different methods to solve LS. eq 
      
      if (kowalski .eq. .true.) then 
C1----Kowalski method
           call findTkow(q,w,ii,JJ,T)
           go to 2233
                                                         
      else
C2----normal method
     
      call makeAV(A,V,q,w,JJ,ii)  !makes and V and most of A 
      
      call DLINCG(n1,A,n1,Ainv,n1) !IMSL, returnsA-inverse
     
      call AmultV(Ainv,V,1,n1,T)  !Returns T= Ainv*V ="T-matix"

      end if
      
 2233 continue     
      

      

      if(propag .eq. .true.) then
C------------- "minimal relativety" eq. 
              factor = starM*2.d0*q(n1)
      else
C-------------"Thomson" eq. 
              factor = dsqrt(starM**2+q(n1)**2)*2.d0*q(n1)
      end if
 
      
 
              
C---------- S-matrix 
              
      S=+(1.d0,0.d0)  - (0.d0,1.d0)*T *factor
   
      
            
C-----------phase shift and inelasticity 
            
      pha=0.5d0*(datan2(dimag(S),dble(S)))*180.d0/pi
      el=dsqrt( dble(S*conjg(S)) )
           
C        test if S is  unitare (=1)
C        write(*,*)'?one?=', S*conjg(S)
      return
      end


      
      

C**   Calculate/Returns A and V matrixes                          **
C**   (spinn singlet andspinn triplet)                            **
C**   q(N) mesh points, q(N+1)="Lab energy" in MEV                **
C**   w(N) weights                                                ** 
C**   JJ = "J" total angular momentum                             **      
C**   ii is used in BONN B potensial V(ii) 1<=ii<=6               **   
C*******************************************************************
      subroutine makeAV(A,V,q,w,JJ,ii)
      include "constants.h" 
      complex*16 A(N+1,N+1), V(N+1),u(N+1)
      real*8 q(N+1),w(N)
      complex*16 cv(6)
      
      call findu(q,w,u)           
       
      do 510,i=1,(N+1)
      do 520,j=1,(N+1)
!      call bonnv( q(i),q(j),JJ,ii,cv)
      call potential( q(i),q(j),JJ,ii,cv)
      
      A(i,j)=-u(j)*cv(ii) 
  520 continue
  510 continue
      do 530,i=1,N+1
      A(i,i)=A(i,i)+(1.d0,0.d0)
!      call bonnv( q(i),q(N+1),JJ,ii,cv)
      call potential( q(i),q(N+1),JJ,ii,cv)
      V(i)=cv(ii)
      
  530 continue
      return
      end
      

 

      
C**   Calculates Big A anv V matrixses ===> "AA" and "VV"         **      
C**                   ( Coubled Channels )                        **
C**   q(N) mesh points, q(N+1)="Lab energy" in MEV                **
C**   w(N) weights                                                **      
C**   JJ = "J" total angular momentum                             **
C**   ii is used in BONN B potensial V(ii) 1<=ii<=6               **
C*******************************************************************
      subroutine makeMM(AA,VV,q,w,JJ)
      include "constants.h" 
      complex*16 AA(2*(N+1),2*(N+1)),  VV(2*(N+1),2),u(N+1)
      real*8 q(N+1),w(N)
      complex*16 cv(6)

      
      n1=N+1
      n2=2*(N+1)
  
      
      call findu(q,w,u)
      
      do 620,i=1,n1
      do 630,j=1,n1
      
!      call bonnv( q(i),q(j),JJ,6,cv) !6: means coubled chanels in bonnpot
      call potential( q(i),q(j),JJ,6,cv)
      
      AA(i,j)=-u(j)*cv(3) 
      AA(i+n1,j+n1)=-u(j)*cv(4) 
      AA(i+n1,j)=-u(j)*cv(6)  
      AA(i,j+n1)=-u(j)*cv(5)
      
  630 continue
  620 continue
           
      do 650, i=1,n1
      
      AA(i,i)=AA(i,i)+(1.d0,0.d0)
      AA(i+n1,i+n1)=AA(i+n1,i+n1)+(1.d0,0.d0) 

!      call bonnv( q(i),q(n1),JJ,6,cv)
      call potential(q(i),q(n1),JJ,6,cv)
      
      VV(i,1)=cv(3)
      VV(i+n1,2)=cv(4) 
      VV(i+n1,1)=cv(6) 
      VV(i,2)=cv(5) 
      
  650 continue    
      return
      end


C**   Calls bonn potensial and returns a bonn potensial...         **
C**   ....for the speified system.                                 **
C**                                                                **
C**   Returns "cv" (see souroutine bonn)                           ** 
C**                                                                ** 
C**   JJ= totaly angular momentum                                  **
C**   q1 and q2 are impulses in momentum space                     **
C**   ii is a parameter in bonn potensial: V(ii)                   **
C********************************************************************
      subroutine bonnv(q1,q2,JJ,ii,cv) 
      IMPLICIT REAL*8 (A-H,O-Z)
      include "constants.h"
      complex*16 cv(6)
      logical propag
      common/cprop/propag
C
      COMMON /CRDWRT/ KREAD,KWRITE,KPUNCH,KDA(9)
C
C
C        ARGUMENTS AND VALUES FOR THE POTENTIAL SUBROUTINE
C
      COMMON /CPOT/   V(6),XMEV,YMEV
      COMMON /CSTATE/ J,HEFORM,SING,TRIP,COUP,ENDEP,LABEL
      LOGICAL HEFORM,SING,TRIP,COUP,ENDEP
C        THIS HAS BEEN THE END OF THE ARGUMENTS AND VALUES OF THE
C        POTENTIAL SUBROUTINE
C
C
      INTEGER LSJ/'LSJ'/,HEL/'HEL'/
      LOGICAL SWITCH/.FALSE./
c
c
      if(propag .eq. .true.)then
C-----------minimal relativety      
            open (unit=5,file='dbonn.d')
      else
C-----------"Thompson"
            open (unit=5,file='dbonntomp.d')
      end if
      
      open (unit=6,file='potbonn.d')
C
C        SET THE FOLLOWING PARAMETERS ONCE FOR EVER
C
      KREAD=5
      KWRITE=6
C
C        THE FOLLOWING STATEMENT MEANS, THAT THE POTENTIAL
C        WILL BE GIVEN IN THE LSJ-FORMALISM
      HEFORM=.FALSE.
C
      NAME=LSJ
      IF (HEFORM) NAME=HEL
C
C        THE FOLLOWING THREE STATEMENTS MEAN, THAT FOR EACH J THE
C        POTENTIAL WILL BE GIVEN FOR THE STATES SINGLET,
C        UNCOUPLED TRIPLET AND THE FOUR COUPLED TRIPLET CASES
      
      SING=.FALSE. 
      TRIP=.FALSE. 
      COUP=.FALSE.
      
      if(ii .EQ. 1) SING=.TRUE.
      
      if(ii .EQ. 3 .AND. JJ .EQ. 0) TRIP=.TRUE. ! these two are swithed !!!!!
      if(ii .EQ. 3 .AND. JJ .EQ. 0) COUP=.TRUE. !????????????????? Ask why!!!!!!!!!!!
      
      if(ii .EQ. 2 .AND. JJ .NE. 0) TRIP=.TRUE. 
      if(ii .EQ. 3 .AND. JJ .NE. 0) COUP=.TRUE.
      if(ii .EQ. 4 .OR. ii .EQ. 5 .OR. ii .EQ. 6) COUP=.TRUE.
      
      J=JJ
      XMEV=q1   
      YMEV=q2
      
      call bonn

      do 700, i=1,6

      cv(i)=V(i)*pi/2.d0*(1.d0,0.d0) ! pi/2 is a normalization factor

      
C**** This is done in potensial. have to use use diferent dbonn.d files      
C-----minimal relativety     
C      cv(i)=(V(i)*pi/2.d0 , 0.d0) ! pi/2 is a normalization factor
C-----Thompson  (change from "minimal relativety"        
C      if(propag .eq. .false.) then
C          cv(i)=cv(i)/redM*dsqrt(dsqrt((redM**2+q1**2)*(redM**2+q2**2)))
C      end if
      
      
  700 continue  
      return
      end








      
C********************************************************************
C**   Calls nna13 potensial and returns a hi energy bonn potensial **
C**   ....for the speified system.                                 **
C**                                                                **
C**   Returns "cv" (see souroutine nna13)                          ** 
C**                                                                ** 
C**   JJ= totaly angular momentum                                  **
C**   q1 and q2 are impulses in momentum space                     **
C**   ii is a parameter in nna13 potensial: V(ii)                  **
C********************************************************************
      subroutine nna13v(q1,q2,JJ,ii,cvout) 
      IMPLICIT REAL*8 (A-H,O-Z), complex*16 (c)
      include "constants.h"
      complex*16 cvout(6)
      logical propag,inmedium,vinmedium
      common/cprop/propag
      common /cinmedium/ starM,qF,inmedium,vinmedium
C
      COMMON /CRDWRT/ KREAD,KWRITE,KPUNCH,KDA(9)
C
C
C        ARGUMENTS AND VALUES FOR THE POTENTIAL SUBROUTINE
C
!      common /cptsc/  cqx(97)
      common /cpts/   qx(97),dx,nx1,ix,iy
      common /cpotc/  cv(6),cxmev,cymev
      common /cstate/ j,heform,sing,trip,coup,endep,label
      common /cpoted/ q0qmev,qfmev,pmev,uanspp,wsnspp,ucnspp,udnspp,
     &                znrl,zrel,smev,noced 
      logical heform,sing,trip,coup,endep 
C        THIS HAS BEEN THE END OF THE ARGUMENTS AND VALUES OF THE
C        POTENTIAL SUBROUTINE
C
C
      INTEGER LSJ/'LSJ'/,HEL/'HEL'/
      LOGICAL SWITCH/.FALSE./
c
c
!      if(inmedium .eq. .true..and. vinmedium.eq..true.)then
      if(vinmedium.eq..true.)then
              open (unit=5,file='dd85in.d')            
      else
              open (unit=5,file='dd85.d')
      end if
      open (unit=6,file='potbonnhi.d')
C
C        SET THE FOLLOWING PARAMETERS ONCE FOR EVER
C
      KREAD=5
      KWRITE=6
C
C        THE FOLLOWING STATEMENT MEANS, THAT THE POTENTIAL
C        WILL BE GIVEN IN THE LSJ-FORMALISM
      HEFORM=.FALSE.
C
      NAME=LSJ
      IF (HEFORM) NAME=HEL
C
C        THE FOLLOWING THREE STATEMENTS MEAN, THAT FOR EACH J THE
C        POTENTIAL WILL BE GIVEN FOR THE STATES SINGLET,
C        UNCOUPLED TRIPLET AND THE FOUR COUPLED TRIPLET CASES
      
      SING=.FALSE. 
      TRIP=.FALSE. 
      COUP=.FALSE.
      
      if(ii .EQ. 1) SING=.TRUE.
      
      if(ii .EQ. 3 .AND. JJ .EQ. 0) TRIP=.TRUE. ! these two are swithed !!!!!
      if(ii .EQ. 3 .AND. JJ .EQ. 0) COUP=.TRUE. !????????????????? Ask why!!!!!!!!!!!
      
      if(ii .EQ. 2 .AND. JJ .NE. 0) TRIP=.TRUE. 
      if(ii .EQ. 3 .AND. JJ .NE. 0) COUP=.TRUE.
      if(ii .EQ. 4 .OR. ii .EQ. 5 .OR. ii .EQ. 6) COUP=.TRUE.
      J=JJ
      nx1=N+1



      
C-----Very slow way to find ix and iy ...
      do 785,i=1,N+1
          if( q1.eq.qx(i) )  then
                  ix=i
                  go to 786
          end if
  785 continue
  786 continue         
      do 787,i=1,N+1
          if( q2.eq.qx(i) )then
                  iy=i
                  go to 788
          end if
  787 continue
  788 continue 
      

      
      cxmev=q1*(1.d0,0.d0)   !qx(ix)   
      cymev=q2*(1.d0,0.d0)   !qx(iy)
!      write(*,*)'914',ix,iy,cxmev,cymev,q0qmev,nx1,qx(33)**2 
      call nna13
!      write(*,*)'949',cv(1)
      do 789, i=1,6
      

      cvout(i)=cv(i)*pi/2.d0 !  pi/2 is a normalization factor

      
C**** This is done in potensial. have to use use diferent dbonn.d files      
C-----minimal relativety     
C      cv(i)=(V(i)*pi/2.d0 , 0.d0) ! pi/2 is a normalization factor
C-----Thompson  (change from "minimal relativety"        
C      if(propag .eq. .false.) then
C          cv(i)=cv(i)/redM*dsqrt(dsqrt((redM**2+q1**2)*(redM**2+q2**2)))
C      end if
      
      
  789 continue  
      return
      end

      
C**********************************************************************






      
      
C====>  all spesial Kovalski methods are comming now            <====== 





      


C**   calculate T from Kovalski's equations                       **
C*******************************************************************      
      subroutine findTkow(q,w,ii,JJ,T)
      include "constants.h"
      real*8 q(N+1),w(N), u(N)
      complex*16 T
      real*8 starM,xint,qF,factor
      logical inmedium, propag,vinmedium
      complex*16 B(N,N),V((N+1),(N+1)),Ts(N)
      complex*16 Binv(N,N),x(N),zum, cv(6)
      common /cinmedium/ starM,qF,inmedium,vinmedium
      
      common/cprop/propag  
      
      n1=N+1

      if(propag .eq. .true.) then
C-------------"minimal relativety"
               factor=starM*q(n1)         
               fac=1.d0
      else
C--------------Thompson
               factor=dsqrt(starM**2+q(n1)**2)*q(n1)
               fac=2.d0*dsqrt(starM**2+q(n1)**2) 
     &              /(dsqrt(starM**2+q(n1)**2)+dsqrt(starM**2+q(i)**2)) 
      end if

      
      
C-----integral term in-medium calculations
      if(inmedium .eq. .true.) then         
           xint=-(1.d0,0.d0)*factor/pi
     &          *dlog(dabs((qF+q(n1))/(qF-q(n1))))
      else
           xint=0.d0
      end if
      
C-----Calculates the potential V            
      do 5510,i=1,(N+1)
      do 5520,j=1,(N+1)
!      call bonnv( q(i),q(j),JJ,ii,cv)
      call potential( q(i),q(j),JJ,ii,cv)
      V(i,j)=cv(ii)
 5520 continue
 5510 continue 
     
      

C-----Calculates u and x      
      if(propag .eq. .true.) then      
C----------"minimal relativety"
           do 5513,i=1,N 
           u(i)=2.d0*starM/pi*w(i)/(q(n1)**2-q(i)**2)
           x(i)=V(i,n1)/V(n1,n1)
 5513      continue
      else
C-------------Thompson
           do 5514,i=1,N
           u(i)=(dsqrt(starM**2+q(n1)**2)+dsqrt(starM**2+q(i)**2))
     &           /pi*w(i)/(q(n1)**2-q(i)**2)
           x(i)=V(i,n1)/V(n1,n1)
 5514      continue
      end if
           
           
C-----Calculates B
      do 5511,i=1,N
      do 5521,j=1,N
      B(i,j)=-( V(i,j)- x(i)*V(n1,j)  )*u(j)*q(j)**2
 5521 continue
 5511 continue
      
      do 5530,i=1,N
      B(i,i)=B(i,i)+(1.d0,0.d0)                                    
 5530 continue
            
            
C-----Calculate the inverse of B  
      call DLINCG(N,B,N,Binv,N) !IMSL, returns B-inverse       


C-----Matrix multiplication of B whith x     
      do 5400,i=1,N
      zum=(0.d0,0.d0)
      do 5420,j=1,N
      zum=zum+ Binv(i,j)*x(j)
 5420 continue
      Ts(i)=zum
 5400 continue
            
      

      
      
C-----Calculate the on-shell T from eq(4.52) 
      zum=(0.d0,0.d0)

      if(propag .eq. .true.) then
C---------"minimal relativety"
           do 5440,i=1,N 
           zum=zum+u(i)*( q(i)**2*V(n1,i)*Ts(i) -  q(n1)**2*V(n1,n1) )
 5440      continue     
      else
C----------Thompson               
           do 5441,i=1,N
           zum=zum+u(i)*( q(i)**2*V(n1,i)*Ts(i) -  q(n1)**2*V(n1,n1)
     &         /(dsqrt(starM**2+q(n1)**2)+dsqrt(starM**2+q(i)**2))
     &         *2.d0*dsqrt(starM**2+q(n1)**2) )

 5441      continue 
      end if
      
      T=V(N+1,N+1) /( (1.d0,0.d0)-zum 
     & +(0.d0,1.d0)*factor*V(n1,n1)-xint*V(n1,n1))
          
      return 
      end

      

      
C*****Coublet Channels*************************************************
      
      subroutine findCTkow(q,w,JJ,T)
      include "constants.h"
      real*8 q(N+1),w(N), u(N)
      complex*16 T(2,2)
      real*8 redM
      complex*16 B(2*N,2*N),V(2*(N+1),2*(N+1)),Ts(2*N,2)
      complex*16 Binv(2*N,2*N),Vinv(2,2),zum1,zum2,x,y,A(2,2)
      complex*16 zum11,zum12,zum21,zum22
      real*8 bv(6)
      common/cmass/redM
      
      n1=N+1
      n2=2*(N+1)
      
C-----Calculates B and x            
      do 6508,i=1,n1
      do 6509,j=1,n1
!      call bonnv( q(i),q(j),JJ,6,bv)
      V(i,j)=bv(3)*(1.d0,0.d0 )
      V(i+n1,j+n1)=bv(4)*(1.d0,0.d0 )
      V(i+n1,j)=bv(6)*(1.d0,0.d0 )
      V(i,j+n1)=bv(5)*(1.d0,0.d0 )
 6509 continue
 6508 continue 
      do 6513,i=1,N 
      u(i)=2.d0*redM/pi*w(i)/(q(n1)**2-q(i)**2)
 6513 continue

C-----inverse of the 2x2-matrix VV(k,k)
      x=1.d0/( V(n1,n1)*V(n2,n2)-V(n1,n2)*V(n2,n1)  )
      
      Vinv(1,1)= x*V(n2,n2)
      Vinv(1,2)=-x*V(n1,n2)
      Vinv(2,1)=-x*V(n2,n1)
      Vinv(2,2)= x*V(n1,n1)
      
      do 6511,i=1,N
      do 6521,j=1,N
       
      y=u(j)*q(j)**2 
      
      B(i,j)      =-( V(i,j)      -V(i,n1)   *Vinv(1,1)*V(n1,j)     )*y
      B(i+n1,j)   =-( V(i+n1,j)   -V(i+n1,n2)*Vinv(2,1)*V(n1,j)     )*y
      B(i,j+n1)   =-( V(i,j+n1)   -V(i,n1)   *Vinv(1,2)*V(n2,j+n1)  )*y
      B(i+n1,j+n1)=-( V(i+n1,j+n1)-V(i+n1,n2)*Vinv(2,2)*V(n2,j+n1)  )*y
 6521 continue
 6511 continue
      
      do 6530,i=1,2*n
      B(i,i)=B(i,i)+(1.d0,0.d0)                                    
!      B(i,2*N)=B(i,2*N)+(1.d0,0.d0)
 6530 continue
            
C-----Calculate the inverse of B  
      call DLINCG(2*N,B,2*N,Binv,2*N) !IMSL, returns B-inverse       

C-----Matrix multiplication   (Binv x V x Vinv)  
      do 6400,i=1,2*N
      zum1=(0.d0,0.d0)
      zum2=(0.d0,0.d0)
      do 6420,j=1,2*N
      zum1=zum1 + Binv(i,j) *V(j,n1)*Vinv(1,1)
     &          + Binv(i,j) *V(j,n2)*Vinv(2,1)
      zum2=zum2 + Binv(i,j) *V(j,n1)*Vinv(1,2)
     $          + Binv(i,j) *V(j,n2)*Vinv(2,2)
 6420 continue
      Ts(i,1)=zum1
      Ts(i,2)=zum2
 6400 continue
            
C-----Calculate the on-shell T from eq(4.52) 
      zum11=(0.d0,0.d0)
      zum12=(0.d0,0.d0)
      zum21=(0.d0,0.d0)
      zum22=(0.d0,0.d0)
      do 6440,i=1,N
      zum11=zum11+u(i)*( q(i)**2*V(n1,i)*Ts(i,1)- q(n1)**2*V(n1,n1) )
      zum12=zum11+u(i)*( q(i)**2*V(n1,i)*Ts(i,2)- q(n1)**2*V(n1,n2) ) 
      zum21=zum11+u(i)*( q(i)**2*V(n2,i)*Ts(i,1)- q(n1)**2*V(n2,n1) ) 
      zum22=zum11+u(i)*( q(i)**2*V(n2,i)*Ts(i,2)- q(n1)**2*V(n2,n2) ) 
 6440 continue
      
      A(1,1)=(1.d0,0.d0)-2.d0*redM/pi*zum11
     $              +(0.d0,1.d0)*q(N+1)*redM*V(n1,n1)
      A(1,2)=(0.d0,0.d0)-2.d0*redM/pi*zum12
     $              +(0.d0,1.d0)*q(N+1)*redM*V(n1,n2)
      A(2,1)=(0.d0,0.d0)-2.d0*redM/pi*zum21
     $              +(0.d0,1.d0)*q(N+1)*redM*V(n2,n1)
      A(2,2)=(1.d0,0.d0)-2.d0*redM/pi*zum22
     $              +(0.d0,1.d0)*q(N+1)*redM*V(n2,n2)

C-----inverse of A
      Vinv(1,1)= 1.d0/( A(1,1)*A(2,2)-V(n1,n2)*V(n2,n1) )* A(2,2)
      Vinv(1,2)=-1.d0/( A(1,1)*A(2,2)-V(n1,n2)*V(n2,n1) )* A(1,2)
      Vinv(2,1)=-1.d0/( A(1,1)*A(2,2)-V(n1,n2)*V(n2,n1) )* A(2,1)
      Vinv(2,2)= 1.d0/( A(1,1)*A(2,2)-V(n1,n2)*V(n2,n1) )* A(1,1)

C-----matrix multiplication
      T(1,1)=Vinv(1,1)*V(n1,n1)+Vinv(1,2)*V(n2,n1)
      T(1,2)=Vinv(1,1)*V(n1,n2)+Vinv(1,2)*V(n2,n2)
      T(2,1)=Vinv(2,1)*V(n1,n1)+Vinv(2,2)*V(n2,n1)
      T(2,2)=Vinv(2,1)*V(n1,n2)+Vinv(2,2)*V(n2,n2)

          
      return 
      end


      


      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  gauss.f: Points and weights for Gaussian quadrature                 c
c                                                                      c
c  taken from: "Projects in Computational Physics" by Landau and Paez  c
c              copyrighted by John Wiley and Sons, New York            c
c                                                                      c
c  written by: Oregon State University Nuclear Theory Group            c
c              Guangliang He & Rubin H. Landau                         c
c              code copyrighted by RH Landau                           c
c  supported by: US National Science Foundation, Northwest Alliance    c
c                for Computational Science and Engineering (NACSE),    c
c                US Department of Energy                               c
c                                                                      c
c  comment: error message occurs if subroutine called without a main   c
c  comment: this file has to reside in the same directory as integ.c   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     rescale rescales the gauss-legendre grid points and weights
c
c     npts     number of points
c     job = 0  rescalling uniformly between (a,b)
c           1  for integral (0,b) with 50% points inside (0, ab/(a+b))
c           2  for integral (a,inf) with 50% inside (a,b+2a)
c     x, w     output grid points and weights.
c
      subroutine gauleg(npts,job,a,b,x,w)
      integer npts,job,m,i,j
      real*8 x(npts),w(npts),a,b,xi
      real*8 t,t1,pp,p1,p2,p3,aj
      real*8 eps,pi,zero,two,one,half,quarter
      parameter (pi = 3.14159265358979323846264338328, eps = 3.0E-14)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      parameter (half=0.5d0,quarter=0.25d0)
c
c *** FIRST EXECTUABLE *************************************************
c
 
      m=(npts+1)/2
      do 1020 i=1,m
         t=cos(pi*(i-quarter)/(npts+half))
 1000    continue
         p1=one
         p2=zero
         aj=zero
         do 1010 j=1,npts
            p3=p2
            p2=p1
            aj=aj+one
            p1=((two*aj-one)*t*p2-(aj-one)*p3)/aj
 1010    continue
         pp=npts*(t*p1-p2)/(t*t-one)
         t1=t
         t=t1-p1/pp
c
         if(abs(t-t1).gt.eps) goto 1000
c
         x(i)=-t
         x(npts+1-i)=t
         w(i)=two/((one-t*t)*pp*pp)
         w(npts+1-i)=w(i)
 1020 continue
c
c rescale the grid points
      if (job.eq.0) then
c     scale to (a,b) uniformly
         do 1030 i=1,npts
            x(i)=x(i)*(b-a)/two+(b+a)/two
            w(i)=w(i)*(b-a)/two
 1030    continue
      elseif (job.eq.1) then
c scale to (0,b) with 50% points inside (0,ab/(a+b))
         do 1040 i=1,npts
            xi=x(i)
            x(i)=a*b*(one+xi)/(b+a-(b-a)*xi)
            w(i)=w(i)*two*a*b*b/((b+a-(b-a)*xi)*(b+a-(b-a)*xi))
 1040    continue
      elseif (job.eq.2) then
c scale to (a,inf) with 50% points inside (a,b+2a)
         do 1050 i=1,npts
            xi=x(i)
            x(i)=(b*xi+b+a+a)/(one-xi)
            w(i)=w(i)*two*(a+b)/((one-xi)*(one-xi))
 1050    continue
      else
         pause 'Wrong value of job'
      endif
c
      return
      end    

      

C*    calculates the analytical phaseshift and the analytical      *
C*    inelasticity using SE in coordinate space                    *
C*******************************************************************
      subroutine anabox(k,pha,inel) !Analytic box-potesial
      include "constants.h"
      complex*16 V,E,k1,k2,xi,box
      real*8 k,a,pha,inel
            
      V=(-14.d0,-1.d0) ! V is the strength of the box potential
      a=2.7d0          ! a is the width of the box potential

      E=k**2/RedM*(1.d0,0.d0)
      k1=(RedM*(E-V))**0.5d0
      k2=(RedM*E)**0.5d0
      xi=(k2*sin(k1*a/chbar)/cos(k1*a/chbar)/k1)*(0.d0,1.d0)
      box=((0.d0,-1.d0)/2.d0*log(((1.d0,0.d0)+xi)/((1.d0,0.d0)-xi))
     & -k2*a/chbar)

3603  if(dble(box) .LE. -pi/2.d0)then
              box=box+(pi,0.d0)
              go to 3603
      end if

      pha=dble(box)*180.d0/pi
      inel=dexp(-2.d0*dimag(box))
      return
      end





C*    The partial potential in momentum space are calculated         *
C*     from complex box potential in coordinate space.               *
C*    "V" is the strengt, and "a" is the width of the box potential  *
C*    in coordinate space                                            *
C*    This potential is explaned in my theses                        *
C*********************************************************************
      subroutine cbox(q1,q2,vb)
      include "constants.h"
      real*8 a,q1,q2
      complex*16 V,vb
      
      V=(-14.d0,-0.d0) ! V is the strength of the box potential
      a=2.7d0          ! a is the width of the box potential

      if(q1.eq. q2) then
              vb=V/(2.d0*q1**2)*(a/chbar-sin(2.d0*q1*a/chbar)/(q1*2.d0))

      else
              vb=V/2.d0/q1/q2*(sin((q1-q2)*a/chbar)/(q1-q2)
     &             -sin((q2+q1)*a/chbar)/(q2+q1))

      end if

      return
      end



C*    The partial potential in momentum space are calculated         *
C*     from complex box potential in coordinate space.               *
C*    This potential is explaned in my theses                        *
C*********************************************************************
      subroutine cpara(k,km,Vkmk)
      include "constants.h"
      real*8 k,km,a,b,c,pionIM
      complex*16 Vkmk,Va,Vb

      Va=(-10.463d0,0.d0)
      Vb=(-1650.6d0,0.d0)
      Vc=(6484.3d0,0.d0)
      a=1.d0
      b=4.d0
      c=7.d0
      pionIM=0.7d0

      Vkmk=       Va*DLOG(((km+k)**2+(pionIM*a*chbar)**2)/
     &                    ((km-k)**2+(pionIM*a*chbar)**2))
      Vkmk=Vkmk + Vb*DLOG(((km+k)**2+(pionIM*b*chbar)**2)/
     &                    ((km-k)**2+(pionIM*b*chbar)**2))
      Vkmk=Vkmk + Vc*DLOG(((km+k)**2+(pionIM*c*chbar)**2)/
     &                    ((km-k)**2+(pionIM*c*chbar)**2))
      Vkmk=Vkmk/(chbar*pionIM*4.d0*k*km)
      return
      end






C***** Calculates the bound state energy ****************************






C*    Calculate the bound state energy for single and coubled channels * 
C***********************************************************************
      subroutine boundE(ii,JJ)
      integer ii,JJ
      write(*,*)'1599',ii,JJ
      if((JJ.gt.0) .and. (ii.gt.1))then
            call boundcc(ii,JJ)
      else  
            call bounds(ii,JJ)     
      end if
      
      return 
      end        
       


C*********************************************************************
C*    Calculate the bound state energy for uncoubled channels        *   
C*    Uses lsj-formalism                                             *
C*    Natural units                                                  * 
C*    Uses "constants.h"                                             *
C*    N is the number of meshpoints (Gauss-Legendre)                 *
C*********************************************************************
      subroutine bounds(nV,JJ) 
      implicit real*8 (a-h,o-z)
      include "constants.h"
      complex*16 E(N),H(N,N),bv(6)
      real*8 q(N),w(N)
      common/cmass/redM 
      write(*,*)'1623hei'
      const=1000.d0 ! half the meshpoints will be in (0,C)
      call gauleg(N,2,0.d0,const,q,w) !Gauss-Legendre (0-inf)

      do 10, j=1,N
      do 20, i=1,N
      call potential(q(i),q(j),JJ,nV,bv)

      H(i,j)= 2.d0/pi*w(j)*q(j)**2*(bv(1))*(1.d0,0.d0)
 20   continue
 10   continue
      do 30,i=1,N
      H(i,i)=H(i,i)+ q(i)**2/redM *(1.d0,0.d0)
   30 continue

      
C----Finding eigenvalues of H
      call DEVLCG(N,H,N,E)
!      write(*,*)'smalest eigenvalue:',dble(E(N))
      do 40, i=1,N
      write(*,*) E(i)
   40 continue
      
      return
      end





C*********************************************************************
C*    Calculate the bound state energy for coubled channels          *   
C*    Uses lsj-formalism                                             *
C*    Natural units                                                  * 
C*    Uses "constants.h"                                             *
C*    N is the number of meshpoints (Gauss-Legendre)                 *
C*********************************************************************
      subroutine boundcc(nV,JJ)
      implicit real*8 (a-h,o-z)
      include "constants.h"
      complex*16 E(2*N),H(2*N,2*N),bv(6)
      real*8 q(N),w(N)
      common/cmass/redM 

      const=1000.d0 ! half the meshpoints will be in (0,C)
      call gauleg(N,2,0.d0,const,q,w) !Gauss-Legendre (0-inf)

      do 10, j=1,N
      do 20, i=1,N
      call potential(q(i),q(j),JJ,nV,bv)

      H(i,j)= 2.d0/pi*w(j)*q(j)**2*(bv(3))*(1.d0,0.d0)
      H(i+N,j+N)= 2.d0/pi*w(j)*q(j)**2*(bv(4))*(1.d0,0.d0)
      H(i+N,j)= 2.d0/pi*w(j)*q(j)**2*(bv(6))*(1.d0,0.d0)
      H(i,j+N)= 2.d0/pi*w(j)*q(j)**2*(bv(5))*(1.d0,0.d0)
 20   continue
 10   continue
      do 30 i=1,N
      H(i,i)=H(i,i)+ q(i)**2/redM*(1.d0,0.d0)
      H(i+N,i+N)=H(i+N,i+N)+ q(i)**2/redM*(1.d0,0.d0)
 30   continue

C----finds the eigenvalues of H
      call DEVLCG(2*N,H,2*N,E)
!      write(*,)'Negative values are poosible bound state energies:'
      do 40,i=1,2*N
      write(*,*) E(i)
   40 continue
      
      return
      end
