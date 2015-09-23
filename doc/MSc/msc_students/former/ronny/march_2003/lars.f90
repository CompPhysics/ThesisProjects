!*****************************************************************************
!    
!     THIS PROGRAM READS IN FERMI ENERGIES FOR NUCLEONS AND HYPERONS
!     TOGETHER WITH THE TOTAL ENERGY PER BARYON.
!     IT CALCULATES THE COMPOSITION OF DENS MATTER
!     AND THE  FREE ENERGY PER BARYON.
!     
!
!            Coded May 1999 by :
!-----------------------------------------------------------------------
!  	Lars Engvik				
!  	Department of Physics		Telephone: (47) 2285 6459
!  	University of Oslo		Telefax   :	+ (47) 22856422
!  	P.O. Box 1048 Blindern	        e-mail : lars.engvik@fys.uio.no
!  	N-0316 Oslo, Norway			
!
!-----------------------------------------------------------------------
!
!     THE INPUT DATA IS READ IN THE SUBROUTINE INITIALIZE (MODULE eos)
!
!     WARNING:
!     The cascade particles are not fully implemented in the program
!    
!
!
!      The program is solves the coupled equations in the following way:
!
!       1. The baryon density is CHOSEN and relevant EOS is stored in
!          a new set of arrays with reduced dimensionalities
!       2. Initial hyperon density and fractions of sigmas and 
!          lambda are  either zero or 
!          taken from previous baryon density.
!          (this determines also the nucleon density)
!           CAUTION: The cascade particles are not fully implemented
!                    in the program. 
!          The various chemical potentials are calculated. 
!       3. Initial (new) proton fraction is set
!           and proton and neutron chemical potentials are calculated
!       4. From the chemical potentials of nucleons the lepton densities
!          are determined   
!       5. if net charge of matter is positive (negative)
!            the proton fraction is reduced (increased) 
!       6. step 3 to 5  is repeated until net charge is close to zero
!       7. If chemical potential of the sigma- is too low (high)
!          the hyperon fraction is increased (decreased)
!       8. step 3 to 7 is repeated until 
!            chemical mu_s- similar to  2*mu_n-u_p
!       9  If chemical potential of the lambda is too low (high)
!          the lambda fraction is increased (decreased)
!      10. step 3 to 9 is repeated until 
!            chemical mu_lam similar to  mu_n
!      11 Similar tests for sigma_0 and sigma_+ 
!         This has so far not been tested for finite fractions!!!!!!!
!      
!        program is easily modified to include cascade 
!
!      WARNING:
!      Program will only find solution if 
!      sigma_ is the first hyperon to appear.
!       
!
       MODULE masses
!
!      Masses of electron, muon, nucleon, lambda, sigma+0- and cascade0-
!
        DOUBLE PRECISION, PUBLIC, PARAMETER  :: wel=0.511
        DOUBLE PRECISION, PUBLIC, PARAMETER  :: wmu=105.7
        DOUBLE PRECISION, PUBLIC, PARAMETER  :: wnuc=938.926
        DOUBLE PRECISION, PUBLIC, PARAMETER  :: wlam=1116.
        DOUBLE PRECISION, PUBLIC, PARAMETER  :: wsmin=1197.
        DOUBLE PRECISION, PUBLIC, PARAMETER  :: wsnul=1193.
        DOUBLE PRECISION, PUBLIC, PARAMETER  :: wsplu=1189.
        DOUBLE PRECISION, PUBLIC, PARAMETER  :: wxinul=1315.
        DOUBLE PRECISION, PUBLIC, PARAMETER  :: wximin=1321.
       END MODULE masses
!
       MODULE constants
!
!       hbar times speed of light in MeV 
!
        DOUBLE PRECISION, PUBLIC, PARAMETER  :: hc=197.328
!
!       conversion factor for units MeV/c^2 to  kg
!
        DOUBLE PRECISION, PUBLIC, PARAMETER  :: wkg=1.78268E-30
        DOUBLE PRECISION, PUBLIC, PARAMETER  :: pi=3.141592653589793d0
        DOUBLE PRECISION, PUBLIC, PARAMETER  :: zero=0.
        DOUBLE PRECISION, PUBLIC, PARAMETER  :: one=1.
        DOUBLE PRECISION, PUBLIC, PARAMETER  :: one_=-1.
       END MODULE constants

       MODULE fermilevels
           DOUBLE PRECISION  ::pf_el,pf_mu
           DOUBLE PRECISION  ::pf_neu,pf_pro,pf_smin,pf_snul,pf_splu
           DOUBLE PRECISION  ::pf_lam,pf_ximin,pf_xinul
         END MODULE fermilevels

       MODULE densities
           DOUBLE PRECISION  ::rho_el,rho_mu,rho_nuc,rho_hyp
           DOUBLE PRECISION  ::rho_smin,rho_snul,rho_splu
           DOUBLE PRECISION  ::rho_lam,rho_ximin,rho_xinul
         END MODULE densities
       MODULE chemical
           DOUBLE PRECISION  ::chem_el,chem_mu,chem_n,chem_p
           DOUBLE PRECISION  ::chem_smin,chem_snul,chem_splu
           DOUBLE PRECISION  ::chem_lam,chem_ximin,chem_xinul
         END MODULE chemical


!      Interpolation routines mainly used in module EOS 
!
       MODULE interpolation
       CONTAINS

           SUBROUTINE interp1d(x,fx,np,x1,fx1)
             IMPLICIT NONE
             INTEGER :: npol,np
             DOUBLE PRECISION :: x(1),fx(1),x1,fx1
             IF(np==1) THEN
                fx1=fx(1)
             ELSE
                npol=MIN(4,np)
                CALL interp1d_tested(x,fx,np,npol,x1,fx1)
             ENDIF
           END SUBROUTINE interp1d


           SUBROUTINE interp1d_tested(x,fx,np,npol,x1,fx1)
!----------------------------------------------------------------------------- 
!     input arrays: x(np),fx(np)   ->function f(x)  
!           value  x1 the point 
!     output value   fx1           -> f(x1)
!     mp: x(mp)< x1 <x(mp+1) 
!     mp=1 if x1 < x(1)
!     mp=np if x1 > x(np)
!     npol= # of interpolation points
!-----------------------------------------------------------------------------
             IMPLICIT NONE
             INTEGER :: np,npol
             DOUBLE PRECISION :: x(np),fx(np),x1,fx1,df
             INTEGER :: mp,l,i
             CALL LOCATE(x,NP,x1,MP)
             L=MIN(MAX(MP-(npol-1)/2,1),NP+1-npol)
             IF(x1.LT.x(1).OR.x1.GT.x(np))THEN
                IF(mp.EQ.0.OR.mp.EQ.np)THEN
                   IF(mp.EQ.0)THEN
                      WRITE(6,*)'warning x1= ',x1,' and  x(',1,')= ',x(1)
                   ELSE
                      WRITE(6,*)'warning x1= ',x1,' and  x(',np,')= ',x(np)
                   ENDIF
                   WRITE(6,*)np,npol
                   DO i=1,np
                      WRITE(6,*)i,x(i),fx(i)
                   ENDDO
                ENDIF
             ENDIF
             CALL POLINT(x(L),fx(L),npol,x1,fx1,df)
           END SUBROUTINE interp1d_tested
         SUBROUTINE locate(xx,n,x,j)
             IMPLICIT NONE
             INTEGER :: n,j
             DOUBLE PRECISION :: xx(n),x
             INTEGER :: jl,ju,jm
             jl=0
             ju=n+1
10           IF(ju-jl.GT.1) THEN 
                jm=(ju+jl)/2
                IF((xx(n).GT.xx(1)).EQV.(x.GT.xx(jm))) THEN           
                   jl=jm
                ELSE
                   ju=jm
                ENDIF
                GOTO 10
             ENDIF
             j=jl
           END SUBROUTINE locate

           
           SUBROUTINE polint(xa,ya,n,x,y,dy)
             IMPLICIT NONE
             INTEGER :: n
             DOUBLE PRECISION :: xa(n),ya(n),x,y,dy
             INTEGER :: ns,i, m
             INTEGER, PARAMETER :: NMAX=101
             DOUBLE PRECISION :: c(NMAX),d(NMAX)
             DOUBLE PRECISION :: dif, dift, ho, hp, w, den 
             ns=1
             dif=ABS(x-xa(1))
             DO i=1,n
                dift = ABS(x-xa(i))
                IF(dift.LT.dif) THEN 
                   ns=i
                   dif=dift
                ENDIF
                c(i)=ya(i)
                d(i)=ya(i)
             ENDDO
             y=ya(ns)
             ns=ns-1
             DO m=1,n-1
                DO i=1,n-m
                   ho=xa(i)-x
                   hp=xa(i+m)-x
                   w=c(i+1)-d(i)
                   den = ho - hp
                   IF(den.EQ.0.d0) THEN
                      WRITE(6,*)'failure in polint'
                      WRITE(6,*)xa
                      WRITE(6,*)ya
                      WRITE(6,*)x
                      WRITE(6,*)y
                   ENDIF
                   den=w/den
                   d(i)=hp*den
                   c(i)=ho*den
                ENDDO
                IF (2*ns.LT.n-m) THEN
                   dy=c(ns+1)
                ELSE 
                   dy=d(ns)
                   ns=ns-1
                ENDIF
                y=y+dy
             ENDDO
           END SUBROUTINE polint

       END MODULE interpolation


       MODULE eos
!
!      The contributions from the various hyperons are grouped 
!      according to when the hyperons  are expected to appear.
!     
!
         INTEGER :: ndens,nasym,n_hyp, n_snul, n_splu, n_lam
!
!        Total energy per baryon:
         DOUBLE PRECISION, DIMENSION (15,10,5,5,5,5) :: uoa
!
!
!        Fermi energy contribution from interaction with nucleons (+kinetic e)
!        is read into arrays:
         DOUBLE PRECISION, DIMENSION (15,10,5,5,5,5) :: u_neu,u_pro,y_neu,y_pro
         DOUBLE PRECISION, DIMENSION (15,10,5,5,5,5) :: u_lam
         DOUBLE PRECISION, DIMENSION (15,10,5,5,5,5) :: u_smin,u_snul,u_splu
!
!
!        Fermi energy contribution from interaction with hyperons
!        is read into arrays:
         DOUBLE PRECISION, DIMENSION (15,10,5,5,5,5) ::  uy_neu,uy_pro,uy_smin
         DOUBLE PRECISION, DIMENSION (15,10,5,5,5,5) ::  uy_lam,uy_snul,uy_splu
         DOUBLE PRECISION, DIMENSION (10,5,5,5,5) ::  e_neu,e_pro
         DOUBLE PRECISION, DIMENSION (10,5,5,5,5) :: eoa
         DOUBLE PRECISION, DIMENSION (10,5,5,5,5) :: e_lam
         DOUBLE PRECISION, DIMENSION (10,5,5,5,5) :: e_smin,e_snul,e_splu
         DOUBLE PRECISION, DIMENSION (10,5,5,5) ::  e1_neu,e1_pro
         DOUBLE PRECISION, DIMENSION (10,5,5,5) :: eoa1
         DOUBLE PRECISION, DIMENSION (10,5,5,5) :: e1_lam
         DOUBLE PRECISION, DIMENSION (10,5,5,5) :: e1_smin,e1_snul,e1_splu
         DOUBLE PRECISION, DIMENSION (10,5,5) ::  e2_neu,e2_pro
         DOUBLE PRECISION, DIMENSION (10,5,5) :: eoa2
         DOUBLE PRECISION, DIMENSION (10,5,5) :: e2_lam
         DOUBLE PRECISION, DIMENSION (10,5,5) :: e2_smin,e2_snul,e2_splu
         DOUBLE PRECISION, DIMENSION (10,5) ::  e3_neu,e3_pro
         DOUBLE PRECISION, DIMENSION (10,5) :: eoa3
         DOUBLE PRECISION, DIMENSION (10,5) :: e3_lam
         DOUBLE PRECISION, DIMENSION (10,5) :: e3_smin,e3_snul,e3_splu
         DOUBLE PRECISION, DIMENSION (10) ::  e4_neu,e4_pro
         DOUBLE PRECISION, DIMENSION (10) :: eoa4
         DOUBLE PRECISION, DIMENSION (10) :: e4_lam
         DOUBLE PRECISION, DIMENSION (10) :: e4_smin,e4_snul,e4_splu
         DOUBLE PRECISION, DIMENSION (10) :: asym
         DOUBLE PRECISION, DIMENSION (15) ::  n_bar,e_bar
         DOUBLE PRECISION, DIMENSION (5) ::  pf_hyp,y_hyp,y1_hyp
         CONTAINS
           SUBROUTINE initialize(rhostart,drho,rhomax,dhyp)
             USE constants
             USE chemical
             USE fermilevels
             USE masses
             USE densities
             USE interpolation
             IMPLICIT NONE
             DOUBLE PRECISION, INTENT(OUT) :: rhostart,drho,rhomax,dhyp
             DOUBLE PRECISION ::  bar_dens
             DOUBLE PRECISION ::  uoa_y,uoa_n,x_p
             INTEGER :: i,j,k,l,l_, n_body
             OPEN(5,file='input.dat')          ! contain all data input for EOS
             OPEN(7,file='read.data')
             OPEN(9,file='test')
             n_snul=1
             n_splu=1
             READ(5,*)n_body                    ! test-para.  3-body eos 
!*************************************************************************
!            dhyp: is the INITIAL step size of 
!                  hyperon density and hyperon fractions
!            drho: baryon density step size
!*********************************************************************
             READ(5,*)rhostart,drho,rhomax,dhyp   
             WRITE(9,*)rhostart,drho,rhomax,dhyp
             READ(5,*)ndens,n_hyp,n_lam,nasym
             read(5,*)(y_hyp(i),i=1,n_hyp)
             WRITE(9,*)(y_hyp(i),i=1,n_hyp)
             read(5,*)(y1_hyp(i),i=1,n_lam)
             WRITE(9,*)(y1_hyp(i),i=1,n_lam)
             read(5,*)(asym(i),i=1,nasym)
             WRITE(9,*)(asym(i),i=1,nasym)
             DO i=1,ndens
                DO k=1,n_hyp
                   DO j=1,nasym
                      DO l_=1,n_lam
                         l=n_lam+1-l_
                         READ(5,*)bar_dens,pf_neu,pf_pro, &
                              & pf_smin,pf_lam, pf_snul,pf_splu ,&
                              &  u_neu(i,j,k,l,1,1),uy_neu(i,j,k,l,1,1), &
                              &  u_pro(i,j,k,l,1,1),uy_pro(i,j,k,l,1,1), &
                              & u_smin(i,j,k,l,1,1),uy_smin(i,j,k,l,1,1),&
                              & u_lam(i,j,k,l,1,1),uy_lam(i,j,k,l,1,1), &
                              & u_snul(i,j,k,l,1,1),uy_snul(i,j,k,l,1,1),&
                              & u_splu(i,j,k,l,1,1),uy_splu(i,j,k,l,1,1),&
                              & uoa_n,uoa_y,uoa(i,j,k,l,1,1)

							  u_neu(i,j,k,l,1,1)=0
							  uy_neu(i,j,k,l,1,1)=0
							  u_pro(i,j,k,l,1,1)=0
							  uy_pro(i,j,k,l,1,1)=0
                              u_smin(i,j,k,l,1,1)=0
							  uy_smin(i,j,k,l,1,1)=0
                              u_lam(i,j,k,l,1,1)=0
							  uy_lam(i,j,k,l,1,1)=0
                              u_snul(i,j,k,l,1,1)=0
							  uy_snul(i,j,k,l,1,1)=0
                              u_splu(i,j,k,l,1,1)=0
							  uy_splu(i,j,k,l,1,1)=0
                              uoa_n=0
							  uoa_y=0
							  uoa(i,j,k,l,1,1)=0





!                         write(9,101)&
!                              &bar_dens,pf_neu,pf_pro, &
!                              & pf_smin,pf_lam
!101                      FORMAT(2X,F8.4,F8.4,F8.4,F8.4,F8.4)
                         rho_hyp=bar_dens*y_hyp(k)                   
                         rho_nuc=bar_dens-rho_hyp
                         x_p=(1.-asym(j))/2.
!
!
!*************************************************************
!                        3-body interaction
!*************************************************************
                         if(n_body==3)then
                            call e_three_body(rho_nuc,x_p,chem_n,chem_p,uoa_n) 
                            u_neu(i,j,k,l,1,1)=chem_n+uy_neu(i,j,k,l,1,1)
                            u_pro(i,j,k,l,1,1)=chem_p+uy_pro(i,j,k,l,1,1)
                            u_neu(i,j,k,l,1,1)=chem_n+uy_neu(i,j,k,l,1,1)
                            u_pro(i,j,k,l,1,1)=chem_p+uy_pro(i,j,k,l,1,1)
                         else
                            u_neu(i,j,k,l,1,1)=&
                                 &u_neu(i,j,k,l,1,1)+uy_neu(i,j,k,l,1,1)
                            u_pro(i,j,k,l,1,1)=&
                                 &u_pro(i,j,k,l,1,1)+uy_pro(i,j,k,l,1,1)
                         endif
!
! **********************************************************************
!               Including YY interaction contribution
!***********************************************************************
                         uoa(i,j,k,l,1,1)=uoa_n+uoa_y
                         u_smin(i,j,k,l,1,1)=u_smin(i,j,k,l,1,1)&
                              &+uy_smin(i,j,k,l,1,1)
                         u_lam(i,j,k,l,1,1)=u_lam(i,j,k,l,1,1)&
                              &+uy_lam(i,j,k,l,1,1)
                         u_snul(i,j,k,l,1,1)=u_snul(i,j,k,l,1,1)&
                              &+uy_snul(i,j,k,l,1,1)
                         u_splu(i,j,k,l,1,1)=u_splu(i,j,k,l,1,1)&
                              &+uy_splu(i,j,k,l,1,1)
                      ENDDO
                   ENDDO
                ENDDO
                n_bar(i)=bar_dens
             ENDDO
!
!            correction of data set !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            
             u_snul(8,1,4,1,1,1)=-15.
             u_smin(8,1,4,1,1,1)=-30.
             u_smin(8,1,4,1,1,1)=-40.
             u_pro(8,1,4,1,1,1)=u_pro(8,1,3,1,1,1)+190.
             DO I=1,N_HYP
                WRITE(9,101)Y_HYP(I),U_SMIN(8,1,I,1,1,1)&
                     &,U_SMIN(8,1,I,2,1,1),U_SMIN(8,1,I,3,1,1)&
                     &,U_SMIN(8,2,I,1,1,1),U_SMIN(8,2,I,2,1,1),U_SMIN(8,2,I,3,1,1)&
                     &,U_SMIN(8,3,I,1,1,1),U_SMIN(8,3,I,2,1,1),U_SMIN(8,3,I,3,1,1)
             ENDDO
101          FORMAT(10F9.2)
             CLOSE(5)
             CLOSE(7)
             CLOSE(9)
           END SUBROUTINE initialize

!*****************************************************************
!          The EOS of the Akmal Pandharipande
!          Include a 3-body nucleonic interaction.
!          This is a parameterized version of  
!          H. Heiselberg, M. Hjorth-Jensen 
!          nucl-th/9902033 (Eq. 49)
! 
           subroutine e_three_body(n,x,chem_n,chem_p,e_n)
             IMPLICIT NONE
             DOUBLE PRECISION,  INTENT(IN) :: n,x 
             DOUBLE PRECISION,  INTENT(OUT) :: chem_n,chem_p,e_n 
             DOUBLE PRECISION , PARAMETER:: delta=0.13
             DOUBLE PRECISION , PARAMETER:: gamma=0.6
             DOUBLE PRECISION , PARAMETER:: e_0=15.8
             DOUBLE PRECISION , PARAMETER:: s_0=32.
             DOUBLE PRECISION , PARAMETER:: n_0=0.16
             DOUBLE PRECISION :: u,e_comp,e_sym,fac_u,fac_x,dedu,dedx
             u=n/n_0
             fac_u=(u-2.-delta)/(1.+delta*u)
             fac_x=(1.-2.*x)**2
             e_comp=e_0*u*fac_u
             e_sym=s_0*fac_x*u**gamma
             e_n=e_comp+e_sym
             dedu=e_0*(fac_u+(1.-fac_u*delta)/(1.+delta*u))
             dedu=dedu+e_sym*gamma/u
             dedx=-4.*s_0*u**gamma*(1.-2.*x)
             chem_n=e_n+u*dedu-x*dedx
             chem_p=chem_n+dedx
           END subroutine e_three_body

           SUBROUTINE  set_baryon_dens(barmev)
             USE masses
             USE chemical
             USE densities
             USE constants
             USE interpolation
             IMPLICIT NONE
             DOUBLE PRECISION,  INTENT(IN) :: barmev
             DOUBLE PRECISION :: ee,baryon
             INTEGER :: i,j,k,l,m,n
             baryon=barmev/hc**3
             DO i=1,n_splu
                DO j=1,n_snul
                   DO k=1,n_lam
                      DO l=1,n_hyp
                         DO m=1,nasym
                            DO n=1,ndens
                               e_bar(n)=u_neu(n,m,l,k,j,i)
                            ENDDO
                            CALL interp1d(n_bar,e_bar,ndens,baryon,ee)
                            e_neu(m,l,k,j,i)=ee
                            DO n=1,ndens
                               e_bar(n)=u_pro(n,m,l,k,j,i)
                            ENDDO
                            CALL interp1d(n_bar,e_bar,ndens,baryon,ee)
                            e_pro(m,l,k,j,i)=ee
                            DO n=1,ndens
                               e_bar(n)=u_smin(n,m,l,k,j,i)
                            ENDDO
                            CALL interp1d(n_bar,e_bar,ndens,baryon,ee)
                            e_smin(m,l,k,j,i)=ee
                            DO n=1,ndens
                               e_bar(n)=u_lam(n,m,l,k,j,i)
                            ENDDO
                            CALL interp1d(n_bar,e_bar,ndens,baryon,ee)
                            e_lam(m,l,k,j,i)=ee
                            DO n=1,ndens
                               e_bar(n)=u_snul(n,m,l,k,j,i)
                            ENDDO
                            CALL interp1d(n_bar,e_bar,ndens,baryon,ee)
                            e_snul(m,l,k,j,i)=ee
                            DO n=1,ndens
                               e_bar(n)=u_splu(n,m,l,k,j,i)
                            ENDDO
                            CALL interp1d(n_bar,e_bar,ndens,baryon,ee)
                            e_splu(m,l,k,j,i)=ee
                            DO n=1,ndens
                               e_bar(n)=uoa(n,m,l,k,j,i)
                            ENDDO
                            CALL interp1d(n_bar,e_bar,ndens,baryon,ee)
                            eoa(m,l,k,j,i)=ee
                          ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
           END SUBROUTINE set_baryon_dens

           SUBROUTINE  set_chem_xinul
             USE masses
             USE chemical
             USE densities
             USE constants
             USE fermilevels
             IMPLICIT NONE
             chem_xinul=e_ferm(wxinul,rho_xinul)
           END SUBROUTINE set_chem_xinul

           SUBROUTINE  set_chem_ximin
             USE masses
             USE chemical
             USE densities
             USE constants
             USE fermilevels
             IMPLICIT NONE
             chem_ximin=e_ferm(wximin,rho_ximin)
           END SUBROUTINE set_chem_ximin
 
           

!***********************************************
!         Store the EOS data with fixed
!         baryon density and sigma+ fraction into new arrays
!
!         Output:
!         The chemical potential of Sigma_+ 
!         The average energy per baryon (av_eoa)
!***********************************************

           SUBROUTINE  fix_dens_splu(xasym,av_eoa)
             USE masses
             USE chemical
             USE densities
             USE constants
             USE interpolation
             USE fermilevels
             IMPLICIT NONE
             DOUBLE PRECISION,  INTENT(IN) :: xasym
             DOUBLE PRECISION,  INTENT(OUT) :: av_eoa
             DOUBLE PRECISION :: ee,y_splu,y_snul,y_lam,hypbar
             INTEGER :: j,k,l,m,n
             pf_splu=(rho_splu*3*pi**2)**.333333333/hc
             hypbar=rho_hyp/(rho_nuc+rho_hyp)
             IF(rho_hyp<=0.)THEN
                y_splu=0.
                y_snul=0.
                y_lam=0.
             ELSE
                y_splu=rho_splu/rho_hyp
                y_snul=rho_snul/rho_hyp
                y_lam=rho_lam/rho_hyp
             ENDIF
             DO j=1,n_snul
                DO k=1,n_lam
                   DO l=1,n_hyp
                      DO m=1,nasym
                         DO n=1,n_splu
                            e_bar(n)=e_neu(m,l,k,j,n)
                         ENDDO
                         CALL interp1d(y1_hyp,e_bar,n_splu,y_splu,ee)
                         e1_neu(m,l,k,j)=ee
                         DO n=1,n_splu
                            e_bar(n)=e_pro(m,l,k,j,n)
                         ENDDO
                         CALL interp1d(y1_hyp,e_bar,n_splu,y_splu,ee)
                         e1_pro(m,l,k,j)=ee
                         DO n=1,n_splu
                            e_bar(n)=e_smin(m,l,k,j,n)
                         ENDDO
                         CALL interp1d(y1_hyp,e_bar,n_splu,y_splu,ee)
                         e1_smin(m,l,k,j)=ee
                         DO n=1,n_splu
                            e_bar(n)=e_lam(m,l,k,j,n)
                         ENDDO
                         CALL interp1d(y1_hyp,e_bar,n_splu,y_splu,ee)
                         e1_lam(m,l,k,j)=ee
                         DO n=1,n_splu
                            e_bar(n)=e_snul(m,l,k,j,n)
                         ENDDO
                         CALL interp1d(y1_hyp,e_bar,n_splu,y_splu,ee)
                         e1_snul(m,l,k,j)=ee
                         DO n=1,n_splu
                            e_bar(n)=e_splu(m,l,k,j,n)
                         ENDDO
                         CALL interp1d(y1_hyp,e_bar,n_splu,y_splu,ee)
                         e1_splu(m,l,k,j)=ee
                         DO n=1,n_splu
                            e_bar(n)=eoa(m,l,k,j,n)
                         ENDDO
                         CALL interp1d(y1_hyp,e_bar,n_splu,y_splu,ee)
                         eoa1(m,l,k,j)=ee
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
!  ********************          sigma_null        ***************
             DO k=1,n_lam
                DO l=1,n_hyp
                   DO m=1,nasym
                      DO n=1,n_snul
                         e_bar(n)=e1_splu(m,l,k,n)
                      ENDDO
                      CALL interp1d(y1_hyp,e_bar,n_snul,y_snul,ee)
                      e2_splu(m,l,k)=ee
                      DO n=1,n_snul
                         e_bar(n)=eoa1(m,l,k,n)
                      ENDDO
                      CALL interp1d(y1_hyp,e_bar,n_snul,y_snul,ee)
                      eoa2(m,l,k)=ee
                   ENDDO
                ENDDO
             ENDDO
!  ********************          lambda        ***************
             DO l=1,n_hyp
                DO m=1,nasym
                   DO n=1,n_lam
                      e_bar(n)=e2_splu(m,l,n)
                   ENDDO
                   CALL interp1d(y1_hyp,e_bar,n_lam,y_lam,ee)
                   e3_splu(m,l)=ee
                   DO n=1,n_lam
                      e_bar(n)=eoa2(m,l,n)
                   ENDDO
                   CALL interp1d(y1_hyp,e_bar,n_lam,y_lam,ee)
                   eoa3(m,l)=ee
                ENDDO
             ENDDO
! *************** hyperon densities ***********************
             DO m=1,nasym
                DO n=1,n_hyp
                   e_bar(n)=e3_splu(m,n)
                ENDDO
                CALL interp1d(y_hyp,e_bar,n_hyp,hypbar,ee)
                e4_splu(m)=ee
                DO n=1,n_hyp
                   e_bar(n)=eoa3(m,n)
                ENDDO
                CALL interp1d(y_hyp,e_bar,n_hyp,hypbar,ee)
                eoa4(m)=ee
             ENDDO
             CALL interp1d(asym,e4_splu,nasym,xasym,ee)
             chem_splu=ee+wsplu
             CALL interp1d(asym,eoa4,nasym,xasym,ee)
             av_eoa=ee
           END SUBROUTINE fix_dens_splu

!***********************************************
!         Store the EOS data with fixed
!         - baryon density 
!         - sigma+ and sigma_0 fractions 
!         into new arrays
!
!         Output:
!         The chemical potential of Sigma_0 
!***********************************************
           SUBROUTINE  fix_dens_snul(xasym)
             USE masses
             USE chemical
             USE densities
             USE constants
             USE interpolation
             USE fermilevels
             IMPLICIT NONE
             DOUBLE PRECISION,  INTENT(IN) :: xasym
             DOUBLE PRECISION :: ee,y_snul,y_lam,hypbar
             INTEGER :: k,l,m,n
             hypbar=rho_hyp/(rho_nuc+rho_hyp)
             IF(rho_hyp<=0.)THEN
                y_snul=0.
                y_lam=0.
             ELSE
                y_snul=rho_snul/rho_hyp
                y_lam=rho_lam/rho_hyp
             ENDIF
             pf_snul=(rho_snul*3*pi**2)**.333333333/hc
             DO k=1,n_lam
                DO l=1,n_hyp
                   DO m=1,nasym
                      DO n=1,n_snul
                         e_bar(n)=e1_neu(m,l,k,n)
                      ENDDO
                      CALL interp1d(y1_hyp,e_bar,n_snul,y_snul,ee)
                      e2_neu(m,l,k)=ee
                      DO n=1,n_snul
                         e_bar(n)=e1_pro(m,l,k,n)
                      ENDDO
                      CALL interp1d(y1_hyp,e_bar,n_snul,y_snul,ee)
                      e2_pro(m,l,k)=ee
                      DO n=1,n_snul
                         e_bar(n)=e1_smin(m,l,k,n)
                      ENDDO
                      CALL interp1d(y1_hyp,e_bar,n_snul,y_snul,ee)
                      e2_smin(m,l,k)=ee
                      DO n=1,n_snul
                         e_bar(n)=e1_lam(m,l,k,n)
                      ENDDO
                      CALL interp1d(y1_hyp,e_bar,n_snul,y_snul,ee)
                      e2_lam(m,l,k)=ee
                      DO n=1,n_snul
                         e_bar(n)=e1_snul(m,l,k,n)
                      ENDDO
                      CALL interp1d(y1_hyp,e_bar,n_snul,y_snul,ee)
                      e2_snul(m,l,k)=ee
                   ENDDO
                ENDDO
             ENDDO
             pf_lam=(rho_lam*3*pi**2)**.333333333/hc
             DO l=1,n_hyp
                DO m=1,nasym
                   DO n=1,n_lam
                      e_bar(n)=e2_snul(m,l,n)
                   ENDDO
                   CALL interp1d(y1_hyp,e_bar,n_lam,y_lam,ee)
                   e3_snul(m,l)=ee
                ENDDO
             ENDDO
             pf_smin=(rho_smin*3*pi**2)**.333333333/hc
             DO m=1,nasym
                DO n=1,n_hyp
                   e_bar(n)=e3_snul(m,n)
                ENDDO
                CALL interp1d(y_hyp,e_bar,n_hyp,hypbar,ee)
                e4_snul(m)=ee
             ENDDO
             CALL interp1d(asym,e4_snul,nasym,xasym,ee)
             chem_snul=ee+wsnul
           END SUBROUTINE fix_dens_snul

!***********************************************
!         Store the EOS data with fixed
!         - baryon density 
!         - sigma+ sigma_0 and lambda fractions 
!         into new arrays
!
!         Output:
!         The chemical potential of lambda 
!***********************************************
           SUBROUTINE  fix_dens_lam(xasym)
             USE masses
             USE chemical
             USE densities
             USE constants
             USE interpolation
             USE fermilevels
             IMPLICIT NONE
             DOUBLE PRECISION,  INTENT(IN) :: xasym
             DOUBLE PRECISION :: ee,y_lam,hypbar
             INTEGER :: l,m,n
             hypbar=rho_hyp/(rho_nuc+rho_hyp)
             IF(rho_hyp<=0.)THEN
                y_lam=0.
             ELSE
                y_lam=rho_lam/rho_hyp
             ENDIF
             pf_lam=(rho_lam*3*pi**2)**.333333333/hc
             DO l=1,n_hyp
                DO m=1,nasym
                   DO n=1,n_lam
                      e_bar(n)=e2_neu(m,l,n)
                   ENDDO
                   CALL interp1d(y1_hyp,e_bar,n_lam,y_lam,ee)
                   e3_neu(m,l)=ee
                   DO n=1,n_lam
                      e_bar(n)=e2_pro(m,l,n)
                   ENDDO
                   CALL interp1d(y1_hyp,e_bar,n_lam,y_lam,ee)
                   e3_pro(m,l)=ee
                   DO n=1,n_lam
                      e_bar(n)=e2_smin(m,l,n)
                   ENDDO
                   CALL interp1d(y1_hyp,e_bar,n_lam,y_lam,ee)
                   e3_smin(m,l)=ee
                   DO n=1,n_hyp
                      e_bar(n)=e2_lam(m,l,n)
                   ENDDO
                   CALL interp1d(y1_hyp,e_bar,n_lam,y_lam,ee)
                   e3_lam(m,l)=ee
                   DO n=1,n_lam
                      e_bar(n)=eoa2(m,l,n)
                   ENDDO
                   CALL interp1d(y1_hyp,e_bar,n_lam,y_lam,ee)
                   eoa3(m,l)=ee
                ENDDO
             ENDDO
             pf_smin=(rho_smin*3*pi**2)**.333333333/hc
             DO m=1,nasym
                DO n=1,n_hyp
                   e_bar(n)=e3_lam(m,n)
                ENDDO
                CALL interp1d(y_hyp,e_bar,n_hyp,hypbar,ee)
                e4_smin(m)=ee
             ENDDO
             CALL interp1d(asym,e4_smin,nasym,xasym,ee)
             chem_lam=ee+wlam
           END SUBROUTINE fix_dens_lam

!***********************************************
!         Store the EOS data with fixed
!         - baryon density 
!         - sigma+, sigma_0, lambda fractions
!         - hyperon density  
!         into new arrays
!
!         Output:
!         The chemical potential of Sigma_- 
!***********************************************
           SUBROUTINE  fix_y_hyp(xasym)
             USE masses
             USE chemical
             USE densities
             USE constants
             USE interpolation
             USE fermilevels
             IMPLICIT NONE
             DOUBLE PRECISION,  INTENT(IN) :: xasym
             DOUBLE PRECISION :: ee,hypbar
             INTEGER :: m,n
             rho_hyp=rho_smin+rho_lam+rho_splu+rho_snul !+rho_ximin+rho_xinul
             hypbar=rho_hyp/(rho_nuc+rho_hyp)
             DO m=1,nasym
                DO n=1,n_hyp
                   e_bar(n)=e3_neu(m,n)
                ENDDO
                CALL interp1d(y_hyp,e_bar,n_hyp,hypbar,ee)
                e4_neu(m)=ee
                DO n=1,n_hyp
                   e_bar(n)=e3_pro(m,n)
                ENDDO
                CALL interp1d(y_hyp,e_bar,n_hyp,hypbar,ee)
                e4_pro(m)=ee
                DO n=1,n_hyp
                   e_bar(n)=e3_smin(m,n)
                ENDDO
                CALL interp1d(y_hyp,e_bar,n_hyp,hypbar,ee)
                e4_smin(m)=ee
                DO n=1,n_hyp
                   e_bar(n)=eoa3(m,n)
                ENDDO
                CALL interp1d(y_hyp,e_bar,n_hyp,hypbar,ee)
                eoa4(m)=ee
             ENDDO
             CALL interp1d(asym,e4_smin,nasym,xasym,ee)
             chem_smin=ee+wsmin     ! e_ferm(wsmin,rho_smin) 
          END SUBROUTINE fix_y_hyp


!***********************************************
!         Fixed
!         - baryon density 
!         - sigma+, sigma_0, lambda fractions
!         - hyperon density  
!         - asymmetry parameter: xasym= (n_n-n_p)/(n_B-n_Y))
!
!         Output:
!         The chemical potential of neutrons and protons
!***********************************************
           SUBROUTINE chem_nuc(xasym)
             USE masses
             USE chemical
             USE densities
             USE interpolation
             USE constants
             USE fermilevels
             IMPLICIT NONE
             DOUBLE PRECISION :: xasym
             DOUBLE PRECISION  :: ee,rho_n,rho_p
             rho_n=rho_nuc*(1.0+xasym)/2.
             rho_p=rho_nuc*(1.0-xasym)/2.
             pf_neu=(3*pi**2*rho_n)**.333333333/hc
             pf_pro=(3*pi**2*rho_n)**.333333333/hc
             CALL interp1d(asym,e4_neu,nasym,xasym,ee)
             chem_n=ee+wnuc
             CALL interp1d(asym,e4_pro,nasym,xasym,ee)
             chem_p=ee+wnuc
           END SUBROUTINE chem_nuc

           FUNCTION e_ferm(mass,rho)
             USE constants
             IMPLICIT NONE
             DOUBLE PRECISION,  INTENT(IN) :: mass,rho
             DOUBLE PRECISION :: e_ferm
             DOUBLE PRECISION :: pf
             pf=(3.*pi**2*rho)**0.33333333333
             IF (mass>200.)THEN
                e_ferm=mass+pf**2/2./mass
             ELSE
                e_ferm=(mass**2+pf**2)**.5
             ENDIF
           END FUNCTION e_ferm


           SUBROUTINE  write_contributions(barmev,xasym)
             USE masses
             USE chemical
             USE densities
             USE constants
             USE interpolation
             IMPLICIT NONE
             DOUBLE PRECISION,  INTENT(IN) :: barmev,xasym
             DOUBLE PRECISION :: ee,baryon,rho_n,rho_p
             DOUBLE PRECISION ::  u1y_neu,u1y_pro,u1y_smin
             DOUBLE PRECISION ::  u1y_lam,u1y_snul,u1y_splu
             DOUBLE PRECISION ::  u1n_neu,u1n_pro,u1n_smin
             DOUBLE PRECISION ::  u1n_lam,u1n_snul,u1n_splu
             INTEGER :: i,j,k,l,m,n
             baryon=barmev/hc**3
             rho_n=rho_nuc*(1.+xasym)/2.
             rho_p=rho_nuc*(1.-xasym)/2.
             ee=e_ferm(wnuc,rho_n)             
             u1n_neu=chem_n -ee
             ee=e_ferm(wnuc,rho_p)             
             u1n_pro=chem_p-ee
             ee=e_ferm(wsmin,rho_smin)             
             u1n_smin=chem_smin-ee
             ee=e_ferm(wlam,rho_lam)             
             u1n_lam=chem_lam-ee
             ee=e_ferm(wsnul,rho_snul)             
             u1n_snul=chem_snul-ee
             ee=e_ferm(wsplu,rho_splu)             
             u1n_splu=chem_splu-ee
             DO i=1,n_splu
                DO j=1,n_snul
                   DO k=1,n_lam
                      DO l=1,n_hyp
                         DO m=1,nasym
                            DO n=1,ndens
                               e_bar(n)=uy_neu(n,m,l,k,j,i)
                            ENDDO
                            CALL interp1d(n_bar,e_bar,ndens,baryon,ee)
                            e_neu(m,l,k,j,i)=ee
                            DO n=1,ndens
                               e_bar(n)=uy_pro(n,m,l,k,j,i)
                            ENDDO
                            CALL interp1d(n_bar,e_bar,ndens,baryon,ee)
                            e_pro(m,l,k,j,i)=ee
                            DO n=1,ndens
                               e_bar(n)=uy_smin(n,m,l,k,j,i)
                            ENDDO
                            CALL interp1d(n_bar,e_bar,ndens,baryon,ee)
                            e_smin(m,l,k,j,i)=ee
                            DO n=1,ndens
                               e_bar(n)=uy_lam(n,m,l,k,j,i)
                            ENDDO
                            CALL interp1d(n_bar,e_bar,ndens,baryon,ee)
                            e_lam(m,l,k,j,i)=ee
                            DO n=1,ndens
                               e_bar(n)=uy_snul(n,m,l,k,j,i)
                            ENDDO
                            CALL interp1d(n_bar,e_bar,ndens,baryon,ee)
                            e_snul(m,l,k,j,i)=ee
                            DO n=1,ndens
                               e_bar(n)=uy_splu(n,m,l,k,j,i)
                            ENDDO
                            CALL interp1d(n_bar,e_bar,ndens,baryon,ee)
                            e_splu(m,l,k,j,i)=ee
                          ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
             call fix_dens_splu(xasym,ee)
             call fix_dens_snul(xasym)
             call fix_dens_lam(xasym)
             call fix_y_hyp(xasym)
             call chem_nuc(xasym)
             u1y_neu=chem_n -wnuc
             u1y_pro=chem_p -wnuc
             u1y_smin=chem_smin -wsmin
             u1y_lam=chem_lam   -wlam
             u1y_snul=chem_snul -wsnul
             u1y_splu=chem_splu -wsplu
             write(17,102)baryon,u1n_neu-u1y_neu,u1n_pro-u1y_pro&
                  &,u1n_smin-u1y_smin,u1n_lam-u1y_lam,u1n_snul-u1y_snul&
                  &,u1n_splu-u1y_splu
             write(18,102)baryon,u1y_neu,u1y_pro&
                  &,u1y_smin,u1y_lam,u1y_snul&
                  &,u1y_splu
102          FORMAT(F9.3,11F9.2)

           END SUBROUTINE write_contributions

       END MODULE eos


       PROGRAM  hyperons
         USE constants
         USE densities
         USE eos
         IMPLICIT NONE
         DOUBLE PRECISION  ::rho,drho,rhomax,dhyp,xasym
         OPEN(6,file='densplus.dat')
		 OPEN(21,file='dens.dat')
         OPEN(4,file='chem.dat')
         OPEN(16,file='eoa.dat')
         OPEN(17,file='u_nuc.dat')
         OPEN(18,file='u_y.dat')
         CALL initialize(rho,drho,rhomax,dhyp)
         rho_hyp=0.
         rho=rho*hc**3     
         drho=drho*hc**3     
         dhyp=dhyp*hc**3     
         rho_nuc=rho
         xasym=0.9
         ENDLESS: DO 
            IF(rho/hc**3 > rhomax) THEN
               EXIT ENDLESS
            ENDIF
            CALL set_baryon_dens(rho)
            CALL composition(rho,dhyp,xasym)
            rho_nuc=rho_nuc+drho
            rho=rho+drho
         ENDDO ENDLESS
         CLOSE(6)
         CLOSE(4)
         CLOSE(16)
         CLOSE(17)
         CLOSE(18)
		 CLOSE(21)

       CONTAINS
         SUBROUTINE composition(rho,dhyp,xasym)
           USE masses
           USE chemical
           USE densities
           USE fermilevels
           USE eos
           IMPLICIT NONE
           DOUBLE PRECISION  ::rho,dhyp,xasym,step,xacc
           DOUBLE PRECISION  ::av_eoa,free_e
           DOUBLE PRECISION  ::rho_n,rho_p, f,e_el,e_mu,ed_el,ed_mu
           INTEGER :: i
           xacc=0.00001d0
           step=dhyp
           DO i=1,5
              CALL with_xi_null(step,xacc,xasym,av_eoa)
              step=step/10.
!*******************************************************
!   Improving the accuracy for each loop
!*******************************************************
           ENDDO
           rho_n=rho_nuc*(1.+xasym)/2.
           rho_p=rho_nuc*(1.-xasym)/2.
           f=hc**3
           if(rho_el==0.)then 
              e_el=0.
              e_mu=0.
           elseif(rho_mu==0.)then
              ed_el=e_dens(wel,pf_el)
              e_el=ed_el/rho_el
              e_mu=.0
           else
              ed_el=e_dens(wel,pf_el)
              e_el=ed_el/rho_el
              ed_mu=e_dens(wmu,pf_mu)
              e_mu=ed_mu/rho_mu
           endif
           av_eoa=av_eoa+(rho_nuc*wnuc)/rho
           av_eoa=av_eoa+(rho_smin*wsmin)/rho
           av_eoa=av_eoa+(rho_lam*wlam)/rho
           av_eoa=av_eoa+(rho_snul*wsnul)/rho
           av_eoa=av_eoa+(rho_splu*wsplu)/rho
           av_eoa=av_eoa-wnuc
           free_e=av_eoa+(ed_el+ed_mu)/rho
           WRITE(16,102)rho/f,e_el,e_mu,av_eoa,free_e
           WRITE(6,100)rho/f,rho_n/f,rho_p/f,rho_el/f,rho_mu/f,rho_smin/f,rho_lam/f,rho_snul/f
		   WRITE(21,100)rho/f,rho_n/f,rho_p/f,rho_el/f,rho_mu/f,rho_smin/f,rho_lam/f,rho_snul/f,0.0d0
           WRITE(4,101)rho/f,chem_n,chem_p,chem_el,chem_smin, chem_lam,chem_snul,chem_splu
           call write_contributions(rho,xasym)
100        FORMAT(11F8.4)
101        FORMAT(F9.3,10F9.1)
102        FORMAT(F9.3,11F9.2)
         END SUBROUTINE composition
         

!***************************************************************************
!       Warning:
!       Cascade particles are not fully implemented in the program
!       The results are valid only if rho_xinul =0 
!***************************************************************************

         SUBROUTINE with_xi_null(step,xacc,xasym,av_eoa)
           USE chemical
           USE densities
           USE constants
           IMPLICIT NONE
           DOUBLE PRECISION,  INTENT(OUT) :: xasym,av_eoa 
           DOUBLE PRECISION  :: step,xacc,sign
           CALL with_xi_minus(step,xacc,xasym,av_eoa)
           CALL set_chem_xinul
           IF(chem_n>chem_xinul.OR.rho_xinul>0)THEN
              IF(chem_n<chem_xinul) THEN
                 sign=-1.
              ELSE
                 sign=1.
              ENDIF
              ENDLESS: DO 
                 rho_smin=rho_smin-step*sign
                 rho_xinul=rho_xinul+step*sign
                 IF(rho_xinul <= 0.) THEN
                    rho_smin=rho_smin+rho_xinul
                    rho_xinul=0.
                    CALL with_xi_minus(step,xacc,xasym,av_eoa)
                    CALL set_chem_xinul
                    EXIT ENDLESS
                 ENDIF
                 CALL with_xi_minus(step,xacc,xasym,av_eoa)
                 CALL set_chem_xinul
                 IF(sign*(chem_n-chem_xinul)<0) THEN           
                    EXIT ENDLESS
                 ENDIF
              ENDDO ENDLESS
           ENDIF
         END SUBROUTINE with_xi_null

!***************************************************************************
!       Warning:
!       Cascade particles are not fully implemented in the program
!       The results are valid only if rho_ximin =0 
!***************************************************************************
         SUBROUTINE with_xi_minus(step,xacc,xasym,av_eoa)
           USE chemical
           USE densities
           USE constants
           IMPLICIT NONE
           DOUBLE PRECISION,  INTENT(OUT) :: xasym,av_eoa 
           DOUBLE PRECISION  :: step,xacc,sign
           CALL with_sigma_plus(step,xacc,xasym,av_eoa)
           CALL set_chem_ximin
           IF(chem_n+chem_el>chem_ximin.OR.rho_ximin>0)THEN
              IF(chem_n+chem_el<chem_ximin) THEN
                 sign=-1.
              ELSE
                 sign=1.
              ENDIF
              ENDLESS: DO 
                 rho_smin=rho_smin-step*sign
                 rho_ximin=rho_ximin+step*sign
                 IF(rho_ximin <= 0.) THEN
                    rho_smin=rho_smin+rho_ximin
                    rho_ximin=0.
                    CALL with_sigma_plus(step,xacc,xasym,av_eoa)
                    CALL set_chem_ximin
                    EXIT ENDLESS
                 ENDIF
                 CALL with_sigma_plus(step,xacc,xasym,av_eoa)
                 CALL set_chem_ximin
                 IF(sign*(chem_n+chem_el-chem_ximin)<0) THEN           
                    EXIT ENDLESS
                 ENDIF
              ENDDO ENDLESS
           ENDIF
         END SUBROUTINE with_xi_minus


         SUBROUTINE with_sigma_plus(step,xacc,xasym,av_eoa)
!***************************************************************************
!       Determine the Sigma_+ density  after lambda, Sigma_0 densities 
!       and  total hyperon density have been fixed.
!***************************************************************************
           USE chemical
           USE densities
           USE constants
           IMPLICIT NONE
           DOUBLE PRECISION,  INTENT(OUT) :: xasym,av_eoa 
           DOUBLE PRECISION  :: step,xacc,sign
           CALL fix_dens_splu(xasym,av_eoa)
           CALL with_sigma_null(step,xacc,xasym)
           IF(chem_p>chem_splu.OR.rho_splu>0)THEN
              IF(chem_p<chem_splu) THEN
                 sign=-1.
              ELSE
                 sign=1.
              ENDIF
              ENDLESS: DO 
                 rho_smin=rho_smin-step*sign
                 rho_splu=rho_splu+step*sign
                 IF(rho_splu <= 0.) THEN
                    rho_smin=rho_smin+rho_splu
                    rho_splu=0.
                    CALL fix_dens_splu(xasym,av_eoa)
                    CALL with_sigma_null(step,xacc,xasym)
                    EXIT ENDLESS
                 ENDIF
                 CALL fix_dens_splu(xasym,av_eoa)
                 CALL with_sigma_null(step,xacc,xasym)
                 IF(sign*(chem_p-chem_splu)<0) THEN           
                    EXIT ENDLESS
                 ENDIF
              ENDDO ENDLESS
           ENDIF
         END SUBROUTINE with_sigma_plus


         SUBROUTINE with_sigma_null(step,xacc,xasym)
!***************************************************************************
!       Determine the Sigma_0 density  after the lambda density
!       and total hyperon density have been fixed.
!***************************************************************************
           USE chemical
           USE densities
           USE constants
           IMPLICIT NONE
           DOUBLE PRECISION  :: step,xacc,xasym,sign
           CALL fix_dens_snul(xasym)
           CALL with_lambda(step,xacc,xasym)
           IF(chem_n>chem_snul.OR.rho_snul>0)THEN
              IF(chem_n<chem_snul) THEN
                 sign=-1.
              ELSE
                 sign=1.
              ENDIF
              ENDLESS: DO 
                 rho_smin=rho_smin-step*sign
                 rho_snul=rho_snul+step*sign
                 IF(rho_snul <= 0.) THEN
                    rho_smin=rho_smin+rho_snul
                    rho_snul=0.
                    CALL fix_dens_snul(xasym)
                    CALL with_lambda(step,xacc,xasym)
                    EXIT ENDLESS
                 ENDIF
                 CALL fix_dens_snul(xasym)
                 CALL with_lambda(step,xacc,xasym)
                 IF(sign*(chem_n-chem_snul)<0) THEN           
                    EXIT ENDLESS
                 ENDIF
              ENDDO ENDLESS
           ENDIF
         END SUBROUTINE with_sigma_null


         SUBROUTINE with_lambda(step,xacc,xasym)
!***************************************************************************
!       Determine the lambda density  after 
!       hyperon density have been fixed.
!***************************************************************************
           USE chemical
           USE densities
           USE constants
           IMPLICIT NONE
           DOUBLE PRECISION  :: step,xacc,xasym,sign
           CALL fix_dens_lam(xasym)
           CALL with_sigma_minus(step,xacc,xasym)
           IF(chem_n>chem_lam.OR.rho_lam>0)THEN
              IF(chem_n<chem_lam) THEN
                 sign=-1.
              ELSE
                 sign=1.
              ENDIF
              ENDLESS: DO 
                 rho_smin=rho_smin-step*sign
                 rho_lam=rho_lam+step*sign
                 IF(rho_lam <= 0.) THEN
                    rho_smin=rho_smin+rho_lam
                    rho_lam=0.
                    CALL fix_dens_lam(xasym)
                    CALL with_sigma_minus(step,xacc,xasym)
                    EXIT ENDLESS
                 ENDIF
                 CALL fix_dens_lam(xasym)
                 CALL with_sigma_minus(step,xacc,xasym)
                 IF(sign*(chem_n-chem_lam)<0) THEN           
                    EXIT ENDLESS
                 ENDIF
              ENDDO ENDLESS
           ENDIF
         END SUBROUTINE with_lambda

         SUBROUTINE with_sigma_minus(step,xacc,xasym)
!***************************************************************************
!       Determine the sigma_ density  after 
!       proton fraction have been fixed.
!***************************************************************************
           USE chemical
           USE densities
           USE constants
           IMPLICIT NONE
           DOUBLE PRECISION  :: step,xacc,xasym,sign
           CALL fix_y_hyp(xasym)
           CALL root(zero,one,xacc,xasym)
           IF(chem_n+chem_el>chem_smin.OR.rho_smin>0)THEN
              IF(chem_n+chem_el<chem_smin) THEN
                 sign=-1.
              ELSE
                 sign=1.
              ENDIF
              ENDLESS: DO 
                 rho_nuc=rho_nuc-step*sign
                 rho_smin=rho_smin+step*sign
                 rho_hyp=rho_hyp+step*sign
                 IF(rho_smin <= 0.) THEN
                    rho_nuc=rho_nuc+rho_smin
                    rho_smin=0.
                    rho_hyp=rho_hyp-rho_smin
                    CALL fix_y_hyp(xasym)
                    CALL root(zero,one,xacc,xasym)
                    EXIT ENDLESS
                 ENDIF
                 CALL fix_y_hyp(xasym)
                 CALL root(zero,one,xacc,xasym)
                 IF(sign*(chem_n+chem_el-chem_smin)<0) THEN           
                    EXIT ENDLESS
                 ENDIF
              ENDDO ENDLESS
           ENDIF
         END SUBROUTINE with_sigma_minus

         SUBROUTINE rootx(x1,x2,xacc,rtbis)
           IMPLICIT NONE
           DOUBLE PRECISION,  INTENT(OUT) :: rtbis
           DOUBLE PRECISION,  INTENT(IN) ::  x1,x2,xacc
           DOUBLE PRECISION ::  dx
           DX=beta_eq(x2)
           DX=x1+x2+xacc
           rtbis=x2
         END SUBROUTINE rootx


!**************************************************************************
!        Find the asymmetry parameter (rtbis) 
!        where matter has zero net charge (beta_eq=0)
!        rtbis fixes the proton and neutron densities
!**************************************************************************
         SUBROUTINE root(x1,x2,xacc,rtbis)
           IMPLICIT NONE
           INTEGER, PARAMETER  :: JMAX=100
           INTEGER :: j
           DOUBLE PRECISION,  INTENT(OUT) :: rtbis
           DOUBLE PRECISION,  INTENT(IN) ::  x1,x2,xacc
           DOUBLE PRECISION ::  dx,f,fmid,xmid
           fmid=beta_eq(x2)
           f=beta_eq(x1)
           IF(f*fmid.GE.0.d0) WRITE(6,*) 'root must be bracketed in rtbis'
           IF(f.LT.0.d0)THEN
              rtbis=x1
              dx=x2-x1
           ELSE
              rtbis=x2
              dx=x1-x2
           ENDIF
           j=1
           endless: DO 
              dx=dx*.5d0
              xmid=rtbis+dx
              fmid=beta_eq(xmid)
              IF(fmid.LE.0.d0)rtbis=xmid
              IF(ABS(dx).LT.xacc .OR. fmid.EQ.0.d0)THEN 
                 EXIT endless
              ELSEIF(j >= jmax)THEN
                 WRITE(6,*) 'too many bisections in rtbis'
                 EXIT endless
              ENDIF
              j=j+1   
           ENDDO endless
           IF(j==jmax) WRITE(6,*) 'too many bisections in rtbis'
         END SUBROUTINE root

         FUNCTION beta_eq(x)
!******************************************************************
!        Determine the total charge excess in hyperonic matter.
!
!        Output:
!        Chemical potential of electrons
!        lepton  densities
!******************************************************************
           USE masses
           USE densities
           USE chemical
           USE fermilevels
           IMPLICIT NONE
           DOUBLE PRECISION,  INTENT(IN) :: x
           DOUBLE PRECISION :: beta_eq
           DOUBLE PRECISION :: uh2,wm2,we2
           DOUBLE PRECISION :: rho_n,rho_p
           rho_n=rho_nuc*(1.+x)/2.
           rho_p=rho_nuc*(1.-x)/2.
           CALL chem_nuc(x)
           chem_el=chem_n-chem_p
           uh2=chem_el**2
           wm2=wmu**2
           we2=wel**2
           pf_mu=0.d0
           IF(uh2.GT.wm2)pf_mu=SQRT(uh2-wm2)*chem_el/ABS(chem_el) !
           pf_el=0.d0
           IF(uh2.GT.we2)pf_el=SQRT(uh2-we2)*chem_el/ABS(chem_el)
!                                                       ^
!                mu_e <0 => positrons <=> rho_e<0      _|  
!
           rho_el=pf_el**3/(3.d0*pi**2)
           rho_mu=pf_mu**3/(3.d0*pi**2)
!           beta_eq=rho_splu+rho_p-rho_el-rho_mu-rho_smin-rho_ximin
           beta_eq=rho_splu+rho_p-rho_el-rho_mu-rho_smin-rho_ximin
           RETURN
         END FUNCTION beta_eq


           FUNCTION fermilevel(mass,chem)
             IMPLICIT NONE
             DOUBLE PRECISION,  INTENT(IN) :: mass,chem
             DOUBLE PRECISION :: fermilevel
             IF (chem < mass )THEN
               fermilevel=0.
             ELSE
                fermilevel=(chem**2-mass**2)**0.5
             ENDIF
           END FUNCTION fermilevel

           
           FUNCTION e_dens(mass,pf)
             IMPLICIT NONE
             DOUBLE PRECISION,  INTENT(IN) :: mass,pf
             DOUBLE PRECISION :: e_dens
             DOUBLE PRECISION :: t,t2,term1,term2
             IF (pf<= 0.)THEN
                e_dens=0.
             ELSE
                t=pf/mass 
                t2=t*t
                term2=LOG(t+(t2+1.)**.5)       
                term1=(2.*t2+1.)*T*(t2+1.)**.5
                e_dens=mass**4*(term1-term2)/(8.*pi**2)
             ENDIF
           END FUNCTION e_dens


           
           
       END PROGRAM hyperons


