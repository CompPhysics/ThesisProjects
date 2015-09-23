! This subroutine gauleg, originates from the f90.lib fortran library. 
! It uses Gauss-Legendre polynomials to numerically estimate integrals.
!  
!      This routine calculates gauss-legendre mesh points and weights
!      input: 
!      x1   : lower limit of the integration interval
!      x2   : upper limit ---------- "" -------------
!      n    : the desired number of mesh points
!      output :
!      x     : gauss-legendre mesh points on the interval (x1,x2)   
!      w     : the corresponding weights

SUBROUTINE gauleg(x1,x2,x,w,n)
	IMPLICIT NONE

	INTEGER :: i, j, m, n
	DOUBLE PRECISION :: eps, x1, x2, x, w 
	DIMENSION :: x(n), w(n) 
	PARAMETER (eps=3.D-14)
	DOUBLE PRECISION :: p1,p2,p3,pp,xl,xm,z,z1

	m  = (n + 1) / 2
	xm = 0.5d0 * (x2 + x1)
	xl = 0.5d0 * (x2 - x1)

	DO i = 1, m
		z1 = 0.
		z  = COS(3.141592654d0 * (i - .25d0) / (n + .5d0))

			DO WHILE(ABS(z - z1) > EPS)
				p1 = 1.
				p2 = 0.

				DO j = 1, n
				   p3 = p2
				   p2 = p1
				   p1 = ((2. * j - 1.) * z * p2 - (j - 1.) * p3) / j
				ENDDO 

				pp = n * (z * p1 - p2) / (z * z - 1.)
				z1 = z
				z  = z - p1 / pp
			ENDDO

		x(i)		 = xm-xl*z
		x(n + 1 - i) = xm+xl*z
		w(i)		 = 2.* xl / ((1. - z * z) * (pp**2))
		w(n + 1 - i) = w(i)
	ENDDO
END SUBROUTINE gauleg

!***********************************************************************************************
!  This subroutine transforms the meshpoints and weights for the neutrino x, so it ranges from 0
!  to infinity, instead of from -1 to 1.
!
SUBROUTINE rel_mesh(x, w, n)
	IMPLICIT NONE

	INTEGER :: i, n
	DOUBLE PRECISION :: c, pi_over_4,cb, hbarc, nuc_mass, x, w
	DIMENSION :: x(n), w(n)

	PARAMETER(c			= 0.75)

	pi_over_4 = ACOS(-1.) / 4.
	hbarc	  = (3.151 *(10**(-26))) / (1.602*(10**(-19))) ! eV
	nuc_mass  = (1.6606*(10**(-27))) / (1.602*(10**(-19))) ! eV / m^2 (c^2?)

	! Estimate integrals
	CALL gauleg(-1.D0, 1.d0, x, w, n)

	! Transform mesh-points & weights
	DO i=1,n
		x(i) = TAN(pi_over_4 * (x(i) + 1.)) * c
		w(i) = (pi_over_4 * c) / (COS(pi_over_4 * (x(i) + 1.))**2) * w(i)
	ENDDO

	!this_array%mesh_point=hbarc*this_array%mesh_point
	!this_array%weight_point=hbarc*nuc_mass*this_array%weight_point
	!this_array%sp_energy=this_array%mesh_point**2
END SUBROUTINE rel_mesh





!***********************************************************************************************
!This module contains subroutines to manipulate	the indata to the program.
!Allocate, read, write,	convert, calculate
!


MODULE particle_density
	TYPE, PUBLIC ::	general_array
	      INTEGER :: n_data
	      DOUBLE PRECISION,	DIMENSION(:), POINTER :: total_baryon, neutron,	&
		    proton, electron, muon, sigma_minus, lambda, sigma0, sigma_plus, &
		    epsilon_nu_npe,epsilon_nu_npm,&			
			epsilon_nu_lpe,epsilon_nu_lpe2,epsilon_nu_lpe3,epsilon_nu_lpetot,&		    
			epsilon_nu_lpm,&			
			epsilon_nu_smne,epsilon_nu_smne2,epsilon_nu_smne3,epsilon_nu_smnetot,&
			epsilon_nu_smnm,epsilon_nu_smnm2,epsilon_nu_smnm3,epsilon_nu_smnmtot,&
			epsilon_nu_smle,epsilon_nu_smlm,epsilon_nu_sms0e,&
		    epsilon_nu_sms0m, epsilon_tot1
	    END	TYPE general_array

	TYPE (general_array), PUBLIC ::	density_distribution, fermilevels_array, &
					energies_array,	chem_array, epsilon_array

CONTAINS



!***********************************************************************************************
!This subroutine allocates room	for densities, fermilevels and energies
!

	SUBROUTINE allocate_dd_array(this_array)
		IMPLICIT NONE
		INTEGER	 :: i
		TYPE(general_array), INTENT(INOUT) :: this_array
		INTEGER	:: n
		n = this_array%n_data
		IF (ASSOCIATED(this_array%total_baryon)) DEALLOCATE(this_array%total_baryon)
		IF (ASSOCIATED(this_array%neutron)) DEALLOCATE(this_array%neutron)
		IF (ASSOCIATED(this_array%proton)) DEALLOCATE(this_array%proton)
		IF (ASSOCIATED(this_array%electron)) DEALLOCATE(this_array%electron)
		IF (ASSOCIATED(this_array%muon)) DEALLOCATE(this_array%muon)
		IF (ASSOCIATED(this_array%sigma_minus))	DEALLOCATE(this_array%sigma_minus)
		IF (ASSOCIATED(this_array%lambda)) DEALLOCATE(this_array%lambda)
		IF (ASSOCIATED(this_array%sigma0)) DEALLOCATE(this_array%sigma0)
		IF (ASSOCIATED(this_array%sigma_plus)) DEALLOCATE(this_array%sigma_plus)

! This allocates room for the particle densities

		ALLOCATE(this_array%total_baryon(n))
		ALLOCATE(this_array%neutron(n))
		ALLOCATE(this_array%proton(n))
		ALLOCATE(this_array%electron(n))
		ALLOCATE(this_array%muon(n))
		ALLOCATE(this_array%sigma_minus(n))
		ALLOCATE(this_array%lambda(n))
		ALLOCATE(this_array%sigma0(n))
		ALLOCATE(this_array%sigma_plus(n))

		DO i = 1,this_array%n_data
		   this_array%total_baryon(i) =	0
		   this_array%neutron(i) = 0
		   this_array%proton(i)=0
		   this_array%electron(i)=0
		   this_array%muon(i)=0
		   this_array%sigma_minus(i)=0
		   this_array%lambda(i)=0
		   this_array%sigma0(i)=0
		   this_array%sigma_plus(i)=0
		ENDDO

	END SUBROUTINE allocate_dd_array


!********************************************************************************************
!This subroutine reads the data	from the dens.dat file to a new	set of variables
!

	SUBROUTINE read_dd_data(this_array)
		IMPLICIT NONE
		INTEGER	 :: i
		TYPE(general_array), INTENT(INOUT) :: this_array
		CHARACTER (LEN = 100), POINTER ::densities
		DOUBLE PRECISION :: tot1, ne1, p1, e_min1, mu_min1, sigma_min1,	&
				    lambda1, sigma_null1, sigma_pluss1


	DO i=1,	this_array%n_data
		READ(1,*)tot1, ne1, p1,	e_min1,	mu_min1, sigma_min1, lambda1,&
			 sigma_null1, sigma_pluss1

		this_array%total_baryon(i)=tot1
		this_array%neutron(i)=ne1
		this_array%proton(i)=p1
		this_array%electron(i)=e_min1
		this_array%muon(i)=mu_min1
		this_array%sigma_minus(i)=sigma_min1
		this_array%lambda(i)=lambda1
		this_array%sigma0(i)=sigma_null1
		this_array%sigma_plus(i)=sigma_pluss1
	ENDDO

	END SUBROUTINE read_dd_data

!********************************************************************************************
!This subroutine writes	the data from the dens.dat file	to a new file densout.dat
!

	SUBROUTINE write_dd_data(this_array)
		IMPLICIT NONE
		INTEGER	 :: i
		TYPE(general_array), INTENT(IN)	:: this_array
!	    CHARACTER (LEN = 100), POINTER :: densities
		DOUBLE PRECISION :: tot1, ne1, p1, e_min1, mu_min1, sigma_min1,	&
				lambda1, sigma_null1, sigma_pluss1
		OPEN(UNIT=2, FILE='densout.dat')


		DO i= 1,this_array%n_data


		  WRITE(2,'(9(d12.6,2x))') this_array%total_baryon(i),	this_array%neutron(i), &
			      this_array%proton(i),this_array%electron(i), &
			      this_array%muon(i), this_array%sigma_minus(i),&
			      this_array%lambda(i),this_array%sigma0(i),&
			      this_array%sigma_plus(i)
		ENDDO

	 END SUBROUTINE	write_dd_data

!****************************************************************************************************
! This subroutine converts the densities to find the fermilevels, and subsequently writes the 
! fermilevels to a file "fermilevels.dat".
!

	 SUBROUTINE fermilevel_convertion
	    IMPLICIT NONE
	    INTEGER  ::	i
!	    TYPE(general_array)this_array
	    DOUBLE PRECISION ::	pi2
	   pi2 = acos(-1.0d0)**2.0d0

		OPEN(UNIT=5,FILE="fermilevels.dat")

	   DO i	= 1,fermilevels_array%n_data
	      fermilevels_array%neutron(i) = (density_distribution%neutron(i)*3*pi2)**(1./3)
	      fermilevels_array%proton(i) = (density_distribution%proton(i)*3*pi2)**(1./3)
	      fermilevels_array%electron(i) = (density_distribution%electron(i)*3*pi2)**(1./3)
	      fermilevels_array%muon(i)	= (density_distribution%muon(i)*3*pi2)**(1./3)
	      fermilevels_array%sigma_minus(i) = (density_distribution%sigma_minus(i)*3*pi2)**(1./3)
	      fermilevels_array%lambda(i) = (density_distribution%lambda(i)*3*pi2)**(1./3)
	      fermilevels_array%sigma0(i) = (density_distribution%sigma0(i)*3*pi2)**(1./3)
	      fermilevels_array%sigma_plus(i) =	(density_distribution%sigma_plus(i)*3*pi2)**(1./3)

	   WRITE (5,'(8(d12.6,2x))')fermilevels_array%neutron(i),fermilevels_array%proton(i)&
                ,fermilevels_array%electron(i),&
		fermilevels_array%muon(i),fermilevels_array%sigma_minus(i),fermilevels_array%lambda(i),&
		fermilevels_array%sigma0(i),fermilevels_array%sigma_plus(i)
	   ENDDO

	 END SUBROUTINE	fermilevel_convertion
!****************************************************************************************************
! Now we must read from	the chem.dat file to get the value for the electron/muon we need in the	
! calculations to find the emissivity.

	SUBROUTINE read_chem_data(this_array)
		IMPLICIT NONE
		INTEGER	 :: i

		TYPE(general_array), INTENT(INOUT) :: this_array
		DOUBLE PRECISION :: ne1chem, p1chem, e_min1chem,  sigma_min1chem, &
				    lambda1chem, sigma_null1chem, sigma_pluss1chem

		OPEN(UNIT=3, FILE='chem.dat')


	DO i=1,	this_array%n_data
		READ(3,*)ne1chem, p1chem, e_min1chem,  sigma_min1chem, &
			 lambda1chem, sigma_null1chem, sigma_pluss1chem
		this_array%neutron(i)=ne1chem
		this_array%proton(i)=p1chem
		this_array%electron(i)=e_min1chem
		this_array%sigma_minus(i)=sigma_min1chem
		this_array%lambda(i)=lambda1chem
		this_array%sigma0(i)=sigma_null1chem
		this_array%sigma_plus(i)=sigma_pluss1chem
	ENDDO

	END SUBROUTINE read_chem_data

!***********************************************************************************************
! This subroutine allocates room	for the	emissivities.
!

	SUBROUTINE allocate_eps_array(this_array)
		IMPLICIT NONE
		TYPE(general_array), INTENT(INOUT) :: this_array
		INTEGER	:: n, i
		n = this_array%n_data
		IF (ASSOCIATED(this_array%epsilon_nu_npe)) DEALLOCATE(this_array%epsilon_nu_npe)
		IF (ASSOCIATED(this_array%epsilon_nu_npm)) DEALLOCATE(this_array%epsilon_nu_npm)
		IF (ASSOCIATED(this_array%epsilon_nu_lpe)) DEALLOCATE(this_array%epsilon_nu_lpe)
				IF (ASSOCIATED(this_array%epsilon_nu_lpe2)) DEALLOCATE(this_array%epsilon_nu_lpe2)
				IF (ASSOCIATED(this_array%epsilon_nu_lpe3)) DEALLOCATE(this_array%epsilon_nu_lpe3)
				IF (ASSOCIATED(this_array%epsilon_nu_lpetot)) DEALLOCATE(this_array%epsilon_nu_lpetot)
		IF (ASSOCIATED(this_array%epsilon_nu_lpm)) DEALLOCATE(this_array%epsilon_nu_lpm)
		IF (ASSOCIATED(this_array%epsilon_nu_smne)) DEALLOCATE(this_array%epsilon_nu_smne)
				IF (ASSOCIATED(this_array%epsilon_nu_smne2)) DEALLOCATE(this_array%epsilon_nu_smne2)
				IF (ASSOCIATED(this_array%epsilon_nu_smne3)) DEALLOCATE(this_array%epsilon_nu_smne3)
				IF (ASSOCIATED(this_array%epsilon_nu_smnetot)) DEALLOCATE(this_array%epsilon_nu_smnetot)
		IF (ASSOCIATED(this_array%epsilon_nu_smnm)) DEALLOCATE(this_array%epsilon_nu_smnm)
				IF (ASSOCIATED(this_array%epsilon_nu_smnm2)) DEALLOCATE(this_array%epsilon_nu_smnm2)
				IF (ASSOCIATED(this_array%epsilon_nu_smnm3)) DEALLOCATE(this_array%epsilon_nu_smnm3)
				IF (ASSOCIATED(this_array%epsilon_nu_smnmtot)) DEALLOCATE(this_array%epsilon_nu_smnmtot)
		IF (ASSOCIATED(this_array%epsilon_nu_smle)) DEALLOCATE(this_array%epsilon_nu_smle)
		IF (ASSOCIATED(this_array%epsilon_nu_smlm)) DEALLOCATE(this_array%epsilon_nu_smlm)
		IF (ASSOCIATED(this_array%epsilon_nu_sms0e)) DEALLOCATE(this_array%epsilon_nu_sms0e)
		IF (ASSOCIATED(this_array%epsilon_nu_sms0m)) DEALLOCATE(this_array%epsilon_nu_sms0m)
		IF (ASSOCIATED(this_array%epsilon_tot1)) DEALLOCATE(this_array%epsilon_tot1)

! This allocates room for the particle densities

		ALLOCATE(this_array%epsilon_nu_npe(n))
		ALLOCATE(this_array%epsilon_nu_npm(n))
		ALLOCATE(this_array%epsilon_nu_lpe(n))
				ALLOCATE(this_array%epsilon_nu_lpe2(n))
				ALLOCATE(this_array%epsilon_nu_lpe3(n))
				ALLOCATE(this_array%epsilon_nu_lpetot(n))

		ALLOCATE(this_array%epsilon_nu_lpm(n))
		ALLOCATE(this_array%epsilon_nu_smne(n))
				ALLOCATE(this_array%epsilon_nu_smne2(n))
				ALLOCATE(this_array%epsilon_nu_smne3(n))
				ALLOCATE(this_array%epsilon_nu_smnetot(n))

		ALLOCATE(this_array%epsilon_nu_smnm(n))
				ALLOCATE(this_array%epsilon_nu_smnm2(n))
				ALLOCATE(this_array%epsilon_nu_smnm3(n))
				ALLOCATE(this_array%epsilon_nu_smnmtot(n))

		ALLOCATE(this_array%epsilon_nu_smle(n))
		ALLOCATE(this_array%epsilon_nu_smlm(n))
		ALLOCATE(this_array%epsilon_nu_sms0e(n))
		ALLOCATE(this_array%epsilon_nu_sms0m(n))
		ALLOCATE(this_array%epsilon_tot1(n))

		DO i = 1,this_array%n_data
		   this_array%epsilon_nu_npe(i)	= 0
		   this_array%epsilon_nu_npm(i)	= 0
		   this_array%epsilon_nu_lpe(i)=0
		   		 this_array%epsilon_nu_lpetot(i)=0
 		   this_array%epsilon_nu_lpm(i)=0
		   this_array%epsilon_nu_smne(i)=0
		   		this_array%epsilon_nu_smnetot(i)=0
		   this_array%epsilon_nu_smnm(i)=0
		   	    this_array%epsilon_nu_smnmtot(i)=0
		   this_array%epsilon_nu_smle(i)=0
		   this_array%epsilon_nu_smlm(i)=0 	
		   this_array%epsilon_nu_sms0e(i)=0
		   this_array%epsilon_nu_sms0m(i)=0
		   this_array%epsilon_tot1(i)=0
		ENDDO

	END SUBROUTINE allocate_eps_array

!****************************************************************************************************
! These are the functions designed to determine the emissivities for the different processes. 
! They are put through the gauleg-routine to get results from the integrations to put in the 
! loops below.
! The formulae are in units of (kBT)**6*1/(2**7*pi**5)(1+g)**2*G/2, with g=1
! x=x(baryon2), y=x(neutrino). The parameter beta is k_B*T in MeV.
!



	DOUBLE PRECISION FUNCTION eps_nu_lpe(x,	y,i)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) ::	x, y
	INTEGER, INTENT(IN) :: i
        DOUBLE PRECISION :: ml, mp, me, beta
        PARAMETER (ml=1115.684d0)  !MeV/c**2 = (*1.78268d-30 kg)
        PARAMETER (mp=938.27231d0) !Mev/c**2 = (1.672649d-27 kg)
        PARAMETER (me=0.51099906d0)!Mev/c**2 = (9.10953d-31 kg)
        PARAMETER (beta=0.08621)

! Correct calculations put in only for those processes that run.
! Starting with lpe, Ev1:

		eps_nu_lpe = (0.5*(ml**2-mp**2-me**2)*(((x/beta)+chem_array%lambda(i))*(y/beta))&
                     -((x/beta+chem_array%lambda(i))*y/beta)**2&
                     -1/3*(((x/beta+chem_array%lambda(i))**2-ml**2)*(y/beta)**2)&
                     *((0.756-1.477)**2)*(1-0.231**2)&
					 *y**2*(1+exp(x))**(-1.0)&
					 *((log(1+exp(fermilevels_array%proton(i)+y-x))&
					 -log(exp(fermilevels_array%proton(i)+1))&
					 -log(0.5*(1+exp(y-x))))/(exp(y-x)-1)))&
					 !done, now Ev1*
					 +(0.5*(ml**2-mp**2-me**2)*(((x/beta)+chem_array%lambda(i))*(y/beta))&
                     -((x/beta+chem_array%lambda(i))*y/beta)**2&
                     -1/3*(((x/beta+chem_array%lambda(i))**2-ml**2)*(y/beta)**2)&
                     *((0.756-1.477)**2)*(1-0.231**2)&
					 *y**2*(1+exp(x))**(-1.0)&
					 *((log(1+exp(fermilevels_array%proton(i)-y-x))&
					 -log(exp(fermilevels_array%proton(i)+1))&
					 -log(0.5*(1+exp(-y-x))))/(exp(-y-x)-1)))


	END FUNCTION eps_nu_lpe


! Continuing with Ev2:
	DOUBLE PRECISION FUNCTION eps_nu_lpe2(x,	y,i)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) ::	x, y
	INTEGER, INTENT(IN) :: i
        DOUBLE PRECISION :: ml, mp, me, beta
        PARAMETER (ml=1115.684d0)  !MeV/c**2 = (*1.78268d-30 kg)
        PARAMETER (mp=938.27231d0) !Mev/c**2 = (1.672649d-27 kg)
        PARAMETER (me=0.51099906d0)!Mev/c**2 = (9.10953d-31 kg)
        PARAMETER (beta=0.0862)

		eps_nu_lpe2 = (0.5*(ml**2-mp**2-me**2)*(((x/beta)+chem_array%proton(i))*(y/beta))&
                     -((x/beta+chem_array%proton(i))*y/beta)**2&
                     -1/3*(((x/beta+chem_array%proton(i))**2-ml**2)*(y/beta)**2)&
                     *((0.595226-0.816497)**2)*(0.053361)&
					 *y**2*(1+exp(-x))**(-1.0)&
					 *((log(1+exp(fermilevels_array%lambda(i)))&
					 +log(1+exp(x+y))&
					 -log(2.0*(exp(fermilevels_array%lambda(i))+exp(x+y))))/(exp(y+x)-1)))&
					 ! done, now Ev2*
				     +(0.5*(ml**2-mp**2-me**2)*(((x/beta)+chem_array%proton(i))*(y/beta))&
                     -((x/beta+chem_array%proton(i))*y/beta)**2&
                     -1/3*(((x/beta+chem_array%proton(i))**2-ml**2)*(y/beta)**2)&
                     *((0.595226-0.816497)**2)*(0.053361)&
					 *y**2*(1+exp(-x))**(-1.0)&
					 *((log(1+exp(fermilevels_array%lambda(i)))&
					 +log(1+exp(x-y))&
					 -log(2.0*(exp(fermilevels_array%lambda(i))+exp(x-y))))/(exp(x-y)-1)))
					 !done

	END FUNCTION eps_nu_lpe2


! Continuing with Ev3:
	DOUBLE PRECISION FUNCTION eps_nu_lpe3(x,	y,i)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) ::	x, y
	INTEGER, INTENT(IN) :: i
        DOUBLE PRECISION :: ml, mp, me, beta
        PARAMETER (ml=1115.684d0)  !MeV/c**2 = (*1.78268d-30 kg)
        PARAMETER (mp=938.27231d0) !Mev/c**2 = (1.672649d-27 kg)
        PARAMETER (me=0.51099906d0)!Mev/c**2 = (9.10953d-31 kg)
        PARAMETER (beta=0.0862)

		eps_nu_lpe3 = 0.053361*2.5935*ml*mp*(x/beta**2)*(y**3)*((1+exp(-x))**-1)&
					 *(log(1+exp(fermilevels_array%lambda(i)))&
					 +log(1+exp(x+y))&
					 -log(2.0*(exp(fermilevels_array%lambda(i))+exp(x+y))))/(exp(y+x)-1)&
					 !done, now for Ev3*
					 +0.053361*2.5935*ml*mp*(x/beta**2)*(y**3)*((1+exp(-x))**-1)&
					 *(log(1+exp(fermilevels_array%lambda(i)))&
					 +log(1+exp(x-y))&
					 -log(2.0*(exp(fermilevels_array%lambda(i))+exp(x-y))))/(exp(x-y)-1)
					 ! done
				
	END FUNCTION eps_nu_lpe3






! Continuing with smne, Ev1:
	DOUBLE PRECISION FUNCTION eps_nu_smne(x,y,i)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) ::	x, y
	INTEGER, INTENT(IN) :: i
        DOUBLE PRECISION :: mn, msm, me, beta
        PARAMETER (mn=939.56563d0) !Mev/c**2 = (1.67470d-27  kg)
        PARAMETER (msm=1197.436d0) !Mev/c**2 = (*1.78268d-30 kg)
        PARAMETER (me=0.51099906d0)!Mev/c**2 = (9.10953d-31 kg)
        PARAMETER (beta=0.0862)

		eps_nu_smne = (0.5*(msm**2-mn**2-me**2)*(((x/beta)+chem_array%sigma_minus(i))*(y/beta))&
                     -((x/beta+chem_array%sigma_minus(i))*y/beta)**2&
                     -1/3*(((x/beta+chem_array%sigma_minus(i))**2-msm**2)*(y/beta)**2)&
                     *((0.756-1.477)**2)*(1-0.231**2)&
					 *y**2*(1+exp(x))**(-1.0)&
					 *((log(1+exp(fermilevels_array%neutron(i)+y-x))&
					 -log(exp(fermilevels_array%neutron(i)+1))&
					 -log(0.5*(1+exp(y-x))))/(exp(y-x)-1)))&
					 !done, now Ev1*
					 +(0.5*(msm**2-mn**2-me**2)*(((x/beta)+chem_array%sigma_minus(i))*(y/beta))&
                     -((x/beta+chem_array%sigma_minus(i))*y/beta)**2&
                     -1/3*(((x/beta+chem_array%sigma_minus(i))**2-msm**2)*(y/beta)**2)&
                     *((0.756-1.477)**2)*(1-0.231**2)&
					 *y**2*(1+exp(x))**(-1.0)&
					 *((log(1+exp(fermilevels_array%neutron(i)-y-x))&
					 -log(exp(fermilevels_array%neutron(i)+1))&
					 -log(0.5*(1+exp(-y-x))))/(exp(-y-x)-1)))
					 !done

	END FUNCTION eps_nu_smne


	DOUBLE PRECISION FUNCTION eps_nu_smne2(x,y,i)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) ::	x, y
	INTEGER, INTENT(IN) :: i
        DOUBLE PRECISION :: mn, msm, me, beta
        PARAMETER (mn=939.56563d0) !Mev/c**2 = (1.67470d-27  kg)
        PARAMETER (msm=1197.436d0) !Mev/c**2 = (*1.78268d-30 kg)
        PARAMETER (me=0.51099906d0)!Mev/c**2 = (9.10953d-31 kg)
        PARAMETER (beta=0.0862)

		eps_nu_smne2 = (0.5*(msm**2-mn**2-me**2)*(((x/beta)+chem_array%neutron(i))*(y/beta))&
                     -((x/beta+chem_array%neutron(i))*y/beta)**2&
                     -1/3*(((x/beta+chem_array%neutron(i))**2-mn**2)*(y/beta)**2)&
                     *((0.595226-0.816497)**2)*(0.053361)&
					 *y**2*(1+exp(-x))**(-1.0)&
					 *((log(1+exp(fermilevels_array%sigma_minus(i)))&
					 +log(1+exp(x+y))&
					 -log(2.0*(exp(fermilevels_array%sigma_minus(i))+exp(x+y))))/(exp(y+x)-1)))&
					 !done
					 +(0.5*(msm**2-mn**2-me**2)*(((x/beta)+chem_array%neutron(i))*(y/beta))&
                     -((x/beta+chem_array%neutron(i))*y/beta)**2&
                     -1/3*(((x/beta+chem_array%neutron(i))**2-mn**2)*(y/beta)**2)&
                     *((0.595226-0.816497)**2)*(0.053361)&
					 *y**2*(1+exp(-x))**(-1.0)&
					 *((log(1+exp(fermilevels_array%sigma_minus(i)))&
					 +log(1+exp(x-y))&
					 -log(2.0*(exp(fermilevels_array%sigma_minus(i))+exp(x-y))))/(exp(x-y)-1)))

	END FUNCTION eps_nu_smne2


	DOUBLE PRECISION FUNCTION eps_nu_smne3(x,y,i)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) ::	x, y
	INTEGER, INTENT(IN) :: i
        DOUBLE PRECISION :: mn, msm, me, beta
        PARAMETER (mn=939.56563d0) !Mev/c**2 = (1.67470d-27  kg)
        PARAMETER (msm=1197.436d0) !Mev/c**2 = (*1.78268d-30 kg)
        PARAMETER (me=0.51099906d0)!Mev/c**2 = (9.10953d-31 kg)
        PARAMETER (beta=0.0862)

		eps_nu_smne3 = 0.053361*(-(1-0.07784))*msm*mn*x*(y**3)*((1+exp(-x))**-1)&
					 *(log(1+exp(fermilevels_array%sigma_minus(i)))&
					 +log(1+exp(x+y))&
					 -log(2.0*(exp(fermilevels_array%sigma_minus(i))+exp(x+y))))/(exp(y+x)-1)&
					 ! done, now *
					 +0.053361*(-(1-0.07784))*msm*mn*x*(y**3)*((1+exp(-x))**-1)&
					 *(log(1+exp(fermilevels_array%sigma_minus(i)))&
					 +log(1+exp(x-y))&
					 -log(2.0*(exp(fermilevels_array%sigma_minus(i))+exp(x-y))))/(exp(x-y)-1)
					 ! done
	END FUNCTION eps_nu_smne3



	DOUBLE PRECISION FUNCTION eps_nu_smnm(x,y,i)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) ::	x, y
	INTEGER, INTENT(IN) :: i
        DOUBLE PRECISION :: mn, msm, mm, beta
		PARAMETER (mn=939.56563d0) !Mev/c**2 = (1.67470d-27  kg)
        PARAMETER (msm=1197.436d0) !Mev/c**2 = (*1.78268d-30 kg)
		PARAMETER (mm=105.658389d0)!MeV/c**2 = (9.10953d-31 kg)
        PARAMETER (beta=0.0862)
		eps_nu_smnm = (0.5*(msm**2-mn**2-mm**2)*(((x/beta)+chem_array%sigma_minus(i))*(y/beta))&
                     -((x/beta+chem_array%sigma_minus(i))*y/beta)**2&
                     -1/3*(((x/beta+chem_array%sigma_minus(i))**2-msm**2)*(y/beta)**2)&
                     *((0.756-1.477)**2)*(1-0.231**2)&
					 *y**2*(1+exp(x))**(-1.0)&
					 *((log(1+exp(fermilevels_array%neutron(i)+y-x))&
					 -log(exp(fermilevels_array%neutron(i)+1))&
					 -log(0.5*(1+exp(y-x))))/(exp(y-x)-1)))&
					 !done
					 +(0.5*(msm**2-mn**2-mm**2)*(((x/beta)+chem_array%sigma_minus(i))*(y/beta))&
                     -((x/beta+chem_array%sigma_minus(i))*y/beta)**2&
                     -1/3*(((x/beta+chem_array%sigma_minus(i))**2-msm**2)*(y/beta)**2)&
                     *((0.756-1.477)**2)*(1-0.231**2)&
					 *y**2*(1+exp(x))**(-1.0)&
					 *((log(1+exp(fermilevels_array%neutron(i)-y-x))&
					 -log(exp(fermilevels_array%neutron(i)+1))&
					 -log(0.5*(1+exp(-y-x))))/(exp(-y-x)-1)))
					 !done

	END FUNCTION eps_nu_smnm


	DOUBLE PRECISION FUNCTION eps_nu_smnm2(x,y,i)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) ::	x, y
	INTEGER, INTENT(IN) :: i
        DOUBLE PRECISION :: mn, msm, mm, beta
		PARAMETER (mn=939.56563d0) !Mev/c**2 = (1.67470d-27  kg)
        PARAMETER (msm=1197.436d0) !Mev/c**2 = (*1.78268d-30 kg)
		PARAMETER (mm=105.658389d0)!MeV/c**2 = (9.10953d-31 kg)
        PARAMETER (beta=0.0862)

		eps_nu_smnm2 = (0.5*(msm**2-mn**2-mm**2)*(((x/beta)+chem_array%neutron(i))*(y/beta))&
                     -((x/beta+chem_array%neutron(i))*y/beta)**2&
                     -1/3*(((x/beta+chem_array%neutron(i))**2-mn**2)*(y/beta)**2)&
                     *((0.595226-0.816497)**2)*(0.053361)&
					 *y**2*(1+exp(-x))**(-1.0)&
					 *((log(1+exp(fermilevels_array%sigma_minus(i)))&
					 +log(1+exp(x+y))&
					 -log(2.0*(exp(fermilevels_array%sigma_minus(i))+exp(x+y))))/(exp(y+x)-1)))&
					 !done
					 +(0.5*(msm**2-mn**2-mm**2)*(((x/beta)+chem_array%neutron(i))*(y/beta))&
                     -((x/beta+chem_array%neutron(i))*y/beta)**2&
                     -1/3*(((x/beta+chem_array%neutron(i))**2-mn**2)*(y/beta)**2)&
                     *((0.595226-0.816497)**2)*(0.053361)&
					 *y**2*(1+exp(-x))**(-1.0)&
					 *((log(1+exp(fermilevels_array%sigma_minus(i)))&
					 +log(1+exp(x-y))&
					 -log(2.0*(exp(fermilevels_array%sigma_minus(i))+exp(x-y))))/(exp(x-y)-1)))
					 ! done

	END FUNCTION eps_nu_smnm2


	DOUBLE PRECISION FUNCTION eps_nu_smnm3(x,y,i)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) ::	x, y
	INTEGER, INTENT(IN) :: i
        DOUBLE PRECISION :: mn, msm, mm, beta
		PARAMETER (mn=939.56563d0) !Mev/c**2 = (1.67470d-27  kg)
        PARAMETER (msm=1197.436d0) !Mev/c**2 = (*1.78268d-30 kg)
		PARAMETER (mm=105.658389d0)!MeV/c**2 = (9.10953d-31 kg)
        PARAMETER (beta=0.0862)

		eps_nu_smnm3 = 0.053361*(-(1-0.07784))*msm*mn*(x/beta**2)*(y**3)*((1+exp(-x))**-1)&
					 *(log(1+exp(fermilevels_array%sigma_minus(i)))&
					 +log(1+exp(x+y))&
					 -log(2.0*(exp(fermilevels_array%sigma_minus(i))+exp(x+y))))/(exp(y+x)-1)&
					 ! done, now *
					 +0.053361*(-(1-0.07784))*msm*mn*(x/beta**2)*(y**3)*((1+exp(-x))**-1)&
					 *(log(1+exp(fermilevels_array%sigma_minus(i)))&
					 +log(1+exp(x-y))&
					 -log(2.0*(exp(fermilevels_array%sigma_minus(i))+exp(x-y))))/(exp(x-y)-1)
						! done

	END FUNCTION eps_nu_smnm3
	



!****************************************************************************************************
! This subroutine tests the triangle inequality to see which processes can go, and
! calculate the	neutrino emissivity epsilon_nu
! The formulae are in units of (k_B*T)**6*1/(2**7*pi**5)(1+g)**2*G/2, with g=1
! 
!
	 SUBROUTINE test_fermilevel
	 IMPLICIT NONE
	 DOUBLE	PRECISION::root23, step
	 INTEGER :: i, nump, ix, iy, iz, iw
	 DOUBLE	PRECISION, ALLOCATABLE,	DIMENSION(:) ::	mp_neutron, wgt_neutron
	 DOUBLE	PRECISION, ALLOCATABLE,	DIMENSION(:) ::	mp_proton, wgt_proton
	 DOUBLE	PRECISION, ALLOCATABLE,	DIMENSION(:) ::	mp_electron, wgt_electron
	 DOUBLE	PRECISION, ALLOCATABLE,	DIMENSION(:) ::	mp_muon, wgt_muon
	 DOUBLE	PRECISION, ALLOCATABLE,	DIMENSION(:) ::	mp_sigma_minus,	wgt_sigma_minus
	 DOUBLE	PRECISION, ALLOCATABLE,	DIMENSION(:) ::	mp_lambda, wgt_lambda
	 DOUBLE	PRECISION, ALLOCATABLE,	DIMENSION(:) ::	mp_sigma0, wgt_sigma0
	 DOUBLE	PRECISION, ALLOCATABLE,	DIMENSION(:) ::	mp_sigma_pluss,	wgt_sigma_pluss

	 DOUBLE	PRECISION, ALLOCATABLE,	DIMENSION(:) ::	mp_neutron_y, wgt_neutron_y
	 DOUBLE	PRECISION, ALLOCATABLE,	DIMENSION(:) ::	mp_proton_y, wgt_proton_y
	 DOUBLE	PRECISION, ALLOCATABLE,	DIMENSION(:) ::	mp_electron_y, wgt_electron_y
	 DOUBLE	PRECISION, ALLOCATABLE,	DIMENSION(:) ::	mp_muon_y, wgt_muon_y
	 DOUBLE	PRECISION, ALLOCATABLE,	DIMENSION(:) ::	mp_sigma_minus_y,	wgt_sigma_minus_y
	 DOUBLE	PRECISION, ALLOCATABLE,	DIMENSION(:) ::	mp_lambda_y, wgt_lambda_y
	 DOUBLE	PRECISION, ALLOCATABLE,	DIMENSION(:) ::	mp_sigma0_y, wgt_sigma0_y
	 DOUBLE	PRECISION, ALLOCATABLE,	DIMENSION(:) ::	mp_sigma_pluss_y,	wgt_sigma_pluss_y

	 root23	= (2.0d0/3.0d0)**(1.0d0/2.0d0)

!		Define function	interfaces for each process.. (10)

!		Allocate arrays...
		WRITE(*,*) ' Read in number of mesh points [nump]'
		READ(*,*) nump

		! Allocate mesh-points and weights for baryon2 (x)
		ALLOCATE ( mp_neutron(nump), wgt_neutron(nump) )
		ALLOCATE ( mp_proton(nump), wgt_proton(nump) )
		ALLOCATE ( mp_electron(nump), wgt_electron(nump) )
		ALLOCATE ( mp_muon(nump), wgt_muon(nump) )
		ALLOCATE ( mp_sigma_minus(nump), wgt_sigma_minus(nump) )
		ALLOCATE ( mp_lambda(nump), wgt_lambda(nump) )
		ALLOCATE ( mp_sigma0(nump), wgt_sigma0(nump) )
		ALLOCATE ( mp_sigma_pluss(nump), wgt_sigma_pluss(nump) )

		! Allocate mesh-points and weigths for neutrinos (y)
		ALLOCATE ( mp_neutron_y(nump), wgt_neutron_y(nump) )
		ALLOCATE ( mp_proton_y(nump), wgt_proton_y(nump) )
		ALLOCATE ( mp_electron_y(nump), wgt_electron_y(nump) )
		ALLOCATE ( mp_muon_y(nump), wgt_muon_y(nump) )
		ALLOCATE ( mp_sigma_minus_y(nump), wgt_sigma_minus_y(nump) )
		ALLOCATE ( mp_lambda_y(nump), wgt_lambda_y(nump) )
		ALLOCATE ( mp_sigma0_y(nump), wgt_sigma0_y(nump) )
		ALLOCATE ( mp_sigma_pluss_y(nump), wgt_sigma_pluss_y(nump) )

	 DO i =	1, epsilon_array%n_data
!*************************************************************************************
!		We now call the gauleg routine and set up mesh points and weights for all 
!       particles.

		CALL gauleg(0.0d0,	fermilevels_array%neutron(i), mp_neutron, wgt_neutron, nump)
		CALL gauleg(0.0d0,	fermilevels_array%proton(i), mp_proton,	wgt_proton, nump)
		CALL gauleg(0.0d0,	fermilevels_array%electron(i), mp_electron, wgt_electron, nump)
		CALL gauleg(0.0d0,	fermilevels_array%muon(i), mp_muon, wgt_muon, nump)
		CALL gauleg(0.0d0,	fermilevels_array%sigma_minus(i), mp_sigma_minus, wgt_sigma_minus, nump)
		CALL gauleg(0.0d0,	fermilevels_array%lambda(i), mp_lambda,	wgt_lambda, nump)
		CALL gauleg(0.0d0,	fermilevels_array%sigma0(i), mp_sigma0,	wgt_sigma0, nump)
		CALL gauleg(0.0d0,	fermilevels_array%sigma_plus(i), mp_sigma_pluss, wgt_sigma_pluss, nump)

		CALL gauleg(0.0d0,  200.0d0, mp_neutron_y,		wgt_neutron_y,		nump)
		CALL gauleg(0.0d0,  200.0d0, mp_proton_y,		wgt_proton_y,		nump)
		CALL gauleg(0.0d0,  200.0d0, mp_electron_y,	    wgt_electron_y,		nump)
		CALL gauleg(0.0d0,  200.0d0, mp_muon_y,		    wgt_muon_y,			nump)
		CALL gauleg(0.0d0,  200.0d0, mp_sigma_minus_y,	wgt_sigma_minus_y,	nump)
		CALL gauleg(0.0d0,  200.0d0, mp_lambda_y,		wgt_lambda_y,		nump)
		CALL gauleg(0.0d0,  200.0d0, mp_sigma0_y,		wgt_sigma0_y,		nump)
		CALL gauleg(0.0d0,  200.0d0, mp_sigma_pluss_y,	wgt_sigma_pluss_y,	nump)

	! First	we test	the process n >	p e v

	
	! Then the process  n >	p mu v

	 

! Then the process  lambda > p e v
	   IF  (fermilevels_array%electron(i)>0&
	   .AND.  fermilevels_array%proton(i)>0&
	   .AND.  fermilevels_array%lambda(i) >	0&
	   .AND.  fermilevels_array%electron(i)	+ fermilevels_array%proton(i) -	&
	       fermilevels_array%lambda(i) > 0)	THEN

				DO ix=1,nump
					DO iy=1,nump
						epsilon_array%epsilon_nu_lpe(i)	= epsilon_array%epsilon_nu_lpe(i)+&
						eps_nu_lpe(mp_lambda(ix), mp_electron_y(iy), i) * wgt_lambda(ix) *&
                                                wgt_electron_y(iy)
					ENDDO
				ENDDO


	   ENDIF

	   	IF  (fermilevels_array%electron(i)>0&
	   .AND.  fermilevels_array%proton(i)>0&
	   .AND.  fermilevels_array%lambda(i) >	0&
	   .AND.  fermilevels_array%electron(i)	+ fermilevels_array%proton(i) -	&
	       fermilevels_array%lambda(i) > 0)	THEN

				DO ix=1,nump
					DO iy=1,nump
						epsilon_array%epsilon_nu_lpe2(i)	= epsilon_array%epsilon_nu_lpe2(i)+&
						eps_nu_lpe2(mp_lambda(ix), mp_electron_y(iy), i) * wgt_lambda(ix) *&
                                                wgt_electron_y(iy)
					ENDDO
				ENDDO


	   ENDIF

	   IF  (fermilevels_array%electron(i)>0&
	   .AND.  fermilevels_array%proton(i)>0&
	   .AND.  fermilevels_array%lambda(i) >	0&
	   .AND.  fermilevels_array%electron(i)	+ fermilevels_array%proton(i) -	&
	       fermilevels_array%lambda(i) > 0)	THEN

				DO ix=1,nump
					DO iy=1,nump
						epsilon_array%epsilon_nu_lpe3(i)	= epsilon_array%epsilon_nu_lpe3(i)+&
						eps_nu_lpe3(mp_lambda(ix), mp_electron_y(iy), i) * wgt_lambda(ix) *&
                                                wgt_electron_y(iy)
					ENDDO
				ENDDO


	   ENDIF


	  


! Then the process  lambda > p mu v



! Then the process  sigma_minus	> n e v
	   IF  (fermilevels_array%electron(i)>0&
	   .AND.  fermilevels_array%neutron(i)>0&
	   .AND.  fermilevels_array%sigma_minus(i) > 0&
	   .AND.  fermilevels_array%electron(i)	+ fermilevels_array%neutron(i) - &
	      fermilevels_array%sigma_minus(i) > 0) THEN

				DO ix=1,nump
					DO iy=1,nump
						epsilon_array%epsilon_nu_smne(i)= epsilon_array%epsilon_nu_smne(i)+&
						eps_nu_smne(mp_sigma_minus(ix), mp_electron_y(iy), i) * wgt_sigma_minus(ix)*&
                                                wgt_electron_y(iy)
					ENDDO
				ENDDO


	   ENDIF

	   IF  (fermilevels_array%electron(i)>0&
	   .AND.  fermilevels_array%neutron(i)>0&
	   .AND.  fermilevels_array%sigma_minus(i) > 0&
	   .AND.  fermilevels_array%electron(i)	+ fermilevels_array%neutron(i) - &
	      fermilevels_array%sigma_minus(i) > 0) THEN

				DO ix=1,nump
					DO iy=1,nump
						epsilon_array%epsilon_nu_smne2(i)= epsilon_array%epsilon_nu_smne2(i)+&
						eps_nu_smne2(mp_sigma_minus(ix), mp_electron_y(iy), i) * wgt_sigma_minus(ix)*&
                                                wgt_electron_y(iy)
					ENDDO
				ENDDO


	   ENDIF

	   IF  (fermilevels_array%electron(i)>0&
	   .AND.  fermilevels_array%neutron(i)>0&
	   .AND.  fermilevels_array%sigma_minus(i) > 0&
	   .AND.  fermilevels_array%electron(i)	+ fermilevels_array%neutron(i) - &
	      fermilevels_array%sigma_minus(i) > 0) THEN

				DO ix=1,nump
					DO iy=1,nump
						epsilon_array%epsilon_nu_smne3(i)= epsilon_array%epsilon_nu_smne3(i)+&
						eps_nu_smne3(mp_sigma_minus(ix), mp_electron_y(iy), i) * wgt_sigma_minus(ix)*&
                                                wgt_electron_y(iy)
					ENDDO
				ENDDO


	   ENDIF

! Then the process  sigma_minus	> n mu v
	   IF  (fermilevels_array%muon(i)>0&
	   .AND.  fermilevels_array%neutron(i)>0&
	   .AND.  fermilevels_array%sigma_minus(i) > 0&
	   .AND.  fermilevels_array%muon(i) + fermilevels_array%neutron(i) - &
	      fermilevels_array%sigma_minus(i) > 0) THEN

				DO ix=1,nump
					DO iy=1,nump
						epsilon_array%epsilon_nu_smnm(i)= epsilon_array%epsilon_nu_smnm(i)+&
						eps_nu_smnm(mp_sigma_minus(ix), mp_electron_y(iy), i) * wgt_sigma_minus(ix) *&
                                                wgt_electron_y(iy)
					ENDDO
				ENDDO

	   ENDIF

	   	IF  (fermilevels_array%muon(i)>0&
	   .AND.  fermilevels_array%neutron(i)>0&
	   .AND.  fermilevels_array%sigma_minus(i) > 0&
	   .AND.  fermilevels_array%muon(i) + fermilevels_array%neutron(i) - &
	      fermilevels_array%sigma_minus(i) > 0) THEN

				DO ix=1,nump
					DO iy=1,nump
						epsilon_array%epsilon_nu_smnm2(i)= epsilon_array%epsilon_nu_smnm2(i)+&
						eps_nu_smnm2(mp_sigma_minus(ix), mp_electron_y(iy), i) * wgt_sigma_minus(ix) *&
                                                wgt_electron_y(iy)
					ENDDO
				ENDDO

	   ENDIF

	
	   IF  (fermilevels_array%muon(i)>0&
	   .AND.  fermilevels_array%neutron(i)>0&
	   .AND.  fermilevels_array%sigma_minus(i) > 0&
	   .AND.  fermilevels_array%muon(i) + fermilevels_array%neutron(i) - &
	      fermilevels_array%sigma_minus(i) > 0) THEN

				DO ix=1,nump
					DO iy=1,nump
						epsilon_array%epsilon_nu_smnm3(i)= epsilon_array%epsilon_nu_smnm3(i)+&
						eps_nu_smnm3(mp_sigma_minus(ix), mp_electron_y(iy), i) * wgt_sigma_minus(ix) *&
                                                wgt_electron_y(iy)
					ENDDO
				ENDDO

	   ENDIF



! Then the process  sigma_minus	> lambda e v


! Then the process  sigma_minus	> lambda mu v
	   

! Then the process  sigma_minus	> sigma0 e v

! Then the process  sigma_minus	> sigma0 mu v
	   
! Now to calculate the total emissivity	from this more accurate estimate.

	  
	  
	   epsilon_array%epsilon_nu_lpetot(i) = &
					   epsilon_array%epsilon_nu_lpe(i) +&
					   epsilon_array%epsilon_nu_lpe2(i) + &
					   epsilon_array%epsilon_nu_lpe3(i) 


	   epsilon_array%epsilon_nu_smnetot(i) = &
					   epsilon_array%epsilon_nu_smne(i) + &
					   epsilon_array%epsilon_nu_smne2(i) + &
					   epsilon_array%epsilon_nu_smne3(i) 

	   epsilon_array%epsilon_nu_smnmtot(i) = &
					   epsilon_array%epsilon_nu_smnm(i) +&
					   epsilon_array%epsilon_nu_smnm2(i) + &
					   epsilon_array%epsilon_nu_smnm3(i) 

	   epsilon_array%epsilon_tot1(i) = &
					   epsilon_array%epsilon_nu_lpe(i) +&
					   epsilon_array%epsilon_nu_lpe2(i) + &
					   epsilon_array%epsilon_nu_lpe3(i) + &
					   epsilon_array%epsilon_nu_smne(i) +&
					   epsilon_array%epsilon_nu_smne2(i) + &
					   epsilon_array%epsilon_nu_smne3(i) + &
					   epsilon_array%epsilon_nu_smnm(i) +&
					   epsilon_array%epsilon_nu_smnm2(i) + &
					   epsilon_array%epsilon_nu_smnm3(i) 
					   

	   WRITE(4,'(4(d12.6,2x))')&
					epsilon_array%epsilon_nu_lpetot(i), &
					epsilon_array%epsilon_nu_smnetot(i),  &
				    epsilon_array%epsilon_nu_smnmtot(i), &
				    epsilon_array%epsilon_tot1(i)

	   step= .2+0.01*i
	   WRITE(*,'(5(d12.6,2x))')step,&
				    epsilon_array%epsilon_nu_lpetot(i), &
					epsilon_array%epsilon_nu_smnetot(i),  &
				    epsilon_array%epsilon_nu_smnmtot(i), &
				    epsilon_array%epsilon_tot1(i)


	   WRITE(7,'(14(d12.6,2x))')step, &
					epsilon_array%epsilon_nu_lpe(i), &					
					epsilon_array%epsilon_nu_lpe2(i), &					
					epsilon_array%epsilon_nu_lpe3(i), &					
					epsilon_array%epsilon_nu_lpetot(i), &
					epsilon_array%epsilon_nu_smne(i),  &
					epsilon_array%epsilon_nu_smne2(i),  &
					epsilon_array%epsilon_nu_smne3(i),  &
					epsilon_array%epsilon_nu_smnetot(i),  &
				    epsilon_array%epsilon_nu_smnm(i), &
					epsilon_array%epsilon_nu_smnm2(i), &
					epsilon_array%epsilon_nu_smnm3(i), &
					epsilon_array%epsilon_nu_smnmtot(i), &
				    epsilon_array%epsilon_tot1(i)

	 ENDDO









	 END SUBROUTINE	test_fermilevel
!***********************************************************************************************************
!																										   *
! We now reach the main program...                                                                         *
!                                                                                                          *
!***********************************************************************************************************



END MODULE particle_density







PROGRAM	prog

  USE particle_density

  IMPLICIT NONE
  INTEGER :: n
  OPEN(UNIT=1,FILE="dens.dat")
  OPEN(UNIT=4,FILE="epsilon1.dat")
  OPEN(UNIT=7,FILE="plot11.dat")
  READ(1,*) n
  density_distribution%n_data =	n
  fermilevels_array%n_data = n
  energies_array%n_data	= n
  epsilon_array%n_data = n
  chem_array%n_data = n

  CALL allocate_dd_array(density_distribution)

  CALL allocate_dd_array(fermilevels_array)

  CALL allocate_dd_array(energies_array)

  CALL allocate_dd_array(chem_array)

  CALL allocate_eps_array(epsilon_array)

  CALL read_dd_data(density_distribution)

  CALL read_chem_data(chem_array)

  CALL write_dd_data(density_distribution)

  CALL fermilevel_convertion

  CALL test_fermilevel

END PROGRAM prog








!***************************************************************************************************
!
!TRANSITIONS
!
!neutron     > proton  + lepton	+ neutrino
!lambda	     > proton  + lepton	+ neutrino
!sigma_minus > neutron + lepton	+ neutrino
!sigma_minus > lambda  + lepton	+ neutrino
!sigma_minus > sigma0  + lepton	+ neutrino
!
!***************************************************************************************************
