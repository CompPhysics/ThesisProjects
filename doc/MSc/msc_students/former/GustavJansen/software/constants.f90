!    This module contains all constants and declarations 
!    of variables read in by the function read_data. These
!    variables are used by many functions.

MODULE constants
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
  ! min and max isospin projection
  INTEGER, PUBLIC :: itzmin, itzmax, strange_min, strange_max
  ! min and max total two-body angular momentum in lab frame
  INTEGER, PUBLIC :: j_lab_min, j_lab_max
  ! min and max total two-body angular momentum in Rel-CoM frame
  INTEGER, PUBLIC :: jmin, jmax
  ! max value of 2n+l = nlmax in rel and CoM frame. lmax is max relative l
  INTEGER, PUBLIC :: nlmax, lmax, nmax, nlmax_model, lab_lmax, lab_nmax
  ! number of integration points in relative and CoM momenta
  INTEGER , PUBLIC :: n_rel, n_cm_mom, n_total
  ! number of starting energies used to compute the g-matrix
  INTEGER, PUBLIC :: n_startenergy_g
  ! number of starting energies used to compute the Q-box, typically 11
  ! can then compute up to the tenth derivative of the Q-box
  INTEGER, PUBLIC :: n_startenergy_veff
  ! number of excitations (using harmonic oscillator picture) for diagrams
  INTEGER, PUBLIC :: number_homega_exct
  ! A of closed shell core
  INTEGER, PUBLIC :: mass_nucleus
  ! Particle charges in multiples of e
  INTEGER, PUBLIC :: sp_charge(6), coulomb_channel(6)
  ! define whether onebody, twobody or threebody interaction
  CHARACTER (LEN=100), PUBLIC :: n_body_int
  !    type of  renormalization and coulomb
  CHARACTER (LEN=100), PUBLIC :: type_of_renormv, coulomb_included
  !    determine the type of Pauli operator
  CHARACTER (LEN=100), PUBLIC :: pauli_operator
  ! the nn potential, cd-bonn, idaho-a, idaho-b etc, with csb and cib options
  CHARACTER (LEN=100), PUBLIC  :: type_of_nn_pot, csb, cib, type_of_yn_pot, &
    type_of_yy_pot
  ! test for square Pauli operator
  LOGICAL, PUBLIC :: square_calculation
  ! arrays of starting energies
  REAL(KIND(1.0D0)), PUBLIC, ALLOCATABLE :: e_start_g(:), wcn(:)
  ! which starting energy is used to compute the effective interaction
  REAL(KIND(1.0D0)), PUBLIC :: starting_energy
  CHARACTER (LEN=100), PUBLIC  :: keep_originalenergies
  ! oscillator energy and oscillator length
  REAL(KIND(1.0D0)), PUBLIC :: hbar_omega, oscl, cutoff
  ! Individual oscillator lengths
  REAL(KIND(1.0D0)), PUBLIC :: oscli(6)
  REAL(KIND(1.0D0)), PUBLIC :: sp_mass(6)
  REAL(KIND(1.0D0)) , PARAMETER, PUBLIC :: p_mass =   938.918725  !   2006 value of average (m_p+m_n)/2
  REAL(KIND(1.0D0)), PARAMETER, PUBLIC :: theta_rot = 0.0! 0.125
  REAL(KIND(1.0D0)), PARAMETER, PUBLIC :: hbarc = 197.326968    !  2006 value
  REAL(KIND(1.0D0)), PARAMETER, PUBLIC :: hb2ip = hbarc*hbarc/p_mass
  REAL(KIND(1.0D0)), PUBLIC, PARAMETER :: pi = 3.141592741012573
  REAL(KIND(1.0D0)), PUBLIC, PARAMETER :: pi_2 = 1.570796370506287
  REAL(KIND(1.0D0)), PUBLIC, PARAMETER :: pi_4 = 0.7853981852531433
  ! Masses of individual particles (Journal of Physics G - Review of particle
  ! physics: Vol 33 July 2006 )
  REAL(KIND(1.0D0)) , PARAMETER, PUBLIC :: proton_mass = 938.272029  ! (p. 955)
  REAL(KIND(1.0D0)) , PARAMETER, PUBLIC :: neutron_mass = 939.565360 ! (p. 963)
  REAL(KIND(1.0D0)) , PARAMETER, PUBLIC :: lambda_mass = 1115.683    ! (p. 1023)
  REAL(KIND(1.0D0)) , PARAMETER, PUBLIC :: sigmap_mass = 1189.37     ! (p. 1039)
  REAL(KIND(1.0D0)) , PARAMETER, PUBLIC :: sigma0_mass = 1192.642    ! (p. 1041)
  REAL(KIND(1.0D0)) , PARAMETER, PUBLIC :: sigmam_mass = 1197.449    ! (p. 1042)

  CHARACTER (LEN=100), PUBLIC :: include_yn, include_yy
  ! define the rank of the one-body operator, e.g., M1 --> lambda = 1, GT
  ! lambda = 1, E2 lambda=2 and so forth.
  INTEGER, PUBLIC :: lambda
  ! Number of BHF iterations needed for self-consistency 
  INTEGER, PUBLIC :: hf_iterations

  CHARACTER (LEN=100), PUBLIC :: type_of_interaction, order_of_interaction

  ! Total mass of the nucleus in question, including hyperons
  REAL(KIND(1.0D0)) , PUBLIC :: total_mass
  ! NUmber of hyperons
  INTEGER, PUBLIC :: n_lambda, n_sigmap, n_sigma0, n_sigmam

END MODULE constants

