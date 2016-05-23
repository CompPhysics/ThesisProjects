% Globals

% hbar times c in units of MeV
global hbarc
global meq1

% Unit dependent scaling factors 
% energy (original units MeV)
global sc_energy_factor
% momentum (original units 1/fm)
global sc_mom_factor

% Proton and neutron masses and reduced mass, in units of MeV
global m_p
global m_n
global redmass_np
global mass

% Type of potential
global pot

% Symmetric or antisymmetric matter
global symmetric

% Parameters for the different potentials
global box_V0
global box_a
global V2potconst
global V2alpha
global Yamaguchilambda;
global Yamaguchibeta;

% Fermi momentum
global k_fermi

% Vector with momenta
global qvec

% Max and min angular momentum J and isospin projection T_z (All zero imply
% test spin matrice)
global jmin
global jmax
global itzmax
global itzmin

% Number of grid points in center-of-mass, energy (omega) and relative
% momentum grids
global n_cm_points
global n_omega_points
global n_rel_points
global n_theta_points
global n_phi_points

% Standard Gaussian grid points and weights mapped to the interval
% [0,infty] with n_gauss points
global n_gauss
global gausspoints
global gaussweights
global gr1
global grw1
global imgr1
global imgrw1


% Closeness limit for the equality of energies in the self-consistency check
global energytolerance

% Maximum number of iterations in self-consistency loop
global maximum_iterations;

% Degrees of polynomial fits
global n_ImGammapolyfit

% The coefficients of the energy polynomial
global energy_coeff

% The table containing the Clebsh-Gordan coefficients
global CGtable