% Parameters for self-consistent solving of in-medium self-energy Sigma and
% vertex function Gamma using Green's functions formalism

% Unit dependent scaling factors 
% energy (original units MeV). 
%sc_energy_factor = meq1./(hbarc.*hbarc); % For hbarc.^2/2.*mass = 1
sc_energy_factor = 1; % For hbarc = 1
% momentum (original units 1/fm)
%sc_mom_factor = meq1./hbarc; % For hbarc.^2/mass = 1
sc_mom_factor = hbarc; % For hbarc = 1
%sc_mom_factor = 1;

% Symmetric matter
symmetric = 1;

% Use angle-averaged and/or exact Pauli operator
average = 1;
exact = 0;

% Type of potential
%pot = 'cdbonn';
%pot = 'Yamaguchi';
pot = 'V20';
%pot = 'V2';
%pot = 'box';

% Parameters for the different potentials
box_V0 = 20.0;
box_a = 2.2./hbarc;

% Parameters for the potential used by Ramos. 
% In units of MeV (presumably)
V2potconst = [-10.463 105.468 -3187.8 9924.3];
% In units of 1/fm 
V2alpha = hbarc.*(0.7).*[1 2 4 6];

% Parameters for the Yamaguchi potential
Yamaguchilambda = 450;
Yamaguchibeta = 6.0.*sc_energy_factor;

% Current mass to be used
%mass = m_p;
%mass = 1.0;
%mass = (hbarc.^2)./2;
mass = redmass_np;

% Fermi momentum
k_fermi = 1.6.*sc_mom_factor;
%k_fermi = 1.8.*sc_mom_factor;
%k_fermi = 2.*sc_mom_factor;
%k_fermi = 2;

% Momentums (in MeV) used to calculate the energy spectrum
qvec = sc_mom_factor.*(0.01:0.1:6);

% Max and min angular momentum J and isospin projection T_z (All 42 imply
% test spin matrice)
%jmin = 0;
%jmax = 2;
jmin = 42;
jmax = 42;
itzmax = 42;
itzmin = 42;

% Number of grid points in center-of-mass, energy (omega) and relative
% momentum grids
n_cm_points = 32;
n_omega_points = 50;
% NBNB!! For values of n_rel_points>100, the code in V2pot.f90 must be
% changed
n_rel_points = 128;
n_theta_points = 3;
n_phi_points = 3;
n_gauss = 32;
n_imG_int = 32;

% Gaussian mesh for interval [0,infty] with n_gauss points
[gr1,grw1] = gauleg(-1.0,1.0,n_gauss);
gausspoints = hbarc.*tan((pi./4.0).*(1+gr1'));
gaussweights = hbarc.*((pi./4.0).*grw1')./ ... 
   (cos((pi./4.0).*(1+gr1')).*(cos((pi./4.0).*(1+gr1'))));

[imgr1,imgrw1] = gauleg(-1.0,1.0,n_imG_int);

% Closeness limit for the equality of energies in the self-consistency check
energytolerance = 0.1;

% Maximum number of iterations in self-consistency loop
maximum_iterations = 15;

% Degrees of polynomials 
n_ImGammapolyfit = 1;



