% 
% This script plots the radial probability density function of 
% the 2s orbital of a hydrogen-like atom with charge number Z = 1
% and Bohr radius a0 = 1.
%
% Kevin Chu
% Autumn 2001
%

% The charge number and Bohr radius.
Z = 1;
a0 = 1;

% The normalization constant for the wave function.
N = 1/4/sqrt(2*pi) * (Z/a0) ^ 1.5;

% Values for r-coordinate.  There is no phi or theta dependence
% for 2s orbitals.
r = Z/a0 * [0 : .01 : 17];
phi = [0 : 0.01 : 2*pi];
theta = [0 : 0.01 : pi];


% The wave function values.
psi = N * (2 - r) .* exp(-r/2);

% Radial wave function
psi_r = 4 * pi * r.^2 .* psi.^2;

figure(1);
plot(r,psi_r);
grid on;
title('|\Psi(r)|^2');
xlabel('r/a_0');
ylabel('Probability density');
