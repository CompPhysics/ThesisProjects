% 
% This script plots the probability density function of the 1s orbital
% of a hydrogen-like atom with charge number Z and Bohr radius a0.
%
% Kevin Chu
% Autumn 2001
%

% The charge number and Bohr radius.
Z = 1;
a0 = 1;
Z_a0 = Z/a0;

% The normalization constant for the wave function.
N = 1/sqrt(pi) * (Z/a0) ^ 1.5;

% Compute contours of squared wave function values.
% psi = N * exp(-sigma);

theta = 0:0.05:pi;
phi = 0:0.05:2*pi;
[theta,phi]=meshgrid(theta,phi);

C = 0.0001;
c = sprintf('%g *exp(-2*r)- %g',N^2,C);
r = fzero(inline(c,'r'),10);
 
% convert to cartesian coordinates
x = r*sin(theta) .* cos(phi);
y = r*sin(theta) .* sin(phi);
z = r*cos(theta);

% Generate contour plot
figure(1);
hold off;
contour3(x,y,z,25);
c = sprintf('|\\Psi^2| = %g Level Surface',C);
title(c);

% Generate surface plot
figure(2);
hold off;
surf(x,y,z);
title(c);

% Generate the radial probability distribution

% Values for r-coordinate.  There is no phi or theta dependence
% for the 1s orbitals.
r = [0:.01:6];
sigma = Z/a0 * r;

% Compute radial probability density
psi_r = 4 * pi * r.^2 .* (N * exp(-sigma)).^2;

% Plot radial distibution function
figure(3);
hold off;
plot(r,psi_r);
grid on;
title('|\Psi(r)|^2');
xlabel('r/a_0');
ylabel('Probability density');

