% 
% This script plots the probability density function of the 2p_z orbital
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
N = 1/4/sqrt(2*pi) * (Z/a0) ^ 1.5;

% Compute contours of the squared wave function values.
% psi = N*sigma*exp(-sigma/2)*cos(theta) 

theta = 0:0.05:pi;
phi = 0:0.05:2*pi+0.1;

[theta,phi]=meshgrid(theta,phi);

siz = size(theta);
r = zeros(1,siz(2));
C = 0.0001;
for i=1:siz(2)
  c = sprintf('%g*r*r*exp(-r)*cos(%g)^2-%g',N^2,theta(1,i),C);
  r(i) = fzero(inline(c,'r'),10);
end

% eliminate NaNs and negative r values.
index = find(r<0);
r(index) = 0;
index = find(r==NaN);
r(index) = 0;

% convert to cartesian coordinates
x = sin(theta).*cos(phi);
x = x*diag(r);
y = sin(theta).*sin(phi);
y = y*diag(r);
z = cos(theta);
z = z*diag(r);

% Generate contour plot
figure(1);
hold off;
contour3(x,y,z,40);
c = sprintf('|\\Psi|^2 = %g Level Surface',C);
title(c);


% Generate surface plot
figure(2);
hold off;
surf(x,y,z);
title(c);


% Generate the radial probability distribution

% Values for r-coordinate.  There is no phi or theta dependence
% for the 1s orbitals.
r_r = [0:15/1000:15];
sigma_r = Z/a0 * r_r;
theta_r = [-pi:2*pi/1000:pi];

% Compute radial probability density
psi_r = (4*pi/3) * r_r.^2 .* (N * sigma_r .* exp(-sigma_r/2)).^2;

figure(3);
hold off;
plot(r_r,psi_r);
grid on;
title('|\Psi(r)|^2');
xlabel('r/a_0');
ylabel('Probability density');

