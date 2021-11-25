% 
% This script plots the probability density function of the 3d_z^2 orbital
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
N = 1/81/sqrt(6*pi) * (Z/a0)^1.5;

% Compute contours of the squared wave function values.
% psi = N*sigma^2*exp(-sigma/3)*(3*cos(theta)^2-1)

theta = 0:0.05:pi;
phi = 0:0.05:2*pi;

[theta,phi]=meshgrid(theta,phi);

siz_theta = size(theta);
siz_phi = size(phi);
r = zeros(1,siz_theta(2));
C = 0.0001;
for i=1:siz_theta(2)
  c = sprintf('%g*r^4*exp(-2*r/3)*(3*(%g)^2-1)^2-%g',N^2,cos(theta(1,i)),C);
  r(i) = fzero(inline(c,'r'),20);
end

% eliminate NaNs and negative r values.
index = find(r<0);
r(index)=0;
index = find(r==NaN);
r(index)=0;

% convert to cartesian coordinates
x = sin(theta).*cos(phi);
x = x*diag(r);
y = sin(theta).*sin(phi);
y = y*diag(r);
z = cos(theta);
z = z*diag(r);

% Generate a contour plot
figure(1); 
hold off;
contour3(x,y,z,40);
c = sprintf('|\\Psi^2| = %g Level Surface',C);
title(c);

% Generate a surface plot
figure(2);
hold off;
surf(x,y,z);
title(c);


