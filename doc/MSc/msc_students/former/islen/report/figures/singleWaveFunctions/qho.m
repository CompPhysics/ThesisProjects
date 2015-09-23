% 
% This script plots the probability density function for a quantum 
% harmonic oscillator for the first few energy levels.
%
% Kevin Chu
% Autumn 2001
%

% a physical parameter
alpha = 1;

% energy levels
num_levels = 10;
n = 0:num_levels-1;

% The normalization constants for the wave functions.
N = 1./sqrt(2.^n .* cumprod(n+1)./(n+1)) * (alpha/pi)^(.25);

% position mesh
x = [-5:.01:5];

% Compute the |psi^2| values.
psi_sq = zeros(4, length(x));

psi_sq(1,:) = N(1)^2 * exp(-x.^2); 
psi_sq(2,:) = N(2)^2 * 4 * x.^2 .* exp(-x.^2); 
psi_sq(3,:) = N(3)^2 * (4 * x.^2 - 2).^2 .* exp(-x.^2); 
psi_sq(4,:) = N(4)^2 * (8 * x.^3 - 12 * x).^2 .* exp(-x.^2); 
psi_sq(4,:) = N(4)^2 * (8 * x.^3 - 12 * x).^2 .* exp(-x.^2); 

% Compute |psi^2| values for high energy level
N_h = N(9);
psi_sq_high = N_h^2 * (256*x.^8-3584*x.^6+13440*x.^4-13440*x.^2+1680).^2 .* exp(-x.^2); 
%psi_sq_high = N_h^2 * (32*x.^5-160*x.^3+120*x).^2 .* exp(-x.^2); 
%psi_sq_high = N_h^2 * (64*x.^6-480*x.^4+720*x.^2-120).^2 .* exp(-x.^2); 
%psi_sq_high = N_h^2 * (128*x.^7-1334*x.^5+3360*x.^3-1680*x).^2 .* exp(-x.^2); 
%psi_sq_high = N_h^2 * (256*x.^8-3584*x.^6+13440*x.^4-13440*x.^2+1680).^2 .* exp(-x.^2); 

h=figure;
plot(x,psi_sq(1,:),'b');
hold on;
plot(x,psi_sq(2,:)+1,'r');
plot(x,psi_sq(3,:)+2,'g');
plot(x,psi_sq_high+3.5,'c');

% plot some reference lines
plot(x,0.0,'k');
plot(x,1,'k');
plot(x,2,'k');
plot(x,3,'k');
plot(x,3.5,'k');
plot(x,4.5,'k');

% fix up the axes
axis([-4.5 4.5 0 5]);
ca = get(h,'CurrentAxes');
yt = get(ca,'YTick');
ytl = get(ca,'YTickLabel');

% convert to characters
yt = [0 0.5 1 1.5 2 2.5 3.5 4];
ytl = char(ytl);

label_len = length(ytl);
for i=1:4
  ytl(2*i+1,:) = '0  ';
  ytl(2*i,:) = '0.5';
end
set(ca,'YTick',yt);
set(ca,'YTickLabel',ytl);

% add some "..." in between
text(0,3.2,'.')
text(0,3.25,'.')
text(0,3.3,'.')

text(0,4.7,'.')
text(0,4.75,'.')
text(0,4.8,'.')

title('Probability Densities for the First Few Quantum Harmonic Oscillator Energy Eigenfunctions');
