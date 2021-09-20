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
num_levels = 6;
n = 0:num_levels-1;

% The normalization constants for the wave functions.
N = 1./sqrt(2.^n .* cumprod(n+1)./(n+1)) * (alpha/pi)^(.25);

% position mesh
x = [-4:.01:4];

% Compute the |psi^2| values.
psi = zeros(num_levels, length(x));

psi(1,:) = N(1) * exp(-x.^2/2); 
psi(2,:) = N(2) * 2 * x .* exp(-x.^2/2); 
psi(3,:) = N(3) * (4 * x.^2 - 2) .* exp(-x.^2/2); 
psi(4,:) = N(4) * (8 * x.^3 - 12 * x) .* exp(-x.^2/2); 
psi(5,:) = N(5) * (16* x.^4 - 48 * x.^2 + 12) .* exp(-x.^2/2); 
psi(6,:) = N(6) * (32* x.^5 - 160 * x.^3 + 120 * x) .* exp(-x.^2/2); 

%figure(1);

%subplot(1,3,1);
plot(x,0.6*psi(1,:)+0.5,'b');
%hold on;
plot(x,0.6*psi(2,:)+1.5,'r');
%hold on;
plot(x,0.6*psi(3,:)+2.5,'g');

% plot the potential well
plot(x,x.^2/2,'m');


%figure(2);
subplot(1,3,1);
plot(x,0.6*psi(4,:)+0.5,'b');
hold on;
plot(x,0.6*psi(5,:)+1.5,'r');
hold on;
plot(x,0.6*psi(6,:)+2.5,'g');


% plot some reference lines
plot(x,0.5,'k');
plot(x,1.5,'k');
plot(x,2.5,'k');

hold on;
% plot the potential well
plot(x,x.^2/2,'m');

% fix the axes
axis([-2.5 2.5 0 2.25]);
title('\psi');

subplot(1,2,1);
plot(x,1.5*psi(1,:).*psi(1,:)+0.5,'b');
hold on;
plot(x,1.5*psi(2,:).*psi(2,:)+1.5,'r');
hold on;
plot(x,1.5*psi(3,:).*psi(3,:)+2.5,'g');

subplot(1,2,1)
plot(x,1.5*psi(4,:).*psi(4,:)+0.5,'b');
hold on;
plot(x,1.5*psi(5,:).*psi(5,:)+1.5,'r');
hold on;
plot(x,1.5*psi(6,:).*psi(6,:)+2.5,'g');

% plot some reference lines
plot(x,0.5,'k');
plot(x,1.5,'k');
plot(x,2.5,'k');

% plot the potential well
plot(x,x.^2/2,'m');

% fix the axes
axis([-2.5 2.5 0 2.25]);
title('|\psi|^2');

hold off;
