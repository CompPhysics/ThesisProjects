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
num_levels = 4;
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
%psi(5,:) = N(5) * (16* x.^4 - 48 * x.^2 + 12) .* exp(-x.^2/2); 
%psi(6,:) = N(6) * (32* x.^5 - 160 * x.^3 + 120 * x) .* exp(-x.^2/2); 

figure(1);

subplot(1,2,1);
plot(x,0.6*psi(1,:)+0.5,'k','LineWidth',2);
hold on;
plot(x,0.6*psi(2,:)+1.5,'b','LineWidth',2);
hold on;
plot(x,0.6*psi(3,:)+2.5,'r','LineWidth',2);
hold on;
plot(x,0.6*psi(4,:)+3.5,'g','LineWidth',2);
hold on;
axis([0 4 0 4])

% plot the potential well
plot(x,x.^2/2,'m','LineWidth',2);


%figure(2);
%subplot(1,2,1);

%plot(x,0.6*psi(4,:)+0.5,'b');
%hold on;
%plot(x,0.6*psi(5,:)+1.5,'r');
%hold on;
%plot(x,0.6*psi(6,:)+2.5,'g');
%hold on;

% plot some reference lines
plot(x,0.5,'k','LineWidth',2);
hold on;
plot(x,1.5,'k','LineWidth',2);
hold on;
plot(x,2.5,'k','LineWidth',2);
hold on;
plot(x,3.5,'k','LineWidth',2);

hold on;
% plot the potential well
%plot(x,x.^2/2,'m');

% fix the axes
axis([-4.5 4.5 0 4.25]);
title('\psi');

subplot(1,2,2);
plot(x,1.5*psi(1,:).*psi(1,:)+0.5,'k', 'LineWidth',2);
hold on;
plot(x,1.5*psi(2,:).*psi(2,:)+1.5,'b','LineWidth',2);
hold on;
plot(x,1.5*psi(3,:).*psi(3,:)+2.5,'r','LineWidth',2);
hold on;
plot(x,1.5*psi(4,:).*psi(4,:)+3.5,'g','LineWidth',2);
%subplot(1,2,1)

%hold on;
%plot(x,1.5*psi(5,:).*psi(5,:)+1.5,'r');
%hold on;
%plot(x,1.5*psi(6,:).*psi(6,:)+2.5,'g');

% plot some reference lines
hold on;
plot(x,0.5,'k');hold on;
plot(x,1.5,'k');hold on;
plot(x,2.5,'k');hold on;
plot(x,3.5,'k');hold on;

% plot the potential well
plot(x,x.^2/2,'m','LineWidth',2);

% fix the axes
axis([-4.5 4.5 0 4.25]);
title('|\psi|^2');

hold off;
