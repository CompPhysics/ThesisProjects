% Function newenergy, determining the coefficients of the energy polynomial

function [new_energies,coeff] = newenergy(Sigma)
globalvalues

new_energies = (qvec.*qvec)./(2.*mass) + Sigma;

%coeff = lsqcurvefit(@energyfcn,energy_coeff,qvec,new_energies)
coeff = polyfit(qvec.*qvec,new_energies,2)
coefftest = polyfit(qvec,new_energies,4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function energyfcn, must have the form of the sp-energy
function F = energyfcn(x,xdata)
F = x(1).*xdata.*xdata + x(2);