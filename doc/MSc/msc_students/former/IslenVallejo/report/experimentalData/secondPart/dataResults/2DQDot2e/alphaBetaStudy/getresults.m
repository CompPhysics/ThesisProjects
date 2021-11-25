function output = getresults(name)

% Run results script:
eval(name)
%%%importdata(name)
  
% Now, E, energy_cut, blabla, gets assigned,
% but only *locally* within this function. Upon exit,
% these are destroyed.

format long e;

save values.mat;

%  output.nsd = nsd;
%  output.nel = nel;
%  output.potpar = potpar;
%  output.nmc = nmc;
%  output.nth = nth;
%  output.nvp = nvp;
%  output.alpha = alpha;
%  %  output.beta = beta;
%  output.percentAccept = percentAccept;
%  output.energy = energy;
%  output.variance = variance;
%  %  output.error = error;
