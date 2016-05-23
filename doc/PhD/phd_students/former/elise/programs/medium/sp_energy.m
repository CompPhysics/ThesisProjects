% Function sp_energy

function energies = sp_energy(q)
globalvalues

% 
% A = energy_coeff(2);
% B = energy_coeff(1);
% energies = A + B.*q.^2;

A = energy_coeff(3);
B = energy_coeff(2);
C = energy_coeff(1);

energies = A + B.*q.^2 + C.*q.^4;
