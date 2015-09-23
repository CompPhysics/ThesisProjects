% Determining the center-of-mass gridpoints and weights, with appropriate
% mapping to integration interval. Uses external C function gauleg.

function Kpoints = cm_mesh;
globalvalues

if n_cm_points == 1;
    %Kpoints = 1.8.*sc_mom_factor;
    %(0.5.*Kpoints.*0.5.*Kpoints)./mass
    %Kpoints = 1.8.*sc_mom_factor;
    %Kpoints = 2;
    Kpoints = 0.017.*sc_mom_factor;
    %Kpoints = 0;
elseif n_cm_points == 2;
    Kpoints = [-1000,1000];
else
    Kmax = 2.*(qvec(end)+k_fermi);
    step = Kmax./(n_cm_points-1);
    Kpoints = 0:step:Kmax;
    %length(Kpoints)
end