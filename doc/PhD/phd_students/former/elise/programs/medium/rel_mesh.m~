% Determining the center-of-mass gridpoints and weights, with appropriate
% mapping to integration interval. Uses external C function gauleg.

function [kpoints,kweights] = rel_mesh;
globalvalues 
if n_rel_points == 1;
    kpoints = 0.567.*sc_mom_factor;
    kweights = 0;
    elseif n_rel_points == 2;
    kpoints = [-1000,1000];
    kweights = [1,1];
else
    
    % % Create mesh
    [x1,w1] = gauleg(-1.0,1.0,n_rel_points);
    % Map to correct interval
    % x = hbarc.*tan((pi./4.0).*(1+x1));
    % w = hbarc.*((pi./4.0).*w1)./ ... 
    %     (cos((pi./4.0).*(1+x1)).*(cos((pi./4.0).*(1+x1))));
    %     
 
    x = sc_mom_factor.*tan((pi./4.0).*(1+x1));
    w = (sc_mom_factor.*((pi./4.0).*w1))./ ... 
        (cos((pi./4.0).*(1+x1)).*(cos((pi./4.0).*(1+x1))));
    % A = energy_coeff(2);
    % B = energy_coeff(1);
    % mpkt = sqrt(((Omegapoints-2.*A)./(2.*B))-((Kpoints.^2).*0.25))
    % x = mpkt.*((1+x1)./(1-x1));
    % w = (2.*mpkt.*w1)./((1-x1).^2);
    
    kpoints = x';
    kweights = w';
    
end