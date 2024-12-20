% Determining the energy gridpoints and weights, with appropriate
% mapping to integration interval. Uses external C function gauleg.

function Omegapoints = omega_mesh(Omegamin,Omegamax,OSmin,OSmax,dis);
globalvalues
if n_omega_points == 1;
    Omegapoints = 173.3;
    %Omegapoints = 10.*sc_energy_factor;
    %Omegapoints = 2.*sp_energy(k_fermi)
    %Omegapoints = 255.9680
    Omegaweights = 0;
elseif n_omega_points == 2;
    Omegapoints = [-100000,100000];
    Omegaweights = [1,1];
else
    %Omegamax = 10000;
    %Omegamin = -10000;
    if (Omegamin<OSmin) & (OSmin<OSmax) & (OSmax<Omegamax)
        Omegamin = Omegamin-0.1.*abs(Omegamin);
        Omegamax = (Omegamax+0.1.*abs(Omegamax));
        innerpoints = round(n_omega_points./dis);
        outerpoints = n_omega_points-innerpoints;
        smallpoints = round(0.5.*outerpoints);
        largepoints = outerpoints-smallpoints;
        largeOstep = abs((Omegamax-OSmax))./(largepoints-1);
        innerstep = abs((OSmax-OSmin))./(innerpoints+1);
        smallOstep = abs(OSmin-Omegamin)./(smallpoints-1);
        smallO = Omegamin:smallOstep:OSmin;
        innerO = OSmin+innerstep:innerstep:OSmax-innerstep;
        largeO = OSmax:largeOstep:Omegamax;
        Omegapoints = [smallO, innerO, largeO];
    else
        Omx = max(OSmax,Omegamax);
        Omi = min(Omegamin,OSmin);
        Omi = (Omi-0.1.*abs(Omi));
        Omx = (Omx+0.1.*abs(Omx));
        step = abs((Omx-Omi))./(n_omega_points-1);
        Omegapoints = Omi:step:Omx;
    end
    
    %length(Omegapoints)
    %     % Create mesh
    %     [x1,w1] = gauleg(-1.0,1.0,n_omega_points);
    %     % Map to correct interval
    % %     x = sc_energy_factor.*tan((pi./4.0).*(1+x1));
    % %     w = sc_energy_factor.*((pi./4.0).*w1)./ ... 
    % %         (cos((pi./4.0).*(1+x1)).*(cos((pi./4.0).*(1+x1))));
    %     
    %     
    %         x = 1.*tan((pi./4.0).*(1+x1));
    %         w = 1.*((pi./4.0).*w1)./ ... 
    %             (cos((pi./4.0).*(1+x1)).*(cos((pi./4.0).*(1+x1))));
    %      
    %     Omegapoints = x';
    %     Omegaweights = w';
end