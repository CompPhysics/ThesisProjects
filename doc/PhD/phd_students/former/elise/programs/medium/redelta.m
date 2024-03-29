% Function which computes the real parts of DeltaGammaUp and DeltaGammaDown

% function [ReDeltaD,ReDeltaU]= redelta(ReDelta,ImGammacoeff,Omegapoints)
function [ReDeltaD,ReDeltaU]= redelta(ReDelta,ImGamma,Omegapoints)
globalvalues
% Initialize
ReDeltaU = zeros(n_cm_points,n_omega_points,n_rel_points);
ReDeltaD = zeros(n_cm_points,n_omega_points,n_rel_points);

twoefermi = 2*sp_energy(k_fermi)
for K = 1:n_cm_points
    for omega = 1:n_omega_points
        Omega = Omegapoints(omega);
        Ompoints = Omegapoints;
        % Set up grid for integration over energy, midpoint midpointa+4*e_fermi
        midpointa = Omega;
        if Omega<twoefermi
            largeromegagrid = twoefermi + (abs(twoefermi)./2).*tan((pi./4.0).*(1+imgr1'));
            largeromegaweights = ((abs(twoefermi)./2).*((pi./4.0).*imgrw1'))./ ... 
                (cos((pi./4.0).*(1+imgr1')).*(cos((pi./4.0).*(1+imgr1'))));
            
            for k = 1:n_rel_points
                %% NBNB! Bare snarvei foreloepig!!
                
                ImGam = squeeze(ImGamma(K,:,k));
                
                omintegral = -(1./pi).*sum(((imgammapol(Ompoints,ImGam,largeromegagrid))...
                    ./(Omega-largeromegagrid)).*largeromegaweights);
                ReDeltaU(K,omega,k) = ReDelta(K,omega,k) - omintegral;
                ReDeltaD(K,omega,k) = omintegral;
                if ReDeltaU(K,omega,k)==0
                    a = ReDelta(K,omega,k)
                    su = (1./pi).*sum(((imgammapol(Ompoints,ImGam,largeromegagrid))...
                        ./(Omega-largeromegagrid)).*largeromegaweights)
                end
                
            end
            
        elseif Omega>twoefermi
            
            omgr = -twoefermi + (abs(twoefermi)./2).*tan((pi./4.0).*(1+imgr1'));
            omw = ((abs(twoefermi)./2).*((pi./4.0).*imgrw1'))./ ... 
                (cos((pi./4.0).*(1+imgr1')).*(cos((pi./4.0).*(1+imgr1'))));
            
            lessomegagrid = -omgr(1,length(omgr):-1:1);
            lessomegaweights = omw(1,length(omw):-1:1);
            
            for k = 1:n_rel_points
                %% NBNB! Bare snarvei foreloepig!!
                
                ImGam = squeeze(ImGamma(K,:,k));
                
                omintegral =  (1./pi).*sum(((imgammapol(Ompoints,ImGam,lessomegagrid))...
                    ./(Omega-lessomegagrid)).*lessomegaweights);
                ReDeltaD(K,omega,k) = ReDelta(K,omega,k) - omintegral;
                ReDeltaU(K,omega,k) =  omintegral; 
                if ReDeltaD(K,omega,k)==0
                    a = ReDelta(K,omega,k)
                    su = sum(((imgammapol(Ompoints,ImGam,largeromegagrid))...
                        ./(Omega-largeromegagrid)).*largeromegaweights)
                end
            end      
        else
            error('Energien er ikke fornuftig')
        end
    end
end

% if ~any(ReDeltaD)
%     error('Only zeros in ReDeltaD')
% end
% if ~any(ReDeltaU)
%     error('Only zeros in ReDeltaU')
% end
if isnan(ReDeltaD)
    error('NaN in ReDeltaD')
end
if isnan(ReDeltaU)
    error('NaN in ReDeltaU')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function imgammapol, determining the value of ImGamma on the points of
% Omegagrid. Must be changed when degree of polynomial is changed
function imgammapoints = imgammapol(Ompoints,imGamma,Omegagrid)
globalvalues

imgammapoints = interp1(Ompoints,imGamma,Omegagrid,'linear','extrap');
% if imgammapoints == 0
%        imGamma
%     end
