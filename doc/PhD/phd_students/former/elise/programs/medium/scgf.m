% Main script of self-consistent solving of in-medium self-energy Sigma and
% vertex function Gamma using Green's functions formalism

clear

% Declare which names are global 
globalvalues
% Constants and general data from script
scgfconstants
scgfdata


% Determine partial waves if first time scgf is executed, or if jmax
% and/or jmin has been changed. If file spinmatfile does not exist,
% run scripts globalvalues scgfconstants scgfdata generate_CG_table from
% command line. Returns matrix where each row contains 
% the set of spin numbers of each partial wave 
load spinmatfile partjmax partjmin spinmat no_chans
if (partjmax~=jmax | partjmin~=jmin)
    disp('Generating new set of partial waves')
    [spinmat,no_chans] = find_partialwaves;
end
%no_chans

% Set up mesh points and weights for center-of-mass
Kpoints = cm_mesh;
%Kpoints

% Set up mesh points and weights for relative momentum k.
[kpoints,kweights] = rel_mesh;

% Starting coefficients for the single-particle energy
%energy_coeff = [1./(2.*mass), 0]; %% NBNBNB!! TEST OM RIKTIG VERDI OVERALT
% SOM GLOBAL
energy_coeff = [0 1./(2.*mass), 0];
%energy_coeff = [-0.00000000018921   0.00087183713162 -43.97428228667606];

% Compute the initial coefficients for the single-particle energy
% polynomial and the energies corresponding to each momentum in qvec. 
% Initial guess free particles.
new_energies = sp_energy(qvec);
energies = zeros(1,length(qvec));
new_energy_coeff = energy_coeff;

iteration = 1;

%%Loop for self-consistently solving the energy spectrum

while sum(abs(new_energies-energies)>energytolerance)>0
    if iteration > maximum_iterations
        energyvec = sp_energy(qvec);
        figure(1)
        plot(qvec./sc_mom_factor,aSigma)        
        figure(2)
        plot(qvec./sc_mom_factor,areSigmaU)
        figure(3)
        plot(qvec./sc_mom_factor,areSigmaD)
        figure(4)
        plot(qvec./sc_mom_factor,aSigmaV)
        %gammaplot((2./pi).*aReDelta(2,:,1),Omegapoints,3,'omega')
        gammaplot(aGamma(1,5,:),kpoints,5,'k')
        gammaplot(aGamma(2,:,1),Omegapoints,6,'omega')
        figure(7)
        olden = (qvec.*qvec)./(2.*mass) + aSigma;
        newen = new_energy_coeff(1).*(qvec.^4) + new_energy_coeff(2).*(qvec.*qvec)+new_energy_coeff(3);
        plot(qvec./hbarc,olden,'b',qvec./hbarc,newen,'g')
        error('Maximum number of iterations allowed exceeded')
        break
    end
    
    % Initialize for new loop
    energy_coeff = new_energy_coeff
    energies = new_energies;
    aGamma = 0;
    Gamma = 0;
    aReDelta = 0;
    Sigma = zeros(1,length(qvec));
    aSigma = zeros(1,length(qvec));
    SigmaV = 0;
    reSigmaU = 0;
    areSigmaU = zeros(1,length(qvec));
    areSigmaD = zeros(1,length(qvec));
    aSigmaV = zeros(1,length(qvec));
    V = 0;
    e_fermi = sp_energy(k_fermi)
    penvec = sp_energy(sqrt(9.*(qvec.^2)));
    menvec = sp_energy(sqrt((qvec.^2)));
    maxenq = max(sp_energy(qvec));
    minenq = min(sp_energy(qvec));
    Omegamin = min([minenq+menvec,minenq+penvec,minenq+ minenq, ...
            -2.*e_fermi - (2.*abs(e_fermi)./2).*tan((pi./4.0).*(1+imgr1(end)))]);
    Omegamax = max([maxenq+menvec,maxenq+penvec,maxenq+ maxenq, ...
            abs(2.*e_fermi) + (2.*e_fermi./2).*tan((pi./4.0).*(1+imgr1(end)))]);
    OSmin = min([minenq+menvec,minenq+penvec,minenq+ minenq]);
    OSmax = max([maxenq+menvec,maxenq+penvec,maxenq+ maxenq]);
    Omegapoints = omega_mesh(Omegamin,Omegamax,OSmin,OSmax,2);
    if length(Omegapoints)~=n_omega_points
        error('Error in generating Omegapoints')
    end
    
    if average==1
        n_pw = size(spinmat,1)
        %% Gamma matrix using angle-averaged Pauli operator
        for pwcounter = 1:n_pw
            pwcounter
            spinvec = spinmat(pwcounter,:);
            V = potential(kpoints,kpoints,spinvec);
            %     % Compute Gamma(K,Omega,kr) (returned as a three-dimensional matrix),
            %     % Re(DeltaGamma) and the matrix element VL for each partial wave, then sum
            [apwgamma,apwV] = gammavertex(spinvec, ...
                Kpoints,Omegapoints,kpoints,kweights,V);
            %         s = spinmat(pwcounter,2);
            %         t = spinmat(pwcounter,5);
            %         Spinfactor = (1-(-1).^(s+t)).*(2.*s+1).*(2.*t+1);
            Spinfactor = 12;
            aGamma = aGamma + Spinfactor.*apwgamma;
            %gamk = squeeze(aGamma(1,1,:))
            aReDeltaL = real(apwgamma)-apwV;
            aReDelta = aReDelta + Spinfactor.*aReDeltaL;
            %gammaplot(apwV(1,5,:),kpoints,5,'k')
        end
        if ~any(aReDelta)
            error('Only zeros in ReDelta')
        end
        % %     % Find the real parts of DeltaGammaUp and DeltaGammaDown
        [aReDeltaD, aReDeltaU] = redelta(aReDelta,imag(aGamma),Omegapoints);
        %     Dredelk = squeeze(aReDeltaD(1,1,:))'
        %     Uredelk = squeeze(aReDeltaU(1,1,:))'
    end
    
    disp('Finished Gamma calculation, starting on self-energy Sigma')
    
    for qvec_ctr = 1:length(qvec)
        % Compute the self-energy Sigma
        q = qvec(qvec_ctr);
        [aSigmat, aSigmaVt, areSigmaDt areSigmaUt] = ...
            sigma(q,sp_energy(q),aReDeltaD,aReDeltaU,V,Kpoints,Omegapoints,kpoints);
        aSigma(qvec_ctr) = aSigmat;
        aSigmaV(qvec_ctr) = aSigmaVt;
        areSigmaU(qvec_ctr) = areSigmaUt;
        areSigmaD(qvec_ctr) = areSigmaDt;
    end
    
    % Find new single-particle energy
    [new_energies,new_energy_coeff] = newenergy(aSigma); % NBNBNB!!!!ENHETER!!
    iteration = iteration+1
end

energy_coeff = new_energy_coeff
energyvec = sp_energy(qvec);
figure(1)
plot(qvec./sc_mom_factor,aSigma)        
figure(2)
plot(qvec./sc_mom_factor,areSigmaU)
figure(3)
plot(qvec./sc_mom_factor,areSigmaD)
figure(4)
plot(qvec./sc_mom_factor,aSigmaV)
%gammaplot((2./pi).*aReDelta(2,:,1),Omegapoints,3,'omega')
gammaplot(aGamma(1,5,:),kpoints,5,'k')
gammaplot(aGamma(2,:,1),Omegapoints,6,'omega')

figure(8)
olden = (qvec.*qvec)./(2.*mass) + aSigma;
newen = energy_coeff(1).*(qvec.*qvec)+energy_coeff(2);
plot(qvec./hbarc,olden,'b',qvec./hbarc,newen,'g')