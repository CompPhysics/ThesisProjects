% Function computing the Gamma vertex

function [pwgamma, Vl] = gammavertex(spinvec,Kpoints,Omegapoints,kpoints,kweights,Vkernel)

globalvalues

kpn = n_rel_points+1;
for K=1:n_cm_points
    for Om=1:n_omega_points
        Omega = Omegapoints(Om);
        Kcm = Kpoints(K);
        kvec = kpoints;
        kwvec = kweights;
        kp = polemomentum(Omegapoints(Om),Kpoints(K),kvec(end));
        % Test of polemomentum
        %kptest = Energy(Kpoints(K),kp) - Omegapoints(K,O)
        if (~isreal(kp)|kp==0)
            m = length(kvec);
            V = Vkernel;
            if ~isempty(find(isnan(V)))
                V
                error('NaN in V')
            end
       
            for idx=1:m
                                u(idx) = ((kwvec(idx).*kvec(idx).*kvec(idx)).*qpaulibarminus(K,kvec(idx))) ...
                                    ./(Omega - Energy(Kcm,kvec(idx)));
                
            end
            for ai=1:m
                Amat(ai,1:m) = -V(ai,1:m).*u;
            end
            Amatrix = Amat + eye(m);

            % Find the R matrix by left-dividing the V matrix by the A matrix
            % (faster than inverting the A matrix and then multiply)
            R = Amatrix\V;

            for k=1:n_rel_points
                pwgamma(K,Om,k) = R(k,k);
                Vl(K,Om,k) = V(k,k);
            end 
        else 
            [R,V] = rmatrix(kvec,kwvec,Kcm,Omega,kp,spinvec,Vkernel); 
            A = afactor(Kcm,kp);
%             if ((K==1) & (Om==1))
%             R(1,kpn)
%             R(kpn,1)
%             R(kpn,kpn)
%         end
            denom = 1 + (A.^2).*(R(kpn,kpn).^2);
            for k=1:n_rel_points
                pwgamma(K,Om,k) = R(k,k) ...
                    - ((A.^2).*(R(k,kpn).*R(kpn,k).*R(kpn,kpn)))./denom ...
                    - (i.*A.*R(k,kpn).*R(kpn,k))./denom;
          
                Vl(K,Om,k) = V(k,k);
            end
        end
        
        
        if ~isempty(find(isnan(pwgamma)))
            error('NaN in pwgamma')
        end
        if ~any(pwgamma)
            pwgamma
            error('Only zeros in pwgamma')
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function rmatrix returns the real reaction matrix, obtained by matrix inversion
% method 
function [R,V] = rmatrix(kv,kwv,K,Omega,kp,spinvec,Vkernel)
globalvalues
kvec = 0;
kwvec = kwv;
kvec = kv;
if all(kvec-kp)    
    k = [kvec, kp];
else error('AU')
end
m = length(k);

% Initialize
V = zeros(m,m);

% Find last matrix elements V_l(k1,k2)
%Vtest = potential(k,k,spinvec);
V(1:(m-1),1:(m-1)) = Vkernel;
V(m,1:m) = potential(kp,k,spinvec);
V(1:m,m) = potential(k,kp,spinvec);
% if (Vtest-V)~=zeros(m,m)
%     error('V gal')
% end
if ~isempty(find(isnan(V)))
    V
    error('NaN in V')
end

Amatrix = zeros(m,m);
Amat = zeros(m,m);
R = zeros(m,m);
u = zeros(1,m);

% Construct the first m points of the u vector
% NBNBNB!! (2./pi)-faktorer!! 
for idx=1:m-1
    u(idx) = ((kwvec(idx).*kvec(idx).*kvec(idx)).*qpaulibarminus(K,kvec(idx))) ...
        ./(Omega - Energy(K,kvec(idx)));
end
% Add the last point
t = sum((kwvec)./((kvec.*kvec)-(kp.*kp)));
u(m) = -Blim(K,kp).*t;

% Construct the A matrix
for ai=1:m
    Amat(ai,1:m) = -V(ai,1:m).*u;
end
Amatrix = Amat + eye(m);
 
% Find the R matrix by left-dividing the V matrix by the A matrix
% (faster than inverting the A matrix and then multiply)
R = Amatrix\V;

if ~isempty(find(isnan(R)))
    Amatrix
    V
    error('NaN in Rmatrix')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function polemomentum: Finds the relative momentum corresponding to the 
% energy Omega
function  kp = polemomentum(Omega,K,kmax)
globalvalues
% A = energy_coeff(2);
% B = energy_coeff(1);
% A = energy_coeff(3);
% B = energy_coeff(2);
% kp = sqrt(((Omega-2.*A)./(2.*B))-((K.^2).*0.25))
A = energy_coeff(3);
B = energy_coeff(2);
C = energy_coeff(1);
if C==0
    kp = sqrt(((Omega-2.*A)./(2.*B))-((K.^2).*0.25));
else
    if ((Energyminusomega(0,K,Omega)<0) & Energyminusomega(kmax.*2,K,Omega)>0) | ...
            ((Energyminusomega(0,K,Omega)>0) & Energyminusomega(kmax.*2,K,Omega)<0)
        options = optimset;
        kp1 = fzero(@Energyminusomega,[0 kmax.*2],options,K,Omega);
        kp = kp1;
        if (abs((Omega-Energy(K,kp)))>abs(0.1.*Omega)) & isreal(kp)
%             K
%             Omega
%             diff = abs((Omega-Energy(K,kp)))
            disp('kp not good')
        end
    else
        kp = i.*pi;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function afactor
function A = afactor(K,q)
A = (pi.*(q.^2).*qpaulibarplus(K,q))./(abs(Eder(K,q)));
%A = (2.*(q.^2).*qpaulibarplus(K,q))./(abs(Eder(K,q)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Energy, the energy as function of center-of-mass and relative momentum 
function en = Energy(K,q)
globalvalues
% A = energy_coeff(2);
% B = energy_coeff(1);
% en = 2.*A + 2.*B.*((K.^2).*0.25 + (q.^2));
A = energy_coeff(3);
B = energy_coeff(2);
C = energy_coeff(1);
en = 2.*A + 2.*B.*((K.^2).*0.25 + (q.^2)) +  ...
    2.*C.*(((K.^2).*0.25 + (q.^2)).^2 + ((K.^2).*(q.^2).*qpaulibarminus(K,q).^3)./3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Energyminusomega, the energy as function of center-of-mass and relative momentum 
function en = Energyminusomega(q,K,Omega)
globalvalues
A = energy_coeff(3);
B = energy_coeff(2);
C = energy_coeff(1);
en = 2.*A + 2.*B.*((K.^2).*0.25 + (q.^2)) +  ...
    2.*C.*(((K.^2).*0.25 + (q.^2)).^2 + ((K.^2).*(q.^2).*qpaulibarminuskvec(K,q).^3)./3)-Omega;
% en = (2.*C./3).*(q.^6) + 2.*C.*K.*(q.^5) + (0.5.*C.*(K.^2) - 2.*C.*(k_fermi.^2)).*(q.^4) + ...
%     (2.*B.*K + C.*(K.^3)).*(q.^3) + ((C./8).*(K.^4) + 2.*(k_fermi.^4) - C.*(K.^2).*(k_fermi.^2)).*(q.^2)  + ...
%     (2.*A.*K + 0.5.*B.*(K.^2) + (C./8).*(K.^5) - K.*Omega).*q + ...
%     ((C./96).*(K.^6) - (2./3).*C.*(k_fermi.^6) - (C./8).*(K.^4).*(k_fermi.^2) + 0.5.*C.*(K.^2).*(k_fermi.^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Eder, the derivative of the energy 
function der = Eder(K,q)
globalvalues
%der = 4.*energy_coeff(1).*q;
B = energy_coeff(2);
C = energy_coeff(1);
der = 4.*B.*q + 8.*C.*q.*((K.^2).*0.25 + (q.^2)) + (4./3).*((K.^2).*q.*C.*qpaulibarminus(K,q).^3) + ...
    2.*C.*(K.^2).*(q.^2).*qpaulibarminus(K,q).^2.*qpaulibarder(K,q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function qpaulibarplus, the weighted, angleindependent Pauli Operator
% Q_pp+Q_hh
function Q = qpaulibarplus(K,q)
globalvalues
%Q = 1;
if K./2 <= k_fermi
    if ((q >= 0) & (q <= (k_fermi- (K.*0.5))))
        Q = 1;
    elseif ((q > (k_fermi- (K.*0.5))) & (q <= sqrt((k_fermi.^2)- ((K.^2).*0.25))))
        Q = ((k_fermi.^2) - ((K.^2).*0.25) - (q.^2)) ./ (K.*q);
    elseif ((q > sqrt((k_fermi.^2)- ((K.^2).*0.25))) & (q <= (k_fermi+ (K.*0.5))))
        Q = ((q.^2) + ((K.^2).*0.25) - (k_fermi.^2)) ./ (K.*q);
    elseif (q > (k_fermi+ (K.*0.5)))
        Q = 1;
    else error('Paulioperator sint')
    end
elseif K./2 > k_fermi
    if ((q >= 0) & (q < ((K.*0.5)-k_fermi)))
        Q = 1;
    elseif ((q >= (K.*0.5)-k_fermi) & (q < (k_fermi+(K.*0.5))))
        Q = ((q.^2) + ((K.^2).*0.25) - (k_fermi.^2)) ./ (K.*q);
    elseif (q >= (k_fermi+ (K.*0.5)))
        Q = 1;
    else error('Paulioperator sint')
    end
else error('Center-of-mass definitely wrong')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function qpaulibarminuskvec, the weighted, angleindependent Pauli Operator
% Q_pp-Q_hh taking a vector of k-values as input
function Qvec = qpaulibarminuskvec(K,kvec)
globalvalues
%Q = 1;
for ctr=1:length(kvec)
    q = kvec(ctr);
    if q<0
        Q(ctr) = 0
        disp('q mindre enn 0!!')
       
    elseif K./2 <= k_fermi
    if ((q >= 0) & (q < (k_fermi- (K.*0.5))))
        Q = -1;
    elseif ((q >= (k_fermi- (K.*0.5))) & (q < sqrt((k_fermi.^2)- ((K.^2).*0.25))))
        Q = -((k_fermi.^2) - ((K.^2).*0.25) - (q.^2)) ./ (K.*q);
    elseif ((q >= sqrt((k_fermi.^2)- ((K.^2).*0.25))) & (q < (k_fermi+ (K.*0.5))))
        Q = ((q.^2) + ((K.^2).*0.25) - (k_fermi.^2)) ./ (K.*q);
    elseif (q >= (k_fermi+ (K.*0.5)))
        Q = 1;
    else error('Paulioperator sint 1')
    end
elseif K./2 > k_fermi
    if ((q >= 0) & (q < ((K.*0.5)-k_fermi)))
        Q = -1;
    elseif ((q >= (K.*0.5)-k_fermi) & (q < (k_fermi+(K.*0.5))))
        Q = ((q.^2) + ((K.^2).*0.25) - (k_fermi.^2)) ./ (K.*q);
    elseif (q >= (k_fermi+ (K.*0.5)))
        Q = 1;
    else error('Paulioperator sint 2')
    end
 
else error('Center-of-mass definitely wrong')
end
Qvec(ctr) = Q;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function qpaulibarminus, the weighted, angleindependent Pauli Operator
% Q_pp-Q_hh
function Q = qpaulibarminus(K,q)
globalvalues
%Q = 1;
if K./2 <= k_fermi
    if ((q >= 0) & (q < (k_fermi- (K.*0.5))))
        Q = -1;
    elseif ((q >= (k_fermi- (K.*0.5))) & (q < sqrt((k_fermi.^2)- ((K.^2).*0.25))))
        Q = -((k_fermi.^2) - ((K.^2).*0.25) - (q.^2)) ./ (K.*q);
    elseif ((q >= sqrt((k_fermi.^2)- ((K.^2).*0.25))) & (q < (k_fermi+ (K.*0.5))))
        Q = ((q.^2) + ((K.^2).*0.25) - (k_fermi.^2)) ./ (K.*q);
    elseif (q >= (k_fermi+ (K.*0.5)))
        Q = 1;
    else error('Paulioperator sint 1')
    end
elseif K./2 > k_fermi
    if ((q >= 0) & (q < ((K.*0.5)-k_fermi)))
        Q = -1;
    elseif ((q >= (K.*0.5)-k_fermi) & (q < (k_fermi+(K.*0.5))))
        Q = ((q.^2) + ((K.^2).*0.25) - (k_fermi.^2)) ./ (K.*q);
    elseif (q >= (k_fermi+ (K.*0.5)))
        Q = 1;
    else error('Paulioperator sint 2')
    end
else error('Center-of-mass definitely wrong')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function qpaulibarder, the derivative of the weighted, angleindependent Pauli Operator
% Q_pp-Q_hh
function Q = qpaulibarder(K,q)
globalvalues
%Q = 1;
if K./2 <= k_fermi
    if ((q >= 0) & (q < (k_fermi- (K.*0.5))))
        Q = 0;
    elseif ((q >= (k_fermi- (K.*0.5))) & (q < sqrt((k_fermi.^2)- ((K.^2).*0.25))))
        Q = -(((K.^2).*0.25) -(k_fermi.^2) - (q.^2) )./ (K.*q.^2);
    elseif ((q >= sqrt((k_fermi.^2)- ((K.^2).*0.25))) & (q < (k_fermi+ (K.*0.5))))
        Q = ((q.^2) - ((K.^2).*0.25) + (k_fermi.^2)) ./ (K.*q.^2);
    elseif (q >= (k_fermi+ (K.*0.5)))
        Q = 0;
    else error('Paulioperator sint')
    end
elseif K./2 > k_fermi
    if ((q >= 0) & (q < ((K.*0.5)-k_fermi)))
        Q = 0;
    elseif ((q >= (K.*0.5)-k_fermi) & (q < (k_fermi+(K.*0.5))))
        Q = ((q.^2) - ((K.^2).*0.25) + (k_fermi.^2)) ./ (K.*q.^2);
    elseif (q >= (k_fermi+ (K.*0.5)))
        Q = 0;
    else error('Derivert Paulioperator sint')
    end
else error('Center-of-mass definitely wrong')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Blim, the factor converting from denominator Omega-Energy to
% k.^2-kp.^2
function b = Blim(K,kp)
globalvalues
%b = (kp.^2).*qpaulibarminus(K,kp).*(-1./(2.*energy_coeff(1)));
B = energy_coeff(2);
C = energy_coeff(1);
b = (kp.^2).*qpaulibarminus(K,kp).*(-1./(2.*B + C.*(K.^2)+4.*C.*(kp.^2)+(2./3).*C.*(K.^2).*(qpaulibarminus(K,kp).^3)));
   