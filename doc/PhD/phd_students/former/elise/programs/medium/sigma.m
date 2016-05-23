%  Function computing the self energy Sigma

function [Sigma, SigmaV, reSigmaD, reSigmaU] = sigma(q,omega,ReDeltaD,ReDeltaU,V,Kpoints,Omegapoints,kpoints)
globalvalues
%warning off all
SigmaV = sigmaV(q,V);
reSigmaD = resigmaD(q,omega,ReDeltaD,Kpoints,Omegapoints,kpoints);
reSigmaU = resigmaU(q,omega,ReDeltaU,Kpoints,Omegapoints,kpoints);
Sigma = SigmaV+ reSigmaD+ reSigmaU;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function sigmaV
function SigmaV = sigmaV(q,V)
globalvalues
SigmaV=0;
%% NBNBNB!! Ikke ferdig enda!! Mangler andre muligheter enn V2 and V20
% Uniform map from begA = 0.5 *abs(k_fermi-q) to endB = 0.5 *(q+k_fermi)
if strcmp(pot,'V20') ==1
    begA = 0.5.*abs(k_fermi-q);
    endB = 0.5.*(q+k_fermi);
    kreexc = (begA+endB).*0.5 + (endB-begA).*0.5.*gr1';
    krewexc = (endB-begA).*0.5.*grw1';
    Vdiage = diag(potential(kreexc,kreexc,[0]));
    if q<k_fermi
     
        % Uniform map from begA = 0 to endB = 0.5 *(k_fermi-q)
        beA = 0;
        enB = 0.5.*(k_fermi-q);
        kre = (beA+enB).*0.5 + (enB-beA).*0.5.*gr1';
        krew = (enB-beA).*0.5.*grw1';
        Vdiag = diag(potential(kre,kre,0));
        SigmaV = 24.*(sum(kre.*kre.*krew.*Vdiag') + ...
            (1./(2.*q)).*sum(kreexc.*krewexc.*Vdiage'.* ...
            (q.*kreexc + 0.25.*(k_fermi.*k_fermi-q.*q) - kreexc.*kreexc)));
    else
        SigmaV = 24.*(1./(2.*q)).*sum(kreexc.*krewexc.*Vdiage'.* ...
            (q.*kreexc + 0.25.*(k_fermi.*k_fermi-q.*q) - kreexc.*kreexc));
    end
 elseif strcmp(pot,'V2')==1
     for potpart=1:4;
         const = V2potconst(potpart);
         alpha = V2alpha(potpart);
    pSigmaV = (8./(3.*pi)).*((k_fermi.^3)./(alpha.^2)) - ...
        k_fermi./pi + ...
        ((q.^2-k_fermi.^2-alpha.^2)./(4.*pi.*q)).*log((alpha.^2+(k_fermi+q).^2)./(alpha.^2+(k_fermi-q).^2)) + ...
        (alpha./pi).*(atan((k_fermi+q)./alpha) + atan((k_fermi-q)./alpha)); 
    %SigmaV = SigmaV+(1./hbarc).*(1./0.7).*(pi./2).*const.*pSigmaV;
    SigmaV = SigmaV+(1./hbarc).*(1./0.7).*const.*pSigmaV;
   
end
else
    warning('Unknown potential; SigmaV = 0')
    SigmaV = 0;
end
%disp('SigmaV')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function resigmaD
function SigmaD = resigmaD(q,omega,ReDeltaD,Kpoints,Omegapoints,kpoints)
globalvalues

if q<k_fermi
% Uniform map from begA = 0 to endB = 0.5 *(q+k_fermi)
    begA = 0;
    endB = 0.5.*(q+k_fermi);
    Pgrid = (begA+endB).*0.5 + (endB-begA).*0.5.*gr1';
    Pweights = (endB-begA).*0.5.*grw1';
else
    % Uniform map from begA = 0.5 *(q-k_fermi) to endB = 0.5 *(q+k_fermi)
    begA = 0.5.*(q-k_fermi);
    endB = 0.5.*(q+k_fermi);
    Pgrid = (begA+endB).*0.5 + (endB-begA).*0.5.*gr1';
    Pweights = (endB-begA).*0.5.*grw1';
end

Pint = 0;
for P = 1:length(Pgrid)
    if (Pgrid(P)>0)&(Pgrid(P)<=abs((k_fermi-q).*0.5))
        cosmin = -1;
    elseif (Pgrid(P)>abs(0.5.*(k_fermi-q)))&(Pgrid(P)<=(k_fermi+q))
        cosmin = ((q.*q)+(4.*Pgrid(P).*Pgrid(P))-k_fermi.*k_fermi)./(4.*q.*Pgrid(P));
    else  P
        error('P ligger i feil intervall')
    end
    cosmax = 1;
     % Uniform map from cosmin to cosmax
    cosgrid = (cosmin+cosmax).*0.5 + (cosmax-cosmin).*0.5.*gr1';
    cosweights = (cosmax-cosmin).*0.5.*grw1';
    
    %NBNBNB! vekter for P? TENK GJENNOM DETTE INTEGRALET!
    kprime = sqrt(q.^2+(4.*Pgrid(P).*Pgrid(P))-(4.*q.*Pgrid(P).*cosgrid));
    krel = sqrt(q.^2+(Pgrid(P).*Pgrid(P))-(2.*q.*Pgrid(P).*cosgrid));
    cosint = 0;
     interpoldelta = reDeltainterp2(krel,2.*Pgrid(P),omega+sp_energy(kprime), ...
             ReDeltaD,Kpoints,Omegapoints,kpoints);
         cosint = sum(Pgrid(P).^2.*cosweights.*interpoldelta);
    Pint = Pint + cosint.*Pweights(P);
end

if isnan(Pint)
       kprime
    krel
    Pgrid
    Pweights
    error('Pint ikke tall???')
end
%SigmaD = (1./((pi.^2))).*Pint;
SigmaD = Pint;
%disp('SigmaD')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function resigmaU
function SigmaU = resigmaU(q,omega,ReDeltaU,Kpoints,Omegapoints,kpoints)
globalvalues

if q<k_fermi
    %% OBS!!!! TRYKKFEIL I RAMOS???
    begA = 0.5.*(k_fermi-q);
    endB = k_fermi;
    Pgrid = (begA+endB).*0.5 + (endB-begA).*0.5.*gr1';
    Pweights = (endB-begA).*0.5.*grw1';
else
    % Uniform map from begA = 0  to endB = k_fermi
    begA = 0;
    endB = k_fermi;
    Pgrid = (begA+endB).*0.5 + (endB-begA).*0.5.*gr1';
    Pweights = (endB-begA).*0.5.*grw1'; 
end

Pint = 0;
for P = 1:length(Pgrid) 
    if q<k_fermi
        %% OBS!!!! TRYKKFEIL I RAMOS???
        if (Pgrid(P)>=0.5.*(k_fermi-q))&(Pgrid(P)<=(k_fermi+q).*0.5)
            cosmax = ((q.*q)+(4.*Pgrid(P).*Pgrid(P))-k_fermi.*k_fermi)./(4.*q.*Pgrid(P));
        elseif (Pgrid(P)>=abs(0.5.*(k_fermi+q)))&(Pgrid(P)<=k_fermi)
            cosmax = 1;
        else P
            Pgrid(P)
            error('P ligger i feil intervall 1')
        end
    else 
       if ((Pgrid(P)>=0) & (Pgrid(P)<=(-k_fermi+q).*0.5))
            cosmax = 1;
        elseif ((Pgrid(P)>=(0.5.*(-k_fermi+q))) & (Pgrid(P)<=k_fermi))
            cosmax = ((q.*q)+(4.*Pgrid(P).*Pgrid(P))-k_fermi.*k_fermi)./(4.*q.*Pgrid(P));
        else P
            Pgrid(P)
            error('P ligger i feil intervall 2')
        end 
    end
    cosmin= -1;
    % Uniform map from cosmin to cosmax
    cosgrid = (cosmin+cosmax).*0.5 + (cosmax-cosmin).*0.5.*gr1';
    cosweights = (cosmax-cosmin).*0.5.*grw1';
    
    %NBNBNB! vekter for P? TENK GJENNOM DETTE INTEGRALET! 
    kprime = sqrt(q.^2+(4.*Pgrid(P).*Pgrid(P))-(4.*q.*Pgrid(P).*cosgrid));
    krel = sqrt(q.^2+(Pgrid(P).*Pgrid(P))-(2.*q.*Pgrid(P).*cosgrid));
    cosint = 0;
    interpoldelta = reDeltainterp2(krel,2.*Pgrid(P),omega+sp_energy(kprime), ...
             ReDeltaU,Kpoints,Omegapoints,kpoints);
         cosint = sum(Pgrid(P).^2.*cosweights.*interpoldelta);
    Pint = Pint + cosint.*Pweights(P);
end

if isnan(Pint)
    kprime
    krel
    Pgrid
    Pweights
    error('Pint ikke tall???')
end
%SigmaU = -(1./((pi.^2))).*Pint;
SigmaU = -Pint;
%disp('SigmaU')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function reDeltainterp2. Interpolates to find ReDeltaD in the points ov vector krel,
% K, Omega
function redeltapoint = reDeltainterp2(krelvec,Kcm,omegavec,ReDelta,Kpoints,Omegapoints,kpoints)
globalvalues

if ((Kcm<=Kpoints(1)) | (Kcm>=Kpoints(n_cm_points)) |  ...
        (omegavec(1)<=Omegapoints(1)) | (omegavec(end)>=Omegapoints(n_omega_points)) | ... 
        (krelvec(1)<=kpoints(1)) | (krelvec(end)>=kpoints(n_rel_points)))
    redeltapoint = 0;
    Kcm
    Omegavec(n_omega_points)
    omegaen
    Omegavec(1)
    krela
    disp('redelta = 0 fordi K, Omega eller krel for stor eller for liten')
else
    Kcmvec = repmat(Kcm,1,length(krelvec));
    redeltapoint = interp3(Omegapoints,Kpoints,kpoints,ReDelta,omegavec,Kcmvec,krelvec);
    %size(redeltapoint)
end