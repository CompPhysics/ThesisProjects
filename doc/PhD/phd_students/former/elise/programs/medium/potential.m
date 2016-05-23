% Function potential, providing the potential V, matrix elements being
% V(k1,k2) in momentum space for different potentials

function [V, coupledV] = potential(k1,k2,spinvec)
globalvalues

m = length(k1);
n = length(k2);
if strcmp(pot,'Argonne') == 1
    % Construct the V matrix, containing the bessel transform of the 
% potential at all mesh points
    for ki1 = 1:m 
        for ki2 = 1:n 
            vkpw = (2/pi).*besseltransform(k1(ki1),k2(ki2), ...
                spinvec,pot,gausspoints,gaussweights);
            V(ki1,ki2) = vkpw(1,1); % NBNBNB! Hva med overlappende pw'er???
        end
    end
    
% elseif strcmp(pot,'Malfliet-Tjon') == 1
%     for ki1 = 1:m 
%         for ki2 = 1:m 
%             V(ki1,ki2) = 42;
%         end
%     end
%     

elseif strcmp(pot,'V20') == 1
    %tidc = cputime;
    for ki1 = 1:m 
        for ki2 = 1:n 
            V(ki1,ki2) = 0;
            x = (V2alpha.^2+k1(ki1).^2+k2(ki2).^2)./(2.*k1(ki1).*k2(ki2));
            V(ki1,ki2) = sum((1./hbarc).*(1./0.7).*(2./pi).* ...
                V2potconst.*(1./(2.*k1(ki1).*k2(ki2))).*0.5.*real(log((1+x)./(1-x))));
%             for ctr=1:4
%                 x = (V2alpha(ctr).^2+k1(ki1).^2+k2(ki2).^2)./(2.*k1(ki1).*k2(ki2));
%                 vv = (2./pi).*V2potconst(ctr).*(1./(2.*k1(ki1).*k2(ki2))).*0.5.*log((1+x)./(1-x));
%                 V(ki1,ki2) = V(ki1,ki2) + (1./hbarc).*(1./0.7).*real(vv);
%             end
        end
    end
    %vtid = cputime-tidc
    
elseif (strcmp(pot,'V2')==1) % | (strcmp(pot,'V20') == 1)
        %tidc = cputime;
        %Vv = (2./pi).*V2pot(spinvec(1),k1,k2,m,n,V2potconst,V2alpha);
        Vv = V2pot(spinvec(1),k1,k2,m,n,V2potconst,V2alpha);
        %V = (1./hbarc).*(1./0.7).*(pi./2).*Vv(1:m,1:n);
        V = (1./hbarc).*(1./0.7).*Vv(1:m,1:n);
        coupledV = [];
        if ~isempty(find(isnan(V)))
            error('Potential gives NaN')
        end
        %vtid = cputime-tidc
        %tidv = cputime;
%        for ki1 = 1:m 
%         for ki2 = 1:m
%             v = 0;
%             for teller = 1:length(V2potconst)
%                 % vkk is in units of 1/Mev^2 (0.7 is in units of 1./fm)
%                 vkk = V2potconst(teller).*(1./hbarc).*(1./0.7).* ...
%                     yukawa(spinvec(1),k(ki1),k(ki2),V2alpha(teller));
%                 v = v + vkk;
%             end
%             V(ki1,ki2) = v;
%         end
%     end 
%     vtid = cputime-tidv
elseif strcmp(pot,'box') == 1
    V0 = box_V0;
    a = box_a;
     %% NBNBNB!!! MAA REPROGRAMMERES TIL K1,K2-FORMAT!!!
     disp('Only works for k1==k2 for the moment')
    for ki1 = 1:m 
        for ki2 = 1:m 
            if (ki2==ki1)
                V(ki1,ki1) = V0.*(1./(2.*k(ki1).*k(ki1))).* ...
                    ((sin(2.*k(ki1).*a)./(2.*k(ki1)))-a);
            else
                V(ki1,ki2) = (V0./(2.*(k(ki1).*k(ki2)))).* ...
                    ((sin((k(ki1)+k(ki2)).*a)./(k(ki1)+k(ki2))) ...
                    - (sin((k(ki1)-k(ki2)).*a)./(k(ki1)-k(ki2))));
            end
        end
    end
elseif strcmp(pot,'Yamaguchi') == 1
    for ki1 = 1:m 
        for ki2 = 1:n 
            V(ki1,ki2) = -Yamaguchilambda.*yamag0(k1(ki1)).*yamag0(k2(ki2));
        end
    end
    elseif strcmp(pot,'cdbonn') == 1
        l = spinvec(1);
        s = spinvec(2);
        j = spinvec(3);
        tz = spinvec(6);
        if (l==j)
            coup = 0;
            coupledV = [];
        else
            coup = 2;
        end
    for ki1 = 1:m 
        for ki2 = 1:n 
            v = cdbonnpot([j s tz coup],[k1(ki1) k2(ki2)],hbarc);
            %if ((l==j) |  (j==0))
            if (l==j)
                % No coupled channels
                V(ki1,ki2) = v(1,1);
            elseif (l == j-1)
                % l=j-1, coupled to j+1. 
                V(ki1,ki2) = v(1,1);
                coupledV(ki1,ki2,1) = v(1,2);
                coupledV(ki1,ki2,2) = v(2,1);
                coupledV(ki1,ki2,3) = v(2,2);
            elseif (l == j+1)
                % l=j+1, coupled to j-1.
                V(ki1,ki2) = v(2,2);
                coupledV(ki1,ki2,1) = v(2,1);
                coupledV(ki1,ki2,2) = v(1,2);
                coupledV(ki1,ki2,3) = v(2,2);
            else
                error('Noe er galt med spinmat!!')
            end
            % The matrix coupledV has <l|V|l'> as first component, l'=l+2
            % or l'=l-2. Second component <l'|V|l>, and third <l'|V|l'>
        end
    end

else error('Unknown potential')
end
% % Test of potential.
% vpot = diag(V)'
% figure(1)
% plot(k(1:m-1)./hbarc,vpot(1:m-1).*hbarc.^3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Yukawa, computing the analytical expression for the Bessel
% transform of a Yukawa potential
function vkk = yukawa(lspin,k1,k2,alpha)
%((alpha.^2)+(k1.^2)+(k2.^2))./(2.*k1.*k2)
%tidc = cputime;
z = complex((((alpha.^2)+(k1.^2)+(k2.^2))./(2.*k1.*k2)),0.0);
%complextid = cputime-tidc
% if  lspin==0
%     [qjm1 qj zz] = legendre2(1,z,0.0+0.0i);
%     vkk = (2./pi).*(1.0./(2.*k1.*k2).*real(qjm1));
% else
%tidl = cputime;
[qjm1 qj zz] = legendre2(lspin,z,0.0+0.0i);
%legendretid = cputime-tidl
%qj
vkk = (2./pi).*(1./(2.*k1.*k2)).*real(qj);
%vkk = (1./(2.*k1.*k2)).*real(qj);

%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Yamag0, computing the g0 function in the Yuamaguchi potential
function gg = yamag0(q)
globalvalues
gg = 1./((q.^2)+(Yamaguchibeta.^2));


         