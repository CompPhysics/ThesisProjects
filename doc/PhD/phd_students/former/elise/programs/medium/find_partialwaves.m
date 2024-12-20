% Determine the number of channels and their spin quantum numbers in
% partial wave expansion [l s j m t tz] 
% NB!! Only 1000 channels possible, must change numbers 
% in fortran code for larger number of channels 

function [spinmat,no_chans] = find_partialwaves
globalvalues

delete spinmatfile.mat
if strcmp(pot,'box') == 1
    spinmat = [0 1 1 0 1 -1];
    no_chans = 1;
elseif (jmin==42 & jmax==42 & itzmin==42 & itzmax==42)
    spinmat = [0 0 0 0 0 0];
    
    %spinmat = [0 0 0 0 1 -1;0 0 0 0 1 0;0 0 0 0 1 1;0 1 1 -1 0 0;0 1 1 0 0 0; 0 1 1 1 0 0]
    %no_chans = 1;
%     spinmat = [0 0 0 0 1 -1;
%                2 0 2 -2 1 -1;
%                2 0 2 -1 1 -1;
%                2 0 2 0 1 -1;
%                2 0 2 1 1 -1;
%                2 0 2 2 1 -1];
%                4 0 4 -4 1 -1;
%                4 0 4 -3 1 -1;
%                4 0 4 -2 1 -1;
%                4 0 4 -1 1 -1;
%                4 0 4 0 1 -1;
%                4 0 4 1 1 -1;
%                4 0 4 2 1 -1;
%                4 0 4 3 1 -1;
%                4 0 4 4 1 -1];
    no_chans = 1;
else 
    [spinmatr,no_chans] = setup_partialwaves(jmin,jmax);
    spinma = spinmatr(:,1:no_chans);
    spinmat = zeros(no_chans,6);
    if symmetric == 1;
        foundctr=0;
        for pwctr=1:no_chans
            spinvec = spinma(:,pwctr);
            l = spinvec(1);
            s = spinvec(2);
            t = spinvec(5);
            if 0.5.*(1-((-1)^(l+s+t)))==1
                foundctr = foundctr+1;
                spinmat(foundctr,:)=spinma(:,pwctr)';
            end
        end
        no_chans = foundctr
        spinmat = spinmat(1:no_chans,:)
        
    else
        spinmat = spinma';
    end
end

partjmin = jmin;
partjmax = jmax;
save spinmatfile partjmax partjmin spinmat no_chans

