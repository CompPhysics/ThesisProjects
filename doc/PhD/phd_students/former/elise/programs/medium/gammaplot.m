% 

function gammaplot(Gamma,axis,fignum,type)
globalvalues

if strcmp(type,'k')==1
Gama = squeeze(Gamma)';
if size(Gama,1)~=1
    disp('Gamma has too many dimensions for plot')
else
plotlim = round(n_rel_points.*(4./5));
 figure(fignum) 
%  plot(axis(1:plotlim)./hbarc,real(Gama(1:plotlim)).*(hbarc.^3),'b', ...
%      axis(1:plotlim)./hbarc,imag(Gama(1:plotlim)).*(hbarc.^3),'g')
 plot(axis(1:plotlim)./hbarc,real(Gama(1:plotlim)),'b', ...
     axis(1:plotlim)./hbarc,imag(Gama(1:plotlim)),'g')
end
elseif strcmp(type,'omega')==1
    Gama = squeeze(Gamma(:,:,1));
if size(Gama,1)~=1
    disp('Gamma has too many dimensions for plot')
else
%plotstart = round(n_omega_points.*(1./4));
plotstart = 1;
plotend = n_omega_points; %-plotstart;
 figure(fignum) 
 plot(axis(plotstart:plotend),real(Gama(plotstart:plotend)).*(hbarc.^3),'b', ...
     axis(plotstart:plotend),imag(Gama(plotstart:plotend)).*(hbarc.^3),'g')
end
else
    disp('Unknown plot type')
end