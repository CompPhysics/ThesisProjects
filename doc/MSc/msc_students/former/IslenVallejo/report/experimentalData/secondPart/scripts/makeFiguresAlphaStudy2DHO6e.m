clear all;
close all;

titleSimu = 'alphaStudy2DHO6e';
savepath pathIslen.m;
chemin=pathIslen;
curentPath=cd;
%  %  %%%%%% fixed parameters
%  %  sys = 'He';
%  %  nsd = 3;
%  %  nel = 2;
%  %  potpar = 2.0;
%  %  nmc = 1e6;
%  %  nth = 1e5;
%  %  dt = 0.01;
%  %  nvp = 2;
%  %  alpha = 1.6785;
%  %  beta = 0.6;
%  %  enableBlocking  = 'yes';
%  %  enable
%  %  outputFile = 'outputFile.m';

%%%%%% VARIABLES !!!!!!!!!!!!!!!!!!!!!! These variables should equal the ones in 'tunning.m'

%%%% 
%%%%%%%% YIELDS ONLY FOR He
dt_range = 0.1;%%[[0.0001:0.0001:0.001] [0.001:0.001:0.01] [0.01: 0.01: 0.1]];
alpha_range = [0.0:0.05:2.0];%2.675%[1.0:0.0625: 2.5]; 
beta_range  = 1.0;
nmc_range = 1e7;
nth_range = 1e5;


%%%%%%%% YIELDS ONLY FOR Be
%  dt_range = 0.001;%%[[0.0001:0.0001:0.001] [0.001:0.001:0.01] [0.01: 0.01: 0.1]];
%  alpha_range = [2.0:0.05:6.0];%2.675%[1.0:0.0625: 2.5]; 
%  beta_range  = 1.0;%[0.01: 0.02: 1];
%  nmc_range = 1e7;
%  nth_range = 0.1*nmc_range;


nbSimu = length(dt_range)*length(alpha_range)*length(beta_range)*length(nmc_range)*length(nth_range);

%%%%% dt, nmc, nth, alpha, beta, accptRatio, energy, variance, error
RESULTS_matrix=zeros(length(dt_range),length(nmc_range),length(nth_range),length(alpha_range),length(beta_range),9);
cd('/mn/kvant/compphys/islen/RESULTS_MATLAB/');
sim=0;
for index_dt = 1:length(dt_range)
  for index_nmc = 1:length(nmc_range)
    for index_nth = 1:length(nth_range)
      for index_alpha = 1:length(alpha_range)
        for index_beta = 1:length(beta_range)
          sim=sim+1;
          nmc = nmc_range(index_nmc);
          nth = nth_range(index_nth);
          dt = dt_range(index_dt);
          alpha = alpha_range(index_alpha);
          beta = beta_range(index_beta);

          %%%%%% READ THE OUTPUT FILE
          
          
          readFile = sprintf('outputFile_%s%d',titleSimu,sim);
          nFiles = sprintf(' file %d over %d\n',sim,numel(RESULTS_matrix));
          getresults(readFile);
          load values.mat
          
          RESULTS_matrix(index_dt,index_nmc,index_nth,index_alpha,index_beta,1) = dt;
          RESULTS_matrix(index_dt,index_nmc,index_nth,index_alpha,index_beta,2) = nmc;
          RESULTS_matrix(index_dt,index_nmc,index_nth,index_alpha,index_beta,3) = nth;
          RESULTS_matrix(index_dt,index_nmc,index_nth,index_alpha,index_beta,4) = alpha;
          RESULTS_matrix(index_dt,index_nmc,index_nth,index_alpha,index_beta,5) = beta;
          RESULTS_matrix(index_dt,index_nmc,index_nth,index_alpha,index_beta,6) = percentAccept;
          RESULTS_matrix(index_dt,index_nmc,index_nth,index_alpha,index_beta,7) = energy;
          RESULTS_matrix(index_dt,index_nmc,index_nth,index_alpha,index_beta,8) = variance;
          RESULTS_matrix(index_dt,index_nmc,index_nth,index_alpha,index_beta,9) = error;
        end
      end
    end
  end
end

%%%% Change to the current path
cd(curentPath)

%  %  system(goBackInitialPath);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                  PLOT                    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%% 2D PLOTTING PART
clear val;
% val(1,:)=RESULTS_matrix(1,1,1,1,:,5);%% x component of val storing all the values of beta for each index of beta (5)
% val(2,:)=RESULTS_matrix(1,1,1,1,:,7);%% y component of val storing all the values of energy (7) for each index of beta (5)

val(1,:)=RESULTS_matrix(1,1,1,:,1,4);%% x component of val storing all the values of alpha for each index of alpha (4)
val(2,:)=RESULTS_matrix(1,1,1,:,1,7);%% y component of val storing all the values of energy (7) for each index of alpha

x=val(1,:);
y=val(2,:);

%%%%Create basic plot
aHeFig = plot(x,y);
fh = figure(1);
hXLabel = xlabel('\alpha', 'FontSize', 16);
hYLabel = ylabel('Energy, \langle E \rangle (au)', 'FontSize', 16);

%%%%% Adjust line properties (Functional)
set(aHeFig,'LineStyle', '-', 'Marker', '.', ...
  'Color', [0 0 .5], ...
  'LineWidth', 2, 'Marker', 'o', ...
  'MarkerSize', 6 , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , [.75 .75 .75]  );


set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde', 'Fontsize', 36);
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 16          );

%%%% Be
% % set(gca, ...
% %   'Box'         , 'on'     , ...
% %   'TickDir'     , 'in'     , ...
% %   'TickLength'  , [.01 .01] , ...
% %   'XMinorTick'  , 'on'      , ...
% %   'YMinorTick'  , 'on'      , ...
% %   'YGrid'       , 'on'      , ...
% %   'XColor'      , [.3 .3 .3], ...
% %   'YColor'      , [.3 .3 .3], ...
% %   'YTick'       , -20.0:0.5:-15.0, ...
% %   'LineWidth'   , 1         );

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , 10.0:0.8:15., ...
  'XTick'       , 0.4:0.2:2.0, ...
  'LineWidth'   , 1         );



axis([0.4 2.0 10. 15.0]); %% 2DHO2e
ylim([9.5 14.5]);

%  %  axis([2.0 6.0 -20.0 -15.0]);    %% Be
%  %  scaletext(1.2);
 
%  %  %%%% 3D PLOTTING PART. Note that now we need 3 variables (4, 5, 7)
%  clear a b c
%  x3D=1:length(alpha_range);
%  y3D=1:length(beta_range);
%  a(x3D,y3D)=RESULTS_matrix(1,1,1,x3D,y3D,4);
%  b(x3D,y3D)=RESULTS_matrix(1,1,1,x3D,y3D,5);
%  c(x3D,y3D)=RESULTS_matrix(1,1,1,x3D,y3D,7);
%  mesh(a,b,c);
%  %  %  [X,Y] = meshgrid(-3:.125:3);
%  %  %  Z = peaks(X,Y);
%  %  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            MINIMUM ENERGY         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  %  %%% MINIMUM ERROR
%  %  %  error=E0(:,:,:,6);
%  %  %  [Y1,I1] =min(error);[Y2,I2]=min(Y1);[Y3,I3]=min(Y2);
%  %  %  minError=error(I1(1,I2(1,1,I3),I3),I2(1,1,I3),I3)
%  %  %  dt_minError=dt_range(I3)
%  %  %  beta_minError=beta_range(I2(1,1,I3))
%  %  %  alpha_minError=alpha_range(I1(1,I2(1,1,I3),I3))
%  %  
%  %  error=E0(:,:,:,6);
%  %  [Y1,I1] =min(error);
%  %  minError=Y1
%  %  alpha_minError=alpha_range(I1)

%%% MINIMUM ENERGY
%  energy=E0(:,:,:,4);
%  [Y1,I1] =min(energy);[Y2,I2]=min(Y1);[Y3,I3]=min(Y2);
%  minEnergy=energy(I1(1,I2(1,1,I3),I3),I2(1,1,I3),I3)
%  dt_minenergy=dt_range(I3)
%  beta_minenergy=beta_range(I2(1,1,I3))
%  alpha_minenergy=alpha_range(I1(1,I2(1,1,I3),I3))

%%%%% dt, nmc, nth, alpha, beta, accptRatio, energy, variance, error
%  energy(1:length(alpha_range),1:length(beta_range))=RESULTS_matrix(1,1,1,:,:,7);
%  [Y1,I1] =min(energy);%% find optimal alpha for each beta
%  [Y2,I2]=min(Y1);%%% find optimal beta with the optimal alpha
%  minEnergy=energy(I1(I2),I2);
%  alpha_minenergy=alpha_range(I1(I2))
%  beta_minenergy=beta_range(I2)

%  
%  for index_alpha = 1:length(alpha_range)
%    for index_beta = 1:length(beta_range)
%      for index_dt = 1:length(dt_range)
%        fprintf(fout,'%0.5g \t %0.5g \t %0.5g \t %0.10g \t %0.10g \t %0.10g \n',E0(index_alpha,index_beta,index_dt,1),E0(index_alpha,index_beta,index_dt,2),E0(index_alpha,index_beta,index_dt,3),E0(index_alpha,index_beta,index_dt,4),E0(index_alpha,index_beta,index_dt,5),E0(index_alpha,index_beta,index_dt,6));	
%      end
%    end
%  end
