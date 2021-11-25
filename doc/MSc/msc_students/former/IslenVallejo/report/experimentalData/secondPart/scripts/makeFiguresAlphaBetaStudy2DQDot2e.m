clear all;
close all;

titleSimu = 'alphaBetaStudy2DQDot2e';
savepath pathIslen.m;
chemin=pathIslen;
curentPath=cd;


dt_range = 0.01;%%[[0.0001:0.0001:0.001] [0.001:0.001:0.01] [0.01: 0.01: 0.1]];
alpha_range = [0.8:0.02:1.2];
beta_range  = [0.08: 0.01: 0.6];
nmc_range = 1e7;
nth_range = 0.1*nmc_range;

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

   %  val(1,:)=RESULTS_matrix(1,1,1,:,1,4);%% x component of val storing all the values of alpha for each index of alpha (4)
   %  val(2,:)=RESULTS_matrix(1,1,1,:,1,7);%% y component of val storing all the values of energy (7) for each index of alpha
   %  
   %  x=val(1,:);
   %  y=val(2,:);
   %  
   %  %%%%Create basic plot
   %  aHeFig = plot(x,y);
   %  fh = figure(1);
   %  hXLabel = xlabel('\alpha', 'FontSize', 16);
   %  hYLabel = ylabel('Energy, \langle E \rangle (au)', 'FontSize', 16);
   %  
   %  %%%%% Adjust line properties (Functional)
   %  set(aHeFig,'LineStyle', '-', 'Marker', '.', ...
   %    'Color', [0 0 .5], ...
   %    'LineWidth', 2, 'Marker', 'o', ...
   %    'MarkerSize', 6 , ...
   %    'MarkerEdgeColor' , [.2 .2 .2]  , ...
   %    'MarkerFaceColor' , [.75 .75 .75]  );
   %  
   %  
   %  set( gca                       , ...
   %      'FontName'   , 'Helvetica' );
   %  set([hXLabel, hYLabel], ...
   %      'FontName'   , 'AvantGarde', 'Fontsize', 36);
   %  set([hXLabel, hYLabel]  , ...
   %      'FontSize'   , 16          );
   %  
   %  %%%% Be
   %  % % set(gca, ...
   %  % %   'Box'         , 'on'     , ...
   %  % %   'TickDir'     , 'in'     , ...
   %  % %   'TickLength'  , [.01 .01] , ...
   %  % %   'XMinorTick'  , 'on'      , ...
   %  % %   'YMinorTick'  , 'on'      , ...
   %  % %   'YGrid'       , 'on'      , ...
   %  % %   'XColor'      , [.3 .3 .3], ...
   %  % %   'YColor'      , [.3 .3 .3], ...
   %  % %   'YTick'       , -20.0:0.5:-15.0, ...
   %  % %   'LineWidth'   , 1         );
   %  
   %  set(gca, ...
   %    'Box'         , 'on'     , ...
   %    'TickDir'     , 'in'     , ...
   %    'TickLength'  , [.01 .01] , ...
   %    'XMinorTick'  , 'on'      , ...
   %    'YMinorTick'  , 'on'      , ...
   %    'YGrid'       , 'on'      , ...
   %    'XColor'      , [.3 .3 .3], ...
   %    'YColor'      , [.3 .3 .3], ...
   %    'YTick'       , 1.0:0.25:5.5, ...
   %    'XTick'       , 0.4:0.2:2.0, ...
   %    'LineWidth'   , 1         );
   %  
   %  
%  %  %  %  
%  %  %  %  axis([0.4 2.0 1.0 5.5]); %% 2DHO2e
%  %  %  %  ylim([1.95 3]);

 
%%%% 3D PLOTTING PART. Note that now we need 3 variables (4, 5, 7)
clear a b c
x3D=1:length(alpha_range);
y3D=1:length(beta_range);
a(x3D,y3D)=RESULTS_matrix(1,1,1,x3D,y3D,4);
b(x3D,y3D)=RESULTS_matrix(1,1,1,x3D,y3D,5);
c(x3D,y3D)=RESULTS_matrix(1,1,1,x3D,y3D,7);
mesh(a,b,c);
xlabel('\alpha');
ylabel('\beta');
zlabel('Energy, \langle E \rangle (au)');
zlim([3.0 3.2]);
set(gca,'XTick',[0.75:0.25:1.25]);


