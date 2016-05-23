clear all;
close all;

titleSimu = 'tunningDtHe';
savepath pathIslen.m;
chemin=pathIslen;
curentPath=cd;

%%%%%% fixed parameters
sys = 'He';
nsd = 3;
nel = 2;
potpar = 2.0;
lambda  = 1.0;
nmc = 1e7;
nth = 1e6;
dt = 0.001;
nvp = 2;
alpha = 1.85;
beta = 0.35;
enableBlocking  = 'no';
enableCorrelation = 'yes'
outputFile = 'outputFile.m';

%%%%%% VARIABLES !!!!!!!!!!!!!!!!!!!!!! This should appear en makeFigure
dt_range = [[0.001:0.001:0.01] [0.01:0.01:0.1] [0.1: 0.1: 1.0]];
alpha_range = 1.85; 
beta_range  = 0.35;
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

val(1,:)=RESULTS_matrix(:,1,1,1,1,1);%% x component of val storing all the values of dt for each index of dt (1)
val(2,:)=RESULTS_matrix(:,1,1,1,1,7);%% y component of val storing all the values of energy (7) for each index of dt

x=val(1,:);
y=val(2,:);
error = RESULTS_matrix(:,1,1,1,1,9);

%%%%Create basic plot
%  h1 = plot(x,y);%%%%,x,RESULTS_matrix(:,1,1,1,1,6));
%  
%  figure(1);
%  hXLabel = xlabel('Time step, dt');
%  hYLabel = ylabel('Energy, \langle E \rangle (au)');
%  
%  %%%%% Adjust line properties (Functional)
%  set(h1,'LineStyle', '-', 'Marker', '.', ...
%    'Color', [0 0 .5], ...
%    'LineWidth', 2, 'Marker', 'o', ...
%    'MarkerSize', 6 , ...
%    'MarkerEdgeColor' , [.2 .2 .2]  , ...
%    'MarkerFaceColor' , [.75 .75 .75]  );
%  
%  set(h2,'LineStyle', '-', 'Marker', 'o', ...
%    'Color', [0 0 .5], ...
%    'LineWidth', 2, 'Marker', 'o', ...
%    'MarkerSize', 6 , ...
%    'MarkerEdgeColor' , [.3 .3 .3]  , ...
%    'MarkerFaceColor' , [.7 .7 .7]  );
%  
%  
%  
%  set( gca                       , ...
%      'FontName'   , 'Helvetica' );
%  set([hXLabel, hYLabel], ...
%      'FontName'   , 'AvantGarde', 'FontSize', 16);
%  set([hXLabel, hYLabel]  , ...
%      'FontSize'   , 16          );
%  
%  
%  axis([0.1 1 -2.8915 -2.8865]); %% He
%  hold off
%  %  axis([2.0 6.0 -20.0 -15.0]);    %% Be
%  %  scaletext(1.2);
 


%% Plot with errorbars
hE = errorbar(x, y, error);
%%%ylim([-2.891 -2.885]);
ylim([-2.8904 -2.8896])
xlim([0.001 0.01]);

%%% Ajust line properties (functional)
set(hE                            , ...
  'LineStyle'       , 'none'      , ...
  'Marker'          , '.'         , ...
  'Color'           , [.3 .3 .3]  );


%%%% Adjust line properties (esthetics)
set(hE                            , ...
  'LineWidth'       , 1           , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 6           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , [.7 .7 .7]  );


% % % adjust error bar width
hE_c  = ...
    get(hE     , 'Children'    );
errorbarXData          = ...
    get(hE_c(2), 'XData'       );
errorbarXData(4:9:end) = ...
    errorbarXData(1:9:end) - 0.00005;
errorbarXData(7:9:end) = ....
    errorbarXData(1:9:end) - 0.00005;
errorbarXData(5:9:end) = ...
    errorbarXData(1:9:end) + 0.00005;
errorbarXData(8:9:end) = ...
    errorbarXData(1:9:end) + 0.00005;
set(hE_c(2), 'XData', errorbarXData);



%%% Add legend and labels
hXLabel = xlabel('dt');
hYLabel = ylabel('Energy, \langle E \rangle (au)');


%%% Add fonts and axis properties
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 14 );

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'XTick'       , [0.001:0.001:0.01], ...
  'LineWidth'   , 1         );






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

%  %  set(gca, ...
%  %    'Box'         , 'on'     , ...
%  %    'TickDir'     , 'in'     , ...
%  %    'TickLength'  , [.01 .01] , ...
%  %    'XMinorTick'  , 'on'      , ...
%  %    'YMinorTick'  , 'on'      , ...
%  %    'YGrid'       , 'on'      , ...
%  %    'XColor'      , [.3 .3 .3], ...
%  %    'YColor'      , [.3 .3 .3], ...
%  %    'YTick'       , -4.2:0.2:-3.0, ...
%  %    'XTick'       , 1.0:0.5:3.0, ...
%  %    'LineWidth'   , 1         );

