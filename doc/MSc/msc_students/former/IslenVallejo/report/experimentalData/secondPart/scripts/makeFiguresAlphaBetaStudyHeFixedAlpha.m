clear all;
close all;
%%%%opengl software;
titleSimu = 'alphaBetaStudyHe';
savepath pathIslen.m;
chemin=pathIslen;
curentPath=cd;

%%%%%%%% YIELDS ONLY FOR He
dt_range = 0.01;%%[[0.0001:0.0001:0.001] [0.001:0.001:0.01] [0.01: 0.01: 0.1]];
alpha_range = [1.0:0.05:3.0];%2.675%[1.0:0.0625: 2.5]; 
beta_range  = [0.05: 0.05: 1];
nmc_range = 1e7;
nth_range = 0.1*nmc_range;


%%%%%%%% YIELDS ONLY FOR Be
%  dt_range = 0.001;%%[[0.0001:0.0001:0.001] [0.001:0.001:0.01] [0.01: 0.01: 0.1]];
%  alpha_range = [2.0:0.05:6.0];%2.675%[1.0:0.0625: 2.5]; 
%  beta_range  = 1.0;%[0.01: 0.02: 1];
%  nmc_range = 1e7;
%  nth_range = 0.1*nmc_range;


nbSimu = length(dt_range)*length(alpha_range)*length(beta_range)*length(nmc_range)*length(nth_range);

%%%%% dt, nmc, nth, alpha, beta, accptRatio, energy, variance, error
RESULTS_matrix=zeros(length(dt_range),length(nmc_range),length(nth_range),length(alpha_range),length(beta_range),8);

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
%  %  %            cd('/mn/kvant/compphys/islen/RESULTS_MATLAB/');
          
          readFile = sprintf('outputFile_%s%d',titleSimu,sim);
          getresults(readFile);
          load values.mat;
          
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
%  %  clear val;
%  %  %  val(1,:)=RESULTS_matrix(1,1,1,1,:,5);%% x component of val storing all the values of beta for each index of beta (5)
%  %  %  val(2,:)=RESULTS_matrix(1,1,1,1,:,7);%% y component of val storing all the values of energy (7) for each index of beta (5)
%  %  
%  %  val(1,:)=RESULTS_matrix(1,1,1,:,1,4);%% x component of val storing all the values of alpha for each index of alpha (4)
%  %  val(2,:)=RESULTS_matrix(1,1,1,:,1,7);%% y component of val storing all the values of energy (7) for each index of alpha
%  %  
%  %  x=val(1,:);
%  %  y=val(2,:);
%  %  
%  %  % Create basic plot
%  %  aHeFig = plot(x,y);
%  %  
%  %  h = figure(1);
%  %  hXLabel = xlabel('\alpha');
%  %  hYLabel = ylabel('Energy, \langle E \rangle (au)');
%  %  
%  %  %%%%% Adjust line properties (Functional)
%  %  set(aHeFig,'LineStyle', '-', 'Marker', '.', ...
%  %    'Color', [0 0 .5], ...
%  %    'LineWidth', 2, 'Marker', 'o', ...
%  %    'MarkerSize', 6 , ...
%  %    'MarkerEdgeColor' , [.2 .2 .2]  , ...
%  %    'MarkerFaceColor' , [.75 .75 .75]  );
%  %  
%  %  
%  %  set( gca                       , ...
%  %      'FontName'   , 'Helvetica' );
%  %  set([hXLabel, hYLabel], ...
%  %      'FontName'   , 'AvantGarde');
%  %  set([hXLabel, hYLabel]  , ...
%  %      'FontSize'   , 10          );
%  %  
%  %  
%  %  set(gca, ...
%  %    'Box'         , 'on'     , ...
%  %    'TickDir'     , 'in'     , ...
%  %    'TickLength'  , [.01 .01] , ...
%  %    'XMinorTick'  , 'on'      , ...
%  %    'YMinorTick'  , 'on'      , ...
%  %    'YGrid'       , 'on'      , ...
%  %    'XColor'      , [.3 .3 .3], ...
%  %    'YColor'      , [.3 .3 .3], ...
%  %    'YTick'       , -20.0:0.5:-15.0, ...
%  %    'LineWidth'   , 1         );
%  %  
%  %  %  axis([1.0 3.0 -4.0 -3.0]); %% He
%  %  axis([2.0 6.0 -20.0 -15.0]);    %% Be
%  %  scaletext(1.2);
 
%  %  %%%% 3D PLOTTING PART. Note that now we need 3 variables (4, 5, 7)
clear a b c
x3D=1:length(alpha_range);
y3D=1:length(beta_range);
a(x3D,y3D)=RESULTS_matrix(1,1,1,x3D,y3D,4);
b(x3D,y3D)=RESULTS_matrix(1,1,1,x3D,y3D,5);
c(x3D,y3D)=RESULTS_matrix(1,1,1,x3D,y3D,7);
mesh(a,b,c);
xlabel('\alpha');
ylabel('\beta');
zlabel('Energy,  \langle E \rangle (au)');
scaletext(1.2);
%  %  %  [X,Y] = meshgrid(-3:.125:3);
%  %  %  Z = peaks(X,Y);
%  %  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            MINIMUM ENERGY         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% MINIMUM ERROR
%%%%% dt, nmc, nth, alpha, beta, accptRatio, energy, variance, error
error(1:length(alpha_range),1:length(beta_range))=RESULTS_matrix(1,1,1,:,:,9);
[Y1,I1] =min(error);%% find optimal alpha for each beta
[Y2,I2]=min(Y1);%%% find optimal beta with the optimal alpha
minError=error(I1(I2),I2)
alpha_minerror=alpha_range(I1(I2))
beta_minerror=beta_range(I2)

%%% MINIMUM ENERGY
%%%%% dt, nmc, nth, alpha, beta, accptRatio, energy, variance, error
energy(1:length(alpha_range),1:length(beta_range))=RESULTS_matrix(1,1,1,:,:,7);
[Y1,I1] =min(energy);%% find optimal alpha for each beta
[Y2,I2]=min(Y1);%%% find optimal beta with the optimal alpha
minEnergy=energy(I1(I2),I2)
alpha_minenergy=alpha_range(I1(I2))
beta_minenergy=beta_range(I2)



%%%%% Eenrgy matrix with optimal alpha
energyBeta(1:length(beta_range))=RESULTS_matrix(1,1,1,I1(I2),:,7);
errorBeta(1:length(beta_range))=RESULTS_matrix(1,1,1,I1(I2),:,9);
%%%% Plot energy as a function beta for optimal alpha

%%%plot(beta_range,energyBeta)


%%% Create basic plot
plot(beta_range, energyBeta, 'LineWidth', 2, 'LineStyle', '-');
      
hold on;


%%% Plot with errorbars
hE = errorbar(beta_range, energyBeta, errorBeta);
ylim([-2.891 -2.885]);
xlim([0.155 0.6]);

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
    errorbarXData(1:9:end) - 0.005;
errorbarXData(7:9:end) = ....
    errorbarXData(1:9:end) - 0.005;
errorbarXData(5:9:end) = ...
    errorbarXData(1:9:end) + 0.005;
errorbarXData(8:9:end) = ...
    errorbarXData(1:9:end) + 0.005;
set(hE_c(2), 'XData', errorbarXData);



%%% Add legend and labels
hXLabel = xlabel('\beta');
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
  'XTick'       , [0.1:0.08:0.6], ...
  'YTick'       , [-2.89:0.002:-2.88], ...
  'LineWidth'   , 1         );




