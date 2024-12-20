clear all;
close all;


acceptanceDtHe=[
0.0001          100.0;
0.0002          100.0;
0.0003           99.99;
0.0004           99.99;
0.0005           99.99;
0.0006           99.98;
0.0007           99.98;
0.0008           99.98;
0.0009           99.97;
0.001            99.97;
0.002            99.91;
0.003            99.84;
0.004            99.77;
0.005            99.68;
0.006            99.59;
0.007            99.49;
0.008            99.39;
0.009            99.29;
0.01             99.18;
0.02             97.95;
0.03             96.54;
0.04             95.02;
0.05             93.44;
0.06             91.83;
0.07             90.17;
0.08             88.51;
0.09             86.83;
0.1              85.16;
0.2              69.11;
0.3              55.16;
0.4              43.61;
0.5              34.27;
0.6              26.78;
0.7              20.85;
0.8              16.20;
0.9               9.52;
];


dt = acceptanceDtHe(:,1);
E =  acceptanceDtHe(:,2); %%%%%%%%%% ACCEPTANCE



hE= ploterr(dt,E,[],[], 'logx')

 
%%% Ajust line properties (functional)
set(hE(1)                            , ...
  'LineStyle'       , '-'      , ...
  'Marker'          , '.'         , ...
  'Color'           , [0 0 .5]  );

% 
% %%%% Adjust line properties (esthetics)
set(hE(1)                            , ...
  'LineWidth'       , 2           , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 6           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , [.7 .7 .7]  );


%%% Add legend and labels
hXLabel = xlabel('Time step, dt');
hYLabel = ylabel('Accepted moves, %');


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
  'LineWidth'   , 1         );

%%%% limits
%  xlim([9e-5 1]);












%%%%hE=plot(dt,E)

%    plot(dt,E, 'r', 'LineWidth', 2);
%    
%   %%% Ajust line properties (functional)
%  set(hE(1)                            , ...
%    'LineStyle'       , 'none'      , ...
%    'Marker'          , '.'         , ...
%    'Color'           , [.3 .3 .3]  );
%  
%  % 
%  % %%%% Adjust line properties (esthetics)
%  set(hE(1)                            , ...
%    'LineWidth'       , 1           , ...
%    'Marker'          , 'o'         , ...
%    'MarkerSize'      , 6           , ...
%    'MarkerEdgeColor' , [.2 .2 .2]  , ...
%    'MarkerFaceColor' , [.7 .7 .7]  );
%  % 
%  % 
%  % % % % adjust error bar width
%  % hE_c  = ...
%  %     get(hE     , 'Children'    );
%  % errorbarXData          = ...
%  %     get(hE_c(2), 'XData'       );
%  % errorbarXData(4:9:end) = errorbarXData(1:9:end) - hodet;
%  % errorbarXData(7:9:end) = errorbarXData(1:9:end) - hodet;
%  % errorbarXData(5:9:end) = errorbarXData(1:9:end) + hodet;
%  % errorbarXData(8:9:end) = errorbarXData(1:9:end) + hodet;
%  % set(hE_c(2), 'XData', errorbarXData);
%  
%  
%  
%  %%% Add legend and labels
%  hXLabel = xlabel('Time step, dt');
%  hYLabel = ylabel('\% of accepted moves');
%  
%  
%  %%% Add fonts and axis properties
%  %  set( gca                       , ...
%  %     'FontName'   , 'Helvetica' );
%  %  set([hXLabel, hYLabel], ...
%  %     'FontName'   , 'AvantGarde');
%  %  set([hXLabel, hYLabel]  , ...
%  %     'FontSize'   , 14 );
%  %  
%  %  set(gca, ...
%  %   'Box'         , 'on'     , ...
%  %   'TickDir'     , 'in'     , ...
%  %   'TickLength'  , [.02 .02] , ...
%  %   'XMinorTick'  , 'on'      , ...
%  %   'YMinorTick'  , 'on'      , ...
%  %   'YGrid'       , 'on'      , ...
%  %   'XColor'      , [.3 .3 .3], ...
%  %   'YColor'      , [.3 .3 .3], ...
%  %   'XTick'       , xstick, ...
%  %   'YTick'       , ystick);
%  %  %   'LineWidth'   , 1         );
%  %  
%  
%  
