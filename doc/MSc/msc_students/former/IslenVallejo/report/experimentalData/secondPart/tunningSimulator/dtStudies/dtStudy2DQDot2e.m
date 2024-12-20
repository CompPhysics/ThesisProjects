matrix = [
%  0.0001   3.0084397   14e-4      100;
%  0.0002   3.0074879   11e-4      100;
%  0.0003   3.0086106   9e-4       100; 
%  0.0004   3.0081546   8e-4       100;
%  0.0005   3.0082982   7e-4       100; 
% 0.0006   3.0088226   6.2e-4     100;
% 0.0007   3.0086718   6e-4       100;
% 0.0008   3.0086234   5.5e-4     100;
% 0.0009   3.0085548   5.2e-4     100; 

0.001     3.0084587   5e-4      100;  
0.002     3.0078266   3.5e-4    99.99;
0.003     3.0079752   3.2e-4    99.98;
0.004     3.0075287   2.5e-4    99.97; 
0.005     3.0075970   2.4e-4    99.96;  
0.006     3.0077140   2.2e-4    99.95;
0.007     3.0076458   2.0e-4    99.94;
0.008     3.0077500   1.9e-4    99.92;
0.009     3.0080983   2.1e-4    99.91


% 0.01     3.0081488    2.0e-4    99.89;
% 0.02     3.0079872    1.6e-4    99.71;
% 0.03     3.0080861    1.4e-6    99.47;
% 0.04     3.0078870    1.2e-4    99.2;
% 0.05     3.0077203    1.08e-4   98.89;
% 0.06     3.0079831    1.14e-4   98.56;
% 0.07     3.0077987    1.7e-4    98.20;
% 0.08     3.0081676    1.5e-4    97.82;
%  %  %  %  %  0.09     3.0081643    97.43
%  %  %  %  %  
%  %  %  %  %  0.1
%  %  %  %  %  0.2
%  %  %  %  %  0.3
%  %  %  %  %  0.4
%  %  %  %  %  0.5
%  %  %  %  %  0.6
%  %  %  %  %  0.7
%  %  %  %  %  0.8
%  %  %  %  %  0.9
];


dt = matrix(:,1);
E =  matrix(:,2);
err =  matrix(:,3);

ploterr(dt, E,[],err)

% hE= ploterr(dt,E,[],err, 'logx')
% 
%  
% %%% Ajust line properties (functional)
% set(hE(1)                            , ...
%   'LineStyle'       , '-'      , ...
%   'Marker'          , '.'         , ...
%   'Color'           , [0 0 .5]  );
% 
% % 
% % %%%% Adjust line properties (esthetics)
% set(hE(1)                            , ...
%   'LineWidth'       , 2           , ...
%   'Marker'          , 'o'         , ...
%   'MarkerSize'      , 6           , ...
%   'MarkerEdgeColor' , [.2 .2 .2]  , ...
%   'MarkerFaceColor' , [.7 .7 .7]  );
% 
% 
% %%% Add legend and labels
% hXLabel = xlabel('Time step, dt');
% hYLabel = ylabel('Energy, \langle E \rangle (au)');
% 
% 
% %%% Add fonts and axis properties
% set( gca                       , ...
%     'FontName'   , 'Helvetica' );
% set([hXLabel, hYLabel], ...
%     'FontName'   , 'AvantGarde');
% set([hXLabel, hYLabel]  , ...
%     'FontSize'   , 14 );
% 
% set(gca, ...
%   'Box'         , 'on'     , ...
%   'TickDir'     , 'in'     , ...
%   'TickLength'  , [.02 .02] , ...
%   'XMinorTick'  , 'on'      , ...
%   'YMinorTick'  , 'on'      , ...
%   'YGrid'       , 'on'      , ...
%   'XColor'      , [.3 .3 .3], ...
%   'YColor'      , [.3 .3 .3], ...
%   'LineWidth'   , 1         );

%%%% limits
%  %  %  %  xlim([9e-5 1]);


