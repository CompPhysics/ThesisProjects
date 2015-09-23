clear all;
close all;


%%%%%%%%%%%%%% DT STUDY FOR HE ON TITAN %%%%%%%%%%%%%%%%
%%%% Parameters: Monte Carlo cycles: 	10M
%%%%             Processors: 					4 
%%%%						 minblock:					 	1
%%%%						 maxblock:					 20k 
%%%%						 blocksize:					 500


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix = [
%%%%%%%%
2e-3		-2.8900198			4.0e-4;   
3e-3		-2.8901157			3.5e-4;   
4e-3		-2.8901421			3.0e-4;   
5e-3		-2.8901340			2.6e-4;   
6e-3		-2.8902146			2.4e-4;   
7e-3		-2.8902377			2.2e-4;   
8e-3		-2.8902043			2.0e-4;   
9e-3		-2.8901123			2.0e-4; 
1e-2		-2.8901204			1.8e-4;  
%%%%%%%%%  
2e-2		-2.8901928			1.4e-4;  
3e-2		-2.8903048			1.3e-4;  
4e-2		-2.8902038			1.2e-4;  
5e-2		-2.8902791			1.2e-4;  
6e-2		-2.8900703			1.2e-4;  
7e-2		-2.8901753			1.3e-4;  
8e-2		-2.8903643			1.3e-4;  
9e-2		-2.8901008			1.3e-4
% 1e-1		-2.8900953			1.3e-4;  
%%%%%%%%%
];



dt = matrix(:,1);
E =  matrix(:,2);
err =  matrix(:,3);

hE= ploterr(dt,E,[],err, 'logx')

 
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
  'LineWidth'   , 1         );

%%%% limits
xlim([9.99999999e-4 1e-1]);


