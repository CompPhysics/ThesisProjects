%% NOTAS: dtHe = 0.01; dtBe = 0.001

clear all;
close all;
delete outputFile.m;
delete dataInput/inputParameters.m;

titleSimu = 'tunningDt2DQDot2e';

%%%%%% fixed parameters
sys = '2DQDot2e';
nsd = 2;
nel = 2;
potpar = 1.0;
lambda  = 1.0;
nvp = 2;
enableBlocking  = 'yes';
enableCorrelation = 'yes'
outputFile = 'outputFile.m';

%%%%%% VARIABLES !!!!!!!!!!!!!!!!!!!!!! This should appear en makeFigure

dt_range = [0.02:0.01:0.1];
alpha_range = 0.868; 
beta_range  = 0.354;
nmc_range = 1e7;
nth_range = 0.1*nmc_range;

%%%% Number of simulations
nbSimu = length(dt_range)*length(alpha_range)*length(beta_range)*length(nmc_range)*length(nth_range);

sim=0;
for index_dt = 1:length(dt_range)
  dt = dt_range(index_dt);
  dt_text = num2str(dt);

  dirBlock = sprintf('blocks%s',strrep(dt_text,'.','p'));
  command= sprintf('mkdir %s',dirBlock);
  system(command);

  for index_nmc = 1:length(nmc_range)
    for index_nth = 1:length(nth_range)
      for index_alpha = 1:length(alpha_range)
        for index_beta = 1:length(beta_range)
          %%% parameters to be variated
          nmc = nmc_range(index_nmc);
          nth = nth_range(index_nth);

          alpha = alpha_range(index_alpha);
          beta = beta_range(index_beta);



          %% loop over simulation configurations
          sim=sim+1;
          sprintf('SIMU #%d/%d\n',sim,nbSimu)
          f = fopen('dataInput/inputParameters.m','w');
          fprintf(f,'sys = %s\n',sys);
          fprintf(f,'nsd = %d\n',nsd);
          fprintf(f,'nel = %d\n',nel);
          fprintf(f,'potpar = %0.3g\n',potpar);
          fprintf(f,'lambda = %0.3g\n',lambda);
          fprintf(f,'nmc = %e\n',nmc);
          fprintf(f,'nth = %e\n',nth);
          fprintf(f,'dt = %0.4e\n',dt);
          fprintf(f,'nvp = %d\n',nvp);	
          fprintf(f,'alpha = %0.4g\n', alpha);
          fprintf(f,'beta = %0.4g\n', beta);
          fprintf(f,'enableBlocking = %s\n', enableBlocking);
          fprintf(f,'enableCorrelation = %s\n', enableCorrelation);
          %%%fprintf(f,'output = %s\n', outputFile);
          fclose(f);
    
          pause(0.7);	
  
          %%%%%% PERFORM THE SIMULATION
          system('mpirun -n 4 ./vmc.exe');


		   commandMoveBlock = sprintf('mv blocks_rank*.dat %s',dirBlock);
		   system(commandMoveBlock);
		   commandMoveOutdata = sprintf('mv outputFile.m %s',dirBlock);
		   system(commandMoveOutdata);
		   commandMoveAnalyse = sprintf('cp analyseBlock.exe %s',dirBlock);
		   system(commandMoveAnalyse);
        end
      end
    end
  end
end


%%% BLOCKING
%%for index_dt = 1:length(dt_range)
  %%dt = dt_range(index_dt);
  %%dt_text = num2str(dt);
  %%dirBlock = sprintf('blocks%s',strrep(dt_text,'.','p'));

  %%goIn = sprintf('cd %s',dirBlock);
  %%system(goIn);
		   
  %%commandBlocking = sprintf('./analyseBlock.exe 4 1 20000 500');
  %goOut = sprintf('cd ..');
  %%system(goOut);
%%end
