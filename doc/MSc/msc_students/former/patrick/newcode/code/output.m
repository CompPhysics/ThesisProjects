     %%%%%%%%% Input Parameters %%%%%%%%   
d= 2           %%% Dimension of the system
R= 0           %%% Maximum energy level to consider
lambda= 1.00000           %%% Strength of confinement
epsilon= 1e-10           %%% Precision of the HF algorithm


Total_Energy= 3.25331413731549901058% iter #001 


E = nan     % Hartree-Fock energy of the quantum dot


%%%%%%%%   Coefficients of the HF orbitals (one column per particle, one row per single particle state) %%%%%%%%
HForbitals = [1.0000000 0.0000000 
0.0000000 1.0000000 
];


%%%%%%%%   Effective potential (row and column corresponding to single particle states in the order given by singleParticleOrbitals, so sorted by identical angular momentum) %%%%%%%%
U = [1.2533141 0.0000000 
0.0000000 1.2533141 
];


%%%%%%%% Table of states with their quantum numbers, grouped by identical angular momentum %%%%%%%%
NbStates= 2;   % Total number of states
states= sparse(2,3);
states(1,1)=0; % ---state #0 |n=0,ml=0,ms=-1>
states(1,2)=0; % ---state #0 |n=0,ml=0,ms=-1>
states(1,3)=-1; % ---state #0 |n=0,ml=0,ms=-1>
states(2,1)=0; % ---state #1 |n=0,ml=0,ms=1>
states(2,2)=0; % ---state #1 |n=0,ml=0,ms=1>
states(2,3)=1; % ---state #1 |n=0,ml=0,ms=1>


%%%%%%%% Table of couples of states with the same {M,S} %%%%%%%%

StateCouples_0000=sparse(2,1);
%---channel #0 {M=0 S=0}: (0,1) 
StateCouples_0000(1,1)= 0;StateCouples_0000(2,1)= 1;

StateCouples_0001=sparse(2,0);
%---channel #1 {M=0 S=1}: 
