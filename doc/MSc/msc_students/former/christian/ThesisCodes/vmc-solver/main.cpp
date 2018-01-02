#include <iostream>
#include <unittest++/UnitTest++.h>
#include <fstream>
#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/repulsivegaussian.h"
#include "WaveFunctions/twoelectrons.h"
#include "WaveFunctions/manyelectrons.h"
#include "WaveFunctions/manyelectrons_coefficients.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/harmonicoscillatorrepulsive.h"
#include "Hamiltonians/harmonicoscillatorelectrons.h"
#include "Hamiltonians/doubleharmonicoscillator.h"
#include "Hamiltonians/squarewell.h"
#include "Hamiltonians/finiteharmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "VariationMethods/steepestdescent.h"
#include "Math/random.h"
#include <mpi.h>
#include <cassert>

using namespace std;
using namespace UnitTest;


int main(int nargs, char* args[]) {

    int numprocs, my_rank;
    double totalE, totalKE, totalPE, totalVariance, totalAcceptanceRate, finalMeanDistance;
    double timeStart, timeEnd, totalTime;

    // Initialize MPI parallelization
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    timeStart = MPI_Wtime();

    int numberOfDimensions  = 2;
    int numberOfParticles   = 2;
    int numberOfSteps       = (int) 1e6;              // Monte Carlo cycles
    double omega            = 1.;                     // Oscillator frequency.
    double alpha            = 1.;//0.946207;//0.988423;//0.98456;//0.7;          // Variational parameter.         //3D: 0.983904
    double beta             = 0.227294;//0.362383;//0.38697;//0.40691;//2.82843;      // Variational parameter.         //3D: 0.376667
    double gamma            = 2.82843;
    double a                = 0.0043;                 // Hard core boson diameter.
    double stepLength       = 0.5;                    // Metropolis step length.
    double equilibration    = 0.1;                    // Amount of the total steps used for equilibration.
    double dt               = 0.01;                   // Time step for importance sampling.
    double aElectrons       = 1.; //1./3
    double C                = 1.;                     // Norm constant.
    double distToWall       = 2.;
    double V0               = 1.;

    //Parameter 1 : 0.988423
    //Parameter 2 : 0.38697

    //----N=2----
    //omega 0.10
    //Parameter 1: 1
    //Parameter 2: 0.497735

    //Optimal Parameter1: 1
    //Optimal Parameter2: 0.674777

    //Optimal Parameter1: 1
    //Optimal Parameter2: 0.661385

    //omega 1.00
    //Parameter 1: 1
    //Parameter 2: 0.0459693

    //Parameter 1: 1
    //Parameter 2: 0.2199619

    //Optimal Parameter1: 1
    //Optimal Parameter2: 0.222418

    //Parameter 1: 1
    //Parameter 2: 0.227294

    //----N=6----
    //omega 0.10
    //Parameter 1: 1
    //Parameter 2: 0.349079


    //----N=2----
    //omega 0.10
    //Optimal Parameter1: 1
    //Optimal Parameter2: 0.424292

    //Optimal Parameter1: 1
    //Optimal Parameter2: 0.65428

    //omega 1.00
    //Optimal Parameter1: 1
    //Optimal Parameter2: 0.00874488    (-0.03037927)





    vec L(3);
    L.fill(0.);
    L(0) = 1.;
    //L(1) = 5.;

    bool analyticalKinetic  = true;
    bool importanceSampling = true;
    bool repulsion          = true;                   // Switch for interacting system or not. (Coulomb for manybody qdot)
    bool quantumDots        = true;                   // Switch for quantum dot system.
    bool twobodyQD          = false;                  // Switch for twobody quantum dot system. (no Slater)
    bool Jastrow            = true;                   // Switch for Jastrow factor. (manybody qdot)
    bool optimizeParameters = false;                  // Switch for optimizing variational parameters.

    bool saveEnergies       = true;
    bool savePositions      = true;

    bool showProgress       = true;
    bool printToTerminal    = true;

    bool useCoeff 		    = true;				  // Coefficients c_ij = <ψ_i|φ_j> from the double well potential.

    bool doubleWell         = false;
    bool finiteWell         = false;
    bool squareWell         = true;

    bool runTests           = false;

    int numMyCycles = numberOfSteps/numprocs;

//    cout << "  -- Settings -- " << boolalpha << endl;
//    cout << " Analytical Kinetic : " << analyticalKinetic << endl;
//    cout << " Importance Sampling : " << importanceSampling << endl;
//    cout << " Repulsion : " << repulsion << endl;

    // Initiate System
    System* system = new System();
    // RandomUniform creates a random initial state
    system->setL                        (L);
    system->setDoubleWellFlag           (doubleWell);
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, my_rank));
    // Select which Hamiltonian and trial wave function to use (interacting or non-interacting)
    if (repulsion && !quantumDots) {
    system->setHamiltonian              (new HarmonicOscillatorRepulsive(system, alpha, omega, a, gamma, analyticalKinetic));
    system->setWaveFunction             (new RepulsiveGaussian(system, alpha, beta, a));
    }
    if (!repulsion && !quantumDots) {
    system->setHamiltonian              (new HarmonicOscillator(system, alpha, omega, analyticalKinetic));
    system->setWaveFunction             (new SimpleGaussian(system, alpha));
    }
    if (quantumDots) {
        if (doubleWell) {
            system->setHamiltonian      (new DoubleHarmonicOscillator(system, L, alpha, omega, analyticalKinetic, repulsion));
        }
        else if (finiteWell) {
            system->setHamiltonian      (new FiniteHarmonicOscillator(system, distToWall, alpha, omega, analyticalKinetic, repulsion));
        }
        else if (squareWell) {
            system->setSquareWellFlag   (squareWell);
            system->setHamiltonian      (new SquareWell(system, V0, distToWall, alpha, omega, analyticalKinetic, repulsion));
        }
        else {
            system->setHamiltonian      (new HarmonicOscillatorElectrons(system, alpha, omega, analyticalKinetic, repulsion));
        }
        if (twobodyQD) {
            system->setWaveFunction     (new TwoElectrons(system, alpha, beta, omega, aElectrons, C, Jastrow));
        }
        else {
            if (useCoeff) {
                vec constants;
                system->retrieveConstantsFromFile("../diagonalization/PlotAndData/Constants.dat", constants);
                assert(omega == constants(0));
                system->setWaveFunction (new ManyElectronsCoefficients(system, alpha, beta, omega, C, Jastrow));
            }
            else {
                system->setWaveFunction (new ManyElectrons(system, alpha, beta, omega, C, Jastrow));
            }
        }
    }
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->setTimeStep                 (dt);
    system->setMyRank                   (my_rank);
    system->setNumProcs                 (numprocs);
    // Optimize parameters
    if (optimizeParameters) {
        system->optimizeParameters          (system, alpha, beta);
    }
    system->setSaveEnergies             (saveEnergies);
    system->setSavePositions            (savePositions);
    // Start Monte Carlo simulation
    system->runMetropolisSteps          (numMyCycles, importanceSampling, showProgress, printToTerminal);
    // Compute MPI averages etc.
    system->MPI_CleanUp                 (totalE, totalKE, totalPE, totalVariance, totalAcceptanceRate, finalMeanDistance,
                                         timeStart, timeEnd, totalTime, numprocs, numberOfSteps);
    // Merge the files from the nodes into one data file
    system->mergeOutputFiles            (numprocs);

    if (runTests) { return RunAllTests(); }
    else { return 0; }
}

//12: 65.7
//20: 155.868

/*
 -- System info --
 Number of particles  : 50
 Number of dimensions : 3
 Number of Metropolis steps run : 10^6
 Number of equilibration steps  : 10^5

  -- Wave function parameters --
 Number of parameters : 2
 Parameter 1 : 0.5
 Parameter 2 : 2.82843

  -- Results --
 Energy : 127.306
 Variance : 0.191128
 Acceptance Rate : 0.891037

Computation Time : 13577

Optimizing alpha using steepest descent:

 -- System info --
 Number of particles  : 500
 Number of dimensions : 3
 Number of Metropolis steps run : 10^5
 Number of equilibration steps  : 10^4

  -- Wave function parameters --
 Number of parameters : 1
 Parameter 1 : 0.5

  -- Results --
 Energy : 750
 Variance : 7.64674e-06
 Acceptance Rate : 0.911399

Computation Time : 4580.39


*/


//----2D Reg. HO----
//----N=2----
//omega 0.01
//Parameter 1 : 0.898516
//Parameter 2 : 0.0778596

//omega 0.10
//Parameter 1: 0.957548
//Parameter 2: 0.174197

//omega 0.28
//Parameter 1 : 0.980494
//Parameter 2 : 0.23758

//omega 0.50
//Parameter 1 : 0.987217
//Parameter 2 : 0.288605

//omega 1.00
//Parameter 1 : 0.992372
//Parameter 2 : 0.369956

//----N=6----
//omega 0.1
//Parameter 1 : 0.871378
//Parameter 2 : 0.201104

//omega 0.28
//Parameter 1 : 0.916587
//Parameter 2 : 0.298214

//omega 0.50
//Parameter 1 : 0.930489
//Parameter 2 : 0.38086

//omega 1.00
//Parameter 1 : 0.94881
//Parameter 2 : 0.509732

//----N=12----
//omega 0.10
//Parameter 1 : 0.770629
//Parameter 2 : 0.241339

//omega 0.28
//Parameter 1: 0.855837
//Parameter 2: 0.350573

//omega 0.50
//Parameter 1: 0.894635
//Parameter 2: 0.434976

//omega 1.00
//Parameter 1: 0.930342
//Parameter 2: 0.566965


//----3D Reg. HO----
//----N=2----
//omega 0.01
//Parameter 1 : 0.909867
//Parameter 2 : 0.0476745

//omega 0.10
//Parameter 1 : 0.988223
//Parameter 2 : 0.0947109

//omega 0.28
//Parameter 1 : 0.995213
//Parameter 2 : 0.142882

//omega 0.50
//Parameter 1 : 0.996944
//Parameter 2 : 0.184135

//omega 1.00
//Parameter 1 : 0.998441
//Parameter 2 : 0.249221


//----N=8----
//omega 0.10
//Parameter 1: 0.933177
//Parameter 2: 0.123694

//omega 0.28
//Parameter 1: 0.959758
//Parameter 2: 0.187537

//omega 0.50
//Parameter 1: 0.966478
//Parameter 2: 0.244179

//omega 1.00
//Parameter 1: 0.974329
//Parameter 2: 0.333485


//----2D Double HO----
//----N=2----
//omega 0.01
//Parameter 1 : 0.888028
//Parameter 2 : 0.0731664

//omega 0.10
//Parameter 1 : 0.948708
//Parameter 2 : 0.16236

//omega 0.28
//Parameter 1 : 1.00086
//Parameter 2 : 0.20455

//omega 0.50
//Parameter 1 : 1.01715
//Parameter 2 : 0.245296

//omega 1.00
//Parameter 1 : 1.02299
//Parameter 2 : 0.320903


//----N=4----
//omega 0.10
//Parameter 1: 0.990441
//Parameter 2: 0.158706

//omega 0.28
//Parameter 1: 0.696118
//Parameter 2: 0.488251

//omega 0.50
//Parameter 1: 0.798288
//Parameter 2: 0.520581

//omega 1.00
//Parameter 1: 0.955759
//Parameter 2: 0.441197


//----N=12----
//omega 0.10
//Parameter 1: 1.01397
//Parameter 2: 0.20283


//omega 0.28
//Parameter 1: 0.978235
//Parameter 2: 0.278034

//omega 0.50
//Parameter 1: 0.962793
//Parameter 2: 0.369301

//omega 1.00
//Parameter 1: 0.988232
//Parameter 2: 0.443659


//----3D Double HO----
//----N=2----
//omega 0.01
//Parameter 1 : 0.975627
//Parameter 2 : 0.042176

//omega 0.10
//Parameter 1 : 0.981495
//Parameter 2 : 0.0910774

//omega 0.28
//Parameter 1 : 1.00215
//Parameter 2 : 0.124679

//omega 0.50
//Parameter 1 : 1.01017
//Parameter 2 : 0.159175

//omega 1.00
//Parameter 1 : 0.998497
//Parameter 2 : 0.265424

//----N=4----
//omega 0.10
//Parameter 1: 1.022585
//Parameter 2: 0.0819739

//omega 0.28
//Parameter 1: 0.910207
//Parameter 2: 0.202278

//omega 0.50
//Parameter 1: 0.944187
//Parameter 2: 0.229013

//omega 1.00
//Parameter 1: 0.946207
//Parameter 2: 0.362383
