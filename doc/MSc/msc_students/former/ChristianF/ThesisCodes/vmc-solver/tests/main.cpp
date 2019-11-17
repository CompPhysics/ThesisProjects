#include "UnitTest++/UnitTest++.h"
#include <string.h>
#include <stdio.h>
#include "../system.h"
#include <iostream>
#include <fstream>
#include "../system.h"
#include "../sampler.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include "../WaveFunctions/simplegaussian.h"
#include "../WaveFunctions/repulsivegaussian.h"
#include "../WaveFunctions/twoelectrons.h"
#include "../WaveFunctions/manyelectrons.h"
#include "../WaveFunctions/manyelectrons_coefficients.h"
#include "../Hamiltonians/squarewell.h"
#include "../Hamiltonians/hamiltonian.h"
#include "../Hamiltonians/harmonicoscillator.h"
#include "../Hamiltonians/harmonicoscillatorrepulsive.h"
#include "../Hamiltonians/harmonicoscillatorelectrons.h"
#include "../InitialStates/initialstate.h"
#include "../InitialStates/randomuniform.h"
#include "../VariationMethods/steepestdescent.h"
#include "../Math/random.h"
#include <mpi.h>
#include <cassert>
#include "wrapper.h"

using namespace std;
using namespace UnitTest;



SUITE(VMCsolver) {
    TEST(Sanity) {CHECK_EQUAL(1,1);}

    Wrapper* wrapper = new Wrapper();

    TEST(Initialization) {

        // Initialize MPI parallelization
        //MPI_Comm_size(MPI_COMM_WORLD, &wrapper->m_numprocs);
        //MPI_Comm_rank(MPI_COMM_WORLD, &wrapper->m_my_rank);
        wrapper->m_numprocs = 1;
        wrapper->m_my_rank = 0;
        //wrapper->m_timeStart = MPI_Wtime();
        wrapper->m_numberOfDimensions  = 2;
        vec L(wrapper->m_numberOfDimensions);
        L.fill(0.);
        L(0) = 0.;
        wrapper->setL(L);

        wrapper->m_numberOfParticles   = 2;
        wrapper->m_numberOfSteps       = (int) 1e6;
        wrapper->m_omega               = 1.;
        wrapper->m_alpha               = 0.98456;
        wrapper->m_beta                = 0.40691;
        wrapper->m_gamma               = 2.82843;
        wrapper->m_a                   = 0.0043;
        wrapper->m_stepLength          = 0.5;
        wrapper->m_equilibration       = 0.1;
        wrapper->m_dt                  = 0.01;
        wrapper->m_aElectrons          = 1.;
        wrapper->m_C                   = 1.;
        wrapper->m_analyticalKinetic   = false;
        wrapper->m_importanceSampling  = true;
        wrapper->m_repulsion           = true;
        wrapper->m_quantumDots         = true;
        wrapper->m_twobodyQD           = false;
        wrapper->m_Jastrow             = true;
        wrapper->m_optimizeParameters  = false;
        wrapper->m_saveEnergies        = false;
        wrapper->m_savePositions       = false;
        wrapper->m_showProgress        = true;
        wrapper->m_printToTerminal     = true;
        wrapper->m_useCoeff 		   = false;

        System* system = new System();

        system->setInitialState             (new RandomUniform(system, wrapper->m_numberOfDimensions, wrapper->m_numberOfParticles, wrapper->m_my_rank));

        if (wrapper->m_repulsion && !wrapper->m_quantumDots) {
            system->setHamiltonian              (new HarmonicOscillatorRepulsive(system, wrapper->m_omega, wrapper->m_a, wrapper->m_gamma, wrapper->m_analyticalKinetic));
            system->setWaveFunction             (new RepulsiveGaussian(system, wrapper->m_alpha, wrapper->m_beta, wrapper->m_a));
        }
        if (!wrapper->m_repulsion && !wrapper->m_quantumDots) {
            system->setHamiltonian              (new HarmonicOscillator(system, wrapper->m_omega, wrapper->m_analyticalKinetic));
            system->setWaveFunction             (new SimpleGaussian(system, wrapper->m_alpha));
        }
        if (wrapper->m_quantumDots) {
            if (wrapper->m_twobodyQD) {
                system->setHamiltonian      (new HarmonicOscillatorElectrons(system, wrapper->m_omega, wrapper->m_analyticalKinetic, wrapper->m_repulsion));
                system->setWaveFunction     (new TwoElectrons(system, wrapper->m_alpha, wrapper->m_beta, wrapper->m_omega, wrapper->m_aElectrons, wrapper->m_C, wrapper->m_Jastrow));
            }
            else {
                system->setHamiltonian      (new HarmonicOscillatorElectrons(system, wrapper->m_omega, wrapper->m_analyticalKinetic, wrapper->m_repulsion));
                if (wrapper->m_useCoeff) {
                    system->setWaveFunction (new ManyElectronsCoefficients(system, wrapper->m_alpha, wrapper->m_beta, wrapper->m_omega, wrapper->m_C, wrapper->m_Jastrow));
                }
                else {
                    system->setWaveFunction (new ManyElectrons(system, wrapper->m_alpha, wrapper->m_beta, wrapper->m_omega, wrapper->m_C, wrapper->m_Jastrow));
                }
            }
        }
        system->setEquilibrationFraction    (wrapper->m_equilibration);
        system->setStepLength               (wrapper->m_stepLength);
        system->setTimeStep                 (wrapper->m_dt);
        system->setMyRank                   (wrapper->m_my_rank);

        system->setNumProcs                 (wrapper->m_numprocs);

        if (wrapper->m_optimizeParameters) {
            system->optimizeParameters          (system, wrapper->m_alpha, wrapper->m_beta);
        }
        system->setSaveEnergies             (wrapper->m_saveEnergies);
        system->setSavePositions            (wrapper->m_savePositions);

        system->runMetropolisSteps          (wrapper->m_numMyCycles, wrapper->m_importanceSampling, wrapper->m_showProgress, wrapper->m_printToTerminal);

        //system->MPI_CleanUp                 (wrapper->m_totalE, wrapper->m_totalKE, wrapper->m_totalPE, wrapper->m_totalVariance, wrapper->m_totalAcceptanceRate, wrapper->m_finalMeanDistance,
        //                                     wrapper->m_timeStart, wrapper->m_timeEnd, wrapper->m_totalTime, wrapper->m_numprocs, wrapper->m_numberOfSteps);
        // Merge the files from the nodes into one data file
        //system->mergeOutputFiles            (wrapper->m_numprocs);
    }
}



int main(int nargs, char* args[]) {
    //MPI_Init(&nargs, &args);
    return RunAllTests();
}

