#include "steepestdescent.h"
#include "../system.h"
#include "../sampler.h"
#include "../InitialStates/randomuniform.h"
#include "../WaveFunctions/wavefunction.h"
#include "../Hamiltonians/hamiltonian.h"
#include <cmath>
#include <iostream>
#include "mpi.h"

using namespace std;

SteepestDescent::SteepestDescent(System* system, double stepLengthSD)
{
    m_system = system;
    m_stepLengthSD = stepLengthSD;
}

void SteepestDescent::obtainOptimalParameter(std::vector<double> parameters, double tol, int maxIterations,
                                             int numberOfMetropolisSteps, bool importanceSampling){

    int iteration = 0;  // Count iterations so we can force quit after a given number max iterations
    int numberOfParameters = parameters.size();
    double diff = 0;
    std::vector<double> parametersNew = parameters;
    int my_rank = m_system->getMyRank();

    for (int i=0; i < numberOfParameters; i++) {
        m_derivativeAvg.push_back(0);
    }

    do{

        // Set up an initial state with the updated parameters
        m_system->getInitialState()->setupInitialState();
        for (int i=0; i < numberOfParameters; i++) {
            m_system->getWaveFunction()->adjustParameter(parameters[i], i);
        }
        m_system->getHamiltonian()->setAlpha(parameters[0]);
        //m_system->getHamiltonian()->setUpHermitePolynomials();

        // Run Monte Carlo simulation to find expectation values
        m_system->runMetropolisSteps(numberOfMetropolisSteps, importanceSampling, false, false);

        std::vector<double> derivative(numberOfParameters);  //derivative of local energy.
        // Expectation values needed to calculate derivative of local energy:
        double energy = m_system->getSampler()->getEnergy();
        std::vector<double> waveFuncEnergy(numberOfParameters);
        std::vector<double> waveFuncDerivative(numberOfParameters);

        for (int i=0; i < numberOfParameters; i++) {
            waveFuncEnergy[i] = m_system->getSampler()->getWaveFuncEnergyParameters()[i];
            waveFuncDerivative[i] = m_system->getSampler()->getWaveFuncDerivativeParameters()[i];
            derivative[i] = 2*(waveFuncEnergy[i] - energy*waveFuncDerivative[i]);
        }
        //cout << "WFE: " << waveFuncEnergy[0] << endl;
        //cout << "WFD: " << waveFuncDerivative[0] << endl;
        //cout << "d: " << derivative[0] << endl;

        for (int i=0; i < numberOfParameters; i++) {
            MPI_Reduce(&derivative[i], &m_derivativeAvg[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if (my_rank==0) {
                derivative[i] = m_derivativeAvg[i]/m_system->getNumProcs();
            }
            MPI_Bcast(&derivative[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

        // Find new parameters
        diff = 0;
        for (int i=0; i < numberOfParameters; i++) {
            parametersNew[i] = parameters[i] - derivative[i]*m_stepLengthSD;
            diff += abs(parametersNew[i] - parameters[i]);  // MAKE SURE ABS RETURNS FLOAT!!!
        }
        //m_stepLengthSD *= 0.8;

        //parametersNew = parameters - derivative*m_stepLengthSD;
        //cout << parametersNew[0] << endl;
        //int a = 0;
        //while(a < 1)
        parameters = parametersNew;   // Update parameters
        iteration++;
        std::string upLine = "\e[A";

        if (my_rank == 0) {
            cout << "Iterations: " << iteration << endl;
            for (int i=0; i < numberOfParameters; i++) {
                cout << "Parameter " << i+1 << ": " << parameters[i] << endl;
                upLine += "\e[A";
            }
            cout << upLine;             //"\033[F";
        }


    }while(diff > tol && iteration < maxIterations);
    // Loop ends when requested tolerance for optimal parameters has been reached or after max iterations.

    if (my_rank == 0) {
        cout << "Total iterations: " << iteration << endl;
        if (iteration==maxIterations) {
            cout << "Max iterations reached.\n";
            for (int i=0; i < numberOfParameters; i++) {
                cout << "Parameter" << i+1 << " at max iterations: " << parameters[i] << endl;
            }

        }
        else {
            for (int i=0; i < numberOfParameters; i++) {
                cout << "Optimal " << "Parameter" << i+1 << ": " << parameters[i] << endl;
            }
        }
    }
    // Performing large MC simulation with optimal parameter:
//    m_system->getInitialState()->setupInitialState();
//    for (int i=0; i < numberOfParameters; i++) {
//        m_system->getWaveFunction()->adjustParameter(parameters[i], i);
//    }
//    m_system->runMetropolisSteps((int) 1e6, importanceSampling, true, true);
}
