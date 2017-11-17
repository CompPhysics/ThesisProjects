#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using namespace std;

Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::setEnergy(double energy) {
    m_energy = energy;
}

void Sampler::setKineticEnergy(double kineticEnergy) {
    m_kineticEnergy = kineticEnergy;
}

void Sampler::setPotentialEnergy(double potentialEnergy) {
    m_potentialEnergy = potentialEnergy;
}

void Sampler::setVariance(double variance) {
    m_variance = variance;
}

void Sampler::setAcceptanceRate(double acceptanceRate) {
    m_acceptanceRate = acceptanceRate;
}

void Sampler::setMeanDistance(double meanDistance) {
    m_meanDistance = meanDistance;
}

void Sampler::sample(bool acceptedStep) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_stepNumber == 0) {
        m_cumulativeEnergy =             0;
        m_cumulativeSquaredEnergy =      0;
        m_cumulativeAcceptedSteps =      0;
        for (int i=0; i < m_system->getWaveFunction()->getNumberOfParameters(); i++) {
            m_cumulativeWFuncDerivativeParameters.push_back(0);
            m_cumulativeWFuncEnergyParameters.push_back(0);
        }
    }

    // Sample local energy
    std::vector<double> energies = m_system->getHamiltonian()->
                                 computeLocalEnergy(m_system->getParticles());
    double localEnergy = energies[0];
    double kineticEnergy = energies[1];
    double potentialEnergy = energies[2];

    m_cumulativeEnergy  += localEnergy;
    m_cumulativeSquaredEnergy += localEnergy*localEnergy;
    m_cumulativeKineticEnergy += kineticEnergy;
    m_cumulativePotentialEnergy += potentialEnergy;

    std::vector<Particle*> particles = m_system->getParticles();

    // Sample mean distance
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    if (numberOfParticles==2) {
        int counter = 0;
        double avgDistance = 0;
        for (int i=0 ; i < numberOfParticles; i++) {
            double r_ij = 0;
            std::vector<double> r_i = particles[i]->getPosition();
            for (int j=i+1; j < numberOfParticles; j++) {
                std::vector<double> r_j = particles[j]->getPosition();
                for (int k=0; k < numberOfDimensions; k++) {
                    r_ij += (r_i[k]-r_j[k])*(r_i[k]-r_j[k]);
                }
                r_ij = sqrt(r_ij);
                avgDistance += r_ij;
                counter ++;
            }
        }
        avgDistance /= counter;
        m_cumulativeDistance += avgDistance;
//        cout << particles[0]->getPosition()[0] << "   " << particles[0]->getPosition()[1] << endl;
//        cout << particles[1]->getPosition()[0] << "   " << particles[1]->getPosition()[1] << endl;
    }

    // Sample things needed for the steepest descent method:
    double waveFunction = m_system->getWaveFunction()->evaluate(particles);
    std::vector<double> waveFuncDerivativeParameters = m_system->getWaveFunction()->computeDerivativeWrtParameters(particles);
    int numberOfParameters = m_system->getWaveFunction()->getNumberOfParameters();
    for (int i=0; i < numberOfParameters; i++) {
        waveFuncDerivativeParameters[i] /= waveFunction;
        m_cumulativeWFuncDerivativeParameters[i] += waveFuncDerivativeParameters[i];
        m_cumulativeWFuncEnergyParameters[i] += waveFuncDerivativeParameters[i]*localEnergy;
    }

    // Sample whether the step is accepted or not in order to find total acceptance ratio
    m_cumulativeAcceptedSteps += acceptedStep;

    saveToFile(localEnergy);
    m_stepNumber++;
}

void Sampler::printOutputToTerminal() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    double  ct = m_system->getComputationTime();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << "    " << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << " Variance : " << m_variance << endl;
    cout << " Kinetic Energy : " << m_kineticEnergy << endl;
    cout << " Potential Energy : " << m_potentialEnergy << endl;
    cout << " Acceptance Rate : " << m_acceptanceRate << endl;
    if (m_meanDistance != 0) cout << " Mean Distance : " << m_meanDistance << endl;
    cout << endl;
    cout << "Computation Time : " << ct << endl;
    cout << endl;
}

void Sampler::computeAverages() {
    // Compute the averages of the sampled quantities.

    // We only sample after the equilibration steps so we need to average by the number of steps after equilibration.
    m_energy = m_cumulativeEnergy / (double) m_stepNumber;
    m_squaredEnergy = m_cumulativeSquaredEnergy / (double) m_stepNumber;
    m_variance = m_squaredEnergy - m_energy*m_energy;

    m_kineticEnergy = m_cumulativeKineticEnergy / (double) m_stepNumber;
    m_potentialEnergy = m_cumulativePotentialEnergy / (double) m_stepNumber;

    int numberOfParticles = m_system->getWaveFunction()->getNumberOfParameters();
    for (int i=0; i < numberOfParticles; i++) {
        m_waveFuncDerivativeParameters.push_back(m_cumulativeWFuncDerivativeParameters[i] / (double) m_stepNumber);
        m_waveFuncEnergyParameters.push_back(m_cumulativeWFuncEnergyParameters[i] / (double) m_stepNumber);
    }

    m_acceptanceRate = m_cumulativeAcceptedSteps / (double) m_stepNumber;

    m_meanDistance = m_cumulativeDistance / (double) m_stepNumber;
}

void Sampler::saveToFile(double localEnergy){
    if (m_system->getSaveEnergies()){
        //if (m_stepNumber == 0 && m_system->getMyRank() == 0) system("rm energies.dat");

        //FILE *outfile = fopen("energies.dat", "ab");
        //char filename[100];
        //sprintf(filename, "energies%i.dat", m_system->getMyRank());
        //FILE *outfile = fopen(filename, "ab");
        fprintf(m_system->getEnergiesFile(), "%f\n", localEnergy);
        //fclose(outfile);
    }

    if (m_system->getSavePositions()){
        //if (m_stepNumber == 0 && m_system->getMyRank() == 0) system("rm positions.dat");

        FILE *outfile = m_system->getPositionsFile();
        std::vector<Particle*> particles = m_system->getParticles();
        int numberOfParticles = m_system->getNumberOfParticles();
        int numberOfDimensions = m_system->getNumberOfDimensions();

        for (int i=0; i < numberOfParticles; i++){
            std::vector<double> r_i = particles[i]->getPosition();
            for (int j=0; j < numberOfDimensions; j++){
                fprintf(outfile, "%f ", r_i[j]);
            }
            fprintf(outfile, "\n");
        }
        //fclose(outfile);
    }
}
