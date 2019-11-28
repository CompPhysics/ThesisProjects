#include <cassert>
#include <iostream>
#include <cmath>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include "harmonicoscillatorrepulsive.h"

HarmonicOscillatorRepulsive::HarmonicOscillatorRepulsive(System* system, double alpha, double omega,
                                                         double a, double gamma,
                                                         bool analyticalKinetic) :
        Hamiltonian(system, analyticalKinetic, alpha, omega) {
    assert(omega > 0);
    m_omega = omega;
    m_a = a;
    m_gamma = gamma;
}

std::vector<double> HarmonicOscillatorRepulsive::computeLocalEnergy(std::vector<Particle*> particles) {

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double potentialEnergy = 0;

    for (int i=0; i < numberOfParticles; i++){
        double rSquared = 0;
        std::vector<double> r_i = particles[i]->getPosition();
        for (int k=0; k < numberOfDimensions; k++){
            if (k==2) rSquared += m_gamma*m_gamma*r_i[k]*r_i[k];
            else rSquared += r_i[k]*r_i[k];
        }
        potentialEnergy += rSquared;

        for (int j=i+1; j < numberOfParticles; j++){
            double r_ijSquared = 0;
            std::vector<double> r_j = particles[j]->getPosition();
            for (int k=0; k < numberOfDimensions; k++){
                r_ijSquared += (r_i[k] - r_j[k]) * (r_i[k] - r_j[k]);
            }

            if (sqrt(r_ijSquared) <= m_a) potentialEnergy += 1e10;  // Pot. energy is infinite when r_ij<=a
        }
    }
    potentialEnergy *= 0.5*m_omega*m_omega;

    double kineticEnergy = 0;

    if (m_analyticalKinetic == true){
        // Using analytical expression.
        double doubleDerivative = m_system->getWaveFunction()->computeDoubleDerivative(particles);
        kineticEnergy = -0.5*doubleDerivative;
    }
    else{
        // Using numerical diff.
        kineticEnergy = computeKineticEnergy(particles);
    }

    std::vector<double> energies(3);
    energies[0] = kineticEnergy + potentialEnergy;
    energies[1] = kineticEnergy;
    energies[2] = potentialEnergy;

    return energies;
}




