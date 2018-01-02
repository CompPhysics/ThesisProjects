#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <iostream>

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */

    double rSum = 0;
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double alpha = m_parameters[0];

    // Sum up rSquared for all particles
    for (int i=0; i<numberOfParticles; i++){

        std::vector<double> r_i = particles[i]->getPosition();
        double rSquared = 0;

        for (int j=0; j<numberOfDimensions; j++){
            rSquared += r_i[j]*r_i[j];
        }
        rSum += rSquared;
        //double g = exp(-alpha*rSquared);
        //prod = prod*g;

    }

    // Evaluate the non-interacting wave function
    double waveFunction = exp(-alpha*rSum);

    return waveFunction;
    //return 0;
}

std::vector<double> SimpleGaussian::computeDerivative(std::vector<class Particle*> particles){
    //Calculates ∇ψ/ψ for the non-interacting wave function using the analytical expression.

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double alpha = m_parameters[0];
    std::vector<double> derivative(numberOfParticles*numberOfDimensions);

    for (int i=0; i < numberOfParticles; i++){
        for (int j=0; j < numberOfDimensions; j++){
            derivative[i+j] = -2*alpha*particles[i]->getPosition()[j];
        }
    }
    return derivative;
}

double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */
    //Calculates ∇²ψ/ψ for the non-interacting wave function using the analytical expression.

    double doubleDerivative = 0;
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double alpha = m_parameters[0];

    for (int i=0; i<numberOfParticles; i++){
        std::vector<double> r_i = particles[i]->getPosition();
        double rSquared = 0;
        for (int j=0; j<numberOfDimensions; j++){
            rSquared += r_i[j]*r_i[j];
        }
        doubleDerivative += 2*alpha*(2*alpha*rSquared - numberOfDimensions);
    }

    return doubleDerivative;
    //return 0;
}

std::vector<double> SimpleGaussian::computeDerivativeWrtParameters(std::vector<Particle *> particles) {
    // Calculates the derivative w.r.t. alpha for the non-interacting wave function using the analytical expression.
    // psi_i = exp(-alpha*r_i*r_i) -> d(psi_i)/d(alpha) = -(r_i*r_i)*exp(-alpha*r_i*r_i).

    std::vector<double> derivative(1);
    double rSum = 0;
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double alpha = m_parameters[0];

    for (int i=0; i < numberOfParticles; i++) {
        std::vector<double> r_i = particles[i]->getPosition();
        double rSquared = 0;
        for (int j=0; j < numberOfDimensions; j++) {
            rSquared += r_i[j]*r_i[j];
        }
        rSum += rSquared;
    }

    derivative[0] = -rSum*exp(-alpha*rSum);
    return derivative;
}

double SimpleGaussian::computeMetropolisRatio(std::vector<Particle *> particles,
                                            int currentParticle, std::vector<double> positionChange) {
    int numberOfDimensions = m_system->getNumberOfDimensions();

    // Evaluate the wave function for current positions
    double waveFunctionOld = evaluate(particles);

    // Change position to trial state
    for (int i=0; i<numberOfDimensions; i++){
        particles[currentParticle]->adjustPosition(positionChange[i], i);
    }

    // Evaluate the wave function for the trial state
    double waveFunctionNew = evaluate(particles);

    return waveFunctionNew*waveFunctionNew / (waveFunctionOld*waveFunctionOld);
}

