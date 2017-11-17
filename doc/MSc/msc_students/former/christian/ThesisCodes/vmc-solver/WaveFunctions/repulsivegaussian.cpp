#include "repulsivegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <iostream>

RepulsiveGaussian::RepulsiveGaussian(System* system, double alpha, double beta, double a) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_a = a;
}

double RepulsiveGaussian::evaluate(std::vector<class Particle*> particles) {

    double rSum = 0;
    double prod = 1;
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double alpha = m_parameters[0];
    double beta = m_parameters[1];

    // Evaluate the wave function
    for (int i=0; i<numberOfParticles; i++){

        std::vector<double> r_i = particles[i]->getPosition();
        double rSquared = 0;

        // uncorrelated term
        for (int k=0; k<numberOfDimensions; k++){
            if (k==2) rSquared += beta*r_i[k]*r_i[k];
            else rSquared += r_i[k]*r_i[k];
        }
        rSum += rSquared;
        //double g = exp(-alpha * rSquared);
        //prod = prod*g;

        // correlation term
        for (int j=i+1; j<numberOfParticles; j++){

            std::vector<double> r_j = particles[j]->getPosition();
            double r_ij = 0;

            for (int k=0; k<numberOfDimensions; k++){
                r_ij += (r_i[k] - r_j[k]) * (r_i[k] - r_j[k]);
            }
            r_ij = sqrt(r_ij);

            double f = 0;
            if (r_ij > m_a) f = 1 - m_a/r_ij;
            prod *= f;
        }
    }

    double waveFunction = exp(-alpha*rSum)*prod;

    return waveFunction;
    //return 0;
}

std::vector<double> RepulsiveGaussian::computeDerivative(std::vector<class Particle*> particles){
    //Calculates ∇ψ/ψ for the interacting wave function using the analytical expression.

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double alpha = m_parameters[0];
    double beta = m_parameters[1];
    std::vector<double> derivative(numberOfParticles*numberOfDimensions);

    // First term
    for (int i=0; i < numberOfParticles; i++){
        std::vector<double> r_i = particles[i]->getPosition();
        for (int k=0; k < numberOfDimensions; k++){
            if (k==2) derivative[i+k] = -2*alpha*beta*r_i[k];
            else derivative[i+k] = -2*alpha*r_i[k];
        }
        // Second term
        for (int j=i+1; j < numberOfParticles; j++){
            std::vector<double> r_j = particles[j]->getPosition();
            double r_ij = 0;
            //double abs_r_i = 0;
            for (int k=0; k < numberOfDimensions; k++){
                r_ij += (r_i[k] - r_j[k]) * (r_i[k] - r_j[k]);
                //abs_r_i += r_i[k]*r_i[k];
            }
            double r_ij2 = r_ij;
            r_ij = sqrt(r_ij);
            //abs_r_i = sqrt(abs_r_i);
            for (int k=0; k < numberOfDimensions; k++){
                if (r_ij > m_a){
                    derivative[i+k] += m_a*(r_i[k]-r_j[k]) / (r_ij2*(r_ij - m_a));
                    //derivative[i+k] += (r_i[k]-r_j[k])/r_ij*(1./(1 - m_a/r_ij))*m_a/(r_ij2);
                    //derivative[i+k] += (1./(1 - m_a/r_ij))*m_a*r_i[k]/(r_ij*r_ij*abs_r_i);
                }
            }
        }
    }
    return derivative;
}

double RepulsiveGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

    //Calculates ∇²ψ/ψ for the interacting wave function using the analytical expression.
    double doubleDerivative = 0;
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double alpha = m_parameters[0];
    double beta = m_parameters[1];

    for (int i=0; i<numberOfParticles; i++){

        //First term
        std::vector<double> r_i = particles[i]->getPosition();
        //double rSquared = 0;
        for (int k=0; k<numberOfDimensions; k++){
            //if (k==2) rSquared += beta*r_i[k]*r_i[k];
            //else rSquared += r_i[k]*r_i[k];
            if (k==2) doubleDerivative += 2*alpha*beta*(2*alpha*beta*r_i[k]*r_i[k] - 1);
            else doubleDerivative += 2*alpha*(2*alpha*r_i[k]*r_i[k] - 1);
        }

        double corr1sum = 0;
        double corr2sum = 0;

        for (int j=i+1; j < numberOfParticles; j++){
            //if (i != j){
            std::vector<double> r_j = particles[j]->getPosition();
            double r_ij = 0;
            //double abs_r_i = 0;
            for (int k=0; k < numberOfDimensions; k++){
                r_ij += (r_i[k] - r_j[k]) * (r_i[k] - r_j[k]);
                //abs_r_i += r_i[k]*r_i[k];
            }
            double r_ij2 = r_ij;
            r_ij = sqrt(r_ij);
            //abs_r_i = sqrt(abs_r_i);

            for (int k=0; k < numberOfDimensions; k++){
                if (r_ij > m_a){
                    double corr1 = m_a*(r_i[k]-r_j[k])/( r_ij2*(r_ij - m_a) );
                    double corr2 = m_a*(m_a-2*r_ij)/(r_ij*r_ij*(r_ij-m_a)*(r_ij-m_a))
                                    + (numberOfDimensions-1.)/r_ij*m_a/(r_ij*(r_ij-m_a));//(1./(1 - m_a/r_ij))*m_a/(r_ij2);

//                    if (k==2) doubleDerivative += 2*(-2*alpha*beta*r_i[k])*corr1;
//                    else doubleDerivative += 2*(-2*alpha*r_i[k])*corr1;

//                    doubleDerivative += corr1*corr1;
//                    doubleDerivative += corr2;

                    corr1sum += corr1;
                    corr2sum += corr2;
                }
            }
        }
        // Second term
        for (int k=0; k < numberOfDimensions; k++){
            if (k==2) doubleDerivative += 2*(-2*alpha*beta*r_i[k])*corr1sum;
            else doubleDerivative += 2*(-2*alpha*r_i[k])*corr1sum;
        }
        // Third and fourth term
        doubleDerivative += corr1sum*corr1sum;
        doubleDerivative += corr2sum;
    }

    return doubleDerivative;
    //return 0;    
}

std::vector<double> RepulsiveGaussian::computeDerivativeWrtParameters(std::vector<Particle *> particles){
    // Calculates the derivative w.r.t. alpha for the interacting wave function using the analytical expression.
    // Correlation factor has no alpha dependency so from the product rule the analytical expression is simply
    // the analytical expression for the non-interacting case multiplied by the correlation factor.

    std::vector<double> derivative(1);
    double rSum = 0;
    double prod = 1;
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double alpha = m_parameters[0];
    double beta = m_parameters[1];

    for (int i=0; i < numberOfParticles; i++) {
        std::vector<double> r_i = particles[i]->getPosition();
        double rSquared = 0;
        for (int k=0; k<numberOfDimensions; k++){
            if (k==2) rSquared += beta*r_i[k]*r_i[k];
            else rSquared += r_i[k]*r_i[k];
        }
        //double g = exp(-alpha * rSquared);
        //prod = prod*g;
        rSum += rSquared;

        for (int j=i+1; j<numberOfParticles; j++){

            std::vector<double> r_j = particles[j]->getPosition();
            double r_ij = 0;

            for (int k=0; k<numberOfDimensions; k++){
                r_ij += (r_i[k] - r_j[k]) * (r_i[k] - r_j[k]);
            }
            r_ij = sqrt(r_ij);

            double f = 0;
            if (r_ij > m_a) f = 1 - m_a/r_ij;
            prod *= f;
        }
    }

    derivative[0] = -rSum*exp(-alpha*rSum)*prod;
    return derivative;
}

double RepulsiveGaussian::computeMetropolisRatio(std::vector<Particle *> particles,
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




