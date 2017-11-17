#include "manyelectrons.h"
#include <cmath>
#include <cassert>
#include "../InitialStates/randomuniform.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include "../Hamiltonians/hamiltonian.h"
#include <iostream>

using namespace std;

bool energyLevelComparison(const pair<int,double> &a,const pair<int,double> &b)
{
       return a.second<b.second;
}

ManyElectrons::ManyElectrons(System* system, double alpha, double beta, double omega, double C, bool Jastrow) :
        WaveFunction(system) {
    assert(omega > 0);
    m_omega = omega;
    assert(alpha >= 0);
    assert(system->getNumberOfDimensions() > 0 && system->getNumberOfDimensions() <= 3);
    m_numberOfDimensions = system->getNumberOfDimensions();
    //m_alpha = alpha;
    //m_alphaOmega = m_alpha*m_omega;
    m_Jastrow = Jastrow;
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_C = C;
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_halfNumberOfParticles = m_numberOfParticles/2.;

    if (m_numberOfParticles == 1) {
        setUpSlaterDetOneParticle();
    }
    else {
        setUpSlaterDet();
    }

    setUpDistances();
    setUpJastrowMat();
}

double ManyElectrons::evaluate(std::vector<class Particle*> particles) {
    // Evaluates the wave function using brute force.

    mat spinUpSlater;
    mat spinDownSlater;
    if (m_numberOfParticles == 1) {
        spinUpSlater = zeros<mat>(m_numberOfParticles, m_numberOfParticles);
        spinDownSlater = zeros<mat>(m_numberOfParticles, m_numberOfParticles);
        spinUpSlater(0,0) = m_SPWFMat(0,0);
    }
    else {
        spinUpSlater = zeros<mat>(m_halfNumberOfParticles, m_halfNumberOfParticles);
        spinDownSlater = zeros<mat>(m_halfNumberOfParticles, m_halfNumberOfParticles);

        for (int i=0; i < m_halfNumberOfParticles; i++) {
            //std::vector<double> rSpinUp = particles[i]->getPosition();//m_system->getParticles()[i]->getPosition();
            //double xSpinUp = rSpinUp[0];
            //double ySpinUp = rSpinUp[1];
            //std::vector<double> rSpinDown = particles[i+m_halfNumberOfParticles]->getPosition();//m_system->getParticles()[i+m_halfNumberOfParticles]->getPosition();
            //double xSpinDown = rSpinDown[0];
            //double ySpinDown = rSpinDown[1];

            for (int j=0; j < m_halfNumberOfParticles; j++) {
                //int nx = m_quantumNumbers(j, 0);
                //int ny = m_quantumNumbers(j, 1);
                spinUpSlater(i,j) = m_SPWFMat(i,j);//evaluateSingleParticleWF(nx, ny, xSpinUp, ySpinUp);
                spinDownSlater(i,j) = m_SPWFMat(i+m_halfNumberOfParticles,j);//evaluateSingleParticleWF(nx, ny, xSpinDown, ySpinDown);
            }
        }
    }

    double beta = m_parameters[1];
    double exponent = 0;
    if (m_Jastrow) {
        for (int i=0; i < m_numberOfParticles; i++) {
            //std::vector<double> r_i = particles[i]->getPosition();

            for (int j=i+1; j < m_numberOfParticles; j++) {
                //std::vector<double> r_j = particles[j]->getPosition();
                //double r_ij = (r_i[0] - r_j[0])*(r_i[0] - r_j[0]) + (r_i[1] - r_j[1])*(r_i[1] - r_j[1]);
//                double r_ij = m_distances(i,j);//sqrt(r_ij);
//                double denom = 1+beta*r_ij;
//                exponent += m_a(i,j)*r_ij/denom;

                exponent += m_JastrowMat(i,j);
            }
        }
    }

    double waveFunction;
    if (m_numberOfParticles == 1) {
        waveFunction = det(spinUpSlater)*exp(exponent);
    }
    else {
        waveFunction = det(spinDownSlater)*det(spinUpSlater)*exp(exponent);
    }

    return waveFunction;
}

std::vector<double> ManyElectrons::computeDerivative(std::vector<class Particle*> particles) {
    //Calculates ∇ψ/ψ for the wave function.

    int i = m_system->getCurrentParticle();
    int numberOfParticles = m_system->getNumberOfParticles();
    std::vector<double> derivative(numberOfParticles*m_numberOfDimensions);

    for (int d = 0; d < m_numberOfDimensions; d++) {
        derivative[i*m_numberOfDimensions+d] = computeSlaterGradient(i)[d];
    }

//    derivative[i*m_numberOfDimensions] = computeSlaterGradient(i)[0]
//                                        ;//+computeJastrowGradient(particles, i)[0];
//    derivative[i*m_numberOfDimensions+1] = computeSlaterGradient(i)[1]
//                                          ;//+computeJastrowGradient(particles, i)[1];

    if (m_Jastrow) {
        for (int d = 0; d < m_numberOfDimensions; d++) {
            derivative[i*m_numberOfDimensions+d] += m_JastrowGrad(i,d);//computeJastrowGradient(particles, i)[d];
            //derivative[i*m_numberOfDimensions+1] += m_JastrowGrad(i,1);//computeJastrowGradient(particles, i)[1];
        }
    }
    return derivative;
    //return 0;
}

std::vector<double> ManyElectrons::computeSlaterGradient(/*std::vector<class Particle*> particles, */int i) {
    // Computes the gradient of the Slater part of the wave function.
    std::vector<double> slaterGradient(m_numberOfDimensions);
    for (int d = 0; d < m_numberOfDimensions; d++) {
        slaterGradient[d] = 0;
    }
//    slaterGradient[0] = 0;
//    slaterGradient[1] = 0;
    //double x = particles[i]->getPosition()[0];
    //double y = particles[i]->getPosition()[1];

    if (i < m_halfNumberOfParticles) {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            //int nx = m_quantumNumbers(j, 0);
            //int ny = m_quantumNumbers(j, 1);
            vec SPWFGradient = m_SPWFDMat(i,j);//computeSPWFDerivative(nx, ny, x, y);
            for (int d = 0; d < m_numberOfDimensions; d++) {
                slaterGradient[d] += SPWFGradient[d]*m_spinUpSlaterInverse(j,i);
            }
//            slaterGradient[0] += SPWFGradient[0]*m_spinUpSlaterInverse(j,i);
//            slaterGradient[1] += SPWFGradient[1]*m_spinUpSlaterInverse(j,i);
        }
    }
    else {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            //int nx = m_quantumNumbers(j, 0);
            //int ny = m_quantumNumbers(j, 1);
            vec SPWFGradient = m_SPWFDMat(i,j);//computeSPWFDerivative(nx, ny, x, y);
            for (int d = 0; d < m_numberOfDimensions; d++) {
                slaterGradient[d] += SPWFGradient[d]*m_spinDownSlaterInverse(j, i-m_halfNumberOfParticles);
            }
//            slaterGradient[0] += SPWFGradient[0]*m_spinDownSlaterInverse(j, i-m_halfNumberOfParticles);
//            slaterGradient[1] += SPWFGradient[1]*m_spinDownSlaterInverse(j, i-m_halfNumberOfParticles);
        }
    }

    return slaterGradient;

}

std::vector<double> ManyElectrons::computeJastrowGradient(std::vector<class Particle*> particles, int k) {
    // Computes the gradient of the Jastrow part of the wave function.
    std::vector<double> jastrowGradient(m_numberOfDimensions);
    for (int d = 0; d < m_numberOfDimensions; d++) {
        jastrowGradient[d] = 0;
    }

//    jastrowGradient[0] = jastrowGradient[1] = 0;

    double beta = m_parameters[1];
    std::vector<double> r_k = particles[k]->getPosition();

    for (int j=0; j < k; j++) {
        std::vector<double> r_j = particles[j]->getPosition();
        double r_kj = m_distances(k,j);//(r_k[0]-r_j[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_j[1])*(r_k[1]-r_j[1]);
        //r_kj = sqrt(r_kj);
        double denom = 1 + beta*r_kj;
        for (int d = 0; d < m_numberOfDimensions; d++) {
            jastrowGradient[d] += (r_k[d]-r_j[d])/r_kj * m_a(k,j)/(denom*denom);
        }
//        jastrowGradient[0] += (r_k[0]-r_j[0])/r_kj * m_a(k, j)/(denom*denom);
//        jastrowGradient[1] += (r_k[1]-r_j[1])/r_kj * m_a(k, j)/(denom*denom);
    }

    for (int j=k+1; j < m_numberOfParticles; j++) {
        std::vector<double> r_j = particles[j]->getPosition();
        double r_kj = m_distances(k,j);//(r_k[0]-r_j[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_j[1])*(r_k[1]-r_j[1]);
        //r_kj = sqrt(r_kj);
        double denom = 1 + beta*r_kj;
        for (int d = 0; d < m_numberOfDimensions; d++) {
            jastrowGradient[d] += (r_k[d]-r_j[d])/r_kj * m_a(k,j)/(denom*denom);
        }
//        jastrowGradient[0] += (r_k[0]-r_j[0])/r_kj * m_a(k, j)/(denom*denom);
//        jastrowGradient[1] += (r_k[1]-r_j[1])/r_kj * m_a(k, j)/(denom*denom);
    }

    return jastrowGradient;

}

double ManyElectrons::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

    //Calculates ∇²ψ/ψ for the wave function.

    double slaterLaplacian = 0;
    double jastrowLaplacian = 0;
    double crossTerm = 0;

    if (m_numberOfParticles == 1) {
        slaterLaplacian += m_SPWFDDMat(0,0)*m_spinUpSlaterInverse(0,0);
    }
    else {
        for (int i=0; i < m_halfNumberOfParticles; i++) {
            //std::vector<double> r_i = particles[i]->getPosition();
            //double x = r_i[0];
            //double y = r_i[1];

            for (int j=0; j < m_halfNumberOfParticles; j++) {
                //int nx = m_quantumNumbers(j, 0);
                //int ny = m_quantumNumbers(j, 1);
                slaterLaplacian += m_SPWFDDMat(i,j)//computeSPWFDoubleDerivative(nx, ny, x, y)
                                  *m_spinUpSlaterInverse(j,i);
            }
        }
        for (int i=m_halfNumberOfParticles; i < m_numberOfParticles; i++) {
            //std::vector<double> r_i = particles[i]->getPosition();
            //double x = r_i[0];
            //double y = r_i[1];

            for (int j=0; j < m_halfNumberOfParticles; j++) {
                //int nx = m_quantumNumbers(j, 0);
                //int ny = m_quantumNumbers(j, 1);
                slaterLaplacian += m_SPWFDDMat(i,j)//computeSPWFDoubleDerivative(nx, ny, x, y)
                                  *m_spinDownSlaterInverse(j,i-m_halfNumberOfParticles);
            }
        }
    }

    if (m_Jastrow) {
        double beta = m_parameters[1];
        double jastrowSum1 = 0;
        double jastrowSum2 = 0;
        int dim = m_numberOfDimensions;

        for (int k=0; k < m_numberOfParticles; k++) {

//            std::vector<double> r_k = particles[k]->getPosition();

            for (int j=0; j < k; j++) {
//                std::vector<double> r_j = particles[j]->getPosition();
                double r_kj = m_distances(k,j);//(r_k[0]-r_j[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_j[1])*(r_k[1]-r_j[1]);
                //r_kj = sqrt(r_kj);
                double denom_kj = (1+beta*r_kj);
                jastrowSum2 += (dim-1)*m_a(k,j)/(r_kj*denom_kj*denom_kj) - 2*m_a(k,j)*beta/(denom_kj*denom_kj*denom_kj);

                for (int i=0; i < k; i++) {

                    for (int d = 0; d < m_numberOfDimensions; d++) {
                        jastrowSum1 += m_dJastrowMat(k,i,d)*m_dJastrowMat(k,j,d);
                    }

//                    std::vector<double> r_i = particles[i]->getPosition();
//                    double r_ki = m_distances(k,i);//(r_k[0]-r_i[0])*(r_k[0]-r_i[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_i[1]);
//                    //r_ki = sqrt(r_ki);
//                    double denom_ki = (1+beta*r_ki);

//                    //double factor1 = (r_k[0]-r_i[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_j[1]);
//                    double factor1 = 0;
//                    for (int d = 0; d < m_numberOfDimensions; d++) {
//                        factor1 += (r_k[d]-r_i[d])*(r_k[d]-r_j[d]);
//                    }

//                    factor1 /= r_ki*r_kj;
//                    double factor2 = m_a(k,i)/(denom_ki*denom_ki) * m_a(k,j)/(denom_kj*denom_kj);
//                    jastrowSum1 += factor1*factor2;
                }
                for (int i=k+1; i < m_numberOfParticles; i++) {

                    for (int d = 0; d < m_numberOfDimensions; d++) {
                        jastrowSum1 += m_dJastrowMat(k,i,d)*m_dJastrowMat(k,j,d);
                    }

//                    std::vector<double> r_i = particles[i]->getPosition();
//                    double r_ki = m_distances(k,i);//(r_k[0]-r_i[0])*(r_k[0]-r_i[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_i[1]);
//                    //r_ki = sqrt(r_ki);
//                    double denom_ki = (1+beta*r_ki);

//                    //double factor1 = (r_k[0]-r_i[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_j[1]);
//                    double factor1 = 0;
//                    for (int d = 0; d < m_numberOfDimensions; d++) {
//                        factor1 += (r_k[d]-r_i[d])*(r_k[d]-r_j[d]);
//                    }

//                    factor1 /= r_ki*r_kj;
//                    double factor2 = m_a(k,i)/(denom_ki*denom_ki) * m_a(k,j)/(denom_kj*denom_kj);
//                    jastrowSum1 += factor1*factor2;
                }
            }
            for (int j=k+1; j < m_numberOfParticles; j++) {
//                std::vector<double> r_j = particles[j]->getPosition();
//                double r_kj = m_distances(k,j);//(r_k[0]-r_j[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_j[1])*(r_k[1]-r_j[1]);
                //r_kj = sqrt(r_kj);
//                double denom_kj = (1+beta*r_kj);
                //jastrowSum2 += (d-1)*m_a(k,j)/(r_kj*denom_kj*denom_kj) - 2*m_a(k,j)*beta/(denom_kj*denom_kj*denom_kj);

                for (int i=0; i < k; i++) {

                    for (int d = 0; d < m_numberOfDimensions; d++) {
                        jastrowSum1 += m_dJastrowMat(k,i,d)*m_dJastrowMat(k,j,d);
                    }

//                    std::vector<double> r_i = particles[i]->getPosition();
//                    double r_ki = m_distances(k,i);//(r_k[0]-r_i[0])*(r_k[0]-r_i[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_i[1]);
//                    //r_ki = sqrt(r_ki);
//                    double denom_ki = (1+beta*r_ki);

//                    //double factor1 = (r_k[0]-r_i[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_j[1]);
//                    double factor1 = 0;
//                    for (int d = 0; d < m_numberOfDimensions; d++) {
//                        factor1 += (r_k[d]-r_i[d])*(r_k[d]-r_j[d]);
//                    }

//                    factor1 /= r_ki*r_kj;
//                    double factor2 = m_a(k,i)/(denom_ki*denom_ki) * m_a(k,j)/(denom_kj*denom_kj);
//                    jastrowSum1 += factor1*factor2;
                }
                for (int i=k+1; i < m_numberOfParticles; i++) {

                    for (int d = 0; d < m_numberOfDimensions; d++) {
                        jastrowSum1 += m_dJastrowMat(k,i,d)*m_dJastrowMat(k,j,d);
                    }

//                    std::vector<double> r_i = particles[i]->getPosition();
//                    double r_ki = m_distances(k,i);//(r_k[0]-r_i[0])*(r_k[0]-r_i[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_i[1]);
//                    //r_ki = sqrt(r_ki);
//                    double denom_ki = (1+beta*r_ki);

//                    //double factor1 = (r_k[0]-r_i[0])*(r_k[0]-r_j[0]) + (r_k[1]-r_i[1])*(r_k[1]-r_j[1]);
//                    double factor1 = 0;
//                    for (int d = 0; d < m_numberOfDimensions; d++) {
//                        factor1 += (r_k[d]-r_i[d])*(r_k[d]-r_j[d]);
//                    }

//                    factor1 /= r_ki*r_kj;
//                    double factor2 = m_a(k,i)/(denom_ki*denom_ki) * m_a(k,j)/(denom_kj*denom_kj);
//                    jastrowSum1 += factor1*factor2;
                }
            }
        }

        jastrowLaplacian = jastrowSum1 + 2*jastrowSum2;

        for (int d=0; d < m_numberOfDimensions; d++) {
            for (int i=0; i < m_numberOfParticles; i++) {
                crossTerm += computeSlaterGradient(i)[d]*m_JastrowGrad(i,d);//computeJastrowGradient(particles, i)[d];
            }
        }
    }

    double laplacian = slaterLaplacian + jastrowLaplacian + 2*crossTerm;
    return laplacian;
    //return 0;
}

std::vector<double> ManyElectrons::computeDerivativeWrtParameters(std::vector<Particle *> particles){
    // Calculates the derivative w.r.t. alpha and beta for the interacting wave function using the analytical expression.

    std::vector<double> derivative(2);
    double slaterUpAlphaDerivative = 0;
    double slaterDownAlphaDerivative = 0;

    for (int i=0; i < m_halfNumberOfParticles; i++) {
        std::vector<double> rSpinUp = particles[i]->getPosition();
//        double xSpinUp = rSpinUp[0];
//        double ySpinUp = rSpinUp[1];
        std::vector<double> rSpinDown = particles[i+m_halfNumberOfParticles]->getPosition();
//        double xSpinDown = rSpinDown[0];
//        double ySpinDown = rSpinDown[1];

        for (int j=0; j < m_halfNumberOfParticles; j++) {
//            int nx = m_quantumNumbers(j, 0);
//            int ny = m_quantumNumbers(j, 1);

            vec n(m_numberOfDimensions);
            for (int d = 0; d < m_numberOfDimensions; d++) {
                n[d] = m_quantumNumbers(j, d);
            }

            slaterUpAlphaDerivative += m_system->getHamiltonian()->computeSPWFAlphaDerivative(n, rSpinUp, j)
                                       *m_spinUpSlaterInverse(j,i);
            slaterDownAlphaDerivative += m_system->getHamiltonian()->computeSPWFAlphaDerivative(n, rSpinDown, j)
                                         *m_spinDownSlaterInverse(j,i);
        }
    }

    double beta = m_parameters[1];
    //double exponent = 0;
    double betaDerivative = 0;
    //double r_ij = 0;

    for (int i=0; i < m_numberOfParticles; i++) {
        //std::vector<double> r_i = particles[i]->getPosition();

        for (int j=i+1; j < m_numberOfParticles; j++) {
            //std::vector<double> r_j = particles[j]->getPosition();

            //for (int k=0; k < m_numberOfDimensions; k++) {
            //    r_ij += (r_i[k]-r_j[k])*(r_i[k]-r_j[k]);
            //}
            //r_ij = sqrt(r_ij);
            double r_ij = m_distances(i,j);

            //exponent += m_a(i,j)*r_ij/(1. + beta*r_ij);
            double denom = (1+beta*r_ij);
            betaDerivative -= m_a(i,j)*r_ij*r_ij/(denom*denom);
        }
    }

    derivative[0] = (slaterUpAlphaDerivative + slaterDownAlphaDerivative)*evaluate(particles);//*exp(exponent);
    derivative[1] = betaDerivative*evaluate(particles);

    return derivative;
    //return 0;
}

double ManyElectrons::computeMetropolisRatio(std::vector<Particle *> particles,
                                            int currentParticle, std::vector<double> positionChange) {
    // Function for calculating the wave function part of the Metropolis ratio
    // both the Slater part and the Jastrow part.

    //std::vector<double> positionOld = particles[currentParticle]->getPosition();
    m_distancesOld = m_distances;
    m_JastrowMatOld = m_JastrowMat;

    for (int i=0; i<m_numberOfDimensions; i++) {
        particles[currentParticle]->adjustPosition(positionChange[i], i);
    }

    //std::vector<double> positionNew = particles[currentParticle]->getPosition();
    m_system->getWaveFunction()->updateDistances(currentParticle);
    m_system->getWaveFunction()->updateSPWFMat(currentParticle);
    m_system->getWaveFunction()->updateJastrow(currentParticle);

    int i = currentParticle;
    double ratioSlaterDet = 0;

    if (i < m_halfNumberOfParticles) {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            //int nx = m_quantumNumbers(j, 0);
            //int ny = m_quantumNumbers(j, 1);
            ratioSlaterDet += m_spinUpSlaterInverse(j,i)
                             *m_SPWFMat(i,j);//evaluateSingleParticleWF(nx, ny, positionNew[0], positionNew[1]);
        }
    }
    else {
        for (int j=0; j < m_halfNumberOfParticles; j++) {
            //int nx = m_quantumNumbers(j, 0);
            //int ny = m_quantumNumbers(j, 1);
            ratioSlaterDet += m_spinDownSlaterInverse(j, i-m_halfNumberOfParticles)
                             *m_SPWFMat(i,j);//evaluateSingleParticleWF(nx, ny, positionNew[0], positionNew[1]);
        }
    }

//    double beta = m_parameters[1];
    double exponent = 0;

    if (m_Jastrow) {
        for (int j=0; j < i; j++) {
            //double r_ijNew = 0;
            //double r_ijOld = 0;
            //std::vector<double> r_j = particles[j]->getPosition();

            //for (int d=0; d < m_numberOfDimensions; d++) {
                //r_ijNew += (positionNew[d] - r_j[d])
                //          *(positionNew[d] - r_j[d]);
                //r_ijOld += (positionOld[d] - r_j[d])
                //          *(positionOld[d] - r_j[d]);
            //}
            //r_ijNew = sqrt(r_ijNew);
            //r_ijOld = sqrt(r_ijOld);
//            double r_ijNew = m_distances(i,j);
//            double r_ijOld = m_distancesOld(i,j);

//            exponent += m_a(i,j)*r_ijNew / (1. + beta*r_ijNew);
//            exponent -= m_a(i,j)*r_ijOld / (1. + beta*r_ijOld);

            exponent += m_JastrowMat(i,j);
            exponent -= m_JastrowMatOld(i,j);
        }
        for (int j=i+1; j < m_numberOfParticles; j++) {
            //double r_ijNew = 0;
            //double r_ijOld = 0;
            //std::vector<double> r_j = particles[j]->getPosition();

            //for (int d=0; d < m_numberOfDimensions; d++) {
                //r_ijNew += (positionNew[d] - r_j[d])
                //          *(positionNew[d] - r_j[d]);
                //r_ijOld += (positionOld[d] - r_j[d])
                //          *(positionOld[d] - r_j[d]);
            //}
            //r_ijNew = sqrt(r_ijNew);
            //r_ijOld = sqrt(r_ijOld);
//            double r_ijNew = m_distances(i,j);
//            double r_ijOld = m_distancesOld(i,j);

//            exponent += m_a(i,j)*r_ijNew / (1. + beta*r_ijNew);
//            exponent -= m_a(i,j)*r_ijOld / (1. + beta*r_ijOld);
            exponent += m_JastrowMat(i,j);
            exponent -= m_JastrowMatOld(i,j);
        }
    }
    double ratioJastrowFactor = exp(exponent);
    m_ratioSlaterDet = ratioSlaterDet;
    m_metropolisRatio = ratioSlaterDet*ratioJastrowFactor;

    return m_metropolisRatio;

}

void ManyElectrons::setUpSlaterDetOneParticle() {

    // Below m_numberOfParticle instead of m_halfNumberOfParticles. Can't have half of one particle.
    m_quantumNumbers = zeros<mat>(m_numberOfParticles, m_numberOfDimensions);

    //----Square well----
    if (m_system->getSquareWellFlag()) {
        for (int d = 0; d < m_numberOfDimensions; d++) {
            m_quantumNumbers(0, d) += 1;
        }
    }
    //----Square well----

    m_spinUpSlater = zeros<mat>(m_numberOfParticles, m_numberOfParticles);
    m_spinDownSlater = zeros<mat>(m_numberOfParticles, m_numberOfParticles);

    m_SPWFMat = zeros<mat>(m_numberOfParticles, m_numberOfParticles);
    m_SPWFDMat = field<vec>(m_numberOfParticles, m_numberOfParticles);
    m_SPWFDDMat = zeros<mat>(m_numberOfParticles, m_numberOfParticles);

    double alpha = m_parameters[0];
    m_system->getHamiltonian()->setAlpha(alpha);

    m_SPWFDMat(0,0) = zeros<vec>(m_numberOfDimensions);

    std::vector<double> r = m_system->getParticles()[0]->getPosition();
    double r2 = 0;

    vec n(m_numberOfDimensions);

    for (int d = 0; d < m_numberOfDimensions; d++) {
        n[d] = m_quantumNumbers(0, d);
        r2 += r[d]*r[d];
    }

    m_a = zeros<mat>(m_numberOfParticles, m_numberOfParticles);
    if (m_numberOfDimensions == 2) { m_a(0,0) = 1./3; }
    else if (m_numberOfDimensions == 3) { m_a(0,0) = 1./4; }


    double expFactor = exp(-alpha*m_omega*(r2)*0.5);

    m_system->getHamiltonian()->setExpFactor(expFactor);

    m_spinUpSlater(0,0) = m_system->getHamiltonian()->evaluateSingleParticleWF(n, r, 0);

    m_SPWFMat(0,0) = m_spinUpSlater(0,0);

    m_SPWFDMat(0,0) = m_system->getHamiltonian()->computeSPWFDerivative(n, r, 0);

    m_SPWFDDMat(0,0) = m_system->getHamiltonian()->computeSPWFDoubleDerivative(n, r, 0);

    m_spinDownSlater(0,0) = m_spinUpSlater(0,0);

    m_spinUpSlaterInverse = m_spinUpSlater.i();
    m_spinDownSlaterInverse = m_spinDownSlater.i();
}

void ManyElectrons::setUpSlaterDet() {
    // Function for setting up the Slater determinant at the begining of the simulation.

    m_quantumNumbers = zeros<mat>(m_halfNumberOfParticles, m_numberOfDimensions);

    if (m_numberOfDimensions == 1) {
        for (int p = 0; p < m_halfNumberOfParticles; p++) {
            m_quantumNumbers(p, 0) = p;
        }
    }

    else if (m_numberOfDimensions == 2) {
        int n = 0;
        int nx = 0;
        int ny = 0;

        for (int p = 0; p < m_halfNumberOfParticles; p++) {
            m_quantumNumbers(p, 0) = nx;    m_quantumNumbers(p, 1) = ny;
            if (ny == n) {
                n++;
                nx = n;
                ny = 0;
            }
            else {
                nx--;
                ny++;
            }
        }
    }

    else {
        int i = 0;
        int nMax = 1;

        for (int n=1; n<m_halfNumberOfParticles; n++) {
            int numberOfEigstates = (n*n*n + 3*n*n +2*n)/6;
            if (numberOfEigstates >= m_halfNumberOfParticles) {
                nMax = n;
                break;
            }
        }

        for (int nx = 0; nx < nMax; nx++) {
            for (int ny = 0; ny < nMax; ny++) {
                for (int nz = 0; nz < nMax; nz++) {
                    if (nx+ny+nz < nMax && i < m_halfNumberOfParticles) {
                        m_quantumNumbers(i,0) = nx;
                        m_quantumNumbers(i,1) = ny;
                        m_quantumNumbers(i,2) = nz;
                        i++;
                    }
                }
            }
        }
        pair<int, int> mapping;
        vector<pair<int, int>> sortingVector;

        for (int i = 0; i < m_halfNumberOfParticles; i++) {
            int nx = m_quantumNumbers(i,0);
            int ny = m_quantumNumbers(i,1);
            int nz = m_quantumNumbers(i,2);

            mapping = make_pair(i, nx+ny+nz);
            sortingVector.push_back(mapping);
        }

        sort(sortingVector.begin(), sortingVector.end(), energyLevelComparison);
        mat quantumNumbersSorted = zeros<mat>(m_halfNumberOfParticles, m_numberOfDimensions);

        for (int i = 0; i < m_halfNumberOfParticles; i++) {
            quantumNumbersSorted(i,0) = m_quantumNumbers(sortingVector[i].first, 0);
            quantumNumbersSorted(i,1) = m_quantumNumbers(sortingVector[i].first, 1);
            quantumNumbersSorted(i,2) = m_quantumNumbers(sortingVector[i].first, 2);
        }
        m_quantumNumbers = quantumNumbersSorted;
    }

    if (m_system->getDoubleWellFlag() && m_numberOfParticles > 2) {
        mat quantumNumbersDoubleWell = zeros<mat>(m_halfNumberOfParticles, m_numberOfDimensions);
        for (int p = 0; p < m_halfNumberOfParticles; p+=2) {
            for (int d = 0; d < m_numberOfDimensions; d++) {
                quantumNumbersDoubleWell(p, d) = m_quantumNumbers(p/2, d);
                if (p+1 < m_halfNumberOfParticles) {
                    quantumNumbersDoubleWell(p+1, d) = m_quantumNumbers(p/2, d);
                }
            }
        }

        m_quantumNumbers = quantumNumbersDoubleWell;
    }

//    //----Square well----
//    if (m_system->getSquareWellFlag()) {
//        for (int p = 0; p < m_halfNumberOfParticles; p++) {
//            for (int d = 0; d < m_numberOfDimensions; d++) {
//                m_quantumNumbers(p, d) += 1;
//            }
//        }
//    }
//    //----Square well----

    m_a = zeros<mat>(m_numberOfParticles, m_numberOfParticles);
    int half = m_halfNumberOfParticles;
    for (int i=0; i < m_numberOfParticles; i++) {
        for (int j=0; j < m_numberOfParticles; j++) {
            if ( ((i < half) && (j < half)) || ((i >= half) && (j >= half)) ) {
                if (m_numberOfDimensions == 2) { m_a(i,j) = 1./3; }
                else if (m_numberOfDimensions == 3) { m_a(i,j) = 1./4; }
            }
            else {
                if (m_numberOfDimensions == 2) { m_a(i,j) = 1.; }
                else if (m_numberOfDimensions == 3) { m_a(i,j) = 1./2; }
            }
        }
    }

    m_spinUpSlater = zeros<mat>(m_halfNumberOfParticles, m_halfNumberOfParticles);
    m_spinDownSlater = zeros<mat>(m_halfNumberOfParticles, m_halfNumberOfParticles);

    m_SPWFMat = zeros<mat>(m_numberOfParticles, m_halfNumberOfParticles);
    m_SPWFDMat = field<vec>(m_numberOfParticles, m_halfNumberOfParticles);
    m_SPWFDDMat = zeros<mat>(m_numberOfParticles, m_halfNumberOfParticles);

    double alpha = m_parameters[0];
    m_system->getHamiltonian()->setAlpha(alpha);

    for (int i=0; i<m_numberOfParticles; i++) {
        for (int j=0; j<m_halfNumberOfParticles; j++) {
            m_SPWFDMat(i,j) = zeros<vec>(m_numberOfDimensions);
        }
    }

    for (int i=0; i < m_halfNumberOfParticles; i++) {
        std::vector<double> rSpinUp = m_system->getParticles()[i]->getPosition();
//        double xSpinUp = rSpinUp[0];
//        double ySpinUp = rSpinUp[1];
        std::vector<double> rSpinDown = m_system->getParticles()[i+m_halfNumberOfParticles]->getPosition();
//        double xSpinDown = rSpinDown[0];
//        double ySpinDown = rSpinDown[1];

        for (int j=0; j < m_halfNumberOfParticles; j++) {
//            nx = m_quantumNumbers(j, 0);
//            ny = m_quantumNumbers(j, 1);
            vec n(m_numberOfDimensions);
            double r2SpinUp = 0;
            double r2SpinDown = 0;
            for (int d = 0; d < m_numberOfDimensions; d++) {
                n[d] = m_quantumNumbers(j, d);
                r2SpinUp += rSpinUp[d]*rSpinUp[d];
                r2SpinDown += rSpinDown[d]*rSpinDown[d];
            }

            double expFactor = exp(-alpha*m_omega*(r2SpinUp)*0.5);
            m_system->getHamiltonian()->setExpFactor(expFactor);

            m_spinUpSlater(i,j) = m_system->getHamiltonian()->evaluateSingleParticleWF(n, rSpinUp, j);

            m_SPWFMat(i,j) = m_spinUpSlater(i,j);

            m_SPWFDMat(i,j) = m_system->getHamiltonian()->computeSPWFDerivative(n, rSpinUp, j);

            m_SPWFDDMat(i,j) = m_system->getHamiltonian()->computeSPWFDoubleDerivative(n, rSpinUp, j);


            expFactor = exp(-alpha*m_omega*(r2SpinDown)*0.5);
            m_system->getHamiltonian()->setExpFactor(expFactor);

            m_spinDownSlater(i,j) = m_system->getHamiltonian()->evaluateSingleParticleWF(n, rSpinDown, j);

            m_SPWFMat(i+m_halfNumberOfParticles, j) = m_spinDownSlater(i,j);

            m_SPWFDMat(i+m_halfNumberOfParticles,j) = m_system->getHamiltonian()->computeSPWFDerivative(n, rSpinDown, j);

            m_SPWFDDMat(i+m_halfNumberOfParticles,j) = m_system->getHamiltonian()->computeSPWFDoubleDerivative(n, rSpinDown, j);

        }
    }

    //cout << m_spinUpSlater << endl;
    m_spinUpSlaterInverse = m_spinUpSlater.i();
    m_spinDownSlaterInverse = m_spinDownSlater.i();
}

void ManyElectrons::setUpDistances() {

    m_distances = zeros<mat>(m_numberOfParticles, m_numberOfParticles);

    for (int i=0; i<m_numberOfParticles; i++) {
        std::vector<double> r_i = m_system->getParticles()[i]->getPosition();

        for (int j=i+1; j<m_numberOfParticles; j++) {
            std::vector<double> r_j = m_system->getParticles()[j]->getPosition();
            double r_ij = 0;

            for (int d = 0; d < m_numberOfDimensions; d++) {
                r_ij += (r_i[d]-r_j[d])*(r_i[d]-r_j[d]);
            }
            m_distances(i,j) = m_distances(j,i) = sqrt(r_ij);
        }
    }
}

void ManyElectrons::setUpJastrowMat() {

    m_JastrowMat = zeros(m_numberOfParticles, m_numberOfParticles);
    m_dJastrowMat = zeros<cube>(m_numberOfParticles, m_numberOfParticles, m_numberOfDimensions);
    m_JastrowGrad = zeros(m_numberOfParticles, m_numberOfDimensions);
    double beta = m_parameters[1];

    for (int i=0; i<m_numberOfParticles; i++) {
        std::vector<double> r_i = m_system->getParticles()[i]->getPosition();
        for (int j=0; j < i; j++) {
            std::vector<double> r_j = m_system->getParticles()[j]->getPosition();
            double r_ij = m_distances(i,j);
            double denom = 1. + beta*r_ij;

            m_JastrowMat(i,j) = m_a(i,j)*r_ij / denom;
            //m_JastrowMat(j,i) = m_JastrowMat(i,j);

            for (int d = 0; d < m_numberOfDimensions; d++) {
                m_JastrowGrad(i,d) += (r_i[d]-r_j[d])/r_ij * m_a(i, j)/(denom*denom);
                //m_JastrowGrad(i,1) += (r_i[1]-r_j[1])/r_ij * m_a(i, j)/(denom*denom);
            }
        }

        for (int j=i+1; j < m_numberOfParticles; j++) {
            std::vector<double> r_j = m_system->getParticles()[j]->getPosition();
            double r_ij = m_distances(i,j);
            double denom = 1. + beta*r_ij;

            m_JastrowMat(i,j) = m_a(i,j)*r_ij / denom;
            //m_JastrowMat(j,i) = m_JastrowMat(i,j);

            for (int d = 0; d < m_numberOfDimensions; d++) {
                m_dJastrowMat(i,j,d) = (r_i[d]-r_j[d])/r_ij * m_a(i, j)/(denom*denom);
                //m_dJastrowMat(i,j,1) = (r_i[1]-r_j[1])/r_ij * m_a(i, j)/(denom*denom);
                m_dJastrowMat(j,i,d) = -m_dJastrowMat(i,j,d);
                //m_dJastrowMat(j,i,1) = -m_dJastrowMat(i,j,1);

                m_JastrowGrad(i,d) += m_dJastrowMat(i,j,d);
                //m_JastrowGrad(i,1) += m_dJastrowMat(i,j,1);
            }
        }
    }
}

void ManyElectrons::updateSlaterDet(int currentParticle) {
    // Function for updating the Slater determinant after every accepted metropolis step.
    int i = currentParticle;
    //std::vector<double> r_i = m_system->getParticles()[i]->getPosition();

    if (i < m_halfNumberOfParticles) {
        mat spinUpSlaterInverseOld = m_spinUpSlaterInverse;

        for (int j=0; j < i; j++) {
            double sum = 0;

            for (int l=0; l <m_halfNumberOfParticles; l++) {
                //int nx = m_quantumNumbers(l, 0);
                //int ny = m_quantumNumbers(l, 1);
                sum += m_SPWFMat(i,l)//evaluateSingleParticleWF(nx, ny, r_i[0], r_i[1])
                      *spinUpSlaterInverseOld(l,j);
            }
            for (int k=0; k < m_halfNumberOfParticles; k++) {
                m_spinUpSlaterInverse(k,j) = spinUpSlaterInverseOld(k,j)
                                            -(sum/m_ratioSlaterDet)*spinUpSlaterInverseOld(k,i);
            }
        }
        for (int j=i+1; j < m_halfNumberOfParticles; j++) {
            double sum = 0;

            for (int l=0; l <m_halfNumberOfParticles; l++) {
                //int nx = m_quantumNumbers(l, 0);
                //int ny = m_quantumNumbers(l, 1);
                sum += m_SPWFMat(i,l)//evaluateSingleParticleWF(nx, ny, r_i[0], r_i[1])
                      *spinUpSlaterInverseOld(l,j);
            }
            for (int k=0; k < m_halfNumberOfParticles; k++) {
                m_spinUpSlaterInverse(k,j) = spinUpSlaterInverseOld(k,j)
                                            -(sum/m_ratioSlaterDet)*spinUpSlaterInverseOld(k,i);
            }
        }
        for (int k=0; k < m_halfNumberOfParticles; k++) {
            m_spinUpSlaterInverse(k,i) = spinUpSlaterInverseOld(k,i)/m_ratioSlaterDet;
        }
    }
    else {
        double iHalf = i-m_halfNumberOfParticles;
        mat spinDownSlaterInverseOld = m_spinDownSlaterInverse;

        for (int j=0; j < iHalf; j++) {
            double sum = 0;

            for (int l=0; l < m_halfNumberOfParticles; l++) {
                //int nx = m_quantumNumbers(l, 0);
                //int ny = m_quantumNumbers(l, 1);
                sum += m_SPWFMat(i,l)//evaluateSingleParticleWF(nx, ny, r_i[0], r_i[1])
                      *spinDownSlaterInverseOld(l,j);
            }
            for (int k=0; k < m_halfNumberOfParticles; k++) {
                m_spinDownSlaterInverse(k,j) = spinDownSlaterInverseOld(k,j)
                                              -(sum/m_ratioSlaterDet)
                                               *spinDownSlaterInverseOld(k, iHalf);
            }
        }
        for (int j=iHalf+1; j < m_halfNumberOfParticles; j++) {
            double sum = 0;

            for (int l=0; l < m_halfNumberOfParticles; l++) {
                //int nx = m_quantumNumbers(l, 0);
                //int ny = m_quantumNumbers(l, 1);
                sum += m_SPWFMat(i,l)//evaluateSingleParticleWF(nx, ny, r_i[0], r_i[1])
                      *spinDownSlaterInverseOld(l,j);
            }
            for (int k=0; k < m_halfNumberOfParticles; k++) {
                m_spinDownSlaterInverse(k,j) = spinDownSlaterInverseOld(k,j)
                                              -(sum/m_ratioSlaterDet)
                                               *spinDownSlaterInverseOld(k, iHalf);
            }
        }
        for (int k=0; k < m_halfNumberOfParticles; k++) {
            m_spinDownSlaterInverse(k, iHalf) = spinDownSlaterInverseOld(k, iHalf)/m_ratioSlaterDet;
        }
    }
}

void ManyElectrons::updateDistances(int currentParticle) {
    // Function for updating the distances between electrons.
    int i = currentParticle;
    std::vector<double> r_i = m_system->getParticles()[i]->getPosition();

    for (int j=0; j<i; j++) {
        std::vector<double> r_j = m_system->getParticles()[j]->getPosition();
        double r_ij = 0;

        for (int d = 0; d < m_numberOfDimensions; d++) {
            r_ij += (r_i[d]-r_j[d])*(r_i[d]-r_j[d]);
        }
        m_distances(i,j) = m_distances(j,i) = sqrt(r_ij);
    }

    for (int j=i+1; j<m_numberOfParticles; j++) {
        std::vector<double> r_j = m_system->getParticles()[j]->getPosition();
        double r_ij = 0;

        for (int d = 0; d < m_numberOfDimensions; d++) {
            r_ij += (r_i[d]-r_j[d])*(r_i[d]-r_j[d]);
        }
        m_distances(i,j) = m_distances(j,i) = sqrt(r_ij);
    }

}

void ManyElectrons::updateSPWFMat(int currentParticle) {

    int i = currentParticle;
    std::vector<double> r_i = m_system->getParticles()[i]->getPosition();
    double alpha = m_parameters[0];

    double r2 = 0;
    for (int d = 0; d < m_numberOfDimensions; d++) {
        r2 += r_i[d]*r_i[d];
    }

    double expFactor = exp(-alpha*m_omega*r2*0.5);
    m_system->getHamiltonian()->setExpFactor(expFactor);

    if (m_numberOfParticles == 1) {
        vec n(m_numberOfDimensions);
        for (int d = 0; d < m_numberOfDimensions; d++) {
            n[d] = m_quantumNumbers(0, d);
        }

        m_SPWFMat(0,0) = m_system->getHamiltonian()->evaluateSingleParticleWF(n, r_i, 0);
        //cout << m_SPWFMat << "    " << r_i[0] << endl;

        m_SPWFDMat(0,0) = m_system->getHamiltonian()->computeSPWFDerivative(n, r_i, 0);

        m_SPWFDDMat(0,0) = m_system->getHamiltonian()->computeSPWFDoubleDerivative(n, r_i, 0);
    }
//    else if (m_system->getDoubleWellFlag()) {
//        for (int j=0; j<m_halfNumberOfParticles; j++) {
//    //        int nx = m_quantumNumbers(j, 0);
//    //        int ny = m_quantumNumbers(j, 1);
//            vec L = m_system->getHamiltonian()->getWellDistance();
//            vec n(m_numberOfDimensions);

//            int denom = 1;
//            for (int d = 0; d < m_numberOfDimensions; d++) {
//                if (L(d) != 0) { denom = 1; }
//            }

//            for (int d = 0; d < m_numberOfDimensions; d++) {
//                n[d] = m_quantumNumbers(j/denom, d);
//            }

//            m_SPWFMat(i,j) = m_system->getHamiltonian()->evaluateSingleParticleWF(n, r_i, j);

//            m_SPWFDMat(i,j) = m_system->getHamiltonian()->computeSPWFDerivative(n, r_i, j);

//            m_SPWFDDMat(i,j) = m_system->getHamiltonian()->computeSPWFDoubleDerivative(n, r_i, j);
//        }
//    }
    else {
        for (int j=0; j<m_halfNumberOfParticles; j++) {
    //        int nx = m_quantumNumbers(j, 0);
    //        int ny = m_quantumNumbers(j, 1);
            vec n(m_numberOfDimensions);
            for (int d = 0; d < m_numberOfDimensions; d++) {
                n[d] = m_quantumNumbers(j, d);
            }

            m_SPWFMat(i,j) = m_system->getHamiltonian()->evaluateSingleParticleWF(n, r_i, j);

            m_SPWFDMat(i,j) = m_system->getHamiltonian()->computeSPWFDerivative(n, r_i, j);

            m_SPWFDDMat(i,j) = m_system->getHamiltonian()->computeSPWFDoubleDerivative(n, r_i, j);
        }
    }

}

void ManyElectrons::updateJastrow(int currentParticle) {

    int p = currentParticle;
    std::vector<double> r_p = m_system->getParticles()[p]->getPosition();
    double beta = m_parameters[1];
    m_dJastrowMatOld = m_dJastrowMat;

    for (int j=0; j<p; j++) {
        std::vector<double> r_j = m_system->getParticles()[j]->getPosition();
        double r_pj = m_distances(p,j);
        double denom = 1. + beta*r_pj;

        m_JastrowMat(p,j) = m_a(p,j)*r_pj / denom;
        m_JastrowMat(j,p) = m_JastrowMat(p,j);

        for (int d = 0; d < m_numberOfDimensions; d++) {
            m_dJastrowMat(p,j,d) = (r_p[d]-r_j[d])/r_pj * m_a(p, j)/(denom*denom);
            //m_dJastrowMat(p,j,1) = (r_p[1]-r_j[1])/r_pj * m_a(p, j)/(denom*denom);
            m_dJastrowMat(j,p,d) = -m_dJastrowMat(p,j,d);
            //m_dJastrowMat(j,p,1) = -m_dJastrowMat(p,j,1);
        }
    }
    for (int j=p+1; j<m_numberOfParticles; j++) {
        std::vector<double> r_j = m_system->getParticles()[j]->getPosition();
        double r_pj = m_distances(p,j);
        double denom = 1. + beta*r_pj;

        m_JastrowMat(p,j) = m_a(p,j)*r_pj / denom;
        m_JastrowMat(j,p) = m_JastrowMat(p,j);

        for (int d = 0; d < m_numberOfDimensions; d++) {
            m_dJastrowMat(p,j,d) = (r_p[d]-r_j[d])/r_pj * m_a(p, j)/(denom*denom);
            //m_dJastrowMat(p,j,1) = (r_p[1]-r_j[1])/r_pj * m_a(p, j)/(denom*denom);
            m_dJastrowMat(j,p,d) = -m_dJastrowMat(p,j,d);
            //m_dJastrowMat(j,p,1) = -m_dJastrowMat(p,j,1);
        }
    }

    m_JastrowGradOld = m_JastrowGrad;

    for (int d = 0; d < m_numberOfDimensions; d++) {
        m_JastrowGrad(p, d) = 0;
        //m_JastrowGrad(p, 1) = 0;

        for (int j=0; j<p; j++) {
            m_JastrowGrad(p, d) += m_dJastrowMat(p,j,d);
            //m_JastrowGrad(p, 1) += m_JastrowMat(p,j,1);
        }
        for (int j=p+1; j<m_numberOfParticles; j++) {
            m_JastrowGrad(p, d) += m_dJastrowMat(p,j,d);
            //m_JastrowGrad(p, 1) += m_JastrowMat(p,j,1);
        }
        for (int i=0; i<p; i++) {
            m_JastrowGrad(i, d) = m_JastrowGradOld(i,d) - m_dJastrowMatOld(i,p,d) + m_dJastrowMat(i,p,d);
            //m_JastrowGrad(i, 1) = m_JastrowGradOld(i,1) - m_JastrowMatOld(i,p,1) + m_JastrowMat(i,p,1);
        }
        for (int i=p+1; i<m_numberOfParticles; i++) {
            m_JastrowGrad(i, d) = m_JastrowGradOld(i,d) - m_dJastrowMatOld(i,p,d) + m_dJastrowMat(i,p,d);
            //m_JastrowGrad(i, 1) = m_JastrowGradOld(i,1) - m_JastrowMatOld(i,p,1) + m_JastrowMat(i,p,1);
        }
    }
}
