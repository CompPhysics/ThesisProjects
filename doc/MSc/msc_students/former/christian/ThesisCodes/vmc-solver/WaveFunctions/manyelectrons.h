#ifndef PROJECT2_MANYELECTRONS_H
#define PROJECT2_MANYELECTRONS_H
#include "wavefunction.h"
#include <armadillo>
using namespace arma;

class ManyElectrons : public WaveFunction {
public:
    ManyElectrons(class System* system, double alpha, double beta, double omega, double C, bool Jastrow);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    double computeMetropolisRatio(std::vector<Particle *> particles, int currentParticle,
                                  std::vector<double> positionChange);
    std::vector<double> computeDerivative(std::vector<class Particle*> particles);
    std::vector<double> computeDerivativeWrtParameters(std::vector<Particle *> particles);
    std::vector<double> computeSlaterGradient(/*std::vector<Particle *> particles, */int i);
    std::vector<double> computeJastrowGradient(std::vector<Particle *> particles, int i);
    void setUpSlaterDet();
    void setUpSlaterDetOneParticle();
    void setUpDistances();
    void setUpJastrowMat();
    void updateSlaterDet(int currentParticle);
    void updateDistances(int currentParticle);
    void updateSPWFMat(int currentParticle);
    void updateJastrow(int currentParticle);

private:
    int m_numberOfParticles = 0;
    int m_halfNumberOfParticles = 0;
    int m_numberOfDimensions = 0;
    int m_k = 0;
    double m_omega = 0;
    double m_alpha = 0;
    double m_alphaOmega = 0;
    double m_C = 0;
    double m_metropolisRatio = 0;
    double m_ratioSlaterDet = 0;
    mat m_quantumNumbers;
    mat m_spinUpSlater;
    mat m_spinDownSlater;
    mat m_spinUpSlaterInverse;
    mat m_spinDownSlaterInverse;
    mat m_distances;
    mat m_distancesOld;
    mat m_SPWFMat;
    field<vec> m_SPWFDMat;
    mat m_SPWFDDMat;
    double m_cDeterminant;
    mat m_JastrowMat;
    mat m_JastrowMatOld;
    cube m_dJastrowMat;
    cube m_dJastrowMatOld;
    mat m_JastrowGrad;
    mat m_JastrowGradOld;
    mat m_a;
    bool m_Jastrow = false;
};

#endif // PROJECT2_MANYELECTRONS_H
