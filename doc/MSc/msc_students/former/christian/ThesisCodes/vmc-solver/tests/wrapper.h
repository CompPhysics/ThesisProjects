#ifndef WRAPPER_H 
#define WRAPPER_H 
#include <armadillo> 
using namespace arma; 

class Wrapper {
public:
    Wrapper();

    vec m_L;
    void setL(vec);

    int m_numprocs, m_my_rank;
    double m_totalE, m_totalKE, m_totalPE, m_totalVariance, m_totalAcceptanceRate, m_finalMeanDistance;
    double m_timeStart, m_timeEnd, m_totalTime;
    int m_numMyCycles = m_numberOfSteps/m_numprocs;

    int m_numberOfDimensions;
    int m_numberOfParticles;
    int m_numberOfSteps;
    double m_omega;
    double m_alpha;
    double m_beta;
    double m_gamma;
    double m_a;
    double m_stepLength;
    double m_equilibration;
    double m_dt;
    double m_aElectrons;
    double m_C;
    bool m_analyticalKinetic;
    bool m_importanceSampling;
    bool m_repulsion;
    bool m_quantumDots;
    bool m_twobodyQD;
    bool m_Jastrow;
    bool m_optimizeParameters;
    bool m_saveEnergies;
    bool m_savePositions;
    bool m_showProgress;
    bool m_printToTerminal;
    bool m_useCoeff; 

};

#endif // WRAPPER_H
