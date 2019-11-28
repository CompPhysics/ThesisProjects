#ifndef WRAPPER_H
#define WRAPPER_H

#include <armadillo>
using namespace arma;

class Wrapper
{
public:
    Wrapper();

    int m_N;
    double m_posMin;
    double m_posMax;
    double m_omega_r;
    int m_nMax;
    int m_nPrimeMax;
    int m_numberOfDimensions;
    double m_h;

    bool m_harmOscPotential;

    vec m_L;
    void setL(vec);

    int m_numberOfEigstates;


};

#endif // WRAPPER_H
