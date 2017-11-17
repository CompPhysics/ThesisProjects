#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
using namespace arma;


class WaveFunction
{
public:
    WaveFunction(class System* system, double omega);
    virtual vec harmonicOscillatorBasis(mat r, int n) = 0;
    virtual vec potential (vec r, double L) = 0;
    vec computeHermitePolynomial(int nValue, vec position);

protected:
    double          m_omega = 0;
    class System*   m_system = nullptr;

};

#endif // WAVEFUNCTION_H
