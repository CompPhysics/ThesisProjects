#ifndef FINITEWELL_H
#define FINITEWELL_H

#include "wavefunction.h"

class FiniteWell : public WaveFunction {
public:
    FiniteWell(class System* system, double omega, double distanceToWall);
    vec harmonicOscillatorBasis(mat x, int n);
    vec potential (vec r, double L);
private:
    double m_distanceToWall = 0;
};


#endif // FINITEWELL_H
