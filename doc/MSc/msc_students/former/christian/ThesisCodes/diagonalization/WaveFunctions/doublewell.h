#ifndef DOUBLEWELL_H
#define DOUBLEWELL_H

#include "wavefunction.h"

class DoubleWell : public WaveFunction {
public:
    DoubleWell(class System* system, double omega);
    vec harmonicOscillatorBasis(mat x, int n);
    vec potential (vec r, double L);
};

#endif // DOUBLEWELL_H
