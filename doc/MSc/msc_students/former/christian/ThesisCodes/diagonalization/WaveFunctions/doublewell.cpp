#include "doublewell.h"
#include "../system.h"
#include "../Math/factorial.h"
#include <cmath>

using namespace std;

DoubleWell::DoubleWell(System *system, double omega)
    : WaveFunction(system, omega) {

}

vec DoubleWell::harmonicOscillatorBasis(mat x, int n) {

    double nFac = factorial(n);

    double n2 = pow(2., n);
    double pi4 = pow(M_PI, -0.25);
    double omega4 = pow(m_omega, 0.25);
    double constant = omega4*pi4/sqrt(nFac*n2);
    vec xAbs2 = x%x;

    vec wavefunc = exp(-0.5*m_omega*xAbs2);

    vec phi = constant*wavefunc%computeHermitePolynomial(n, x);

    return phi;
}

vec DoubleWell::potential (vec r, double L) {
    return 0.5*m_omega*m_omega*(r%r - 2*abs(r)*L + L*L);
}
