#include "finitewell.h"
#include "../system.h"
#include "../Math/factorial.h"
#include <cmath>

FiniteWell::FiniteWell(System *system, double omega, double distanceToWall)
    : WaveFunction(system, omega) {

    m_distanceToWall = distanceToWall;
}

vec FiniteWell::harmonicOscillatorBasis(mat x, int n) {

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

vec FiniteWell::potential (vec r, double L) {
    int N = m_system->getN();
    vec V(N+1);
    double distanceToWall2 = m_distanceToWall*m_distanceToWall;

    for (int i = 0; i < N+1; i++) {
        double r_i2 = r[i]*r[i];
        if (r_i2 > distanceToWall2) { V[i] = distanceToWall2; }
        else { V[i] = r_i2; }
    }

    return 0.5*m_omega*m_omega*V;
}
