#include "squarewell.h"
#include "../system.h"
#include "../Math/factorial.h"
#include <cmath>

using namespace std;

SquareWell::SquareWell(System *system, double omega, double V0, double distanceToWall)
    : WaveFunction(system, omega) {

    m_V0 = V0;
    m_distanceToWall = distanceToWall;
}

vec SquareWell::harmonicOscillatorBasis(mat x, int n) {

    double nFac = factorial(n);

    double n2 = pow(2., n);
    double pi4 = pow(M_PI, -0.25);
    double omega4 = pow(m_omega, 0.25);
    double constant = omega4*pi4/sqrt(nFac*n2);

    vec xAbs2 = x%x;

    vec wavefunc = exp(-0.5*m_omega*xAbs2);

    vec phi = constant*wavefunc%computeHermitePolynomial(n, x);

//    mat eigvals = m_system->getEigvals();
//    double E = eigvals.col(0)[n];
//    double k = sqrt(2*E);
//    double kPrime = sqrt(2*(E-m_V0));
//    double alpha = sqrt(2*(m_V0-E));

//    double A;
//    double B;
//    double F;
//    double G;
//    double H;
//    double I;

//    if (n%2 == 0) {
//        A = 0.;
//        B = 1.;
//        F = 0.;
//        G = B*cos(k*m_distanceToWall)/exp(-alpha*m_distanceToWall);
//        H = G;
//        I = 0.;
//    }

//    else {
//        A = 1.;
//        B = 0.;
//        F = 0.;
//        G = -A*sin(k*m_distanceToWall)/exp(-alpha*m_distanceToWall);
//        H = -G;
//        I = 0.;
//    }

//    int N = m_system->getN();
//    vec phi(N-1);
//    for (int i = 0; i < N-1; i++) {
//        if (x[i] < -m_distanceToWall) {
//            if (E > m_V0) {
//                phi[i] = sin(kPrime*x[i]) + cos(kPrime*x[i]);
//            }
//            else {
//                phi[i] = F*exp(-alpha*x[i]) + G*exp(alpha*x[i]);
//            }
//        }
//        else if (x[i] > m_distanceToWall) {
//            if (E > m_V0) {
//                phi[i] = sin(kPrime*x[i]) + cos(kPrime*x[i]);
//            }
//            else {
//                phi[i] = H*exp(-alpha*x[i]) + I*exp(alpha*x[i]);
//            }
//        }
//        else {
//            phi[i] = A*sin(k*x[i]) + B*cos(k*x[i]);
//        }
//    }

    return phi;
}

vec SquareWell::potential (vec r, double L) {
    int N = m_system->getN();
    vec V(N+1);

    for (int i = 0; i < N+1; i++) {
        if (abs(r[i]) >= m_distanceToWall) { V[i] = m_V0; }
        else { V[i] = 0; }
    }

    return V;
}
