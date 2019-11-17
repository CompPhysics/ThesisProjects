#include "wavefunction.h"

WaveFunction::WaveFunction(System* system, double omega) {
    m_system = system;
    m_omega = omega;
}

vec WaveFunction::computeHermitePolynomial(int nValue, vec position) {
    // Computes Hermite polynomials.
    double omegaSqrt = sqrt(m_omega);

    vec factor = 2*omegaSqrt*position;

    vec HermitePolynomialPP = zeros(position.size());        // H_{n-2}
    vec HermitePolynomialP = ones(position.size());          // H_{n-1}
    vec HermitePolynomial = HermitePolynomialP;              // H_n

    for (int n=1; n <= nValue; n++) {
        HermitePolynomial = factor%HermitePolynomialP - 2*(n-1)*HermitePolynomialPP;
        HermitePolynomialPP = HermitePolynomialP;
        HermitePolynomialP = HermitePolynomial;
    }

    return HermitePolynomial;
}

//    vec x = position*omegaSqrt;

//    if (nValue == 0) {
//        return ones(x.size());
//    }
//    else if (nValue == 1) {
//        return 2*x;
//    }
//    else if (nValue == 2) {
//        return 4*x%x - 2;
//    }
//    else if (nValue == 3) {
//        return 8*x%x%x - 12*x;
//    }
//    else if (nValue == 4) {
//        return 16*x%x%x%x - 48*x%x + 12;
//    }
//    else if (nValue == 5) {
//        return 32*x%x%x%x%x -160*x%x%x +120*x;
//    }
//    else if (nValue == 6) {
//        return 64*x%x%x%x%x%x - 480*x%x%x%x + 720*x%x - 120;
//    }
//    else if (nValue == 7) {
//        return 128*x%x%x%x%x%x%x - 1344*x%x%x%x%x + 3360*x%x%x - 1680;
//    }
//    else if (nValue == 8) {
//        return 256*x%x%x%x%x%x%x%x - 3584*x%x%x%x%x%x + 13440*x%x%x%x - 13440*x%x + 1680;
//    }
