#include <iostream>
#include <armadillo>
#include <time.h>
#include <cassert>
#include "../system.h"
#include "../Math/factorial.h"
#include "../WaveFunctions/wavefunction.h"
#include "../WaveFunctions/doublewell.h"
#include "UnitTest++/UnitTest++.h"
#include "wrapper.h"

using namespace std;
using namespace arma;
using namespace UnitTest;

SUITE(Diagonalization) {
    // Check that factorial function is performing:
    TEST(Factorial) {
        CHECK_EQUAL(87178291200, factorial(14));
    }
    Wrapper* wrapper = new Wrapper();

    TEST(Initialization) {
        wrapper->m_N                    = 1000;
        wrapper->m_posMin               = -10;
        wrapper->m_posMax               = 10;
        wrapper->m_omega_r              = 0.5;
        wrapper->m_nMax                 = 3;
        wrapper->m_nPrimeMax            = 2;
        wrapper->m_numberOfDimensions   = 3;
        wrapper->m_h                    = (wrapper->m_posMax-wrapper->m_posMin)/wrapper->m_N;

        wrapper->m_harmOscPotential     = true;

        vec L(wrapper->m_numberOfDimensions);
        L.fill(0.);
        L(0) = 5.;
        wrapper->setL(L);

        if (wrapper->m_numberOfDimensions == 2) {
            //numberOfEigstates = int(0.5*(nMax+1)*(nMax+2));
            wrapper->m_numberOfEigstates = int(0.5*(wrapper->m_nMax)*(wrapper->m_nMax+1));
        }

        else if (wrapper->m_numberOfDimensions == 3) {
            //numberOfEigstates = int((nMax+1)*(nMax+2)*(nMax+3)/6.);
            wrapper->m_numberOfEigstates = int((wrapper->m_nMax)*(wrapper->m_nMax+1)*(wrapper->m_nMax+2)/6.);
        }
        else { wrapper->m_numberOfEigstates = wrapper->m_nMax; }

        assert(wrapper->m_nPrimeMax <= wrapper->m_numberOfEigstates);
    }

    TEST(McTestFace) {
        //Set up the vector x and the matrix A:
        mat r = zeros(wrapper->m_N+1,wrapper->m_numberOfDimensions);
        vec rAbs = zeros(wrapper->m_N+1);

        for (int d = 0; d < wrapper->m_numberOfDimensions; d++) {
            r.col(d) = wrapper->m_posMin + linspace(0, wrapper->m_N, wrapper->m_N+1)*wrapper->m_h;
            rAbs += r.col(d)%r.col(d);
        }

        rAbs = sqrt(rAbs);

        cube diagMat(wrapper->m_N-1, wrapper->m_N-1, wrapper->m_numberOfDimensions);
        cube eigvecs(wrapper->m_N-1, wrapper->m_N-1, wrapper->m_numberOfDimensions);
        mat  eigvals(wrapper->m_N-1, wrapper->m_numberOfDimensions);

        mat saveEigenvector         = ones(wrapper->m_N-1, wrapper->m_numberOfEigstates);
        cube saveSepEigenvector		= zeros(wrapper->m_N-1, wrapper->m_numberOfEigstates, wrapper->m_numberOfDimensions);
        mat SavePositionvector      = zeros(wrapper->m_N-1, wrapper->m_numberOfDimensions+1);

        for (int d = 0; d < wrapper->m_numberOfDimensions; d++) {
            SavePositionvector.col(d)   = r.col(d).subvec(1, wrapper->m_N-1);    //Saves the y vector for output.
        }

        SavePositionvector.col(wrapper->m_numberOfDimensions) = rAbs.subvec(1, wrapper->m_N-1);    //Saves the r vector for output.

        //Init system
        System* system = new System(wrapper->m_omega_r, wrapper->m_numberOfDimensions, wrapper->m_h, wrapper->m_N);

        system->setWaveFunction(new DoubleWell(system, wrapper->m_omega_r));

        system->diagonalizeMatrix(r, wrapper->m_L, wrapper->m_N, diagMat);
        system->findEigenstate(eigvals, eigvecs, diagMat,
                               saveEigenvector, saveSepEigenvector,
                               wrapper->m_numberOfEigstates, wrapper->m_nMax);

        // Check to see if eigenvalues correspond to the first values for the given potential.
        if (wrapper->m_harmOscPotential) {
            CHECK_CLOSE(0.5, eigvals(0,0), 0.0001);
        }
        // Check that quantum numbers are correct
        CHECK_EQUAL(wrapper->m_nPrimeMax, system->getQuantumNumbers()(wrapper->m_numberOfEigstates - 1, 0));

        if (wrapper->m_omega_r == 0.5 && wrapper->m_L(0) == 5) {
            // BEGIN harmOsc: Check that the potential returns correct value
            vec tmp_resPot(3);
            tmp_resPot(0) = 6.25; tmp_resPot(1) = 6.25; tmp_resPot(2) = 1.5625;

            int tmp_L = 5.;

            vec tmp_r(3);
            tmp_r(0) = -10; tmp_r(1) = 0; tmp_r(2) = 2.5;
            CHECK_ARRAY_EQUAL(system->getWaveFunction()->potential(tmp_r, tmp_L), tmp_resPot, 3);
            tmp_L = 0.;
            tmp_resPot(0) = 25; tmp_resPot(1) = 0; tmp_resPot(2) = 1.5625;
            CHECK_ARRAY_EQUAL(system->getWaveFunction()->potential(tmp_r, tmp_L), tmp_resPot, 3);
            // END harmOsc

            // BEGIN wavefunc: Check that the wavefunction returns correct value for a given quantum number
            double tmp_resWave1 = 0.0045717;
            double tmp_resWave2 = -0.3826216;
            vec tmp_r2(1);
            tmp_r2(0) = -4.44; //"Random" coordinate.
            CHECK_CLOSE(system->getWaveFunction()->harmonicOscillatorBasis(tmp_r2, 0)(0), tmp_resWave1, 0.0001);
            CHECK_CLOSE(system->getWaveFunction()->harmonicOscillatorBasis(tmp_r2, 5)(0), tmp_resWave2, 0.0001);
            // END wavefunc
        }




        cout << endl << "eigvals, Armadillo:" << endl;
        int displayVals = 15;
        for (int i = 0; i < displayVals; ++i) {
            for (int d = 0; d < wrapper->m_numberOfDimensions; d++) {
                cout << i+1 << ": E"<< d <<": " << eigvals.col(d)(i);
            }
            cout << endl;
        }
        cout << endl;
        cout << "Computation time (sec):" << endl;
        cout << "Aramadillo: " << system->getComputationTime() << endl;

    }
}

int main() {
    return RunAllTests();
}
