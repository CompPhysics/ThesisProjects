#include "system.h"
#include "WaveFunctions/wavefunction.h"


System::System(double omega, int numberOfDimensions, double h, int N, bool createSupPos) {
    m_omega = omega;
    setNumberOfDimensions(numberOfDimensions);
    setN(N);
    setStepLength(h);
    m_createSupPos = createSupPos;
}

void System::diagonalizeMatrix(mat r, vec L, int N, cube &diagMat, mat &savePotential) {
    double Constant = 1./(2*m_h*m_h);
    mat V(N+1, m_numberOfDimensions);
    for (int d = 0; d < m_numberOfDimensions; d++) {
        V.col(d) = m_waveFunction->potential(r.col(d), L(d));
        diagMat.slice(d).diag(0)  =  2.*Constant + V.col(d).subvec(1,N-1);     //Set d_i elements in A
        diagMat.slice(d).diag(1)  = -1.*Constant*ones(N-2);               //Set e_i elements in A
        diagMat.slice(d).diag(-1) = diagMat.slice(d).diag(1);                         //Set e_i elements in A
    }

    savePotential = V;
    return;
}

void System::findEigenstate(mat &eigvals, cube eigvecs, cube diagMat,
                            mat &saveEigenvector,
                            cube &saveSepEigenvector,
                            int numberOfEigstates, int nMax) {
    clock_t start1, finish1;
    start1 = clock();

    // Finding eigenvalues and eigenvectors using armadillo:
    for (int d = 0; d < m_numberOfDimensions; d++) {
        vec eigvalsTemp = eigvals.col(d);
        eig_sym(eigvalsTemp, eigvecs.slice(d),  diagMat.slice(d));
        eigvals.col(d) = eigvalsTemp;
    }

    cube eigVecsTemp = zeros(m_N-1, m_N-1/*numberOfEigstates*/, m_numberOfDimensions);
    m_qNumbers = zeros(numberOfEigstates, m_numberOfDimensions);
    m_numberOfEigstates = numberOfEigstates;

    if (m_numberOfDimensions == 1) {
        for (int i = 0; i < numberOfEigstates; i++) {
            m_qNumbers(i, 0) = i;
        }
        int d = 0;
        eigVecsTemp.slice(d) = eigvecs.slice(d).submat(0,0,m_N-2,m_N-2/*numberOfEigstates-1*/);
        saveEigenvector = eigVecsTemp.slice(d);
        saveSepEigenvector.slice(d) = eigVecsTemp.slice(d);
    }


    else if (m_numberOfDimensions == 2) {
        int nx 	= 0;
        int ny 	= 0;
        int n 	= 0;
        int i 	= 0;

        for (int d = 0; d < m_numberOfDimensions; d++) {
            eigVecsTemp.slice(d) = eigvecs.slice(d);//.submat(0,0,m_N-2,m_N-2/*numberOfEigstates-1*/);
            saveSepEigenvector.slice(d) = eigVecsTemp.slice(d);
        }

        while (n < nMax) {
            m_qNumbers(i, 0) = nx;
            m_qNumbers(i, 1) = ny;
            saveEigenvector.col(i) = eigVecsTemp.slice(0).col(nx)%eigVecsTemp.slice(1).col(ny);

            cout << "i: " << i << " nx: " << nx << " ny: " << ny << " n: " << n << endl;

            if (ny == n) {
                n++;
                nx = n;
                ny = 0;
            }
            else {
                nx--;
                ny++;
            }

            i++;
        }
    }

    else if (m_numberOfDimensions == 3) {

        int i = 0;
        for (int nx = 0; nx < nMax; nx++) {
            for (int ny = 0; ny < nMax; ny++) {
                for (int nz = 0; nz < nMax; nz++) {
                    if (nx+ny+nz < nMax) {
                        m_qNumbers(i,0) = nx;
                        m_qNumbers(i,1) = ny;
                        m_qNumbers(i,2) = nz;
                        i++;
                    }
                }
            }
        }

//        m_qNumbers(0, 0) = 0; m_qNumbers(0, 1) = 0; m_qNumbers(0, 2) = 0;

//        if (numberOfEigstates >= 4) {
//            m_qNumbers(1, 0) = 1; m_qNumbers(1, 1) = 0; m_qNumbers(1, 2) = 0;
//            m_qNumbers(2, 0) = 0; m_qNumbers(2, 1) = 1; m_qNumbers(2, 2) = 0;
//            m_qNumbers(3, 0) = 0; m_qNumbers(3, 1) = 0; m_qNumbers(3, 2) = 1;
//        }

//        if (numberOfEigstates >= 10) {
//            m_qNumbers(4, 0) = 2; m_qNumbers(4, 1) = 0; m_qNumbers(4, 2) = 0;
//            m_qNumbers(5, 0) = 1; m_qNumbers(5, 1) = 1; m_qNumbers(5, 2) = 0;
//            m_qNumbers(6, 0) = 1; m_qNumbers(6, 1) = 0; m_qNumbers(6, 2) = 1;
//            m_qNumbers(7, 0) = 0; m_qNumbers(7, 1) = 2; m_qNumbers(7, 2) = 0;
//            m_qNumbers(8, 0) = 0; m_qNumbers(8, 1) = 1; m_qNumbers(8, 2) = 1;
//            m_qNumbers(9, 0) = 0; m_qNumbers(9, 1) = 0; m_qNumbers(9, 2) = 2;
//        }

        for (int d = 0; d < m_numberOfDimensions; d++) {
            eigVecsTemp.slice(d) = eigvecs.slice(d).submat(0,0,m_N-2,m_N-2/*numberOfEigstates-1*/);
            saveSepEigenvector.slice(d) = eigVecsTemp.slice(d);
        }


        for (int i = 0; i < numberOfEigstates; i++) {
            int nx = m_qNumbers(i, 0);
            int ny = m_qNumbers(i, 1);
            int nz = m_qNumbers(i, 2);

            saveEigenvector.col(i) = eigVecsTemp.slice(0).col(nx)
                                    %eigVecsTemp.slice(1).col(ny)
                                    %eigVecsTemp.slice(2).col(nz);
        }
    }

    else {
        cout << "Number of dimensions must be 1, 2 or 3." << endl;
    }

//    mat temp = zeros(numberOfEigstates, m_numberOfDimensions);
//    for (int p = 0; p < numberOfEigstates; p+=2) {
//        for (int d = 0; d < m_numberOfDimensions; d++) {
//            temp(p, d) = m_qNumbers(p/2, d);
//            if (p+1 < numberOfEigstates) {
//                temp(p+1, d) = m_qNumbers(p/2, d);
//            }
//        }
//    }
//    m_qNumbers = temp;
//    cout << m_qNumbers << endl;



    m_psi = saveSepEigenvector;

    m_eigvals = eigvals;

    finish1 = clock();
    m_computationTime = ((finish1-start1)/(double) CLOCKS_PER_SEC);

    return;
}

void System::findCoefficients(int nMax, int nPrimeMax, vec x, mat &C, int currentDim){
    cout << "Finding coefficients for dimension " << currentDim+1 << " of " <<  m_numberOfDimensions << endl;
    cout.flush();
    std::string upLine = "\033[F";
    for	(int nPrime = 0; nPrime < nPrimeMax; nPrime++) {
        cout << "nPrime = " << nPrime << " of " << nPrimeMax-1 << endl;
        for (int nx = 0; nx < nMax; nx++) {
            cout << "[" << int(double(nx)/nMax * 100.0) << " %]\r";
            cout.flush();
            double innerprod = 0;
            for (int i = 0; i < m_N-1; i++) {
                innerprod += m_psi.slice(currentDim).col(nPrime)(i)*m_waveFunction->harmonicOscillatorBasis(x, nx)(i);
            }
            C(nx, nPrime) = innerprod;
        }
        cout << upLine;
        //nPrime++;   //Only need even nPrimes due to double well (degeneracy = 2).
    }
    cout << upLine;

//    for	(int nPrime = 0; nPrime < nPrimeMax; nPrime++) {
//        for (int nx = 0; nx < nMax; nx+=2) {
//            C(nx, nPrime) = C(nx/2, nPrime);
//            if (nx+1 < nMax) {
//                C(nx+1, nPrime) = C(nx/2, nPrime);
//            }
//        }
//    }

    C *= m_h;

//    cout << endl << endl << endl;
//    double test = 0;
//    for (int i = 0; i < m_N-1; i++) {
//        test += (m_psi.slice(currentDim).col(0)(i) - C(0,0)*m_waveFunction->harmonicOscillatorBasis(x, 0)(i))
//               *(m_psi.slice(currentDim).col(0)(i) - C(0,0)*m_waveFunction->harmonicOscillatorBasis(x, 0)(i));
//    }
//    cout << sqrt(test) << endl;

//    double p = 0;
//    for (int nx = 0; nx < nMax; nx++) {
//        p += C(nx, 0)*C(nx, 0);
//    }
//    cout << p << endl;

//    double a = 1;
//    while (a > 0) {
//        double b = 0;
//    }
}

mat System::findSuperPos(mat r, int nMax, int nPrimeMax, cube &supPosSep, cube &saveC) {


//    for (int i = 0; i < m_N-10; i++) {
//    cout << m_psi.slice(0).col(0)(i) << "    " << r(i) << endl;
//    }

//    int i = 1;
//    while (i > 0) {
//        int a = 0;
//    }

    mat rCut = zeros(m_N-1, m_numberOfDimensions);

    for (int d=0; d < m_numberOfDimensions; d++) {
        rCut.col(d) = r.col(d).subvec(1,m_N-1);
    }

    //rCut.col(1) = zeros(m_N-1);
    //rCut.col(2) = zeros(m_N-1);

    cube C = zeros(nMax, nPrimeMax, m_numberOfDimensions);
    mat supPos = zeros(m_N-1, nPrimeMax);

    //This is wrong!:
    for (int d = 0; d < m_numberOfDimensions; d++) {
        findCoefficients(nMax, nPrimeMax, rCut.col(d), C.slice(d), d);
        //C.slice(d) = C.slice(d)%C.slice(d)/dot(C.slice(d),C.slice(d));
        //saveC %= C.slice(d);
    }
    saveC = C;

//    for (int nPrime = 0; nPrime < nPrimeMax; nPrime++) {
//        for (int n = 0; n < nMax; n++) {
//            vec plusTerm = ones(m_N-1);
//            for (int d = 0; d < m_numberOfDimensions; d++) {
//                plusTerm %= C(n, nPrime, d)*m_waveFunction->harmonicOscillatorBasis(rCut.col(d), n);
//            }
//            supPos.col(nPrime) += plusTerm;
//        }
//        nPrime++;   //Only need even nPrimes due to double well (degeneracy = 2).
//    }

    if (m_createSupPos) {
    if (m_numberOfDimensions == 1) {
        for (int nPrime = 0; nPrime < nPrimeMax; nPrime++) {
            for (int i = 0; i < m_numberOfEigstates; i++) {
                int nx = m_qNumbers(i, 0);

                vec plusTerm = C(nx, nPrime, 0)*m_waveFunction->harmonicOscillatorBasis(rCut.col(0), nx);

                supPos.col(nPrime) += plusTerm;
                supPosSep.slice(0).col(nPrime) += plusTerm;
            }
        }

//        for (int nPrime = 0; nPrime < nPrimeMax; nPrime++) {
//            for (int d = 0; d < m_numberOfDimensions; d++) {
//                vec plusTerm = zeros(m_N-1);
//                for (int n = 0; n < nMax; n++) {
//                    plusTerm += C(n, nPrime, d)*m_waveFunction->harmonicOscillatorBasis(rCut.col(d), n);
//                }
//                supPos.col(nPrime) %= plusTerm;
//                supPosSep.slice(d).col(nPrime) = plusTerm;
//            }
//            //nPrime++;   //Only need even nPrimes due to double well (degeneracy = 2).
//        }
    }

    else if (m_numberOfDimensions == 2) {

        for (int nPrime = 0; nPrime < nPrimeMax; nPrime++) {

            for (int i = 0; i < m_numberOfEigstates; i++) {
                int nx = m_qNumbers(i, 0);
                int ny = m_qNumbers(i, 1);

                vec plusTermX = C(nx, nPrime, 0)*m_waveFunction->harmonicOscillatorBasis(rCut.col(0), nx);
                vec plusTermY = C(ny, nPrime, 1)*m_waveFunction->harmonicOscillatorBasis(rCut.col(1), ny);

                supPos.col(nPrime) += plusTermX%plusTermY;
                //supPosSep.slice(0).col(nPrime) += plusTermX;
                //supPosSep.slice(1).col(nPrime) += plusTermY;
            }

            for (int n = 0; n < nMax; n++) {

                vec plusTermXSep = C(n, nPrime, 0)*m_waveFunction->harmonicOscillatorBasis(rCut.col(0), n);
                vec plusTermYSep = C(n, nPrime, 1)*m_waveFunction->harmonicOscillatorBasis(rCut.col(1), n);

                supPosSep.slice(0).col(nPrime) += plusTermXSep;
                supPosSep.slice(1).col(nPrime) += plusTermYSep;

            }
        }

//        int nPrime = 0;
//        int nx 	= 0;
//        int ny 	= 0;
//        int n 	= 0;
//        int i 	= 0;

//        vec plusTermX = zeros(m_N-1);
//        vec plusTermY = zeros(m_N-1);
//        vec plusTerm = zeros(m_N-1);
//        while (n < nMax) {
//            plusTermX = C(nx, nPrime, 0)*m_waveFunction->harmonicOscillatorBasis(rCut.col(0), nx);
//            plusTermY = C(ny, nPrime, 1)*m_waveFunction->harmonicOscillatorBasis(rCut.col(1), ny);
//            plusTerm += plusTermX%plusTermY;
//            supPosSep.slice(0).col(nPrime) += plusTermX;
//            supPosSep.slice(1).col(nPrime) += plusTermY;

//            cout << "i: " << i << " nx: " << nx << " ny: " << ny << " n: " << n << endl;

//            if (ny == n) {
//                n++;
//                nx = n;
//                ny = 0;
//            }
//            else {
//                nx--;
//                ny++;
//            }

//            i++;
//        }
//        //supPosSep.slice(0).col(nPrime) /= sqrt(dot(supPosSep.slice(0).col(nPrime),supPosSep.slice(0).col(nPrime)));
//        //supPosSep.slice(1).col(nPrime) /= sqrt(dot(supPosSep.slice(1).col(nPrime),supPosSep.slice(1).col(nPrime)));
//        supPos.col(nPrime) = plusTerm;
    }

    else if (m_numberOfDimensions == 3) {

        for (int nPrime = 0; nPrime < nPrimeMax; nPrime++) {
            for (int i = 0; i < m_numberOfEigstates; i++) {
                int nx = m_qNumbers(i, 0);
                int ny = m_qNumbers(i, 1);
                int nz = m_qNumbers(i, 2);

                vec plusTermX = C(nx, nPrime, 0)*m_waveFunction->harmonicOscillatorBasis(rCut.col(0), nx);
                vec plusTermY = C(ny, nPrime, 1)*m_waveFunction->harmonicOscillatorBasis(rCut.col(1), ny);
                vec plusTermZ = C(nz, nPrime, 2)*m_waveFunction->harmonicOscillatorBasis(rCut.col(2), nz);

                supPos.col(nPrime) += plusTermX%plusTermY%plusTermZ;
                //supPosSep.slice(0).col(nPrime) += plusTermX;
                //supPosSep.slice(1).col(nPrime) += plusTermY;
                //supPosSep.slice(2).col(nPrime) += plusTermZ;
            }

            for (int n = 0; n < nMax; n++) {

                vec plusTermXSep = C(n, nPrime, 0)*m_waveFunction->harmonicOscillatorBasis(rCut.col(0), n);
                vec plusTermYSep = C(n, nPrime, 1)*m_waveFunction->harmonicOscillatorBasis(rCut.col(1), n);
                vec plusTermZSep = C(n, nPrime, 2)*m_waveFunction->harmonicOscillatorBasis(rCut.col(2), n);

                supPosSep.slice(0).col(nPrime) += plusTermXSep;
                supPosSep.slice(1).col(nPrime) += plusTermYSep;
                supPosSep.slice(2).col(nPrime) += plusTermZSep;

            }
        }
    }

    else {
        cout << "Number of dimensions must be 1, 2 or 3." << endl;
    }
    }
    return supPos;

}


void System::setWaveFunction(WaveFunction *waveFunction) { m_waveFunction = waveFunction; }
void System::setN(double N) {m_N = N; }
void System::setStepLength(double h) { m_h = h; }
void System::setNumberOfDimensions(int numberOfDimensions) { m_numberOfDimensions = numberOfDimensions; }
void System::setNumberOfParticles(int numberOfParticles) { m_numberOfParticles = numberOfParticles; }
void System::setComputationTime(double computationTime) { m_computationTime = computationTime; }
