#include "quantumdots.h"
#include "Coulomb_Functions.hpp"

using namespace std;

bool energyComparison(const pair<int,double> &a,const pair<int,double> &b)
{
       return a.second<b.second;
}

QuantumDots::QuantumDots(int numberOfShells, int numberOfParticles, int numberOfDimensions, double omega)
{
    m_numberOfShells = numberOfShells;
    m_numberOfParticles = numberOfParticles;
    m_numberOfDimensions = numberOfDimensions;
    m_omega = omega;

    setUpBasisStates();


}

void QuantumDots::setUpBasisStates() {

    QuantumState quantumState;
    vector<int> oddNumberedShells;
    int n;
    int m;
    int sm = -1;
    double s = 0.5;

    for (int shell = 1; shell <= m_numberOfShells; shell++) {
        if (shell%2 != 0) oddNumberedShells.push_back(shell);
    }

    for (int i = 0; i < oddNumberedShells.size(); i++) {
        // Negative m (and m=0)
        for (int j = 0; j <= m_numberOfShells - oddNumberedShells[i]; j++) {
            n = i;
            m = -j;
            quantumState.setQuantumNumbers(n, m, sm, s);
            m_quantumStates.push_back(quantumState);
        }

        // Positive m
        for (int j = 1; j <= m_numberOfShells - oddNumberedShells[i]; j++) {
            n = i;
            m = j;
            quantumState.setQuantumNumbers(n, m, sm, s);
            m_quantumStates.push_back(quantumState);
        }
    }

    // Sort quantum states
    pair<int, double> mapping;
    vector<pair<int, double>> sortingVector;

    for (int i = 0; i < m_quantumStates.size(); i++) {
        vec quantumNumbers = m_quantumStates[i].getQuantumNumbers();
        int n = quantumNumbers[0];
        int m = quantumNumbers[1];
        double energyOfState = m_omega*(2.*n + abs(m) + 0.5*m_numberOfDimensions);
        //cout << energyOfState << endl;

        mapping = make_pair(i, energyOfState);
        sortingVector.push_back(mapping);
    }

    sort(sortingVector.begin(), sortingVector.end(), energyComparison);
    vector<QuantumState> sortedStates;

    for (int i = 0; i < sortingVector.size(); i++) {
        sortedStates.push_back(m_quantumStates[sortingVector[i].first]);
        m_quantumStates[sortingVector[i].first].flipSpin();
        sortedStates.push_back(m_quantumStates[sortingVector[i].first]);
    }

    m_quantumStates = sortedStates;
    m_numberOfStates = m_quantumStates.size();

    calculateUnperturbedEnergiesHO();

    int dim = m_numberOfStates;
    create4DMatrix(m_twoBodyMatrix, dim, dim, dim, dim);

    fillTwoBodyMatrix();


//    for (int i = 0; i < m_quantumStates.size(); i++) {
//        cout << m_quantumStates[i].getQuantumNumbers() << endl;
//    }

}

void QuantumDots::calculateUnperturbedEnergiesHO() {
    m_unperturbedEnergiesHO = zeros(m_numberOfStates);
    for(int i = 0; i < m_numberOfStates; i++) {
        vec quantumNumbers = m_quantumStates[i].getQuantumNumbers();
        int n = quantumNumbers[0];
        int m = quantumNumbers[1];
        m_unperturbedEnergiesHO(i) = m_omega*(2.*n + abs(m) + 0.5*m_numberOfDimensions);
    }
}

void QuantumDots::create4DMatrix(double**** &matrix, int dim1, int dim2, int dim3, int dim4) {

    matrix = new double***[dim1];

    for(int i=0; i < dim1; i++) {
        matrix[i] = new double**[dim2];
        for (int j=0; j < dim2; j++) {
            matrix[i][j] = new double*[dim3];
            for(int k = 0; k < dim3; k++) {
                matrix[i][j][k] = new double[dim4];
            }
        }
    }
    for(int i=0; i < dim1; i++) {
        for (int j=0; j < dim2; j++) {
            for(int k = 0; k < dim3; k++) {
                for(int l = 0; l < dim4; l++) {
                    matrix[i][j][k][l] = 0;
                }
            }
        }
    }
}

void QuantumDots::fillTwoBodyMatrix() {
    // Precalculate two body elements.
    double**** twoBodyMatrixPol;
    int dim = m_numberOfStates;
    create4DMatrix(twoBodyMatrixPol, dim, dim, dim, dim);

    for(int i = 0; i < m_numberOfStates; i++) {
        vec quantumNumbersAlpha = m_quantumStates[i].getQuantumNumbers();
        int nAlpha = quantumNumbersAlpha[0];
        int mAlpha = quantumNumbersAlpha[1];
        int smAlpha = quantumNumbersAlpha[2];

        for(int j = 0; j < m_numberOfStates; j++) {
            vec quantumNumbersBeta = m_quantumStates[j].getQuantumNumbers();
            int nBeta = quantumNumbersBeta[0];
            int mBeta = quantumNumbersBeta[1];
            int smBeta = quantumNumbersBeta[2];

            for(int k = 0; k < m_numberOfStates; k++) {
                vec quantumNumbersGamma = m_quantumStates[k].getQuantumNumbers();
                int nGamma = quantumNumbersGamma[0];
                int mGamma = quantumNumbersGamma[1];
                int smGamma = quantumNumbersGamma[2];

                for(int l = 0; l < m_numberOfStates; l++) {
                    vec quantumNumbersDelta = m_quantumStates[l].getQuantumNumbers();
                    int nDelta = quantumNumbersDelta[0];
                    int mDelta = quantumNumbersDelta[1];
                    int smDelta = quantumNumbersDelta[2];

                    if ((smAlpha == smBeta && smGamma == smDelta)){
                        twoBodyMatrixPol[i][k][j][l] = Coulomb_HO(m_omega, nAlpha, mAlpha, nGamma, mGamma, nBeta, mBeta,  nDelta, mDelta);
                    }
                    if ((smAlpha == smDelta && smGamma == smBeta)){
                        twoBodyMatrixPol[i][k][l][j] = Coulomb_HO(m_omega, nAlpha, mAlpha, nGamma, mGamma, nDelta, mDelta, nBeta, mBeta);
                    }
                }
            }
        }
    }

    double frac = 1./4;

    for(int i = 0; i < m_numberOfStates; i++) {

        for(int j = 0; j < m_numberOfStates; j++) {

            for(int k = 0; k < m_numberOfStates; k++) {

                for(int l = 0; l < m_numberOfStates; l++) {

                    m_twoBodyMatrix[i][k][j][l] = /*frac**/(twoBodyMatrixPol[i][k][j][l]/* + twoBodyMatrixPol[i][k-2][j-2][l]
                                                       +twoBodyMatrixPol[i][k-2][j][l-2] + twoBodyMatrixPol[i-2][k][j-2][l]
                                                       +twoBodyMatrixPol[i-2][k][j][l-2] + twoBodyMatrixPol[i-2][k-2][j-2][l-2]*/);
                    //cout << m_twoBodyMatrix[i][j][k][l] << endl;
                }
            }
        }
    }

    int n = m_quantumStates[5].getQuantumNumbers()[0];
    int m = m_quantumStates[5].getQuantumNumbers()[1];

    int ny = (n - m + abs(m))/2.;
    int nx = m + ny;

    cout << nx << " " << ny << endl;





    int i = 4;
    int j = 4;
    int k = 5;
    int l = 5;
    m_twoBodyMatrix[i][k][j][l] = frac*(twoBodyMatrixPol[i][k][j][l] + twoBodyMatrixPol[i][k-2][j-2][l]
                                       +twoBodyMatrixPol[i][k-2][j][l-2] + twoBodyMatrixPol[i-2][k][j-2][l]
                                       +twoBodyMatrixPol[i-2][k][j][l-2] + twoBodyMatrixPol[i-2][k-2][j-2][l-2]);
    cout << m_twoBodyMatrix[i][k][j][l] << endl;
}

mat QuantumDots::computeDensityMatrix() {
    mat densityMatrix = zeros(m_numberOfStates, m_numberOfStates);

    for (int k = 0; k < m_numberOfStates; k++) {
            for (int l = 0; l < m_numberOfStates; l++) {
                for (int i=0; i < m_numberOfParticles; i++) {
                        densityMatrix(k, l) += m_coefficients(k,i)*m_coefficients(l,i);
                }
            }
        }
    return densityMatrix;
}

void QuantumDots::computeHartreeFockMatrix(mat densityMatrix) {
    m_HartreeFockMatrix = zeros(m_numberOfStates,m_numberOfStates);
    double FockElement = 0;

    for(int i = 0; i < m_numberOfStates; i++) {
        vec quantumNumbers = m_quantumStates[i].getQuantumNumbers();
        //int nAlpha = quantumNumbers[0];
        //int mAlpha = quantumNumbers[1];
        int smAlpha = quantumNumbers[2];

        for(int j = 0; j < m_numberOfStates; j++) {
            vec quantumNumbersBeta = m_quantumStates[j].getQuantumNumbers();
            //int nBeta = quantumNumbersBeta[0];
            //int mBeta = quantumNumbersBeta[1];
            int smBeta = quantumNumbersBeta[2];

            for(int k = 0; k < m_numberOfStates; k++) {
                vec quantumNumbersGamma = m_quantumStates[k].getQuantumNumbers();
                //int nGamma = quantumNumbersGamma[0];
                //int mGamma = quantumNumbersGamma[1];
                int smGamma = quantumNumbersGamma[2];

                for(int l = 0; l < m_numberOfStates; l++) {
                    vec quantumNumbersDelta = m_quantumStates[l].getQuantumNumbers();
                    //int nDelta = quantumNumbersDelta[0];
                    //int mDelta = quantumNumbersDelta[1];
                    int smDelta = quantumNumbersDelta[2];

                    double TBME = 0;
                    double tbme1 = 0;
                    double tbme2 = 0;

                    if ((smAlpha == smBeta && smGamma == smDelta)){
                        tbme1 = m_twoBodyMatrix[i][k][j][l];
                    }
                    if ((smAlpha == smDelta && smGamma == smBeta)){
                        tbme2 = m_twoBodyMatrix[i][k][l][j];
                    }

                    TBME = tbme1 - tbme2;
                    FockElement += densityMatrix(k,l)*TBME;
                }
            }
            if (i == j) {
                m_HartreeFockMatrix(i, i) += m_unperturbedEnergiesHO(i);
            }
            m_HartreeFockMatrix(i, j) += FockElement;
            FockElement = 0;
        }
    }
}

double QuantumDots::computeHartreeFockEnergyDifference() {

    return (accu(abs(m_eigvals - m_eigvalsOld))) / (double)m_numberOfStates;
}


void QuantumDots::computeHartreeFockEnergy(mat densityMatrix) {

    int FermiLevel = m_numberOfParticles;
    double selfConsistentFieldIterations = 0;
    double exchangePart = 0;
    double singleParticleEnergies = 0;

    for(int f = 0; f < FermiLevel; f++){
        singleParticleEnergies += m_eigvals(f);
    }

    for(int i = 0; i < m_numberOfStates; i++) {
        vec quantumNumbers = m_quantumStates[i].getQuantumNumbers();
        //int nAlpha = quantumNumbers[0];
        //int mAlpha = quantumNumbers[1];
        int smAlpha = quantumNumbers[2];

        for(int j = 0; j < m_numberOfStates; j++) {
            vec quantumNumbersBeta = m_quantumStates[j].getQuantumNumbers();
            //int nBeta = quantumNumbersBeta[0];
            //int mBeta = quantumNumbersBeta[1];
            int smBeta = quantumNumbersBeta[2];

            for(int k = 0; k < m_numberOfStates; k++) {
                vec quantumNumbersGamma = m_quantumStates[k].getQuantumNumbers();
                //int nGamma = quantumNumbersGamma[0];
                //int mGamma = quantumNumbersGamma[1];
                int smGamma = quantumNumbersGamma[2];

                for(int l = 0; l < m_numberOfStates; l++) {
                    vec quantumNumbersDelta = m_quantumStates[l].getQuantumNumbers();
                    //int nDelta = quantumNumbersDelta[0];
                    //int mDelta = quantumNumbersDelta[1];
                    int smDelta = quantumNumbersDelta[2];

                    double TBME = 0;
                    double tbme1 = 0;
                    double tbme2 = 0;

                    if ((smAlpha == smBeta && smGamma == smDelta)){
                        tbme1 = m_twoBodyMatrix[i][k][j][l];
                    }
                    if ((smAlpha == smDelta && smGamma == smBeta)){
                        tbme2 = m_twoBodyMatrix[i][k][l][j];
                    }

                    TBME = tbme1 - tbme2;
                    selfConsistentFieldIterations = densityMatrix(i,j)*densityMatrix(k,l)*TBME;
                    exchangePart += selfConsistentFieldIterations;
                }
            }
        }
    }

    double HartreeFockEnergy = singleParticleEnergies - 0.5*exchangePart;
    // Uncoment for debug
    //cout << "SPEnergies " << SingleParticleEnergies << endl;
    //cout << "Exchange " << ExchangePart << endl;
    cout << "===================================================================" << endl;
    cout << setprecision(12);
    cout << "Number of electrons: " << m_numberOfParticles << endl;
    cout << "Number of shells: " << m_numberOfShells << endl;
    cout << "Omega: " << m_omega << endl;
    cout << "Total energy: " << HartreeFockEnergy << endl;
    //writeToFile(HF_Energy, NumberOfParticles, m_EnergyCutOff, homega);
}


void QuantumDots::runHartreeFock() {

    m_coefficients = eye(m_numberOfStates, m_numberOfStates);

    double difference = 1.; //dummy value to handle first iteration
    double epsilon = 10e-8;

    m_eigvalsOld = zeros(m_numberOfStates);

    int i = 0;
    while (epsilon < difference && i < 1000){
        mat xDensityMatrix = computeDensityMatrix();
        computeHartreeFockMatrix(xDensityMatrix);
        eig_sym(m_eigvals, m_eigvecs, m_HartreeFockMatrix);

        m_coefficients = m_eigvecs;
        difference = computeHartreeFockEnergyDifference();
        m_eigvalsOld = m_eigvals;
        i++;

    }

    mat yDensityMatrix = computeDensityMatrix();
    computeHartreeFockEnergy(yDensityMatrix);
    cout << "Number of iterations: " << i << endl;
    cube saveC = zeros(m_numberOfStates, m_numberOfStates, 2);
    saveC.slice(0) = m_coefficients;
    saveC.save("../Hartree-Fock/Data/CoefficientsHF.dat", arma_ascii);
}

