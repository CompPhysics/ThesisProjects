#ifndef SYSTEM_H
#define SYSTEM_H

#include <armadillo>

using namespace arma;

class System {
public:
    System(double omega, int numberOfDimensions, double h, int N, bool createSupPos);
    void diagonalizeMatrix(mat r, vec L, int N, cube &diagMat, mat &savePotential);
    void findEigenstate(mat &eigvals, cube eigvecs, cube diagMat, mat &saveEigenvector, cube &saveSeparateEigenvector, int numberOfEigstates, int nMax);
    void findCoefficients(int nMax, int nPrimeMax, vec x, mat &C, int currentDim);
    mat findSuperPos(mat r, int nMax, int nPrimeMax, cube &supPosSep, cube &saveC);

    void setN(double N);
    void setStepLength(double h);
    void setNumberOfDimensions(int numberOfDimensions);
    void setNumberOfParticles(int numberOfParticles);
    void setComputationTime(double computationTime);
    void setWaveFunction(class WaveFunction* waveFunction);

    class WaveFunction* getWaveFunction()   { return m_waveFunction; }
    double getN()                             { return m_N; }
    double getStepLength()                  { return m_h; }
    int getNumberOfDimensions()             { return m_numberOfDimensions; }
    int getNumberOfParticles()              { return m_numberOfParticles; }
    double getComputationTime()             { return m_computationTime; }
    mat getEigvals()                        { return m_eigvals; }

private:
    int                     m_N                      = 0;
    int                     m_numberOfEigstates      = 0;
    double                  m_h                      = 0;
    int                     m_numberOfDimensions     = 0;
    int                     m_numberOfParticles      = 0;
    double                  m_computationTime        = 0;
    double                  m_omega                  = 0;
    mat						m_qNumbers;
    mat                     m_eigvals;
    cube                    m_psi;
    class WaveFunction*     m_waveFunction           = nullptr;
    bool                    m_createSupPos           = false;
};


#endif // SYSTEM_H
