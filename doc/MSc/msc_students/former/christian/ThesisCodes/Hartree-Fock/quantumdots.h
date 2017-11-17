#ifndef QUANTUMDOTS_H
#define QUANTUMDOTS_H
#include <armadillo>
#include "quantumstate.h"

using namespace arma;


class QuantumDots
{
public:
    QuantumDots(int numberOfShells, int numberOfParticles, int numberOfDimensions, double omega);
    void setUpBasisStates();
    void calculateUnperturbedEnergiesHO();
    void create4DMatrix(double ****&matrix, int dim1, int dim2, int dim3, int dim4);
    void fillTwoBodyMatrix();

    mat computeDensityMatrix();
    void computeHartreeFockMatrix(mat densityMatrix);
    double computeHartreeFockEnergyDifference();
    void computeHartreeFockEnergy(mat densityMatrix);
    void runHartreeFock();

private:
    int m_numberOfShells = 0;
    int m_numberOfParticles = 0;
    int m_numberOfDimensions = 0;
    int m_numberOfStates = 0;
    double m_omega = 0;
    double**** m_twoBodyMatrix;
    vec m_unperturbedEnergiesHO;
    vec m_eigvals;
    vec m_eigvalsOld;
    mat m_eigvecs;
    mat m_coefficients;
    mat m_HartreeFockMatrix;

    std::vector<QuantumState> m_quantumStates;
};

#endif // QUANTUMDOTS_H
