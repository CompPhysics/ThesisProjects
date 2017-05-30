#ifndef QUANTUMSTATE_H
#define QUANTUMSTATE_H
#include <armadillo>
using namespace arma;

class QuantumState
{
public:
    QuantumState();
    void setQuantumNumbers(int n, int m, int sm, double s);
    void flipSpin();
    vec getQuantumNumbers() {return {(double)m_n, (double)m_m, (double)m_sm, m_s} ;}

private:
    //Quantum numbers n, m, sm (spin projection), s.
    int m_n;
    int m_m;
    int m_sm;
    double m_s;
};

#endif // QUANTUMSTATE_H
