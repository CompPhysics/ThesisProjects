#ifndef PROJECT2_HAMILTONIAN_H
#define PROJECT2_HAMILTONIAN_H
#include <vector>
#include <armadillo>
#include "../HermitePolynomials/hermitepolynomials.h"

using namespace arma;

class Hamiltonian {
public:
    Hamiltonian(class System* system, bool analyticalKinetic, double alpha, double omega);
    double computeKineticEnergy(std::vector<class Particle*> particles);
    void setUpHermitePolynomials();
    virtual std::vector<double> computeLocalEnergy(std::vector<class Particle*> particles) = 0;// { particles = particles; }
    virtual double evaluateSingleParticleWF(vec n, std::vector<double> r, int j) { return n[0]*r[0]*j; }
    virtual std::vector<double> computeSPWFDerivative(vec n, std::vector<double> r, int j) { j = j; n = n; return r; }
    virtual double computeSPWFDoubleDerivative(vec n, std::vector<double> r, int j) { return n[0]*r[0]*j; }
    virtual double computeSPWFAlphaDerivative(vec n, std::vector<double> r, int j) { return n[0]*r[0]*j; }
    double getAnalytic(){ return m_analyticalKinetic; }
    void setExpFactor(double expFactor) { m_expFactor = expFactor; }
    void setAlpha(double alpha) { m_alpha =  alpha; }
    virtual vec getWellDistance()   { return zeros(1); }
    HermitePolynomials** getHermitePolynomials() {return m_hermitePolynomials;}
    HermitePolynomials** getHermitePolynomialsDerivative() {return m_hermitePolynomialsDerivative;}
    HermitePolynomials** getHermitePolynomialsDoubleDerivative() {return m_hermitePolynomialsDoubleDerivative;}

    virtual double computeHermitePolynomial(int nValue, double position) { return nValue*position; }
    virtual double computeHermitePolynomialDerivative(int nValue, double position) { return nValue*position; }
    virtual double computeHermitePolynomialDoubleDerivative(int nValue, double position) { return nValue*position; }
    virtual double computeHermitePolynomialAlphaDerivative(int nValue, double position) { return nValue*position; }

protected:
    class System* m_system = nullptr;
    bool m_analyticalKinetic = false;
    double m_expFactor = 0;
    double m_alpha = 0;
    double m_omega = 0;

    HermitePolynomials** m_hermitePolynomials;
    HermitePolynomials** m_hermitePolynomialsDerivative;
    HermitePolynomials** m_hermitePolynomialsDoubleDerivative;
};

#endif // PROJECT2_HAMILTONIAN_H
