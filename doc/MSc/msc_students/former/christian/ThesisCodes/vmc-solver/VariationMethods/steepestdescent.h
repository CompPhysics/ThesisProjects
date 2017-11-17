#ifndef PROJECT2_STEEPESTDESCENT_H
#define PROJECT2_STEEPESTDESCENT_H
#include "iostream"
#include <vector>

class SteepestDescent {
public:
    SteepestDescent(class System* system, double stepLengthSD);
    void obtainOptimalParameter(std::vector<double> parameters, double tol,
                                int maxIterations, int numberOfMetropolisSteps, bool importanceSampling);

private:
    double m_stepLengthSD = 0;
    std::vector<double> m_derivativeAvg;
    class System* m_system = nullptr;
};

#endif // PROJECT2_STEEPESTDESCENT_H
