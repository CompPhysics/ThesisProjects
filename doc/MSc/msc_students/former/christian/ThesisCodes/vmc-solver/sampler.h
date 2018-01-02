#ifndef PROJECT2_SAMPLER_H
#define PROJECT2_SAMPLER_H
#include <vector>

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void setEnergy(double energy);
    void setKineticEnergy(double kineticEnergy);
    void setPotentialEnergy(double potentialEnergy);
    void setVariance(double variance);
    void setAcceptanceRate(double acceptanceRate);
    void setMeanDistance(double meanDistance);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    void saveToFile(double localEnergy);
    double getEnergy()                   { return m_energy; }
    double getKineticEnergy()            { return m_kineticEnergy; }
    double getPotentialEnergy()          { return m_potentialEnergy; }
    double getMeanDistance()             { return m_meanDistance; }
    double getVariance()                 { return m_variance; }
    double getAcceptanceRate()           { return m_acceptanceRate; }
    std::vector<double> getWaveFuncDerivativeParameters()  { return m_waveFuncDerivativeParameters; }
    std::vector<double> getWaveFuncEnergyParameters()      { return m_waveFuncEnergyParameters; }

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    int     m_cumulativeAcceptedSteps = 0;
    double  m_energy = 0;
    double  m_squaredEnergy = 0;
    double  m_kineticEnergy = 0;
    double  m_potentialEnergy = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeSquaredEnergy = 0;
    double  m_cumulativeKineticEnergy = 0;
    double  m_cumulativePotentialEnergy = 0;
    double  m_variance = 0;
    double  m_acceptanceRate = 0;
    std::vector<double>  m_waveFuncDerivativeParameters = std::vector<double>();
    std::vector<double>  m_waveFuncEnergyParameters = std::vector<double>();
    std::vector<double>  m_cumulativeWFuncDerivativeParameters = std::vector<double>();
    std::vector<double>  m_cumulativeWFuncEnergyParameters = std::vector<double>();
    double  m_cumulativeDistance = 0;
    double  m_meanDistance = 0;
    class System* m_system = nullptr;
};

#endif // PROJECT2_SAMPLER_H
