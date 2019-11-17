#ifndef PROJECT2_SYSTEM_H
#define PROJECT2_SYSTEM_H
#include <vector>
#include <armadillo>
#include <iostream>

using namespace arma;

class System {
public:
    bool metropolisStep             (int currentParticle);
    bool metropolisStepImpSampling  (int currentParticle);
    void runMetropolisSteps         (int numberOfMetropolisSteps, bool importanceSampling,
                                     bool showProgress, bool printToTerminal);
    void optimizeParameters         (System* system, double alpha, double beta);
    void MPI_CleanUp                (double &totalE, double &totalKE, double &totalPE,
                                     double &totalVariance, double &totalAcceptanceRate,
                                     double &finalMeanDistance, double &timeStart,
                                     double &timeEnd, double &totalTime, int numprocs, int numberOfSteps);
    void mergeOutputFiles           (int numprocs);
    double calculateGreensFunction  (std::vector<double> positionOld, std::vector<double> positionNew);
    std::vector<double> quantumForce();
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setTimeStep                (double dt);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    void setNumberOfMetropolisSteps (int numberOfMetropolisSteps);
    void setCurrentParticle         (int currentParticle);
    void setMyRank                  (int my_rank);
    void setNumProcs                (int numprocs);
    void setComputationTime         (double computationTime);
    void setSaveEnergies            (bool saveEnergies);
    void setSavePositions           (bool savePositions);
    void retrieveCoefficientsFromFile			(std::string fileName, cube &loadCoefficients);
    void retrieveConstantsFromFile              (std::string fileName, vec &loadConstants);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    class InitialState*             getInitialState()   { return m_initialState; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    int getCurrentParticle()             { return m_currentParticle; }
    int getMyRank()                     { return m_my_rank; }
    int getNumProcs()                   { return m_numprocs; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }
    double getComputationTime()         { return m_computationTime; }
    bool getSaveEnergies()              { return m_saveEnergies; }
    FILE* getEnergiesFile()             { return m_outfileE; }
    bool getSavePositions()             { return m_savePositions; }
    FILE* getPositionsFile()            { return m_outfileP; }

    void setDoubleWellFlag(bool doubleWell)                { m_doubleWell = doubleWell; }
    bool getDoubleWellFlag()                               { return m_doubleWell; }

    void setSquareWellFlag(bool squareWell)                { m_squareWell = squareWell; }
    bool getSquareWellFlag()                               { return m_squareWell; }

    void setL(vec L)                                       { m_L = L; }
    vec getL()                                             { return m_L; }

private:
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    int                             m_currentParticle = 0;
    int                             m_my_rank = 0;
    int                             m_numprocs = 0;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;
    double                          m_dt = 0.01;
    double                          m_computationTime = 0;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
    bool                            m_printToTerminal = false;
    bool                            m_saveEnergies = false;
    FILE*                           m_outfileE;
    bool                            m_savePositions = false;
    FILE*                           m_outfileP;
    bool                            m_doubleWell = false;
    bool                            m_squareWell = false;
    vec                             m_L;
};

#endif // PROJECT2_SYSTEM_H
