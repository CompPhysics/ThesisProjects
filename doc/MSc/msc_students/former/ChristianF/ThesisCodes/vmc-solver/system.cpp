#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "VariationMethods/steepestdescent.h"
#include "Math/random.h"
#include <iostream>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <time.h>
#include <mpi.h>

using namespace std;
using namespace arma;

bool System::metropolisStep(int currentParticle) {
    /* Perform the actual Metropolis step: Take the current particle and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */

    // Change position of current particle by a random amount creating a trial state
    setCurrentParticle(currentParticle);
    std::vector<double> positionChange(m_numberOfDimensions);

    for (int i=0; i<m_numberOfDimensions; i++){
        positionChange[i] = (Random::nextDouble()*2-1)*m_stepLength;
    }

    // Metropolis ratio
    double qratio = m_waveFunction->computeMetropolisRatio(m_particles, currentParticle, positionChange);

    // Check if trial state is accepted
    if (Random::nextDouble() <= qratio){
        m_waveFunction->updateSlaterDet(currentParticle);
        return true;
    }

    for (int i=0; i<m_numberOfDimensions; i++){
        // If trial state is not accepted, revert to old position for chosen particle (revert to old state)
        m_particles[currentParticle]->adjustPosition(-positionChange[i], i);
        m_waveFunction->updateDistances(currentParticle);
        m_waveFunction->updateSPWFMat(currentParticle);
        m_waveFunction->updateJastrow(currentParticle);
    }

    return false;
}

bool System::metropolisStepImpSampling(int currentParticle){

//    // Choose a random particle to change the position of
//    int randomParticle = Random::nextInt(m_numberOfParticles);
//    setRandomParticle(randomParticle);
    setCurrentParticle(currentParticle);

    std::vector<double> positionChange(m_numberOfDimensions);
    double D = 0.5;

    // Keep old position for Greens function
    std::vector<double> positionOld = m_particles[currentParticle/*randomParticle*/]->getPosition();

    // Change position of current particle
    for (int i=0; i < m_numberOfDimensions; i++){
        positionChange[i] = Random::nextGaussian(0., sqrt(m_dt)) + D*m_dt*quantumForce()[i];
    }

    double qratio = m_waveFunction->computeMetropolisRatio(m_particles, currentParticle/*randomParticle*/, positionChange);

    // Keep new position for Greens function
    std::vector<double> positionNew = m_particles[currentParticle/*randomParticle*/]->getPosition();

    // Evaluate Greens functions and find Metropolis-Hastings ratio:
    double GreensFunctionNew = calculateGreensFunction(positionNew, positionOld);

    for (int i=0; i < m_numberOfDimensions; i++){
        m_particles[currentParticle/*randomParticle*/]->adjustPosition(-positionChange[i], i);
    }

    double GreensFunctionOld = calculateGreensFunction(positionOld, positionNew);

    qratio *= GreensFunctionNew / GreensFunctionOld;

    // If move is accepted give the random particle the new position, otherwise keep the old position
    if (Random::nextDouble() <= qratio){
        for (int i=0; i<m_numberOfDimensions; i++){
            m_particles[currentParticle/*randomParticle*/]->adjustPosition(positionChange[i], i);
        }
        m_waveFunction->updateSlaterDet(currentParticle/*randomParticle*/);
        return true;
    }

    m_waveFunction->updateDistances(currentParticle/*randomParticle*/);
    m_waveFunction->updateSPWFMat(currentParticle/*randomParticle*/);
    m_waveFunction->updateJastrow(currentParticle/*randomParticle*/);

    return false;
}

std::vector<double> System::quantumForce(){
    // Calculate the quantum (drift) force. computeDerivative gives the gradient of the wave function

    std::vector<double> qForce;

    for (int i=0; i < m_numberOfDimensions; i++){
        qForce.push_back(
             2*m_waveFunction->computeDerivative(m_particles)[m_currentParticle*m_numberOfDimensions + i]);
    }
    return qForce;
}

double System::calculateGreensFunction(std::vector<double> positionOld, std::vector<double> positionNew){
    // Calculate the Greens function using the drift force from the quantumForce function

    double GreensFunction = 0;
    double D = 0.5;

    for (int i=0; i < m_numberOfDimensions; i++){
        double tmp = positionNew[i] - positionOld[i] - D*m_dt*quantumForce()[i];
        GreensFunction += tmp*tmp;
    }

    GreensFunction = exp(-GreensFunction/4*D*m_dt);
    return GreensFunction;
}

void System::runMetropolisSteps(int numberOfMetropolisSteps, bool importanceSampling,
                                bool showProgress, bool printToTerminal) {
    // Initialize Monte Carlo simulation
    //m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    // Find CPU time
    clock_t start, finish;
    start = clock();
    int percent = numberOfMetropolisSteps/100;
    int progress = 0;
    double equilibrationSteps = m_equilibrationFraction*numberOfMetropolisSteps;

    for (int i = 0; i < numberOfMetropolisSteps; i++) {
        // Update progress
        if (showProgress && m_my_rank == 0) {
            if (i%percent==0){
                progress += 1;
                cout << progress << "%" << "\n\033[F";
            }
        }
        for (int currentParticle = 0; currentParticle < m_numberOfParticles; currentParticle++){
            bool acceptedStep;
            // Ask Metropolis step functions whether the step was accepted or not
            if (importanceSampling)  acceptedStep = metropolisStepImpSampling(currentParticle);
            else                     acceptedStep = metropolisStep(currentParticle);

            /* Here you should sample the energy (and maybe other things using
             * the m_sampler instance of the Sampler class. Make sure, though,
             * to only begin sampling after you have let the system equilibrate
             * for a while. You may handle this using the fraction of steps which
             * are equilibration steps; m_equilibrationFraction.
             */

            if (i>=equilibrationSteps){
                // Sample energy etc.
                m_sampler->sample(acceptedStep);
            }
        }
    }

    finish = clock();
    m_computationTime = (finish-start)/ (double) CLOCKS_PER_SEC;
    m_printToTerminal = printToTerminal;//if (printToTerminal){;} //m_sampler->printOutputToTerminal();

    // Compute final expectation values
    m_sampler->computeAverages();
}

void System::optimizeParameters(System* system, double alpha, double beta) {
    // Steepest descent:
    int maxIterations             = 1000;
    int numberOfStepsSD           = (int) 1e3;
    double stepLengthSD           = 0.01;
    std::vector<double> initialParameters(2);
    initialParameters[0] = alpha;
    initialParameters[1] = beta;
    //double initialAlpha           = alpha;//0.7;
    double tol                    = 1e-5;//0.001;
    bool importanceSamplingSD     = false;
    //std::string parameterAlpha    = "alpha";

    if (m_my_rank == 0) {
        cout << "Optimizing parameters using steepest descent:" << endl;
    }
    SteepestDescent* steepestDescent = new SteepestDescent(system, stepLengthSD);
    steepestDescent->obtainOptimalParameter(initialParameters, tol, maxIterations,
                                            numberOfStepsSD, importanceSamplingSD);

//    maxIterations             = 100;
//    numberOfStepsSD           = (int) 1e5;
//    stepLengthSD              = 0.01;
//    double initialBeta        = beta;//1.142;//0.505;
//    tol                       = 1e-6;//0.001;
//    importanceSamplingSD      = false;
//    std::string parameterBeta = "beta";

//    if (m_my_rank == 0) {
//        cout << "Optimizing beta using steepest descent:" << endl;
//        SteepestDescent* steepestDescent = new SteepestDescent(system, stepLengthSD);
//        steepestDescent->obtainOptimalParameter(initialBeta, parameterBeta, tol, maxIterations,
//                                            numberOfStepsSD, importanceSamplingSD);
//    }
}

void System::MPI_CleanUp(double &totalE, double &totalKE, double &totalPE,
                         double &totalVariance, double &totalAcceptanceRate,
                         double &finalMeanDistance, double &timeStart, double &timeEnd,
                         double &totalTime, int numprocs, int numberOfSteps) {
    double e = m_sampler->getEnergy();
    double e_k = m_sampler->getKineticEnergy();
    double e_p = m_sampler->getPotentialEnergy();
    double variance = m_sampler->getVariance();
    double acceptanceRate = m_sampler->getAcceptanceRate();
    double mean_distance = m_sampler->getMeanDistance();

    MPI_Reduce(&e, &totalE, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&e_k, &totalKE, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&e_p, &totalPE, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&variance, &totalVariance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&acceptanceRate, &totalAcceptanceRate, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mean_distance, &finalMeanDistance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    timeEnd = MPI_Wtime();
    totalTime = timeEnd-timeStart;

    if (m_my_rank == 0){
        totalE /= numprocs;
        totalKE /= numprocs;
        totalPE /= numprocs;
        totalVariance /= numprocs;
        totalAcceptanceRate /= numprocs;
        finalMeanDistance /= numprocs;
        setNumberOfMetropolisSteps(numberOfSteps);
        setComputationTime(totalTime);
        m_sampler->setEnergy(totalE);
        m_sampler->setKineticEnergy(totalKE);
        m_sampler->setPotentialEnergy(totalPE);
        m_sampler->setVariance(totalVariance);
        m_sampler->setAcceptanceRate(totalAcceptanceRate);
        m_sampler->setMeanDistance(finalMeanDistance);
        if (m_printToTerminal) {m_sampler->printOutputToTerminal();}
    }
    if (m_saveEnergies) fclose(m_outfileE);
    if (m_savePositions) fclose(m_outfileP);
    MPI_Finalize();
}

void System::mergeOutputFiles(int numprocs) {
    if (m_saveEnergies){
        std::ofstream outfile("../vmc-solver/AnalyseData/Data/energies.dat", std::ios_base::binary);
        for (int i=0; i < numprocs; i++){
            char nodeFileName[100];
            sprintf(nodeFileName, "energiesNode%i.dat", i);
            std::ifstream nodeFile(nodeFileName, std::ios_base::binary);

            outfile << nodeFile.rdbuf();
            nodeFile.close();
            //remove(nodeFileName);
        }
        outfile.close();
    }

    if (m_savePositions){
        std::ofstream outfile("../vmc-solver/AnalyseData/Data/positions.dat", std::ios_base::binary);
        for (int i=0; i < numprocs; i++){
            char nodeFileName[100];
            sprintf(nodeFileName, "positionsNode%i.dat", i);
            std::ifstream nodeFile(nodeFileName, std::ios_base::binary);

            outfile << nodeFile.rdbuf();
            nodeFile.close();
            //remove(nodeFileName);
        }
        outfile.close();
    }

//    sleep(1);
//    for (int i=0; i < numprocs; i++){
//        char nodeFileNameE[100];
//        sprintf(nodeFileNameE, "energiesNode%i.dat", i);
//        char nodeFileNameP[100];
//        sprintf(nodeFileNameP, "positionsNode%i.dat", i);
//        remove(nodeFileNameE);
//        remove(nodeFileNameP);
//    }
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setTimeStep(double dt) {
    assert(dt > 0);
    m_dt = dt;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
    m_particles = initialState->getParticles();
}

void System::setNumberOfMetropolisSteps(int numberOfMetropolisSteps) {
    m_numberOfMetropolisSteps = numberOfMetropolisSteps;
}

void System::setCurrentParticle(int currentParticle) {
    m_currentParticle = currentParticle;
}

void System::setMyRank(int my_rank) {
    m_my_rank = my_rank;
}

void System::setNumProcs(int numprocs) {
    m_numprocs = numprocs;
}

void System::setComputationTime(double computationTime) {
    m_computationTime = computationTime;
}

void System::setSaveEnergies(bool saveEnergies) {
    m_saveEnergies = saveEnergies;
    if (saveEnergies){
        if (m_my_rank == 0){
            system("rm ../vmc-solver/AnalyseData/Data/energies.dat");
        }
        char outfileName[100];
        //char removeCommand[100];
        sprintf(outfileName, "energiesNode%i.dat", m_my_rank);
        remove(outfileName);
        //sprintf(removeCommand, "rm %s", outfileName);
        //system(removeCommand);
        m_outfileE = fopen(outfileName, "ab");
    }
}

void System::setSavePositions(bool savePositions) {
    m_savePositions = savePositions;
    if (savePositions){
        if (m_my_rank == 0){
            system("rm ../vmc-solver/AnalyseData/Data/positions.dat");
        }
        char outfileName[100];
        //char removeCommand[100];
        sprintf(outfileName, "positionsNode%i.dat", m_my_rank);
        remove(outfileName);
        //sprintf(removeCommand, "rm %s", outfileName);
        //system(removeCommand);
        m_outfileP = fopen(outfileName, "ab");
    }
}


void System::retrieveCoefficientsFromFile(string fileName, cube &loadCoefficients) {
    loadCoefficients.load(fileName, arma_ascii);
    return;
}

void System::retrieveConstantsFromFile(string fileName, vec &loadConstants) {
    loadConstants.load(fileName, raw_ascii);
    return;
}
