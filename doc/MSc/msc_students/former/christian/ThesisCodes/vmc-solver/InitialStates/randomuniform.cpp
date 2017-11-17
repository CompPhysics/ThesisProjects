#include "randomuniform.h"
#include <iostream>
#include <cassert>
#include "../Math/random.h"
#include "../particle.h"
#include "../system.h"

using namespace std;

RandomUniform::RandomUniform(System*    system,
                             int        numberOfDimensions,
                             int        numberOfParticles, int my_rank)  :
        InitialState(system) {
    assert(numberOfDimensions > 0 && numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    m_numberOfParticles  = numberOfParticles;
    m_my_rank            = my_rank;

    /* The Initial State class is in charge of everything to do with the
     * initialization of the system; this includes determining the number of
     * particles and the number of dimensions used. To make sure everything
     * works as intended, this information is passed to the system here.
     */
    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
    setupInitialState();
}

void RandomUniform::setupInitialState() {

    long idum = -1-m_my_rank;
    Random::setSeed(idum);

    // Create random positions for all particles:
    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();

        for (int j=0; j < m_numberOfDimensions; j++) {
            /* This is where you should actually place the particles in
             * some positions, according to some rule. Since this class is
             * called random uniform, they should be placed randomly according
             * to a uniform distribution here. However, later you will write
             * more sub-classes of the InitialState class in which the
             * particles are placed in other configurations.
             *
             * Note: For now, the particles are simply placed in positions
             * according to their index in the particles list (this is
             * obviously NOT a good idea).
             */

            double pos = pow(-1, i%2)*Random::nextDouble();

            // Ensure that starting position of particles are reasonably close to the center of the well when a double well is used:
            vec L = m_system->getL();
            bool doubleWellFlag = m_system->getDoubleWellFlag();
            if (L[j] != 0 && doubleWellFlag) pos *= L[j];
            position.push_back(pos);


            //position.push_back(i);
        }
        // Create particles:
        m_particles.push_back(new Particle());
        m_particles.at(i)->setNumberOfDimensions(m_numberOfDimensions);
        m_particles.at(i)->setPosition(position);
    }
}
