#ifndef PROJECT2_INITIALSTATE_H
#define PROJECT2_INITIALSTATE_H
#include <vector>

class InitialState {
public:
    InitialState(class System* system);
    virtual void setupInitialState() = 0;
    std::vector<class Particle*> getParticles() { return m_particles; }

protected:
    class System* m_system = nullptr;
    std::vector<Particle*> m_particles;// = std::vector<Particle*>();
    int m_numberOfDimensions = 0;
    int m_numberOfParticles  = 0;
    int m_my_rank            = 0;
};

#endif // PROJECT2_INITIALSTATE_H
