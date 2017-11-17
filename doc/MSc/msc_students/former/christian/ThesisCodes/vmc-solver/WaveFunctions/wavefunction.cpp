#include "wavefunction.h"


WaveFunction::WaveFunction(System* system) {
    m_system = system;
}

void WaveFunction::adjustParameter(double parameter, int parameterNumber){
    // Adjust variational parameter for steepest descent
    m_parameters[parameterNumber] = parameter;
}
