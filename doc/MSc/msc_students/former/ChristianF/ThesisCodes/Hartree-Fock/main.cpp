#include <iostream>
#include "quantumdots.h"

using namespace std;

int main()
{
    int numberOfShells = 2;
    int numberOfParticles = 2;
    int numberOfDimensions = 2;     //Only 2D for now.
    double omega = 1.;

    QuantumDots quantumDots(numberOfShells, numberOfParticles, numberOfDimensions, omega);
    quantumDots.runHartreeFock();
}

