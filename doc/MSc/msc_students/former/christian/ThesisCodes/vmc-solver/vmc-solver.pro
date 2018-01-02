TEMPLATE = app
CONFIG  += console c++11
CONFIG  -= app_bundle
CONFIG  -= qt
QMAKE_CXXFLAGS_RELEASE -= -O
QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2

QMAKE_CXXFLAGS_RELEASE += -O3

QMAKE_CXXFLAGS_RELEASE += -ffast-math

SOURCES += main.cpp \
    system.cpp \
    Hamiltonians/hamiltonian.cpp \
    Hamiltonians/harmonicoscillator.cpp \
    particle.cpp \
    WaveFunctions/wavefunction.cpp \
    InitialStates/initialstate.cpp \
    InitialStates/randomuniform.cpp \
    Math/random.cpp \
    Math/factorial.cpp \
    sampler.cpp \
    WaveFunctions/simplegaussian.cpp \
    Hamiltonians/harmonicoscillatorrepulsive.cpp \
    WaveFunctions/repulsivegaussian.cpp \
    VariationMethods/steepestdescent.cpp \
    VariationMethods/conjugategradient.cpp \
    Hamiltonians/harmonicoscillatorelectrons.cpp \
    WaveFunctions/twoelectrons.cpp \
    WaveFunctions/manyelectrons.cpp \
    BasisFunctions/basisfunctions.cpp \
    WaveFunctions/manyelectrons_coefficients.cpp \
    Hamiltonians/squarewell.cpp \
    Hamiltonians/finiteharmonicoscillator.cpp \
    Hamiltonians/doubleharmonicoscillator.cpp \
    HermitePolynomials/hermitepolynomials.cpp
    DMC/walker.cpp

LIBS += -llapack -lblas -larmadillo -lUnitTest++

HEADERS += \
    system.h \
    Hamiltonians/hamiltonian.h \
    Hamiltonians/harmonicoscillator.h \
    particle.h \
    WaveFunctions/wavefunction.h \
    InitialStates/initialstate.h \
    InitialStates/randomuniform.h \
    Math/random.h \
    Math/factorial.h \
    sampler.h \
    WaveFunctions/simplegaussian.h \
    Hamiltonians/harmonicoscillatorrepulsive.h \
    WaveFunctions/repulsivegaussian.h \
    VariationMethods/steepestdescent.h \
    VariationMethods/conjugategradient.h \
    Hamiltonians/harmonicoscillatorelectrons.h \
    WaveFunctions/twoelectrons.h \
    WaveFunctions/manyelectrons.h \
    BasisFunctions/basisfunctions.h \
    WaveFunctions/manyelectrons_coefficients.h \
    Hamiltonians/squarewell.h \
    Hamiltonians/finiteharmonicoscillator.h \
    Hamiltonians/doubleharmonicoscillator.h \
    HermitePolynomials/hermitepolynomials.h
    DMC/walker.h

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
