TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS_RELEASE += -O3

QMAKE_CXXFLAGS_RELEASE += -ffast-math

SOURCES += main.cpp \
    system.cpp \
    WaveFunctions/wavefunction.cpp \
    initialize_electrons.cpp \
    WaveFunctions/doublewell.cpp \
    Math/factorial.cpp \
    WaveFunctions/finitewell.cpp \
    WaveFunctions/squarewell.cpp

LIBS += -llapack -lopenblas -larmadillo

HEADERS += system.h \
    WaveFunctions/wavefunction.h \
    WaveFunctions/doublewell.h \
    Math/factorial.h \
    WaveFunctions/finitewell.h \
    WaveFunctions/squarewell.h
