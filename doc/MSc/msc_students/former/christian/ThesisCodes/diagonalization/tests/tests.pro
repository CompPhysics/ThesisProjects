TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    ../system.cpp \
    ../WaveFunctions/wavefunction.cpp \
    ../initialize_electrons.cpp \
    ../WaveFunctions/doublewell.cpp \
    ../Math/factorial.cpp \
    wrapper.cpp

LIBS += -llapack -lblas -larmadillo -lUnitTest++

HEADERS += ../system.h \
    ../WaveFunctions/wavefunction.h \
    ../WaveFunctions/doublewell.h \
    ../Math/factorial.h \
    wrapper.h
