TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    quantumdots.cpp \
    quantumstate.cpp \
    Coulomb_Functions.cpp

LIBS += -llapack -lopenblas -larmadillo

HEADERS += \
    quantumdots.h \
    quantumstate.h \
    Coulomb_Functions.hpp

