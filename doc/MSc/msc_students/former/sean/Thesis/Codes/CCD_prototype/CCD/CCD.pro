TEMPLATE = app
CONFIG  += console
CONFIG  -= app_bundle
CONFIG  -= qt
CONFIG(release, debug|release){CONFIG += optimize_full}
QMAKE_CFLAGS = -m64
QMAKE_LFLAGS = -m64 -fopenmp -lgfortran
QMAKE_LFLAGS -= -O1
QMAKE_LFLAGS += -O3
QMAKE_CXXFLAGS = -m64 -O3 -fopenmp -std=c++11

SOURCES += Systems/NLO2opt/constants.f90 \
    Systems/NLO2opt/spin_isospin_operators.f90 \
    Systems/NLO2opt/ang_mom_module.f90 \
    Systems/NLO2opt/chiral_nnlo_opt_simple.f90 \
    Systems/NLO2opt/chipot_f90_wrapper.f90

SOURCES += main.cpp \
    Systems/heg.cpp \
    Systems/mp.cpp \
    Systems/pm.cpp \
    master.cpp \
    Systems/system.cpp \
    makeampmat.cpp \
    makeintmat.cpp \
    diagrams.cpp \
    Systems/chiral.cpp

HEADERS += \
    Systems/heg.h \
    Systems/mp.h \
    Systems/pm.h \
    master.h \
    Systems/system.h \
    makeampmat.h \
    makeintmat.h \
    diagrams.h \
    Systems/chiral.h

#FORTRAN_SOURCES += \


#INCLUDEPATH += /usr/include/openmpi-x86_64
#QMAKE_CXX = mpicxx #/usr/lib64/openmpi/bin/mpicxx #mpicxx
#QMAKE_CXX += -O3
#QMAKE_CXX_RELEASE = $$QMAKE_CXX
#QMAKE_CXX_DEBUG = $$QMAKE_CXX
