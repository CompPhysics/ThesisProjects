TEMPLATE = app
CONFIG  += console c++11
CONFIG  -= app_bundle
CONFIG  -= qt
QMAKE_CXXFLAGS += -std=c++11

SOURCES += main.cpp \
    testclass.cpp

HEADERS += \
    testclass.h

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp

#INCLUDEPATH += /usr/local/include/openmpi
#INCLUDEPATH += /usr/local/lib
#QMAKE_CXX = /usr/local/bin/mpicxx #/usr/lib64/openmpi/bin/mpicxx #mpicxx
#QMAKE_CXX_RELEASE = $$QMAKE_CXX
#QMAKE_CXX_DEBUG = $$QMAKE_CXX
#QMAKE_LINK = $$QMAKE_CXX
#QMAKE_CC = /usr/local/bin/mpicc #/usr/lib64/openmpi/bin/mpicc #mpicc

#QMAKE_CFLAGS += $$system(mpicc --showme:compile)
#QMAKE_LFLAGS += $$system(mpicxx --showme:link)
#QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
#QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
