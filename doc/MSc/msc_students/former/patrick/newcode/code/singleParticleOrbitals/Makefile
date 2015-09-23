project.app: singleParticleOrbitals.o test1.o
	g++ -g singleParticleOrbitals.o test1.o -o project.app

singleParticleOrbitals.o: singleParticleOrbitals.h singleParticleOrbitals.cpp
	g++   -D DEBUGCONSTRUCTOR=1 -D COMMENT=1 -D CHECK=1 -c singleParticleOrbitals.cpp

test1.o: singleParticleOrbitals.h test1.cpp
	g++  -c test1.cpp