# shortcut for the compiler
#FOR PRODUCTION RUN:
#      - SET COMMENT and CHECK to zero in order to reduce computation time since they are usually placed in big FOR loops.
#      - remove the option -pg used with "gprof executable_name" to look at efficiency of the code
# level 1 of comment: DEBUGCONSTRUCTOR
# level 2 of comment: COMMENT
# level 3 of comment: CHECK
cc = g++ -D DEBUGCONSTRUCTOR=1 -D COMMENT=1 -D CHECK=1 -D CHECK2=0
copt = -O2
lib = -I ../singleParticleOrbitals -I ../complex  -I ../ConfigFile -I ../blitz-0.9


default: project.app

project.app: test01.o Complex.o singleParticleOrbitals.o ConfigFile.o singleParticleWaveFunctions.o
	$(cc) $(copt) -o project.app  test01.o Complex.o singleParticleOrbitals.o ConfigFile.o singleParticleWaveFunctions.o


Complex.o: ../complex/Complex.h ../complex/Complex.cpp 
	$(cc) $(copt) $(lib) -c ../complex/Complex.cpp

singleParticleWaveFunctions.o: ../singleParticleOrbitals/singleParticleOrbitals.h singleParticleWaveFunctions.h singleParticleWaveFunctions.cpp
	$(cc) $(copt) -c $(lib) ../singleParticleWaveFunctions/singleParticleWaveFunctions.cpp

singleParticleOrbitals.o: ../singleParticleOrbitals/singleParticleOrbitals.cpp ../singleParticleOrbitals/singleParticleOrbitals.h
	$(cc) $(copt) -c $(lib) ../singleParticleOrbitals/singleParticleOrbitals.cpp

ConfigFile.o: ../ConfigFile/ConfigFile.h ../ConfigFile/ConfigFile.cpp
	$(cc) $(copt) -c $(lib) ../ConfigFile/ConfigFile.cpp

test01.o: test01.cpp ../complex/Complex.h  ../singleParticleOrbitals/singleParticleOrbitals.h ../CoulombMatrix/CoulombMatrix.h ../ConfigFile/ConfigFile.h 
	$(cc) $(copt) $(lib) -c  test01.cpp

clean: *.o *.app
	@ rm -f *.o *.app
