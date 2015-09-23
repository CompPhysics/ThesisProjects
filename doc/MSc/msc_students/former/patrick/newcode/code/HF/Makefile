# shortcut for the compiler
#FOR PRODUCTION RUN:
#      - SET COMMENT and CHECK to zero in order to reduce computation time since they are usually placed in big FOR loops.
#      - remove the option -pg used with "gprof executable_name" to look at efficiency of the code
# level 1 of comment: DEBUGCONSTRUCTOR
# level 2 of comment: COMMENT
# level 3 of comment: CHECK
cc = g++ -pg -D DEBUGCONSTRUCTOR=1 -D COMMENT=1 -D CHECK=1 -D CHECK2=0
copt = -O2
lib = -I ../singleParticleOrbitals -I ../complex  -I ../CoulombMatrix/ -I ../blitz-0.9 -I ../ConfigFile -I ../singleParticleWaveFunctions -I ../singleParticleEnergies -I ../HF -I ../lib
#lib = -I lib -I singleParticleOrbitals -I complex


default: project.app #CVShtml

CVShtml: cvsHTML_script.py
	cvsHTML_script.py

project.app: main.o Complex.o singleParticleOrbitals.o CoulombMatrix.o ConfigFile.o singleParticleEnergies.o singleParticleWaveFunctions.o HFalgorithm.o lib.o
	$(cc) $(copt) -o project.app  main.o Complex.o singleParticleOrbitals.o CoulombMatrix.o ConfigFile.o singleParticleEnergies.o singleParticleWaveFunctions.o HFalgorithm.o lib.o

Complex.o: ../complex/Complex.h ../complex/Complex.cpp 
	$(cc) $(copt) $(lib) -c ../complex/Complex.cpp  

HFalgorithm.o: ../CoulombMatrix/CoulombMatrix.h ../singleParticleOrbitals/singleParticleOrbitals.h ../singleParticleWaveFunctions/singleParticleWaveFunctions.h ../singleParticleEnergies/singleParticleEnergies.h ../HF/HFalgorithm.h ../HF/HFalgorithm.cpp ../lib/lib.h
	$(cc) $(copt) -c $(lib) ../HF/HFalgorithm.cpp

CoulombMatrix.o: ../CoulombMatrix/CoulombMatrix.cpp ../CoulombMatrix/CoulombMatrix.h ../singleParticleOrbitals/singleParticleOrbitals.h
	$(cc) $(copt) -c $(lib) ../CoulombMatrix/CoulombMatrix.cpp

singleParticleEnergies.o: ../singleParticleOrbitals/singleParticleOrbitals.h ../singleParticleWaveFunctions/singleParticleWaveFunctions.h ../singleParticleEnergies/singleParticleEnergies.h ../singleParticleEnergies/singleParticleEnergies.cpp
	$(cc) $(copt) -c $(lib) ../singleParticleEnergies/singleParticleEnergies.cpp

singleParticleWaveFunctions.o: ../singleParticleOrbitals/singleParticleOrbitals.h ../singleParticleWaveFunctions/singleParticleWaveFunctions.h ../singleParticleWaveFunctions/singleParticleWaveFunctions.cpp
	$(cc) $(copt) -c $(lib) ../singleParticleWaveFunctions/singleParticleWaveFunctions.cpp

singleParticleOrbitals.o: ../singleParticleOrbitals/singleParticleOrbitals.cpp ../singleParticleOrbitals/singleParticleOrbitals.h
	$(cc) $(copt) -c $(lib) ../singleParticleOrbitals/singleParticleOrbitals.cpp

ConfigFile.o: ../ConfigFile/ConfigFile.h ../ConfigFile/ConfigFile.cpp
	$(cc) $(copt) -c $(lib) ../ConfigFile/ConfigFile.cpp

lib.o: ../lib/lib.h ../lib/lib.cpp
	$(cc) $(copt) -c $(lib) ../lib/lib.cpp

main.o: ../HF/main.cpp ../complex/Complex.h  ../singleParticleOrbitals/singleParticleOrbitals.h ../CoulombMatrix/CoulombMatrix.h ../ConfigFile/ConfigFile.h  ../singleParticleWaveFunctions/singleParticleWaveFunctions.h ../singleParticleEnergies/singleParticleEnergies.h ../lib/lib.h
	$(cc) $(copt) $(lib) -c  ../HF/main.cpp

clean: *.o *.app
	@ rm -f *.o *.app
