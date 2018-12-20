# shortcut for the compiler
#FOR PRODUCTION RUN:
#      - SET COMMENT and CHECK to zero in order to reduce computation time since they are usually placed in big FOR loops.
#      - remove the option -pg used with "gprof executable_name" to look at efficiency of the code
cc = g++ -D DEBUGCONSTRUCTOR=1 -D COMMENT=1 -D CHECK=1
copt = -O3
lib = -I ../singleParticleOrbitals -I ../CoulombMatrix
#lib = -I ../singleParticleOrbitals -I ../CoulombMatrix -I ../lib

default: project.app

project.app: test1.o CoulombMatrix.o singleParticleOrbitals.o
	$(cc) $(copt) -o project.app -pg test1.o CoulombMatrix.o singleParticleOrbitals.o

#lib.o: ../lib/lib.cpp ../lib/lib.h
#	$(cc) $(copt) $(lib)-I ../lib -c ../lib/lib.cpp

CoulombMatrix.o: CoulombMatrix.cpp  CoulombMatrix.h ../singleParticleOrbitals/singleParticleOrbitals.h
	$(cc) $(copt)  $(lib) -c CoulombMatrix.cpp  

singleParticleOrbitals.o: ../singleParticleOrbitals/singleParticleOrbitals.cpp ../singleParticleOrbitals/singleParticleOrbitals.h
	$(cc) $(copt) $(lib) -c ../singleParticleOrbitals/singleParticleOrbitals.cpp

test1.o: test1.cpp ../CoulombMatrix/CoulombMatrix.h ../singleParticleOrbitals/singleParticleOrbitals.h
	#$(cc) $(copt) -c -I ../lib -I ../CoulombMatrix -I ../singleParticleOrbitals test1.cpp
	$(cc) $(copt) -c $(lib) test1.cpp

clean:
	rm -f *.o
