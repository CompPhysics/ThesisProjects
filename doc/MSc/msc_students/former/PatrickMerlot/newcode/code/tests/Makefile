# shortcut for the compiler
cc = g++
copt = -O3

default: project.app  SolvingEigValProblem.app SolvingEigValProblem2.app

project.app: tests.o State_Vector.o Complex.o All_States.o lib.o
	$(cc) $(copt) -o project.app tests.o State_Vector.o Complex.o lib.o

SolvingEigValProblem2.app: lib.o SolvingEigValProblem2.o
	$(cc) $(copt) -o SolvingEigValProblem2.app SolvingEigValProblem2.o lib.o

SolvingEigValProblem2.o: SolvingEigValProblem2.cpp
	$(cc) $(copt) -c -I ../lib/ SolvingEigValProblem2.cpp

SolvingEigValProblem.app: lib.o SolvingEigValProblem.o
	$(cc) $(copt) -o SolvingEigValProblem.app SolvingEigValProblem.o lib.o

SolvingEigValProblem.o: SolvingEigValProblem.cpp
	$(cc) $(copt) -c -I ../lib SolvingEigValProblem.cpp

All_States.o: ../complex/Complex.h ../StateVector/State_Vector.h ../AllStates/All_States.h
	$(cc) $(copt) -I ../complex -I ../StateVector -I ../AllStates -c ../AllStates/All_States.cpp

State_Vector.o: ../complex/Complex.h ../StateVector/State_Vector.h ../StateVector/State_Vector.cpp 
	$(cc) $(copt) -I ../complex -c ../StateVector/State_Vector.cpp

Complex.o: ../complex/Complex.h ../complex/Complex.cpp 
	$(cc) $(copt) -I ../complex -c ../complex/Complex.cpp  

lib.o: ../lib/lib.cpp ../lib/lib.h
	$(cc) $(copt) -c -I ../lib/ ../lib/lib.cpp

tests.o: tests.cpp ../complex/Complex.h ../StateVector/State_Vector.h ../AllStates/All_States.h ../lib/lib.h
	$(cc) $(copt) -c -I ../lib/ -I ../complex -I ../StateVector -I ../AllStates -c tests.cpp

clean: *.o *.app out.m
	@ rm -f *.o *.app out.m