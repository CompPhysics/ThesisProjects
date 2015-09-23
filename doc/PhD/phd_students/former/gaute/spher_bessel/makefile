# general makefile for F90 - choose PROG = name of given program 

# here we define compiler options, libraries and the target

F90 = lf95  -O3 --tpp --staticlink

# PROG = sphericalbessel  

# here we make the executable file 

sphericalbes.exe : beschb.o sphericalbessel.o bessjy.o chebev.o sphbes.o

	${F90} -o sphericalbes.exe beschb.o sphericalbessel.o bessjy.o sphbes.o chebev.o 

beschb.o : beschb.f90
	${F90} -c beschb.f90 

sphericalbessel.o : sphericalbessel.f90
		  ${F90} -c -free sphericalbessel.f90

bessjy.o : bessjy.f90
	 ${F90} -c bessjy.f90

chebev.o : chebev.f90
	 ${F90} -c chebev.f90

sphbes.o : sphbes.f90
	 ${F90} -c sphbes.f90