project.app: Complex.o test1.o
	g++ Complex.o test1.o -o project.app

Complex.o: Complex.h Complex.cpp
	g++  -c Complex.cpp

test1.o: Complex.h test1.cpp
	g++  -c test1.cpp