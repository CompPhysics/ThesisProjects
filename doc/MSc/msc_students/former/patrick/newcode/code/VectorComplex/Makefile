default: project.app

project.app: test2.o Vector_Complex.o Complex.o
	g++ -o project.app test2.o Vector_Complex.o Complex.o

Vector_Complex.o: ../dir1/Complex.h Vector_Complex.h Vector_Complex.cpp 
	g++ -I ../dir1 -c Vector_Complex.cpp

Complex.o: ../dir1/Complex.h ../dir1/Complex.cpp 
	g++ -I ../dir1 -c ../dir1/Complex.cpp  

test2.o: test2.cpp ../dir1/Complex.h Vector_Complex.h 
	g++ -I ../dir1 -c test2.cpp 

clean:
	rm -f *.o
