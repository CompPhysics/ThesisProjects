# shortcut for the compiler
#FOR PRODUCTION RUN:
#      - SET COMMENT and CHECK to zero in order to reduce computation time since they are usually placed in big FOR loops.
#      - remove the option -pg used with "gprof executable_name" to look at efficiency of the code
# level 1 of comment: DEBUGCONSTRUCTOR
# level 2 of comment: COMMENT
# level 3 of comment: CHECK
#cc = g++ -pg -D DEBUGCONSTRUCTOR=1 -D COMMENT=1 -D CHECK=1 -D CHECK2=0
cc = g++
copt = -O3
lib = -I ~/master/C++/BLITZ_test/blitz-0.9

default: example1.app

example1.app: example1.o
	$(cc) $(copt) -o example1.app  example1.o

example1.o: example1.cpp
	$(cc) $(copt) $(lib) -c  example1.cpp

clean:
	rm -f *.o *.app
