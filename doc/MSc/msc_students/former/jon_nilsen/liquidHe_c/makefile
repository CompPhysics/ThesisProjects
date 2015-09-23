# General makefile for c - choose PROG =   name of given program

# Here we define compiler option, libraries and the  target
CC= gcc -O3 -pg -Wall
PROG= liquidHelium
# this is the math library in C
LIB = -lm lib.o
# Here we make the executable file 
${PROG} :	${PROG}.o
		${CC} ${PROG}.o ${LIB} -o ${PROG}

# whereas here we create the object file

${PROG}.o :	${PROG}.c
		${CC} -c ${PROG}.c

