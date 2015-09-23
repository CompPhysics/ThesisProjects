PROG=statistics
COMP=c++ -finline-functions -Wall #-pg #-g #-pg mpic++
DEBUG=#-D_DEBUG_
OPT=-O3

${PROG}:		${PROG}.o 
			${COMP} ${DEBUG} ${OPT} ${PROG}.o  -o statistics

${PROG}.o:		${PROG}.cpp 
			${COMP} ${DEBUG} ${OPT} -c ${PROG}.cpp

clear:
			rm -f *.o

clear~:
			rm -f *~

clearAll:	
		rm -f statistics
		make clear
		make clear~
		rm -f ./#* ./.#*
