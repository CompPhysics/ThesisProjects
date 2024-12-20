set terminal lua #fulldoc \
#solid originreset plotsize 9,6.5
set output 'plotExecTimeHePyCpp.tex'

set format x "%g"
set format y "%g"

unset key
#set label 1 '$Energy, au$' #at 5.5, 0.0193014, 0 left norotate back nopoint offset character 0, 0, 0
#set label 2 '$\sigma$' at 7.73607, 0.108213, 0 left norotate back nopoint offset character 0, 0, 0

set xlabel 'Monte Carlo cycles'  
set xrange [10000 : 110000] noreverse nowriteback
set mxtics 5

set ylabel 'Execution time, s'
set yrange [0 : 190] noreverse nowriteback
set mytics 5


plot 'execTimeHePython.data' using 1:2 lt 3 lw 5 w lp, \
'execTimeHeCppNew.data' using 1:2 lt 2 lw 5 w lp

