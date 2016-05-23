set terminal lua #fulldoc \
#solid originreset plotsize 9,6.5
set output 'plotDiffTimePyCpp.tex'

set format x "%g"
set format y "%g"

unset key
set label 1 '$t = 18\times 10^{-4}mc - 0.59$' at 60000, 80, 0 left norotate back nopoint offset character 0, 0, 0
#set label 2 '$\sigma$' at 7.73607, 0.108213, 0 left norotate back nopoint offset character 0, 0, 0


set xlabel 'Monte Carlo cycles'  
set xrange [10000 : 110000] noreverse nowriteback
set mxtics 5

set ylabel 'time$_{Py}$ - time$_{C++}$'
set yrange [0 : 180] noreverse nowriteback
set mytics 5



plot 'relativeTimePyCpp.data'using 1:4 lt 3 lw 5 w lp

