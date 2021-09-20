set terminal lua fulldoc \
solid originreset plotsize 9,6.5
set output 'plotEnergiesDt001.tex'

set format x "%g"
set format y "%g"

unset key
#set label 1 '' #at 5.5, 0.0193014, 0 left norotate back nopoint offset character 0, 0, 0
#set label 2 '$\sigma$' at 7.73607, 0.108213, 0 left norotate back nopoint offset character 0, 0, 0

set xlabel 'Variantional parameter, $\beta$'  
set xrange [0.01: 0.09] noreverse nowriteback
set mxtics 5

set ylabel 'Energy $\langle E \rangle$, au'
set yrange [-2.85 : -2.81] noreverse nowriteback
set mytics 5


#plot 'dt001.data' using 2:3 lt 2 lw 3 w lp
plot 'dt001.data' using 2:3:5 lt 2 lw 2 w yerrorbars
