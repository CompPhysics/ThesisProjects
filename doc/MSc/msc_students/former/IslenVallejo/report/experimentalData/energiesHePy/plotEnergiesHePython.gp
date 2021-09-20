set terminal lua #fulldoc \
solid originreset #plotsize 9,6.5
set output 'plotEnergiesHePython.tex'

set format x "%g"
set format y "%g"

unset key
#set label 1 '' #at 5.5, 0.0193014, 0 left norotate back nopoint offset character 0, 0, 0
#set label 2 '$\sigma$' at 7.73607, 0.108213, 0 left norotate back nopoint offset character 0, 0, 0

set xlabel 'Variantional parameter, $\beta$'  
set xrange [0.1: 0.6] noreverse nowriteback
set mxtics 10

set mxtics 10
set ylabel 'Energy $\langle E \rangle$, au'
set yrange [-2.88 : -2.83] noreverse nowriteback
set mytics 5


#plot 'energiesHePy.data' using 1:2 lt 2 lw 3 w lp


set multiplot
plot 'energiesHePy.data' using 1:2:4 lw 1 w yerrorbars, 'energiesHePy.data' using 1:2 lt 4 lw 4 w l

set origin 0.4,0.4
set size 0.55,0.5
set xrange[0.4:0.6]
set xtics 0.4
set xtics font ",20"
unset xlabel
set yrange[-2.878:-2.87]
set ytics 0.01
unset ylabel
unset label
replot
set nomultiplot

