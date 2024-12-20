reset
set terminal lua #fulldoc \
#solid originreset plotsize 9,6.5
set output 'plotEnergiesHePlainCppdBeta003.tex'

set format x "%g"
set format y "%g"

unset key
#set label 1 '' #at 5.5, 0.0193014, 0 left norotate back nopoint offset character 0, 0, 0
#set label 2 '$\sigma$' at 7.73607, 0.108213, 0 left norotate back nopoint offset character 0, 0, 0

set xlabel 'Variantional parameter, $\beta$'  
#set xrange [0.01: 0.09] noreverse nowriteback
set mxtics 10

set ylabel 'Energy $\langle E \rangle$, au'
#set yrange [-2.85 : -2.81] noreverse nowriteback
set mytics 10

#set style line  100 linetype 1 pointtype  2 pointsize 1


set multiplot
plot 'resultsHe10MdBeta003.dat' using 1:2:4 lw 1 w yerrorbars, 'resultsHe10MdBeta003.dat' using 1:2 lt 4 lw 4 w l

set origin 0.4,0.4
set size 0.55,0.5
set xrange[0.35:0.9]
set xtics 0.4
set xtics font ",20"
unset xlabel
set yrange[-2.878:-2.87]
set ytics 0.01
unset ylabel
unset label
replot
set nomultiplot
