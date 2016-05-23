reset
set terminal lua fulldoc \
solid originreset plotsize 9,6.5
set output 'plotExecTimeMixedVsCppAlone.tex'

set format x "%g"
set format y "%g"

unset key
#set label 1 '' #at 5.5, 0.0193014, 0 left norotate back nopoint offset character 0, 0, 0
#set label 2 '$\sigma$' at 7.73607, 0.108213, 0 left norotate back nopoint offset character 0, 0, 0

set xlabel 'Monte Carlo cycles'  
set xrange [10000 : 110000] noreverse nowriteback
set mxtics 5

set ylabel 'Execution time, s'
set yrange [0 : 2.8] noreverse nowriteback
set mytics 5

set key top left

plot 'execTimeMixedHeSinThermOpt1.data' using 1:2 lt 1 lw 3 w lp t 'Mixed -O1', 'execTimeCppAloneOpt1.data' using 2:1 lt 2 lw 3 w lp t 'C++ -O1','execTimeMixedHeSinThermOpt2.data' using 1:2 lt 3 lw 3 w lp t 'Mixed -O2','execTimeCppAloneOpt2.data' using 2:1 lt 4 lw 5 w lp t 'C++ -O2','execTimeMixedHeSinThermOpt3.data' using 1:2 lt 5 lw 3 w lp t 'Mixed -O3','execTimeCppAloneOpt3.data' using 2:1 lt 6 lw 3 w lp t 'C++ -O3'












        
