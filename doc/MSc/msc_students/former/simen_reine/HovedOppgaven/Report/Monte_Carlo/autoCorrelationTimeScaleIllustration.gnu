set out 'Monte_Carlo/autoCorrelationTimeScaleIllustration.tex'
set terminal pstex
f(x) = exp(-x) + 0.2*exp(-3*(x-5)*(x-5))+exp(-3*(x-7.3)*(x-7.3))
set xlabel 'x'
set label 'A' at 0.9,0.1
set label 'B' at 7.15,0.1
plot [x=0:10] f(x) title '$\Psi(x)$'
set nolabel
