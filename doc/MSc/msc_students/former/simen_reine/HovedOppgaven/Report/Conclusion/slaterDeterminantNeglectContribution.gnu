set out 'Conclusion/slaterDeterminantNeglectContribution.tex'
set terminal pstex
set xlabel 'r'
Z=10
f(r) = exp(-Z*r)
g(r) = r*r*exp(-Z*r/3)
h(r) = exp(Z*(r-4))
x(r) = (4-r)*(4-r)*exp(Z*(r-4)/3)
set xlabel 'r'
set label 'A' at 0.1,0.1
set label 'B' at 0.7,0.1
set label 'C' at 3.3,0.1
set label 'D' at 3.9,0.1
set out '.trash'
plot [r=0:4] f(r) title '$\Psi_{1s}^1(r)$' 
replot g(r) title '$\Psi_{3d}^1(r)$' 
replot x(r) title '$\Psi_{3d}^2(r)$'
!rm -f .trash
set out 'Conclusion/slaterDeterminantNeglectContribution.tex'
replot h(r) title '$\Psi_{1s}^2(r)$'
set nolabel

