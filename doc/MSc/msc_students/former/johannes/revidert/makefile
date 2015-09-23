# makefile for pdflatex

#epstopdf *.eps
PDFL=pdflatex
BIB=bibtex
infinitematter.pdf:	infinitematter.tex frontpage.tex notations.tex introduction.tex chapter2.tex chapter3.tex chapter4.tex chapter5.tex chapter6.tex chapter7.tex chapter8.tex chapter9.tex appendixA.tex appendixB.tex appendixC.tex appendixD.tex
	${PDFL} infinitematter && ${PDFL} infinitematter && ${PDFL} infinitematter
#infinitematter.bbl:	infinitematter.tex
#	${BIB} infinitematter #infinitematter.pdf:	infinitematter.tex
#	${PDFL} infinitematter  ${BIB} infinitematter && 
#	${PDFL} infinitematter
#PDF: infinitematter.tex
#	kpdf infinitematter.pdf

make clean:
	rm *.eps *.aux *.log
