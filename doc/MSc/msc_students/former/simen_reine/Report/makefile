#hovedfag.ps: hovedfag.dvi 
#	dvips -o hovedfag.ps hovedfag.dvi
#	gv hovedfag.ps&

hovedfag.ps: hovedfag.tex Appendix/appendix.tex Atom_Physics/atom_physics.tex Basic_QM/basic_qm.tex Conclusion/conclusion.tex Implementation/implementation.tex Introduction/introduction.tex Monte_Carlo/monte_carlo.tex Preface/preface.tex Results/results.tex References/references.bib
	latex hovedfag.tex
	bibtex hovedfag
	latex hovedfag.tex
	latex hovedfag.tex
	dvips -o hovedfag.ps hovedfag.dvi


clear:		
		rm -f *.dvi  *.bbl *.blg *.log *.toc *.aux
		rm -f Appendix/appendix.aux
		rm -f Atom_Physics/atom_physics.aux
		rm -f Basic_QM/basic_qm.aux
		rm -f Conclusion/conclusion.aux
		rm -f Implementation/implementation.aux
		rm -f Introduction/introduction.aux 
		rm -f Monte_Carlo/monte_carlo.aux
		rm -f Preface/preface.aux 
		rm -f Results/results.aux
		touch hovedfag.tex

clear~:		
		rm -f *~
		rm -f Appendix/*~
		rm -f Atom_Physics/*~
		rm -f Basic_QM/*~
		rm -f Conclusion/*~
		rm -f Implementation/*~
		rm -f Introduction/*~
		rm -f Monte_Carlo/*~
		rm -f Preface/*~
		rm -f Results/*~

clearAll:		
		make clear
		make clear~
		rm -f *.ps
		rm -f ./#* ./.#*
		rm -f -r .xvpics
