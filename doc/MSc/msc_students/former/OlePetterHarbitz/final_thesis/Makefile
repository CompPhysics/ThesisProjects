TEXFILES=refs.bib thesis.tex chapter1.tex chapter2.tex chapter3.tex chapter4.tex chapter5.tex 1*.tex orb_occ_1*.tex *1*_0p3-0d3.eps goldstone.eps nn-potential.eps goldstoneLD.eps goldstoneTD.eps C14E2.eps C16E2.eps qbox1.eps qbox2.eps

dvi: thesis.dvi
ps: thesis.ps
pdf: thesis.pdf
zip: thesis.tgz
thesis.dvi: $(TEXFILES)
	latex thesis.tex && bibtex thesis && latex thesis.tex && latex thesis.tex
thesis.ps: thesis.dvi
	dvips thesis.dvi
thesis.pdf: thesis.ps
	ps2pdf thesis.ps
clean:
	rm thesis.{aux,bbl,blg,log,dvi,pdf,ps,toc}
zip: tar vfzc thesis.tgz refs.bib thesis.tex chapter1.tex chapter2.tex chapter3.tex chapter4.tex chapter5.tex 1*.tex orb_occ_1*.tex *1*_0p3-0d3.eps goldstone.eps nn-potential.eps goldstoneLD.eps goldstoneTD.eps C14E2.eps C16E2.eps
