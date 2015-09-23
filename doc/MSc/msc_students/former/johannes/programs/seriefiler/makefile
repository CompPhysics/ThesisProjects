XLF= ifort -O3 -static  -C -u -g -traceback -warn all  #--tpp --staticlink
#lf95
twoholes.exe :   complex-modules.o twoholes.o complex-onebody3.o complex-interaction.o complex-library.o  complex-f77lib.o complex-gmatrix.o /mn/tid/kjerne-s1/johannre/fysikk/subatom/kjernefys/master/programmer/seriefiler/naglib_int/e01baf.o /mn/tid/kjerne-s1/johannre/fysikk/subatom/kjernefys/master/programmer/seriefiler/naglib_int/e02baf.o /mn/tid/kjerne-s1/johannre/fysikk/subatom/kjernefys/master/programmer/seriefiler/naglib_int/e02bbf.o /mn/tid/kjerne-s1/johannre/fysikk/subatom/kjernefys/master/programmer/seriefiler/naglib_int/p01aaf.o /mn/tid/kjerne-s1/johannre/fysikk/subatom/kjernefys/master/programmer/seriefiler/naglib_int/e01acf.o cmplx-bessel.o gsm-contour.o n3lo.o

	       ${XLF} -o  twoholes.exe  complex-modules.o twoholes.o complex-onebody3.o complex-interaction.o complex-library.o  complex-f77lib.o complex-gmatrix.o /mn/tid/kjerne-s1/johannre/fysikk/subatom/kjernefys/master/programmer/seriefiler/naglib_int/e01baf.o /mn/tid/kjerne-s1/johannre/fysikk/subatom/kjernefys/master/programmer/seriefiler/naglib_int/e02baf.o /mn/tid/kjerne-s1/johannre/fysikk/subatom/kjernefys/master/programmer/seriefiler/naglib_int/e02bbf.o /mn/tid/kjerne-s1/johannre/fysikk/subatom/kjernefys/master/programmer/seriefiler/naglib_int/p01aaf.o /mn/tid/kjerne-s1/johannre/fysikk/subatom/kjernefys/master/programmer/seriefiler/naglib_int/e01acf.o cmplx-bessel.o gsm-contour.o n3lo.o



complex-modules.o :   complex-modules.f90
	${XLF} -c  complex-modules.f90

twoholes.o :   twoholes.f90
	${XLF} -c twoholes.f90

complex-onebody3.o :   complex-onebody3.f90
	${XLF} -c complex-onebody3.f90

complex-interaction.o :   complex-interaction.f90
	${XLF} -c complex-interaction.f90

complex-library.o :   complex-library.f90
	${XLF} -c complex-library.f90 

complex-f77lib.o : complex-f77lib.f
	${XLF} -c  complex-f77lib.f

complex-gmatrix.o : complex-gmatrix.f90
	${XLF}  -c complex-gmatrix.f90

nag_interpol.o :  nag_interpol.f
	${XLF} -c  nag_interpol.f

e01acf.o :  e01acf.f
	${XLF} -c  /mn/tid/kjerne-s1/johannre/fysikk/subatom/kjernefys/master/programmer/seriefiler/naglib_int/e01acf.f

e01baf.o :  e01baf.f
	${XLF} -c  e01baf.f

e02baf.o :  e02baf.f  
	${XLF} -c e02baf.f  

e02bbf.o : e02bbf.f
	${XLF} -c e02bbf.f

p01aaf.o :  p01aaf.f 
	${XLF} -c p01aaf.f 

cmplx-bessel.o : cmplx-bessel.f
	${XLF} -c cmplx-bessel.f

gsm-contour.o : gsm-contour.f90
	${XLF} -c  gsm-contour.f90

n3lo.o : n3lo.f
	${XLF} -c n3lo.f 
