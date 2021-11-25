XLF= lf95 -O3 --tpp --staticlink


prog_vmche4.exe :   vmche4.o f90lib.o 
	       ${XLF} -o  prog_vmche4.exe  vmche4.o f90lib.o

f90lib.o :   f90lib.f90
		${XLF} -c  f90lib.f90


vmche4.o :   vmche4.f90
		${XLF} -c vmche4.f90
