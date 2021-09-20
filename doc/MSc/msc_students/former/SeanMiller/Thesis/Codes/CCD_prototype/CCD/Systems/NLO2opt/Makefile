CC= c++
FC= gfortran
FFlags=
CFLAGS= -std=c++11 -O3
LIBS= -lgfortran
EXECUTABLES = main.exe

O_FILES = chi_test.o chipot_cpp_wrapper.o chipot_f90_wrapper.o constants.o chiral_nnlo_opt_simple.o ang_mom_module.o spin_isospin_operators.o InfMatterSPBasis.cpp SPBasis.cpp util.cpp

$(EXECUTABLES): $(O_FILES)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm *.o *~ *# $(EXECUTABLES)
