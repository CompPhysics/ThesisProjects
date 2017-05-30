#ifndef chipot_cpp_wrapper_hpp_
#define chipot_cpp_wrapper_hpp_
#include <complex>
#include "SPBasis.hpp"

std::complex<double> chipot_cpp_wrapper(SPBasis *basis, double density, int p, int q, int r, int s);

void chipot_regression_test();

#endif
