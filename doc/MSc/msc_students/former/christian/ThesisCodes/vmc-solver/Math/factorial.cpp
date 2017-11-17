#include "factorial.h"
#include <cassert>
#include <stdio.h>

double factorial(int n) {
    assert(n>=0);
    double fact = 1.;
    for (int i = 2; i <= n; i++) {
        fact *= i;
    }
    return fact;
}

