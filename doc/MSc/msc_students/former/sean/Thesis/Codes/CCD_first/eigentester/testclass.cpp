#include "testclass.h"

#include <omp.h>

using namespace std;

testClass::testClass()
{
}

void testClass::testFunc(){

    for (int i=0; i<1e4; i++){
        places[i] = (double) i;
    }
}


int testClass::testy( int id ){
    int length = D1.size();
    int returnVal;
    int n = omp_get_max_threads();

    #pragma for num_threads(n)
    for (int i=0; i<length; i++){
        if (id == D1[i]){ returnVal = D2[i]; }
    }

    return returnVal;
}
