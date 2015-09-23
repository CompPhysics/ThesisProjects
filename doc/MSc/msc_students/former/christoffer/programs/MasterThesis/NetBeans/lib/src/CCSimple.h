/* 
 * File:   CCSimple.h
 * Author: toffyrn
 *
 * Created on 27. desember 2011, 09:10
 */

#ifndef CCSIMPLE_H
#define	CCSIMPLE_H

#include "System.h"
#include <armadillo>
#include <iomanip>

namespace toffyrn
{
    namespace libUNK
    {

        class CCSimple
        {
        public:
            static double testSolver(System * sys, int maxiter, int debugFlags = DEBUG_NONE);
            static const int DEBUG_NONE = 0;
            static const int DEBUG_T1 = 1 << 1;
            static const int DEBUG_T2 = 1 << 2;
            static const int DEBUG_T_ALL = DEBUG_T1 | DEBUG_T2;
            static const int DEBUG_I1 = 1 << 3;
            static const int DEBUG_I2 = 1 << 4;
            static const int DEBUG_I3 = 1 << 5;
            static const int DEBUG_I4 = 1 << 6;
            static const int DEBUG_I5 = 1 << 7;
            static const int DEBUG_I6 = 1 << 8;
            static const int DEBUG_I7 = 1 << 9;
            static const int DEBUG_I8 = 1 << 10;
            static const int DEBUG_I9 = 1 << 11;
            static const int DEBUG_I10 = 1 << 12;
            static const int DEBUG_I11 = 1 << 13;
            static const int DEBUG_I_T1 = DEBUG_I1 | DEBUG_I2 | DEBUG_I3 | DEBUG_I4 | DEBUG_I5;
            static const int DEBUG_I_T2 = DEBUG_I6 | DEBUG_I7 | DEBUG_I8 | DEBUG_I9 | DEBUG_I10 | DEBUG_I11;
            static const int DEBUG_I_ALL = DEBUG_I_T1 | DEBUG_I_T2;
            static const int DEBUG_ENERGY = 1 << 14;
            static const int DEBUG_INTERACTIONS = 1 << 15;
            static const int DEBUG_D1 = 1 << 16;
            static const int DEBUG_D2 = 1 << 17;
            static const int DEBUG_ALL = DEBUG_T_ALL | DEBUG_I_ALL | DEBUG_ENERGY | DEBUG_INTERACTIONS | DEBUG_D1 | DEBUG_D2;
        private:

        };

    }
}

#endif	/* CCSIMPLE_H */

