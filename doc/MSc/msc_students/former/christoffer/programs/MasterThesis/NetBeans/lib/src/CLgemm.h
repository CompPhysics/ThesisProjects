/* 
 * File:   CLgemm.h
 * Author: toffyrn
 *
 * Created on 15. februar 2012, 12:49
 */

#ifndef CLGEMM_H
#define	CLGEMM_H

#define __CL_ENABLE_EXCEPTIONS
#include "GEMM.h"
#include <CL/cl.hpp>
#include <clAmdBlas.h>

namespace toffyrn
{
    namespace libUNK
    {

        class CLgemm : public GEMM
        {
        public:
            CLgemm();
            CLgemm(const CLgemm& orig);
            virtual ~CLgemm();
            virtual void dgemm(
                    arma::mat &res,
                    arma::mat const &left,
                    arma::mat const &right,
                    double alpha = 1,
                    double beta = 0,
                    bool transL = false,
                    bool transR = false);
        private:
            cl::Context createSomeContext();
            cl::Context context;
            cl::CommandQueue queue;

        };
    }
}

#endif	/* CLGEMM_H */

