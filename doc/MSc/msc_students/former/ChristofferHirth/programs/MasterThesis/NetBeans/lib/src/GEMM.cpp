/* 
 * File:   GEMM.cpp
 * Author: toffyrn
 * 
 * Created on 24. januar 2012, 13:43
 */

#include "GEMM.h"
#include <iostream>

using namespace toffyrn::libUNK;

GEMM::GEMM()
{
    tot_time = 0;
}

GEMM::GEMM(const GEMM& orig)
{
}

GEMM::~GEMM()
{
}

void GEMM::dgemm(
        arma::mat &res,
        arma::mat const &left,
        arma::mat const &right,
        double alpha,
        double beta,
        bool transL,
        bool transR)
{
    timer.tic();
    if (transL == false && transR == false)
        res = alpha * left * right + beta * res;
    else if (transL == true && transR == false)
        res = alpha * trans(left) * right + beta * res;
    else if (transL == false && transR == true)
        res = alpha * left * trans(right) + beta * res;
    else if (transL == true && transR == true)
        res = alpha * trans(left) * trans(right) + beta * res;
    tot_time += timer.toc();

    return;
}

