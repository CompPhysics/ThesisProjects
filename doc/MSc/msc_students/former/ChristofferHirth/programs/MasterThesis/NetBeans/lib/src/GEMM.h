/* 
 * File:   GEMM.h
 * Author: toffyrn
 *
 * Created on 24. januar 2012, 13:43
 */

#ifndef GEMM_H
#define	GEMM_H

#include <armadillo>

namespace toffyrn
{
    namespace libUNK
    {

        class GEMM
        {
        public:
            GEMM();
            GEMM(const GEMM& orig);
            virtual ~GEMM();

            /**
             * @brief Do matrix multiplication: 
             * res = alpha * left * right + beta * res
             * 
             * @param res Result, will be changed
             * @param left Left matrix
             * @param right Right matrix
             * @param alpha Scale the matrix product
             * @param beta Add beta * res to the result.
             * @param transL Transpose left matrix?
             * @param transR Transpose right matrix?
             */
            virtual void dgemm(
                    arma::mat &res,
                    arma::mat const &left,
                    arma::mat const &right,
                    double alpha = 1,
                    double beta = 0,
                    bool transL = false,
                    bool transR = false);

            virtual void operator()(
                    arma::mat &res,
                    arma::mat const &left,
                    arma::mat const &right,
                    double alpha = 1,
                    double beta = 0,
                    bool transL = false,
                    bool transR = false)
            {
                dgemm(res, left, right, alpha, beta, transL, transR);
            }

            virtual arma::mat dgemm(
                    arma::mat const &left,
                    arma::mat const &right,
                    double alpha = 1.0,
                    bool transL = false,
                    bool transR = false)
            {
                int rows = left.n_rows;
                int cols = right.n_cols;
                arma::mat res = arma::zeros<arma::mat > (rows, cols);

                double beta = 0.0; //result matrix is undefined before multiplication.
                dgemm(res, left, right, alpha, beta, transL, transR);
                return res;
            }

            /**
             * @brief Total time spent doing matrix multiplication.
             * @return GEMM::tot_time
             */
            virtual double get_tot_time() const
            {
                return tot_time;
            }
        protected:
            /** Timer to see how many seconds used on matmult */
            arma::wall_clock timer;
            /** Seconds used accumulates into this variable */
            double tot_time;

        private:

        };
    }
}

#endif	/* GEMM_H */

