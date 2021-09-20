/* 
 * File:   HF.h
 * Author: toffyrn
 *
 * Created on 31. januar 2012, 09:57
 */

#ifndef HF_H
#define	HF_H

#include "System.h"

namespace toffyrn
{
    namespace libUNK
    {

        class HF
        {
        public:
            HF(double precision = 0.001, int max_iter = 20);
            virtual ~HF();
            void set_system(System * sys);
            std::pair<double, arma::mat> solve_ground_state_energy(int debug_flags = DEBUG_NONE, std::ostream &debug_stream = std::cout);


            ////FLAGS FOR DEBUGGING
            static const int DEBUG_NONE = 0;
            static const int DEBUG_ENERGY = 1 << 14;
            static const int DEBUG_INTERACTIONS = 1 << 15;
            static const int DEBUG_ALL = DEBUG_ENERGY | DEBUG_INTERACTIONS;
            //END OF DEBUG FLAGS
        private:
            System * sys;
            int max_iter;
            double precision;

            double energy(arma::mat const &C_inner, arma::mat const &h0);
            arma::mat hf_matrix(arma::mat const &C_inner, arma::mat const &h0);
            void sort_eigvecval(arma::mat &eigvec, arma::vec &eigval);

        };

    }
}

#endif	/* HF_H */

