/* 
 * File:   HOBasis.h
 * Author: toffyrn
 *
 * Created on 29. desember 2011, 15:23
 */

#ifndef HOBASIS_H
#define	HOBASIS_H

#include "Basis.h"

namespace toffyrn
{
    namespace libUNK
    {

        class HOBasis : public Basis
        {
        public:
            HOBasis(int filledR, int totalR);
            virtual arma::ivec stateMap(int p) const;
            virtual int lambda_neg(int p, int q) const;
            virtual int lambda_1p(int p) const;
            virtual int dim_lmd_1p() const;

            virtual HOBasis * clone() const
            {
                return new HOBasis(*this);
            }
        protected:
            virtual void initStateMap(int nStates);
            virtual void initMappings();
        private:
            int shells_tot;
            int shells_filled;
            arma::imat statemapMat;
        };

    }
}


#endif	/* HOBASIS_H */

