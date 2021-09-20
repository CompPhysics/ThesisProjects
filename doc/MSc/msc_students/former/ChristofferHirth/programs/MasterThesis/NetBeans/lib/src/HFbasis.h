/* 
 * File:   HFbasis.h
 * Author: toffyrn
 *
 * Created on 16. februar 2012, 17:21
 */

#ifndef HFBASIS_H
#define	HFBASIS_H

#include "System.h"

namespace toffyrn
{
    namespace libUNK
    {

        class HFbasis : public System
        {
        public:
            HFbasis(System const * originalSys, arma::mat const &Coeff);
            virtual ~HFbasis();
            virtual double f_elem(std::size_t p, std::size_t q) const;
            virtual double v_elem(
                    std::size_t p, std::size_t q,
                    std::size_t r, std::size_t s) const;
        protected:
            virtual void transformElements(arma::mat const &C);
            double transformOneElem(int p, int q, int r, int s, arma::mat const &C);
        private:

        };

    }
}
#endif	/* HFBASIS_H */

