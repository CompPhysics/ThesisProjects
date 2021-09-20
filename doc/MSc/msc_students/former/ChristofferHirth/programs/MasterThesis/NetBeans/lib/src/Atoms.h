/* 
 * File:   Atoms.h
 * Author: chrishir
 *
 * Created on 7. desember 2011, 12:49
 */

#ifndef ATOMS_H
#define	ATOMS_H

#include "System.h"
#include <armadillo>

namespace toffyrn
{
    namespace libUNK
    {

        class Atoms : public System
        {
        public:
            Atoms(int numElec = 2);
            virtual ~Atoms();

            virtual double f_elem(std::size_t p, std::size_t q) const;
            virtual double v_elem(
                    std::size_t p, std::size_t q,
                    std::size_t r, std::size_t s) const;

            double spatial_integral(int p, int q, int r, int s) const;
        private:
            double Z ;

        };


    }
}


#endif	/* ATOMS_H */

