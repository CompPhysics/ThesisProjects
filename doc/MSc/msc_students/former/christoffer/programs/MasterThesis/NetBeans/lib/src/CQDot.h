/* 
 * File:   CQDot.h
 * Author: toffyrn
 *
 * Created on 19. desember 2011, 15:05
 */

#ifndef CQDOT_H
#define	CQDOT_H

#include "System.h"
#include <algorithm>
#include <armadillo>
#include "HOBasis.h"
#include <iostream>
#include <map>

namespace toffyrn
{
    namespace libUNK
    {

        class CQDot : public System
        {
        public:
            CQDot(int filledR, int shellsR, std::string fname, bool nmms_file = false);
            virtual ~CQDot();
            virtual double f_elem(std::size_t p, std::size_t q) const;
            virtual double v_elem(
                    std::size_t p, std::size_t q,
                    std::size_t r, std::size_t s) const;

        private:
            std::string filename;
            void fillMatrixElements(bool nmms_file);
            void putOneElement(int p, int q, int r, int s, double element);
            void readNMMSfile(std::ifstream &file);

        };

        class Generate
        {
        public:
            Generate();
            void genFile(int shells, std::string fname);

        private:
            arma::vec logfacVec;
            arma::vec lgammaVec;
            arma::mat lgprd3_compon;

            // compute a part of the Coulomb matrix element: <12|V|34>
            double coulomb(const int, const int, const int, const int, const int, const int, const int, const int) const;

            //compute (-1)^k
            int minusPower(const int) const;

            // computes log(n!)
            double LogFac(const int, bool calc);
            double LogFac(const int) const;

            // Computes the first ratio in the Asinimovas expression
            double LogRatio1(const int, const int, const int, const int) const;

            // Computes the 2nd ratio in the Asinimovas expression
            double LogRatio2(const int) const;

            // computes first product of indices in the Anisimovas/Matulis expression
            double Product1(const int, const int, const int, const int, const int, const int, const int, const int) const;

            // Computes the log of the 2nd product in the Asinimovas expression
            double LogProduct2(const int, const int, const int, const int, const int, const int, const int, const int, const int, const int, const int, const int) const;

            // Computes the log of the 3rd product in the Asinimovas expression
            double LogProduct3(const int, const int, const int, const int, const int, const int, const int, const int) const;

            // The function lgamma() computes the logarithm of the gamma function of real argument x
            double lgamma(int xTimes2) const;
            double lgamma(int xTimes2, bool calc);
        };
    }
}

#endif	/* CQDOT_H */

