/* 
 * File:   System.h
 * Author: chrishir
 *
 * Created on 6. desember 2011, 14:08
 */

#ifndef SYSTEM_H
#define	SYSTEM_H

#include <armadillo>
#include "Basis.h"

namespace toffyrn
{
    namespace libUNK
    {

        /**
         * @brief Class representing a physical system.
         * 
         * Provides the matrix elements f_pq and v_pqrs as well as the Basis
         * these elements belong to.
         * 
         * The simplest subclass must implement the two methods "f_elem()" 
         * and "v_elem()", as well as setting System::basis to a sensible Basis,
         * before calling fillMatrixElements() in the constructor.
         */
        class System
        {
        public:
            /**
             * @brief Default constructor.
             * 
             * Sets basis == NULL.
             */
            System();

            /**
             * @brief Copy constructor.
             * 
             * @param sys Object to copy from.
             */
            System(const System& sys);

            /**
             * @brief Destructor. basis is deleted.
             */
            virtual ~System();
            
            /**
             * @brief Assignment operator.
             * @param sys System to assign from.
             * @return Itself.
             */
            System& operator= (const System& sys);

            /**
             * @brief Get f-matrix element.
             * @param p Index of the outgoing state.
             * @param q Index of the incoming state.
             * @return <p|f|q>
             */
            virtual double f_elem(std::size_t p, std::size_t q) const = 0;

            /**
             * @brief Get v-matrix element, antisymmetrized.
             * @param p Index of the left outgoing state.
             * @param q Index of the right outgoing state.
             * @param r Index of the left incoming state.
             * @param s Index of the right incoming state.
             * @return <pq||rs>
             */
            virtual double v_elem(
                    std::size_t p, std::size_t q,
                    std::size_t r, std::size_t s) const = 0;

            /**
             * @brief Get the address of f_hh matrix.
             * @return f_hh
             */
            virtual arma::mat const * get_f_hh() const
            {
                return &f_hh;
            }

            /**
             * @brief Get the address of f_ph matrix.
             * @return f_ph
             */
            virtual arma::mat const * get_f_ph() const
            {
                return &f_ph;
            }

            /**
             * @brief Get the address of f_pp matrix.
             * @return f_pp
             */
            virtual arma::mat const * get_f_pp() const
            {
                return &f_pp;
            }

            /**
             * @brief Get the address of v_hhhh matrix
             * \f$ \langle \mu||\mu' \rangle_{(\lambda)} \f$,
             * stored as ".at(lmd)(mu,mu')".
             * @sa @ref configurations
             * @return v_hhhh
             */
            virtual std::vector<arma::mat> const * get_v_hhhh() const
            {
                return &v_hhhh;
            }

            /**
             * @brief Get the address of v_phhh matrix 
             * \f$ \langle \nu||\mu \rangle_{(\lambda)} \f$, 
             * stored as ".at(lmd)(nu,mu)".
             * @sa @ref configurations
             * @return v_phhh
             */
            virtual std::vector<arma::mat> const * get_v_phhh() const
            {
                return &v_phhh;
            }

            /**
             * @brief Get the address of v_pphh matrix
             * \f$ \langle \xi||\mu \rangle_{(\lambda)} \f$,
             * stored as ".at(lmd)(xi,mu)".
             * @sa @ref configurations
             * @return pphh
             */
            virtual std::vector<arma::mat> const * get_v_pphh() const
            {
                return &v_pphh;
            }

            /**
             * @brief Get the address of v_phph matrix
             * \f$ \langle \nu||\nu' \rangle_{(\lambda)} \f$,
             *  stored as ".at(lmd)(nu,nu')".
             * @sa @ref configurations
             * @return v_phph
             */
            virtual std::vector<arma::mat> const * get_v_phph() const
            {
                return &v_phph;
            }

            /**
             * @brief Get the address of v_ppph matrix
             * \f$ \langle \xi||\nu \rangle_{(\lambda)} \f$,
             *  stored as ".at(lmd)(xim, nu)".
             * @sa @ref configurations
             * @return v_ppph
             */
            virtual std::vector<arma::mat> const * get_v_ppph() const
            {
                return &v_ppph;
            }

            /**
             * @brief Get the address of v_pppp matrix
             * \f$ \langle \xi||\xi' \rangle_{(\lambda)} \f$,
             *  stored as ".at(lmd)(xi, xi')".
             * @sa @ref configurations
             * @return v_pppp
             */
            virtual std::vector<arma::mat> const * get_v_pppp() const
            {
                return &v_pppp;
            }

            /**
             * @brief Return the Basis this system is in.
             * @return A const pointer to correcsponding Basis.
             */
            virtual Basis const * get_basis() const
            {
                return basis;
            }



        protected:
            /** \f$ f_{lm} \f$ */
            arma::mat f_hh;
            /** \f$ f_{dl} \f$ */
            arma::mat f_ph;
            /** \f$ f_{de} \f$ */
            arma::mat f_pp;
            /** \f$ \langle \mu||\mu' \rangle_{(\lambda)} \f$ */
            std::vector<arma::mat> v_hhhh;
            /** \f$ \langle \nu||\mu \rangle_{(\lambda)} \f$ */
            std::vector<arma::mat> v_phhh;
            /** \f$ \langle \xi||\mu \rangle_{(\lambda)} \f$ */
            std::vector<arma::mat> v_pphh;
            /** \f$ \langle \nu||\nu' \rangle_{(\lambda)} \f$ */
            std::vector<arma::mat> v_phph;
            /** \f$ \langle \xi||\nu \rangle_{(\lambda)} \f$ */
            std::vector<arma::mat> v_ppph;
            /** \f$ \langle \xi||\xi' \rangle_{(\lambda)} \f$ */
            std::vector<arma::mat> v_pppp;
            /** The Basis we are working in */
            Basis * basis;

            /**
             * @brief Fill all matrix elements: f_hh, f_ph, f_pp, v_hhhh,
             * v_phhh, v_pphh, v_phph, v_ppph, v_pppp.
             * 
             * The default implementation uses the supplied Basis, and the 
             * two functions f_elem() and v_elem().
             */
            virtual void fillMatrixElements();


        private:

        };

    }
}

#endif	/* SYSTEM_H */

