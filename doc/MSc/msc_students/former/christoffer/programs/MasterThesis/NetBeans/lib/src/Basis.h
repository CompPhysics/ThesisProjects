/* 
 * File:   Basis.h
 * Author: chrishir
 *
 * Created on 7. desember 2011, 10:49
 */

#ifndef BASIS_H
#define	BASIS_H

#include <armadillo>

namespace toffyrn
{
    namespace libUNK
    {

        /**
         * @brief Representing a basis, with a number of occupied states, 
         * holestates, and a number unoccupied states, particlestates.
         * 
         * In addition a mapping from two quantum numbers into configurations 
         * and channels is represented.
         */
        class Basis
        {
        public:
            /**
             * @brief Empty constructor. Sets this to an empty basis (0 states).
             */
            Basis();

            /**
             * @brief Constructor initializing a basis with nP (virtual) particlestates, 
             * and nH (occupied) holestates.
             * 
             * This baseclass implementation calls initMappings() to set up mappings.
             * 
             * @param nH Number of holestates, the actual number of particles in a system.
             * @param nP Numper of virtual particlestates.
             */
            Basis(int nH, int nP);

            /**
             * @brief Destructor. No objects needs explicit cleaning.
             */
            virtual ~Basis();

            /**
             * @brief Clone method. Replaces copy constructor, to ensure copying correct class type.
             * @return A copy of this "Basis" object.
             */
            virtual Basis * clone() const
            {
                return new Basis(*this);
            }

            /**
             * @brief Getter for nH, the number of holestates.
             * @return The number of holestates, the actual number of particles in a system.
             */
            int get_nH() const
            {
                return nH;
            }

            /**
             * @brief Getter for nP, the number of virtual particlestates.
             * @return The number of virtual particlestates in the basis.
             */
            int get_nP() const
            {
                return nP;
            }

            /**
             * @brief Getter for general mapping 
             * \f$(p,q) \rightarrow (\lambda,\pi)\f$.
             * @sa map_pq_lmdPI
             * @return unsigned \f$(2\times (n_H+n_P)^2)\f$ matrix map_pq_lmdPI.
             */
            arma::umat const * get_map_pq_lmdPI() const
            {
                return &map_pq_lmdPI;
            }

            /**
             * @brief Getter for general mapping 
             * \f$ (\lambda,\pi) \rightarrow (p,q) \f$.
             * @sa map_lmdPI_pq
             * @return map_lmdPI_pq, a STL vector filled with an unsigned 
             * vector for each channel \f$\lambda\f$.
             */
            std::vector<arma::uvec> const * get_map_lmdPI_pq() const
            {
                return &map_lmdPI_pq;
            }

            /**
             * @brief Getter for hole-hole mapping
             * \f$(l,m) \rightarrow (\lambda,\mu)\f$.
             * @sa map_lm_lmdMU
             * @return unsigned matrix \f$(2\times n_H^2)\f$ map_lm_lmdMU.
             */
            arma::umat const * get_map_lm_lmdMU() const
            {
                return &map_lm_lmdMU;
            }

            /**
             * @brief Getter for hole-hole mapping
             * \f$(\lambda,\mu) \rightarrow (l,m)\f$.
             * @sa map_lmdMU_lm
             * @return map_lmdMU_lm, a STL vector filled with an unsigned 
             * vector for each channel \f$\lambda\f$.
             */
            std::vector<arma::uvec> const * get_map_lmdMU_lm() const
            {
                return &map_lmdMU_lm;
            }

            /**
             * @brief Getter for particle-hole mapping
             * \f$(d,l) \rightarrow (\lambda,\nu)\f$.
             * @sa map_dl_lmdNU
             * @return unsigned matrix \f$(2\times n_H n_P)\f$ map_dl_lmdNU.
             */
            arma::umat const * get_map_dl_lmdNU() const
            {
                return &map_dl_lmdNU;
            }

            /**
             * @brief Getter for particle-hole mapping
             * \f$(\lambda,\nu) \rightarrow (d,l)\f$.
             * @sa map_lmdNU_dl
             * @return map_lmdNU_dl, a STL vector filled with an unsigned
             * vector for each channel \f$\lambda\f$.
             */
            std::vector<arma::uvec> const * get_map_lmdNU_dl() const
            {
                return &map_lmdNU_dl;
            }

            /**
             * @brief Getter for particle-particle mapping
             * \f$(d,e) \rightarrow (\lambda,\xi)\f$.
             * @sa map_de_lmdXI
             * @return unsigned matrix \f$(2\times n_P^2)\f$ map_de_lmdXI.
             */
            arma::umat const * get_map_de_lmdXI() const
            {
                return &map_de_lmdXI;
            }

            /**
             * @brief Getter for particle-particle mapping
             * \f$(\lambda,\xi) \rightarrow (d,e)\f$.
             * @sa map_lmdXI_de
             * @return map_lmdXI_de, a STL vector filled with an unsigned 
             * vector for each channel \f$\lambda\f$.
             */
            std::vector<arma::uvec> const * get_map_lmdXI_de() const
            {
                return &map_lmdXI_de;
            }

            /**
             * @brief The number of channels existing for two-particle configurations; \f$ dim(\lambda)\f$.
             * @return \f$ dim(\lambda)\f$
             */
            virtual size_t dim_lmd_2p() const
            {
                return map_lmdPI_pq.size();
            }

            /**
             * @brief Get the quantum numbers corresponding to a specific state.
             * 
             * Each state within a basis is labelled by one integer from 
             * 0 up to the total number of states minus one.
             * This function allows a specific basis to translate this into 
             * a set of quantum numbers, usually three/four for 2D/3D problems.
             * 
             * The default implementation is empty, returning a zero-sized vector.
             * 
             * @param state A integer within [0,nStates-1].
             * @return A vector of integers, each labelling one quantum number.
             */
            virtual arma::ivec stateMap(int state) const
            {
                arma::ivec map = arma::zeros<arma::ivec > (0);
                return map;
            }

            /**
             * @brief Returning the "inverse" channel \f$ pq^{-1} \f$.
             * Dealing with channels for additive quantum numbers, here
             * q is subtracted instead of added.
             * 
             * The dimensionality of the inverse channels should be equal to 
             * the two-particle channels \f$ dim(\lambda) \f$.
             * 
             * This default Basis class operates with one channel (0).
             * 
             * @param p State for particle one.
             * @param q State for particle two, inverse.
             * @return The "inverse" channel for \f$ pq^{-1} \f$.
             */
            virtual int lambda_neg(int p, int q_neg) const
            {
                return 0;
            }

            /**
             * @brief Get the channel corresponding to a single particlestate.
             * @return The single-particle channel \f$ \lambda_{1p} \f$.
             */
            virtual int lambda_1p(int p) const
            {
                return 0;
            }

            /**
             * @brief Dimensionality of \f$ dim(\lambda_{1p}) \f$ -- the number 
             * transition-channels through one-particle interactions.
             * @return \f$ dim(\lambda_{1p}) \f$.
             */
            virtual int dim_lmd_1p() const
            {
                return 0;
            }

        protected:
            /** The number of holestates in our basis. 
             * (This is the actual number of particles.)
             */
            int nH;
            /** The number of particlestates in our basis. 
             * (This is the unoccupied/virtual states.)
             */
            int nP;

            /**
             * Containing the mapping from two general states \f$(p,q)\f$ into
             * \f$(\lambda,\pi)\f$.
             * The first row contains \f$\lambda\f$ and the second contains 
             * \f$\pi\f$, corresponding to row \f$(p + q * (n_P+n_H))\f$.
             */
            arma::umat map_pq_lmdPI;

            /** 
             * Containing the mapping from a channel \f$\lambda\f$ and
             * a general configuration \f$\pi\f$ into two general states
             * \f$(p,q)\f$.
             * An uvec is contained in each element of this vector, so
             * that ".at(lmd)(pi)" has the value \f$(p+q*(n_P+n_H))\f$.
             */
            std::vector<arma::uvec> map_lmdPI_pq;

            /**
             * Containing the mapping from two hole states \f$(l,m)\f$ into
             * \f$(\lambda,\mu)\f$.
             * The first row contains \f$\lambda\f$ and the second contains 
             * \f$\mu\f$, corresponding to row \f$(l + m * n_H)\f$.
             */
            arma::umat map_lm_lmdMU;

            /** 
             * Containing the mapping from a channel \f$\lambda\f$ and
             * a hole-hole configuration \f$\mu\f$ into two hole states
             * \f$(l,m)\f$.
             * An uvec is contained in each element of this vector, so
             * that ".at(lmd)(mu)" has the value \f$(l+m*n_H)\f$.
             */
            std::vector<arma::uvec> map_lmdMU_lm;

            /**
             * Containing the mapping from one particle one hole state \f$(d,l)\f$
             * into \f$(\lambda,\nu)\f$.
             * The first row contains \f$\lambda\f$ and the second contains \f$\nu\f$,
             * corresponding to row \f$(d + l* n_P)\f$.
             */
            arma::umat map_dl_lmdNU;

            /** 
             * Containing the mapping from a channel \f$\lambda\f$ and
             * a particle-hole configuration \f$\nu\f$ into one particle state, and
             * one hole state \f$(d,l)\f$.
             * An uvec is contained in each element of this vector, so
             * that ".at(lmd)(nu)" has the value \f$(d+l*n_P)\f$.
             */
            std::vector<arma::uvec> map_lmdNU_dl;

            /** 
             * Containing the mapping from two particle states \f$(d,e)\f$ into
             * \f$(\lambda,\xi)\f$.
             * The first row contains \f$\lambda\f$ values and the second 
             * contains \f$\xi\f$, corresponding to the row \f$(d + e * n_P)\f$.
             */
            arma::umat map_de_lmdXI;

            /** 
             * Containing the mapping from a channel \f$\lambda\f$ and
             * a particle-particle configuration \f$\xi\f$ into two particle states
             * \f$(d,e)\f$.
             * An uvec is contained in each element of this vector, so
             * that ".at(lmd)(xi)" has the value \f$(d+e*n_P)\f$.
             */
            std::vector<arma::uvec> map_lmdXI_de;

            /**
             * @brief Initialize all mappings.
             * 
             * This default baseclass implementation creates only one channel
             * \f$(\lambda)\f$, and maps as follow:
             * @li Free configurations \f$ \pi = p+q*(n_H+n_P) \f$.
             * @li Hole-hole configurations \f$ \mu = l+m*n_H \f$.
             * @li Particle-hole configurations \f$ \nu = d+l*n_P \f$.
             * @li Particle-particle configurations \f$ \xi = d+e*n_P \f$.
             * 
             * Resulting in the need of storing the complete interaction matrices.
             * Although not memory efficient, it is general, and can be used on any system.
             */
            virtual void initMappings();

        private:

        };
    }
}

#endif	/* BASIS_H */


/**
 * @page configurations Configurations and channels
 * Info about tp configurations.
 * 
 * //TODO: More details and info about these transition-channels and two-particle configurations.
 */
