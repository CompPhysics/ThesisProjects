// Copyright (c) 2015-2016, Justin Gage Lietz
// All rights reserved.

#ifndef SPBasis_hpp_
#define SPBasis_hpp_

#include <iomanip>

/**
 * This SPBasis class is the interface for solver methods (like CCD)
 * to interact with any arbitrary basis. The goal is to have
 * the minimum necessary member variables and functions
 * without exposing the individual basis implementation.
 * It would be nice to standardize this so that new bases
 * can be added and swapped in easily. The basic idea in mind is that
 * for a solver that uses up to 2-body forces, almost all of the two-body
 * matrix elements (TMBEs) <p,q|v|r,s> will be zero. One way to store
 * the TBMEs to reduce memory requirements, is to take advantage of the
 * known symmetries of the interaction, and not store the matrix elements
 * that break conservation laws. For example, if v conserves spin projection
 * (ms), then ms_p + ms_q = ms_r + ms_s. So any combinations p,q,r,s that
 * do not obey this have <p,q|v|r,s> = 0, and we do not store them. What we
 * do, is to store only blocks of the full v matrix that obey conservation
 * laws. Abstractly, we can say that for single particle (s.p.) states p,q
 * with good quantum numbers, adding their quantum numbers p + q defines a
 * set of two-body quantum numbers, for which r + s must equal to get a
 * non-zero TBME. So the non-zero blocks of the v matrix are defined by
 * two-body quantum numbers p + q. I will do a full write-up about this
 * in a PDF at some point and include it in the folder.
 */

struct PQBundle{
	std::size_t p;
	std::size_t q;

	bool operator <(const PQBundle other) const {
		if ( p < other.p ){
			return 1;
		} else if ( p == other.p ){
			return q < other.q;
		}
		return 0;
	}
};

struct PQRBundle{
	std::size_t p;
	std::size_t q;
	std::size_t r;

	bool operator <(const PQRBundle other) const {
		if ( p < other.p ){
			return 1;
		} else if ( p == other.p ){
			if ( q < other.q ){
				return 1;
			} else if ( q == other.q ){
				return r < other.r;
			}
		}
		return 0;
	}
};

/**
 * ChannelDims is a small struct that carries row and col dimensions
 * for a given symmetry block.
 */
struct ChannelDims{
	std::size_t ppDim;
	std::size_t hhDim;
	std::size_t hpDim;
	std::size_t phDim;
};

/**
 * maps are a small structs that carry the single particle
 * indices for a given symmetry block.
 */
struct ChannelMaps{
	// These store the sp indices that are in a given channel
	PQBundle * ppMap;
	PQBundle * hhMap;
	PQBundle * hpMap;
	PQBundle * phMap;
};

struct ThreeBodyChannelDims{
	std::size_t phhDim;
	std::size_t hppDim;
};
struct ThreeBodyChannelMaps{
	PQRBundle * phhMap;
	PQRBundle * hppMap;
};


/**
 * The SPBasis class is mostly an abstract class. Some member functions
 * are implemented, but most of the basis functionality is implemented
 * in a basis specific class that inherits from SPBasis. Look at
 * InfMatterSPBasis for an example of how I construct a basis.
 */
class SPBasis{
	public:
		// All bases need these members to interact with CCD
		std::size_t basisIndicator; //cheat parameter for now?
		std::size_t nSpstates;
		std::size_t nParticles;
		std::size_t nChannels;
		int ** indexMap;
		double * spEnergy;
		ChannelDims * chanDims;
		ChannelDims * chanModDims;
		ChannelMaps * chanMaps;
		ChannelMaps * chanModMaps;
		ThreeBodyChannelDims * threeBodyChanDims;
		ThreeBodyChannelMaps * threeBodyChanMaps;
		// delete these. Move destructor into SPBasis???
		std::size_t * channelIndices;
		std::size_t * inverseChannelIndices;
		std::size_t * modChannelIndices;
		std::size_t * inverseModChannelIndices;
		std::size_t * spChannelIndices;
		std::size_t * inverseSPChannelIndices;

		std::size_t maxChannelSize_h_hpp;
		std::size_t maxChannelSize_phh_p;
		std::size_t maxChannelSize_hhhh;
		std::size_t maxChannelSize_hhhp;
		std::size_t maxChannelSize_hhpp;
		std::size_t maxChannelSize_hpph;
		std::size_t maxChannelSize_hphp;
		std::size_t maxChannelSize_hppp;
		std::size_t maxChannelSize_pppp;

		/** SPBasis needs a constructor to set up everything.
		 * This constructor does nothing, as most information
		 * as most of this information is basis specific. At
		 * the very least, nParticles, nSpstates, nChannels
		 * indexMap and spEnergy should be constructed here.
		 */
		SPBasis() {}

		/**
		 * Might remove this at some point. All it does is
		 * call setUpChannelValues(), Dims(), and Maps().
		 */
		virtual void setUpTwoStateChannels() = 0;

		/**
		 * Basis specific. Sets up a basis specific struct
		 * that contains the values of the 2-body quantum
		 * number for all of the channels. Right now this is
		 * not used, but I can see a use for having the
		 * 2-body quantum numbers stored somewhere.
		 */
		virtual void setUpChannelValues() = 0;


		/**
		 * This function determines a mapping from s.p. indices
		 * into a regular channel. |p,q> has symmetry p+q. This
		 * is necessary for setUpChannelDims() and ...Maps().
		 * @param p is the first index of a ket state |p,q>
		 * @param q is the second index of a ket state |p,q>
		 */
		virtual std::size_t TBchanIndexFunction(std::size_t p, std::size_t q) = 0;

		/**
		 * These functions determine a mapping from s.p. indices
		 * into a modified channel, or modChannel. |p,q^{-1}> has
		 * symmetry p-q. This is necessary for setUpChannelDims()
		 * and ...Maps().
		 * @param p is the first index of a ket state |p,q^{-1}>
		 * @param q is the second index of a ket state |p,q^{-1}>
		 */
		virtual std::size_t TBmodChanIndexFunction(std::size_t p, std::size_t q) = 0;

		/**
		 * These functions determine a mapping from s.p. indices
		 * into a three body channel. <p,q|v|r,s> to <p|v|q^{-1},r,s>
		 * has symmetry r + s - q. This is necessary for setUpChannelDims()
		 * and ...Maps().
		 * @param q_inv is the first index of a ket state |q^{-1},r,s>
		 * @param r is the second index of a ket state |q^{-1},r,s>
		 * @param s is the third index of a ket state |q^{-1},r,s>
		 */
		virtual std::size_t spIndex_from3Body(std::size_t q_inv, std::size_t r, std::size_t s) = 0;

		/**
		 * This functions returns 1 if a mapping from s.p. indices exists,
		 * 0 if it does not.
		 * @param q_inv is the first index of a ket state |q^{-1},r,s>
		 * @param r is the second index of a ket state |q^{-1},r,s>
		 * @param s is the third index of a ket state |q^{-1},r,s>
		 */
		virtual int spIndexExists_from3Body(std::size_t q_inv, std::size_t r, std::size_t s) = 0;

		// free all the memory.
		virtual void deallocate() = 0;
		// quite verbose print
		virtual void printBasis() = 0;
		// Two-Body Matrix Elements are dependent on a basis
		// so this has to be carried around
		// Default is ANTI-SYMMETRIC
		virtual double calcTBME(std::size_t p, std::size_t q, std::size_t r, std::size_t s) = 0;

		// Non-antisymmetric TBMEs
		virtual double calc_TBME_not_antisym(std::size_t p, std::size_t q, std::size_t r, std::size_t s) = 0;

		// Now are the functions that are implemented in SPBasis.cpp.
		// Should be the same across all bases. Would be nice to move
		// as much functionality as possible here, so that creating
		// new bases requires less work.

		/**
		 * Takes the default s.p. basis energies, and makes them the
		 * normal ordered s.p. energies since CCD assumes the Hamilonian
		 * is normal ordered.
		 * e_i_normalordered = e_i + sum_j <i,j|v|i,j>
		 * e_a_normalordered = e_a + sum_j <a,j|v|a,j>
		 */
		void rotateSpEnergiesToNormalOrdered();

		/**
		 * Calculates the reference energy. Assumes that the
		 * s.p. energies are *NOT* normal ordered.
		 * E_ref = sum_i e_i + 0.5*sum_i,j <i,j|v|i,j>
		 */
		double referenceEnergy();

		/**
		 * Generates chanDims, chanModDims and threeBodyChanDims
		 * by using a basis specific implementation of the
		 * index functions.
		 */
		void setUpChannelDims();

		/**
		 * Generates chanMaps, chanModMaps and threeBodyChanMaps
		 * by using a basis specific implementation of the
		 * index functions. Needs setUpChannelDims to be ran
		 * first
		 */
		void setUpChannelMaps();

		/**
		 * If the channelIndices have been permuted, this
		 * sets up the corresponding inverseChannelIndices
		 * array, and rotates the chanDims accordingly.
		 */
		void setUpInverseChannelIndicesAndPermuteDims();
};

#endif /*ABSTRACTSPBASIS_HPP*/
