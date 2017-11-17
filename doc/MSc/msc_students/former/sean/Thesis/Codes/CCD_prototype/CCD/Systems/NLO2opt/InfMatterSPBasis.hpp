// Copyright (c) 2015-2016, Justin Gage Lietz
// All rights reserved.

#ifndef InfMatterBasis_hpp_
#define InfMatterBasis_hpp_

#include "SPBasis.hpp"
#include <unordered_map>

struct ChannelBundle{
	int chanNx;
	int chanNy;
	int chanNz;
	int chanSz;
	int chanTz;
};

class InfMatterSPBasis: public SPBasis{
	public:
		// Physical constants and relevant basis quantities
		const double hbarc = 197.3269788; // MeVfm
		//double m_neutronc2 = 939.565378; // MeV // accurate value
		const double m_neutronc2 = 939.565; // MeV Gaute's value
		//double m_protonc2 = 938.272046; // MeV
		//double m_symmc2 = (m_neutronc2 + m_protonc2)*0.5;
		double massc2;
		double boxLength;
		ChannelBundle * chanValue;
		ChannelBundle * chanModValue;

		// density of the matter in fm^-3
		double density;

		// Number of particles species, 1 for neutrons, 2 for protons and neutrons
		int tzMax;

		// kind of hacked parameter, keeps track of how many channels
		int chanTzMax;
		// same but for the mod Channels
		int modChanTzMax;

		// highest single nx, ny, or nz value in basis
		int nMax;
		int nParticleMax;

		// number of non-empty fermi spheres in the calculations
		// number of non-empty fermi spheres to fill with particles
		std::size_t nShells;
		std::size_t nParticleShells;


		// max sp energy
		// max sp energy for holes
		std::size_t EMax;
		std::size_t fermiEnergy;


		// hash to map from quantum numbers -> sp index
		// need to use the keyGen_spHash first to generate key
		// feed 5 quantum numbers into keyGen to get a key,
		// then use the key to get the sp index
		std::unordered_map<std::size_t, std::size_t> sp_hash;


		InfMatterSPBasis(std::size_t basisIndicator, double densityIn, std::size_t tzMaxIn, std::size_t nShellsIn, std::size_t nParticleShellsIn);
		void generateIndexMap();
		void generateBasis();
		void setUpTwoStateChannels();
		void setUpChannelValues();

		void printBasis();
		void deallocate();
		int spinExchangeMtxEle(int i, int j, int k, int l);
		int kron_del(int i, int j);

		std::size_t keyGen_spHash(int nx, int ny, int nz, int sz, int tz);
		std::size_t TBchanIndexFunction(std::size_t p, std::size_t q);
		std::size_t TBmodChanIndexFunction(std::size_t p, std::size_t q);
		std::size_t spIndex_from3Body(std::size_t q_inv, std::size_t r, std::size_t s);
		int spIndexExists_from3Body(std::size_t q_inv, std::size_t r, std::size_t s);

		double calcTBME(std::size_t p, std::size_t q, std::size_t r, std::size_t s);
		double calc_TBME_not_antisym(std::size_t p, std::size_t q, std::size_t r, std::size_t s);
};

#endif /* INFMATTER_HPP */
