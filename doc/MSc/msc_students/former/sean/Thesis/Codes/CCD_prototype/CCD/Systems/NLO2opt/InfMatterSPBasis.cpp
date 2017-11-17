// Copyright (c) 2014-2016, Justin Gage Lietz
// All rights reserved.

// Need to write a PDF about this biz.

#include "InfMatterSPBasis.hpp"
#include "util.hpp"
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cmath>

InfMatterSPBasis::InfMatterSPBasis(std::size_t basisIndicatorIn, double densityIn, std::size_t tzMaxIn,
		std::size_t nShellsIn, std::size_t nParticleShellsIn){

	this->basisIndicator = basisIndicatorIn;
	this->density = densityIn;
	this->tzMax = tzMaxIn;
	this->nShells = nShellsIn;
	this->nParticleShells = nParticleShellsIn;
	if(this->tzMax == 1){
		this->chanTzMax = -2;
		this->modChanTzMax = 0;
	} else if(this->tzMax == 2){
		this->chanTzMax = 2;
		this->modChanTzMax = 2;
	} else {
		this->chanTzMax = 8000;
		this->modChanTzMax = 8000;
		std::cout << "wrong tzMax" << std::endl;
	}

	// Primary goal: Find the number of sp states
	// also find EMax and nMax
	this->EMax = 0;
	for(std::size_t i = 1; i < this->nShells; this->EMax++){
		if(isASumOfThreeSquares(this->EMax+1)){
			i++;
		}
	}

	this->nMax = smallestSquareRootAtLeast(this->EMax);

	this->nSpstates = 0;
	int nMaxActual = 0;
	for( int nx = -this->nMax; nx <= this->nMax; nx++){
		for( int ny = -this->nMax; ny <= this->nMax; ny++){
			for( int nz = -this->nMax; nz <= this->nMax; nz++){
				if( (std::size_t) nx*nx + ny*ny + nz*nz <= this->EMax){
					if(nx > nMaxActual){
						nMaxActual = nx;
					}
					for( int sz = -1; sz <= 1; sz = sz+2){
						for( int tz = -1; tz <= this->tzMax - 1; tz = tz+2){
							this->nSpstates++;
						} // end tz loop
					} // end sz loop
				}
			} // end nz loop
		} // end ny loop
	} // end nx loop

	this->nMax = nMaxActual;

	// now find the number of particles
	this->fermiEnergy = 0;
	for(std::size_t i = 1; i < this->nParticleShells; this->fermiEnergy++){
		if(isASumOfThreeSquares(this->fermiEnergy+1)){
			i++;
		}
	}

	this->nParticleMax = smallestSquareRootAtLeast(this->fermiEnergy);

	this->nParticles = 0;
	for( int nx = -this->nParticleMax; nx <= this->nParticleMax; nx++){
		for( int ny = -this->nParticleMax; ny <= this->nParticleMax; ny++){
			for( int nz = -this->nParticleMax; nz <= this->nParticleMax; nz++){
				if(nx > nMaxActual){
					nMaxActual = nx;
				}
				if( (std::size_t) nx*nx + ny*ny + nz*nz <= this->fermiEnergy){
					for( int sz = -1; sz <= 1; sz = sz+2){
						for( int tz = -1; tz <= this->tzMax - 1; tz = tz+2){
							this->nParticles++;
						} // end tz loop
					} // end sz loop
				}

			} // end nz loop
		} // end ny loop
	} // end nx loop
	this->nParticleMax = nMaxActual;

	// Now build relevant sp structures
	this->generateIndexMap();
	this->generateBasis();
}


void InfMatterSPBasis::generateIndexMap(){

	this->indexMap = new int* [this->nSpstates];
	for(std::size_t i = 0; i < this->nSpstates; i++){
		this->indexMap[i] = new int[5];
	}

	std::size_t *shells = new std::size_t[this->EMax + 2];
	std::size_t E;

	for( E = 0; E <= this->EMax + 1; E++){
		shells[E] = 0;
	}

	//Determine how many single particle states are in each shell.
	for( int nx = -this->nMax; nx <= this->nMax; nx++){
		for( int ny = -this->nMax; ny <= this->nMax; ny++){
			for( int nz = -this->nMax; nz <= this->nMax; nz++){
				E = nx*nx + ny*ny + nz*nz;
				if( E <= this->EMax){
					for( int sz = -1; sz <= 1; sz = sz+2){
						for( int tz = -1; tz <= this->tzMax - 1; tz = tz+2){
							shells[E + 1]++;
						} // end tz loop
					} // end sz loop
				} // end if
			} // end nz loop
		} // end ny loop
	} // end nx loop

	//Determine offsets into index array for each energy level
	for( E = 1; E <= this->EMax; E++){
		shells[E] += shells[E - 1];
	}

	// Construct the index map
	// and inversely a sp index hash
	int hashKey;
	for( int nx = -this->nMax; nx <= this->nMax; nx++){
		for( int ny = -this->nMax; ny <= this->nMax; ny++){
			for( int nz = -this->nMax; nz <= this->nMax; nz++){
				E = nx*nx + ny*ny + nz*nz;
				if( E <= this->EMax){
					for( int sz = -1; sz <= 1; sz = sz+2){
						for( int tz = -1; tz <= this->tzMax - 1; tz = tz+2){
							this->indexMap[shells[E]][0] = nx;
							this->indexMap[shells[E]][1] = ny;
							this->indexMap[shells[E]][2] = nz;
							this->indexMap[shells[E]][3] = sz;
							this->indexMap[shells[E]][4] = tz;
							hashKey = this->keyGen_spHash(nx,ny,nz,sz,tz);
							sp_hash[hashKey] = shells[E];
							shells[E]++;
						} // end tz loop
					} // end sz loop
				} // end if
			} // end nz loop
		} // end ny loop
	} // end nx loop
	delete[] shells;

} // end generateIndexMap


// requires generate index map first
void InfMatterSPBasis::generateBasis(){
	this->spEnergy = new double[this->nSpstates];
	if( this->tzMax == 1){
		this->massc2 = this->m_neutronc2;
	} else if( this->tzMax == 2){
		//massc2 = m_symmc2;
		this->massc2 = this->m_neutronc2;
	} else {
		std::cout << "tzMax = 1 or 2 ONLY" << std::endl;
	}
	double prefactor = this->hbarc * this->hbarc/(2.*this->massc2);
	this->boxLength= pow(nParticles/density, 1./3.);
	//double E;
	for(std::size_t i = 0; i < this->nSpstates; i++){
		this->spEnergy[i] = (prefactor*4*M_PI*M_PI/(this->boxLength*this->boxLength) ) *
			(this->indexMap[i][0]*this->indexMap[i][0] +
			 this->indexMap[i][1]*this->indexMap[i][1] +
			 this->indexMap[i][2]*this->indexMap[i][2]);
	}
} // end generateBasis

void InfMatterSPBasis::setUpTwoStateChannels(){
	this->setUpChannelValues();
	this->setUpChannelDims();
	this->setUpChannelMaps();
}

// This is not really used currently. Keep this around for debugging maybe.
void InfMatterSPBasis::setUpChannelValues(){
	int channelNmax = 2*this->nMax;
	int channelTzMax = this->chanTzMax;
	int modChannelTzMax = this->modChanTzMax;
	int channelSzMax = 2;

	// first count how many unique channels there are
	// could potentially nuke these loops and just overallocate
	// upper bound is known

	// Before we multiply these numbers together, we cast them to
	// size_t to avoid overflowing int.
	this->nChannels = ((std::size_t)(2*channelNmax+1))
		* ((std::size_t)(2*channelNmax+1))
		* ((std::size_t)(2*channelNmax+1))
		* ((std::size_t)(3*(2+channelTzMax/2)));
	this->chanValue = new ChannelBundle[this->nChannels];
	this->chanModValue = new ChannelBundle[this->nChannels];

	std::size_t channelCount = 0;
	for(int ichanNx = -channelNmax; ichanNx <= channelNmax; ichanNx++){
		for(int ichanNy = -channelNmax; ichanNy <= channelNmax; ichanNy++){
			for(int ichanNz = -channelNmax; ichanNz <= channelNmax; ichanNz++){
				for(int ichanSz = -channelSzMax; ichanSz <= channelSzMax; ichanSz = ichanSz +2){
					for(int ichanTz = -2; ichanTz <= channelTzMax; ichanTz = ichanTz +2){
						this->chanValue[channelCount].chanNx = ichanNx;
						this->chanValue[channelCount].chanNy = ichanNy;
						this->chanValue[channelCount].chanNz = ichanNz;
						this->chanValue[channelCount].chanSz = ichanSz;
						this->chanValue[channelCount].chanTz = ichanTz;

						channelCount++;
					}
				}
			}
		}
	}

	channelCount = 0;
	for(int ichanNx = -channelNmax; ichanNx <= channelNmax; ichanNx++){
		for(int ichanNy = -channelNmax; ichanNy <= channelNmax; ichanNy++){
			for(int ichanNz = -channelNmax; ichanNz <= channelNmax; ichanNz++){
				for(int ichanSz = -channelSzMax; ichanSz <= channelSzMax; ichanSz = ichanSz +2){
					for(int imodChanTz = -modChannelTzMax; imodChanTz <= modChannelTzMax; imodChanTz = imodChanTz +2){
						this->chanModValue[channelCount].chanNx = ichanNx;
						this->chanModValue[channelCount].chanNy = ichanNy;
						this->chanModValue[channelCount].chanNz = ichanNz;
						this->chanModValue[channelCount].chanSz = ichanSz;
						this->chanModValue[channelCount].chanTz = imodChanTz;

						channelCount++;
					}
				}
			}
		}
	}
}


void InfMatterSPBasis::printBasis(){
	std::cout << "mass*c^2: " << this->massc2 << "MeV" << std::endl;
	std::cout << "L = " << this->boxLength << "fm V= " << this->boxLength*this->boxLength*this->boxLength << "fm^3" << std::endl;
	std::cout << "k_fermi " << std::setprecision(16) <<
		pow((6.*M_PI*M_PI*this->density/(2.*this->tzMax)),1./3.) << "fm^-1" << std::endl;
	std::cout << "SPBasis:" << std::endl;
	std::cout << "i nx ny nz sz tz E" << std::endl;
	for(std::size_t p = 0; p < this->nSpstates; p++ ) {
		std::cout << p << " " << this->indexMap[p][0] << " " << this->indexMap[p][1] << " " << this->indexMap[p][2] << " " << this->indexMap[p][3] << " " << this->indexMap[p][4] << " " << this->spEnergy[p] << std::endl;
	}
} // end printBasis

void InfMatterSPBasis::deallocate(){

	for(std::size_t i = 0; i < this->nSpstates; i++){
		delete [] this->indexMap[i];
	}
	delete [] this->indexMap;
	delete [] this->spEnergy;
	delete [] this->chanValue;
	delete [] this->chanDims;
	delete [] this->chanMaps;
	delete [] this->chanModValue;
	delete [] this->chanModDims;
	delete [] this->chanModMaps;
} // end deallocate


double InfMatterSPBasis::calcTBME(std::size_t p, std::size_t q, std::size_t r, std::size_t s){
	return this->calc_TBME_not_antisym(p,q,r,s) -
		this->calc_TBME_not_antisym(p,q,s,r);
}


double InfMatterSPBasis::calc_TBME_not_antisym(std::size_t p, std::size_t q, std::size_t r, std::size_t s){
	double vout = 0.;
	int *qi = this->indexMap[p];
	int *qj = this->indexMap[q];
	int *qk = this->indexMap[r];
	int *ql = this->indexMap[s];

	int initialMomx, initialMomy, initialMomz;
	int finalMomx, finalMomy, finalMomz;
	// Momentum Conservation Checks.
	initialMomx = qi[0] + qj[0];
	initialMomy = qi[1] + qj[1];
	initialMomz = qi[2] + qj[2];
	finalMomx = qk[0] + ql[0];
	finalMomy = qk[1] + ql[1];
	finalMomz = qk[2] + ql[2];

	// maybe add a conservation of spin if here.
	if( initialMomx == finalMomx && initialMomy == finalMomy && initialMomz == finalMomz ){
		double L = pow(this->nParticles/this->density,1./3.);
		double V_R, V_T, V_S;
		double V_0R, V_0T, V_0S;
		double kappa_R, kappa_T, kappa_S;
		double ki[3];
		double kj[3];
		double kk[3];
		double kl[3];
		double relMomBra[3];
		double relMomKet[3];
		double relMomTransf[3];
		double qSquared, spinEx, isoSpinEx;
		double IsIt, PsIt, PsPt, IsPt;
		V_0R = 200; //MeV
		V_0T = 178; //MeV
		V_0S = 91.85; //MeV
		kappa_R = 1.487; //fm^-2
		kappa_T = 0.639; //fm^-2
		kappa_S = 0.465; //fm^-2

		qSquared = 0.;
		for( int i = 0; i < 3; i++){
			ki[i] = 2*M_PI*qi[i]/L;
			kj[i] = 2*M_PI*qj[i]/L;
			kk[i] = 2*M_PI*qk[i]/L;
			kl[i] = 2*M_PI*ql[i]/L;
			relMomBra[i] = 0.5*(ki[i] - kj[i]);
			relMomKet[i] = 0.5*(kk[i] - kl[i]);
			relMomTransf[i] = relMomBra[i] - relMomKet[i];
			qSquared += relMomTransf[i]*relMomTransf[i];
		}

		V_R = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5)*exp(-qSquared/(4*kappa_R));
		V_T = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5)*exp(-qSquared/(4*kappa_T));
		V_S = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5)*exp(-qSquared/(4*kappa_S));

		spinEx = spinExchangeMtxEle(qi[3],qj[3],qk[3],ql[3]);
		isoSpinEx = spinExchangeMtxEle(qi[4],qj[4],qk[4],ql[4]);

		// 4 terms, IsIt, PsIt, PsPt, IsPt
		// identity spin, identity isospin
		IsIt = kron_del(qi[3],qk[3])*kron_del(qj[3],ql[3])*kron_del(qi[4],qk[4])*kron_del(qj[4],ql[4]);
		// Exchange spin, identity isospin
		PsIt = spinEx*kron_del(qi[4],qk[4])*kron_del(qj[4],ql[4]);
		// Exchange spin, Exchange isospin
		PsPt = spinEx*isoSpinEx;
		// identity spin, Exchange isospin
		IsPt = kron_del(qi[3],qk[3])*kron_del(qj[3],ql[3])*isoSpinEx;

		vout = 0.5*(V_R + 0.5*V_T + 0.5*V_S)*IsIt
			+ 0.25*(V_T - V_S)*PsIt
			- 0.5*(V_R + 0.5*V_T + 0.5*V_S)*PsPt
			- 0.25*(V_T - V_S)*IsPt;

	} // end Momentum Conversation Check

	return vout;
}


int InfMatterSPBasis::spinExchangeMtxEle(int i, int j, int k, int l){
	if( i == l && j == k ){
		return 1;
	} else {
		return 0;
	}
} // end spinEx



int InfMatterSPBasis::kron_del(int i, int j){
	if(i != j){
		return 0;
	}
	return 1;
} // end kron_del


std::size_t InfMatterSPBasis::TBchanIndexFunction(std::size_t p, std::size_t q){
	int chanNx = this->indexMap[p][0] + this->indexMap[q][0];
	int chanNy = this->indexMap[p][1] + this->indexMap[q][1];
	int chanNz = this->indexMap[p][2] + this->indexMap[q][2];
	int chanSz = this->indexMap[p][3] + this->indexMap[q][3];
	int chanTz = this->indexMap[p][4] + this->indexMap[q][4];

	std::size_t index;
	std::size_t TBnMax = 4*this->nMax + 1;
	std::size_t TBsMax = 3;
	std::size_t TBtMax = 2+this->chanTzMax/2;// 1 or 3 pnm vs snm

	chanNx += TBnMax/2;
	chanNy += TBnMax/2;
	chanNz += TBnMax/2;
	chanSz = (chanSz+2)/2;
	chanTz = (chanTz+2)/2;  // 0 or 0,1,2
	index = ((std::size_t)chanNx)*TBnMax*TBnMax*TBsMax*TBtMax
		+ ((std::size_t)chanNy)*TBnMax*TBsMax*TBtMax
		+ ((std::size_t)chanNz)*TBsMax*TBtMax
		+ ((std::size_t)chanSz)*TBtMax
		+ ((std::size_t)chanTz);
	//return index;
	return this->inverseChannelIndices[index];
} // end indexMap


std::size_t InfMatterSPBasis::TBmodChanIndexFunction(std::size_t p, std::size_t q){
	int modChanNx = this->indexMap[p][0] - this->indexMap[q][0];
	int modChanNy = this->indexMap[p][1] - this->indexMap[q][1];
	int modChanNz = this->indexMap[p][2] - this->indexMap[q][2];
	int modChanSz = this->indexMap[p][3] - this->indexMap[q][3];
	int modChanTz = this->indexMap[p][4] - this->indexMap[q][4];

	std::size_t index;
	std::size_t TBnMax = 4*this->nMax + 1;
	std::size_t TBsMax = 3;
	std::size_t TBtMax = 1 + this->modChanTzMax; // 1 or 3

	modChanNx += TBnMax/2;
	modChanNy += TBnMax/2;
	modChanNz += TBnMax/2;
	modChanSz = (modChanSz+2)/2;
	modChanTz = (modChanTz+this->modChanTzMax)/2; // 0 or 0,1,2
	index = ((std::size_t)modChanNx)*TBnMax*TBnMax*TBsMax*TBtMax
		+ ((std::size_t)modChanNy)*TBnMax*TBsMax*TBtMax
		+ ((std::size_t)modChanNz)*TBsMax*TBtMax
		+ ((std::size_t)modChanSz)*TBtMax
		+ ((std::size_t)modChanTz);
	// inverseChannelIndices is the indentity if a load
	// balancer has not yet been ran.
	return this->inverseModChannelIndices[index];
} // end TBmodChanIndexFunction

// take in 5 quantum numbers, maps this to a unique SP label.
// THEN this label is a key to to hashed sp state index
std::size_t InfMatterSPBasis::keyGen_spHash(int nx, int ny, int nz, int sz, int tz){
	std::size_t nRange = 2*this->nMax + 1;
	std::size_t sRange = 2;
	std::size_t tzRange = this->tzMax;
	std::size_t index;
	// want all inputs to start at index 0
	nx += nRange/2;
	ny += nRange/2;
	nz += nRange/2;
	sz = (sz + 1)/2; // {-1,1} -> {0,1}
tz = (tz + tzRange)/2; // {-1} -> {0} or {-1,1} -> {0,1}
index = ((std::size_t)nx)*nRange*nRange*sRange*tzRange +
	((std::size_t)ny)*nRange*sRange*tzRange +
	((std::size_t)nz)*sRange*tzRange +
	((std::size_t)sz)*tzRange +
	((std::size_t)tz);
return index;
} // end keyGen_spHash

// Finds p given q^-1, r ,s
// returns -1 if p is not in the sp basis
std::size_t InfMatterSPBasis::spIndex_from3Body(std::size_t q_inv, std::size_t r, std::size_t s){
	int nx, ny, nz, sz, tz;
	nx = this->indexMap[r][0] + this->indexMap[s][0] - this->indexMap[q_inv][0];
	ny = this->indexMap[r][1] + this->indexMap[s][1] - this->indexMap[q_inv][1];
	nz = this->indexMap[r][2] + this->indexMap[s][2] - this->indexMap[q_inv][2];
	sz = this->indexMap[r][3] + this->indexMap[s][3] - this->indexMap[q_inv][3];
	tz = this->indexMap[r][4] + this->indexMap[s][4] - this->indexMap[q_inv][4];

	std::size_t key = this->keyGen_spHash(nx,ny,nz,sz,tz);

	return this->inverseSPChannelIndices[ sp_hash.at(key) ];
} // end spIndex_from3Body

int InfMatterSPBasis::spIndexExists_from3Body(std::size_t q_inv, std::size_t r, std::size_t s){
	int nx, ny, nz, sz, tz;
	nx = this->indexMap[r][0] + this->indexMap[s][0] - this->indexMap[q_inv][0];
	ny = this->indexMap[r][1] + this->indexMap[s][1] - this->indexMap[q_inv][1];
	nz = this->indexMap[r][2] + this->indexMap[s][2] - this->indexMap[q_inv][2];
	sz = this->indexMap[r][3] + this->indexMap[s][3] - this->indexMap[q_inv][3];
	tz = this->indexMap[r][4] + this->indexMap[s][4] - this->indexMap[q_inv][4];

	std::size_t key = this->keyGen_spHash(nx,ny,nz,sz,tz);

	if( sp_hash.count(key) == 0 ){
		return 0;
	}

	// return (this->sp_hash.count(key) != 0);

	return ( TBchanIndexFunction(r,s) == TBchanIndexFunction(q_inv,sp_hash[key]) );
} // end spIndex_from3Body
