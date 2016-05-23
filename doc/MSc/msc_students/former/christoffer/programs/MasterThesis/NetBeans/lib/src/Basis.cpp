/* 
 * File:   Basis.cpp
 * Author: chrishir
 * 
 * Created on 7. desember 2011, 10:49
 */

#include "Basis.h"

using namespace arma;
using namespace toffyrn::libUNK;

Basis::Basis()
{
    this->nH = 0;
    this->nP = 0;
}

Basis::Basis(int nH, int nP)
{
    //Store number of states.
    this->nH = nH;
    this->nP = nP;

    //Initialize default mapping
    initMappings();
}

Basis::~Basis()
{
}

void Basis::initMappings()
{
    //General mapping
    //pi = p+q*nTot
    int nTot = nH + nP;
    map_lmdPI_pq.clear();
    map_lmdPI_pq.push_back(linspace<uvec > (0, nTot * nTot - 1, nTot * nTot));
    map_pq_lmdPI = umat(2, nTot * nTot);
    map_pq_lmdPI.row(0) = zeros<urowvec > (nTot * nTot);
    map_pq_lmdPI.row(1) = linspace<urowvec > (0, nTot * nTot - 1, nTot * nTot);

    //holehole mapping
    //mu = l+m*nH
    map_lmdMU_lm.clear();
    map_lmdMU_lm.push_back(linspace<uvec > (0, nH * nH - 1, nH * nH));
    map_lm_lmdMU = umat(2, nH * nH);
    map_lm_lmdMU.row(0) = zeros<urowvec > (nH * nH);
    map_lm_lmdMU.row(1) = linspace<urowvec > (0, nH * nH - 1, nH * nH);

    //Particle-hole mapping
    //nu = d + l* nP
    map_lmdNU_dl.clear();
    map_lmdNU_dl.push_back(linspace<uvec > (0, nH * nP - 1, nH * nP));
    map_dl_lmdNU = umat(2, nP * nH);
    map_dl_lmdNU.row(0) = zeros<urowvec > (nP * nH);
    map_dl_lmdNU.row(1) = linspace<urowvec > (0, nH * nP - 1, nH * nP);

    //Particle-particle mapping
    //xi = d + e *nP
    map_lmdXI_de.clear();
    map_lmdXI_de.push_back(linspace<uvec > (0, nP * nP - 1, nP * nP));
    map_de_lmdXI = umat(2, nP * nP);
    map_de_lmdXI.row(0) = zeros<urowvec > (nP * nP);
    map_de_lmdXI.row(1) = linspace<urowvec > (0, nP * nP - 1, nP * nP);
}

