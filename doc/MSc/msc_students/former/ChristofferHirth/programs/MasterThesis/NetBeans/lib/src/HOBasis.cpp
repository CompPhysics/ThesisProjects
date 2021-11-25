/**
 * @file   HOBasis.cpp
 * @author toffyrn
 * 
 * Created on 29. desember 2011, 15:23
 */

#include "HOBasis.h"

using namespace toffyrn::libUNK;
using namespace arma;
using namespace std;

HOBasis::HOBasis(int filledR, int totalR)
{
    //Store the number of shells
    shells_filled = filledR;
    shells_tot = totalR;

    //Store number of states
    this->nH = filledR * (filledR + 1);
    int nTOT = totalR * (totalR + 1);
    this->nP = nTOT - nH;

    //Initialize statemap
    initStateMap(nTOT);
    //Initialize mappings
    initMappings();
}

void HOBasis::initStateMap(int nStates)
{
    //Allocate size
    statemapMat = zeros<imat > (3, nStates);

    //p -> [n,m,ms]
    //
    //up -> 0, do -> 1
    //
    //0 -> [0,0,up]
    //1 -> [0,0,do]
    //2 -> [0,-1,up]
    //3 -> [0,-1,do]
    //4 -> [0,+1,up]
    //5 -> [0,+1,do]
    //6 -> [0,-2,up]
    //7 -> [0,-2,do]
    //8 -> [1,0,up]
    //9 -> [1,0,do]
    //10-> [0,2,up]
    //11-> [0,2,do]
    // ...
    int R = 0;
    int m = 0;
    int ms = 0;
    for (int state = 0; state < nStates; state++)
    {
        int n = (R - abs(m)) / 2;
        //Save the info
        ivec3 res;
        res << n << m << ms;
        statemapMat.col(state) = res;

        //Update for next state
        if (ms == 1)
        {
            ms = 0;
            if (m == R)
            {
                R++;
                m = -R;
            } else
                m += 2;
        } else
            ms++;
    }
}

void HOBasis::initMappings()
{
    //The maximum M value (m1 + m2)
    int maxM = 2 * (shells_tot - 1);
    //The total number of M values (0,...,totM-1) -> (-maxM,...,0,...,maxM)
    int totM = 2 * maxM + 1;
    //The total number of Ms values (0,1,2) -> (1,0,-1) 
    int totMs = 3;
    //Number lmd values (0,1,...,dimLmd-1) where lmd = M + Ms * totM
    int dimLmd = totMs * totM;
    //total number of states.
    int nTOT = nH + nP;

    //Allocate/clear mapping containers
    map_pq_lmdPI = zeros<umat > (2, nTOT * nTOT);
    map_lmdPI_pq.clear();
    map_de_lmdXI = zeros<umat > (2, nP * nP);
    map_lmdXI_de.clear();
    map_dl_lmdNU = zeros<umat > (2, nP * nH);
    map_lmdNU_dl.clear();
    map_lm_lmdMU = zeros<umat > (2, nH * nH);
    map_lmdMU_lm.clear();

    //Create mappings
    //(size is unknown, therefore using stl Vectors temporary.)
    vector < vector<int > > stlVec_pq(dimLmd);
    vector < vector<int > > stlVec_de(dimLmd);
    vector < vector<int > > stlVec_dl(dimLmd);
    vector < vector<int > > stlVec_lm(dimLmd);

    //Loop over all tp states (p,q)
    for (int pq = 0; pq < nTOT * nTOT; pq++)
    {
        int p = pq % nTOT;
        int q = pq / nTOT;

        ivec3 p_nmm_s = stateMap(p);
        ivec3 q_nmm_s = stateMap(q);
        int M = p_nmm_s(1) + q_nmm_s(1) + maxM; //0,1,...,totM
        int Ms = p_nmm_s(2) + q_nmm_s(2); //0,1,2
        int lmd = M + Ms * totM;

        //Fill the general mapping
        stlVec_pq.at(lmd).push_back(pq);
        int pi = stlVec_pq.at(lmd).size() - 1;
        map_pq_lmdPI(0, pq) = lmd;
        map_pq_lmdPI(1, pq) = pi;

        //Particle particle state (d,e)
        if (p >= nH && q >= nH)
        {
            int d = p - nH;
            int e = q - nH;
            int de = d + e *nP;
            stlVec_de.at(lmd).push_back(de);
            int xi = stlVec_de.at(lmd).size() - 1;
            map_de_lmdXI(0, de) = lmd;
            map_de_lmdXI(1, de) = xi;
        }

        //Particle hole state (d,l)
        if (p >= nH && q < nH)
        {
            int d = p - nH;
            int l = q;
            int dl = d + l * nP;
            stlVec_dl.at(lmd).push_back(dl);
            int nu = stlVec_dl.at(lmd).size() - 1;
            map_dl_lmdNU(0, dl) = lmd;
            map_dl_lmdNU(1, dl) = nu;
        }

        //Hole hole state (l,m)
        if (p < nH && q < nH)
        {
            int l = p;
            int m = q;
            int lm = l + m * nH;
            stlVec_lm.at(lmd).push_back(lm);
            int mu = stlVec_lm.at(lmd).size() - 1;
            map_lm_lmdMU(0, lm) = lmd;
            map_lm_lmdMU(1, lm) = mu;
        }
    }

    //Copy mappings from stl containers to armadillo
    for (int lmd = 0; lmd < dimLmd; lmd++)
    {
        map_lmdPI_pq.push_back(conv_to<uvec>::from(stlVec_pq.at(lmd)));
        map_lmdMU_lm.push_back(conv_to<uvec>::from(stlVec_lm.at(lmd)));
        map_lmdNU_dl.push_back(conv_to<uvec>::from(stlVec_dl.at(lmd)));
        map_lmdXI_de.push_back(conv_to<uvec>::from(stlVec_de.at(lmd)));
    }

}

ivec HOBasis::stateMap(int p) const
{
    return statemapMat.col(p);
}

int HOBasis::lambda_neg(int p, int q) const
{
    //The maximum M value (m1 + m2)
    int maxM = 2 * (shells_tot - 1);
    //The total number of M values (0,...,totM-1) -> (-maxM,...,0,...,maxM)
    int totM = 2 * maxM + 1;

    ivec3 p_nmm_s = stateMap(p);
    ivec3 q_nmm_s = stateMap(q);

    int M = p_nmm_s(1) - q_nmm_s(1) + maxM; //0,1,...,totM
    int Ms = p_nmm_s(2) - q_nmm_s(2) + 1; //0,1,2
    int lmd = M + Ms * totM;

    return lmd;
}

int HOBasis::lambda_1p(int p) const
{
    //The maximum M value for one particle
    int maxM = shells_tot - 1;
    //The total number of M values (0,...,totM-1) -> (-maxM,...,0,...,maxM)
    int totM = 2 * maxM + 1;

    ivec3 p_nmm_s = stateMap(p);

    int M = p_nmm_s(1) + maxM; //0,1,...,totM
    int Ms = p_nmm_s(2); //0,1

    int lmd = M + Ms * totM;

    return lmd;
}

int HOBasis::dim_lmd_1p() const
{
    //The maximum M value for one particle
    int maxM = shells_tot - 1;
    //The total number of M values (0,...,totM-1) -> (-maxM,...,0,...,maxM)
    int totM = 2 * maxM + 1;

    //One particle can have two spin states (up or down), and totM angular momentum states.
    return 2 * totM;
}

