/* 
 * File:   System.cpp
 * Author: chrishir
 * 
 * Created on 6. desember 2011, 14:08
 */

#include <vector>

#include "System.h"

using namespace std;
using namespace toffyrn::libUNK;
using namespace arma;

System::System()
{
    basis = NULL;
}

System::~System()
{
    if (basis != NULL)
        delete basis;
}

System::System(const System& sys)
{
    //Copy the basis
    this->basis = NULL;
    if (sys.basis != NULL)
        this->basis = sys.basis->clone();

    //Copy all other members.
    this->f_hh = sys.f_hh;
    this->f_ph = sys.f_ph;
    this->f_pp = sys.f_pp;
    this->v_hhhh = sys.v_hhhh;
    this->v_phhh = sys.v_phhh;
    this->v_phph = sys.v_phph;
    this->v_pphh = sys.v_pphh;
    this->v_ppph = sys.v_ppph;
    this->v_pppp = sys.v_pppp;
}

System& System::operator =(const System& sys)
{
    //Check for selfassignments
    if (this == &sys)
        return *this;

    //TODO: Should we not delete basis before setting it to NULL?
    //Copy the basis
    this->basis = NULL;
    if (sys.basis != NULL)
        this->basis = sys.basis->clone();

    //Copy all other members.
    this->f_hh = sys.f_hh;
    this->f_ph = sys.f_ph;
    this->f_pp = sys.f_pp;
    this->v_hhhh = sys.v_hhhh;
    this->v_phhh = sys.v_phhh;
    this->v_phph = sys.v_phph;
    this->v_pphh = sys.v_pphh;
    this->v_ppph = sys.v_ppph;
    this->v_pppp = sys.v_pppp;
}

void System::fillMatrixElements()
{
    //dimensions
    int nH = basis->get_nH();
    int nP = basis->get_nP();
    int nT = nH + nP;
    //Mappings
    vector<uvec> const * mapHH = basis->get_map_lmdMU_lm();
    vector<uvec> const * mapPH = basis->get_map_lmdNU_dl();
    vector<uvec> const * mapPP = basis->get_map_lmdXI_de();

    //Filling v_hhhh -> <lm||no> -> <mu||mu'>(lmd)
    v_hhhh.clear();
    for (size_t lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        size_t dimMU = mapHH->at(lmd).size();
        mat v_mumu(dimMU, dimMU);

        for (size_t mu1 = 0; mu1 < dimMU; mu1++)
            for (size_t mu2 = 0; mu2 < dimMU; mu2++)
            {
                int lm = mapHH->at(lmd)(mu1);
                int l = lm % nH;
                int m = lm / nH;
                int no = mapHH->at(lmd)(mu2);
                int n = no % nH;
                int o = no / nH;

                v_mumu(mu1, mu2) = v_elem(l, m, n, o);
            }

        v_hhhh.push_back(v_mumu);
    }

    //Filling v_phhh -> <dl||mn> -> <nu||mu>(lmd)
    v_phhh.clear();
    for (size_t lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        size_t dimMU = mapHH->at(lmd).size();
        size_t dimNU = mapPH->at(lmd).size();
        mat v_numu(dimNU, dimMU);

        for (size_t nu = 0; nu < dimNU; nu++)
            for (size_t mu = 0; mu < dimMU; mu++)
            {
                int dl = mapPH->at(lmd)(nu);
                int d = dl % nP;
                int l = dl / nP;
                int mn = mapHH->at(lmd)(mu);
                int m = mn % nH;
                int n = mn / nH;

                v_numu(nu, mu) = v_elem(d + nH, l, m, n);
            }

        v_phhh.push_back(v_numu);
    }

    //Filling v_pphh -> <de||lm> -> <xi||mu>(lmd)
    v_pphh.clear();
    for (size_t lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        size_t dimXI = mapPP->at(lmd).size();
        size_t dimMU = mapHH->at(lmd).size();
        mat v_ximu(dimXI, dimMU);

        for (size_t xi = 0; xi < dimXI; xi++)
            for (size_t mu = 0; mu < dimMU; mu++)
            {
                int de = mapPP->at(lmd)(xi);
                int d = de % nP;
                int e = de / nP;
                int lm = mapHH->at(lmd)(mu);
                int l = lm % nH;
                int m = lm / nH;

                v_ximu(xi, mu) = v_elem(d + nH, e + nH, l, m);
            }

        v_pphh.push_back(v_ximu);
    }


    //Filling v_phph -> <dl||em> -> <nu||nu'>(lmd)
    v_phph.clear();
    for (size_t lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        size_t dimNU = mapPH->at(lmd).size();
        mat v_nunu(dimNU, dimNU);

        for (size_t nu1 = 0; nu1 < dimNU; nu1++)
            for (size_t nu2 = 0; nu2 < dimNU; nu2++)
            {
                int dl = mapPH->at(lmd)(nu1);
                int d = dl % nP;
                int l = dl / nP;
                int em = mapPH->at(lmd)(nu2);
                int e = em % nP;
                int m = em / nP;

                v_nunu(nu1, nu2) = v_elem(d + nH, l, e + nH, m);
            }

        v_phph.push_back(v_nunu);
    }

    //Filling v_ppph -> <de||fl> -> <xi||nu>(lmd)
    v_ppph.clear();
    for (size_t lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        size_t dimXI = mapPP->at(lmd).size();
        size_t dimNU = mapPH->at(lmd).size();
        mat v_xinu(dimXI, dimNU);

        for (size_t xi = 0; xi < dimXI; xi++)
            for (size_t nu = 0; nu < dimNU; nu++)
            {
                int de = mapPP->at(lmd)(xi);
                int d = de % nP;
                int e = de / nP;
                int fl = mapPH->at(lmd)(nu);
                int f = fl % nP;
                int l = fl / nP;

                v_xinu(xi, nu) = v_elem(d + nH, e + nH, f + nH, l);
            }

        v_ppph.push_back(v_xinu);
    }

    //Filling v_pppp -> <de||fg> -> <xi||xi'>(lmd)
    v_pppp.clear();
    for (size_t lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        size_t dimXI = mapPP->at(lmd).size();
        mat v_xixi(dimXI, dimXI);

        for (size_t xi1 = 0; xi1 < dimXI; xi1++)
            for (size_t xi2 = 0; xi2 < dimXI; xi2++)
            {
                int de = mapPP->at(lmd)(xi1);
                int d = de % nP;
                int e = de / nP;
                int fg = mapPP->at(lmd)(xi2);
                int f = fg % nP;
                int g = fg / nP;

                v_xixi(xi1, xi2) = v_elem(d + nH, e + nH, f + nH, g + nH);
            }

        v_pppp.push_back(v_xixi);
    }

    //Filling f_matrices
    f_hh = mat(nH, nH);
    for (size_t l = 0; l < nH; l++)
        for (size_t m = 0; m < nH; m++)
            f_hh(l, m) = f_elem(l, m);
    f_ph = mat(nP, nH);
    for (size_t d = 0; d < nP; d++)
        for (size_t l = 0; l < nH; l++)
            f_ph(d, l) = f_elem(d + nH, l);
    f_pp = mat(nP, nP);
    for (size_t d = 0; d < nP; d++)
        for (size_t e = 0; e < nP; e++)
            f_pp(d, e) = f_elem(d + nH, e + nH);
}

