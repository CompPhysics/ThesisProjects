/* 
 * File:   HFbasis.cpp
 * Author: toffyrn
 * 
 * Created on 16. februar 2012, 17:21
 */

#include "HFbasis.h"
#include "System.h"

using namespace toffyrn::libUNK;
using namespace arma;
using namespace std;

HFbasis::HFbasis(System const * originalSys, arma::mat const &Coeff)
{
    this->basis = originalSys->get_basis()->clone();
    this->f_hh = *originalSys->get_f_hh();
    this->f_ph = *originalSys->get_f_ph();
    this->f_pp = *originalSys->get_f_pp();
    this->v_hhhh = *originalSys->get_v_hhhh();
    this->v_phhh = *originalSys->get_v_phhh();
    this->v_pphh = *originalSys->get_v_pphh();
    this->v_phph = *originalSys->get_v_phph();
    this->v_ppph = *originalSys->get_v_ppph();
    this->v_pppp = *originalSys->get_v_pppp();
    transformElements(Coeff);
}

HFbasis::~HFbasis()
{
}

double HFbasis::f_elem(std::size_t p, std::size_t q) const
{
    int nH = basis->get_nH();
    int nP = basis->get_nP();

    if (p < nH && q < nH)
        return f_hh(p, q);
    else if (p >= nH && q < nH)
        return f_ph(p - nH, q);
    else if (p < nH && q >= nH)
        return f_ph(q - nH, p);
    else if (p >= nH && q >= nH)
        return f_pp(p - nH, q - nH);
    else
        throw std::string("Unreachable element in f_elem.");
}

double HFbasis::v_elem(
        std::size_t p, std::size_t q,
        std::size_t r, std::size_t s) const
{
    //Needed mappings & constants
    int nH = basis->get_nH();
    int nP = basis->get_nP();
    umat const * mapHH = basis->get_map_lm_lmdMU();
    umat const * mapPH = basis->get_map_dl_lmdNU();
    umat const * mapPP = basis->get_map_de_lmdXI();

    //Element value and phasefactor
    double elem = 0;
    double phase = 1.0;

    //Is there 0, 1 or 2 particle states in bra or ket?
    int pStatesLeft = 0;
    if (p >= nH)
        pStatesLeft++;
    if (q >= nH)
        pStatesLeft++;
    int pStatesRight = 0;
    if (r >= nH)
        pStatesRight++;
    if (s >= nH)
        pStatesRight++;

    //Gather particle states to the left.
    if (pStatesRight > pStatesLeft)
    {
        //Swap p,q <-> r,s
        int tmp = p;
        p = r;
        r = tmp;
        tmp = q;
        q = s;
        s = tmp;
        //Left <-> Right
        tmp = pStatesRight;
        pStatesRight = pStatesLeft;
        pStatesLeft = tmp;
    }
    if (q > p)
    {
        int tmp = p;
        p = q;
        q = tmp;
        phase *= -1;
    }
    if (s > r)
    {
        int tmp = r;
        r = s;
        s = tmp;
        phase *= -1;
    }



    if (pStatesLeft == 0 && pStatesRight == 0)
    { //hhhh
        int pq = p + q * nH;
        int lmd_pq = (*mapHH)(0, pq);
        int mu_pq = (*mapHH)(1, pq);
        int rs = r + s * nH;
        int lmd_rs = (*mapHH)(0, rs);
        int mu_rs = (*mapHH)(1, rs);
        if (lmd_pq == lmd_rs)
            elem = v_hhhh.at(lmd_pq)(mu_pq, mu_rs);
    } else if (pStatesLeft == 1 && pStatesRight == 0)
    { //phhh
        int pq = (p - nH) + q * nP;
        int lmd_pq = (*mapPH)(0, pq);
        int nu_pq = (*mapPH)(1, pq);
        int rs = r + s * nH;
        int lmd_rs = (*mapHH)(0, rs);
        int mu_rs = (*mapHH)(1, rs);
        if (lmd_pq == lmd_rs)
            elem = v_phhh.at(lmd_pq)(nu_pq, mu_rs);
    } else if (pStatesLeft == 1 && pStatesRight == 1)
    { //phph
        int pq = (p - nH) + q * nP;
        int lmd_pq = (*mapPH)(0, pq);
        int nu_pq = (*mapPH)(1, pq);
        int rs = (r - nH) + s * nP;
        int lmd_rs = (*mapPH)(0, rs);
        int nu_rs = (*mapPH)(1, rs);
        if (lmd_pq == lmd_rs)
            elem = v_phph.at(lmd_pq)(nu_pq, nu_rs);
    } else if (pStatesLeft == 2 && pStatesRight == 0)
    { //pphh
        int pq = (p - nH) + (q - nH) * nP;
        int lmd_pq = (*mapPP)(0, pq);
        int xi_pq = (*mapPP)(1, pq);
        int rs = r + s * nH;
        int lmd_rs = (*mapHH)(0, rs);
        int mu_rs = (*mapHH)(1, rs);
        if (lmd_pq == lmd_rs)
            elem = v_pphh.at(lmd_pq)(xi_pq, mu_rs);
    } else if (pStatesLeft == 2 && pStatesRight == 1)
    { //ppph
        int pq = (p - nH) + (q - nH) * nP;
        int lmd_pq = (*mapPP)(0, pq);
        int xi_pq = (*mapPP)(1, pq);
        int rs = (r - nH) + s * nP;
        int lmd_rs = (*mapPH)(0, rs);
        int nu_rs = (*mapPH)(1, rs);
        if (lmd_pq == lmd_rs)
            elem = v_ppph.at(lmd_pq)(xi_pq, nu_rs);
    } else if (pStatesLeft == 2 && pStatesRight == 2)
    { //pppp
        int pq = (p - nH) + (q - nH) * nP;
        int lmd_pq = (*mapPP)(0, pq);
        int xi_pq = (*mapPP)(1, pq);
        int rs = (r - nH) + (s - nH) * nP;
        int lmd_rs = (*mapPP)(0, rs);
        int xi_rs = (*mapPP)(1, rs);
        if (lmd_pq == lmd_rs)
            elem = v_pppp.at(lmd_pq)(xi_pq, xi_rs);
    }

    return phase * elem;
}

double HFbasis::transformOneElem(int p, int q, int r, int s, mat const &C)
{
    int nP = basis->get_nP();
    int nH = basis->get_nH();
    int nTot = nH + nP;
    vector<uvec> const * ppMap = basis->get_map_lmdXI_de();
    vector<uvec> const * phMap = basis->get_map_lmdNU_dl();
    vector<uvec> const * hhMap = basis->get_map_lmdMU_lm();
    int dimLMD = basis->dim_lmd_2p();

    double elm = 0;
    for (int lmd = 0; lmd < dimLMD; lmd++)
    {
        int dimMU = hhMap->at(lmd).size();
        int dimNU = phMap->at(lmd).size();
        int dimXI = ppMap->at(lmd).size();

        //hhhh
        for (int mu_ab = 0; mu_ab < dimMU; mu_ab++)
            for (int mu_gd = 0; mu_gd < dimMU; mu_gd++)
            {
                int ab = hhMap->at(lmd)(mu_ab);
                int alpha = ab % nH;
                int beta = ab / nH;
                int gd = hhMap->at(lmd)(mu_gd);
                int gamma = gd % nH;
                int delta = gd / nH;
                double cur_elem = v_hhhh.at(lmd)(mu_ab, mu_gd);
                elm += C(p, alpha) * C(q, beta) * C(r, gamma) * C(s, delta) * cur_elem;
            }
        //phhh
        for (int nu_ab = 0; nu_ab < dimNU; nu_ab++)
            for (int mu_gd = 0; mu_gd < dimMU; mu_gd++)
            {
                int ab = phMap->at(lmd)(nu_ab);
                int alpha = ab % nP + nH;
                int beta = ab / nP;
                int gd = hhMap->at(lmd)(mu_gd);
                int gamma = gd % nH;
                int delta = gd / nH;
                double cur_elem = v_phhh.at(lmd)(nu_ab, mu_gd);
                elm += C(p, alpha) * C(q, beta) * C(r, gamma) * C(s, delta) * cur_elem;
                elm -= C(p, beta) * C(q, alpha) * C(r, gamma) * C(s, delta) * cur_elem;
                elm += C(p, gamma) * C(q, delta) * C(r, alpha) * C(s, beta) * cur_elem;
                elm -= C(p, gamma) * C(q, delta) * C(r, beta) * C(s, alpha) * cur_elem;
            }
        //pphh
        for (int xi_ab = 0; xi_ab < dimXI; xi_ab++)
            for (int mu_gd = 0; mu_gd < dimMU; mu_gd++)
            {
                int ab = ppMap->at(lmd)(xi_ab);
                int alpha = ab % nP + nH;
                int beta = ab / nP + nH;
                int gd = hhMap->at(lmd)(mu_gd);
                int gamma = gd % nH;
                int delta = gd / nH;
                double cur_elem = v_pphh.at(lmd)(xi_ab, mu_gd);
                elm += C(p, alpha) * C(q, beta) * C(r, gamma) * C(s, delta) * cur_elem;
                elm += C(p, gamma) * C(q, delta) * C(r, alpha) * C(s, beta) * cur_elem;
            }
        //phph
        for (int nu_ab = 0; nu_ab < dimNU; nu_ab++)
            for (int nu_gd = 0; nu_gd < dimNU; nu_gd++)
            {
                int ab = phMap->at(lmd)(nu_ab);
                int alpha = ab % nP + nH;
                int beta = ab / nP;
                int gd = phMap->at(lmd)(nu_gd);
                int gamma = gd % nP + nH;
                int delta = gd / nP;
                double cur_elem = v_phph.at(lmd)(nu_ab, nu_gd);
                elm += C(p, alpha) * C(q, beta) * C(r, gamma) * C(s, delta) * cur_elem;
                elm += C(p, beta) * C(q, alpha) * C(r, delta) * C(s, gamma) * cur_elem;
                elm -= C(p, beta) * C(q, alpha) * C(r, gamma) * C(s, delta) * cur_elem;
                elm -= C(p, alpha) * C(q, beta) * C(r, delta) * C(s, gamma) * cur_elem;
            }
        //ppph
        for (int xi_ab = 0; xi_ab < dimXI; xi_ab++)
            for (int nu_gd = 0; nu_gd < dimNU; nu_gd++)
            {
                int ab = ppMap->at(lmd)(xi_ab);
                int alpha = ab % nP + nH;
                int beta = ab / nP + nH;
                int gd = phMap->at(lmd)(nu_gd);
                int gamma = gd % nP + nH;
                int delta = gd / nP;
                double cur_elem = v_ppph.at(lmd)(xi_ab, nu_gd);
                elm += C(p, alpha) * C(q, beta) * C(r, gamma) * C(s, delta) * cur_elem;
                elm -= C(p, alpha) * C(q, beta) * C(r, delta) * C(s, gamma) * cur_elem;
                elm += C(p, gamma) * C(q, delta) * C(r, alpha) * C(s, beta) * cur_elem;
                elm -= C(p, delta) * C(q, gamma) * C(r, alpha) * C(s, beta) * cur_elem;
            }
        //pppp
        for (int xi_ab = 0; xi_ab < dimXI; xi_ab++)
            for (int xi_gd = 0; xi_gd < dimXI; xi_gd++)
            {
                int ab = ppMap->at(lmd)(xi_ab);
                int alpha = ab % nP + nH;
                int beta = ab / nP + nH;
                int gd = ppMap->at(lmd)(xi_gd);
                int gamma = gd % nP + nH;
                int delta = gd / nP + nH;
                double cur_elem = v_pppp.at(lmd)(xi_ab, xi_gd);
                elm += C(p, alpha) * C(q, beta) * C(r, gamma) * C(s, delta) * cur_elem;
            }
    }

    return elm;
}

void HFbasis::transformElements(mat const &C)
{
    int nH = basis->get_nH();
    int nP = basis->get_nP();
    int nTot = nH + nP;
    int dimLMD2 = basis->dim_lmd_2p();
    vector<uvec> const * mapHH = basis->get_map_lmdMU_lm();
    vector<uvec> const * mapPH = basis->get_map_lmdNU_dl();
    vector<uvec> const * mapPP = basis->get_map_lmdXI_de();

    //Get the singleparticle interactions.
    mat h0(nTot, nTot);
    span hSpan(0, nH - 1);
    span pSpan(nH, nTot - 1);
    h0(hSpan, hSpan) = f_hh;
    h0(pSpan, hSpan) = f_ph;
    h0(hSpan, pSpan) = trans(h0(pSpan, hSpan));
    h0(pSpan, pSpan) = f_pp;
    for (int p = 0; p < nTot; p++)
        for (int q = 0; q < nTot; q++)
        {
            double u_pq = 0;
            for (int i = 0; i < nH; i++)
                u_pq += v_elem(p, i, q, i);
            h0(p, q) -= u_pq;
        }
    //Transform h0
    h0 = C * h0 * C.t();

    //TODO: Transform 2particle interactions and verify results.
    cout << "2p elem\n";
    vector<mat> new_hhhh = v_hhhh;
    vector<mat> new_phhh = v_phhh;
    vector<mat> new_pphh = v_pphh;
    vector<mat> new_phph = v_phph;
    vector<mat> new_ppph = v_ppph;
    vector<mat> new_pppp = v_pppp;
    for (int lmd2 = 0; lmd2 < dimLMD2; lmd2++)
    {
        int dimMU = mapHH->at(lmd2).size();
        int dimNU = mapPH->at(lmd2).size();
        int dimXI = mapPP->at(lmd2).size();

        //hhhh
        for (int mu_ij = 0; mu_ij < dimMU; mu_ij++)
            for (int mu_kl = 0; mu_kl < dimMU; mu_kl++)
            {
                int ij = mapHH->at(lmd2)(mu_ij);
                int i = ij % nH;
                int j = ij / nH;
                int kl = mapHH->at(lmd2)(mu_kl);
                int k = kl % nH;
                int l = kl / nH;
                new_hhhh.at(lmd2)(mu_ij, mu_kl) = transformOneElem(i, j, k, l, C);
            }
        //phhh
        for (int nu_aj = 0; nu_aj < dimNU; nu_aj++)
            for (int mu_kl = 0; mu_kl < dimMU; mu_kl++)
            {
                int aj = mapPH->at(lmd2)(nu_aj);
                int a = aj % nP;
                int j = aj / nP;
                int kl = mapHH->at(lmd2)(mu_kl);
                int k = kl % nH;
                int l = kl / nH;
                new_phhh.at(lmd2)(nu_aj, mu_kl) = transformOneElem(a + nH, j, k, l, C);
            }
        //pphh
        for (int xi_ab = 0; xi_ab < dimXI; xi_ab++)
            for (int mu_kl = 0; mu_kl < dimMU; mu_kl++)
            {
                int ab = mapPP->at(lmd2)(xi_ab);
                int a = ab % nP;
                int b = ab / nP;
                int kl = mapHH->at(lmd2)(mu_kl);
                int k = kl % nH;
                int l = kl / nH;
                new_pphh.at(lmd2)(xi_ab, mu_kl) = transformOneElem(a + nH, b + nH, k, l, C);
            }
        //phph
        for (int nu_aj = 0; nu_aj < dimNU; nu_aj++)
            for (int nu_cl = 0; nu_cl < dimNU; nu_cl++)
            {
                int aj = mapPH->at(lmd2)(nu_aj);
                int a = aj % nP;
                int j = aj / nP;
                int cl = mapPH->at(lmd2)(nu_cl);
                int c = cl % nP;
                int l = cl / nP;
                new_phph.at(lmd2)(nu_aj, nu_cl) = transformOneElem(a + nH, j, c + nH, l, C);
            }
        //ppph
        for (int xi_ab = 0; xi_ab < dimXI; xi_ab++)
            for (int nu_cl = 0; nu_cl < dimNU; nu_cl++)
            {
                int ab = mapPP->at(lmd2)(xi_ab);
                int a = ab % nP;
                int b = ab / nP;
                int cl = mapPH->at(lmd2)(nu_cl);
                int c = cl % nP;
                int l = cl / nP;
                new_ppph.at(lmd2)(xi_ab, nu_cl) = transformOneElem(a + nH, b + nH, c + nH, l, C);
            }
        //pppp
        for (int xi_ab = 0; xi_ab < dimXI; xi_ab++)
            for (int xi_cd = 0; xi_cd < dimXI; xi_cd++)
            {
                int ab = mapPP->at(lmd2)(xi_ab);
                int a = ab % nP;
                int b = ab / nP;
                int cd = mapPP->at(lmd2)(xi_cd);
                int c = cd % nP;
                int d = cd / nP;
                new_pppp.at(lmd2)(xi_ab, xi_cd) = transformOneElem(a + nH, b + nH, c + nH, d + nH, C);
            }
    }
    cout << "EO2p elem\n";
    v_hhhh = new_hhhh;
    v_phhh = new_phhh;
    v_pphh = new_pphh;
    v_phph = new_phph;
    v_ppph = new_ppph;
    v_pppp = new_pppp;


    //Find f_hh
    for (int i = 0; i < nH; i++)
        for (int j = 0; j < nH; j++)
        {
            double u_ij = 0;
            for (int k = 0; k < nH; k++)
                u_ij += v_elem(i, k, j, k);
            f_hh(i, j) = h0(i, j) + u_ij;
        }
    //Find f_ph 
    for (int a = 0; a < nP; a++)
        for (int j = 0; j < nH; j++)
        {
            double u_aj = 0;
            for (int k = 0; k < nH; k++)
                u_aj += v_elem(a + nH, k, j, k);
            f_ph(a, j) = h0(a + nH, j) + u_aj;
        }
    //Find f_pp
    for (int a = 0; a < nP; a++)
        for (int b = 0; b < nP; b++)
        {
            double u_ab = 0;
            for (int k = 0; k < nH; k++)
                u_ab += v_elem(a + nH, k, b + nH, k);
            f_pp(a, b) = h0(a + nH, b + nH) + u_ab;
        }
    //Ensure matrices are symmetric //Some methods may diverge otherwise.
    f_hh = symmatl(f_hh);
    f_pp = symmatl(f_pp);

    return;
}


