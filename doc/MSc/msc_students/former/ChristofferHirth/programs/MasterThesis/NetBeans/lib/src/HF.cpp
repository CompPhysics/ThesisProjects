/**
 * @file   HF.cpp
 * @author toffyrn
 * 
 * Created on 31. januar 2012, 09:57
 */

#include <iomanip>

#include "HF.h"

using namespace toffyrn::libUNK;
using namespace arma;
using namespace std;

HF::HF(double precision, int max_iter)
{
    this->max_iter = max_iter;
    this->precision = precision;
}

HF::~HF()
{
}

void HF::set_system(System* sys)
{
    this->sys = sys;
}

pair<double, mat> HF::solve_ground_state_energy(int debug_flags, std::ostream& debug_stream)
{
    //SP states to use
    int numHOLEstates = sys->get_basis()->get_nH();
    int numPARTstates = sys->get_basis()->get_nP();
    int numSPstates = numPARTstates + numHOLEstates;

    //Get the singleparticle interactions.
    mat h0(numSPstates, numSPstates);
    span hSpan(0, numHOLEstates - 1);
    span pSpan(numHOLEstates, numSPstates - 1);
    h0(hSpan, hSpan) = *sys->get_f_hh();
    h0(pSpan, hSpan) = *sys->get_f_ph();
    h0(hSpan, pSpan) = trans(h0(pSpan, hSpan));
    h0(pSpan, pSpan) = *sys->get_f_pp();
    for (int p = 0; p < numSPstates; p++)
        for (int q = 0; q < numSPstates; q++)
        {
            double u_pq = 0;
            for (int i = 0; i < numHOLEstates; i++)
                u_pq += sys->v_elem(p, i, q, i);
            h0(p, q) -= u_pq;
        }

    //Setting up The coefficient matrix
    mat C = eye< mat > (numSPstates, numSPstates);
    //And derived types
    mat C_holeXall = C.submat(span(0, numHOLEstates - 1), span::all);
    mat C_inner = C_holeXall.t() * C_holeXall;

    //Store energy
    double E_old = energy(C_inner, h0);

    for (int iteration = 0; iteration < max_iter; iteration++)
    {
        //Build hartree fock matrix
        mat H_hf = hf_matrix(C_inner, h0);

        //Finding eigenvalues/vectors
        vec eigval;
        mat eigvec;
        eig_sym(eigval, eigvec, H_hf);

        //Sorting values/vectors
        sort_eigvecval(eigvec, eigval);

        C = trans(eigvec);

        C_holeXall = C.submat(span(0, numHOLEstates - 1), span::all);
        C_inner = C_holeXall.t() * C_holeXall;

        //E_ref    
        double E_new = energy(C_inner, h0);
        cout << std::setprecision(10) << E_new << endl;
        if (fabs(E_new - E_old) < precision)
        {
            E_old = E_new;
            break;
        }
        E_old = E_new;

    }

    //Pair up the energy and Coefficients.
    return make_pair<double, mat > (E_old, C);
}

void HF::sort_eigvecval(mat &eigvec, vec &eigval)
{
    for (int idx = 1; idx < eigval.n_elem; idx++)
    {
        for (int search = idx; search > 0; search--)
            if (eigval(search) < eigval(search - 1))
            {
                double buffer = eigval(search - 1);
                eigval(search - 1) = eigval(search);
                eigval(search) = buffer;

                vec buffer_vec = eigvec.col(search - 1);
                eigvec.col(search - 1) = eigvec.col(search);
                eigvec.col(search) = buffer_vec;
            }
    }

    return;
}

mat HF::hf_matrix(mat const &C_inner, mat const &h0)
{
    Basis const * basis = sys->get_basis();
    vector<uvec> const * ppMap = basis->get_map_lmdXI_de();
    vector<uvec> const * phMap = basis->get_map_lmdNU_dl();
    vector<uvec> const * hhMap = basis->get_map_lmdMU_lm();
    int dimLMD = basis->dim_lmd_2p();
    int nH = sys->get_basis()->get_nH();
    int nP = sys->get_basis()->get_nP();

    mat H_hf = h0;
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
                double elem = sys->get_v_hhhh()->at(lmd)(mu_ab, mu_gd);
                H_hf(alpha, gamma) += C_inner(beta, delta) * elem;
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
                double elem = sys->get_v_phhh()->at(lmd)(nu_ab, mu_gd);
                H_hf(alpha, gamma) += C_inner(beta, delta) * elem;
                H_hf(gamma, alpha) += C_inner(beta, delta) * elem;
                H_hf(beta, gamma) -= C_inner(alpha, delta) * elem;
                H_hf(gamma, beta) -= C_inner(delta, alpha) * elem;
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
                double elem = sys->get_v_pphh()->at(lmd)(xi_ab, mu_gd);
                H_hf(alpha, gamma) += C_inner(beta, delta) * elem;
                H_hf(gamma, alpha) += C_inner(delta, beta) * elem;
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
                double elem = sys->get_v_phph()->at(lmd)(nu_ab, nu_gd);
                H_hf(alpha, gamma) += C_inner(beta, delta) * elem;
                H_hf(alpha, delta) -= C_inner(beta, gamma) * elem;
                H_hf(beta, gamma) -= C_inner(alpha, delta) * elem;
                H_hf(beta, delta) += C_inner(alpha, gamma) * elem;
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
                double elem = sys->get_v_ppph()->at(lmd)(xi_ab, nu_gd);
                H_hf(alpha, gamma) += C_inner(beta, delta) * elem;
                H_hf(alpha, delta) -= C_inner(beta, gamma) * elem;
                H_hf(gamma, alpha) += C_inner(delta, beta) * elem;
                H_hf(delta, alpha) -= C_inner(gamma, beta) * elem;
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
                double elem = sys->get_v_pppp()->at(lmd)(xi_ab, xi_gd);
                H_hf(alpha, gamma) += C_inner(beta, delta) * elem;
            }
    }

    return H_hf;
}

double HF::energy(mat const &C_inner, mat const &h0)
{
    Basis const * basis = sys->get_basis();
    vector<uvec> const * ppMap = basis->get_map_lmdXI_de();
    vector<uvec> const * phMap = basis->get_map_lmdNU_dl();
    vector<uvec> const * hhMap = basis->get_map_lmdMU_lm();
    int dimLMD = basis->dim_lmd_2p();
    int nH = sys->get_basis()->get_nH();
    int nP = sys->get_basis()->get_nP();

    double E_ref = accu(C_inner % h0);
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
                double elem = sys->get_v_hhhh()->at(lmd)(mu_ab, mu_gd);
                E_ref += 0.5 * C_inner(alpha, gamma) * C_inner(beta, delta) * elem;
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
                double elem = sys->get_v_phhh()->at(lmd)(nu_ab, mu_gd);
                E_ref += 0.5 * C_inner(alpha, gamma) * C_inner(beta, delta) * elem;
                E_ref += 0.5 * C_inner(gamma, alpha) * C_inner(beta, delta) * elem;
                E_ref -= 0.5 * C_inner(beta, gamma) * C_inner(alpha, delta) * elem;
                E_ref -= 0.5 * C_inner(gamma, beta) * C_inner(delta, alpha) * elem;
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
                double elem = sys->get_v_pphh()->at(lmd)(xi_ab, mu_gd);
                E_ref += 0.5 * C_inner(alpha, gamma) * C_inner(beta, delta) * elem;
                E_ref += 0.5 * C_inner(gamma, alpha) * C_inner(delta, beta) * elem;
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
                double elem = sys->get_v_phph()->at(lmd)(nu_ab, nu_gd);
                E_ref += 0.5 * C_inner(alpha, gamma) * C_inner(beta, delta) * elem;
                E_ref -= 0.5 * C_inner(alpha, delta) * C_inner(beta, gamma) * elem;
                E_ref -= 0.5 * C_inner(beta, gamma) * C_inner(alpha, delta) * elem;
                E_ref += 0.5 * C_inner(beta, delta) * C_inner(alpha, gamma) * elem;
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
                double elem = sys->get_v_ppph()->at(lmd)(xi_ab, nu_gd);
                E_ref += 0.5 * C_inner(alpha, gamma) * C_inner(beta, delta) * elem;
                E_ref -= 0.5 * C_inner(alpha, delta) * C_inner(beta, gamma) * elem;
                E_ref += 0.5 * C_inner(gamma, alpha) * C_inner(delta, beta) * elem;
                E_ref -= 0.5 * C_inner(delta, alpha) * C_inner(gamma, beta) * elem;
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
                double elem = sys->get_v_pppp()->at(lmd)(xi_ab, xi_gd);
                E_ref += 0.5 * C_inner(alpha, gamma) * C_inner(beta, delta) * elem;
            }
    }

    return E_ref;
}
