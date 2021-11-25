/* 
 * File:   CC.cpp
 * Author: chrishir
 * 
 * Created on 6. desember 2011, 14:46
 */

#include "CC.h"
#include "GEMM.h"
#include <iomanip>

using namespace toffyrn::libUNK;
using namespace std;
using namespace arma;

CC::CC(double precision, int max_iter)
{
    //Matrix multiplicator (standard)
    mult = new GEMM;

    this->max_iter = max_iter;
    this->precision = precision;
    reset_timing_info();
}

CC::~CC()
{
    if (mult != NULL)
        delete(mult);
}

void CC::set_system(System* sys)
{
    this->sys = sys;
}

double CC::solve_ground_state_energy(int debug_flags, ostream &debug_stream)
{
    Basis const * basis = sys->get_basis();
    vector<uvec> const * ppMap = basis->get_map_lmdXI_de();
    vector<uvec> const * hhMap = basis->get_map_lmdMU_lm();
    umat const * ppMapInv = basis->get_map_de_lmdXI();
    umat const * hhMapInv = basis->get_map_lm_lmdMU();
    int dimLMD2 = basis->dim_lmd_2p();
    int nP = basis->get_nP();
    int nH = basis->get_nH();

    //Prepare additional mappings needed for efficient solutions.
    init_additional_mappings();

    //Print all interaction elements for DEBUG
    if (debug_flags & DEBUG_INTERACTIONS)
    {
        cout << "f_hh:\n" << *sys->get_f_hh() << endl;
        cout << "f_ph:\n" << *sys->get_f_ph() << endl;
        cout << "f_pp:\n" << *sys->get_f_pp() << endl;
        for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
            debug_stream << "v_hhhh - lmd=" << lmd << endl << sys->get_v_hhhh()->at(lmd) << endl;
        for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
            debug_stream << "v_phhh - lmd=" << lmd << endl << sys->get_v_phhh()->at(lmd) << endl;
        for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
            debug_stream << "v_phph - lmd=" << lmd << endl << sys->get_v_phph()->at(lmd) << endl;
        for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
            debug_stream << "v_pphh - lmd=" << lmd << endl << sys->get_v_pphh()->at(lmd) << endl;
        for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
            debug_stream << "v_ppph - lmd=" << lmd << endl << sys->get_v_ppph()->at(lmd) << endl;
        for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
            debug_stream << "v_pppp - lmd=" << lmd << endl << sys->get_v_pppp()->at(lmd) << endl;
    }

    //Find reference energy.
    double E_ref = 0;
    for (int k = 0; k < nH; k++)
        E_ref += (*sys->get_f_hh())(k, k);
    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        int dimMU = hhMap->at(lmd).size();
        for (int mu = 0; mu < dimMU; mu++)
            E_ref -= 0.5 * sys->get_v_hhhh()->at(lmd)(mu, mu);
    }
    if (debug_flags & DEBUG_ENERGY)
        debug_stream << std::setprecision(10) << "E_ref: " << E_ref << endl;

    //Prepare amplitudes.
    mat t1_old = zeros<mat > (sys->get_basis()->get_nP(), sys->get_basis()->get_nH());
    vector<mat> t2_old;
    {
        Basis const * basis = sys->get_basis();
        vector<uvec> const * mapPP = basis->get_map_lmdXI_de();
        vector<uvec> const * mapHH = basis->get_map_lmdMU_lm();

        for (size_t lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
        {
            size_t dimXI = mapPP->at(lmd).size();
            size_t dimMU = mapHH->at(lmd).size();

            t2_old.push_back(zeros<mat > (dimXI, dimMU));
        }
    }

    double E_ccsd_old = 0;
    for (int iter = 0; iter < max_iter; iter++)
    {
        if (debug_flags != 0)
            debug_stream << "============ ITERATION " << iter << " ==============\n";

        wall_clock time;
        mat i1 = c_i1(t1_old);
        if (debug_flags & DEBUG_I1)
            debug_stream << "i1\n" << i1 << endl;
        //
        time.tic();
        mat i2 = c_i2(t1_old);
        timing_info(2) += time.toc();
        if (debug_flags & DEBUG_I2)
            debug_stream << "i2\n" << i2 << endl;
        //
        mat i3 = c_i3(i2, t1_old, t2_old);
        if (debug_flags & DEBUG_I3)
            debug_stream << "i3\n" << i3 << endl;
        //
        vector<mat> i5 = c_i5(t1_old);
        if (debug_flags & DEBUG_I5)
            for (int idx = 0; idx < i5.size(); idx++)
                debug_stream << "i5 - lmd=" << idx << endl << i5.at(idx) << endl;
        //
        vector<mat> i4 = c_i4(i5, t1_old);
        if (debug_flags & DEBUG_I4)
            for (int idx = 0; idx < i4.size(); idx++)
                debug_stream << "i4 - lmd=" << idx << endl << i4.at(idx) << endl;
        //
        time.tic();
        mat d1 = c_d1(i1, i3);
        timing_info(14) += time.toc();
        if (debug_flags & DEBUG_D1)
            debug_stream << "d1\n" << d1 << endl;
        //
        mat t1_new = c_t1(i1, i2, i3, i4, d1, t1_old, t2_old);
        if (debug_flags & DEBUG_T1)
            debug_stream << "t1\n" << t1_new << endl;
        //
        i1 = c_i1(t1_new);
        if (debug_flags & DEBUG_I1)
            debug_stream << "i1\n" << i1 << endl;
        //
        time.tic();
        i2 = c_i2(t1_new);
        timing_info(2) += time.toc();
        if (debug_flags & DEBUG_I2)
            debug_stream << "i2\n" << i2 << endl;
        //
        i3 = c_i3(i2, t1_new, t2_old);
        if (debug_flags & DEBUG_I3)
            debug_stream << "i3\n" << i3 << endl;
        //
        i5 = c_i5(t1_new);
        if (debug_flags & DEBUG_I5)
            for (int idx = 0; idx < i5.size(); idx++)
                debug_stream << "i5 - lmd=" << idx << endl << i5.at(idx) << endl;
        //
        i4 = c_i4(i5, t1_new);
        if (debug_flags & DEBUG_I4)
            for (int idx = 0; idx < i4.size(); idx++)
                debug_stream << "i4 - lmd=" << idx << endl << i4.at(idx) << endl;
        //
        vector<mat> i6 = c_i6(i5, t1_new, t2_old);
        if (debug_flags & DEBUG_I6)
            for (int idx = 0; idx < i6.size(); idx++)
                debug_stream << "i6 - lmd=" << idx << endl << i6.at(idx) << endl;
        //
        mat i7 = c_i7(i1, i2, t1_new, t2_old);
        if (debug_flags & DEBUG_I7)
            debug_stream << "i7\n" << i7 << endl;
        //
        vector<mat> i9 = c_i9(t1_new);
        if (debug_flags & DEBUG_I9)
            for (int idx = 0; idx < i9.size(); idx++)
                debug_stream << "i9 - lmd=" << idx << endl << i9.at(idx) << endl;
        //
        vector<mat> i8 = c_i8(i4, i9, t1_new, t2_old);
        if (debug_flags & DEBUG_I8)
            for (int idx = 0; idx < i9.size(); idx++)
                debug_stream << "i8 - lmd=" << idx << endl << i8.at(idx) << endl;
        //
        time.tic();
        vector<mat> i10 = c_i10(i6, i9, t1_new, t2_old);
        timing_info(10) += time.toc();
        if (debug_flags & DEBUG_I10)
            for (int idx = 0; idx < i10.size(); idx++)
                debug_stream << "i10 - lmd=" << idx << endl << i10.at(idx) << endl;
        //
        vector<mat> i11 = c_i11(t1_new);
        if (debug_flags & DEBUG_I11)
            for (int idx = 0; idx < i11.size(); idx++)
                debug_stream << "i11 - lmd=" << idx << endl << i11.at(idx) << endl;
        //
        time.tic();
        vector<mat> d2 = c_d2(i3, i6, i7);
        timing_info(15) += time.toc();
        if (debug_flags & DEBUG_D2)
            for (int idx = 0; idx < d2.size(); idx++)
                debug_stream << "d2 - lmd=" << idx << endl << d2.at(idx) << endl;
        //
        vector<mat> t2_new = c_t2(i3, i6, i7, i8, i10, i11, d2, t1_new, t2_old);
        if (debug_flags & DEBUG_T2)
            for (int idx = 0; idx < t2_new.size(); idx++)
                debug_stream << "t2 - lmd=" << idx << endl << t2_new.at(idx) << endl;

        //
        t1_old = t1_new;
        t2_old = t2_new;


        //Calculate energy!
        time.tic();
        double E_ccsd = 0;
        E_ccsd += accu((*sys->get_f_ph()) % t1_new);
        for (int lmd2 = 0; lmd2 < dimLMD2; lmd2++)
        {
            E_ccsd += 0.25 * accu(sys->get_v_pphh()->at(lmd2) % t2_new.at(lmd2));

            int dimXI = ppMap->at(lmd2).size();
            int dimMU = hhMap->at(lmd2).size();

            for (int xi_ab = 0; xi_ab < dimXI; xi_ab++)
                for (int mu_ij = 0; mu_ij < dimMU; mu_ij++)
                {
                    int ab = ppMap->at(lmd2)(xi_ab);
                    int a = ab % nP;
                    int b = ab / nP;
                    int ij = hhMap->at(lmd2)(mu_ij);
                    int i = ij % nH;
                    int j = ij / nH;

                    E_ccsd += 0.5 * sys->get_v_pphh()->at(lmd2)(xi_ab, mu_ij) * t1_new(a, i) * t1_new(b, j);
                }
        }
        timing_info(0) += time.toc();
        if (debug_flags & DEBUG_ENERGY)
            debug_stream << std::setprecision(10) << E_ccsd << endl;

        if (fabs(E_ccsd - E_ccsd_old) < precision)
        {
            E_ccsd_old = E_ccsd;
            break;
        }
        E_ccsd_old = E_ccsd;

    }

    return E_ccsd_old + E_ref;
}

mat CC::c_i1(mat const &t1_old)
{
    //#1 f_ad
    wall_clock timer;
    timer.tic();
    mat i1 = *sys->get_f_pp();
    timing_info(1)(0) += timer.toc();

    //#2 <de||al> t^e_l 
    timer.tic();
    c_i1_t2(i1, t1_old);
    timing_info(1)(1) += timer.toc();

    return i1;
}

void CC::c_i1_t2(arma::mat &i1, arma::mat const &t1_old)
{
    //Needed constructs
    Basis const * basis = sys->get_basis();
    vector<uvec> const * mapPP = basis->get_map_lmdXI_de();
    vector<uvec> const * mapPH = basis->get_map_lmdNU_dl();
    int nP = basis->get_nP();
    int dimLMD2 = basis->dim_lmd_2p();

    //Temporary vector(!) t1_el1
    vector<vec> t1_el1;
    for (int lmd = 0; lmd < dimLMD2; lmd++)
    {
        int dimPH1 = map_p_h1.at(lmd).size();
        t1_el1.push_back(zeros<vec > (dimPH1));

        //Fill t1
        for (int idx_el1 = 0; idx_el1 < dimPH1; idx_el1++)
        {
            int el = map_p_h1.at(lmd)(idx_el1);
            int e = el % nP;
            int l = el / nP;
            t1_el1.at(lmd)(idx_el1) = t1_old(e, l);
        }
    }

    //Matmult
    vector<vec> i1_ad1;
    for (int lmd = 0; lmd < dimLMD2; lmd++)
        i1_ad1.push_back(mult->dgemm(v_ad1_el1.at(lmd), t1_el1.at(lmd)));

    //Fill i1
    for (int lmd = 0; lmd < dimLMD2; lmd++)
    {
        int dimPP1 = map_p_p1.at(lmd).size();
        for (int idx_ad1 = 0; idx_ad1 < dimPP1; idx_ad1++)
        {
            int ad = map_p_p1.at(lmd)(idx_ad1);
            int a = ad % nP;
            int d = ad / nP;
            i1(a, d) += i1_ad1.at(lmd)(idx_ad1);
        }
    }


    return;
}

mat CC::c_i2(arma::mat const &t1_old)
{
    //Needed objects
    Basis const * basis = sys->get_basis();
    vector<uvec> const * mapPP = basis->get_map_lmdXI_de();
    vector<uvec> const * mapHH = basis->get_map_lmdMU_lm();
    int nH = basis->get_nH();
    int nP = basis->get_nP();

    //#1 f_dl
    mat i2 = *sys->get_f_ph();

    //#2 <de||lm> t^e_m
    for (size_t lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        size_t dimXI = mapPP->at(lmd).size();
        size_t dimMU = mapHH->at(lmd).size();

        for (size_t xi = 0; xi < dimXI; xi++)
            for (size_t mu = 0; mu < dimMU; mu++)
            {
                int de = mapPP->at(lmd)(xi);
                int d = de % nP;
                int e = de / nP;
                int lm = mapHH->at(lmd)(mu);
                int l = lm % nH;
                int m = lm / nH;

                i2(d, l) += sys->get_v_pphh()->at(lmd)(xi, mu) * t1_old(e, m);
            }
    }

    return i2;
}

mat CC::c_i3(arma::mat const &i2, arma::mat const &t1_old, std::vector<arma::mat> const &t2_old)
{
    //Needed objects
    Basis const * basis = sys->get_basis();
    vector<uvec> const * mapPP = basis->get_map_lmdXI_de();
    vector<uvec> const * mapPH = basis->get_map_lmdNU_dl();
    vector<uvec> const * mapHH = basis->get_map_lmdMU_lm();
    umat const * mapHHinv = basis->get_map_lm_lmdMU();
    int nH = basis->get_nH();
    int nP = basis->get_nP();

    //#1 f_li
    wall_clock time;
    time.tic();
    mat i3 = *sys->get_f_hh();
    timing_info(3)(0) += time.toc();

    //#2  <di||ml> t^d_m
    time.tic();
    for (size_t lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        size_t dimNU = mapPH->at(lmd).size();
        size_t dimMU = mapHH->at(lmd).size();

        for (size_t nu = 0; nu < dimNU; nu++)
            for (size_t mu = 0; mu < dimMU; mu++)
            {
                int di = mapPH->at(lmd)(nu);
                int d = di % nP;
                int i = di / nP;
                int ml = mapHH->at(lmd)(mu);
                int m = ml % nH;
                int l = ml / nH;

                i3(l, i) += sys->get_v_phhh()->at(lmd)(nu, mu) * t1_old(d, m);
            }
    }
    timing_info(3)(1) += time.toc();

    //#3 0.5 <de||lm> t^de_im
    time.tic();
    c_i3_t3(i3, t2_old);
    timing_info(3)(2) += time.toc();

    //#4 i2_dl * t^d_i
    time.tic();
    for (size_t l = 0; l < nH; l++)
        for (size_t i = 0; i < nH; i++)
            for (size_t d = 0; d < nP; d++)
                i3(l, i) += i2(d, l) * t1_old(d, i);
    timing_info(3)(3) += time.toc();

    return i3;
}

void CC::c_i3_t3(mat &i3, vector<mat> const &t2_old)
{
    //Needed objects
    Basis const * basis = sys->get_basis();
    vector<uvec> const * mapPP = basis->get_map_lmdXI_de();
    vector<uvec> const * mapPH = basis->get_map_lmdNU_dl();
    vector<uvec> const * mapHH = basis->get_map_lmdMU_lm();
    umat const * mapHHinv = basis->get_map_lm_lmdMU();
    int nH = basis->get_nH();
    int nP = basis->get_nP();

    vector<mat> const * v_pphh = sys->get_v_pphh();
    for (int l = 0; l < nH; l++)
        for (int i = 0; i < nH; i++)
        {
            double i3_li = 0;
            for (int m = 0; m < nH; m++)
            {
                int lm = l + m * nH;
                int lmd_lm = (*mapHHinv)(0, lm);
                int mu_lm = (*mapHHinv)(1, lm);
                int im = i + m * nH;
                int lmd_im = (*mapHHinv)(0, im);
                int mu_im = (*mapHHinv)(1, im);

                if (lmd_lm == lmd_im)
                    i3_li += dot(v_pphh->at(lmd_lm)(span::all, mu_lm), t2_old.at(lmd_lm)(span::all, mu_im));
            }
            i3(l, i) += 0.5 * i3_li;
        }


    //    for (size_t lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    //    {
    //        size_t dimXI = mapPP->at(lmd).size();
    //
    //        for (size_t xi = 0; xi < dimXI; xi++)
    //            for (int l = 0; l < nH; l++)
    //                for (int m = 0; m < nH; m++)
    //                    for (int i = 0; i < nH; i++)
    //                    {
    //                        int lm = l + m * nH;
    //                        int lmd_lm = (*mapHHinv)(0, lm);
    //                        int mu_lm = (*mapHHinv)(1, lm);
    //                        int im = i + m * nH;
    //                        int lmd_im = (*mapHHinv)(0, im);
    //                        int mu_im = (*mapHHinv)(1, im);
    //
    //                        if (lmd_lm == lmd_im && lmd_lm == lmd)
    //                            i3(l, i) += 0.5 * sys->get_v_pphh()->at(lmd)(xi, mu_lm) * t2_old.at(lmd)(xi, mu_im);
    //                    }
    //    }

    return;
}

std::vector<arma::mat> CC::c_i4(
        std::vector<arma::mat> const &i5,
        arma::mat const &t1_old)
{
    //Needed objects
    Basis const * basis = sys->get_basis();
    umat const * mapPPinv = basis->get_map_de_lmdXI();
    umat const * mapPHinv = basis->get_map_dl_lmdNU();
    umat const * mapHHinv = basis->get_map_lm_lmdMU();
    int nH = basis->get_nH();
    int nP = basis->get_nP();

    //#1 i5^di_lm
    wall_clock time;
    time.tic();
    std::vector<arma::mat> i4 = i5;
    timing_info(4)(0) += time.toc();

    //#2 0.5 <ed||lm> t^e_i
    time.tic();
    c_i5_t2(i4, t1_old); //t2 of i4 and i5 are equal.
    timing_info(4)(1) += time.toc();

    return i4;
}

std::vector<arma::mat> CC::c_i5(arma::mat const &t1_old)
{
    //Needed objects
    Basis const * basis = sys->get_basis();
    umat const * mapPPinv = basis->get_map_de_lmdXI();
    umat const * mapPHinv = basis->get_map_dl_lmdNU();
    umat const * mapHHinv = basis->get_map_lm_lmdMU();
    int nH = basis->get_nH();
    int nP = basis->get_nP();

    //#1 - <di||lm>
    wall_clock time;
    time.tic();
    std::vector<arma::mat> i5 = *sys->get_v_phhh();
    for (size_t idx = 0; idx < i5.size(); idx++)
        i5.at(idx) = -i5.at(idx);
    timing_info(5)(0) += time.toc();

    //#2 0.5 <ed||lm> t^e_i
    time.tic();
    c_i5_t2(i5, t1_old);
    timing_info(5)(1) += time.toc();

    return i5;
}

void CC::c_i5_t2(vector<mat> &i5, mat const &t1_old)
{
    //Needed objects
    Basis const * basis = sys->get_basis();
    vector<uvec> const * mapPP = basis->get_map_lmdXI_de();
    vector<uvec> const * mapPH = basis->get_map_lmdNU_dl();
    vector<uvec> const * mapHH = basis->get_map_lmdMU_lm();
    umat const * mapPPinv = basis->get_map_de_lmdXI();
    umat const * mapPHinv = basis->get_map_dl_lmdNU();
    umat const * mapHHinv = basis->get_map_lm_lmdMU();
    int nH = basis->get_nH();
    int nP = basis->get_nP();

    for (int lmd2 = 0; lmd2 < basis->dim_lmd_2p(); lmd2++)
    {
        int dimXI = mapPP->at(lmd2).size();
        int dimMU = mapHH->at(lmd2).size();

        for (int xi_ed = 0; xi_ed < dimXI; xi_ed++)
        {
            int ed = mapPP->at(lmd2)(xi_ed);
            int e = ed % nP;
            int d = ed / nP;

            for (int mu_lm = 0; mu_lm < dimMU; mu_lm++)
            {
                int lmd1 = map_p_inv(0, e);
                int dimH = map_h.at(lmd1).size();

                for (int idx_i = 0; idx_i < dimH; idx_i++)
                {
                    int i = map_h.at(lmd1)(idx_i);
                    int di = d + i * nP;
                    int lmd_di = (*mapPHinv)(0, di);
                    int nu_di = (*mapPHinv)(1, di);

                    if (lmd_di != lmd2)
                        throw string("Wrong lmd_di in c_i5_t2.");

                    i5.at(lmd2)(nu_di, mu_lm) += 0.5 *
                            sys->get_v_pphh()->at(lmd2)(xi_ed, mu_lm) * t1_old(e, i);
                }
            }
        }
    }

    return;
}

mat CC::c_d1(arma::mat const &i1, arma::mat const &i3)
{
    int nP = sys->get_basis()->get_nP();
    int nH = sys->get_basis()->get_nH();

    mat d1(nP, nH);
    for (int a = 0; a < nP; a++)
        for (int i = 0; i < nH; i++)
            d1(a, i) = -i1(a, a) + i3(i, i);

    return d1;
}

vector<mat> CC::c_i6(
        std::vector<arma::mat> const &i5,
        arma::mat const &t1_old,
        std::vector<arma::mat> const &t2_old)
{
    Basis const * basis = sys->get_basis();

    //#1 <lm||ij>
    wall_clock timer;
    timer.tic();
    vector<mat> i6 = *sys->get_v_hhhh();
    timing_info(6)(0) += timer.toc();

    //#2 0.5 <de||lm> t^de_ij
    // (lmd) i6^mu1_mu2 += 0.5 <xi||mu1> t^xi_mu2
    timer.tic();
    for (size_t lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
        mult->dgemm(i6.at(lmd), sys->get_v_pphh()->at(lmd), t2_old.at(lmd), 0.5, 1.0, true, false);
    timing_info(6)(1) += timer.toc();

    timer.tic();
    c_i6_t3(i6, i5, t1_old);
    timing_info(6)(2) += timer.toc();


    return i6;
}

void CC::c_i6_t3(vector<mat> &i6, vector<mat> const &i5, mat const &t1_old)
{
    //Needed objects
    Basis const * basis = sys->get_basis();
    vector<uvec> const * mapPH = basis->get_map_lmdNU_dl();
    umat const * mapPHinv = basis->get_map_dl_lmdNU();
    umat const * mapHHinv = basis->get_map_lm_lmdMU();
    int nH = basis->get_nH();
    int nP = basis->get_nP();
    int dimLMD2 = basis->dim_lmd_2p();

    //#3 i5^di_lm t^d_j
    for (int lmd = 0; lmd < dimLMD2; lmd++)
    {
        int dimNU = mapPH->at(lmd).size();

        if (i5.at(lmd).n_elem != 0)
        {
            for (int nu_di = 0; nu_di < dimNU; nu_di++)
            {
                int di = mapPH->at(lmd)(nu_di);
                int d = di % nP;
                int i = di / nP;

                vec i5_ALL_di = trans(i5.at(lmd).row(nu_di));
                for (int j = 0; j < nH; j++)
                {
                    int ij = i + j * nH;
                    int lmd_ij = (*mapHHinv)(0, ij);
                    int mu_ij = (*mapHHinv)(1, ij);
                    if (lmd_ij == lmd)
                        i6.at(lmd).col(mu_ij) += t1_old(d, j) * i5_ALL_di;

                    int ji = j + i * nH;
                    int lmd_ji = (*mapHHinv)(0, ji);
                    int mu_ji = (*mapHHinv)(1, ji);
                    if (lmd_ji == lmd)
                        i6.at(lmd).col(mu_ji) -= t1_old(d, j) * i5_ALL_di;
                }
            }
        }
    }


    return;
}

mat CC::c_i7(
        mat const &i1,
        mat const &i2,
        mat const &t1_old,
        vector<mat> const &t2_old)
{
    //#1 i1_bd
    wall_clock timer;
    timer.tic();
    mat i7 = i1;
    timing_info(7)(0) += timer.toc();

    //#2 i2_dl t^b_l
    timer.tic();
    mult->dgemm(i7, t1_old, i2, -1, 1, false, true);
    timing_info(7)(1) += timer.toc();

    //#3 0.5 <de||lm> t^eb_lm
    timer.tic();
    c_i7_t3(i7, t2_old);
    timing_info(7)(2) += timer.toc();

    return i7;
}

void CC::c_i7_t3(mat &i7, vector<mat> const &t2_old)
{
    Basis const * basis = sys->get_basis();
    vector<uvec> const * mapHH = basis->get_map_lmdMU_lm();
    vector<uvec> const * mapPP = basis->get_map_lmdXI_de();
    umat const * mapPPinv = basis->get_map_de_lmdXI();
    umat const * mapHHinv = basis->get_map_lm_lmdMU();
    int nP = basis->get_nP();
    int nH = basis->get_nH();
    int dimLMD1 = basis->dim_lmd_1p();

    //Temporary storage t2_lme1_b
    vector<mat> t2_lme1_b;
    for (int lmd1 = 0; lmd1 < dimLMD1; lmd1++)
    {
        int dimP = map_p.at(lmd1).size();
        int dimHHP1 = map_h_h_p1.at(lmd1).size();
        t2_lme1_b.push_back(zeros<mat > (dimHHP1, dimP));

        for (int idx_lme1 = 0; idx_lme1 < dimHHP1; idx_lme1++)
            for (int idx_b = 0; idx_b < dimP; idx_b++)
            {
                int b = map_p.at(lmd1)(idx_b);
                int lme = map_h_h_p1.at(lmd1)(idx_lme1);
                int l = lme % nH;
                int me = lme / nH;
                int m = me % nH;
                int e = me / nH;

                int eb = e + b * nP;
                int lmd2_eb = (*mapPPinv)(0, eb);
                int xi_eb = (*mapPPinv)(1, eb);
                int lm = l + m * nH;
                int lmd2_lm = (*mapHHinv)(0, lm);
                int mu_lm = (*mapHHinv)(1, lm);
                if (lmd2_eb != lmd2_lm)
                    throw std::string("Mismatching lmd's when filling t2_lme1_b in c_i7_t3");

                t2_lme1_b.at(lmd1)(idx_lme1, idx_b) = t2_old.at(lmd2_lm)(xi_eb, mu_lm);
            }
    }

    //Matmult
    vector<mat> i7_d_b;
    for (int lmd1 = 0; lmd1 < dimLMD1; lmd1++)
        i7_d_b.push_back(mult->dgemm(v_d_lme1.at(lmd1), t2_lme1_b.at(lmd1), 0.5));

    //Fill back into i7
    for (int lmd1 = 0; lmd1 < dimLMD1; lmd1++)
    {
        int dimP = map_p.at(lmd1).size();
        for (int idx_d = 0; idx_d < dimP; idx_d++)
            for (int idx_b = 0; idx_b < dimP; idx_b++)
            {
                int d = map_p.at(lmd1)(idx_d);
                int b = map_p.at(lmd1)(idx_b);
                i7(b, d) += i7_d_b.at(lmd1)(idx_d, idx_b);
            }
    }

    return;
}

vector<mat> CC::c_i9(arma::mat const &t1_old)
{
    Basis const * basis = sys->get_basis();
    umat const * mapPPinv = basis->get_map_de_lmdXI();
    umat const * mapPHinv = basis->get_map_dl_lmdNU();
    vector<uvec> const * mapPH = basis->get_map_lmdNU_dl();
    int nP = basis->get_nP();
    int nH = basis->get_nH();

    //#1 -<bl||dj>
    wall_clock timer;
    timer.tic();
    vector<mat> i9 = *sys->get_v_phph();
    for (int idx = 0; idx < i9.size(); idx++)
        i9.at(idx) = -i9.at(idx);
    timing_info(9)(0) += timer.toc();

    //#2 0.5 <ed||bl>t^e_j
    timer.tic();
    c_i9_t2(i9, t1_old);
    timing_info(9)(1) += timer.toc();

    return i9;
}

void CC::c_i9_t2(vector<mat> &i9, mat const &t1_old)
{
    Basis const * basis = sys->get_basis();
    umat const * mapPPinv = basis->get_map_de_lmdXI();
    umat const * mapPHinv = basis->get_map_dl_lmdNU();
    vector<uvec> const * mapPH = basis->get_map_lmdNU_dl();
    int nP = basis->get_nP();
    int nH = basis->get_nH();

    vector<mat> v_bl_ed;
    for (int lmd2 = 0; lmd2 < basis->dim_lmd_2p(); lmd2++)
        v_bl_ed.push_back(trans(sys->get_v_ppph()->at(lmd2)));

    for (int e = 0; e < nP; e++)
        for (int d = 0; d < nP; d++)
        {
            int ed = e + d * nP;
            int lmd_ed = (*mapPPinv)(0, ed);
            int xi_ed = (*mapPPinv)(1, ed);
            mat v_all_ed = v_bl_ed.at(lmd_ed)(span::all, xi_ed);

            for (int j = 0; j < nH; j++)
            {
                int dj = d + j * nP;
                int lmd_dj = (*mapPHinv)(0, dj);
                int nu_dj = (*mapPHinv)(1, dj);
                if (lmd_ed == lmd_dj)
                    i9.at(lmd_dj)(span::all, nu_dj) += 0.5 * t1_old(e, j) * v_all_ed;
            }
        }

    return;
}

std::vector<arma::mat> CC::c_i8(
        std::vector<arma::mat> const &i4,
        std::vector<arma::mat> const &i9,
        arma::mat const &t1_old,
        std::vector<arma::mat> const &t2_old)
{
    //#1 i9^bl_dj
    wall_clock time;
    time.tic();
    vector<mat> i8 = i9;
    timing_info(8)(0) += time.toc();

    //#2 0.5 <ed||bl> t^e_j
    time.tic();
    c_i8_t2(i8, t1_old);
    timing_info(8)(1) += time.toc();

    //#3 i4^dj_lm t^b_m
    time.tic();
    c_i8_t3(i8, i4, t1_old);
    timing_info(8)(2) += time.toc();

    //#4 0.5 <de||lm> t^eb_mj
    time.tic();
    c_i8_t4(i8, t2_old);
    timing_info(8)(3) += time.toc();

    return i8;
}

void CC::c_i8_t2(vector<mat> &i8, mat const &t1_old)
{
    Basis const * basis = sys->get_basis();
    vector<uvec> const * mapHH = basis->get_map_lmdMU_lm();
    vector<uvec> const * mapPH = basis->get_map_lmdNU_dl();
    vector<uvec> const * mapPP = basis->get_map_lmdXI_de();
    umat const * mapPPinv = basis->get_map_de_lmdXI();
    umat const * mapPHinv = basis->get_map_dl_lmdNU();
    umat const * mapHHinv = basis->get_map_lm_lmdMU();
    int nP = basis->get_nP();
    int nH = basis->get_nH();
    int dimLMD2 = basis->dim_lmd_2p();

    for (int lmd2 = 0; lmd2 < dimLMD2; lmd2++)
    {
        int dimXI = mapPP->at(lmd2).size();
        int dimNU = mapPH->at(lmd2).size();
        mat const &v_ppph = trans(sys->get_v_ppph()->at(lmd2));

        for (int xi_ed = 0; xi_ed < dimXI; xi_ed++)
            for (int nu_dj = 0; nu_dj < dimNU; nu_dj++)
            {
                int ed = mapPP->at(lmd2)(xi_ed);
                int e = ed % nP;
                int d1 = ed / nP;
                int dj = mapPH->at(lmd2)(nu_dj);
                int d2 = dj % nP;
                int j = dj / nP;

                if (d1 == d2)
                    i8.at(lmd2).col(nu_dj) += 0.5 * v_ppph.col(xi_ed) * t1_old(e, j);
            }
    }

    return;
}

void CC::c_i8_t3(vector<mat> &i8, vector<mat> const &i4, mat const &t1_old)
{
    Basis const * basis = sys->get_basis();
    vector<uvec> const * mapHH = basis->get_map_lmdMU_lm();
    vector<uvec> const * mapPH = basis->get_map_lmdNU_dl();
    vector<uvec> const * mapPP = basis->get_map_lmdXI_de();
    umat const * mapPPinv = basis->get_map_de_lmdXI();
    umat const * mapPHinv = basis->get_map_dl_lmdNU();
    umat const * mapHHinv = basis->get_map_lm_lmdMU();
    int nP = basis->get_nP();
    int nH = basis->get_nH();


    for (int b = 0; b < nP; b++)
        for (int l = 0; l < nH; l++)
            for (int m = 0; m < nH; m++)
            {
                int bl = b + l * nP;
                int lmd_bl = (*mapPHinv)(0, bl);
                int nu_bl = (*mapPHinv)(1, bl);
                int lm = l + m * nH;
                int lmd_lm = (*mapHHinv)(0, lm);
                int mu_lm = (*mapHHinv)(1, lm);

                if (lmd_bl == lmd_lm)
                    i8.at(lmd_bl)(nu_bl, span::all) += t1_old(b, m) * trans(i4.at(lmd_lm)(span::all, mu_lm));
            }


    return;
}

void CC::c_i8_t4(vector<mat> &i8, vector<mat> const &t2_old)
{
    Basis const * basis = sys->get_basis();
    vector<uvec> const * mapHH = basis->get_map_lmdMU_lm();
    vector<uvec> const * mapPH = basis->get_map_lmdNU_dl();
    vector<uvec> const * mapPP = basis->get_map_lmdXI_de();
    umat const * mapPPinv = basis->get_map_de_lmdXI();
    umat const * mapPHinv = basis->get_map_dl_lmdNU();
    umat const * mapHHinv = basis->get_map_lm_lmdMU();
    int nP = basis->get_nP();
    int nH = basis->get_nH();

    //Temporary storage 
    vector<mat> i8_dl_bj;
    vector<mat> t_em_bj;
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_2p(); lmd1++)
    {
        int dimP_H1 = map_p_h1.at(lmd1).size();
        int dimP1_H = map_p1_h.at(lmd1).size();
        t_em_bj.push_back(zeros<mat > (dimP_H1, dimP1_H));
    }
    //Fill storage
    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        int dimMU = mapHH->at(lmd).size();
        int dimXI = mapPP->at(lmd).size();
        for (int xi_eb = 0; xi_eb < dimXI; xi_eb++)
            for (int mu_mj = 0; mu_mj < dimMU; mu_mj++)
            {
                int eb = mapPP->at(lmd)(xi_eb);
                int e = eb % nP;
                int b = eb / nP;
                int mj = mapHH->at(lmd)(mu_mj);
                int m = mj % nH;
                int j = mj / nH;

                int bj = b + j * nP;
                int em = e + m * nP;
                int lmd1_b1_j = map_p1_h_inv(0, bj);
                int b1_j = map_p1_h_inv(1, bj);
                int lmd1_e_m1 = map_p_h1_inv(0, em);
                int e_m1 = map_p_h1_inv(1, em);

                if (lmd1_e_m1 != lmd1_b1_j)
                    cerr << "WRONG lmd1's!!!\n";
                t_em_bj.at(lmd1_b1_j)(e_m1, b1_j) = t2_old.at(lmd)(xi_eb, mu_mj);
            }
    }

    //Matrix product.
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_2p(); lmd1++)
        i8_dl_bj.push_back(mult->dgemm(v_d1l_em1.at(lmd1), t_em_bj.at(lmd1), 0.5));

    //Map i8 back
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_2p(); lmd1++)
    {
        int dimP1_H = map_p1_h.at(lmd1).size();
        for (int d1_l = 0; d1_l < dimP1_H; d1_l++)
            for (int b1_j = 0; b1_j < dimP1_H; b1_j++)
            {
                int dl = map_p1_h.at(lmd1)(d1_l);
                int d = dl % nP;
                int l = dl / nP;
                int bj = map_p1_h.at(lmd1)(b1_j);
                int b = bj % nP;
                int j = bj / nP;

                int bl = b + l *nP;
                int lmd_bl = (*mapPHinv)(0, bl);
                int nu_bl = (*mapPHinv)(1, bl);
                int dj = d + j* nP;
                int lmd_dj = (*mapPHinv)(0, dj);
                int nu_dj = (*mapPHinv)(1, dj);

                if (lmd_dj != lmd_bl)
                    cerr << "WRONG lmd's!!!!\n";
                i8.at(lmd_dj)(nu_bl, nu_dj) += i8_dl_bj.at(lmd1)(d1_l, b1_j);
            }
    }

    return;
}

std::vector<arma::mat> CC::c_i10(
        std::vector<arma::mat> const &i6,
        std::vector<arma::mat> const &i9,
        arma::mat const &t1_old,
        std::vector<arma::mat> const &t2_old)
{
    Basis const * basis = sys->get_basis();
    vector<uvec> const * mapPH = basis->get_map_lmdNU_dl();
    umat const * mapPHinv = basis->get_map_dl_lmdNU();
    vector<uvec> const * mapHH = basis->get_map_lmdMU_lm();
    umat const * mapHHinv = basis->get_map_lm_lmdMU();
    int nP = basis->get_nP();
    int nH = basis->get_nH();

    //#1 <al||ij> 
    vector<mat> i10 = *sys->get_v_phhh();

    //#2 0.5 <de||al> t^de_ij
    vector<mat> const &v_ppph = *sys->get_v_ppph();
    for (int lmd = 0; lmd < i10.size(); lmd++)
        mult->dgemm(i10.at(lmd), v_ppph.at(lmd), t2_old.at(lmd), 0.5, 1.0, true, false);

    //#3 i9^al_di t^d_j
    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        int dimNU = mapPH->at(lmd).size();
        int dimMU = mapHH->at(lmd).size();

        for (int nu_al = 0; nu_al < dimNU; nu_al++)
            for (int mu_ij = 0; mu_ij < dimMU; mu_ij++)
            {
                int ij = mapHH->at(lmd)(mu_ij);
                int i = ij % nH;
                int j = ij / nH;

                for (int d = 0; d < nP; d++)
                {
                    int di = d + i * nP;
                    int lmd_di = (*mapPHinv)(0, di);
                    int nu_di = (*mapPHinv)(1, di);

                    if (lmd == lmd_di)
                        i10.at(lmd)(nu_al, mu_ij) += i9.at(lmd)(nu_al, nu_di) * t1_old(d, j);
                }
            }
    }

    //#4 -i9^al_dj t^d_i
    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        int dimNU = mapPH->at(lmd).size();
        int dimMU = mapHH->at(lmd).size();

        for (int nu_al = 0; nu_al < dimNU; nu_al++)
            for (int mu_ij = 0; mu_ij < dimMU; mu_ij++)
            {
                int ij = mapHH->at(lmd)(mu_ij);
                int i = ij % nH;
                int j = ij / nH;

                for (int d = 0; d < nP; d++)
                {
                    int dj = d + j * nP;
                    int lmd_dj = (*mapPHinv)(0, dj);
                    int nu_dj = (*mapPHinv)(1, dj);

                    if (lmd == lmd_dj)
                        i10.at(lmd)(nu_al, mu_ij) -= i9.at(lmd)(nu_al, nu_dj) * t1_old(d, i);
                }
            }
    }

    //#5 0.5 i6^lm_ij t^a_m
    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        int dimNU = mapPH->at(lmd).size();
        int dimMU = mapHH->at(lmd).size();

        for (int nu_al = 0; nu_al < dimNU; nu_al++)
            for (int mu_ij = 0; mu_ij < dimMU; mu_ij++)
            {
                int al = mapPH->at(lmd)(nu_al);
                int a = al % nP;
                int l = al / nP;

                for (int m = 0; m < nH; m++)
                {
                    int lm = l + m * nH;
                    int lmd_lm = (*mapHHinv)(0, lm);
                    int mu_lm = (*mapHHinv)(1, lm);

                    if (lmd == lmd_lm)
                        i10.at(lmd)(nu_al, mu_ij) += 0.5 * i6.at(lmd)(mu_lm, mu_ij) * t1_old(a, m);
                }
            }
    }
    //i10^al_ij

    return i10;
}

std::vector<arma::mat> CC::c_i11(arma::mat const &t1_old)
{
    //#1 <ab||dj>
    wall_clock time;
    time.tic();
    vector<mat> i11 = *sys->get_v_ppph();
    timing_info(11)(0) += time.toc();

    //#2  0.5 <ab||de> t^e_j
    time.tic();
    c_i11_t2(i11, t1_old);
    timing_info(11)(1) += time.toc();

    return i11;
}

void CC::c_i11_t2(vector<mat> &i11, mat const &t1_old)
{
    Basis const * basis = sys->get_basis();
    vector<uvec> const * mapPP = basis->get_map_lmdXI_de();
    umat const * mapPPinv = basis->get_map_de_lmdXI();
    vector<uvec> const * mapPH = basis->get_map_lmdNU_dl();
    umat const * mapPHinv = basis->get_map_dl_lmdNU();
    int nP = basis->get_nP();
    int nH = basis->get_nH();
    int dimLMD1 = basis->dim_lmd_1p();

    //Rewriting t1
    vector<mat> t1_e_j;
    for (int lmd1 = 0; lmd1 < dimLMD1; lmd1++)
    {
        int dimP = map_p.at(lmd1).size();
        int dimH = map_h.at(lmd1).size();
        t1_e_j.push_back(zeros<mat > (dimP, dimH));

        for (int e_idx = 0; e_idx < dimP; e_idx++)
            for (int j_idx = 0; j_idx < dimH; j_idx++)
            {
                int e = map_p.at(lmd1)(e_idx);
                int j = map_h.at(lmd1)(j_idx);

                t1_e_j.at(lmd1)(e_idx, j_idx) = t1_old(e, j);
            }
    }

    //Matrix mult.
    vector<mat> i11_abd1_j;
    for (int lmd1 = 0; lmd1 < dimLMD1; lmd1++)
        i11_abd1_j.push_back(mult->dgemm(v_abd1_e.at(lmd1), t1_e_j.at(lmd1), 0.5));

    //Add terms to i11
    for (int lmd1 = 0; lmd1 < dimLMD1; lmd1++)
    {
        int dimABD1 = map_p_p_p1.at(lmd1).size();
        int dimJ = map_h.at(lmd1).size();

        for (int abd_idx = 0; abd_idx < dimABD1; abd_idx++)
            for (int j_idx = 0; j_idx < dimJ; j_idx++)
            {
                int abd = map_p_p_p1.at(lmd1)(abd_idx);
                int a = abd % nP;
                int bd = abd / nP;
                int b = bd % nP;
                int d = bd / nP;
                int j = map_h.at(lmd1)(j_idx);

                int ab = a + b * nP;
                int dj = d + j * nP;

                int lmd2_ab = (*mapPPinv)(0, ab);
                int xi_ab = (*mapPPinv)(1, ab);
                int lmd2_dj = (*mapPHinv)(0, dj);
                int nu_dj = (*mapPHinv)(1, dj);
                if (lmd2_ab != lmd2_dj)
                    throw string("Mismatch between lmd2_ab and lmd2_dj in i11 term 2.");
                i11.at(lmd2_ab)(xi_ab, nu_dj) += i11_abd1_j.at(lmd1)(abd_idx, j_idx);
            }
    }

    return;
}

std::vector<arma::mat> CC::c_d2(
        arma::mat const &i3,
        std::vector<arma::mat> const &i6,
        arma::mat const &i7)
{
    Basis const * basis = sys->get_basis();
    vector<uvec> const * mapPP = basis->get_map_lmdXI_de();
    vector<uvec> const * mapHH = basis->get_map_lmdMU_lm();
    int nP = basis->get_nP();
    int nH = basis->get_nH();

    // -0.5<ab||ab> + i3_ii + i3_jj 
    // -0.5 i6^ij_ij - i7_bb - i7_aa
    std::vector<arma::mat> d2;
    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        int dimXI = mapPP->at(lmd).size();
        int dimMU = mapHH->at(lmd).size();
        mat ximu(dimXI, dimMU);

        for (int xi = 0; xi < dimXI; xi++)
            for (int mu = 0; mu < dimMU; mu++)
            {
                //-0.5 <ab||ab>
                ximu(xi, mu) = -0.5 * sys->get_v_pppp()->at(lmd)(xi, xi);

                // + i3_ii + i3_jj
                int ij = mapHH->at(lmd)(mu);
                int i = ij % nH;
                int j = ij / nH;
                ximu(xi, mu) += i3(i, i) + i3(j, j);

                // - 0.5 i6^ij_ij
                ximu(xi, mu) += -0.5 * i6.at(lmd)(mu, mu);

                // - i7_aa - i7_bb
                int ab = mapPP->at(lmd)(xi);
                int a = ab % nP;
                int b = ab / nP;
                ximu(xi, mu) += -i7(a, a) - i7(b, b);
            }

        d2.push_back(ximu);
    }

    return d2;
}

arma::mat CC::c_t1(
        arma::mat const &i1,
        arma::mat const &i2,
        arma::mat const &i3,
        std::vector<arma::mat> const &i4,
        arma::mat const &d1,
        arma::mat const &t1_old,
        std::vector<arma::mat> const &t2_old)
{
    //#1 f_ai
    wall_clock time;
    time.tic();
    arma::mat t1_new = *sys->get_f_ph();
    timing_info(12)(0) += time.toc();

    //#2 <la||di> t^d_l ->  - <al||di> t^d_l
    time.tic();
    c_t1_t2(t1_new, t1_old);
    timing_info(12)(1) += time.toc();

    //#3 0.5 <de||al> t^de_il
    time.tic();
    c_t1_t3(t1_new, t2_old);
    timing_info(12)(2) += time.toc();

    //#4 i1_ad t^d_i
    time.tic();
    mult->dgemm(t1_new, i1, t1_old, 1.0, 1.0);
    timing_info(12)(3) += time.toc();

    //#5 - i3_li t^a_l
    time.tic();
    mult->dgemm(t1_new, t1_old, i3, -1.0, 1.0);
    timing_info(12)(4) += time.toc();

    //#6 0.5 i4^di_lm t^da_lm
    time.tic();
    c_t1_t6(t1_new, i4, t2_old);
    timing_info(12)(5) += time.toc();

    //#7 i2_dl t^ad_il 
    time.tic();
    c_t1_t7(t1_new, i2, t2_old);
    timing_info(12)(6) += time.toc();

    // Add and divide by denominator.
    t1_new = t1_new / d1 + t1_old;

    return t1_new;
}

void CC::c_t1_t2(mat &t1_new, mat const &t1_old)
{
    Basis const * basis = sys->get_basis();
    vector <uvec> const * mapPH = basis->get_map_lmdNU_dl();
    int nP = basis->get_nP();
    int dimLMD2 = basis->dim_lmd_2p();

    for (int lmd2 = 0; lmd2 < dimLMD2; lmd2++)
    {
        int dimNU = mapPH->at(lmd2).size();
        for (int nu_al = 0; nu_al < dimNU; nu_al++)
            for (int nu_di = 0; nu_di < dimNU; nu_di++)
            {
                int al = mapPH->at(lmd2)(nu_al);
                int a = al % nP;
                int l = al / nP;
                int di = mapPH->at(lmd2)(nu_di);
                int d = di % nP;
                int i = di / nP;
                t1_new(a, i) -= sys->get_v_phph()->at(lmd2)(nu_al, nu_di) * t1_old(d, l);
            }
    }

    return;
}

void CC::c_t1_t3(mat& t1_new, vector<mat> const &t2_old)
{
    Basis const * basis = sys->get_basis();
    umat const * ppMapInv = basis->get_map_de_lmdXI();
    umat const * phMapInv = basis->get_map_dl_lmdNU();
    umat const * hhMapInv = basis->get_map_lm_lmdMU();
    int nP = basis->get_nP();
    int nH = basis->get_nH();
    vector<mat> const * v_ppph = sys->get_v_ppph();

    for (int a = 0; a < nP; a++)
        for (int i = 0; i < nH; i++)
        {
            double tai = 0;
            for (int l = 0; l < nH; l++)
            {
                int al = a + l * nP;
                int lmd_al = (*phMapInv)(0, al);
                int nu_al = (*phMapInv)(1, al);
                int il = i + l *nH;
                int lmd_il = (*hhMapInv)(0, il);
                int mu_il = (*hhMapInv)(1, il);

                if (lmd_al == lmd_il)
                    tai += dot(v_ppph->at(lmd_al).col(nu_al), t2_old.at(lmd_al).col(mu_il));
            }
            t1_new(a, i) += 0.5 * tai;
        }

    return;
}

void CC::c_t1_t6(
        mat &t1_new,
        vector<mat> const &i4,
        vector<mat> const &t2_old)
{
    Basis const * basis = sys->get_basis();
    umat const * ppMapInv = basis->get_map_de_lmdXI();
    umat const * phMapInv = basis->get_map_dl_lmdNU();
    umat const * hhMapInv = basis->get_map_lm_lmdMU();
    int nP = basis->get_nP();
    int nH = basis->get_nH();

    for (int i = 0; i < nH; i++)
        for (int a = 0; a < nP; a++)
        {
            double tai = 0;
            for (int d = 0; d < nP; d++)
            {
                int di = d + i * nP;
                int lmd_di = (*phMapInv)(0, di);
                int nu_di = (*phMapInv)(1, di);
                int da = d + a * nP;
                int lmd_da = (*ppMapInv)(0, da);
                int xi_da = (*ppMapInv)(1, da);

                if (lmd_da == lmd_di)
                    tai += dot(i4.at(lmd_da)(nu_di, span::all), t2_old.at(lmd_da)(xi_da, span::all));
            }
            t1_new(a, i) += 0.5 * tai;
        }

    return;
}

void CC::c_t1_t7(mat &t1_new, mat const &i2, vector< mat> const &t2_old)
{
    Basis const * basis = sys->get_basis();
    vector<uvec> const * mapPP = basis->get_map_lmdXI_de();
    vector<uvec> const * mapHH = basis->get_map_lmdMU_lm();
    int nP = basis->get_nP();
    int nH = basis->get_nH();
    int dimLMD2 = basis->dim_lmd_2p();

    for (int lmd2 = 0; lmd2 < dimLMD2; lmd2++)
    {
        int dimXI = mapPP->at(lmd2).size();
        int dimMU = mapHH->at(lmd2).size();

        for (int xi_ad = 0; xi_ad < dimXI; xi_ad++)
            for (int mu_il = 0; mu_il < dimMU; mu_il++)
            {
                int ad = mapPP->at(lmd2)(xi_ad);
                int a = ad % nP;
                int d = ad / nP;
                int il = mapHH->at(lmd2)(mu_il);
                int i = il % nH;
                int l = il / nH;

                t1_new(a, i) += i2(d, l) * t2_old.at(lmd2)(xi_ad, mu_il);
            }
    }

    return;
}

std::vector<arma::mat> CC::c_t2(
        arma::mat const &i3,
        std::vector<arma::mat> const &i6,
        arma::mat const &i7,
        std::vector<arma::mat> const &i8,
        std::vector<arma::mat> const &i10,
        std::vector<arma::mat> const &i11,
        std::vector<arma::mat> const &d2,
        arma::mat const &t1_old,
        std::vector<arma::mat> const &t2_old)
{
    Basis const * basis = sys->get_basis();
    umat const * ppMapInv = basis->get_map_de_lmdXI();
    umat const * phMapInv = basis->get_map_dl_lmdNU();
    umat const * hhMapInv = basis->get_map_lm_lmdMU();
    vector<uvec> const * ppMap = basis->get_map_lmdXI_de();
    vector<uvec> const * phMap = basis->get_map_lmdNU_dl();
    vector<uvec> const * hhMap = basis->get_map_lmdMU_lm();
    int nP = basis->get_nP();
    int nH = basis->get_nH();

    //#1 <ab||ij> 
    wall_clock timer;
    timer.tic();
    vector<mat> t2_new = *sys->get_v_pphh();
    timing_info(13)(0) += timer.toc();
   
    //#2 0.5 <ab||de> t^de_ij
    timer.tic();
    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
        mult->dgemm(t2_new.at(lmd), sys->get_v_pppp()->at(lmd), t2_old.at(lmd), 0.5, 1);
    timing_info(13)(1) += timer.toc();

    //#3 - i3_li t^ab_lj - i3_lj t^ab_il
    timer.tic();
    c_t2_t3(t2_new, i3, t2_old);
    timing_info(13)(2) += timer.toc();


    //#4 0.5 i6^lm_ij t^ab_lm
    timer.tic();
    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
        mult->dgemm(t2_new.at(lmd), t2_old.at(lmd), i6.at(lmd), 0.5, 1.0);
    timing_info(13)(3) += timer.toc();

    //#5 + i7_bd t^ad_ij + i7_ad t^db_ij
    timer.tic();
    c_t2_t5(t2_new, i7, t2_old);
    timing_info(13)(4) += timer.toc();

    //#6 + P_ij P_ab i8^bl_dj t^ad_il
    timer.tic();
    c_t2_t6(t2_new, i8, t2_old);
    timing_info(13)(5) += timer.toc();

    //#7 - P_ab i10^al_ij t^b_l
    timer.tic();
    c_t2_t7(t2_new, i10, t1_old);
    timing_info(13)(6) += timer.toc();

    //#8 + P_ij i11^ab_dj t^d_i
    timer.tic();
    c_t2_t8(t2_new, i11, t1_old);
    timing_info(13)(7) += timer.toc();

    // Add and divide by denominator.    
    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
        t2_new.at(lmd) = t2_new.at(lmd) / d2.at(lmd) + t2_old.at(lmd);

    return t2_new;
}

void CC::c_t2_t3(vector<mat> &t2_new, mat const &i3, vector<mat> const &t2_old)
{
    Basis const * basis = sys->get_basis();
    umat const * ppMapInv = basis->get_map_de_lmdXI();
    umat const * phMapInv = basis->get_map_dl_lmdNU();
    umat const * hhMapInv = basis->get_map_lm_lmdMU();
    vector<uvec> const * ppMap = basis->get_map_lmdXI_de();
    vector<uvec> const * phMap = basis->get_map_lmdNU_dl();
    vector<uvec> const * hhMap = basis->get_map_lmdMU_lm();
    int nP = basis->get_nP();
    int nH = basis->get_nH();

    //Temporary storage
    vector<mat> t2old_abj1_l;
    vector<mat> t2old_abi1_l;
    vector<mat> i3_l_i;
    vector<mat> &i3_l_j = i3_l_i;
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        int dimPPH1 = map_p_p_h1.at(lmd1).size();
        int dimH = map_h.at(lmd1).size();
        t2old_abj1_l.push_back(zeros<mat > (dimPPH1, dimH));
        t2old_abi1_l.push_back(zeros<mat > (dimPPH1, dimH));
        i3_l_i.push_back(zeros<mat > (dimH, dimH));
    }
    //Fill t2old
    for (int lmd2 = 0; lmd2 < basis->dim_lmd_2p(); lmd2++)
    {
        int dimXI = ppMap->at(lmd2).size();
        int dimMU = hhMap->at(lmd2).size();

        for (int xi_ab = 0; xi_ab < dimXI; xi_ab++)
            for (int mu_lj = 0; mu_lj < dimMU; mu_lj++)
            {
                int ab = ppMap->at(lmd2)(xi_ab);
                int a = ab % nP;
                int b = ab / nP;
                int lj = hhMap->at(lmd2)(mu_lj);
                int l = lj % nH;
                int j = lj / nH;

                int abj = a + (b + j * nP) * nP;
                int lmd_abj1 = map_p_p_h1_inv(0, abj);
                int idx_abj1 = map_p_p_h1_inv(1, abj);
                int lmd_l = map_h_inv(0, l);
                int idx_l = map_h_inv(1, l);
                if (lmd_abj1 != lmd_l)
                    throw string("Mismathcing lmd's in c_t2_t3, when filling t2old.");

                t2old_abj1_l.at(lmd_l)(idx_abj1, idx_l) = t2_old.at(lmd2)(xi_ab, mu_lj);
            }
    }
    for (int lmd2 = 0; lmd2 < basis->dim_lmd_2p(); lmd2++)
    {
        int dimXI = ppMap->at(lmd2).size();
        int dimMU = hhMap->at(lmd2).size();

        for (int xi_ab = 0; xi_ab < dimXI; xi_ab++)
            for (int mu_il = 0; mu_il < dimMU; mu_il++)
            {
                int ab = ppMap->at(lmd2)(xi_ab);
                int a = ab % nP;
                int b = ab / nP;
                int il = hhMap->at(lmd2)(mu_il);
                int i = il % nH;
                int l = il / nH;

                int abi = a + (b + i * nP) * nP;
                int lmd_abi1 = map_p_p_h1_inv(0, abi);
                int idx_abi1 = map_p_p_h1_inv(1, abi);
                int lmd_l = map_h_inv(0, l);
                int idx_l = map_h_inv(1, l);
                if (lmd_abi1 != lmd_l)
                    throw string("Mismathcing lmd's in c_t2_t3, when filling t2old.");

                t2old_abi1_l.at(lmd_l)(idx_abi1, idx_l) = t2_old.at(lmd2)(xi_ab, mu_il);
            }
    }
    //Fill i3_l_i
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        int dimH = map_h.at(lmd1).size();
        for (int idx_l = 0; idx_l < dimH; idx_l++)
            for (int idx_i = 0; idx_i < dimH; idx_i++)
            {
                int l = map_h.at(lmd1)(idx_l);
                int i = map_h.at(lmd1)(idx_i);
                i3_l_i.at(lmd1)(idx_l, idx_i) = i3(l, i);
            }
    }

    //Matmult
    vector<mat> t2new_abj1_i;
    vector<mat> t2new_abi1_j;
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        //        t2new_abj1_i.push_back(t2old_abj1_l.at(lmd1) * i3_l_i.at(lmd1));
        //        t2new_abi1_j.push_back(t2old_abi1_l.at(lmd1) * i3_l_j.at(lmd1));
        t2new_abj1_i.push_back(mult->dgemm(t2old_abj1_l.at(lmd1), i3_l_i.at(lmd1)));
        t2new_abi1_j.push_back(mult->dgemm(t2old_abi1_l.at(lmd1), i3_l_j.at(lmd1)));
    }

    //Fill t2_new
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        int dimPPH1 = map_p_p_h1.at(lmd1).size();
        int dimH = map_h.at(lmd1).size();

        for (int idx_abj1 = 0; idx_abj1 < dimPPH1; idx_abj1++)
            for (int idx_i = 0; idx_i < dimH; idx_i++)
            {
                int abj = map_p_p_h1.at(lmd1)(idx_abj1);
                int a = abj % nP;
                int bj = abj / nP;
                int b = bj % nP;
                int j = bj / nP;
                int i = map_h.at(lmd1)(idx_i);

                int ab = a + b * nP;
                int lmd_ab = (*ppMapInv)(0, ab);
                int xi_ab = (*ppMapInv)(1, ab);
                int ij = i + j * nH;
                int lmd_ij = (*hhMapInv)(0, ij);
                int mu_ij = (*hhMapInv)(1, ij);
                int ji = j + i * nH;
                int lmd_ji = (*hhMapInv)(0, ji);
                int mu_ji = (*hhMapInv)(1, ji);

                if (lmd_ab != lmd_ij)
                    throw string("Mismathcing lmd's in c_t2_t3, when filling t2new.");
                t2_new.at(lmd_ab)(xi_ab, mu_ij) -= t2new_abj1_i.at(lmd1)(idx_abj1, idx_i);
                if (lmd_ab != lmd_ji)
                    throw string("Mismathcing lmd's in c_t2_t3, when filling t2new.");
                t2_new.at(lmd_ab)(xi_ab, mu_ji) -= t2new_abi1_j.at(lmd1)(idx_abj1, idx_i);
            }
    }


    return;
}

void CC::c_t2_t5(vector<mat> &t2_new, mat const &i7, vector<mat> const &t2_old)
{
    Basis const * basis = sys->get_basis();
    umat const * ppMapInv = basis->get_map_de_lmdXI();
    umat const * phMapInv = basis->get_map_dl_lmdNU();
    umat const * hhMapInv = basis->get_map_lm_lmdMU();
    vector<uvec> const * ppMap = basis->get_map_lmdXI_de();
    vector<uvec> const * phMap = basis->get_map_lmdNU_dl();
    vector<uvec> const * hhMap = basis->get_map_lmdMU_lm();
    int nP = basis->get_nP();
    int nH = basis->get_nH();

    //Temporary storage
    vector<mat> i7_b_d;
    vector<mat> t2old_d_ija1;
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        int dimP = map_p.at(lmd1).size();
        int dimHHP1 = map_h_h_p1.at(lmd1).size();

        i7_b_d.push_back(zeros<mat > (dimP, dimP));
        t2old_d_ija1.push_back(zeros<mat > (dimP, dimHHP1));
    }

    //Fill t2old_d_ija1
    for (int lmd2 = 0; lmd2 < basis->dim_lmd_2p(); lmd2++)
    {
        int dimXI = ppMap->at(lmd2).size();
        int dimMU = hhMap->at(lmd2).size();

        for (int xi_ad = 0; xi_ad < dimXI; xi_ad++)
            for (int mu_ij = 0; mu_ij < dimMU; mu_ij++)
            {
                int ad = ppMap->at(lmd2)(xi_ad);
                int a = ad % nP;
                int d = ad / nP;
                int ij = hhMap->at(lmd2)(mu_ij);
                int i = ij % nH;
                int j = ij / nH;

                int ija = i + (j + a * nH) * nH;
                int lmd_ija1 = map_h_h_p1_inv(0, ija);
                int idx_ija1 = map_h_h_p1_inv(1, ija);
                int lmd_d = map_p_inv(0, d);
                int idx_d = map_p_inv(1, d);
                if (lmd_ija1 != lmd_d)
                    throw string("Mismathcing lmd's in c_t2_t5, when filling t2old.");
                t2old_d_ija1.at(lmd_d)(idx_d, idx_ija1) = t2_old.at(lmd2)(xi_ad, mu_ij);
            }
    }

    //Fill i7_b_d
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        int dimP = map_p.at(lmd1).size();
        for (int idx_b = 0; idx_b < dimP; idx_b++)
            for (int idx_d = 0; idx_d < dimP; idx_d++)
            {
                int b = map_p.at(lmd1)(idx_b);
                int d = map_p.at(lmd1)(idx_d);

                i7_b_d.at(lmd1)(idx_b, idx_d) = i7(b, d);
            }
    }

    //MAtmult
    vector<mat> t2new_b_ija1;
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
        t2new_b_ija1.push_back(mult->dgemm(i7_b_d.at(lmd1), t2old_d_ija1.at(lmd1)));

    //Map back into t2_new
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        int dimP = map_p.at(lmd1).size();
        int dimHHP1 = map_h_h_p1.at(lmd1).size();

        for (int idx_b = 0; idx_b < dimP; idx_b++)
            for (int idx_ija1 = 0; idx_ija1 < dimHHP1; idx_ija1++)
            {
                int b = map_p.at(lmd1)(idx_b);
                int ija = map_h_h_p1.at(lmd1)(idx_ija1);
                int i = ija % nH;
                int ja = ija / nH;
                int j = ja % nH;
                int a = ja / nH;

                int ab = a + b * nP;
                int lmd2_ab = (*ppMapInv)(0, ab);
                int xi_ab = (*ppMapInv)(1, ab);
                int ba = b + a * nP;
                int lmd2_ba = (*ppMapInv)(0, ba);
                int xi_ba = (*ppMapInv)(1, ba);
                int ij = i + j *nH;
                int lmd2_ij = (*hhMapInv)(0, ij);
                int mu_ij = (*hhMapInv)(1, ij);
                double t2abij_elem = t2new_b_ija1.at(lmd1)(idx_b, idx_ija1);

                if (lmd2_ab != lmd2_ij)
                    throw string("Mismathcing lmd's in c_t2_t5, when filling t2new(1).");
                t2_new.at(lmd2_ab)(xi_ab, mu_ij) += t2abij_elem;

                if (lmd2_ba != lmd2_ij)
                    throw string("Mismathcing lmd's in c_t2_t5, when filling t2new(2).");
                t2_new.at(lmd2_ba)(xi_ba, mu_ij) -= t2abij_elem;
            }
    }

    return;
}

void CC::c_t2_t6(vector<mat> &t2_new, vector<mat> const &i8, vector<mat> const &t2_old)
{
    Basis const * basis = sys->get_basis();
    umat const * ppMapInv = basis->get_map_de_lmdXI();
    umat const * phMapInv = basis->get_map_dl_lmdNU();
    umat const * hhMapInv = basis->get_map_lm_lmdMU();
    vector<uvec> const * ppMap = basis->get_map_lmdXI_de();
    vector<uvec> const * phMap = basis->get_map_lmdNU_dl();
    vector<uvec> const * hhMap = basis->get_map_lmdMU_lm();
    int nP = basis->get_nP();
    int nH = basis->get_nH();

    //Temporary storage
    vector<mat> i8_bj1_dl1;
    vector<mat> t2old_dl1_a1i;
    vector<mat> t2new_bj1_a1i;
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_2p(); lmd1++)
    {
        int dimPH1 = map_p_h1.at(lmd1).size();
        int dimP1H = map_p1_h.at(lmd1).size();

        t2old_dl1_a1i.push_back(zeros<mat > (dimPH1, dimP1H));
        i8_bj1_dl1.push_back(zeros<mat > (dimPH1, dimPH1));
    }
    //Copy into i8 and t2
    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        int dimMU = hhMap->at(lmd).size();
        int dimNU = phMap->at(lmd).size();
        int dimXI = ppMap->at(lmd).size();

        for (int nu_bl = 0; nu_bl < dimNU; nu_bl++)
            for (int nu_dj = 0; nu_dj < dimNU; nu_dj++)
            {
                int bl = phMap->at(lmd)(nu_bl);
                int b = bl % nP;
                int l = bl / nP;
                int dj = phMap->at(lmd)(nu_dj);
                int d = dj % nP;
                int j = dj / nP;

                int bj = b + j * nP;
                int lmd_b_j1 = map_p_h1_inv(0, bj);
                int b_j1 = map_p_h1_inv(1, bj);
                int dl = d + l * nP;
                int lmd_d_l1 = map_p_h1_inv(0, dl);
                int d_l1 = map_p_h1_inv(1, dl);
                if (lmd_b_j1 != lmd_d_l1)
                    throw string("Mismathcing lmd's in c_t2_t6, when filling i8.");

                i8_bj1_dl1.at(lmd_b_j1)(b_j1, d_l1) = i8.at(lmd)(nu_bl, nu_dj);
            }

        for (int xi_da = 0; xi_da < dimXI; xi_da++)
            for (int mu_li = 0; mu_li < dimMU; mu_li++)
            {
                int da = ppMap->at(lmd)(xi_da);
                int d = da % nP;
                int a = da / nP;
                int li = hhMap->at(lmd)(mu_li);
                int l = li % nH;
                int i = li / nH;

                int dl = d + l * nP;
                int lmd_d_l1 = map_p_h1_inv(0, dl);
                int d_l1 = map_p_h1_inv(1, dl);
                int ai = a + i * nP;
                int lmd_a1_i = map_p1_h_inv(0, ai);
                int a1_i = map_p1_h_inv(1, ai);
                if (lmd_d_l1 != lmd_a1_i)
                    throw string("Mismathcing lmd's in c_t2_t6, when filling t2_old.");

                t2old_dl1_a1i.at(lmd_d_l1)(d_l1, a1_i) = t2_old.at(lmd)(xi_da, mu_li);
            }
    }

    //Matrix mult 
    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
        t2new_bj1_a1i.push_back(mult->dgemm(i8_bj1_dl1.at(lmd), t2old_dl1_a1i.at(lmd)));

    //Map back into t2_new.
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_2p(); lmd1++)
    {
        int dim_PH1 = map_p_h1.at(lmd1).size();
        int dim_P1H = map_p1_h.at(lmd1).size();

        for (int b_j1 = 0; b_j1 < dim_PH1; b_j1++)
            for (int a1_i = 0; a1_i < dim_P1H; a1_i++)
            {
                int bj = map_p_h1.at(lmd1)(b_j1);
                int b = bj % nP;
                int j = bj / nP;
                int ai = map_p1_h.at(lmd1)(a1_i);
                int a = ai % nP;
                int i = ai / nP;

                int ab = a + b * nP;
                int lmd_ab = (*ppMapInv)(0, ab);
                int xi_ab = (*ppMapInv)(1, ab);
                int ij = i + j * nH;
                int lmd_ij = (*hhMapInv)(0, ij);
                int mu_ij = (*hhMapInv)(1, ij);
                int ba = b + a * nP;
                int lmd_ba = (*ppMapInv)(0, ba);
                int xi_ba = (*ppMapInv)(1, ba);
                int ji = j + i *nH;
                int lmd_ji = (*hhMapInv)(0, ji);
                int mu_ji = (*hhMapInv)(1, ji);
                double t2_bj1_a1i_elem = t2new_bj1_a1i.at(lmd1)(b_j1, a1_i);

                if (lmd_ab != lmd_ij)
                    throw string("Mismathcing lmd's in c_t2_t6, when filling t2_new(1).");
                t2_new.at(lmd_ab)(xi_ab, mu_ij) += t2_bj1_a1i_elem;

                if (lmd_ab != lmd_ji)
                    throw string("Mismathcing lmd's in c_t2_t6, when filling t2_new(2).");
                t2_new.at(lmd_ab)(xi_ab, mu_ji) -= t2_bj1_a1i_elem;

                if (lmd_ba != lmd_ij)
                    throw string("Mismathcing lmd's in c_t2_t6, when filling t2_new(3).");
                t2_new.at(lmd_ba) (xi_ba, mu_ij) -= t2_bj1_a1i_elem;

                if (lmd_ba != lmd_ji)
                    throw string("Mismathcing lmd's in c_t2_t6, when filling t2_new(4).");
                t2_new.at(lmd_ba)(xi_ba, mu_ji) += t2_bj1_a1i_elem;
            }
    }

    return;
}

void CC::c_t2_t7(vector<mat> &t2_new, vector<mat> const &i10, mat const &t1_old)
{
    Basis const * basis = sys->get_basis();
    umat const * ppMapInv = basis->get_map_de_lmdXI();
    umat const * phMapInv = basis->get_map_dl_lmdNU();
    umat const * hhMapInv = basis->get_map_lm_lmdMU();
    vector<uvec> const * ppMap = basis->get_map_lmdXI_de();
    vector<uvec> const * phMap = basis->get_map_lmdNU_dl();
    vector<uvec> const * hhMap = basis->get_map_lmdMU_lm();
    int nP = basis->get_nP();
    int nH = basis->get_nH();

    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        int dimXI = ppMap->at(lmd).size();

        for (int xi_ab = 0; xi_ab < dimXI; xi_ab++)
        {
            int ab = ppMap->at(lmd)(xi_ab);
            int a = ab % nP;
            int b = ab / nP;
            int ba = b + a * nP;
            if ((*ppMapInv)(0, ba) != lmd)
                throw string("Mismathing lmd's for ab and ba in c_t2_t7.");
            int xi_ba = (*ppMapInv)(1, ba);

            for (int l = 0; l < nH; l++)
            {
                int al = a + l * nP;
                int lmd_al = (*phMapInv)(0, al);
                int nu_al = (*phMapInv)(1, al);
                if (lmd_al == lmd)
                {
                    arma::mat tmp = t1_old(b, l) * i10.at(lmd)(nu_al, span::all);
                    t2_new.at(lmd)(xi_ab, span::all) -= tmp;
                    t2_new.at(lmd)(xi_ba, span::all) += tmp;
                }
            }
        }
    }

    return;
}

void CC::c_t2_t8(vector<mat> &t2_new, vector<mat> const &i11, mat const &t1_old)
{
    Basis const * basis = sys->get_basis();
    umat const * ppMapInv = basis->get_map_de_lmdXI();
    umat const * phMapInv = basis->get_map_dl_lmdNU();
    umat const * hhMapInv = basis->get_map_lm_lmdMU();
    vector<uvec> const * ppMap = basis->get_map_lmdXI_de();
    vector<uvec> const * phMap = basis->get_map_lmdNU_dl();
    vector<uvec> const * hhMap = basis->get_map_lmdMU_lm();
    int nP = basis->get_nP();
    int nH = basis->get_nH();

    //Temporary storage
    vector<mat> i11_abj1_d;
    vector<mat> t1_d_i;
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        int dimPPH1 = map_p_p_h1.at(lmd1).size();
        int dimP = map_p.at(lmd1).size();
        int dimH = map_h.at(lmd1).size();

        i11_abj1_d.push_back(zeros<mat > (dimPPH1, dimP));
        t1_d_i.push_back(zeros<mat > (dimP, dimH));
    }

    //Fill i11_abj1_d
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        int dimPPH1 = map_p_p_h1.at(lmd1).size();
        int dimP = map_p.at(lmd1).size();

        mat &i11_abj1_d_lmd1 = i11_abj1_d.at(lmd1);

        for (int idx_abj1 = 0; idx_abj1 < dimPPH1; idx_abj1++)
        {
            int abj = map_p_p_h1.at(lmd1)(idx_abj1);
            int a = abj % nP;
            int bj = abj / nP;
            int b = bj % nP;
            int j = bj / nP;
            int ab = a + b * nP;
            int lmd2_ab = (*ppMapInv)(0, ab);
            int xi_ab = (*ppMapInv)(1, ab);

            mat const &i11_lmd2 = i11.at(lmd2_ab);

            for (int idx_d = 0; idx_d < dimP; idx_d++)
            {
                int d = map_p.at(lmd1)(idx_d);
                int dj = d + j * nP;
                int lmd2_dj = (*phMapInv)(0, dj);
                int nu_dj = (*phMapInv)(1, dj);

                if (lmd2_ab != lmd2_dj)
                    throw string("Mismathcing lmd's in c_t2_t8, when filling i11.");
                i11_abj1_d_lmd1(idx_abj1, idx_d) = i11_lmd2(xi_ab, nu_dj);
            }
        }
    }

    //Fill t1_d_i
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        int dimP = map_p.at(lmd1).size();
        int dimH = map_h.at(lmd1).size();

        for (int idx_d = 0; idx_d < dimP; idx_d++)
            for (int idx_i = 0; idx_i < dimH; idx_i++)
            {
                int d = map_p.at(lmd1)(idx_d);
                int i = map_h.at(lmd1)(idx_i);
                t1_d_i.at(lmd1)(idx_d, idx_i) = t1_old(d, i);
            }
    }

    //Matrix multiplication
    vector<mat> t2new_abj1_i;
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
        t2new_abj1_i.push_back(mult->dgemm(i11_abj1_d.at(lmd1), t1_d_i.at(lmd1)));

    //Fill t2_new
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        int dimPPH1 = map_p_p_h1.at(lmd1).size();
        int dimH = map_h.at(lmd1).size();

        for (int idx_abj1 = 0; idx_abj1 < dimPPH1; idx_abj1++)
            for (int idx_i = 0; idx_i < dimH; idx_i++)
            {
                int abj = map_p_p_h1.at(lmd1)(idx_abj1);
                int a = abj % nP;
                int bj = abj / nP;
                int b = bj % nP;
                int j = bj / nP;
                int i = map_h.at(lmd1)(idx_i);

                int ab = a + b * nP;
                int lmd2_ab = (*ppMapInv)(0, ab);
                int xi_ab = (*ppMapInv)(1, ab);
                int ij = i + j * nH;
                int lmd2_ij = (*hhMapInv)(0, ij);
                int mu_ij = (*hhMapInv)(1, ij);
                int ji = j + i * nH;
                int lmd2_ji = (*hhMapInv)(0, ji);
                int mu_ji = (*hhMapInv)(1, ji);
                double t2_elem = t2new_abj1_i.at(lmd1)(idx_abj1, idx_i);

                if (lmd2_ab != lmd2_ij)
                    throw string("Mismathcing lmd's in c_t2_t8, when filling t2new(1).");
                t2_new.at(lmd2_ab)(xi_ab, mu_ij) += t2_elem;

                if (lmd2_ab != lmd2_ji)
                    throw string("Mismathcing lmd's in c_t2_t8, when filling t2new(2).");
                t2_new.at(lmd2_ab)(xi_ab, mu_ji) -= t2_elem;
            }
    }




    //    for (int a = 0; a < nP; a++)
    //        for (int b = 0; b < nP; b++)
    //        {
    //            int ab = a + b * nP;
    //            int lmd_ab = (*ppMapInv)(0, ab);
    //            int xi_ab = (*ppMapInv)(1, ab);
    //
    //            for (int i = 0; i < nH; i++)
    //                for (int j = 0; j < nH; j++)
    //                {
    //                    int ij = i + j * nH;
    //                    int lmd_ij = (*hhMapInv)(0, ij);
    //                    int mu_ij = (*hhMapInv)(1, ij);
    //
    //                    if (lmd_ab == lmd_ij)
    //                        for (int d = 0; d < nP; d++)
    //                        {
    //                            //+ i11^ab_dj t^d_i
    //                            int dj = d + j * nP;
    //                            int lmd_dj = (*phMapInv)(0, dj);
    //                            int nu_dj = (*phMapInv)(1, dj);
    //                            if (lmd_ab == lmd_dj)
    //                                t2_new.at(lmd_ab)(xi_ab, mu_ij) += i11.at(lmd_ab)(xi_ab, nu_dj) * t1_old(d, i);
    //
    //                            //- i11^ab_di t^d_j
    //                            int di = d + i *nP;
    //                            int lmd_di = (*phMapInv)(0, di);
    //                            int nu_di = (*phMapInv)(1, di);
    //                            if (lmd_ab == lmd_di)
    //                                t2_new.at(lmd_ab)(xi_ab, mu_ij) -= i11.at(lmd_ab)(xi_ab, nu_di) * t1_old(d, j);
    //                        }
    //                }
    //        }

    return;
}

void CC::reset_timing_info()
{
    //16 different timings
    this->timing_info = arma::field<arma::vec > (16);

    //Energy is counted as one term
    timing_info(0) = zeros<arma::vec > (1);

    //Intermediates
    timing_info(1) = zeros<arma::vec > (2);
    timing_info(2) = zeros<arma::vec > (1);
    timing_info(3) = zeros<arma::vec > (4);
    timing_info(4) = zeros<arma::vec > (2);
    timing_info(5) = zeros<arma::vec > (2);
    timing_info(6) = zeros<arma::vec > (3);
    timing_info(7) = zeros<arma::vec > (3);
    timing_info(8) = zeros<arma::vec > (4);
    timing_info(9) = zeros<arma::vec > (2);
    timing_info(10) = zeros<arma::vec > (1);
    timing_info(11) = zeros<arma::vec > (2);

    //T1  &  T2
    timing_info(12) = zeros<arma::vec > (7);
    timing_info(13) = zeros<arma::vec > (8);

    //D1  &  D2
    timing_info(14) = zeros<arma::vec > (1);
    timing_info(15) = zeros<arma::vec > (1);
}

void CC::init_additional_mappings()
{
    Basis const * basis = sys->get_basis();
    vector<uvec> const * mapHH = basis->get_map_lmdMU_lm();
    vector<uvec> const * mapPH = basis->get_map_lmdNU_dl();
    vector<uvec> const * mapPP = basis->get_map_lmdXI_de();
    umat const * mapPPinv = basis->get_map_de_lmdXI();
    umat const * mapPHinv = basis->get_map_dl_lmdNU();
    umat const * mapHHinv = basis->get_map_lm_lmdMU();
    int nP = basis->get_nP();
    int nH = basis->get_nH();

    //Create mappings pp^{-1}  (assuming lmd^-1 has same dimension as lmd)
    vector <vector<int > > p_p1_vec(basis->dim_lmd_2p());
    map_p_p1.clear();
    map_p_p1_inv = zeros<umat > (2, nP * nP);
    for (int pp = 0; pp < nP * nP; pp++)
    {
        int p1 = pp % nP;
        int p2 = pp / nP;

        int lmd_p_p1 = basis->lambda_neg(p1 + nH, p2 + nH);
        map_p_p1_inv(0, pp) = lmd_p_p1;
        map_p_p1_inv(1, pp) = p_p1_vec.at(lmd_p_p1).size();
        p_p1_vec.at(lmd_p_p1).push_back(pp);
    }
    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        uvec p_p1_uvec = zeros<uvec > (p_p1_vec.at(lmd).size());
        for (int idx = 0; idx < p_p1_uvec.size(); idx++)
            p_p1_uvec(idx) = p_p1_vec.at(lmd).at(idx);
        map_p_p1.push_back(p_p1_uvec);
    }

    //Create mappings ph^{-1} and p^{-1}h (assuming lmd^-1 has same dimension as lmd)
    vector <vector<int > > p_h1_vec(basis->dim_lmd_2p());
    map_p_h1.clear();
    map_p_h1_inv = zeros<umat > (2, nP * nH);
    vector <vector<int > > p1_h_vec(basis->dim_lmd_2p());
    map_p1_h.clear();
    map_p1_h_inv = zeros<umat > (2, nP * nH);
    for (int ph = 0; ph < nP * nH; ph++)
    {
        int p = ph % nP;
        int h = ph / nP;

        //p_h1
        int lmd_p_h1 = basis->lambda_neg(p + nH, h);
        map_p_h1_inv(0, ph) = lmd_p_h1;
        map_p_h1_inv(1, ph) = p_h1_vec.at(lmd_p_h1).size();
        p_h1_vec.at(lmd_p_h1).push_back(ph);

        //p1_h
        int lmd_p1_h = basis->lambda_neg(h, p + nH);
        map_p1_h_inv(0, ph) = lmd_p1_h;
        map_p1_h_inv(1, ph) = p1_h_vec.at(lmd_p1_h).size();
        p1_h_vec.at(lmd_p1_h).push_back(ph);
    }
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_2p(); lmd1++)
    {
        uvec p_h1_uvec = zeros<uvec > (p_h1_vec.at(lmd1).size());
        for (int nu1 = 0; nu1 < p_h1_uvec.size(); nu1++)
            p_h1_uvec(nu1) = p_h1_vec.at(lmd1).at(nu1);
        map_p_h1.push_back(p_h1_uvec);

        uvec p1_h_uvec = zeros<uvec > (p1_h_vec.at(lmd1).size());
        for (int nu1 = 0; nu1 < p1_h_uvec.size(); nu1++)
            p1_h_uvec(nu1) = p1_h_vec.at(lmd1).at(nu1);
        map_p1_h.push_back(p1_h_uvec);
    }

    //mapping for holes j
    vector<vector<int > > h_vec(basis->dim_lmd_1p());
    map_h.clear();
    map_h_inv = zeros<umat > (2, nH);
    for (int j = 0; j < nH; j++)
    {
        int lmd1 = basis->lambda_1p(j);
        map_h_inv(0, j) = lmd1;
        map_h_inv(1, j) = h_vec.at(lmd1).size();
        h_vec.at(lmd1).push_back(j);
    }
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        uvec uvecH(h_vec.at(lmd1).size());
        for (int idx = 0; idx < uvecH.size(); idx++)
            uvecH(idx) = h_vec.at(lmd1).at(idx);
        map_h.push_back(uvecH);
    }

    //Mapping for particles e
    vector<vector<int > > p_vec(basis->dim_lmd_1p());
    map_p.clear();
    map_p_inv = zeros<umat > (2, nP);
    for (int e = 0; e < nP; e++)
    {
        int lmd1 = basis->lambda_1p(e + nH);
        map_p_inv(0, e) = lmd1;
        map_p_inv(1, e) = p_vec.at(lmd1).size();
        p_vec.at(lmd1).push_back(e);
    }
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        uvec uvecP(p_vec.at(lmd1).size());
        for (int idx = 0; idx < uvecP.size(); idx++)
            uvecP(idx) = p_vec.at(lmd1).at(idx);
        map_p.push_back(uvecP);
    }

    //Create mappings p_p_p1 and h for <ab|x|dj> -> <abd1|x|j>
    //Also <ab|x|de> -> <abd1|x|e>
    //
    //Mappings for particles abd, relating them to same channel.
    vector<vector<int > > ppp1_vec(basis->dim_lmd_1p());
    map_p_p_p1.clear();
    map_p_p_p1_inv = 999999999 * ones<umat > (2, nP * nP * nP);
    for (int lmd2 = 0; lmd2 < basis->dim_lmd_2p(); lmd2++)
    {
        int dimXI = mapPP->at(lmd2).size();
        int dimNU = mapPH->at(lmd2).size();

        for (int xi_ab = 0; xi_ab < dimXI; xi_ab++)
        {
            int ab = mapPP->at(lmd2)(xi_ab);
            int a = ab % nP;
            int b = ab / nP;

            for (int nu_dj = 0; nu_dj < dimNU; nu_dj++)
            {
                int dj = mapPH->at(lmd2)(nu_dj);
                int d = dj % nP;
                int j = dj / nP;
                int abd = a + (b + d * nP) * nP;
                int lmd1 = map_h_inv(0, j);

                if (map_p_p_p1_inv(0, abd) != 999999999 || map_p_p_p1_inv(1, abd) != 999999999)
                {
                    if (map_p_p_p1_inv(0, abd) != lmd1)
                        cerr << "ERROR: MAPPING ALREADY DEFINED TO A DIFFERENT CHANNEL!\n";
                } else
                {
                    map_p_p_p1_inv(0, abd) = lmd1;
                    map_p_p_p1_inv(1, abd) = ppp1_vec.at(lmd1).size();
                    ppp1_vec.at(lmd1).push_back(abd);
                }
            }

            for (int xi_de = 0; xi_de < dimXI; xi_de++)
            {
                int de = mapPP->at(lmd2)(xi_de);
                int d = de % nP;
                int e = de / nP;
                int abd = a + (b + d * nP) * nP;
                int lmd1 = map_p_inv(0, e);

                if (map_p_p_p1_inv(0, abd) != 999999999 || map_p_p_p1_inv(1, abd) != 999999999)
                {
                    if (map_p_p_p1_inv(0, abd) != lmd1)
                        cerr << "ERROR: MAPPING ALREADY DEFINED TO A DIFFERENT CHANNEL!\n";
                } else
                {
                    map_p_p_p1_inv(0, abd) = lmd1;
                    map_p_p_p1_inv(1, abd) = ppp1_vec.at(lmd1).size();
                    ppp1_vec.at(lmd1).push_back(abd);
                }
            }
        }
    }
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        uvec uvecPPP1(ppp1_vec.at(lmd1).size());
        for (int idx = 0; idx < uvecPPP1.size(); idx++)
            uvecPPP1(idx) = ppp1_vec.at(lmd1).at(idx);
        map_p_p_p1.push_back(uvecPPP1);
    }

    //Mapping <abj^-1|x|l>
    //
    //p_p_h1
    vector<vector<int > > pph1_vec(basis->dim_lmd_1p());
    map_p_p_h1.clear();
    map_p_p_h1_inv = 999999999 * ones<umat > (2, nP * nP * nH);
    for (int lmd2 = 0; lmd2 < basis->dim_lmd_2p(); lmd2++)
    {
        int dimXI = mapPP->at(lmd2).size();
        int dimMU = mapHH->at(lmd2).size();
        int dimNU = mapPH->at(lmd2).size();

        for (int xi_ab = 0; xi_ab < dimXI; xi_ab++)
        {
            int ab = mapPP->at(lmd2)(xi_ab);
            int a = ab % nP;
            int b = ab / nP;

            for (int mu_jl = 0; mu_jl < dimMU; mu_jl++)
            {
                int jl = mapHH->at(lmd2)(mu_jl);
                int j = jl % nH;
                int l = jl / nH;
                int abj = a + (b + j * nP) * nP;
                int lmd1 = map_h_inv(0, l);

                if (map_p_p_h1_inv(0, abj) != 999999999 || map_p_p_h1_inv(1, abj) != 999999999)
                {
                    if (map_p_p_h1_inv(0, abj) != lmd1)
                        cerr << "ERROR: MAPPING ALREADY DEFINED TO A DIFFERENT CHANNEL!\n";
                } else
                {
                    map_p_p_h1_inv(0, abj) = lmd1;
                    map_p_p_h1_inv(1, abj) = pph1_vec.at(lmd1).size();
                    pph1_vec.at(lmd1).push_back(abj);
                }
            }

            for (int nu_ej = 0; nu_ej < dimNU; nu_ej++)
            {
                int ej = mapPH->at(lmd2)(nu_ej);
                int e = ej % nP;
                int j = ej / nP;
                int abj = a + (b + j * nP) * nP;
                int lmd1 = map_p_inv(0, e);

                if (map_p_p_h1_inv(0, abj) != 999999999 || map_p_p_h1_inv(1, abj) != 999999999)
                {
                    if (map_p_p_h1_inv(0, abj) != lmd1)
                        cerr << "ERROR: MAPPING ALREADY DEFINED TO A DIFFERENT CHANNEL!\n";
                } else
                {
                    map_p_p_h1_inv(0, abj) = lmd1;
                    map_p_p_h1_inv(1, abj) = pph1_vec.at(lmd1).size();
                    pph1_vec.at(lmd1).push_back(abj);
                }
            }
        }
    }
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        uvec uvecPPH1(pph1_vec.at(lmd1).size());
        for (int idx = 0; idx < uvecPPH1.size(); idx++)
            uvecPPH1(idx) = pph1_vec.at(lmd1).at(idx);
        map_p_p_h1.push_back(uvecPPH1);
    }

    //Mapping h_h_p1 for <ad|t2|ij> -> <d|t2|ija^-1>
    vector<vector<int > > hhp1_vec(basis->dim_lmd_1p());
    map_h_h_p1.clear();
    map_h_h_p1_inv = 999999999 * ones<umat > (2, nP * nH * nH);
    for (int lmd2 = 0; lmd2 < basis->dim_lmd_2p(); lmd2++)
    {
        int dimXI = mapPP->at(lmd2).size();
        int dimMU = mapHH->at(lmd2).size();
        int dimNU = mapPH->at(lmd2).size();

        for (int mu_ij = 0; mu_ij < dimMU; mu_ij++)
        {
            int ij = mapHH->at(lmd2)(mu_ij);
            int i = ij % nH;
            int j = ij / nH;

            for (int xi_ad = 0; xi_ad < dimXI; xi_ad++)
            {
                int ad = mapPP->at(lmd2)(xi_ad);
                int a = ad % nP;
                int d = ad / nP;

                int ija = i + (j + a * nH) * nH;
                int lmd1 = map_p_inv(0, d);

                if (map_h_h_p1_inv(0, ija) != 999999999 || map_h_h_p1_inv(1, ija) != 999999999)
                {
                    if (map_h_h_p1_inv(0, ija) != lmd1)
                        cerr << "ERROR: MAPPING ALREADY DEFINED TO A DIFFERENT CHANNEL!\n";
                } else
                {
                    map_h_h_p1_inv(0, ija) = lmd1;
                    map_h_h_p1_inv(1, ija) = hhp1_vec.at(lmd1).size();
                    hhp1_vec.at(lmd1).push_back(ija);
                }
            }

            for (int nu_al = 0; nu_al < dimNU; nu_al++)
            {
                int al = mapPH->at(lmd2)(nu_al);
                int a = al % nP;
                int l = al / nP;
                int ija = i + (j + a * nH) * nH;
                int lmd1 = map_h_inv(0, l);

                if (map_h_h_p1_inv(0, ija) != 999999999 || map_h_h_p1_inv(1, ija) != 999999999)
                {
                    if (map_h_h_p1_inv(0, ija) != lmd1)
                        cerr << "ERROR: MAPPING ALREADY DEFINED TO A DIFFERENT CHANNEL!\n";
                } else
                {
                    map_h_h_p1_inv(0, ija) = lmd1;
                    map_h_h_p1_inv(1, ija) = hhp1_vec.at(lmd1).size();
                    hhp1_vec.at(lmd1).push_back(ija);
                }
            }
        }
    }
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        uvec uvecHHP1(hhp1_vec.at(lmd1).size());
        for (int idx = 0; idx < uvecHHP1.size(); idx++)
            uvecHHP1(idx) = hhp1_vec.at(lmd1).at(idx);
        map_h_h_p1.push_back(uvecHHP1);
    }


    //Map interaction elements <d1_l||e_m1>
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_2p(); lmd1++)
    {
        int dimP_H1 = map_p_h1.at(lmd1).size();
        int dimP1_H = map_p1_h.at(lmd1).size();

        v_d1l_em1.push_back(zeros<mat > (dimP1_H, dimP_H1));
    }
    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        int dimMU = mapHH->at(lmd).size();
        int dimXI = mapPP->at(lmd).size();
        for (int xi_ed = 0; xi_ed < dimXI; xi_ed++)
            for (int mu_ml = 0; mu_ml < dimMU; mu_ml++)
            {
                int ed = mapPP->at(lmd)(xi_ed);
                int e = ed % nP;
                int d = ed / nP;
                int ml = mapHH->at(lmd)(mu_ml);
                int m = ml % nH;
                int l = ml / nH;

                int dl = d + l * nP;
                int em = e + m * nP;
                int lmd1_d1_l = map_p1_h_inv(0, dl);
                int d1_l = map_p1_h_inv(1, dl);
                int lmd1_e_m1 = map_p_h1_inv(0, em);
                int e_m1 = map_p_h1_inv(1, em);

                if (lmd1_e_m1 != lmd1_d1_l)
                    cerr << "WRONG lmd1's!!!\n";
                v_d1l_em1.at(lmd1_d1_l)(d1_l, e_m1) = sys->get_v_pphh()->at(lmd)(xi_ed, mu_ml);
            }
    }

    //Interaction elements <a_b_d1||e>
    v_abd1_e.clear();
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        int dimE = map_p.at(lmd1).size();
        int dimABD1 = map_p_p_p1.at(lmd1).size();
        //        cout << "lmd1: " << lmd1 << "  SIZES: " << dimABD1 << "," << dimE << endl;
        v_abd1_e.push_back(zeros<mat > (dimABD1, dimE));
    }
    for (int lmd2 = 0; lmd2 < basis->dim_lmd_2p(); lmd2++)
    {
        int dimXI = mapPP->at(lmd2).size();

        for (int xi_ab = 0; xi_ab < dimXI; xi_ab++)
            for (int xi_de = 0; xi_de < dimXI; xi_de++)
            {
                int ab = mapPP->at(lmd2)(xi_ab);
                int a = ab % nP;
                int b = ab / nP;
                int de = mapPP->at(lmd2)(xi_de);
                int d = de % nP;
                int e = de / nP;

                int abd = a + (b + d * nP) * nP;
                int lmd1 = map_p_inv(0, e);
                int e_idx = map_p_inv(1, e);
                int abd_idx = map_p_p_p1_inv(1, abd);
                if (map_p_p_p1_inv(0, abd) != lmd1)
                    cerr << "lmd1's does not match!\n";
                v_abd1_e.at(lmd1)(abd_idx, e_idx) = sys->get_v_pppp()->at(lmd2)(xi_ab, xi_de);
            }
    }

    //Interaction elements <ad1||el1>
    v_ad1_el1.clear();
    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        int dimPP1 = map_p_p1.at(lmd).size();
        int dimPH1 = map_p_h1.at(lmd).size();
        v_ad1_el1.push_back(zeros<mat > (dimPP1, dimPH1));
    }
    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        int dimXI = mapPP->at(lmd).size();
        int dimNU = mapPH->at(lmd).size();

        for (int xi_de = 0; xi_de < dimXI; xi_de++)
            for (int nu_al = 0; nu_al < dimNU; nu_al++)
            {
                int de = mapPP->at(lmd)(xi_de);
                int d = de % nP;
                int e = de / nP;
                int al = mapPH->at(lmd)(nu_al);
                int a = al % nP;
                int l = al / nP;

                int ad = a + d * nP;
                int lmd_ad1 = map_p_p1_inv(0, ad);
                int idx_ad1 = map_p_p1_inv(1, ad);
                int el = e + l * nP;
                int lmd_el1 = map_p_h1_inv(0, el);
                int idx_el1 = map_p_h1_inv(1, el);
                if (lmd_ad1 != lmd_el1)
                    throw std::string("Mismatching lmd's when creating additional matrix elements <ad1||el1>.");

                double elm = sys->get_v_ppph()->at(lmd)(xi_de, nu_al);
                v_ad1_el1.at(lmd_ad1)(idx_ad1, idx_el1) = elm;
            }
    }

    //Interaction elements <d||lme1>
    v_d_lme1.clear();
    for (int lmd1 = 0; lmd1 < basis->dim_lmd_1p(); lmd1++)
    {
        int dimP = map_p.at(lmd1).size();
        int dimHHP1 = map_h_h_p1.at(lmd1).size();
        v_d_lme1.push_back(zeros<mat > (dimP, dimHHP1));

        for (int idx_d = 0; idx_d < dimP; idx_d++)
            for (int idx_lme1 = 0; idx_lme1 < dimHHP1; idx_lme1++)
            {
                int d = map_p.at(lmd1)(idx_d);
                int lme = map_h_h_p1.at(lmd1)(idx_lme1);
                int l = lme % nH;
                int me = lme / nH;
                int m = me % nH;
                int e = me / nH;

                int de = d + e * nP;
                int lmd2_de = (*mapPPinv)(0, de);
                int xi_de = (*mapPPinv)(1, de);
                int lm = l + m * nH;
                int lmd2_lm = (*mapHHinv)(0, lm);
                int mu_lm = (*mapHHinv)(1, lm);
                if (lmd2_de != lmd2_lm)
                    throw std::string("Mismatching lmd's when creating additional matrix elements <d||lme1>");

                double elm = sys->get_v_pphh()->at(lmd2_de)(xi_de, mu_lm);
                v_d_lme1.at(lmd1)(idx_d, idx_lme1) = elm;
            }
    }

    return;
}
