/* 
 * File:   main.cpp
 * Author: chrishir
 *
 * Created on 6. desember 2011, 14:01
 */

#include <cstdlib>
#include "System.h"
#include "Basis.h"
#include "HF.h"
#include "Atoms.h"
#include "CC.h"
#include "CQDot.h"
#include "CCSimple.h"
#include "CLgemm.h"
#include "HFbasis.h"


using namespace std;
using namespace toffyrn::libUNK;
using namespace arma;

/*
 * 
 */
int main(int argc, char** argv)
{
    bool testMATMULT = false;

    if (testMATMULT)
    {
        int m = atoi(argv[1]);
        int n = atoi(argv[2]);
        int p = atoi(argv[3]);
        int times = atoi(argv[4]);
        cout << "m: " << m << "  n: " << n << "  p: " << p << "  times: " << times << endl;
        CLgemm mult;

        mat C_amd = ones<mat > (m, n);
        mat C_arm = ones<mat > (m, n);
        mat A(m, p);
        A.randn();
        mat B(p, n);
        B.randn();

        wall_clock timer;
        timer.tic();
        for (int it = 0; it < times; it++)
            C_arm += 0.5 * A * B;
        double timespent = timer.toc();
        cout << "Time spent ARMA: " << timespent << endl;

        timer.tic();
        for (int it = 0; it < times; it++)
            mult.dgemm(C_amd, A, B, 0.5, 1);
        timespent = timer.toc();
        cout << "Time spent O_CL: " << timespent << endl;

        mat diff = abs(C_amd - C_arm);
        cout << "MAXDIFF: " << max(max(diff)) << endl;

    } else
    {

        //    Generate gen;
        //    gen.genFile(8, "tpElem8.dat");

        //sys
        //CQDot sys(3, 14, "tpElem15.dat", false);
        CQDot sys(3, 18, "tp_w=1_R=20_nmms.simen.dat", true);
        //Atoms sys;

        //        HF hf(0.00000001, 100);
        //        hf.set_system(&sys);
        //        pair<double, mat> energy_coeff = hf.solve_ground_state_energy();
        //        cout << std::setprecision(10) << energy_coeff.first << endl;
        //        HFbasis hfb(&sys, energy_coeff.second);

        System * chosen_sys = &sys;
        //System * chosen_sys = &hfb;

        CC ccsd(1e-7, 100);
        ccsd.set_system(chosen_sys);
        ccsd.set_mult(new CLgemm());

        if (argc == 2 && strcmp(argv[1], "test") == 0)
        {
            CCSimple ccsi;
            cout << std::setprecision(10) << ccsi.testSolver(chosen_sys, 40, CCSimple::DEBUG_ENERGY) << endl;
        } else
        {
            cout << std::setprecision(10) << ccsd.solve_ground_state_energy(CC::DEBUG_ENERGY) << endl;
            //cout << ccsd.get_timing_info() << endl;
        }

        cout << "Time spent in GEMM class: " << ccsd.get_mult()->get_tot_time() << endl;
    }

    return 0;
}

