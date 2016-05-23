/**
 * @file   CCSimple.cpp
 * @author toffyrn
 * 
 * Created on 27. desember 2011, 09:10
 */

#include "CCSimple.h"

using namespace arma;
using namespace std;

double toffyrn::libUNK::CCSimple::testSolver(System* sys, int maxiter, int debugFlags)
{
    int nP = sys->get_basis()->get_nP();
    int nH = sys->get_basis()->get_nH();
    int nT = nH + nP;

    //Find reference energy.
    double E_ref = 0;
    for (int k = 0; k < nH; k++)
    {
        E_ref += sys->f_elem(k, k);
        for (int l = 0; l < nH; l++)
            E_ref -= 0.5 * sys->v_elem(k, l, k, l);
    }
    if (debugFlags & DEBUG_ENERGY)
        cout << "E_ref: " << E_ref << endl;

    //Set up amplitudes
    //t_i^a
    mat t1_new = zeros<mat > (nH, nP);
    mat t1_old = t1_new;
    //t_ij^ab
    field<mat> t2_new(nH, nH);
    for (int i = 0; i < nH; i++)
        for (int j = 0; j < nH; j++)
            t2_new(i, j) = zeros<mat > (nP, nP);
    field<mat> t2_old = t2_new;

    //Find the "initial" energy
    double E_ccsd = 0;
    for (int i = 0; i < nH; i++)
        for (int a = 0; a < nP; a++)
            E_ccsd += sys->f_elem(i, a + nH) * t1_new(i, a);
    for (int i = 0; i < nH; i++)
        for (int j = 0; j < nH; j++)
            for (int a = 0; a < nP; a++)
                for (int b = 0; b < nP; b++)
                {
                    E_ccsd += 0.25 * sys->v_elem(i, j, a + nH, b + nH) * t2_new(i, j)(a, b);
                    E_ccsd += 0.5 * sys->v_elem(i, j, a + nH, b + nH) * t1_new(i, a) * t1_new(j, b);
                }
    if (debugFlags & DEBUG_ENERGY)
        cout << "E_ccsd: " << E_ccsd << endl;

    double E_new = E_ref;
    double E_old = E_new;

    //Setting up arrays for interactions
    field<mat> hhhh(nH, nH);
    for (int i = 0; i < nH; i++)
        for (int j = 0; j < nH; j++)
        {
            hhhh(i, j) = zeros<mat > (nH, nH);
            for (int k = 0; k < nH; k++)
                for (int l = 0; l < nH; l++)
                    hhhh(i, j)(k, l) = sys->v_elem(i, j, k, l);
        }
    if (debugFlags & DEBUG_INTERACTIONS)
        cout << "hhhh:\n" << hhhh << endl;
    //
    field<mat> phhh(nP, nH);
    for (int a = 0; a < nP; a++)
        for (int i = 0; i < nH; i++)
        {
            phhh(a, i) = zeros<mat > (nH, nH);
            for (int j = 0; j < nH; j++)
                for (int k = 0; k < nH; k++)
                    phhh(a, i)(j, k) = sys->v_elem(a + nH, i, j, k);
        }
    if (debugFlags & DEBUG_INTERACTIONS)
        cout << "phhh:\n" << phhh << endl;
    //
    field<mat> pphh(nP, nP);
    for (int a = 0; a < nP; a++)
        for (int b = 0; b < nP; b++)
        {
            pphh(a, b) = zeros<mat > (nH, nH);
            for (int j = 0; j < nH; j++)
                for (int k = 0; k < nH; k++)
                    pphh(a, b)(j, k) = sys->v_elem(a + nH, b + nH, j, k);
        }
    if (debugFlags & DEBUG_INTERACTIONS)
        cout << "pphh:\n" << pphh << endl;
    //
    field<mat> phph(nP, nH);
    for (int a = 0; a < nP; a++)
        for (int i = 0; i < nH; i++)
        {
            phph(a, i) = zeros<mat > (nP, nH);
            for (int b = 0; b < nP; b++)
                for (int k = 0; k < nH; k++)
                    phph(a, i)(b, k) = sys->v_elem(a + nH, i, b + nH, k);
        }
    if (debugFlags & DEBUG_INTERACTIONS)
        cout << "phph:\n" << phph << endl;
    //
    field<mat> ppph(nP, nP);
    for (int a = 0; a < nP; a++)
        for (int b = 0; b < nP; b++)
        {
            ppph(a, b) = zeros<mat > (nP, nH);
            for (int c = 0; c < nP; c++)
                for (int k = 0; k < nH; k++)
                    ppph(a, b)(c, k) = sys->v_elem(a + nH, b + nH, c + nH, k);
        }
    if (debugFlags & DEBUG_INTERACTIONS)
        cout << "ppph:\n" << ppph << endl;
    //
    field<mat> pppp(nP, nP);
    for (int a = 0; a < nP; a++)
        for (int b = 0; b < nP; b++)
        {
            pppp(a, b) = zeros<mat > (nP, nP);
            for (int c = 0; c < nP; c++)
                for (int d = 0; d < nP; d++)
                    pppp(a, b)(c, d) = sys->v_elem(a + nH, b + nH, c + nH, d + nH);
        }
    if (debugFlags & DEBUG_INTERACTIONS)
        cout << "pppp:\n" << pppp << endl;
    //
    mat pp(nP, nP);
    for (int a = 0; a < nP; a++)
        for (int b = 0; b < nP; b++)
            pp(a, b) = sys->f_elem(a + nH, b + nH);
    if (debugFlags & DEBUG_INTERACTIONS)
        cout << "pp:\n" << pp << endl;
    //
    mat ph(nP, nH);
    for (int a = 0; a < nP; a++)
        for (int i = 0; i < nH; i++)
            ph(a, i) = sys->f_elem(a + nH, i);
    if (debugFlags & DEBUG_INTERACTIONS)
        cout << "ph:\n" << ph << endl;
    //
    mat hp(nH, nP);
    for (int a = 0; a < nP; a++)
        for (int i = 0; i < nH; i++)
            hp(i, a) = sys->f_elem(i, a + nH);
    if (debugFlags & DEBUG_INTERACTIONS)
        cout << "hp:\n" << hp << endl;
    //
    mat hh(nH, nH);
    for (int i = 0; i < nH; i++)
        for (int j = 0; j < nH; j++)
            hh(i, j) = sys->f_elem(i, j);
    if (debugFlags & DEBUG_INTERACTIONS)
        cout << "hh:\n" << hh << endl;

    //Iterate
    for (int iter = 0; iter < maxiter; iter++)
    {
        if (debugFlags != DEBUG_NONE)
            cout << "============ ITERATION " << iter << " ==============\n";

        //Save last iteration
        t1_old = t1_new;
        t2_old = t2_new;

        //Intermediates
        mat i1(nP, nP);
        for (int a = 0; a < nP; a++)
            for (int d = 0; d < nP; d++)
            {
                i1(a, d) = pp(a, d);
                for (int l = 0; l < nH; l++)
                    for (int e = 0; e < nP; e++)
                        i1(a, d) += ppph(d, e)(a, l) * t1_old(l, e);
            }
        if (debugFlags & DEBUG_I1)
            cout << "i1:\n" << i1 << endl;
        //
        mat i2(nH, nP);
        for (int l = 0; l < nH; l++)
            for (int d = 0; d < nP; d++)
            {
                i2(l, d) = hp(l, d);
                for (int m = 0; m < nH; m++)
                    for (int e = 0; e < nP; e++)
                        i2(l, d) += pphh(d, e)(l, m) * t1_old(m, e);
            }
        if (debugFlags & DEBUG_I2)
            cout << "i2:\n" << i2 << endl;
        //
        mat i3(nH, nH);
        for (int l = 0; l < nH; l++)
            for (int i = 0; i < nH; i++)
            {
                i3(l, i) = hh(l, i);
                for (int m = 0; m < nH; m++)
                    for (int d = 0; d < nP; d++)
                        i3(l, i) += phhh(d, i)(m, l) * t1_old(m, d);
                for (int m = 0; m < nH; m++)
                    for (int d = 0; d < nP; d++)
                        for (int e = 0; e < nP; e++)
                            i3(l, i) += 0.5 * pphh(d, e)(l, m) * t2_old(i, m)(d, e);
                for (int d = 0; d < nP; d++)
                    i3(l, i) += i2(l, d) * t1_old(i, d);
            }
        if (debugFlags & DEBUG_I3)
            cout << "i3:\n" << i3 << endl;
        //
        field<mat> i5(nH, nP);
        for (int i = 0; i < nH; i++)
            for (int d = 0; d < nP; d++)
            {
                i5(i, d) = zeros<mat > (nH, nH);
                for (int l = 0; l < nH; l++)
                    for (int m = 0; m < nH; m++)
                    {
                        i5(i, d)(l, m) = -phhh(d, i)(l, m);
                        for (int e = 0; e < nP; e++)
                            i5(i, d)(l, m) += 0.5 * pphh(e, d)(l, m) * t1_old(i, e);
                    }
            }
        if (debugFlags & DEBUG_I5)
            cout << "i5:\n" << i5 << endl;
        //
        field<mat> i4(nH, nP);
        for (int i = 0; i < nH; i++)
            for (int d = 0; d < nP; d++)
            {
                i4(i, d) = zeros<mat > (nH, nH);
                for (int l = 0; l < nH; l++)
                    for (int m = 0; m < nH; m++)
                    {
                        i4(i, d)(l, m) = i5(i, d)(l, m);
                        for (int e = 0; e < nP; e++)
                            i4(i, d)(l, m) += 0.5 * pphh(e, d)(l, m) * t1_old(i, e);
                    }
            }
        if (debugFlags & DEBUG_I4)
            cout << "i4:\n" << i4 << endl;
        //            }
        mat d1(nH, nP);
        for (int i = 0; i < nH; i++)
            for (int a = 0; a < nP; a++)
            {
                //d1(i, a) = phph(a, i)(a, i) - i1(a, a) + i3(i, i);
                d1(i, a) = -i1(a, a) + i3(i, i);
            }
        if (debugFlags & DEBUG_D1)
            cout << "D1:\n" << d1 << endl;
        //
        field<mat> i6(nH, nH);
        for (int i = 0; i < nH; i++)
            for (int j = 0; j < nH; j++)
            {
                i6(i, j) = zeros<mat > (nH, nH);
                for (int l = 0; l < nH; l++)
                    for (int m = 0; m < nH; m++)
                    {
                        i6(i, j)(l, m) = hhhh(l, m)(i, j);
                        for (int d = 0; d < nP; d++)
                            for (int e = 0; e < nP; e++)
                                i6(i, j)(l, m) += 0.5 * pphh(d, e)(l, m) * t2_old(i, j)(d, e);
                        for (int d = 0; d < nP; d++)
                            i6(i, j)(l, m) += i5(i, d)(l, m) * t1_old(j, d) - i5(j, d)(l, m) * t1_old(i, d);
                    }
            }
        if (debugFlags & DEBUG_I6)
            cout << "i6:\n" << i6 << endl;
        //
        mat i7(nP, nP);
        for (int b = 0; b < nP; b++)
            for (int d = 0; d < nP; d++)
            {
                i7(b, d) = i1(b, d);
                for (int l = 0; l < nH; l++)
                    i7(b, d) -= i2(l, d) * t1_old(l, b);
                for (int l = 0; l < nH; l++)
                    for (int m = 0; m < nH; m++)
                        for (int e = 0; e < nP; e++)
                            i7(b, d) += 0.5 * pphh(d, e)(l, m) * t2_old(l, m)(e, b);
            }
        if (debugFlags & DEBUG_I7)
            cout << "i7:\n" << i7 << endl;
        //
        field<mat> i9(nP, nH);
        for (int d = 0; d < nP; d++)
            for (int j = 0; j < nH; j++)
            {
                i9(d, j) = zeros<mat > (nH, nP);
                for (int l = 0; l < nH; l++)
                    for (int b = 0; b < nP; b++)
                    {
                        i9(d, j)(l, b) = -phph(b, l)(d, j);
                        for (int e = 0; e < nP; e++)
                            i9(d, j)(l, b) += 0.5 * ppph(e, d)(b, l) * t1_old(j, e);
                    }
            }
        if (debugFlags & DEBUG_I9)
            cout << "i9\n" << i9 << endl;
        //
        field<mat> i8(nP, nH);
        for (int d = 0; d < nP; d++)
            for (int j = 0; j < nH; j++)
            {
                i8(d, j) = zeros<mat > (nH, nP);
                for (int l = 0; l < nH; l++)
                    for (int b = 0; b < nP; b++)
                    {
                        i8(d, j)(l, b) = i9(d, j)(l, b);
                        for (int e = 0; e < nP; e++)
                            i8(d, j)(l, b) += 0.5 * ppph(e, d)(b, l) * t1_old(j, e);
                        for (int m = 0; m < nH; m++)
                            i8(d, j)(l, b) += i4(j, d)(l, m) * t1_old(m, b);
                        for (int m = 0; m < nH; m++)
                            for (int e = 0; e < nP; e++)
                                i8(d, j)(l, b) += 0.5 * pphh(d, e)(l, m) * t2_old(m, j)(e, b);
                    }
            }
        if (debugFlags & DEBUG_I8)
            cout << "i8\n" << i8 << endl;
        //
        field<mat> i10(nH, nH);
        for (int i = 0; i < nH; i++)
            for (int j = 0; j < nH; j++)
            {
                i10(i, j) = zeros<mat > (nP, nH);
                for (int a = 0; a < nP; a++)
                    for (int l = 0; l < nH; l++)
                    {
                        i10(i, j)(a, l) = phhh(a, l)(i, j);
                        for (int d = 0; d < nP; d++)
                            for (int e = 0; e < nP; e++)
                                i10(i, j)(a, l) += 0.5 * ppph(d, e)(a, l) * t2_old(i, j)(d, e);
                        for (int d = 0; d < nP; d++)
                            i10(i, j)(a, l) += i9(d, i)(l, a) * t1_old(j, d) - i9(d, j)(l, a) * t1_old(i, d);
                        for (int m = 0; m < nH; m++)
                            i10(i, j)(a, l) += 0.5 * i6(i, j)(l, m) * t1_old(m, a);
                    }
            }
        if (debugFlags & DEBUG_I10)
            cout << "i10\n" << i10 << endl;
        //
        field<mat> i11(nP, nH);
        for (int d = 0; d < nP; d++)
            for (int j = 0; j < nH; j++)
            {
                i11(d, j) = zeros<mat > (nP, nP);
                for (int a = 0; a < nP; a++)
                    for (int b = 0; b < nP; b++)
                    {
                        i11(d, j)(a, b) = ppph(a, b)(d, j);
                        for (int e = 0; e < nP; e++)
                            i11(d, j)(a, b) += 0.5 * pppp(a, b)(d, e) * t1_old(j, e);
                    }
            }
        if (debugFlags & DEBUG_I11)
            cout << "i11\n" << i11 << endl;
        //
        field<mat> d2(nH, nH);
        for (int i = 0; i < nH; i++)
            for (int j = 0; j < nH; j++)
            {
                d2(i, j) = zeros<mat > (nP, nP);
                for (int a = 0; a < nP; a++)
                    for (int b = 0; b < nP; b++)
                    {
                        d2(i, j)(a, b) = -0.5 * pppp(a, b)(a, b);
                        d2(i, j)(a, b) += i3(i, i) + i3(j, j);
                        d2(i, j)(a, b) += -0.5 * i6(i, j)(i, j);
                        d2(i, j)(a, b) += -i7(b, b) - i7(a, a);
                        //                        d2(i, j)(a, b) += -i8(b, j)(j, b) - i8(a, j)(j, a) - i8(b, i)(i, b) - i8(a, i)(i, a);
                    }
            }
        if (debugFlags & DEBUG_D2)
            cout << "D2: \n" << d2 << endl;

        //T1 amplitudes
        for (int i = 0; i < nH; i++)
            for (int a = 0; a < nP; a++)
            {
                t1_new(i, a) = ph(a, i);
                for (int l = 0; l < nH; l++)
                    for (int d = 0; d < nP; d++)
                        //if (l != i || d != a)
                        t1_new(i, a) -= phph(a, l)(d, i) * t1_old(l, d);
                for (int l = 0; l < nH; l++)
                    for (int d = 0; d < nP; d++)
                        for (int e = 0; e < nP; e++)
                            t1_new(i, a) += 0.5 * ppph(d, e)(a, l) * t2_old(i, l)(d, e);
                for (int d = 0; d < nP; d++)
                    if (d != a)
                        t1_new(i, a) += i1(a, d) * t1_old(i, d);
                for (int l = 0; l < nH; l++)
                    if (l != i)
                        t1_new(i, a) -= i3(l, i) * t1_old(l, a);
                for (int l = 0; l < nH; l++)
                    for (int m = 0; m < nH; m++)
                        for (int d = 0; d < nP; d++)
                            t1_new(i, a) += 0.5 * i4(i, d)(l, m) * t2_old(l, m)(d, a);
                for (int l = 0; l < nH; l++)
                    for (int d = 0; d < nP; d++)
                        t1_new(i, a) += i2(l, d) * t2_old(i, l)(a, d);
                //Divide
                t1_new(i, a) = t1_new(i, a) / d1(i, a);
            }
        if (debugFlags & DEBUG_T1)
            cout << "T1:\n" << t1_new << endl;


        //Find t2 amplitudes
        for (int i = 0; i < nH; i++)
            for (int j = 0; j < nH; j++)
                for (int a = 0; a < nP; a++)
                    for (int b = 0; b < nP; b++)
                    {
                        t2_new(i, j)(a, b) = pphh(a, b)(i, j);
                        for (int d = 0; d < nP; d++)
                            for (int e = 0; e < nP; e++)
                                if (d != a || e != b)
                                    t2_new(i, j)(a, b) += 0.5 * pppp(a, b)(d, e) * t2_old(i, j)(d, e);
                        //if (i == 0 && j == 0 && a == 0 && b == 0)
                        //cout << "1: t2(0,0)(0,0) " << t2_new(0, 0)(0, 0) << endl;
                        for (int l = 0; l < nH; l++)
                            if (l != i)
                                t2_new(i, j)(a, b) += -i3(l, i) * t2_old(l, j)(a, b);
                        for (int l = 0; l < nH; l++)
                            if (l != j)
                                t2_new(i, j)(a, b) += -i3(l, j) * t2_old(i, l)(a, b);
                        //if (i == 0 && j == 0 && a == 0 && b == 0)
                        //  cout << "2: t2(0,0)(0,0) " << t2_new(0, 0)(0, 0) << endl;
                        for (int l = 0; l < nH; l++)
                            for (int m = 0; m < nH; m++)
                                if (l != i || m != j)
                                    t2_new(i, j)(a, b) += 0.5 * t2_old(l, m)(a, b) * i6(i, j)(l, m);
                        //if (i == 0 && j == 0 && a == 0 && b == 0)
                        //  cout << "3: t2(0,0)(0,0) " << t2_new(0, 0)(0, 0) << endl;
                        for (int d = 0; d < nP; d++)
                            if (d != b)
                                t2_new(i, j)(a, b) += t2_old(i, j)(a, d) * i7(b, d);
                        for (int d = 0; d < nP; d++)
                            if (d != a)
                                t2_new(i, j)(a, b) += t2_old(i, j)(d, b) * i7(a, d);
                        //if (i == 0 && j == 0 && a == 0 && b == 0)
                        //  cout << "4: t2(0,0)(0,0) " << t2_new(0, 0)(0, 0) << endl;
                        for (int l = 0; l < nH; l++)
                            for (int d = 0; d < nP; d++)
                                //if (l != j || d != b)
                                t2_new(i, j)(a, b) += t2_old(i, l)(a, d) * i8(d, j)(l, b);
                        //if (i == 0 && j == 0 && a == 0 && b == 0)
                        //  cout << "5: t2(0,0)(0,0) " << t2_new(0, 0)(0, 0) << endl;
                        for (int l = 0; l < nH; l++)
                            for (int d = 0; d < nP; d++)
                                //if (l != j || d != a)
                                t2_new(i, j)(a, b) += t2_old(i, l)(d, b) * i8(d, j)(l, a);
                        //if (i == 0 && j == 0 && a == 0 && b == 0)
                        //  cout << "5: t2(0,0)(0,0) " << t2_new(0, 0)(0, 0) << endl;
                        for (int l = 0; l < nH; l++)
                            for (int d = 0; d < nP; d++)
                                //if (l != i || d != b)
                                t2_new(i, j)(a, b) += t2_old(l, j)(a, d) * i8(d, i)(l, b);
                        //if (i == 0 && j == 0 && a == 0 && b == 0)
                        //  cout << "5: t2(0,0)(0,0) " << t2_new(0, 0)(0, 0) << endl;
                        for (int l = 0; l < nH; l++)
                            for (int d = 0; d < nP; d++)
                                //if (l != i || d != a)
                                t2_new(i, j)(a, b) += t2_old(l, j)(d, b) * i8(d, i)(l, a);
                        //if (i == 0 && j == 0 && a == 0 && b == 0)
                        //  cout << "5: t2(0,0)(0,0) " << t2_new(0, 0)(0, 0) << endl;
                        for (int l = 0; l < nH; l++)
                            t2_new(i, j)(a, b) += t1_old(l, a) * i10(i, j)(b, l) - t1_old(l, b) * i10(i, j)(a, l);
                        //if (i == 0 && j == 0 && a == 0 && b == 0)
                        //  cout << "6: t2(0,0)(0,0) " << t2_new(0, 0)(0, 0) << endl;
                        for (int d = 0; d < nP; d++)
                            t2_new(i, j)(a, b) += t1_old(i, d) * i11(d, j)(a, b) - t1_old(j, d) * i11(d, i)(a, b);
                        //if (i == 0 && j == 0 && a == 0 && b == 0)
                        //  cout << "7: t2(0,0)(0,0) " << t2_new(0, 0)(0, 0) << endl;

                        //Divide by denominator
                        t2_new(i, j)(a, b) = t2_new(i, j)(a, b) / d2(i, j)(a, b);
                    }
        if (debugFlags & DEBUG_T2)
            cout << "T2:\n" << t2_new << endl;


        E_ccsd = 0;
        for (int i = 0; i < nH; i++)
            for (int a = 0; a < nP; a++)
                E_ccsd += hp(i, a) * t1_new(i, a);
        for (int i = 0; i < nH; i++)
            for (int j = 0; j < nH; j++)
                for (int a = 0; a < nP; a++)
                    for (int b = 0; b < nP; b++)
                    {
                        E_ccsd += 0.25 * pphh(a, b)(i, j) * t2_new(i, j)(a, b);
                        E_ccsd += 0.5 * pphh(a, b)(i, j) * t1_new(i, a) * t1_new(j, b);
                    }
        if (debugFlags & DEBUG_ENERGY)
            cout << "E_ccsd: " << setprecision(10) << E_ccsd << endl;

        //        {
        //            //Test the amplitudes.
        //            double summaT = 0;
        //            for (int i = 0; i < nH; i++)
        //                for (int a = 0; a < nP; a++)
        //                {
        //                    double summa = sys->f_elem(a + nH, i);
        //                    for (int l = 0; l < nH; l++)
        //                    {
        //                        summa -= i3(l, i) * t1_new(l, a);
        //                        for (int d = 0; d < nP; d++)
        //                        {
        //                            summa += sys->v_elem(l, a + nH, d + nH, i) * t1_new(l, d);
        //                            for (int e = 0; e < nP; e++)
        //                                summa += 0.5 * sys->v_elem(a + nH, l, d + nH, e + nH) * t2_new(i, l)(d, e);
        //                            for (int m = 0; m < nH; m++)
        //                                summa += 0.5 * i4(i, d)(l, m) * t2_new(l, m)(d, a);
        //                            summa += i2(l, d) * t2_new(i, l)(a, d);
        //                        }
        //                    }
        //                    for (int d = 0; d < nP; d++)
        //                        summa += i1(a, d) * t1_new(i, d);
        //                    summaT += abs(summa);
        //                }
        //        }

    }

    if (debugFlags & DEBUG_ENERGY)
        cout << "E_ref:  " << E_ref << endl;

    E_ccsd = 0;
    for (int i = 0; i < nH; i++)
        for (int a = 0; a < nP; a++)
            E_ccsd += hp(i, a) * t1_new(i, a);
    for (int i = 0; i < nH; i++)
        for (int j = 0; j < nH; j++)
            for (int a = 0; a < nP; a++)
                for (int b = 0; b < nP; b++)
                {
                    E_ccsd += 0.25 * pphh(a, b)(i, j) * t2_new(i, j)(a, b);
                    E_ccsd += 0.5 * pphh(a, b)(i, j) * t1_new(i, a) * t1_new(j, b);
                }
    if (debugFlags & DEBUG_ENERGY)
    {
        cout << "E_ccsd: " << setprecision(10) << E_ccsd << endl;
        cout << "E_tot:  " << setprecision(10) << E_ref + E_ccsd << endl;
    }

    return E_ref + E_ccsd;
}

