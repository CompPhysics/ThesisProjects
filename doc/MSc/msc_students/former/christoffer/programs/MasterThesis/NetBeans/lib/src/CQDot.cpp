/**
 * @file   CQDot.cpp
 * @author toffyrn
 * 
 * Created on 19. desember 2011, 15:05
 */

#include "CQDot.h"

using namespace toffyrn::libUNK;
using namespace arma;
using namespace std;

CQDot::CQDot(int filledR, int totalR, std::string fname, bool nmms_file)
{
    //Initialize
    basis = new HOBasis(filledR, totalR);
    filename = fname;
    fillMatrixElements(nmms_file);
}

CQDot::~CQDot()
{
}

void Generate::genFile(int shells, std::string fname)
{
    HOBasis basis(shells, shells); //All shells are filled, but it doesn't matter since only statemap is needed.
    std::ofstream file(fname.c_str(), ios::binary);

    int states = shells * (shells + 1);
    for (unsigned short p = 0; p < states; p++)
        for (unsigned short q = 0; q < p; q++)
        {
            char buffer[50];
            double douP = p;
            sprintf(buffer, "\rWriting tpElem: %.2f %%", (100.0 * (douP / states + q / (douP * states))));
            cout << buffer;

            for (unsigned short r = 0; r < states; r++)
                for (unsigned short s = 0; s < r; s++)
                {
                    ivec3 pM = basis.stateMap(p);
                    ivec3 qM = basis.stateMap(q);
                    ivec3 rM = basis.stateMap(r);
                    ivec3 sM = basis.stateMap(s);

                    double directT;
                    if (pM(2) == rM(2) && qM(2) == sM(2))
                        directT = coulomb(pM(0), pM(1), qM(0), qM(1), sM(0), sM(1), rM(0), rM(1));
                    else
                        directT = 0.0;
                    double exchangeT;
                    if (pM(2) == sM(2) && qM(2) == rM(2))
                        exchangeT = coulomb(pM(0), pM(1), qM(0), qM(1), rM(0), rM(1), sM(0), sM(1));
                    else
                        exchangeT = 0.0;

                    double element = directT - exchangeT;

                    if (element != 0.0)
                    {
                        file.write((char*) &p, sizeof (p));
                        file.write((char*) &q, sizeof (q));
                        file.write((char*) &r, sizeof (r));
                        file.write((char*) &s, sizeof (s));
                        file.write((char*) &element, sizeof (element));
                    }
                }
        }

    file.close();
}

void CQDot::fillMatrixElements(bool nmms_file)
{
    int nH = basis->get_nH();
    int nP = basis->get_nP();

    //Prepare interaction matrices.
    v_hhhh.clear();
    v_phhh.clear();
    v_pphh.clear();
    v_phph.clear();
    v_ppph.clear();
    v_pppp.clear();
    for (int lmd = 0; lmd < basis->dim_lmd_2p(); lmd++)
    {
        int dimMU = basis->get_map_lmdMU_lm()->at(lmd).size();
        int dimNU = basis->get_map_lmdNU_dl()->at(lmd).size();
        int dimXI = basis->get_map_lmdXI_de()->at(lmd).size();

        v_hhhh.push_back(zeros<mat > (dimMU, dimMU));
        v_phhh.push_back(zeros<mat > (dimNU, dimMU));
        v_pphh.push_back(zeros<mat > (dimXI, dimMU));
        v_phph.push_back(zeros<mat > (dimNU, dimNU));
        v_ppph.push_back(zeros<mat > (dimXI, dimNU));
        v_pppp.push_back(zeros<mat > (dimXI, dimXI));
    }


    //Open file
    std::ifstream file(filename.c_str(), ios::binary);
    if (!file.is_open())
        throw std::string("Could not open file.");

    if (nmms_file)
        readNMMSfile(file);
    else
    {
        //Loop through file
        while (!file.eof())
        {
            unsigned short p, q, r, s;
            double element;
            file.read((char*) &p, sizeof (p));
            file.read((char*) &q, sizeof (q));
            file.read((char*) &r, sizeof (r));
            file.read((char*) &s, sizeof (s));
            file.read((char*) &element, sizeof (element));

            putOneElement(p, q, r, s, element);
            putOneElement(p, q, s, r, -element);
            putOneElement(q, p, s, r, element);
            putOneElement(q, p, r, s, -element);
        }
    }
    //(and close it)
    file.close();

    //Generate f_matrices
    f_hh = zeros<mat > (nH, nH);
    for (int i = 0; i < nH; i++)
        for (int j = 0; j < nH; j++)
            f_hh(i, j) = f_elem(i, j);
    f_ph = zeros<mat > (nP, nH);
    for (int a = 0; a < nP; a++)
        for (int i = 0; i < nH; i++)
            f_ph(a, i) = f_elem(a + nH, i);
    f_pp = zeros<mat > (nP, nP);
    for (int a = 0; a < nP; a++)
        for (int b = 0; b < nP; b++)
            f_pp(a, b) = f_elem(a + nH, b + nH);
}

void CQDot::readNMMSfile(std::ifstream &file)
{
    unsigned short n_a, n_b, n_c, n_d;
    short m_a, m_b, m_c, m_d;
    char ms_a, ms_b, ms_c, ms_d;
    double element;

    //Creating reversed state map
    map<unsigned int, unsigned int> revMap;
    int ntot = basis->get_nH() + basis->get_nP();
    int maxM = 40000; //40000 M-values * 40000 N-values * 2 spin values fits in uint.
    int hMaxM = maxM / 2; //Used to center the m values.
    for (int state = 0; state < ntot; state++)
    {
        ivec3 nmms = basis->stateMap(state);
        unsigned int n = nmms(0);
        int m = nmms(1);
        unsigned int ms = nmms(2);
        unsigned int key = n + (m + (ms + hMaxM) * maxM) * maxM;
        revMap[key] = state;
    }

    while (!file.eof())
    {
        file.read((char*) &n_a, sizeof (n_a));
        file.read((char*) &m_a, sizeof (m_a));
        file.read((char*) &ms_a, sizeof (ms_a));
        file.read((char*) &n_b, sizeof (n_b));
        file.read((char*) &m_b, sizeof (m_b));
        file.read((char*) &ms_b, sizeof (ms_b));
        file.read((char*) &n_c, sizeof (n_c));
        file.read((char*) &m_c, sizeof (m_c));
        file.read((char*) &ms_c, sizeof (ms_c));
        file.read((char*) &n_d, sizeof (n_d));
        file.read((char*) &m_d, sizeof (m_d));
        file.read((char*) &ms_d, sizeof (ms_d));
        file.read((char*) &element, sizeof (element));
        if (file.eof())
            break;

        if (ms_a == -1)
            ms_a = 1;
        else if (ms_a == +1)
            ms_a = 0;
        unsigned int key = n_a + (m_a + (ms_a + hMaxM) * maxM) * maxM;
        map<unsigned int, unsigned int>::iterator found = revMap.find(key);
        if (found == revMap.end())
            continue;
        int p = found->second;

        if (ms_b == -1)
            ms_b = 1;
        else if (ms_b == +1)
            ms_b = 0;
        key = n_b + (m_b + (ms_b + hMaxM) * maxM) * maxM;
        found = revMap.find(key);
        if (found == revMap.end())
            continue;
        int q = found->second;

        if (ms_c == -1)
            ms_c = 1;
        else if (ms_c == +1)
            ms_c = 0;
        key = n_c + (m_c + (ms_c + hMaxM) * maxM) * maxM;
        found = revMap.find(key);
        if (found == revMap.end())
            continue;
        int r = found->second;

        if (ms_d == -1)
            ms_d = 1;
        else if (ms_d == +1)
            ms_d = 0;
        key = n_d + (m_d + (ms_d + hMaxM) * maxM) * maxM;
        found = revMap.find(key);
        if (found == revMap.end())
            continue;
        int s = found->second;

        putOneElement(p, q, r, s, element);
        putOneElement(p, q, s, r, -element);
        putOneElement(q, p, s, r, element);
        putOneElement(q, p, r, s, -element);
    }

    return;
}

void CQDot::putOneElement(int p, int q, int r, int s, double element)
{
    if (max(max(p, q), max(r, s)) >= (basis->get_nH() + basis->get_nP()))
        return;

    int nH = basis->get_nH();
    int nP = basis->get_nP();

    umat const * hhmap = basis->get_map_lm_lmdMU();
    umat const * phmap = basis->get_map_dl_lmdNU();
    umat const * ppmap = basis->get_map_de_lmdXI();

    if (p < nH && q < nH && r < nH && s < nH)
    { //hhhh
        int pq = p + q * nH;
        int lmd_pq = (*hhmap)(0, pq);
        int mu_pq = (*hhmap)(1, pq);
        int rs = r + s* nH;
        int lmd_rs = (*hhmap)(0, rs);
        int mu_rs = (*hhmap)(1, rs);
        if (lmd_pq != lmd_rs)
            throw std::string("Transition between different channels non-zero.");
        v_hhhh.at(lmd_pq)(mu_pq, mu_rs) = element;
    } else
        if (p >= nH && q < nH && r < nH && s < nH)
    { //phhh
        int pq = (p - nH) + q * nP;
        int lmd_pq = (*phmap)(0, pq);
        int nu_pq = (*phmap)(1, pq);
        int rs = r + s * nH;
        int lmd_rs = (*hhmap)(0, rs);
        int mu_rs = (*hhmap)(1, rs);
        if (lmd_pq != lmd_rs)
            throw std::string("Transition between different channels non-zero.");
        v_phhh.at(lmd_pq)(nu_pq, mu_rs) = element;
    } else
        if (p >= nH && q >= nH && r < nH && s < nH)
    { //pphh
        int pq = (p - nH) + (q - nH) * nP;
        int lmd_pq = (*ppmap)(0, pq);
        int xi_pq = (*ppmap)(1, pq);
        int rs = r + s * nH;
        int lmd_rs = (*hhmap)(0, rs);
        int mu_rs = (*hhmap)(1, rs);
        if (lmd_pq != lmd_rs)
            throw std::string("Transition between different channels non-zero.");
        v_pphh.at(lmd_pq)(xi_pq, mu_rs) = element;
    } else
        if (p >= nH && q < nH && r >= nH && s < nH)
    { //phph
        int pq = (p - nH) + q * nP;
        int lmd_pq = (*phmap)(0, pq);
        int nu_pq = (*phmap)(1, pq);
        int rs = (r - nH) + s * nP;
        int lmd_rs = (*phmap)(0, rs);
        int nu_rs = (*phmap)(1, rs);
        if (lmd_pq != lmd_rs)
            throw std::string("Transition between different channels non-zero.");
        v_phph.at(lmd_pq)(nu_pq, nu_rs) = element;
    } else
        if (p >= nH && q >= nH && r >= nH && s < nH)
    { //ppph
        int pq = (p - nH) + (q - nH) * nP;
        int lmd_pq = (*ppmap)(0, pq);
        int xi_pq = (*ppmap)(1, pq);
        int rs = (r - nH) + s * nP;
        int lmd_rs = (*phmap)(0, rs);
        int nu_rs = (*phmap)(1, rs);
        if (lmd_pq != lmd_rs)
            throw std::string("Transition between different channels non-zero.");
        v_ppph.at(lmd_pq)(xi_pq, nu_rs) = element;
    } else
        if (p >= nH && q >= nH && r >= nH && s >= nH)
    { //pppp
        int pq = (p - nH) + (q - nH) * nP;
        int lmd_pq = (*ppmap)(0, pq);
        int xi_pq = (*ppmap)(1, pq);
        int rs = (r - nH) + (s - nH) * nP;
        int lmd_rs = (*ppmap)(0, rs);
        int xi_rs = (*ppmap)(1, rs);
        if (lmd_pq != lmd_rs)
            throw std::string("Transition between different channels non-zero.");
        v_pppp.at(lmd_pq)(xi_pq, xi_rs) = element;
    }
}

double CQDot::f_elem(std::size_t p, std::size_t q) const
{
    double u_pq = 0.0;
    for (int i = 0; i < basis->get_nH(); i++)
        u_pq += v_elem(p, i, q, i);

    if (p != q)
        return u_pq;

    ivec3 pM = basis->stateMap(p);
    int n = pM(0);
    int m = pM(1);

    //TODO: Do not hardcode omega;
    double omega = 1.0;
    return omega * (2 * n + abs(m) + 1) + u_pq;
}

double CQDot::v_elem(
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

double Generate::coulomb(const int n1, const int m1, const int n2, const int m2, const int n3, const int m3, const int n4, const int m4) const
{
    const double * lgv = lgammaVec.memptr();
    const double * lgf = logfacVec.memptr();

    double coulombint = 0.0;
    if (m1 + m2 == m3 + m4)
    {
        double temp;

        int gamma1 = 0;
        int gamma2 = 0;
        int gamma3 = 0;
        int gamma4 = 0;
        int G = 0;
        for (int j1 = 0; j1 <= n1; j1++)
            for (int j2 = 0; j2 <= n2; j2++)
                for (int j3 = 0; j3 <= n3; j3++)
                    for (int j4 = 0; j4 <= n4; j4++)
                    {
                        gamma1 = (int) (j1 + j4 + 0.5 * (abs(m1) + m1) + 0.5 * (abs(m4) - m4));
                        gamma2 = (int) (j2 + j3 + 0.5 * (abs(m2) + m2) + 0.5 * (abs(m3) - m3));
                        gamma3 = (int) (j2 + j3 + 0.5 * (abs(m3) + m3) + 0.5 * (abs(m2) - m2));
                        gamma4 = (int) (j1 + j4 + 0.5 * (abs(m4) + m4) + 0.5 * (abs(m1) - m1));

                        G = gamma1 + gamma2 + gamma3 + gamma4;

                        double Lgratio1 = LogRatio1(j1, j2, j3, j4);
                        double Lgproduct2 = LogProduct2(n1, m1, n2, m2, n3, m3, n4, m4, j1, j2, j3, j4);
                        double Lgratio2 = LogRatio2(G);

                        double lgprd3_sum1 = lgf[gamma1] + lgf[gamma2] + lgf[gamma3] + lgf[gamma4];

                        temp = 0.0;
                        for (int l1 = 0; l1 <= gamma1; l1++)
                        {
                            double lgprd3_sum2 = lgprd3_sum1 - lgf[l1] - lgf[gamma1 - l1];
                            for (int l2 = 0; l2 <= gamma2; l2++)
                            {
                                double lgprd3_sum3 = lgprd3_sum2 - lgf[l2] - lgf[gamma2 - l2];
                                for (int l3 = 0; l3 <= gamma3; l3++)
                                    for (int l4 = 0; l4 <= gamma4; l4++)
                                    {
                                        if ((l1 + l2) == (l3 + l4))
                                        {
                                            int lambda = l1 + l2 + l3 + l4;
                                            double lgprd3 = lgprd3_sum3 - lgf[l3] - lgf[l4] - lgf[gamma3 - l3] - lgf[gamma4 - l4];
                                            if ((gamma2 + gamma3 - l2 - l3) % 2 == 0)
                                                temp += exp(lgprd3 + lgv[2 + lambda] + lgv[G - lambda + 1]);
                                            else
                                                temp -= exp(lgprd3 + lgv[2 + lambda] + lgv[G - lambda + 1]);
                                        }
                                    }
                            }
                        }


                        coulombint += minusPower(j1 + j2 + j3 + j4) * exp(Lgratio1 + Lgproduct2 + Lgratio2) * temp;

                    }
        coulombint *= Product1(n1, m1, n2, m2, n3, m3, n4, m4);
    }
    return coulombint;
}

/* Computes the first ratio in the Asinimovas expression
 */
double Generate::LogRatio1(const int j1, const int j2, const int j3, const int j4) const
{
    double temp = -LogFac(j1) - LogFac(j2) - LogFac(j3) - LogFac(j4);
    return temp;
}

/* Computes the 2nd ratio in the Asinimovas expression
 */
double Generate::LogRatio2(const int G) const
{
    double temp = -1 * (G + 1)*0.5 * log(2);
    return temp;
}

/* Computes the log of the 2nd product in the Asinimovas expression
 */
double Generate::LogProduct2(const int n1, const int m1, const int n2, const int m2, const int n3, const int m3, const int n4, const int m4, const int j1, const int j2, const int j3, const int j4) const
{
    double temp = LogFac(n1 + abs(m1)) + LogFac(n2 + abs(m2)) + LogFac(n3 + abs(m3)) + LogFac(n4 + abs(m4)) - LogFac(n1 - j1) - LogFac(n2 - j2) - LogFac(n3 - j3) - LogFac(n4 - j4) - LogFac(j1 + abs(m1)) - LogFac(j2 + abs(m2)) - LogFac(j3 + abs(m3)) - LogFac(j4 + abs(m4));
    return temp;
}

/* Computes the log of the 3rd product in the Asinimovas expression
 */
double Generate::LogProduct3(const int l1, const int l2, const int l3, const int l4, const int gamma1, const int gamma2, const int gamma3, const int gamma4) const
{
    const double * lgf = logfacVec.memptr();
    double temp = lgf[gamma1] + lgf[gamma2] + lgf[gamma3] + lgf[gamma4] - lgf[l1] - lgf[l2] - lgf[l3] - lgf[l4] - lgf[gamma1 - l1] - lgf[gamma2 - l2] - lgf[gamma3 - l3] - lgf[gamma4 - l4];
    return temp;
}

/* computes (-1)^k
 */
int Generate::minusPower(const int k) const
{
    int temp = abs(k % 2)*-2 + 1; // gives 1 if k is even, gives -1 if k is odd
    return temp;
}

/* computes log(n!)
 */
double Generate::LogFac(const int n, bool calc)
{
    if (n > 170)
    {
        printf("#### Too big integer in LogFac(n)!!!!#### \n\n");
        exit(1);
    }

    if (calc == true)
    {
        logfacVec = vec(n + 1);
        logfacVec(0) = 0;
        logfacVec(1) = 0;
        //(The first two elements are zero)

        for (int elem = 2; elem <= n; elem++)
            logfacVec(elem) = log(elem) + logfacVec(elem - 1);
    }

    return logfacVec(n);
}

/* computes log(n!)
 */
double Generate::LogFac(const int n) const
{
    if (n > 170)
    {
        printf("#### Too big integer in LogFac(n)!!!!#### \n\n");
        exit(1);
    }

    return logfacVec(n);
}

/* computes first product of indices in the Anisimovas/Matulis expression
 */
double Generate::Product1(const int n1, const int ml1, const int n2, const int ml2, const int n3, const int ml3, const int n4, const int ml4) const
{
    double temp = 0;
    temp = LogFac(n1) + LogFac(n2) + LogFac(n3) + LogFac(n4) - LogFac(n1 + abs(ml1)) - LogFac(n2 + abs(ml2)) - LogFac(n3 + abs(ml3)) - LogFac(n4 + abs(ml4));
    temp *= 0.5;
    return exp(temp);
}

/* lgamma.cpp -- log gamma function of real argument.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
//  Returns log(gamma) of real argument.
//  NOTE: Returns 1e308 if argument is 0 or negative.
//  taken on the web [http://www.crbond.com/math.htm]
 */
double Generate::lgamma(int xTimes2) const
{
#if CHECK
    if (x >= 171)
    {
        printf("#### Too big integer in lgamma(x) to give accurate result!!!!#### \n\n");
        exit(1);
    }
#endif

    return lgammaVec(xTimes2);
}// end of lgamma function

double Generate::lgamma(int xTimes2, bool calc)
{
#if CHECK
    if (x >= 171)
    {
        printf("#### Too big integer in lgamma(x) to give accurate result!!!!#### \n\n");
        exit(1);
    }
#endif
    if (calc)
    {
        int xTimes2max = xTimes2;
        lgammaVec = zeros<vec > (xTimes2max + 1);
        for (xTimes2 = 0; xTimes2 <= xTimes2max; xTimes2++)
        {
            double x = xTimes2 / 2.0;

            double x0, x2, xp, gl, gl0;
            int n, k;
            static double a[] = {
                8.333333333333333e-02,
                -2.777777777777778e-03,
                7.936507936507937e-04,
                -5.952380952380952e-04,
                8.417508417508418e-04,
                -1.917526917526918e-03,
                6.410256410256410e-03,
                -2.955065359477124e-02,
                1.796443723688307e-01,
                -1.39243221690590
            };

            x0 = x;
            if (x <= 0.0)
            {
                lgammaVec(xTimes2) = 1e308;
                continue;
            } else if ((x == 1.0) || (x == 2.0))
            {
                lgammaVec(xTimes2) = 0.0;
                continue;
            } else if (x <= 7.0)
            {
                n = (int) (7 - x);
                x0 = x + n;
            }
            x2 = 1.0 / (x0 * x0);
            xp = 2.0 * M_PI;
            gl0 = a[9];
            for (k = 8; k >= 0; k--)
            {
                gl0 = gl0 * x2 + a[k];
            }
            gl = gl0 / x0 + 0.5 * log(xp)+(x0 - 0.5) * log(x0) - x0;
            if (x <= 7.0)
            {
                for (k = 1; k <= n; k++)
                {
                    gl -= log(x0 - 1.0);
                    x0 -= 1.0;
                }
            }
            lgammaVec(xTimes2) = gl;
        }

        xTimes2 = xTimes2max;
    }

    return lgammaVec(xTimes2);
}// end of lgamma function

Generate::Generate()
{
    //Fill the needed vectors for coulomb potential
    LogFac(169, true);
    lgamma(340, true);

    //Fill lgprd3_comp
    lgprd3_compon = zeros<mat > (169, 169);
    for (int l = 0; l < 169; l++)
        for (int gamma = 0; gamma < 169; gamma++)
            if (gamma - l >= 0)
                lgprd3_compon(gamma, l) = logfacVec(gamma) - logfacVec(l) - logfacVec(gamma - l);
}
