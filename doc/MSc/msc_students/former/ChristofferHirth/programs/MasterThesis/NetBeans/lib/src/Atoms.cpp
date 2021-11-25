/* 
 * File:   Atoms.cpp
 * Author: chrishir
 * 
 * Created on 7. desember 2011, 12:49
 */

#include "Atoms.h"

using namespace std;
using namespace arma;
using namespace toffyrn::libUNK;

Atoms::Atoms(int numElec)
{
    //Charge from nucleus defaults to number of electrons.
    Z = numElec;

    //There are only 8 states implemented in spatial_integral().
    basis = new Basis(numElec, 8 - numElec);

    //Fill matrices
    fillMatrixElements();
}

Atoms::~Atoms()
{
}

double Atoms::f_elem(std::size_t p, std::size_t q) const
{
    double u_pq = 0;

    for (int i = 0; i < basis->get_nH(); i++)
        u_pq += v_elem(p, i, q, i);

    if (p == q)
    {
        int n = p / 2 + 1;
        return -pow((double) Z, 2) / (2 * n * n) + u_pq;
    } else
        return u_pq;
}

double Atoms::v_elem(std::size_t p, std::size_t q, std::size_t r, std::size_t s) const
{
    //If r == s the element is zero
    //<pq||rr> = <pq|v|rr> - <pq|v|rr> = 0
    if (r == s)
        return 0;

    //Get spin of particles
    bool pUP = p % 2 == 0;
    bool qUP = q % 2 == 0;
    bool rUP = r % 2 == 0;
    bool sUP = s % 2 == 0;

    //Map p,q,r,s to n quantum number
    //0,1 -> 1
    //2,3 -> 2
    //...
    p = p / 2 + 1;
    q = q / 2 + 1;
    r = r / 2 + 1;
    s = s / 2 + 1;

    //First term of element
    double term1;
    if (pUP == rUP && qUP == sUP)
        term1 = spatial_integral(p, q, r, s);
    else
        term1 = 0;

    //Second term
    double term2;
    if (pUP == sUP && qUP == rUP)
        term2 = spatial_integral(p, q, s, r);
    else
        term2 = 0;

    return term1 - term2;
}

double Atoms::spatial_integral(int p, int q, int r, int s) const
{
    double E;

    if (p == 1 && q == 1 && r == 1 && s == 1)
        E = (5 * ((double) Z)) / 8;
    else if (p == 1 && q == 1 && r == 1 && s == 2)
        E = (4096 * sqrt(2)*((double) Z)) / 64827;
    else if (p == 1 && q == 1 && r == 1 && s == 3)
        E = (1269 * sqrt(3)*((double) Z)) / 50000;
    else if (p == 1 && q == 1 && r == 1 && s == 4)
        E = (416415744 * ((double) Z)) / 15083778125;
    else if (p == 1 && q == 1 && r == 2 && s == 1)
        E = (4096 * sqrt(2)*((double) Z)) / 64827;
    else if (p == 1 && q == 1 && r == 2 && s == 2)
        E = (16 * ((double) Z)) / 729;
    else if (p == 1 && q == 1 && r == 2 && s == 3)
        E = (110592 * sqrt(6)*((double) Z)) / 24137569;
    else if (p == 1 && q == 1 && r == 2 && s == 4)
        E = (98304 * sqrt(2)*((double) Z)) / 19487171;
    else if (p == 1 && q == 1 && r == 3 && s == 1)
        E = (1269 * sqrt(3)*((double) Z)) / 50000;
    else if (p == 1 && q == 1 && r == 3 && s == 2)
        E = (110592 * sqrt(6)*((double) Z)) / 24137569;
    else if (p == 1 && q == 1 && r == 3 && s == 3)
        E = (189 * ((double) Z)) / 32768;
    else if (p == 1 && q == 1 && r == 3 && s == 4)
        E = (1808400384 * sqrt(3)*((double) Z)) / 852891037441;
    else if (p == 1 && q == 1 && r == 4 && s == 1)
        E = (416415744 * ((double) Z)) / 15083778125;
    else if (p == 1 && q == 1 && r == 4 && s == 2)
        E = (98304 * sqrt(2)*((double) Z)) / 19487171;
    else if (p == 1 && q == 1 && r == 4 && s == 3)
        E = (1808400384 * sqrt(3)*((double) Z)) / 852891037441;
    else if (p == 1 && q == 1 && r == 4 && s == 4)
        E = (22848 * ((double) Z)) / 9765625;
    else if (p == 1 && q == 2 && r == 1 && s == 1)
        E = (4096 * sqrt(2)*((double) Z)) / 64827;
    else if (p == 1 && q == 2 && r == 1 && s == 2)
        E = (17 * ((double) Z)) / 81;
    else if (p == 1 && q == 2 && r == 1 && s == 3)
        E = (1555918848 * sqrt(6)*((double) Z)) / 75429903125;
    else if (p == 1 && q == 2 && r == 1 && s == 4)
        E = (288456704 * sqrt(2)*((double) Z)) / 14206147659;
    else if (p == 1 && q == 2 && r == 2 && s == 1)
        E = (16 * ((double) Z)) / 729;
    else if (p == 1 && q == 2 && r == 2 && s == 2)
        E = (512 * sqrt(2)*((double) Z)) / 84375;
    else if (p == 1 && q == 2 && r == 2 && s == 3)
        E = (2160 * sqrt(3)*((double) Z)) / 823543;
    else if (p == 1 && q == 2 && r == 2 && s == 4)
        E = (376832 * ((double) Z)) / 129140163;
    else if (p == 1 && q == 2 && r == 3 && s == 1)
        E = (110592 * sqrt(6)*((double) Z)) / 24137569;
    else if (p == 1 && q == 2 && r == 3 && s == 2)
        E = (29943 * sqrt(3)*((double) Z)) / 13176688;
    else if (p == 1 && q == 2 && r == 3 && s == 3)
        E = (1216512 * sqrt(2)*((double) Z)) / 815730721;
    else if (p == 1 && q == 2 && r == 3 && s == 4)
        E = (423788544 * sqrt(6)*((double) Z)) / 762939453125;
    else if (p == 1 && q == 2 && r == 4 && s == 1)
        E = (98304 * sqrt(2)*((double) Z)) / 19487171;
    else if (p == 1 && q == 2 && r == 4 && s == 2)
        E = (108838912 * ((double) Z)) / 44840334375;
    else if (p == 1 && q == 2 && r == 4 && s == 3)
        E = (10148806656 * sqrt(6)*((double) Z)) / 19073486328125;
    else if (p == 1 && q == 2 && r == 4 && s == 4)
        E = (39 * ((double) Z)) / (32768 * sqrt(2));
    else if (p == 1 && q == 3 && r == 1 && s == 1)
        E = (1269 * sqrt(3)*((double) Z)) / 50000;
    else if (p == 1 && q == 3 && r == 1 && s == 2)
        E = (1555918848 * sqrt(6)*((double) Z)) / 75429903125;
    else if (p == 1 && q == 3 && r == 1 && s == 3)
        E = (815 * ((double) Z)) / 8192;
    else if (p == 1 && q == 3 && r == 1 && s == 4)
        E = (11694850770862080 * sqrt(3)*((double) Z)) / 702392443647273463;
    else if (p == 1 && q == 3 && r == 2 && s == 1)
        E = (110592 * sqrt(6)*((double) Z)) / 24137569;
    else if (p == 1 && q == 3 && r == 2 && s == 2)
        E = (2160 * sqrt(3)*((double) Z)) / 823543;
    else if (p == 1 && q == 3 && r == 2 && s == 3)
        E = (37826560 * sqrt(2)*((double) Z)) / 22024729467;
    else if (p == 1 && q == 3 && r == 2 && s == 4)
        E = (487489536 * sqrt(6)*((double) Z)) / 762939453125;
    else if (p == 1 && q == 3 && r == 3 && s == 1)
        E = (189 * ((double) Z)) / 32768;
    else if (p == 1 && q == 3 && r == 3 && s == 2)
        E = (1216512 * sqrt(2)*((double) Z)) / 815730721;
    else if (p == 1 && q == 3 && r == 3 && s == 3)
        E = (617 * ((double) Z)) / (314928 * sqrt(3));
    else if (p == 1 && q == 3 && r == 3 && s == 4)
        E = (30254432256 * ((double) Z)) / 41426511213649;
    else if (p == 1 && q == 3 && r == 4 && s == 1)
        E = (1808400384 * sqrt(3)*((double) Z)) / 852891037441;
    else if (p == 1 && q == 3 && r == 4 && s == 2)
        E = (10148806656 * sqrt(6)*((double) Z)) / 19073486328125;
    else if (p == 1 && q == 3 && r == 4 && s == 3)
        E = (90581886173184 * ((double) Z)) / 129457847542653125;
    else if (p == 1 && q == 3 && r == 4 && s == 4)
        E = (74450880 * sqrt(3)*((double) Z)) / 285311670611;
    else if (p == 1 && q == 4 && r == 1 && s == 1)
        E = (416415744 * ((double) Z)) / 15083778125;
    else if (p == 1 && q == 4 && r == 1 && s == 2)
        E = (288456704 * sqrt(2)*((double) Z)) / 14206147659;
    else if (p == 1 && q == 4 && r == 1 && s == 3)
        E = (11694850770862080 * sqrt(3)*((double) Z)) / 702392443647273463;
    else if (p == 1 && q == 4 && r == 1 && s == 4)
        E = (22513 * ((double) Z)) / 390625;
    else if (p == 1 && q == 4 && r == 2 && s == 1)
        E = (98304 * sqrt(2)*((double) Z)) / 19487171;
    else if (p == 1 && q == 4 && r == 2 && s == 2)
        E = (376832 * ((double) Z)) / 129140163;
    else if (p == 1 && q == 4 && r == 2 && s == 3)
        E = (487489536 * sqrt(6)*((double) Z)) / 762939453125;
    else if (p == 1 && q == 4 && r == 2 && s == 4)
        E = (5053 * ((double) Z)) / (3538944 * sqrt(2));
    else if (p == 1 && q == 4 && r == 3 && s == 1)
        E = (1808400384 * sqrt(3)*((double) Z)) / 852891037441;
    else if (p == 1 && q == 4 && r == 3 && s == 2)
        E = (423788544 * sqrt(6)*((double) Z)) / 762939453125;
    else if (p == 1 && q == 4 && r == 3 && s == 3)
        E = (30254432256 * ((double) Z)) / 41426511213649;
    else if (p == 1 && q == 4 && r == 3 && s == 4)
        E = (1243165779 * sqrt(3)*((double) Z)) / 4564986729776;
    else if (p == 1 && q == 4 && r == 4 && s == 1)
        E = (22848 * ((double) Z)) / 9765625;
    else if (p == 1 && q == 4 && r == 4 && s == 2)
        E = (39 * ((double) Z)) / (32768 * sqrt(2));
    else if (p == 1 && q == 4 && r == 4 && s == 3)
        E = (74450880 * sqrt(3)*((double) Z)) / 285311670611;
    else if (p == 1 && q == 4 && r == 4 && s == 4)
        E = (1804351488 * ((double) Z)) / 6179146071875;
    else if (p == 2 && q == 1 && r == 1 && s == 1)
        E = (4096 * sqrt(2)*((double) Z)) / 64827;
    else if (p == 2 && q == 1 && r == 1 && s == 2)
        E = (16 * ((double) Z)) / 729;
    else if (p == 2 && q == 1 && r == 1 && s == 3)
        E = (110592 * sqrt(6)*((double) Z)) / 24137569;
    else if (p == 2 && q == 1 && r == 1 && s == 4)
        E = (98304 * sqrt(2)*((double) Z)) / 19487171;
    else if (p == 2 && q == 1 && r == 2 && s == 1)
        E = (17 * ((double) Z)) / 81;
    else if (p == 2 && q == 1 && r == 2 && s == 2)
        E = (512 * sqrt(2)*((double) Z)) / 84375;
    else if (p == 2 && q == 1 && r == 2 && s == 3)
        E = (29943 * sqrt(3)*((double) Z)) / 13176688;
    else if (p == 2 && q == 1 && r == 2 && s == 4)
        E = (108838912 * ((double) Z)) / 44840334375;
    else if (p == 2 && q == 1 && r == 3 && s == 1)
        E = (1555918848 * sqrt(6)*((double) Z)) / 75429903125;
    else if (p == 2 && q == 1 && r == 3 && s == 2)
        E = (2160 * sqrt(3)*((double) Z)) / 823543;
    else if (p == 2 && q == 1 && r == 3 && s == 3)
        E = (1216512 * sqrt(2)*((double) Z)) / 815730721;
    else if (p == 2 && q == 1 && r == 3 && s == 4)
        E = (10148806656 * sqrt(6)*((double) Z)) / 19073486328125;
    else if (p == 2 && q == 1 && r == 4 && s == 1)
        E = (288456704 * sqrt(2)*((double) Z)) / 14206147659;
    else if (p == 2 && q == 1 && r == 4 && s == 2)
        E = (376832 * ((double) Z)) / 129140163;
    else if (p == 2 && q == 1 && r == 4 && s == 3)
        E = (423788544 * sqrt(6)*((double) Z)) / 762939453125;
    else if (p == 2 && q == 1 && r == 4 && s == 4)
        E = (39 * ((double) Z)) / (32768 * sqrt(2));
    else if (p == 2 && q == 2 && r == 1 && s == 1)
        E = (16 * ((double) Z)) / 729;
    else if (p == 2 && q == 2 && r == 1 && s == 2)
        E = (512 * sqrt(2)*((double) Z)) / 84375;
    else if (p == 2 && q == 2 && r == 1 && s == 3)
        E = (2160 * sqrt(3)*((double) Z)) / 823543;
    else if (p == 2 && q == 2 && r == 1 && s == 4)
        E = (376832 * ((double) Z)) / 129140163;
    else if (p == 2 && q == 2 && r == 2 && s == 1)
        E = (512 * sqrt(2)*((double) Z)) / 84375;
    else if (p == 2 && q == 2 && r == 2 && s == 2)
        E = (77 * ((double) Z)) / 512;
    else if (p == 2 && q == 2 && r == 2 && s == 3)
        E = (5870679552.0 * sqrt(6)*((double) Z)) / 669871503125.0;
    else if (p == 2 && q == 2 && r == 2 && s == 4)
        E = (31363072.0 * sqrt(2)*((double) Z)) / 4202539929.0;
    else if (p == 2 && q == 2 && r == 3 && s == 1)
        E = (2160 * sqrt(3)*((double) Z)) / 823543;
    else if (p == 2 && q == 2 && r == 3 && s == 2)
        E = (5870679552.0 * sqrt(6)*((double) Z)) / 669871503125.0;
    else if (p == 2 && q == 2 && r == 3 && s == 3)
        E = (73008 * ((double) Z)) / 9765625;
    else if (p == 2 && q == 2 && r == 3 && s == 4)
        E = (14739259392.0 * sqrt(3)*((double) Z)) / 6131066257801.0;
    else if (p == 2 && q == 2 && r == 4 && s == 1)
        E = (376832 * ((double) Z)) / 129140163;
    else if (p == 2 && q == 2 && r == 4 && s == 2)
        E = (31363072.0 * sqrt(2)*((double) Z)) / 4202539929.0;
    else if (p == 2 && q == 2 && r == 4 && s == 3)
        E = (14739259392.0 * sqrt(3)*((double) Z)) / 6131066257801.0;
    else if (p == 2 && q == 2 && r == 4 && s == 4)
        E = (424 * ((double) Z)) / 177147;
    else if (p == 2 && q == 3 && r == 1 && s == 1)
        E = (110592 * sqrt(6)*((double) Z)) / 24137569;
    else if (p == 2 && q == 3 && r == 1 && s == 2)
        E = (2160 * sqrt(3)*((double) Z)) / 823543;
    else if (p == 2 && q == 3 && r == 1 && s == 3)
        E = (37826560 * sqrt(2)*((double) Z)) / 22024729467;
    else if (p == 2 && q == 3 && r == 1 && s == 4)
        E = (487489536 * sqrt(6)*((double) Z)) / 762939453125;
    else if (p == 2 && q == 3 && r == 2 && s == 1)
        E = (29943 * sqrt(3)*((double) Z)) / 13176688;
    else if (p == 2 && q == 3 && r == 2 && s == 2)
        E = (5870679552 * sqrt(6)*((double) Z)) / 669871503125;
    else if (p == 2 && q == 3 && r == 2 && s == 3)
        E = (32857 * ((double) Z)) / 390625;
    else if (p == 2 && q == 3 && r == 2 && s == 4)
        E = (55508689880137728 * sqrt(3)*((double) Z)) / 5049196699148208943;
    else if (p == 2 && q == 3 && r == 3 && s == 1)
        E = (1216512 * sqrt(2)*((double) Z)) / 815730721;
    else if (p == 2 && q == 3 && r == 3 && s == 2)
        E = (73008 * ((double) Z)) / 9765625;
    else if (p == 2 && q == 3 && r == 3 && s == 3)
        E = (6890942464 * sqrt(2 / 3)*((double) Z)) / 1210689028125;
    else if (p == 2 && q == 3 && r == 3 && s == 4)
        E = (69158928384 * sqrt(2)*((double) Z)) / 34271896307633;
    else if (p == 2 && q == 3 && r == 4 && s == 1)
        E = (423788544 * sqrt(6)*((double) Z)) / 762939453125;
    else if (p == 2 && q == 3 && r == 4 && s == 2)
        E = (14739259392 * sqrt(3)*((double) Z)) / 6131066257801;
    else if (p == 2 && q == 3 && r == 4 && s == 3)
        E = (36645380390912 * sqrt(2)*((double) Z)) / 24984212408264457;
    else if (p == 2 && q == 3 && r == 4 && s == 4)
        E = (145503 * sqrt(3 / 2)*((double) Z)) / 134217728;
    else if (p == 2 && q == 4 && r == 1 && s == 1)
        E = (98304 * sqrt(2)*((double) Z)) / 19487171;
    else if (p == 2 && q == 4 && r == 1 && s == 2)
        E = (376832 * ((double) Z)) / 129140163;
    else if (p == 2 && q == 4 && r == 1 && s == 3)
        E = (487489536 * sqrt(6)*((double) Z)) / 762939453125;
    else if (p == 2 && q == 4 && r == 1 && s == 4)
        E = (5053 * ((double) Z)) / (3538944 * sqrt(2));
    else if (p == 2 && q == 4 && r == 2 && s == 1)
        E = (108838912 * ((double) Z)) / 44840334375;
    else if (p == 2 && q == 4 && r == 2 && s == 2)
        E = (31363072.0 * sqrt(2)*((double) Z)) / 4202539929.0;
    else if (p == 2 && q == 4 && r == 2 && s == 3)
        E = (55508689880137728.0 * sqrt(3)*((double) Z)) / 5049196699148208943.0;
    else if (p == 2 && q == 4 && r == 2 && s == 4)
        E = (4043 * ((double) Z)) / 78732;
    else if (p == 2 && q == 4 && r == 3 && s == 1)
        E = (10148806656.0 * sqrt(6)*((double) Z)) / 19073486328125.0;
    else if (p == 2 && q == 4 && r == 3 && s == 2)
        E = (14739259392.0 * sqrt(3)*((double) Z)) / 6131066257801.0;
    else if (p == 2 && q == 4 && r == 3 && s == 3)
        E = (69158928384.0 * sqrt(2)*((double) Z)) / 34271896307633.0;
    else if (p == 2 && q == 4 && r == 3 && s == 4)
        E = (2496169683.0 * sqrt(3 / 2)*((double) Z)) / 1677721600000.0;
    else if (p == 2 && q == 4 && r == 4 && s == 1)
        E = (39 * ((double) Z)) / (32768 * sqrt(2));
    else if (p == 2 && q == 4 && r == 4 && s == 2)
        E = (424 * ((double) Z)) / 177147;
    else if (p == 2 && q == 4 && r == 4 && s == 3)
        E = (145503 * sqrt(3 / 2)*((double) Z)) / 134217728;
    else if (p == 2 && q == 4 && r == 4 && s == 4)
        E = (21252608 * sqrt(2)*((double) Z)) / 35595703125;
    else if (p == 3 && q == 1 && r == 1 && s == 1)
        E = (1269 * sqrt(3)*((double) Z)) / 50000;
    else if (p == 3 && q == 1 && r == 1 && s == 2)
        E = (110592 * sqrt(6)*((double) Z)) / 24137569;
    else if (p == 3 && q == 1 && r == 1 && s == 3)
        E = (189 * ((double) Z)) / 32768;
    else if (p == 3 && q == 1 && r == 1 && s == 4)
        E = (1808400384 * sqrt(3)*((double) Z)) / 852891037441;
    else if (p == 3 && q == 1 && r == 2 && s == 1)
        E = (1555918848 * sqrt(6)*((double) Z)) / 75429903125;
    else if (p == 3 && q == 1 && r == 2 && s == 2)
        E = (2160 * sqrt(3)*((double) Z)) / 823543;
    else if (p == 3 && q == 1 && r == 2 && s == 3)
        E = (1216512 * sqrt(2)*((double) Z)) / 815730721;
    else if (p == 3 && q == 1 && r == 2 && s == 4)
        E = (10148806656 * sqrt(6)*((double) Z)) / 19073486328125;
    else if (p == 3 && q == 1 && r == 3 && s == 1)
        E = (815 * ((double) Z)) / 8192;
    else if (p == 3 && q == 1 && r == 3 && s == 2)
        E = (37826560 * sqrt(2)*((double) Z)) / 22024729467;
    else if (p == 3 && q == 1 && r == 3 && s == 3)
        E = (617 * ((double) Z)) / (314928 * sqrt(3));
    else if (p == 3 && q == 1 && r == 3 && s == 4)
        E = (90581886173184 * ((double) Z)) / 129457847542653125;
    else if (p == 3 && q == 1 && r == 4 && s == 1)
        E = (11694850770862080 * sqrt(3)*((double) Z)) / 702392443647273463;
    else if (p == 3 && q == 1 && r == 4 && s == 2)
        E = (487489536 * sqrt(6)*((double) Z)) / 762939453125;
    else if (p == 3 && q == 1 && r == 4 && s == 3)
        E = (30254432256 * ((double) Z)) / 41426511213649;
    else if (p == 3 && q == 1 && r == 4 && s == 4)
        E = (74450880 * sqrt(3)*((double) Z)) / 285311670611;
    else if (p == 3 && q == 2 && r == 1 && s == 1)
        E = (110592 * sqrt(6)*((double) Z)) / 24137569;
    else if (p == 3 && q == 2 && r == 1 && s == 2)
        E = (29943 * sqrt(3)*((double) Z)) / 13176688;
    else if (p == 3 && q == 2 && r == 1 && s == 3)
        E = (1216512 * sqrt(2)*((double) Z)) / 815730721;
    else if (p == 3 && q == 2 && r == 1 && s == 4)
        E = (423788544 * sqrt(6)*((double) Z)) / 762939453125;
    else if (p == 3 && q == 2 && r == 2 && s == 1)
        E = (2160 * sqrt(3)*((double) Z)) / 823543;
    else if (p == 3 && q == 2 && r == 2 && s == 2)
        E = (5870679552 * sqrt(6)*((double) Z)) / 669871503125;
    else if (p == 3 && q == 2 && r == 2 && s == 3)
        E = (73008 * ((double) Z)) / 9765625;
    else if (p == 3 && q == 2 && r == 2 && s == 4)
        E = (14739259392 * sqrt(3)*((double) Z)) / 6131066257801;
    else if (p == 3 && q == 2 && r == 3 && s == 1)
        E = (37826560 * sqrt(2)*((double) Z)) / 22024729467;
    else if (p == 3 && q == 2 && r == 3 && s == 2)
        E = (32857 * ((double) Z)) / 390625;
    else if (p == 3 && q == 2 && r == 3 && s == 3)
        E = (6890942464 * sqrt(2 / 3)*((double) Z)) / 1210689028125;
    else if (p == 3 && q == 2 && r == 3 && s == 4)
        E = (36645380390912 * sqrt(2)*((double) Z)) / 24984212408264457;
    else if (p == 3 && q == 2 && r == 4 && s == 1)
        E = (487489536 * sqrt(6)*((double) Z)) / 762939453125;
    else if (p == 3 && q == 2 && r == 4 && s == 2)
        E = (55508689880137728 * sqrt(3)*((double) Z)) / 5049196699148208943;
    else if (p == 3 && q == 2 && r == 4 && s == 3)
        E = (69158928384 * sqrt(2)*((double) Z)) / 34271896307633;
    else if (p == 3 && q == 2 && r == 4 && s == 4)
        E = (145503 * sqrt(3 / 2)*((double) Z)) / 134217728;
    else if (p == 3 && q == 3 && r == 1 && s == 1)
        E = (189 * ((double) Z)) / 32768;
    else if (p == 3 && q == 3 && r == 1 && s == 2)
        E = (1216512 * sqrt(2)*((double) Z)) / 815730721;
    else if (p == 3 && q == 3 && r == 1 && s == 3)
        E = (617 * ((double) Z)) / (314928 * sqrt(3));
    else if (p == 3 && q == 3 && r == 1 && s == 4)
        E = (30254432256 * ((double) Z)) / 41426511213649;
    else if (p == 3 && q == 3 && r == 2 && s == 1)
        E = (1216512 * sqrt(2)*((double) Z)) / 815730721;
    else if (p == 3 && q == 3 && r == 2 && s == 2)
        E = (73008 * ((double) Z)) / 9765625;
    else if (p == 3 && q == 3 && r == 2 && s == 3)
        E = (6890942464 * sqrt(2 / 3)*((double) Z)) / 1210689028125;
    else if (p == 3 && q == 3 && r == 2 && s == 4)
        E = (69158928384 * sqrt(2)*((double) Z)) / 34271896307633;
    else if (p == 3 && q == 3 && r == 3 && s == 1)
        E = (617 * ((double) Z)) / (314928 * sqrt(3));
    else if (p == 3 && q == 3 && r == 3 && s == 2)
        E = (6890942464 * sqrt(2 / 3)*((double) Z)) / 1210689028125;
    else if (p == 3 && q == 3 && r == 3 && s == 3)
        E = (17 * ((double) Z)) / 256;
    else if (p == 3 && q == 3 && r == 3 && s == 4)
        E = (2486755845603328 * ((double) Z)) / (158298797548828125 * sqrt(3));
    else if (p == 3 && q == 3 && r == 4 && s == 1)
        E = (30254432256 * ((double) Z)) / 41426511213649;
    else if (p == 3 && q == 3 && r == 4 && s == 2)
        E = (69158928384 * sqrt(2)*((double) Z)) / 34271896307633;
    else if (p == 3 && q == 3 && r == 4 && s == 3)
        E = (2486755845603328 * ((double) Z)) / (158298797548828125 * sqrt(3));
    else if (p == 3 && q == 3 && r == 4 && s == 4)
        E = (2560158144.0 * ((double) Z)) / 678223072849.0;
    else if (p == 3 && q == 4 && r == 1 && s == 1)
        E = (1808400384.0 * sqrt(3)*((double) Z)) / 852891037441.0;
    else if (p == 3 && q == 4 && r == 1 && s == 2)
        E = (423788544.0 * sqrt(6)*((double) Z)) / 762939453125.0;
    else if (p == 3 && q == 4 && r == 1 && s == 3)
        E = (30254432256.0 * ((double) Z)) / 41426511213649.0;
    else if (p == 3 && q == 4 && r == 1 && s == 4)
        E = (1243165779.0 * sqrt(3)*((double) Z)) / 4564986729776.0;
    else if (p == 3 && q == 4 && r == 2 && s == 1)
        E = (10148806656.0 * sqrt(6)*((double) Z)) / 19073486328125.0;
    else if (p == 3 && q == 4 && r == 2 && s == 2)
        E = (14739259392.0 * sqrt(3)*((double) Z)) / 6131066257801.0;
    else if (p == 3 && q == 4 && r == 2 && s == 3)
        E = (69158928384 * sqrt(2)*((double) Z)) / 34271896307633;
    else if (p == 3 && q == 4 && r == 2 && s == 4)
        E = (2496169683.0 * sqrt(3 / 2)*((double) Z)) / 1677721600000.0;
    else if (p == 3 && q == 4 && r == 3 && s == 1)
        E = (90581886173184.0 * ((double) Z)) / 129457847542653125.0;
    else if (p == 3 && q == 4 && r == 3 && s == 2)
        E = (36645380390912.0 * sqrt(2)*((double) Z)) / 24984212408264457.0;
    else if (p == 3 && q == 4 && r == 3 && s == 3)
        E = (2486755845603328.0 * ((double) Z)) / (158298797548828125.0 * sqrt(3));
    else if (p == 3 && q == 4 && r == 3 && s == 4)
        E = (621550729.0 * ((double) Z)) / 13841287201.0;
    else if (p == 3 && q == 4 && r == 4 && s == 1)
        E = (74450880.0 * sqrt(3)*((double) Z)) / 285311670611.0;
    else if (p == 3 && q == 4 && r == 4 && s == 2)
        E = (145503 * sqrt(3 / 2)*((double) Z)) / 134217728.0;
    else if (p == 3 && q == 4 && r == 4 && s == 3)
        E = (2560158144.0 * ((double) Z)) / 678223072849.0;
    else if (p == 3 && q == 4 && r == 4 && s == 4)
        E = (413631006610176000.0 * sqrt(3)*((double) Z)) / 249430673908303812379.0;
    else if (p == 4 && q == 1 && r == 1 && s == 1)
        E = (416415744 * ((double) Z)) / 15083778125;
    else if (p == 4 && q == 1 && r == 1 && s == 2)
        E = (98304 * sqrt(2)*((double) Z)) / 19487171;
    else if (p == 4 && q == 1 && r == 1 && s == 3)
        E = (1808400384 * sqrt(3)*((double) Z)) / 852891037441;
    else if (p == 4 && q == 1 && r == 1 && s == 4)
        E = (22848 * ((double) Z)) / 9765625;
    else if (p == 4 && q == 1 && r == 2 && s == 1)
        E = (288456704 * sqrt(2)*((double) Z)) / 14206147659;
    else if (p == 4 && q == 1 && r == 2 && s == 2)
        E = (376832 * ((double) Z)) / 129140163;
    else if (p == 4 && q == 1 && r == 2 && s == 3)
        E = (423788544 * sqrt(6)*((double) Z)) / 762939453125;
    else if (p == 4 && q == 1 && r == 2 && s == 4)
        E = (39 * ((double) Z)) / (32768 * sqrt(2));
    else if (p == 4 && q == 1 && r == 3 && s == 1)
        E = (11694850770862080 * sqrt(3)*((double) Z)) / 702392443647273463;
    else if (p == 4 && q == 1 && r == 3 && s == 2)
        E = (487489536 * sqrt(6)*((double) Z)) / 762939453125;
    else if (p == 4 && q == 1 && r == 3 && s == 3)
        E = (30254432256 * ((double) Z)) / 41426511213649;
    else if (p == 4 && q == 1 && r == 3 && s == 4)
        E = (74450880 * sqrt(3)*((double) Z)) / 285311670611;
    else if (p == 4 && q == 1 && r == 4 && s == 1)
        E = (22513 * ((double) Z)) / 390625;
    else if (p == 4 && q == 1 && r == 4 && s == 2)
        E = (5053 * ((double) Z)) / (3538944 * sqrt(2));
    else if (p == 4 && q == 1 && r == 4 && s == 3)
        E = (1243165779 * sqrt(3)*((double) Z)) / 4564986729776;
    else if (p == 4 && q == 1 && r == 4 && s == 4)
        E = (1804351488 * ((double) Z)) / 6179146071875;
    else if (p == 4 && q == 2 && r == 1 && s == 1)
        E = (98304 * sqrt(2)*((double) Z)) / 19487171;
    else if (p == 4 && q == 2 && r == 1 && s == 2)
        E = (108838912 * ((double) Z)) / 44840334375;
    else if (p == 4 && q == 2 && r == 1 && s == 3)
        E = (10148806656.0 * sqrt(6)*((double) Z)) / 19073486328125.0;
    else if (p == 4 && q == 2 && r == 1 && s == 4)
        E = (39 * ((double) Z)) / (32768 * sqrt(2));
    else if (p == 4 && q == 2 && r == 2 && s == 1)
        E = (376832 * ((double) Z)) / 129140163;
    else if (p == 4 && q == 2 && r == 2 && s == 2)
        E = (31363072.0 * sqrt(2)*((double) Z)) / 4202539929.0;
    else if (p == 4 && q == 2 && r == 2 && s == 3)
        E = (14739259392 * sqrt(3)*((double) Z)) / 6131066257801;
    else if (p == 4 && q == 2 && r == 2 && s == 4)
        E = (424 * ((double) Z)) / 177147;
    else if (p == 4 && q == 2 && r == 3 && s == 1)
        E = (487489536 * sqrt(6)*((double) Z)) / 762939453125;
    else if (p == 4 && q == 2 && r == 3 && s == 2)
        E = (55508689880137728 * sqrt(3)*((double) Z)) / 5049196699148208943;
    else if (p == 4 && q == 2 && r == 3 && s == 3)
        E = (69158928384 * sqrt(2)*((double) Z)) / 34271896307633;
    else if (p == 4 && q == 2 && r == 3 && s == 4)
        E = (145503 * sqrt(3 / 2)*((double) Z)) / 134217728;
    else if (p == 4 && q == 2 && r == 4 && s == 1)
        E = (5053 * ((double) Z)) / (3538944 * sqrt(2));
    else if (p == 4 && q == 2 && r == 4 && s == 2)
        E = (4043 * ((double) Z)) / 78732;
    else if (p == 4 && q == 2 && r == 4 && s == 3)
        E = (2496169683.0 * sqrt(3 / 2)*((double) Z)) / 1677721600000.0;
    else if (p == 4 && q == 2 && r == 4 && s == 4)
        E = (21252608.0 * sqrt(2)*((double) Z)) / 35595703125.0;
    else if (p == 4 && q == 3 && r == 1 && s == 1)
        E = (1808400384.0 * sqrt(3)*((double) Z)) / 852891037441.0;
    else if (p == 4 && q == 3 && r == 1 && s == 2)
        E = (10148806656.0 * sqrt(6)*((double) Z)) / 19073486328125.0;
    else if (p == 4 && q == 3 && r == 1 && s == 3)
        E = (90581886173184.0 * ((double) Z)) / 129457847542653125.0;
    else if (p == 4 && q == 3 && r == 1 && s == 4)
        E = (74450880.0 * sqrt(3)*((double) Z)) / 285311670611.0;
    else if (p == 4 && q == 3 && r == 2 && s == 1)
        E = (423788544 * sqrt(6)*((double) Z)) / 762939453125;
    else if (p == 4 && q == 3 && r == 2 && s == 2)
        E = (14739259392 * sqrt(3)*((double) Z)) / 6131066257801;
    else if (p == 4 && q == 3 && r == 2 && s == 3)
        E = (36645380390912 * sqrt(2)*((double) Z)) / 24984212408264457;
    else if (p == 4 && q == 3 && r == 2 && s == 4)
        E = (145503 * sqrt(3 / 2)*((double) Z)) / 134217728;
    else if (p == 4 && q == 3 && r == 3 && s == 1)
        E = (30254432256 * ((double) Z)) / 41426511213649;
    else if (p == 4 && q == 3 && r == 3 && s == 2)
        E = (69158928384 * sqrt(2)*((double) Z)) / 34271896307633;
    else if (p == 4 && q == 3 && r == 3 && s == 3)
        E = (2486755845603328.0 * ((double) Z)) / (158298797548828125.0 * sqrt(3));
    else if (p == 4 && q == 3 && r == 3 && s == 4)
        E = (2560158144.0 * ((double) Z)) / 678223072849.0;
    else if (p == 4 && q == 3 && r == 4 && s == 1)
        E = (1243165779.0 * sqrt(3)*((double) Z)) / 4564986729776.0;
    else if (p == 4 && q == 3 && r == 4 && s == 2)
        E = (2496169683.0 * sqrt(3 / 2)*((double) Z)) / 1677721600000.0;
    else if (p == 4 && q == 3 && r == 4 && s == 3)
        E = (621550729 * ((double) Z)) / 13841287201;
    else if (p == 4 && q == 3 && r == 4 && s == 4)
        E = (413631006610176000.0 * sqrt(3)*((double) Z)) / 249430673908303812379.0;
    else if (p == 4 && q == 4 && r == 1 && s == 1)
        E = (22848 * ((double) Z)) / 9765625;
    else if (p == 4 && q == 4 && r == 1 && s == 2)
        E = (39 * ((double) Z)) / (32768 * sqrt(2));
    else if (p == 4 && q == 4 && r == 1 && s == 3)
        E = (74450880 * sqrt(3)*((double) Z)) / 285311670611;
    else if (p == 4 && q == 4 && r == 1 && s == 4)
        E = (1804351488 * ((double) Z)) / 6179146071875;
    else if (p == 4 && q == 4 && r == 2 && s == 1)
        E = (39 * ((double) Z)) / (32768 * sqrt(2));
    else if (p == 4 && q == 4 && r == 2 && s == 2)
        E = (424 * ((double) Z)) / 177147;
    else if (p == 4 && q == 4 && r == 2 && s == 3)
        E = (145503 * sqrt(3 / 2)*((double) Z)) / 134217728;
    else if (p == 4 && q == 4 && r == 2 && s == 4)
        E = (21252608 * sqrt(2)*((double) Z)) / 35595703125;
    else if (p == 4 && q == 4 && r == 3 && s == 1)
        E = (74450880 * sqrt(3)*((double) Z)) / 285311670611;
    else if (p == 4 && q == 4 && r == 3 && s == 2)
        E = (145503 * sqrt(3 / 2)*((double) Z)) / 134217728;
    else if (p == 4 && q == 4 && r == 3 && s == 3)
        E = (2560158144.0 * ((double) Z)) / 678223072849.0;
    else if (p == 4 && q == 4 && r == 3 && s == 4)
        E = (413631006610176000.0 * sqrt(3)*((double) Z)) / 249430673908303812379.0;
    else if (p == 4 && q == 4 && r == 4 && s == 1)
        E = (1804351488 * ((double) Z)) / 6179146071875;
    else if (p == 4 && q == 4 && r == 4 && s == 2)
        E = (21252608 * sqrt(2)*((double) Z)) / 35595703125;
    else if (p == 4 && q == 4 && r == 4 && s == 3)
        E = (413631006610176000.0 * sqrt(3)*((double) Z)) / 249430673908303812379.0;
    else if (p == 4 && q == 4 && r == 4 && s == 4)
        E = (19541 * ((double) Z)) / 524288;
    else
        throw std::string("TROUBLE with matrix element.");


    return E;
}


