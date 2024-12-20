//This file is maintained by an external python script and should not be edited manually.
#ifndef BASISBANK_H
#define BASISBANK_H
#include <armadillo>
#include <basis.h>
#include <primitive.h>
using namespace std;
using namespace arma;
 
class basisbank{
public:
    basisbank(basis BS);
    basisbank();
    basis bs;
    string basistype;    void add_3_21G_h(vec3 corePos);
    void add_3_21G_he(vec3 corePos);
    void add_3_21G_li(vec3 corePos);
    void add_3_21G_be(vec3 corePos);
    void add_3_21G_b(vec3 corePos);
    void add_3_21G_c(vec3 corePos);
    void add_3_21G_n(vec3 corePos);
    void add_3_21G_o(vec3 corePos);
    void add_3_21G_f(vec3 corePos);
    void add_3_21G_ne(vec3 corePos);
    void add_3_21G_na(vec3 corePos);
    void add_3_21G_mg(vec3 corePos);
    void add_3_21G_al(vec3 corePos);
    void add_3_21G_si(vec3 corePos);
    void add_3_21G_p(vec3 corePos);
    void add_3_21G_s(vec3 corePos);
    void add_3_21G_cl(vec3 corePos);
    void add_3_21G_ar(vec3 corePos);
    void add_3_21G_k(vec3 corePos);
    void add_3_21G_ca(vec3 corePos);
    void add_3_21G_sc(vec3 corePos);
    void add_3_21G_ti(vec3 corePos);
    void add_3_21G_v(vec3 corePos);
    void add_3_21G_cr(vec3 corePos);
    void add_3_21G_mn(vec3 corePos);
    void add_3_21G_fe(vec3 corePos);
    void add_3_21G_co(vec3 corePos);
    void add_3_21G_ni(vec3 corePos);
    void add_3_21G_cu(vec3 corePos);
    void add_3_21G_zn(vec3 corePos);
    void add_3_21G_ga(vec3 corePos);
    void add_3_21G_ge(vec3 corePos);
    void add_3_21G_as(vec3 corePos);
    void add_3_21G_se(vec3 corePos);
    void add_3_21G_br(vec3 corePos);
    void add_6_311G2df2pd_h(vec3 corePos);
    void add_6_311G2df2pd_he(vec3 corePos);
    void add_6_311G2df2pd_li(vec3 corePos);
    void add_6_311G2df2pd_be(vec3 corePos);
    void add_6_311G2df2pd_b(vec3 corePos);
    void add_6_311G2df2pd_c(vec3 corePos);
    void add_6_311G2df2pd_n(vec3 corePos);
    void add_6_311G2df2pd_o(vec3 corePos);
    void add_6_311G2df2pd_f(vec3 corePos);
    void add_6_311G_h(vec3 corePos);
    void add_6_311G_he(vec3 corePos);
    void add_6_311G_li(vec3 corePos);
    void add_6_311G_be(vec3 corePos);
    void add_6_311G_b(vec3 corePos);
    void add_6_311G_c(vec3 corePos);
    void add_6_311G_n(vec3 corePos);
    void add_6_311G_o(vec3 corePos);
    void add_6_311G_f(vec3 corePos);
    void add_6_311G_ne(vec3 corePos);
    void add_6_311G_na(vec3 corePos);
    void add_6_311G_mg(vec3 corePos);
    void add_6_311G_al(vec3 corePos);
    void add_6_311G_si(vec3 corePos);
    void add_6_311G_p(vec3 corePos);
    void add_6_311G_s(vec3 corePos);
    void add_6_311G_cl(vec3 corePos);
    void add_STO_3G_h(vec3 corePos);
    void add_STO_3G_he(vec3 corePos);
    void add_STO_3G_li(vec3 corePos);
    void add_STO_3G_be(vec3 corePos);
    void add_STO_3G_b(vec3 corePos);
    void add_STO_3G_c(vec3 corePos);
    void add_STO_3G_n(vec3 corePos);
    void add_STO_3G_o(vec3 corePos);
    void add_STO_3G_f(vec3 corePos);
    void add_STO_3G_ne(vec3 corePos);
    void add_STO_3G_na(vec3 corePos);
    void add_STO_3G_mg(vec3 corePos);
    void add_STO_3G_al(vec3 corePos);
    void add_STO_3G_si(vec3 corePos);
    void add_STO_3G_p(vec3 corePos);
    void add_STO_3G_s(vec3 corePos);
    void add_STO_3G_cl(vec3 corePos);
    void add_STO_3G_ar(vec3 corePos);
    void add_STO_3G_k(vec3 corePos);
    void add_STO_3G_ca(vec3 corePos);
    void add_STO_3G_sc(vec3 corePos);
    void add_STO_3G_ti(vec3 corePos);
    void add_STO_3G_v(vec3 corePos);
    void add_STO_3G_cr(vec3 corePos);
    void add_STO_3G_mn(vec3 corePos);
    void add_STO_3G_fe(vec3 corePos);
    void add_STO_3G_co(vec3 corePos);
    void add_STO_3G_ni(vec3 corePos);
    void add_STO_3G_cu(vec3 corePos);
    void add_STO_3G_zn(vec3 corePos);
    void add_STO_3G_ga(vec3 corePos);
    void add_STO_3G_ge(vec3 corePos);
    void add_STO_3G_as(vec3 corePos);
    void add_STO_3G_se(vec3 corePos);
    void add_STO_3G_br(vec3 corePos);
    void add_STO_6G_h(vec3 corePos);
    void add_STO_6G_he(vec3 corePos);
    void add_STO_6G_li(vec3 corePos);
    void add_STO_6G_be(vec3 corePos);
    void add_STO_6G_b(vec3 corePos);
    void add_STO_6G_c(vec3 corePos);
    void add_STO_6G_n(vec3 corePos);
    void add_STO_6G_o(vec3 corePos);
    void add_STO_6G_f(vec3 corePos);
    void add_STO_6G_ne(vec3 corePos);
    void add_STO_6G_na(vec3 corePos);
    void add_STO_6G_mg(vec3 corePos);
    void add_STO_6G_al(vec3 corePos);
    void add_STO_6G_si(vec3 corePos);
    void add_STO_6G_p(vec3 corePos);
    void add_STO_6G_s(vec3 corePos);
    void add_STO_6G_cl(vec3 corePos);
    void add_STO_6G_ar(vec3 corePos);
    void add_STO_6G_k(vec3 corePos);
    void add_STO_6G_ca(vec3 corePos);
    void add_STO_6G_sc(vec3 corePos);
    void add_STO_6G_ti(vec3 corePos);
    void add_STO_6G_v(vec3 corePos);
    void add_STO_6G_cr(vec3 corePos);
    void add_STO_6G_mn(vec3 corePos);
    void add_STO_6G_fe(vec3 corePos);
    void add_STO_6G_co(vec3 corePos);
    void add_STO_6G_ni(vec3 corePos);
    void add_STO_6G_cu(vec3 corePos);
    void add_STO_6G_zn(vec3 corePos);
    void add_STO_6G_ga(vec3 corePos);
    void add_STO_6G_ge(vec3 corePos);
    void add_STO_6G_as(vec3 corePos);
    void add_STO_6G_se(vec3 corePos);
    void add_STO_6G_br(vec3 corePos);
};
#endif // BASISBANK_H