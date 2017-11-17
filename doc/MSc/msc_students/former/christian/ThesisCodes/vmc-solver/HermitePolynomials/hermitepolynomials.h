#ifndef HERMITEPOLYNOMIALS_H
#define HERMITEPOLYNOMIALS_H

#include "math.h"


class HermitePolynomials
{
public:
    HermitePolynomials();
    virtual double eval(double x) = 0;
protected:
    double m_alpha = 0;
    double m_omega = 0;
    double m_alphaomega = 0;
};



class HermitePolynomial_0 : public HermitePolynomials {
public:

    HermitePolynomial_0(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_0 : public HermitePolynomials {
public:

    dell_HermitePolynomial_0(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_0 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_0(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 0 ----------------------


class HermitePolynomial_1 : public HermitePolynomials {
public:

    HermitePolynomial_1(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_1 : public HermitePolynomials {
public:

    dell_HermitePolynomial_1(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_1 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_1(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 1 ----------------------


class HermitePolynomial_2 : public HermitePolynomials {
public:

    HermitePolynomial_2(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_2 : public HermitePolynomials {
public:

    dell_HermitePolynomial_2(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_2 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_2(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 2 ----------------------


class HermitePolynomial_3 : public HermitePolynomials {
public:

    HermitePolynomial_3(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_3 : public HermitePolynomials {
public:

    dell_HermitePolynomial_3(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_3 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_3(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 3 ----------------------


class HermitePolynomial_4 : public HermitePolynomials {
public:

    HermitePolynomial_4(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_4 : public HermitePolynomials {
public:

    dell_HermitePolynomial_4(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_4 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_4(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 4 ----------------------


class HermitePolynomial_5 : public HermitePolynomials {
public:

    HermitePolynomial_5(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_5 : public HermitePolynomials {
public:

    dell_HermitePolynomial_5(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_5 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_5(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 5 ----------------------


class HermitePolynomial_6 : public HermitePolynomials {
public:

    HermitePolynomial_6(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_6 : public HermitePolynomials {
public:

    dell_HermitePolynomial_6(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_6 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_6(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 6 ----------------------


class HermitePolynomial_7 : public HermitePolynomials {
public:

    HermitePolynomial_7(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_7 : public HermitePolynomials {
public:

    dell_HermitePolynomial_7(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_7 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_7(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 7 ----------------------


class HermitePolynomial_8 : public HermitePolynomials {
public:

    HermitePolynomial_8(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_8 : public HermitePolynomials {
public:

    dell_HermitePolynomial_8(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_8 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_8(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 8 ----------------------


class HermitePolynomial_9 : public HermitePolynomials {
public:

    HermitePolynomial_9(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_9 : public HermitePolynomials {
public:

    dell_HermitePolynomial_9(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_9 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_9(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 9 ----------------------


class HermitePolynomial_10 : public HermitePolynomials {
public:

    HermitePolynomial_10(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_10 : public HermitePolynomials {
public:

    dell_HermitePolynomial_10(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_10 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_10(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 10 ----------------------


class HermitePolynomial_11 : public HermitePolynomials {
public:

    HermitePolynomial_11(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_11 : public HermitePolynomials {
public:

    dell_HermitePolynomial_11(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_11 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_11(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 11 ----------------------


class HermitePolynomial_12 : public HermitePolynomials {
public:

    HermitePolynomial_12(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_12 : public HermitePolynomials {
public:

    dell_HermitePolynomial_12(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_12 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_12(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 12 ----------------------


class HermitePolynomial_13 : public HermitePolynomials {
public:

    HermitePolynomial_13(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_13 : public HermitePolynomials {
public:

    dell_HermitePolynomial_13(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_13 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_13(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 13 ----------------------


class HermitePolynomial_14 : public HermitePolynomials {
public:

    HermitePolynomial_14(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_14 : public HermitePolynomials {
public:

    dell_HermitePolynomial_14(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_14 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_14(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 14 ----------------------


class HermitePolynomial_15 : public HermitePolynomials {
public:

    HermitePolynomial_15(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_15 : public HermitePolynomials {
public:

    dell_HermitePolynomial_15(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_15 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_15(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 15 ----------------------


class HermitePolynomial_16 : public HermitePolynomials {
public:

    HermitePolynomial_16(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_16 : public HermitePolynomials {
public:

    dell_HermitePolynomial_16(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_16 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_16(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 16 ----------------------


class HermitePolynomial_17 : public HermitePolynomials {
public:

    HermitePolynomial_17(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_17 : public HermitePolynomials {
public:

    dell_HermitePolynomial_17(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_17 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_17(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 17 ----------------------


class HermitePolynomial_18 : public HermitePolynomials {
public:

    HermitePolynomial_18(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_18 : public HermitePolynomials {
public:

    dell_HermitePolynomial_18(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_18 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_18(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 18 ----------------------


class HermitePolynomial_19 : public HermitePolynomials {
public:

    HermitePolynomial_19(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_19 : public HermitePolynomials {
public:

    dell_HermitePolynomial_19(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_19 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_19(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 19 ----------------------


class HermitePolynomial_20 : public HermitePolynomials {
public:

    HermitePolynomial_20(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_20 : public HermitePolynomials {
public:

    dell_HermitePolynomial_20(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_20 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_20(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 20 ----------------------


class HermitePolynomial_21 : public HermitePolynomials {
public:

    HermitePolynomial_21(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_21 : public HermitePolynomials {
public:

    dell_HermitePolynomial_21(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_21 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_21(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 21 ----------------------


class HermitePolynomial_22 : public HermitePolynomials {
public:

    HermitePolynomial_22(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_22 : public HermitePolynomials {
public:

    dell_HermitePolynomial_22(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_22 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_22(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 22 ----------------------


class HermitePolynomial_23 : public HermitePolynomials {
public:

    HermitePolynomial_23(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_23 : public HermitePolynomials {
public:

    dell_HermitePolynomial_23(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_23 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_23(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 23 ----------------------


class HermitePolynomial_24 : public HermitePolynomials {
public:

    HermitePolynomial_24(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_24 : public HermitePolynomials {
public:

    dell_HermitePolynomial_24(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_24 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_24(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 24 ----------------------


class HermitePolynomial_25 : public HermitePolynomials {
public:

    HermitePolynomial_25(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_25 : public HermitePolynomials {
public:

    dell_HermitePolynomial_25(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_25 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_25(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 25 ----------------------


class HermitePolynomial_26 : public HermitePolynomials {
public:

    HermitePolynomial_26(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_26 : public HermitePolynomials {
public:

    dell_HermitePolynomial_26(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_26 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_26(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 26 ----------------------


class HermitePolynomial_27 : public HermitePolynomials {
public:

    HermitePolynomial_27(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_27 : public HermitePolynomials {
public:

    dell_HermitePolynomial_27(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_27 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_27(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 27 ----------------------


class HermitePolynomial_28 : public HermitePolynomials {
public:

    HermitePolynomial_28(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_28 : public HermitePolynomials {
public:

    dell_HermitePolynomial_28(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_28 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_28(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 28 ----------------------


class HermitePolynomial_29 : public HermitePolynomials {
public:

    HermitePolynomial_29(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_29 : public HermitePolynomials {
public:

    dell_HermitePolynomial_29(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_29 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_29(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 29 ----------------------


class HermitePolynomial_30 : public HermitePolynomials {
public:

    HermitePolynomial_30(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_30 : public HermitePolynomials {
public:

    dell_HermitePolynomial_30(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_30 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_30(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 30 ----------------------


class HermitePolynomial_31 : public HermitePolynomials {
public:

    HermitePolynomial_31(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_31 : public HermitePolynomials {
public:

    dell_HermitePolynomial_31(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_31 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_31(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 31 ----------------------


class HermitePolynomial_32 : public HermitePolynomials {
public:

    HermitePolynomial_32(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_32 : public HermitePolynomials {
public:

    dell_HermitePolynomial_32(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_32 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_32(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 32 ----------------------


class HermitePolynomial_33 : public HermitePolynomials {
public:

    HermitePolynomial_33(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_33 : public HermitePolynomials {
public:

    dell_HermitePolynomial_33(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_33 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_33(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 33 ----------------------


class HermitePolynomial_34 : public HermitePolynomials {
public:

    HermitePolynomial_34(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_34 : public HermitePolynomials {
public:

    dell_HermitePolynomial_34(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_34 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_34(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 34 ----------------------


class HermitePolynomial_35 : public HermitePolynomials {
public:

    HermitePolynomial_35(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_35 : public HermitePolynomials {
public:

    dell_HermitePolynomial_35(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_35 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_35(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 35 ----------------------


class HermitePolynomial_36 : public HermitePolynomials {
public:

    HermitePolynomial_36(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_36 : public HermitePolynomials {
public:

    dell_HermitePolynomial_36(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_36 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_36(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 36 ----------------------


class HermitePolynomial_37 : public HermitePolynomials {
public:

    HermitePolynomial_37(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_37 : public HermitePolynomials {
public:

    dell_HermitePolynomial_37(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_37 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_37(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 37 ----------------------


class HermitePolynomial_38 : public HermitePolynomials {
public:

    HermitePolynomial_38(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_38 : public HermitePolynomials {
public:

    dell_HermitePolynomial_38(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_38 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_38(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 38 ----------------------


class HermitePolynomial_39 : public HermitePolynomials {
public:

    HermitePolynomial_39(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_39 : public HermitePolynomials {
public:

    dell_HermitePolynomial_39(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_39 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_39(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 39 ----------------------


class HermitePolynomial_40 : public HermitePolynomials {
public:

    HermitePolynomial_40(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_40 : public HermitePolynomials {
public:

    dell_HermitePolynomial_40(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_40 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_40(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 40 ----------------------


class HermitePolynomial_41 : public HermitePolynomials {
public:

    HermitePolynomial_41(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_41 : public HermitePolynomials {
public:

    dell_HermitePolynomial_41(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_41 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_41(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 41 ----------------------


class HermitePolynomial_42 : public HermitePolynomials {
public:

    HermitePolynomial_42(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_42 : public HermitePolynomials {
public:

    dell_HermitePolynomial_42(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_42 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_42(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 42 ----------------------


class HermitePolynomial_43 : public HermitePolynomials {
public:

    HermitePolynomial_43(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_43 : public HermitePolynomials {
public:

    dell_HermitePolynomial_43(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_43 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_43(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 43 ----------------------


class HermitePolynomial_44 : public HermitePolynomials {
public:

    HermitePolynomial_44(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_44 : public HermitePolynomials {
public:

    dell_HermitePolynomial_44(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_44 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_44(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 44 ----------------------


class HermitePolynomial_45 : public HermitePolynomials {
public:

    HermitePolynomial_45(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_45 : public HermitePolynomials {
public:

    dell_HermitePolynomial_45(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_45 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_45(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 45 ----------------------


class HermitePolynomial_46 : public HermitePolynomials {
public:

    HermitePolynomial_46(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_46 : public HermitePolynomials {
public:

    dell_HermitePolynomial_46(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_46 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_46(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 46 ----------------------


class HermitePolynomial_47 : public HermitePolynomials {
public:

    HermitePolynomial_47(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_47 : public HermitePolynomials {
public:

    dell_HermitePolynomial_47(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_47 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_47(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 47 ----------------------


class HermitePolynomial_48 : public HermitePolynomials {
public:

    HermitePolynomial_48(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_48 : public HermitePolynomials {
public:

    dell_HermitePolynomial_48(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_48 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_48(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 48 ----------------------


class HermitePolynomial_49 : public HermitePolynomials {
public:

    HermitePolynomial_49(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_49 : public HermitePolynomials {
public:

    dell_HermitePolynomial_49(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_49 : public HermitePolynomials {
public:

    lapl_HermitePolynomial_49(double alpha, double omega);
    virtual double eval(double x);

};


//---------------------- END 49 ----------------------


#endif // HERMITEPOLYNOMIALS_H
