from mpmath import *
from sympy.utilities.codegen import codegen
from sympy import *

mp.dps = 15
mp.pretty = True
nMax = 50

printExpressions = True
printConstructors = False
printDeclarations = False
printVectorElements = False

if (printExpressions):
	x, a, w, aw = symbols('x, m_alpha, m_omega, m_alphaomega')

	H = zeros(nMax)
	Hd = zeros(nMax)
	Hdd = zeros(nMax)

	H[0] = 1

	for n in range(1, nMax):
		H[n] = simplify(2*sqrt(a*w)*x*H[n-1] - diff(H[n-1], x))

	for n in range(0, nMax):
		Hd[n] = simplify(diff(H[n], x))

	for n in range(0, nMax):
		Hdd[n] = simplify(diff(H[n], x, 2))

	#for n in range(0, nMax):
	#	print codegen(("HermitePolynomial_%i::eval" %n, H[n]), "c", "file", header=False)[0][1]
	#	print codegen(("dell_HermitePolynomial_%i::eval" %n, Hd[n]), "c", "file", header=False)[0][1]
	#	print codegen(("lapl_HermitePolynomial_%i::eval" %n, Hdd[n]), "c", "file", header=False)[0][1]

	for n in range(0, nMax):
		print("""double HermitePolynomial_%i::eval(double x) {

   double H;
   H = %s;
   return H;

}




double dell_HermitePolynomial_%i::eval(double x) {

   double dell_H;
   dell_H = %s;
   return dell_H;

}




double lapl_HermitePolynomial_%i::eval(double x) {

   double lapl_H;
   lapl_H = %s;
   return lapl_H;

}""" %(n, printing.ccode(H[n]), n, printing.ccode(Hd[n]), n, printing.ccode(Hdd[n])))
	print("""

//---------------------- END %i ----------------------

""" %n)


if (printConstructors):
	for n in range(0, nMax):
		print("""HermitePolynomial_%i::HermitePolynomial_%i(double alpha, double omega)
: HermitePolynomials() { 
    m_alpha = alpha;
    m_omega = omega;
    m_alphaomega = alpha*omega;
}

dell_HermitePolynomial_%i::dell_HermitePolynomial_%i(double alpha, double omega)
: HermitePolynomials() { 
    m_alpha = alpha;
    m_omega = omega;
    m_alphaomega = alpha*omega;
}

lapl_HermitePolynomial_%i::lapl_HermitePolynomial_%i(double alpha, double omega)
: HermitePolynomials() { 
    m_alpha = alpha;
    m_omega = omega;
    m_alphaomega = alpha*omega;
}""" %(n, n, n, n, n, n))
		print("""

//---------------------- END %i ----------------------

""" %n)

if (printDeclarations):
	for n in range(0, nMax):
		print("""class HermitePolynomial_%i : public HermitePolynomials {
public:

    HermitePolynomial_%i(double alpha, double omega);
    virtual double eval(double x);

};


class dell_HermitePolynomial_%i : public HermitePolynomials {
public:

    dell_HermitePolynomial_%i(double alpha, double omega);
    virtual double eval(double x);

};


class lapl_HermitePolynomial_%i : public HermitePolynomials {
public:

    lapl_HermitePolynomial_%i(double alpha, double omega);
    virtual double eval(double x);

};""" %(n, n, n, n, n, n))
		print("""

//---------------------- END %i ----------------------

""" %n)

if (printVectorElements):
	for n in range(0, nMax):
		print("""m_hermitePolynomials[%i] = new HermitePolynomial_%i(m_alpha, m_omega);""" %(n,n))
	print("")

	for n in range(0, nMax):
		print("""m_hermitePolynomialsDerivative[%i] = new dell_HermitePolynomial_%i(m_alpha, m_omega);""" %(n,n))
	print("")
	
	for n in range(0, nMax):
		print("""m_hermitePolynomialsDoubleDerivative[%i] = new lapl_HermitePolynomial_%i(m_alpha, m_omega);""" %(n,n))
	



















