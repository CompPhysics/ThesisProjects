#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <eigen3/Eigen/Dense>
#include <armadillo>


using namespace std;


double singleBody(int i, arma::mat SPB){
    return 0.5*SPB(i,0);
}

int KD(int p, int q){
    return (p==q);  //return 1 if true, 0 if false
}

//Coloumb interaction. The "v..." notation stands for vector, and the arguments are 2d vectors where the first
//argument is the momentum and the second is the spin, i.e. vAlpha = (k_alpha, m_alpha)
double coloumbAS(int q, int p, int r, int s, double L3, int d, arma::mat SPB){

    double returnVal = 0;
    double pi = M_PI;

    bool cond1 =   KD(SPB(q,0)+SPB(p,0), SPB(r,0)+SPB(s,0));
    bool cond2 = ( KD(SPB(q,4),SPB(r,4)) && KD(SPB(p,4),SPB(s,4)) && KD(SPB(q,0),SPB(r,0)) );
    bool cond3 = ( KD(SPB(q,4),SPB(s,4)) && KD(SPB(p,4),SPB(r,4)) && KD(SPB(q,0),SPB(s,0)) );


    if ( cond1 ){
        if ( cond2 ){
            if      (d==2){returnVal += 1/( SPB(q,0) - SPB(r,0) );}
            else if (d==3){returnVal += 1/( ( SPB(q,0) - SPB(r,0) )*( SPB(q,0) - SPB(r,0) ) );}
        }
        if ( cond3 ){
            if      (d==2){returnVal -= 1/( SPB(q,0) - SPB(s,0) );}
            else if (d==3){returnVal -= 1/( ( SPB(q,0) - SPB(s,0) )*( SPB(q,0) - SPB(s,0) ) );}
        }
    }

    return (1/(pi*cbrt(L3)))*returnVal;
}

double calcEHF(double rho, int N, int d, arma::mat SPB){

    double L3 = N/rho;
    double EHF = 0;

    for (int q=0; q<N; q++){
        for (int p=0; p<N; p++){
            if (q==p){EHF += singleBody(q, SPB);}    //this can actually moved up a for-loop and drop the if-test

            for (int k=0; k<N; k++){
                EHF += 0.5*coloumbAS(q, k, p, k, L3, d, SPB);
            }
        }
    }

    return EHF;
}

int main()
{
    int nParticles = 14;
    double pi   = M_PI;
    int N = nParticles;

    arma::mat SPB = arma::zeros<arma::mat>(N, 5); //arma uses (row,col) indexing

    if (N==2){
        SPB(0,0) = 0; SPB(0,1) = 0; SPB(0,2) = 0; SPB(0,3) = 0; SPB(0,4) = 1;
        SPB(1,0) = 0; SPB(1,1) = 0; SPB(1,2) = 0; SPB(1,3) = 0; SPB(1,4) = -1;
    }
    else if (N==14){
        SPB(0,0) = 0; SPB(0,1) = 0; SPB(0,2) = 0; SPB(0,3) = 0; SPB(0,4) = 1;
        SPB(1,0) = 0; SPB(1,1) = 0; SPB(1,2) = 0; SPB(1,3) = 0; SPB(1,4) = -1;
        SPB(2,0) = 1; SPB(2,1) = 1; SPB(2,2) = 0; SPB(2,3) = 0; SPB(2,4) = 1;
        SPB(3,0) = 1; SPB(3,1) = 1; SPB(3,2) = 0; SPB(3,3) = 0; SPB(3,4) = -1;
        SPB(4,0) = 1; SPB(4,1) = 0; SPB(4,2) = 1; SPB(4,3) = 0; SPB(4,4) = 1;
        SPB(5,0) = 1; SPB(5,1) = 0; SPB(5,2) = 1; SPB(5,3) = 0; SPB(5,4) = -1;
        SPB(6,0) = 1; SPB(6,1) = 0; SPB(6,2) = 0; SPB(6,3) = 1; SPB(6,4) = 1;
        SPB(7,0) = 1; SPB(7,1) = 0; SPB(7,2) = 0; SPB(7,3) = 1; SPB(7,4) = -1;
        SPB(8,0) = 1; SPB(8,1) = -1; SPB(8,2) = 0; SPB(8,3) = 0; SPB(8,4) = 1;
        SPB(9,0) = 1; SPB(9,1) = -1; SPB(9,2) = 0; SPB(9,3) = 0; SPB(9,4) = -1;
        SPB(10,0) = 1; SPB(10,1) = 0; SPB(10,2) = -1; SPB(10,3) = 0; SPB(10,4) = 1;
        SPB(11,0) = 1; SPB(11,1) = 0; SPB(11,2) = -1; SPB(11,3) = 0; SPB(11,4) = -1;
        SPB(12,0) = 1; SPB(12,1) = 0; SPB(12,2) = 0; SPB(12,3) = -1; SPB(12,4) = 1;
        SPB(13,0) = 1; SPB(13,1) = 0; SPB(13,2) = 0; SPB(13,3) = -1; SPB(13,4) = -1;
    }

    double upperRho = 0.3;
    int d = 3;
    int numSteps = 100;

    ofstream myfile1;
    ofstream myfile2;
    myfile1.open ("example1.txt");
    myfile2.open ("example2.txt");
    for (int i=1; i<numSteps; i++){
        double dens = i*upperRho/numSteps;
        double rs = cbrt(3/(4*pi*dens));
        double var = calcEHF(dens, N, d, SPB);
        double ana = 0.5*(2.21/(rs*rs) - 0.916/rs);
        myfile1 << dens << "     " << var << endl;
        myfile2 << dens << "     " << ana << endl;
    }
    myfile1.close();
    myfile2.close();

}
