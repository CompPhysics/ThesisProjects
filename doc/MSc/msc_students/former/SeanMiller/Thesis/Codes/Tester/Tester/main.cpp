#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <eigen3/Eigen/Dense>

using namespace std;

int main()
{
    int nParticles = 14;
    int SPB_dim = 10;                       //dimension (dim) of single-particle basis (SPB)
    double pi = M_PI;

    int N = nParticles;

    ofstream SPB ("SPB.dat");
    if (N==2){
        SPB << "0 0 0 0" << endl;       //n^2 nx ny nz
    }
    else if (N==14){
        SPB << "0 0 0 0" << endl;
        SPB << "1 1 0 0" << endl;
        SPB << "1 0 1 0" << endl;
        SPB << "1 0 0 1" << endl;
        SPB << "1 -1 0 0" << endl;
        SPB << "1 0 -1 0" << endl;
        SPB << "1 0 0 -1" << endl;
    }
    else cout << "Unable to open file";
    SPB.close();
    SPB.clear();

    double m    = 0.51;           //MeV
    double hbar = 0.1973*1e-12;   //MeV um
    double L    = 0.0005;         //m

    ofstream OBI ("OBI.dat");     //OBI: one-body interactions
    SPB.open("SPB.dat");
    while ( !SPB.eof() ){
        int nSquare
        SPB >> nSquare;
        double epsilon = 2*hbar*hbar*pi*nSquare/(m*L*L);    //epsilon_i = <i|h|i>
        OBI << epsilon << endl;
    }
}

