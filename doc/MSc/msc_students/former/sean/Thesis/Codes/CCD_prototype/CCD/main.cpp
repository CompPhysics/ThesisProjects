//imported files
//other libraries
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <fstream>
#include <math.h>
#include <chrono>
#include <omp.h>
#include "eigen3/Eigen/Dense"

typedef std::chrono::high_resolution_clock Clock;   //needed for timing

//author made files
#include "makestatespace.h"
#include "master.h"
#include "Systems/system.h"
#include "Systems/heg.h"
#include "Systems/mp.h"
#include "Systems/chiral.h"
#include "makeampmat.h"
#include "makeintmat.h"


using namespace std;

//argv will accept: System (string), Nh (int), Nb (int), rs/rho (double), precision (double), degree of triples (CCD, CCDT-1, CCDT-2, etc), #threads
int main(int argc, char** argv)
{
    Eigen::initParallel();

    double       eps     = 1e-16;              //remember to adjust setprecision in master when changing this
    double       conFac  = 1;                          //convergence factor
    const double pi      = M_PI;


    bool    intermediates = true;                   //turn on/off intermediates in CCD eqs
    //bool    CCDT          = true;                   //turn on/off CCDT-1
    bool    timer         = true;                   //turn on/off timer
    bool    relaxation    = true;                   //turn on/off relaxation when updating amplitudes
    double  alpha         = 0.873;                  //relaxation parameter (I found 0.873 to be best)
    int     threads       = 4;                      //number of threads, default is whatever you put here

    bool    threadsOn;
    if (threads > 1){threadsOn = true;}
    else{threadsOn = false;}

    if (atoi(argv[7]) > 1){
        threadsOn = true;

        Eigen::initParallel();
        threads = atoi(argv[7]);
        omp_set_dynamic(0);
        omp_set_num_threads(threads);
        Eigen::setNbThreads(threads);
        int n = Eigen::nbThreads( );
        std::cout << "OMP threads: " << n << std::endl;
    }
    else{
        omp_set_dynamic(0);
        omp_set_num_threads(threads);
            std::cout << "Threading is turned off, meaning:" << std::endl;
            std::cout << "- Eigen won't run parallel matrix products" << std::endl;
            std::cout << "- Permutations will be performed in serial" << std::endl;
            std::cout << "- Diagrams T5b and T5c will be run in serial" << std::endl;

    }

    if (argc != 8 && argc != 1){
        std::cout << "Failure: Missmatch command line arguments" << std::endl;
        std::cout << "The command line arguments are: Model, #particles, #shells, rs/rho, desired precision (down to 1e-16), which CC approximation" << std::endl;
        std::cout << "The CC models are: 0 -> full CCD" << std::endl;
        std::cout << "                   1 -> CCDT-1" << std::endl;
        std::cout << "                   2 -> CCDT-2" << std::endl;
        std::cout << "                   3 -> full CCDT" << std::endl;
        std::cout << "If you require any other set of CC diagrams, edit them manually in master.cpp (in function Iterator)" << std::endl;

        return 0;
    }

    //we use natural units
    int Nh; int Nb;
    if (argc==8){//user defined size
        Nh = atoi(argv[2]);				//number of particles
        Nb = atoi(argv[3]);				//number of closed-shells (n^2=0, n^2=1, n^2=2, etc... For NB=2 is min for N=14)
        eps = atof(argv[5]);
    }
    else{        //default size
        Nh = 14;                        //number of particles
        Nb = 6;                         //number of closed-shells (n^2=0, n^2=1, n^2=2, etc... For NB=2 is min for N=14)
    }
    double  rs;     //Wigner Seitz radius
    double  rho;    //Density
    double  L3;     //Box volume
    double  L2;
    double  L1;
    double  m;      //Particle mass

    Master* master = new Master;
    master->setSize(Nh, Nb);

    //set system and physical parameters
    if (argc != 8){
        std::cout << "No arguments given, running default setup" << std::endl;
        std::cout << "Default setup: HEG for Nh=14, Nb=3, rs=1.0, 1e-16 precision, CCDT" << std::endl;

        /*m          = 1;             //Electron mass [MeV?]
        rs         = 1;
        double  rb = 1.;            //Bohr radius [MeV^-1]
        double  r1 = pow(rs*rb, 3);
        L3         = 4.*pi*Nh*r1/3.;
        L2         = pow(L3, 2./3.);
        L1         = pow(L3, 1./3.);
        master->setSystem(new HEG(master, m, L3, L2, L1))*/;
        m   = 939.5653;              //Neutron mass [MeV]
        rho = 0.2;//7.9999998211860657E-002;
        L3  = double(Nh)/rho;
        L2  = pow(L3, 2./3.);
        L1  = pow(L3, 1./3.);
        master->setSystem(new MP(master, m, L3, L2, L1));
    }
    else if (std::string(argv[1]) == "HEG"){
        m          = 1;             //Electron mass [MeV?]
        rs         = atof(argv[4]);
        double  rb = 1.;            //Bohr radius [MeV^-1]
        double  r1 = pow(rs*rb, 3);
        L3         = 4.*pi*Nh*r1/3.;
        L2         = pow(L3, 2./3.);
        L1         = pow(L3, 1./3.);
        master->setSystem(new HEG(master, m, L3, L2, L1));
    }
    else if (std::string(argv[1]) == "MP"){
        m   = 939.565;              //Neutron mass [MeV]
        rho = atof(argv[4]);
        L3  = double(Nh)/rho;
        L2  = pow(L3, 2./3.);
        L1  = pow(L3, 1./3.);
        master->setSystem(new MP(master, m, L3, L2, L1));
    }
    else{
        std::cout << "Failure: You need to submit a model, e.g. HEG or MP" << std::endl;
        std::cout << "Which model is specified in the first command line argument" << std::endl;

        return 0;
    }

    bool CCDT;

    master->setIntermediates(intermediates);
    master->setRelaxation(relaxation, alpha);
    master->setTimer(timer);
    master->setThreads_forMaster(threadsOn, threads/*atoi(argv[7])*/);
    if (argc == 8){
        master->setCCType(atoi(argv[6]));
        if (atoi(argv[6])==0){
            master->setTriples(false);
            CCDT = false;
        }
        else{
            master->setTriples(true);
            CCDT = true;
        }
    }
    else{
        master->setCCType(3);
        master->setTriples(true);
        CCDT = true;
    }

    cout << "C++ code" << endl;

    double EHF; double Eref; double Ecc;
    Eref = master->CC_Eref();
    EHF = master->CC_E_HF();

    auto t1 = Clock::now();

    //typedef MakeIntMat::variable_type variable_type;

    bool run_CC = true;

    if (run_CC){
        if (CCDT){
            double ECCDT = master->CC_master(eps, conFac);
            std::cout << std::endl;
            std::cout << "Correlation energy" << std::endl;
            std::cout << "Delta ECCDT:    "<< std::right << std::setw(21) << ECCDT << std::endl;
            Ecc = ECCDT;
        }
        else{
            double ECCD = master->CC_master(eps, conFac);
            std::cout << std::endl;
            std::cout << "Correlation energy" << std::endl;
            std::cout << "Delta ECCD:     "<< std::right << std::setw(21) << ECCD << std::endl;
            Ecc = ECCD;
        }
        std::cout << "E_corr/A:       " << std::fixed << std::setprecision (16) << std::right << std::setw(21) << (Ecc)/Nh << std::endl;
    }

    auto t2 = Clock::now();

    std::cout << std::endl;
    std::cout << "Hartree-Fock energy" << std::endl;
    std::cout << "E_HF:           " << std::fixed << std::setprecision (16) << std::right << std::setw(21) << EHF << std::endl;
    std::cout << "E_HF/A:         " << std::fixed << std::setprecision (16) << std::right << std::setw(21) << (EHF)/Nh << std::endl;
    std::cout << std::endl;
    std::cout << "Reference energy" << std::endl;
    std::cout << "Eref:           " << std::fixed << std::setprecision (16) << std::right << std::setw(21) << Eref << std::endl;
    std::cout << "E_ref/A:        " << std::fixed << std::setprecision (16) << std::right << std::setw(21) << (Eref)/Nh << std::endl;
    std::cout << std::endl;
    std::cout << "Total energy" << std::endl;
    std::cout << "E/A:            " << std::fixed << std::setprecision (16) << std::right << std::setw(21) << (Eref+Ecc)/Nh << std::endl;
    std::cout << std::endl;

    if (intermediates){
        std::cout << "Total time used: "
                  << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
                  << " seconds, with intermediates ON" << std::endl;
    }
    else{
        std::cout << "Total time used: "
                  << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
                  << " seconds, with intermediates OFF" << std::endl;
    }

    return 0;

}
