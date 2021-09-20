#ifndef MASTER_H
#define MASTER_H

#include <Systems/system.h>
#include "eigen3/Eigen/Dense"


class Master
{
public:
    class System*     m_system             = nullptr;
    class MakeAmpMat* m_ampClass           = nullptr;
    class MakeIntMat* m_intClass           = nullptr;
    class Diagrams*   m_diagrams           = nullptr;
    int               m_Nh                 = 0;        //number of particles
    int               m_Nb                 = 0;
    int               m_Ns                 = 0;
    int               m_Np                 = 0;
    bool              m_triplesOn          = false;
    int               m_CC_type            = 0;
    bool              m_intermediatesOn    = false;
    bool              m_timerOn            = false;
    bool              m_relaxation         = false;
    bool              m_threadsOn          = false;
    int               m_threads            = 1;
    double            m_alpha              = 1;

    int countCC_iters = 0;

    //this is the type to be used in the interactions and amplitudes
    //set to std::complex if you wish to work with the chiral model, otherwise use double
    typedef double variable_type;
    typedef Eigen::Matrix<variable_type, Eigen::Dynamic, Eigen::Dynamic> MatrixX;

    typedef Eigen::Matrix<unsigned long int, Eigen::Dynamic, Eigen::Dynamic> MatrixXuli;

    void setSize();
    void setSystem(class System* system);
    void setThreads_forMaster(bool argument, int num);
    void setThreads();
    void setSize(int Nh, int Nb);
    void setIntermediates(bool argument);
    void setTriples(bool argument);
    void setTimer(bool argument);
    void setRelaxation(bool argument, double alpha);
    void setClasses();
    void setCCType(int type);

    double CC_Eref();
    double CC_E_HF();
    double CC_master(double eps, double conFac);
    variable_type Iterator(double eps, double conFac, variable_type E_MBPT2);
};

#endif // MASTER_H
