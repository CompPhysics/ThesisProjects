#ifndef SYSTEM_H
#define SYSTEM_H

#include <master.h>
#include <complex>
//#include "../makeintmat.h"
#include "eigen3/Eigen/Dense"

class System
{
public:
    typedef double variable_type;
    System(class Master* master);
    int                    m_Nh;        //number of particles
    int                    m_Nb;        //number of closed-shells
    int                    m_Ns;
    int                    m_dk;
    const double    m_hbarc = 197.32697188;  //[Mev fm] (/c)
    Eigen::VectorXi below_fermi;
    Eigen::VectorXi above_fermi;
    Eigen::MatrixXi m_states;                         //each column is a state where the rows are quantum numbers of that states
    class Master*   m_master        = nullptr;
    virtual void    makeStateSpace  ()                              = 0;
    virtual variable_type  assym           (int p, int q, int r, int s)    = 0;
    virtual variable_type  assym_single    (int p, int q)                  = 0;
    virtual variable_type  h0              (int p)                         = 0;
    virtual variable_type  f               (int p)                         = 0;
    virtual int     kUnique1        (int p, int s1)                                                             = 0;
    virtual int     kUnique2        (int p, int q, int s1, int s2)                                              = 0;
    virtual int     kUnique3        (int p, int q, int r, int s1, int s2, int s3)                               = 0;
    virtual int     kUnique4        (int p, int q, int r, int s, int s1, int s2, int s3, int s4)                = 0;
    virtual int     kUnique5        (int p, int q, int r, int s, int t, int s1, int s2, int s3, int s4, int s5) = 0;
};

#endif // SYSTEM_H
