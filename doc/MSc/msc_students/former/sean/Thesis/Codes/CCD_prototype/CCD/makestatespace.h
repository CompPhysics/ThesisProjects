#ifndef MAKESTATESPACE_H
#define MAKESTATESPACE_H

#include <eigen3/Eigen/Dense>

class MakeStateSpace
{
public:
    MakeStateSpace();
    int             Ns;
    Eigen::VectorXi below_fermi;
    Eigen::VectorXi above_fermi;
    Eigen::MatrixXi states;
};

#endif // MAKESTATESPACE_H
