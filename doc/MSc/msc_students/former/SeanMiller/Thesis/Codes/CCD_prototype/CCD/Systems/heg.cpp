#include "heg.h"
#include <iostream>

using namespace std;

HEG::HEG(Master* master, double m, double L3, double L2, double L1) : System(master) /* Homogenous Electron Gas */
{
    m_Nh = master->m_Nh;
    m_Nb = master->m_Nb;
    m_dk = 2*m_Nb + 1;
    m_master = master;

    m_m  = m;
    m_L3 = L3;
    m_L2 = L2;
    m_L1 = L1;
    makeStateSpace();
}

void HEG::makeStateSpace(){
    m_states.conservativeResize(m_states.rows()+5, Eigen::NoChange);    //sets the format of states
    //start for-loops
    for (int n2=0; n2<m_Nb+1; n2++){
        for (int nx=-n2; nx<n2+1; nx++){
            for (int ny=-n2; ny<n2+1; ny++){
                for (int nz=-n2; nz<n2+1; nz++){
                    if (nx*nx + ny*ny + nz*nz == n2){
                        m_states.conservativeResize(Eigen::NoChange, m_states.cols()+2);
                        m_states.col(m_states.cols()-2) << n2,nx,ny,nz, 1;
                        m_states.col(m_states.cols()-1) << n2,nx,ny,nz,-1;
                    }
                }
            }
        }
    }

    m_master->m_Ns = m_states.cols();
    below_fermi = Eigen::VectorXi::LinSpaced(m_Nh,0,m_Nh);
    above_fermi = Eigen::VectorXi::LinSpaced(m_Ns,m_Nh,m_Ns);
}

int HEG::kUnique1(int k, int s1){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::VectorXi mom = s1*kk;

    /*int val = 0;
    for (int i = 0; i<mom.rows();i++){
        if (val < mom(i)){
            val = mom(i);
        }
    }

    int dk = 2*val + 1;*/

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

//I think using eigen here is a bit over-the-top for such a function, but whatevs~
int HEG::kUnique2(int k, int p, int s1, int s2){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::Vector4i kp( m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p) );
    Eigen::VectorXi mom = s1*kk + s2*kp;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

int HEG::kUnique3(int k, int p, int q, int s1, int s2, int s3){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::Vector4i kp( m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p) );
    Eigen::Vector4i kq( m_states(1,q), m_states(2,q), m_states(3,q), m_states(4,q) );
    Eigen::VectorXi mom = s1*kk + s2*kp + s3*kq;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

int HEG::kUnique4(int k, int p, int q, int s, int s1, int s2, int s3, int s4){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::Vector4i kp( m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p) );
    Eigen::Vector4i kq( m_states(1,q), m_states(2,q), m_states(3,q), m_states(4,q) );
    Eigen::Vector4i ks( m_states(1,s), m_states(2,s), m_states(3,s), m_states(4,s) );
    Eigen::VectorXi mom = s1*kk + s2*kp + s3*kq + s4*ks;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

int HEG::kUnique5(int k, int p, int q, int s, int t, int s1, int s2, int s3, int s4, int s5){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::Vector4i kp( m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p) );
    Eigen::Vector4i kq( m_states(1,q), m_states(2,q), m_states(3,q), m_states(4,q) );
    Eigen::Vector4i ks( m_states(1,s), m_states(2,s), m_states(3,s), m_states(4,s) );
    Eigen::Vector4i kt( m_states(1,t), m_states(2,t), m_states(3,t), m_states(4,t) );
    Eigen::VectorXi mom = s1*kk + s2*kp + s3*kq + s4*ks + s5*kt;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

System::variable_type HEG::f(int p){
    variable_type returnVal = h0(p);
    for (int i=0; i<m_Nh; i++){
        returnVal += assym_single(p, i);
    };
    return returnVal;
}

System::variable_type HEG::h0(int p){
    double energy = m_states(0,p);
    return energy*2*pi*pi/(m_m*m_L2);
}

System::variable_type HEG::assym(int p, int q, int r, int s){

    Eigen::Vector3i kp( m_states(1,p), m_states(2,p), m_states(3,p) );
    Eigen::Vector3i kq( m_states(1,q), m_states(2,q), m_states(3,q) );
    Eigen::Vector3i kr( m_states(1,r), m_states(2,r), m_states(3,r) );
    Eigen::Vector3i ks( m_states(1,s), m_states(2,s), m_states(3,s) );
    int sp = m_states(4,p);
    int sq = m_states(4,q);
    int sr = m_states(4,r);
    int ss = m_states(4,s);

    if ( vecDelta(kp+kq, kr+ks) ){

        double returnVal = 0;
        if ( vecDelta( kp, kr) == 0 ){
            returnVal += (sp==sr)*(sq==ss)/( (double)(kr-kp).squaredNorm() );
        }
        if ( vecDelta( kp, ks) == 0){
            returnVal -= (sp==ss)*(sq==sr)/( (double)(ks-kp).squaredNorm() );
        }
        return returnVal/(m_L1*pi);
    }
    else{
        return 0;
    }
}

//can test this func by copying assym, and set int r = p, int s = q
System::variable_type HEG::assym_single(int p, int q){

    Eigen::Vector3i kp( m_states(1,p), m_states(2,p), m_states(3,p) );
    Eigen::Vector3i kq( m_states(1,q), m_states(2,q), m_states(3,q) );
    int sp = m_states(4,p);
    int sq = m_states(4,q);

    double returnVal = 0;
    if ( vecDelta( kp, kq) == 0){
        returnVal -= (sp==sq)/( (double)(kq-kp).squaredNorm() );
    }
    return returnVal/(m_L1*pi);
}

bool HEG::vecDelta(Eigen::VectorXi v1, Eigen::VectorXi v2){
    int dim1 = v1.rows();
    int dim2 = v2.rows();

    if (dim1 != dim2){
        cout << "dimensional error" << endl;
        return 0;
    }
    else{
        bool returnInt = 1;
        for (int i =0; i<dim1; i++){
            returnInt *= ( v1(i)==v2(i) );
        }
        return returnInt;
    }
}
