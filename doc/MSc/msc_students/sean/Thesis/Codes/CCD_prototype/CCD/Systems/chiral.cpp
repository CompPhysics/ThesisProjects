#include "chiral.h"
#include <iostream>

using namespace std;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

extern "C" {
    // get real and imaginary components of matrix element.
    void chipot_f90_wrapper_(double *matel_real, double *matel_im,
                            int *Nparticles, double *rho,
                            int *ps, int *pt, int *ppx, int *ppy, int* ppz,
                            int *qs, int *qt, int *qpx, int *qpy, int *qpz,
                            int *rs, int *rt, int *rpx, int *rpy, int *rpz,
                            int *ss, int *st, int *spx, int *spy, int *spz);
}

CHIRAL::CHIRAL(class Master* master, double m, double L3, double L2, double L1) : System(master) /* Chiral Potential */
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

void CHIRAL::makeStateSpace(){
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
                    } //end of nx^2 + ny^2 + nz^2 == n^2
                }//end of nz-loop
            }//end of ny-loop
        }//end of nx-loop
    }//end of n2-loop

    m_master->m_Ns = m_states.cols();
    below_fermi = Eigen::VectorXi::LinSpaced(m_Nh,0,m_Nh);
    above_fermi = Eigen::VectorXi::LinSpaced(m_Ns,m_Nh,m_Ns);
}

//I think using eigen here is a bit over-the-top for such a function, but whatevs~
int CHIRAL::kUnique1(int k, int s1){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::VectorXi mom = s1*kk;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

int CHIRAL::kUnique2(int k, int p, int s1, int s2){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::Vector4i kp( m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p) );
    Eigen::VectorXi mom = s1*kk + s2*kp;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

int CHIRAL::kUnique3(int k, int p, int q, int s1, int s2, int s3){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::Vector4i kp( m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p) );
    Eigen::Vector4i kq( m_states(1,q), m_states(2,q), m_states(3,q), m_states(4,q) );
    Eigen::VectorXi mom = s1*kk + s2*kp + s3*kq;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

int CHIRAL::kUnique4(int k, int p, int q, int s, int s1, int s2, int s3, int s4){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::Vector4i kp( m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p) );
    Eigen::Vector4i kq( m_states(1,q), m_states(2,q), m_states(3,q), m_states(4,q) );
    Eigen::Vector4i ks( m_states(1,s), m_states(2,s), m_states(3,s), m_states(4,s) );
    Eigen::VectorXi mom = s1*kk + s2*kp + s3*kq + s4*ks;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

int CHIRAL::kUnique5(int k, int p, int q, int s, int t, int s1, int s2, int s3, int s4, int s5){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::Vector4i kp( m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p) );
    Eigen::Vector4i kq( m_states(1,q), m_states(2,q), m_states(3,q), m_states(4,q) );
    Eigen::Vector4i ks( m_states(1,s), m_states(2,s), m_states(3,s), m_states(4,s) );
    Eigen::Vector4i kt( m_states(1,t), m_states(2,t), m_states(3,t), m_states(4,t) );
    Eigen::VectorXi mom = s1*kk + s2*kp + s3*kq + s4*ks + s5*kt;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

System::variable_type CHIRAL::f(int p){
    variable_type returnVal = h0(p);
    for (int i=0; i<m_Nh; i++){
        //if (i!=p){returnVal += assym_single(p, i);}
        returnVal += assym_single(p, i);
    };
    //cout << returnVal << endl;
    return returnVal;
}

System::variable_type CHIRAL::h0(int p){
    variable_type energy = (variable_type)m_states(0,p);
   // cout << energy*2.*pi*pi*m_hbarc*m_hbarc/(m_m*m_L2) << endl;
    return 2.*energy*pi*pi*m_hbarc*m_hbarc/(m_m*m_L2);
}

System::variable_type CHIRAL::assym(int p, int q, int r, int s){
    Eigen::Vector3i kp( m_states(1,p), m_states(2,p), m_states(3,p) );
    Eigen::Vector3i kq( m_states(1,q), m_states(2,q), m_states(3,q) );
    Eigen::Vector3i kr( m_states(1,r), m_states(2,r), m_states(3,r) );
    Eigen::Vector3i ks( m_states(1,s), m_states(2,s), m_states(3,s) );
    int sp = m_states(4,p); int tp = 1;//= m_states(5,p);
    int sq = m_states(4,q); int tq = 1;//= m_states(5,q);
    int sr = m_states(4,r); int tr = 1;//= m_states(5,r);
    int ss = m_states(4,s); int ts = 1;//= m_states(5,s);

    //these tests should be already performed through k_unique
    if ( vecDelta(kp+kq, kr+ks) == 0){ return 0;}   //momentum conservation
    if ( sp+sq != sr+ss){ return 0; }               //spin conservation
    if ( tp+tq != tr+ts){ return 0; }               //isospin conservation

    double matel_real   = 0;
    double matel_im     = 0;
    double rho          = m_Nh/(m_L3);
    int isospin         = 1;
    variable_type result;

    /*std::cout << kp(0) << " " << kp(1) << " " << kp(2) << " " << sp << std::endl;
    std::cout << kq(0) << " " << kq(1) << " " << kq(2) << " " << sq << std::endl;
    std::cout << std::endl;*/

    //Morten thinks +1=neutrons in the fortran code
    chipot_f90_wrapper_(&matel_real, &matel_im,
                        &m_Nh, &rho,
                        &m_states(4,p), &isospin, &m_states(1,p), &m_states(2,p), &m_states(3,p),
                        &m_states(4,q), &isospin, &m_states(1,q), &m_states(2,q), &m_states(3,q),
                        &m_states(4,r), &isospin, &m_states(1,r), &m_states(2,r), &m_states(3,r),
                        &m_states(4,s), &isospin, &m_states(1,s), &m_states(2,s), &m_states(3,s));

    double matel_real1   = 0;
    double matel_im1     = 0;
    variable_type result1;
    chipot_f90_wrapper_(&matel_real1, &matel_im1,
                        &m_Nh, &rho,
                        &m_states(4,p), &isospin, &m_states(1,p), &m_states(2,p), &m_states(3,p),
                        &m_states(4,q), &isospin, &m_states(1,q), &m_states(2,q), &m_states(3,q),
                        &m_states(4,s), &isospin, &m_states(1,s), &m_states(2,s), &m_states(3,s),
                        &m_states(4,r), &isospin, &m_states(1,r), &m_states(2,r), &m_states(3,r));

    /*int nxp = kp(0); int nyp = kp(1); int nzp = kp(2);
    int nxq = kq(0); int nyq = kq(1); int nzq = kq(2);
    int nxr = kr(0); int nyr = kr(1); int nzr = kr(2);
    int nxs = ks(0); int nys = ks(1); int nzs = ks(2);

    chipot_f90_wrapper_(&matel_real, &matel_im,
                            &m_Nh, &rho,
                            &sp, &isospin, &nxp, &nyp, &nzp,
                            &sq, &isospin, &nxq, &nyq, &nzq,
                            &sr, &isospin, &nxr, &nyr, &nzr,
                            &ss, &isospin, &nxs, &nys, &nzs);
    chipot_f90_wrapper_(&matel_real1, &matel_im1,
                            &m_Nh, &rho,
                            &sp, &isospin, &nxp, &nyp, &nzp,
                            &sq, &isospin, &nxq, &nyq, &nzq,
                            &ss, &isospin, &nxs, &nys, &nzs,
                            &sr, &isospin, &nxr, &nyr, &nzr);*/


    //function as in wrapper
    /*chipot_f90_wrapper_(&matel_real, &matel_im,
                        &Nparticles, &rho,
                        &ps, &pt, &ppx, &ppy, &ppz,
                        &qs, &qt, &qpx, &qpy, &qpz,
                        &rs, &rt, &rpx, &rpy, &rpz,
                        &ss, &st, &spx, &spy, &spz);*/
    //double real = sgn(matel_real)*pow(sgn(matel_real)*matel_real, 1./3);
    //double imag = sgn(matel_im)*pow(sgn(matel_im)*matel_im, 1./3);


    //std::cout << real << " " << matel_real << std::endl;

    //result.real(real);
    //result.imag(imag);

    /*result.real(matel_real);
    result.imag(matel_im);

    result1.real(matel_real1);
    result1.imag(matel_im1);*/

    return (result)*m_hbarc*m_hbarc/(2.*m_m);
}

System::variable_type CHIRAL::assym_single(int p, int q){
    /*if (m_states(1,p) == m_states(1,q)){
    std::cout << "i: "<< m_states(1,p) << " " << m_states(2,p) << " " << m_states(3,p) << std::endl;
    std::cout << "j: "<< m_states(1,q) << " " << m_states(2,q) << " " << m_states(3,q) << std::endl;
    std::cout << assym(p,q,p,q) << std::endl;
    std::cout << std::endl;
    }*/
    return assym(p,q,p,q);
}

bool CHIRAL::vecDelta(Eigen::VectorXi v1, Eigen::VectorXi v2){
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
