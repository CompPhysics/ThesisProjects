#include "mp.h"
#include <iostream>

using namespace std;

MP::MP(Master* master, double m, double L3, double L2, double L1) : System(master) /* Minnesota Potential */
{
    m_Nh = master->m_Nh;
    m_Nb = master->m_Nb;
    m_dk = 2*m_Nb + 1;
    m_master = master;

    m_m  = m;
    m_L3 = L3;
    m_L2 = L2;
    m_L1 = L1;
    V_0R_fac = V_0R*pow(pi/kappa_R, 1.5)/m_L3;
    V_0S_fac = -V_0S*pow(pi/kappa_S, 1.5)/m_L3; //the minus signs come from morten's code
    V_0T_fac = -V_0T*pow(pi/kappa_T, 1.5)/m_L3;
    piOverL  = pi/m_L1;
    makeStateSpace();
}

void MP::makeStateSpace(){
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
int MP::kUnique1(int k, int s1){
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
int MP::kUnique2(int k, int p, int s1, int s2){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::Vector4i kp( m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p) );
    Eigen::VectorXi mom = s1*kk + s2*kp;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

int MP::kUnique3(int k, int p, int q, int s1, int s2, int s3){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::Vector4i kp( m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p) );
    Eigen::Vector4i kq( m_states(1,q), m_states(2,q), m_states(3,q), m_states(4,q) );
    Eigen::VectorXi mom = s1*kk + s2*kp + s3*kq;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

int MP::kUnique4(int k, int p, int q, int s, int s1, int s2, int s3, int s4){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::Vector4i kp( m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p) );
    Eigen::Vector4i kq( m_states(1,q), m_states(2,q), m_states(3,q), m_states(4,q) );
    Eigen::Vector4i ks( m_states(1,s), m_states(2,s), m_states(3,s), m_states(4,s) );
    Eigen::VectorXi mom = s1*kk + s2*kp + s3*kq + s4*ks;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

int MP::kUnique5(int k, int p, int q, int s, int t, int s1, int s2, int s3, int s4, int s5){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::Vector4i kp( m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p) );
    Eigen::Vector4i kq( m_states(1,q), m_states(2,q), m_states(3,q), m_states(4,q) );
    Eigen::Vector4i ks( m_states(1,s), m_states(2,s), m_states(3,s), m_states(4,s) );
    Eigen::Vector4i kt( m_states(1,t), m_states(2,t), m_states(3,t), m_states(4,t) );
    Eigen::VectorXi mom = s1*kk + s2*kp + s3*kq + s4*ks + s5*kt;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}
/*int MP::kUnique1(int k, int s1){
    Eigen::Matrix<int, 5, 1> kk;
    kk << m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k), m_states(5,k);
    Eigen::VectorXi mom = s1*kk;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk + mom(4)*m_dk*m_dk*m_dk*m_dk;
    return kuni;
}

int MP::kUnique2(int k, int p, int s1, int s2){
    Eigen::Matrix<int, 5, 1> kk;
    Eigen::Matrix<int, 5, 1> kp;
    kk << m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k), m_states(5,k);
    kp << m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p), m_states(5,p);
    Eigen::VectorXi mom = s1*kk+s2*kp;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk + mom(4)*m_dk*m_dk*m_dk*m_dk;
    return kuni;
}

int MP::kUnique3(int k, int p, int q, int s1, int s2, int s3){
    Eigen::Matrix<int, 5, 1> kk;
    Eigen::Matrix<int, 5, 1> kp;
    Eigen::Matrix<int, 5, 1> kq;
    kk << m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k), m_states(5,k);
    kp << m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p), m_states(5,p);
    kq << m_states(1,q), m_states(2,q), m_states(3,q), m_states(4,q), m_states(5,q);
    Eigen::VectorXi mom = s1*kk+s2*kp+s3*kq;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk + mom(4)*m_dk*m_dk*m_dk*m_dk;
    return kuni;
}


int MP::kUnique4(int k, int p, int q, int s, int s1, int s2, int s3, int s4){
    Eigen::Matrix<int, 5, 1> kk;
    Eigen::Matrix<int, 5, 1> kp;
    Eigen::Matrix<int, 5, 1> kq;
    Eigen::Matrix<int, 5, 1> ks;
    kk << m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k), m_states(5,k);
    kp << m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p), m_states(5,p);
    kq << m_states(1,q), m_states(2,q), m_states(3,q), m_states(4,q), m_states(5,q);
    ks << m_states(1,s), m_states(2,s), m_states(3,s), m_states(4,s), m_states(5,s);
    Eigen::VectorXi mom = s1*kk+s2*kp+s3*kq+s4*ks;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk + mom(4)*m_dk*m_dk*m_dk*m_dk;
    return kuni;
}

int MP::kUnique5(int k, int p, int q, int s, int t, int s1, int s2, int s3, int s4, int s5){
    Eigen::Matrix<int, 5, 1> kk;
    Eigen::Matrix<int, 5, 1> kp;
    Eigen::Matrix<int, 5, 1> kq;
    Eigen::Matrix<int, 5, 1> ks;
    Eigen::Matrix<int, 5, 1> kt;
    kk << m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k), m_states(5,k);
    kp << m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p), m_states(5,p);
    kq << m_states(1,q), m_states(2,q), m_states(3,q), m_states(4,q), m_states(5,q);
    ks << m_states(1,s), m_states(2,s), m_states(3,s), m_states(4,s), m_states(5,s);
    kt << m_states(1,t), m_states(2,t), m_states(3,t), m_states(4,t), m_states(5,t);
    Eigen::VectorXi mom = s1*kk+s2*kp+s3*kq+s4*ks+s5*kt;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk + mom(4)*m_dk*m_dk*m_dk*m_dk;
    return kuni;
}*/

System::variable_type MP::f(int p){
    variable_type returnVal = h0(p);
    for (int i=0; i<m_Nh; i++){
        returnVal += assym_single(p, i);
    };
    return returnVal;
}

System::variable_type MP::h0(int p){
    double energy = m_states(0,p);
    return energy*2*pi*pi*m_hbarc*m_hbarc/(m_m*m_L2);
}

/*
double MP::assym(int p, int q, int r, int s){
    Eigen::Vector3i kp( m_states(1,p), m_states(2,p), m_states(3,p) );
    Eigen::Vector3i kq( m_states(1,q), m_states(2,q), m_states(3,q) );
    Eigen::Vector3i kr( m_states(1,r), m_states(2,r), m_states(3,r) );
    Eigen::Vector3i ks( m_states(1,s), m_states(2,s), m_states(3,s) );
    int sp = m_states(4,p); int tp = m_states(5,p);
    int sq = m_states(4,q); int tq = m_states(5,q);
    int sr = m_states(4,r); int tr = m_states(5,r);
    int ss = m_states(4,s); int ts = m_states(5,s);

    //these tests should be already performed through k_unique
    if ( vecDelta(kp+kq, kr+ks) == 0){ return 0;}   //momentum conservation
    if ( sp+sq != sr+ss){ return 0; }               //spin conservation
    if ( tp+tq != tr+ts){ return 0; }               //isospin conservation

    //I'm not certain if this indexing is correct
    //I've dropped the factor of 0.5 since I think it might be a misprint
    //piOverL is needed for momentum, not just quantum numbers
    double q2_dir = piOverL*((kp-kq-kr+ks).adjoint()*(kp-kq-kr+ks))(0,0);
    double q2_ex  = piOverL*((kp-kq+kr-ks).adjoint()*(kp-kq+kr-ks))(0,0);

    double V_1R = V_0R_fac*exp(-q2_dir/(4*kappa_R));
    double V_1T = V_0T_fac*exp(-q2_dir/(4*kappa_T));
    double V_1S = V_0S_fac*exp(-q2_dir/(4*kappa_S));

    double V_2R = V_0R_fac*exp(-q2_ex /(4*kappa_R));
    double V_2T = V_0T_fac*exp(-q2_ex /(4*kappa_T));
    double V_2S = V_0S_fac*exp(-q2_ex /(4*kappa_S));

    double returnVal = 0;

    //copying morten's code, didn't help much
    bool alpha1 = (sp==ss && sq==sr);                       bool beta1 = (sp==sr && sq==ss);
    bool alpha2 = (tp==ts && tq==tr);                       bool beta2 = (tp==tr && tq==ts);
    bool alpha3 = (sp==sr && sq==ss && tp==tr && tq==ts);   bool beta3 = (sp==ss && sq==sr && tp==ts && tq==tr);
    bool alpha4 = alpha1*(tp==tr && tq==ts);                bool beta4 = beta1*(tp==ts && tq==tr);
    bool alpha5 = alpha1*alpha2;                            bool beta5 = beta1*beta2;
    bool alpha6 = alpha2*(sp==sr && sq==ss);                bool beta6 = beta2*(sp==ss && sq==sr);

    returnVal += V_1R*(alpha3-alpha5) + 0.5*V_1T*(alpha3+alpha4-alpha5-alpha6) + 0.5*V_1S*(alpha3-alpha4-alpha5+alpha6);
    returnVal -= V_2R*(beta3-beta5) + 0.5*V_2T*(beta3+beta4-beta5-beta6) + 0.5*V_2S*(beta3-beta4-beta5+beta6);

    return 0.5*returnVal;
}*/


// TEST FUNCTIONS START
double MP::assym_test(int i, int j, int r, int s){
      double VRfactor = V_0R/(m_L3)*pow(M_PI/kappa_R,1.5);
      double VTfactor = -V_0T/(m_L3)*pow(M_PI/kappa_T,1.5);
      double VSfactor = -V_0S/(m_L3)*pow(M_PI/kappa_S,1.5);

      double V_R1, V_T1, V_S1, V_R2, V_T2, V_S2;
      double kX1, kY1, kZ1, kX2, kY2, kZ2;
      double qSquared1, spinEx1, isoSpinEx1, qSquared2, spinEx2, isoSpinEx2;
      double IsIt1, PsIt1, PsPt1, IsPt1, IsIt2, PsIt2, PsPt2, IsPt2;

      if(m_states(1,i) + m_states(1,j) != m_states(1,r) + m_states(1,s)){ return 0.0; }
      if(m_states(2,i) + m_states(2,j) != m_states(2,r) + m_states(2,s)){ return 0.0; }
      if(m_states(3,i) + m_states(3,j) != m_states(3,r) + m_states(3,s)){ return 0.0; }
      if(m_states(4,i) + m_states(4,j) != m_states(4,r) + m_states(4,s)){ return 0.0; }
      //if(m_states(5,i) + m_states(5,j) != m_states(5,r) + m_states(5,s)){ return 0.0; }

      kX1 = piOverL * (m_states(1,i) - m_states(1,j) - m_states(1,r) + m_states(1,s));
      kY1 = piOverL * (m_states(2,i) - m_states(2,j) - m_states(2,r) + m_states(2,s));
      kZ1 = piOverL * (m_states(3,i) - m_states(3,j) - m_states(3,r) + m_states(3,s));

      kX2 = piOverL * (m_states(1,i) - m_states(1,j) - m_states(1,s) + m_states(1,r));
      kY2 = piOverL * (m_states(2,i) - m_states(2,j) - m_states(2,s) + m_states(2,r));
      kZ2 = piOverL * (m_states(3,i) - m_states(3,j) - m_states(3,s) + m_states(3,r));

      qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
      qSquared2 = kX2 * kX2 + kY2 * kY2 + kZ2 * kZ2;

      V_R1 = VRfactor * exp(-qSquared1/(4*kappa_R));
      V_T1 = VTfactor * exp(-qSquared1/(4*kappa_T));
      V_S1 = VSfactor * exp(-qSquared1/(4*kappa_S));

      V_R2 = VRfactor * exp(-qSquared2/(4*kappa_R));
      V_T2 = VTfactor * exp(-qSquared2/(4*kappa_T));
      V_S2 = VSfactor * exp(-qSquared2/(4*kappa_S));

      spinEx1 = spinExchangeTerm(m_states(4,i), m_states(4,j), m_states(4,r), m_states(4,s));
      isoSpinEx1 = 1;//spinExchangeTerm(m_states(5,i), m_states(5,j), m_states(5,r), m_states(5,s));

      spinEx2 = spinExchangeTerm(m_states(4,i), m_states(4,j), m_states(4,s), m_states(4,r));
      isoSpinEx2 = 1;//spinExchangeTerm(m_states(5,i), m_states(5,j), m_states(5,s), m_states(5,r));

      IsIt1 = kroneckerDelta(m_states(4,i), m_states(4,r)) * kroneckerDelta(m_states(4,j), m_states(4,s));/* *
              kroneckerDelta(m_states(5,i), m_states(5,r)) * kroneckerDelta(m_states(5,j), m_states(5,s));*/
      PsIt1 = spinEx1;// * kroneckerDelta(m_states(5,i), m_states(5,r)) * kroneckerDelta(m_states(5,j), m_states(5,s));
      PsPt1 = spinEx1 * isoSpinEx1;
      IsPt1 = kroneckerDelta(m_states(4,i), m_states(4,r))*kroneckerDelta(m_states(4,j), m_states(4,s)) * isoSpinEx1;

      IsIt2 = kroneckerDelta(m_states(4,i), m_states(4,s)) * kroneckerDelta(m_states(4,j), m_states(4,r));/* *
              kroneckerDelta(m_states(5,i), m_states(5,s)) * kroneckerDelta(m_states(5,j), m_states(5,r));*/
      PsIt2 = spinEx2;// * kroneckerDelta(m_states(5,i), m_states(5,s)) * kroneckerDelta(m_states(5,j), m_states(5,r));
      PsPt2 = spinEx2 * isoSpinEx2;
      IsPt2 = kroneckerDelta(m_states(4,i), m_states(4,s)) * kroneckerDelta(m_states(4,j), m_states(4,r)) * isoSpinEx2;

      return 0.5 * (V_R1 + 0.5*V_T1 + 0.5*V_S1) * IsIt1 +
              0.25 * (V_T1 - V_S1) * PsIt1 -
              0.5 * (V_R1 + 0.5*V_T1 + 0.5*V_S1) * PsPt1 -
              0.25 * (V_T1 - V_S1) * IsPt1 -
              0.5 * (V_R2 + 0.5*V_T2 + 0.5*V_S2) * IsIt2 -
              0.25 * (V_T2 - V_S2) * PsIt2 +
              0.5 * (V_R2 + 0.5*V_T2 + 0.5*V_S2) * PsPt2 +
              0.25 * (V_T2 - V_S2) * IsPt2;
}

int MP::kroneckerDelta(const int &i, const int &j) {
  if(i != j){ return 0; }
  return 1;
}

int MP::spinExchangeTerm(const int &i, const int &j, const int &k, const int &l) {
  if(i == l && j == k){ return 1; }
  else{ return 0; }
}
//TEST FUNCTIONS END

System::variable_type MP::assym(int p, int q, int r, int s){
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

    double q2_dir = piOverL*piOverL*(kp-kq-kr+ks).squaredNorm();   // <pq||rs> (direct)
    double q2_ex = piOverL*piOverL*(kp-kq+kr-ks).squaredNorm();    // <pq||sr> (exchange)

   // cout << q2_ex << endl;

    double VR_dir = V_0R_fac*exp(-q2_dir/(4.0*kappa_R));
    double VT_dir = V_0T_fac*exp(-q2_dir/(4.0*kappa_T));
    double VS_dir = V_0S_fac*exp(-q2_dir/(4.0*kappa_S));

    double VR_ex = V_0R_fac*exp(-q2_ex/(4.0*kappa_R));
    double VT_ex = V_0T_fac*exp(-q2_ex/(4.0*kappa_T));
    double VS_ex = V_0S_fac*exp(-q2_ex/(4.0*kappa_S));


    bool Ps_dir = (sp==ss)*(sq==sr);    //exchange spins
    bool Cs_dir = (sp==sr)*(sq==ss);    //direct spins
    bool Pt_dir = (tp==ts)*(tq==tr);    //exchange isospins
    bool Ct_dir = (tp==tr)*(tq==ts);    //direct isospins

    bool Ps_ex = (sp==sr)*(sq==ss);    //exchange spins
    bool Cs_ex = (sp==ss)*(sq==sr);    //direct spins
    bool Pt_ex = (tp==tr)*(tq==ts);    //exchange isospins
    bool Ct_ex = (tp==ts)*(tq==tr);    //direct isospins

    double returnVal = 0;

    returnVal += (VR_dir + 0.5*VT_dir + 0.5*VS_dir)*Cs_dir*Ct_dir
               + 0.5*(VT_dir - VS_dir)*Ps_dir*Ct_dir
               - (VR_dir + 0.5*VT_dir + 0.5*VS_dir)*Ps_dir*Pt_dir
               - 0.5*(VT_dir - VS_dir)*Cs_dir*Pt_dir;

    returnVal += -(VR_ex + 0.5*VT_ex + 0.5*VS_ex)*Cs_ex*Ct_ex
               - 0.5*(VT_ex - VS_ex)*Ps_ex*Ct_ex
               + (VR_ex + 0.5*VT_ex + 0.5*VS_ex)*Ps_ex*Pt_ex
               + 0.5*(VT_ex - VS_ex)*Cs_ex*Pt_ex;

    //use to compare with morten's code
    /*double temp = assym_test(p,q,r,s);
    double diff = abs(temp - 0.5*returnVal);
    if (diff > 1e-14){std::cout << diff << std::endl;}*/

    return 0.5*returnVal;
    //return assym_test(p,q,r,s);
}

System::variable_type MP::assym_single(int p, int q){
    return assym(p,q,p,q);
}

//can do a test by copying assym, and set int s = p, int r = q
//worked for Nh=Nb=2 and Nh=14, Nb=3
/*double MP::assym_single(int p, int q){
    Eigen::Vector3i kp( m_states(1,p), m_states(2,p), m_states(3,p) ); //r = p
    Eigen::Vector3i kq( m_states(1,q), m_states(2,q), m_states(3,q) ); //s = q
    int sp = m_states(4,p); int tp = m_states(5,p);
    int sq = m_states(4,q); int tq = m_states(5,q);

    //I'm not certain if this indexing is correct
    //I've dropped the factor of 0.5 since I think it might be a misprint
    //piOverL is needed for momentum, not just quantum numbers
    //double q2_dir = 0;
    double q2_ex  = piOverL*((kp-kq).adjoint()*(kp-kq))(0,0);

    double V_1R = V_0R_fac;
    double V_1T = V_0T_fac;
    double V_1S = V_0S_fac;

    double V_2R = V_0R_fac*exp(-q2_ex /kappa_R);
    double V_2T = V_0T_fac*exp(-q2_ex /kappa_T);
    double V_2S = V_0S_fac*exp(-q2_ex /kappa_S);

    double returnVal = 0;

    //copying morten's code, didn't help much
    bool alpha1 = (sp==sq);       bool beta1 = true;
    bool alpha2 = (tp==tq);       bool beta2 = true;
    bool alpha3 = true;           bool beta3 = (sp==sq && tp==tq);
    bool alpha4 = alpha1;         bool beta4 = beta1*(tp==tq);
    bool alpha5 = alpha1*alpha2;  bool beta5 = beta1*beta2;
    bool alpha6 = alpha2;         bool beta6 = beta2*(sp==sq);

    returnVal += V_1R*(alpha3-alpha5) + 0.5*V_1T*(alpha3+alpha4-alpha5-alpha6) + 0.5*V_1S*(alpha3-alpha4-alpha5+alpha6);
    returnVal -= V_2R*(beta3-beta5) + 0.5*V_2T*(beta3+beta4-beta5-beta6) + 0.5*V_2S*(beta3-beta4-beta5+beta6);

    return 0.5*returnVal;
}*/

bool MP::vecDelta(Eigen::VectorXi v1, Eigen::VectorXi v2){
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
