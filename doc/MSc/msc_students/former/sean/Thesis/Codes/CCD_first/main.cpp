#include <iostream>
#include <math.h>
#include <cmath>
#include <eigen3/Eigen/Dense>

using namespace std;

//TEMPLATES
template <typename T> //credit: http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes#12399290
vector<size_t> sort_indexes(const vector<T> &v) {

    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}

//CLASSES

//CLASS HEG
class HEG{
private:
    Eigen::VectorXi below_fermi = Eigen::VectorXi::LinSpaced(m_Nh,0,m_Nh);
    Eigen::VectorXi above_fermi;
    Eigen::MatrixXi states;                         //each column is a state where the rows are quantum numbers of that states

    bool vecDelta(Eigen::VectorXi v1, Eigen::VectorXi v2){
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

    double  pi        =   M_PI;
    double  m_m       =   0;                          //electron mass [MeV]
    double  m_L3      =   0;                          //box volume
    double  m_L2      =   0;
    double  m_L1      =   0;

public:

    int            m_Ns                 = 0;
    int            m_Nh                 = 0;        //number of particles
    int            m_Nb                 = 0;        //number of closed-shells (n^2=0, n^2=1, n^2=2, etc... For NB=2 can at max have N=14)
    void           setNumOfHoles(int Nh)    {m_Nh = Nh;}
    void           setNumOfShells(int Nb)   {m_Nb = Nb;}

    void setConstants(double m, double L3, double L2, double L1){
        m_m  = m;
        m_L3 = L3;
        m_L2 = L2;
        m_L1 = L1;
    }

    void makeStateSpace(){
        states.conservativeResize(states.rows()+5, Eigen::NoChange);    //sets the format of states

        //start for-loops
        for (int n2=0; n2<m_Nb; n2++){
            for (int nx=-n2; nx<n2+1; nx++){
                for (int ny=-n2; ny<n2+1; ny++){
                    for (int nz=-n2; nz<n2+1; nz++){
                        if (nx*nx + ny*ny + nz*nz == n2){
                            states.conservativeResize(Eigen::NoChange, states.cols()+2);
                            states.col(states.cols()-2) << n2,nx,ny,nz, 1;
                            states.col(states.cols()-1) << n2,nx,ny,nz,-1;
                        }
                    }
                }
            }
        }

        //end for-loops
        m_Ns = states.cols();
        above_fermi = Eigen::VectorXi::LinSpaced(m_Ns,m_Nh,m_Ns);
    }

    //I think using eigen here is a bit over-the-top for such a function, but whatevs~
    int kUnique2(int k, int p){
        Eigen::VectorXi kk = states.col(k);
        Eigen::VectorXi kp = states.col(p);
        Eigen::VectorXi mom = kk+kp;

        int val = 0;
        for (int i = 0; i<mom.rows();i++){
            if (val < mom(i)){
                val = mom(i);
            }
        }

        int dk = 2*val + 1;
        int kuni = mom(1) + mom(2)*dk + mom(3)*dk*dk + mom(4)*dk*dk*dk;
        return kuni;
    }

    double f(int p){
        double returnVal = h0(p);
        for (int i=0; i<below_fermi.rows(); i++){
            returnVal += assym(p,below_fermi(i),p,below_fermi(i));
        };
        return returnVal;
    }

    double h0(int p){
        double energy = states(p,0);
        return energy*2*pi*pi/(m_m*m_L2);
    }

    //Eigen::MatrixXi* assym(int kp, int kq, int kr, int ks){
    double assym(int p, int q, int r, int s){
        int sp = states(p,4);
        Eigen::Vector3i kp( states(p,1), states(p,2), states(p,3) );
        int sq = states(q,4);
        Eigen::Vector3i kq( states(q,1), states(q,2), states(q,3) );
        int sr = states(r,4);
        Eigen::Vector3i kr( states(r,1), states(r,2), states(r,3) );
        int ss = states(s,4);
        Eigen::Vector3i ks( states(s,1), states(s,2), states(s,3) );

        if ( vecDelta(kp+kq, kr+ks) ){
            double returnVal = 0;
            if ( vecDelta( kp, kr) ){
                returnVal += (sp==sr)*(sq==ss)/( (kr-kp).squaredNorm() );
            }
            if ( vecDelta( kp, ks) ){
                returnVal -= (sp==ss)*(sq==sr)/( (ks-kp).squaredNorm() );
            }
            return returnVal/(m_L1*pi);
        }
        else{
            return 0;
        }
    }

};

//CLASS DIAGRAMS
class Diagrams{
private:
public:
    Eigen::MatrixXi D1(){}
    /*Eigen::MatrixXi D2(){}
    Eigen::MatrixXi D3(){}
    Eigen::MatrixXi D4(){}
    Eigen::MatrixXi D5(){}
    Eigen::MatrixXi D6(){}
    Eigen::MatrixXi D7(){}
    Eigen::MatrixXi D8(){}
    Eigen::MatrixXi D9(){}*/
};

//CLASS MAKEINTMAT
class MakeIntMat{
private:
    int m_Nh = 0;
    int m_Ns = 0;
    //index-keepers
    Eigen::VectorXi identify_pppp;
    Eigen::VectorXi identify_hhpp;
    Eigen::VectorXi identify_fock;
    //memory storage
    Eigen::MatrixXi intMat_pppp;
    Eigen::MatrixXi intMat_hhpp;
    Eigen::MatrixXi fockMat;
    Eigen::Matrix<int, 3, (m_Ns-m_Nh)*(m_Ns-m_Nh)>   blockArrays_pp;
    Eigen::Matrix<int, 3, m_Nh*m_Nh>             blockArrays_hh;
    bool contractor(int i, int j){ return i==j; } //contracts repeated elements to a single edit
public:
    class HEG*     m_system             = nullptr;
    void setNumOfStates(){m_Ns =

    vector<Eigen::MatrixXf*> makeBlockMat(){
        Eigen::Matrix<int, 3, (m_Ns-m_Nh)*(m_Ns-m_Nh)>   blockArrays_pp_temp;
        Eigen::Matrix<int, 3, m_Nh*m_Nh>             blockArrays_hh_temp;

        //make index-arrays
        int index = 0;
        for (a=m_Nh; a<m_Ns; a++){
            for (b=m_Nh; b<m_Ns; b++){
                blockArrays_pp_temp.col(index) = system::kUnique2(a,b),a,b;
            }
        }
        index = 0;
        for (i=0; i<m_Nh; i++){
            for (j=0; j<m_Nh; j++){
                blockArrays_hh_temp.col(index) = system::kUnique2(i,j),i,j;
            }
        }
        //end of making index-arrays

        //sorting
        index = 0;
        Eigen::VectorXi tempVec_pp = blockArrays_pp_temp.row(0);
        //credit: http://stackoverflow.com/questions/26094379/typecasting-eigenvectorxd-to-stdvector#26094708
        vector<int> sortVec_pp(tempVec_pp.data(), tempVec_pp.data() + blockArrays_pp_temp.cols()); //this vector now contains all k_unique

        //index i is the sorted indices of temp_vec
        for (auto i: sort_indexes(sortVec_pp)) {
            blockArrays_pp.col(index) << blockArrays_pp_temp.col(i) << endl;
            index += 1;
        }
        //end of sorting

        //blocking
        std::vector<int>::iterator it;
        it=std::unique_copy (sortVec_pp.begin(), it, sortVec_pp.begin(), contractor);
        sortVec_pp.resize( std::distance(sortVec_pp.begin(),it) ); //now sortvec_pp has no repeated k_unique values

        std::vector<Eigen::MatrixXf*> returnBlocks;

        int range_upper = 0;
        int range_lower = 0;

        for (int l=0; l<sortVec_pp.size(); l++){
            val = sortVec_pp[l];
            for (int h=0; h<m_Ns-m_Nh; h++){
                if ( val == blockArrays_pp(0,h) ){
                    range_upper += 1;
                }
            }
            returnBlocks.push_back( makeBlock(range_lower, range_upper) );//NEED MATRIX GENERATOR HERE
        }
        //end of blocking

        return returnBlocks;
    }

    Eigen::MatrixXf* makeBlock(int range_lower, int range_upper){
        int dim = range_upper - range_lower;
        Eigen::MatrixXf<double, dim,dim> returnMat;
        for (int i = range_lower; i<range_upper; i++){
            for (int j = range_lower; j<range_upper; j++){
                returnMat(i,j) = system.assym(blockArrays_pp(1,i),blockArrays_pp(2,i),blockArrays_pp(1,j),blockArrays_pp(2,j));
            }
        }
        return returnMat*;
    }

};

//CLASS MAKEAMPMAT
class MakeAmpMat{
private:
    //index-keepers
    Eigen::VectorXi identify_amp;   //same form as identify_hhpp
    //memory storage
    Eigen::MatrixXi T2;
public:
    Eigen::MatrixXi* makeMat();
};

//CLASS MASTER
class Master{
public:
    class HEG*     m_system             = nullptr;
    int            m_Nh                 = 0;        //number of particles
    int            m_Nb                 = 0;        //number of closed-shells (n^2=0, n^2=1, n^2=2, etc... For NB=2 can at max have N=14)
    int            m_Ns                 = 0;
    void           setNumOfHoles(int Nh)    {m_Nh = Nh;}
    void           setNumOfShells(int Nb)   {m_Nb = Nb;}
    void           setNumOfStates()         {m_Ns = m_system->m_Ns;}

    void setSystem(class HEG* system){
        m_system = system;
    }

    double Iterator(double eps, double conFac){
        double ECCD;
        double ECCD_old;
        while (conFac < eps){

            ECCD_old =
            conFac = abs(ECCD - ECCD_old);
        }
        return ECCD;
    }
};

//MAIN
int main()
{
    double  pi      =   M_PI;
    int     Nh      =   2;							//number of particles
    int     Nb      =   4;							//number of closed-shells (n^2=0, n^2=1, n^2=2, etc... For NB=2 can at max have N=14)
    double  rs      =   1.0;                        //Wigner Seitz radius
    double  rb      =   1;                          //Bohr radius [MeV^-1]
    double  m       =   1;                          //electron mass [MeV]
    double  L3      =   4*pi*Nh*rs/3;               //box volume
    double  L2      =   pow(L3, 2/3);
    double  L1      =   pow(L3, 1/3);

    double  eps     =   1e-10;
    double  conFac  =   1;                          //convergence factor

    HEG*    heg    = new HEG;
    heg->setConstants(m, L3, L2, L1);

    Master* master = new Master;
    master->setNumOfHoles(Nh);
    master->setNumOfShells(Nb);
    master->setSystem(heg);

    master->Iterator(eps, conFac);


}

