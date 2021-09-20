#ifndef MAKEINTMAT_H
#define MAKEINTMAT_H

#include "eigen3/Eigen/Dense"
#include <Systems/system.h>
#include <Systems/heg.h>
#include <Systems/mp.h>
#include <iostream>
#include <map>
#include <complex>
#include <unordered_map>
#include <sparsepp/spp.h>

class MakeIntMat
{
private:

    typedef Eigen::Matrix<short int, Eigen::Dynamic, Eigen::Dynamic>         MatrixXsi;

    //index-keepers
    Eigen::VectorXi identify_pppp;
    Eigen::VectorXi identify_hhpp;
    Eigen::VectorXi identify_fock;
    //memory storage
    MatrixXsi intMat_pppp;
    MatrixXsi intMat_hhpp;
    MatrixXsi fockMat;

    int m_printLength = 15;

    //bool contractor(int i, int j){ return i==j; } //contracts repeated elements to a single edit
public:

    //this is the type to be used in the interactions and amplitudes
    //set to std::complex<double> if you wish to work with the chiral model, otherwise use double
    typedef double variable_type;
    typedef Eigen::Matrix<variable_type, Eigen::Dynamic, Eigen::Dynamic> MatrixX;

    int m_Nh        = 0;
    int m_Ns        = 0;
    int m_Np        = 0;
    int m_Nh2       = 0;
    int m_Nh3       = 0;
    int m_Nh4       = 0;
    int m_Nh5       = 0;
    int m_NhNs      = 0;
    int m_NhNs2     = 0;
    unsigned long int m_Nh2Ns     = 0;
    unsigned long int m_Nh3Ns     = 0;
    unsigned long int m_Nh3Ns2    = 0;

    //Eigen matrix with unsigned long int as type
    typedef Eigen::Matrix<unsigned long int, Eigen::Dynamic, Eigen::Dynamic> MatrixXuli;
    //typedef Eigen::Matrix<short int, Eigen::Dynamic, Eigen::Dynamic>         MatrixXsi;


    /* Due to the momentum-conservation relation (kp+kq = kr+ks for <pq||rs>),
     * and due to alignment of matrices, we need several blockArrays, and corresponding
     * sortVecs and indexHolders. We'd like the blockArray generation to be pretty general,
     * otherwise the generation of blockArrays would require a nightmare of if-tests and whatnot.
     * So, since the CC diagrams require many different alignments, giving various rewritings of
     * the momentum relation, I've used the convention _(p/m)_(p/h), where in the first, p/m stands for
     * plus or minus between p/h in the next part. I.e. _ppm_hpp means total momentum kh+kp-kp.
     * Probably there's a better way, but I can't muck about with a little problem for long.
     */

    //blockArrays hold quantum numbers
    MatrixXsi blockArrays_p_h;
    MatrixXsi blockArrays_p_p;
    MatrixXsi blockArrays_pp_hh;
    MatrixXsi blockArrays_pp_hp;
    MatrixXsi blockArrays_pp_ph;
    MatrixXsi blockArrays_pp_pp;
    MatrixXsi blockArrays_hp_s;   //for Vhphp
    MatrixXsi blockArrays_pm_hh;
    MatrixXsi blockArrays_pm_hp;
    MatrixXsi blockArrays_pm_ph;
    MatrixXsi blockArrays_pm_pp;
    MatrixXsi blockArrays_ppp_hhh;
    MatrixXsi blockArrays_ppm_hhh;
    MatrixXsi blockArrays_ppm_hhp;
    MatrixXsi blockArrays_ppm_pph;
    MatrixXsi blockArrays_ppp_ppp;
    MatrixXsi blockArrays_pppm_hhhp;
    MatrixXsi blockArrays_ppmm_hhpp;
    MatrixXsi blockArrays_ppmm_pphh;
    MatrixXsi blockArrays_pppm_ppph;
    MatrixXsi blockArrays_pppmm_hhhpp;
    MatrixXsi blockArrays_pppmm_ppphh;

    //T3 permutations, needed to restack t_ijk^abc matrices
    //same shape as corresponding blockArrays
    MatrixXsi blockArrays_ppp_hhh_Pij;    // ijk <-> jik
    MatrixXsi blockArrays_ppp_hhh_Pik;    // ijk <-> kji
    MatrixXsi blockArrays_ppp_hhh_Pjk;    // ijk <-> ikj
    MatrixXsi blockArrays_ppp_hhh_Pijik;  // ijk <-> kij
    MatrixXsi blockArrays_ppp_hhh_Pijjk;  // ijk <-> jki

    MatrixXsi blockArrays_ppp_ppp_Pab;    // abc <-> bac
    MatrixXsi blockArrays_ppp_ppp_Pac;    // abc <-> cba
    MatrixXsi blockArrays_ppp_ppp_Pbc;    // abc <-> acb
    MatrixXsi blockArrays_ppp_ppp_Pabac;  // abc <-> cab
    MatrixXsi blockArrays_ppp_ppp_Pabbc;  // abc <-> bca

    void makePermutations(); //fills all permutaion matrices

    //sortVec holds all distinct kUnique for each index series
    std::vector<int> sortVec_p_h;
    std::vector<int> sortVec_p_p;
    std::vector<int> sortVec_pp_hh;
    std::vector<int> sortVec_pp_hp;
    std::vector<int> sortVec_pp_ph;
    std::vector<int> sortVec_pp_pp;
    std::vector<int> sortVec_hp_s;      //for Vhphp
    std::vector<int> sortVec_pm_hh;
    std::vector<int> sortVec_pm_hp;
    std::vector<int> sortVec_pm_ph;
    std::vector<int> sortVec_pm_pp;
    std::vector<int> sortVec_ppp_hhh;
    std::vector<int> sortVec_ppm_hhh;
    std::vector<int> sortVec_ppm_hhp;
    std::vector<int> sortVec_ppm_pph;
    std::vector<int> sortVec_ppp_ppp;
    std::vector<int> sortVec_pppm_hhhp;
    std::vector<int> sortVec_ppmm_hhpp;
    std::vector<int> sortVec_ppmm_pphh;
    std::vector<int> sortVec_pppm_ppph;
    std::vector<int> sortVec_pppmm_hhhpp;
    std::vector<int> sortVec_pppmm_ppphh;

    std::vector<int> fullVec_p_h;
    std::vector<int> fullVec_p_p;
    std::vector<int> fullVec_pp_hh;
    std::vector<int> fullVec_pp_hp;
    std::vector<int> fullVec_pp_ph;
    std::vector<int> fullVec_pp_pp;
    std::vector<int> fullVec_hp_s;      //for Vhphp
    std::vector<int> fullVec_pm_hh;
    std::vector<int> fullVec_pm_hp;
    std::vector<int> fullVec_pm_ph;
    std::vector<int> fullVec_pm_pp;
    std::vector<int> fullVec_ppp_hhh;
    std::vector<int> fullVec_ppm_hhh;
    std::vector<int> fullVec_ppm_hhp;
    std::vector<int> fullVec_ppm_pph;
    std::vector<int> fullVec_ppp_ppp;
    std::vector<int> fullVec_pppm_hhhp;
    std::vector<int> fullVec_ppmm_hhpp;
    std::vector<int> fullVec_ppmm_pphh;
    std::vector<int> fullVec_pppm_ppph;
    std::vector<int> fullVec_pppmm_hhhpp;
    std::vector<int> fullVec_pppmm_ppphh;

    //indexHolder holds upper and lower bound of indices for a certain kUnique, same indexing as the corresponding matrices
    MatrixXuli indexHolder_p_h;
    MatrixXuli indexHolder_p_p;
    MatrixXuli indexHolder_pp_hh;
    MatrixXuli indexHolder_pp_hp;
    MatrixXuli indexHolder_pp_ph;
    MatrixXuli indexHolder_pp_pp;
    MatrixXuli indexHolder_hp_s;   //for Vhphp
    MatrixXuli indexHolder_pm_hh;
    MatrixXuli indexHolder_pm_hp;
    MatrixXuli indexHolder_pm_ph;
    MatrixXuli indexHolder_pm_pp;
    MatrixXuli indexHolder_ppp_hhh;
    MatrixXuli indexHolder_ppm_hhh;
    MatrixXuli indexHolder_ppm_hhp;
    MatrixXuli indexHolder_ppm_pph;
    MatrixXuli indexHolder_ppp_ppp;
    MatrixXuli indexHolder_pppm_hhhp;
    MatrixXuli indexHolder_ppmm_hhpp;
    MatrixXuli indexHolder_ppmm_pphh;
    MatrixXuli indexHolder_pppm_ppph;
    MatrixXuli indexHolder_pppmm_hhhpp;
    MatrixXuli indexHolder_pppmm_ppphh;

    // these additional boundsHolders are possibly not necessary after the implementation of make3x1- and make2x2Block
    //these are needed for Qa
    MatrixXuli boundsHolder_hhpp_hh;
    MatrixXuli boundsHolder_hhpp_pp;
    MatrixXuli boundsHolder_hhhppp_hhh;
    MatrixXuli boundsHolder_hhhppp_ppp;

    //these are needed for Qc and Qd
    MatrixXsi indexHolder_h_pph_hpp;
    MatrixXsi indexHolder_h_pph_h;
    MatrixXsi indexHolder_hhp_p_hhp;
    MatrixXsi indexHolder_hhp_p_p;

    //MatrixXsi indexHolder_hp;
    //MatrixXsi indexHolder_ph;
    //MatrixXsi boundsHolder_pppp_pp;

    void setThreads(int numthreads);
    int m_numThreads;

    MakeIntMat();
    System* m_system = nullptr;
    void                            mapper_1(std::vector<int> &sortVecIn, std::vector<int> &fullVecIn, MatrixXsi &blockArraysIn, int i1, int s1);
    void                            mapper_2(std::vector<int> &sortVecIn, std::vector<int> &fullVecIn, MatrixXsi &blockArraysIn, int i1, int i2, int s1, int s2);
    void                            mapper_2_alt();
    void                            mapper_3(std::vector<int> &sortVecIn, std::vector<int> &fullVecIn, MatrixXsi &blockArraysIn, int i1, int i2, int i3, int s1, int s2, int s3);
    void                            mapper_4(std::vector<int> &sortVecIn, std::vector<int> &fullVecIn, MatrixXsi &blockArraysIn, int i1, int i2, int i3, int i4, int s1, int s2, int s3, int s4);
    void                            mapper_5(std::vector<int> &sortVecIn, std::vector<int> &fullVecIn, MatrixXsi &blockArraysIn, int i1, int i2, int i3, int i4, int i5, int s1, int s2, int s3, int s4, int s5);
    void                            mapper_hp();        //special func made for Lc diagram, not necessary IF I fix sign convention for ph and hp

    void                            makeBlockMat(System* system, int Nh, int Ns);
    void                            setTriples(bool argument);
    bool                            m_triplesOn;

    int                             numOfKu;        //number of blocks in V_hh_pp
    int                             numOfKu3;       //number of blocks in V_hhh_ppp
    MatrixX                 makeSquareBlock(MatrixXsi& array, unsigned long int range_lower, unsigned long int range_upper);
    MatrixX                 makeSquareBlock_s(MatrixXsi& array, unsigned long range_lower, unsigned long range_upper);
    MatrixX                 makeRektBlock(MatrixXsi& array1, MatrixXsi& array2, unsigned long int range_lower1, unsigned long int range_upper1, unsigned long int range_lower2, unsigned long int range_upper2);


    //CCD terms
    MatrixX                 I1_makemat(int channel1, int channel2);
    MatrixX                 I2_makemat(int channel1, int channel2);
    MatrixX                 I3_makemat(int channel1, int channel2);
    MatrixX                 I4_makemat(int channel1, int channel2);

    //T3 contributions to T2
    MakeIntMat::MatrixX                 D10b_makemat(int channel1, int channel2);
    MakeIntMat::MatrixX                 D10c_makemat(int channel1, int channel2);

    //linear T2 terms in T3
    MakeIntMat::MatrixX                 T1a_makemat(int channel1, int channel2);
    MakeIntMat::MatrixX                 T1b_makemat(int channel1, int channel2);

    //linear T3 terms in T3 (these already stored directly, and need no "make" function)

    //quadratic T2 terms in T3
    MakeIntMat::MatrixX                 T3b_makemat(int channel1, int channel2);
    MakeIntMat::MatrixX                 T3c_makemat(int channel1, int channel2);
    MakeIntMat::MatrixX                 T3d_makemat(int channel1, int channel2);
    MakeIntMat::MatrixX                 T3e_makemat(int channel1, int channel2);

    //T2*T3 terms in T3
    MakeIntMat::MatrixX                 T5a_makemat(int channel1, int channel2);
    MakeIntMat::MatrixX                 T5b_makemat(int channel1, int channel2);
    MakeIntMat::MatrixX                 T5c_makemat(int channel1, int channel2);
    MakeIntMat::MatrixX                 T5d_makemat(int channel1, int channel2);
    MakeIntMat::MatrixX                 T5e_makemat(int channel1, int channel2);
    MakeIntMat::MatrixX                 T5f_makemat(int channel1, int channel2); //T5f and T5g are identical, but it's nice to think of them seperate
    MakeIntMat::MatrixX                 T5g_makemat(int channel1, int channel2);

    MakeIntMat::MatrixX                 make3x1Block(int ku, int i1, int i2, int i3, int i4);
    MakeIntMat::MatrixX                 make2x2Block(int ku, int i1, int i2, int i3, int i4);
    MakeIntMat::MatrixX                 make2x2Block_alt(int channel);

    void                            makeMatMap_hhhp(MatrixXsi& array1, MatrixXsi& array2, unsigned long int range_lower1, unsigned long int range_upper1, unsigned long int range_lower2, unsigned long int range_upper2);
    void                            makeMatMap_hhpp(MatrixXsi& array1, MatrixXsi& array2, unsigned long int range_lower1, unsigned long int range_upper1, unsigned long int range_lower2, unsigned long int range_upper2);
    void                            makeMatMap_ppph(MatrixXsi& array1, MatrixXsi& array2, unsigned long int range_lower1, unsigned long int range_upper1, unsigned long int range_lower2, unsigned long int range_upper2);

    void                            makeMatVec(MatrixXsi& array1, MatrixXsi& array2, unsigned long int range_lower1, unsigned long int range_upper1, unsigned long int range_lower2, unsigned long int range_upper2);
    unsigned long int               Identity_hhhp(int h1, int h2, int h3, int p1);
    unsigned long int               Identity_hhpp(int h1, int h2, int p1, int p2);
    unsigned long int               Identity_hh(int h1, int h2);
    unsigned long int               Identity_pp(int p1, int p2);
    unsigned long int               Identity_ppph(int p1, int p2, int p3, int h1);
    unsigned long int               Identity_hhhppp(int h1, int h2, int h3, int p1, int p2, int p3);
    unsigned long int               Identity_hhh(int h1, int h2, int h3);
    unsigned long int               Identity_ppp(int p1, int p2, int p3);
    unsigned long int               Identity_hhhhhp(int h1, int h2, int h3, int h4, int h5, int p1); //this is a special function for T3e alone
    spp::sparse_hash_map<unsigned long int, variable_type>           Vhhhp_elements; //needed for T3
    spp::sparse_hash_map<unsigned long int, variable_type>           Vhhpp_elements;
    spp::sparse_hash_map<unsigned long int, variable_type>           Vppph_elements; //needed for T3

    /*std::pair<unsigned long, int> returnId_hhhp(const int i1, const int i2, const int i3, const int i4);
    int returnId_hhpp(const int i1, const int i2, const int i3, const int i4);
    int returnId_hppp(const int i1, const int i2, const int i3, const int i4);*/

    //Eigen::VectorXd                 Vhhpp_vector;
    std::vector<variable_type>             Vhhpp_vector;
    int                             Vhhpp_counter = 0;

    //interaction matrices for CCD
    std::vector<MakeIntMat::MatrixX>   Vhhpp;

    //these are special
    std::vector<MakeIntMat::MatrixX>   Vpppp; //for La
    std::vector<MakeIntMat::MatrixX>   Vhhhh; //for Lb
    std::vector<MakeIntMat::MatrixX>   Vhphp; //for Lc

    std::vector<MatrixXsi>   V_hp_hp;
    std::vector<MatrixXsi>   V_hh_pp;

    //vectors with kUnique for each matrix, indices match ( that is, Vhhhh_i[h] <-> Vhhhh[h] )
    std::vector<int>               Vhhhh_i;
    std::vector<int>               Vhphp_i;
    std::vector<int>               Vhhpp_i;
    std::vector<int>               Vpppp_i;
    std::vector<int>               Vhhhppp_i;
};

#endif // MAKEINTMAT_H
