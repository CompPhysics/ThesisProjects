#ifndef MAKEAMPMAT_H
#define MAKEAMPMAT_H

#include "eigen3/Eigen/Dense"
#include <makeintmat.h>
#include <complex>
#include <Systems/system.h>
#include <Systems/heg.h>
#include <Systems/mp.h>
#include <Systems/pm.h>

class MakeAmpMat
{
public:
    MakeAmpMat();

    typedef MakeIntMat::variable_type   variable_type;
    typedef MakeIntMat::MatrixX         MatrixX;

    void makeFockMaps();
    void emptyFockMaps();
    spp::sparse_hash_map<int, variable_type> FockMap_h;
    spp::sparse_hash_map<int, variable_type> FockMap_p;
    std::vector<MatrixX> denomMat;
    std::vector<MatrixX> denomMat3;

    int m_Counter = 0; //a test counter, used for debugging sessions

    void setThreads(int numthreads);
    int m_numThreads;

    MakeIntMat*  m_intClass = nullptr;
    System*      m_system   = nullptr;
    std::vector<MatrixX>  Amplitudes;

    typedef Eigen::Matrix<unsigned long int, Eigen::Dynamic, Eigen::Dynamic> MatrixXuli;
    typedef Eigen::Matrix<short int, Eigen::Dynamic, Eigen::Dynamic>         MatrixXsi;

    spp::sparse_hash_map<unsigned long int, variable_type>           T2_elements;
    spp::sparse_hash_map<unsigned long int, variable_type>           T2_temp;
    spp::sparse_hash_map<unsigned long int, variable_type>           T2_elements_new;

    //these are no longer in use
    spp::sparse_hash_map<int, variable_type>           T3_elements;
    spp::sparse_hash_map<int, variable_type>           T3_temp;
    spp::sparse_hash_map<int, variable_type>           T3_elements_new;

    //The T5b and T5c use very demanding remappings, so it's far more efficient, both in CPU and memory,
    //to have matrices storing the indices, rather than finding them on the go.
    //These matrices are made in mapper_5 in the intClass, because this is a late fix
    std::vector<MatrixXuli>              T3_T5b_indices;
    std::vector<MatrixXsi>               T3_T5b_indices_signs;
    std::vector<MatrixXuli>              T3_T5c_indices;
    std::vector<MatrixXsi>               T3_T5c_indices_signs;

    int returnId(const int i1, const int i2, const int i3, const int i4);

    spp::sparse_hash_map<unsigned long int, unsigned long int>        T3_elements_I;

    std::vector<variable_type>        T3_elements_A;        //holds T3 amplitudes
    std::vector<variable_type>        T3_elements_A_new;    //holds new T3 amplitudes
    std::vector<variable_type>        denom3_elements;
    std::vector<variable_type>        T3_elements_A_temp;   //holds temporary diagram contributions
    void                              T3_makeMap(MatrixX inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6);
    void                              T3_makeDirectMat();
    MatrixX                           T3_buildDirectMat(int channel, std::vector<variable_type>& T_vec);
    std::vector<MatrixXuli>           T3_directMat;         //holds indices for T3_elements_A to make t_ijk^abc

    MatrixX                 make3x1Block(int ku, int i1, int i2, int i3, int i4, spp::sparse_hash_map<unsigned long int, variable_type> &T_list);
    MatrixX                 make2x2Block(int ku, int i1, int i2, int i3, int i4, spp::sparse_hash_map<unsigned long int, variable_type>& T_list);

    MatrixX                 make3x3Block(int ku, int i1, int i2, int i3, int i4, int i5, int i6, spp::sparse_hash_map<unsigned long int, variable_type>& T_list);
    MatrixX                 make3x3Block_I(int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<variable_type>& T_vec);
    MatrixX                 make3x3Block_I_D10c(int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<variable_type>& T_vec);



    //CCD terms
    MatrixX                 I1_makemat_1(int channel1, int channel2);
    MatrixX                 I1_makemat_2(int channel1, int channel2);
    MatrixX                 I2_makemat_1(int channel1, int channel2);
    MatrixX                 I2_makemat_2(int channel1, int channel2);
    MatrixX                 I3_makemat_1(int channel1, int channel2);
    MatrixX                 I3_makemat_2(int channel1, int channel2);
    MatrixX                 I4_makemat_1(int channel1, int channel2);
    MatrixX                 I4_makemat_2(int channel1, int channel2);

    //T3 contributions to T2
    MatrixX                 D10b_makemat(int channel1, int channel2);
    MatrixX                 D10c_makemat(int channel1, int channel2);

    //linear T2 terms in T3
    MatrixX                 T1a_makemat(int channel1, int channel2);
    MatrixX                 T1b_makemat(int channel1, int channel2);

    //linear T3 terms in T3
    MatrixX                 T2c_makemat(int channel1, int channel2);
    MatrixX                 T2d_makemat(int channel1, int channel2);
    MatrixX                 T2e_makemat(int channel1, int channel2);


    //quadratic T2 terms in T3
    MatrixX                 T3b_makemat_1(int channel1, int channel2);
    MatrixX                 T3b_makemat_2(int channel1, int channel2);
    MatrixX                 T3b_makemat_3(int channel1, int channel2);  //The convention makemat_1, makemat_2, and makemat_3 comes from the diagrams
    MatrixX                 T3c_makemat_1(int channel1, int channel2);  //there are three matrices in these, so you need three constructors
    MatrixX                 T3c_makemat_2(int channel1, int channel2);  //a pain, trust me, I know
    MatrixX                 T3c_makemat_3(int channel1, int channel2);
    MatrixX                 T3d_makemat_1(int channel1, int channel2);
    MatrixX                 T3d_makemat_2(int channel1, int channel2);
    MatrixX                 T3d_makemat_3(int channel1, int channel2);
    MatrixX                 T3e_makemat_1(int channel1, int channel2);
    MatrixX                 T3e_makemat_2(int channel1, int channel2);
    MatrixX                 T3e_makemat_3(int channel1, int channel2);

    //T2*T3 terms in T3
    MatrixX                 T5a_makemat_1(int channel1, int channel2);
    MatrixX                 T5a_makemat_2(int channel1, int channel2);
    MatrixX                 T5b_makemat_1(int channel1, int channel2);
    MatrixX                 T5b_makemat_2(int channel1, int channel2);
    MatrixXuli              T5b_makemat_2_I(int channel1, int channel2); //index, not double
    MatrixXsi               T5b_makemat_2_I_signs(int channel1, int channel2); //signs of indices, not double
    MatrixX                 T5c_makemat_1(int channel1, int channel2);
    MatrixX                 T5c_makemat_2(int channel1, int channel2);
    MatrixXuli              T5c_makemat_2_I(int channel1, int channel2); //index, not double
    MatrixXsi               T5c_makemat_2_I_signs(int channel1, int channel2); //signs of indices, not double
    MatrixX                 T5d_makemat_1(int channel1, int channel2);
    MatrixX                 T5d_makemat_2(int channel1, int channel2);
    MatrixX                 T5e_makemat_1(int channel1, int channel2);
    MatrixX                 T5e_makemat_2(int channel1, int channel2);
    MatrixX                 T5f_makemat_1(int channel1, int channel2);
    MatrixX                 T5f_makemat_2(int channel1, int channel2);
    MatrixX                 T5g_makemat_1(int channel1, int channel2);
    MatrixX                 T5g_makemat_2(int channel1, int channel2);

    //CCD terms
    void                            I1_inverse(MatrixX& inMat, int channel1, int channel2);
    void                            I2_inverse(MatrixX& inMat, int channel1, int channel2);
    void                            I3_inverse(MatrixX inMat, int channel1, int channel2);
    void                            I4_inverse(MatrixX inMat, int channel1, int channel2);

    //T3 contributions to T2
    void                            D10b_inverse(MatrixX& inMat, int channel1, int channel2);
    void                            D10c_inverse(MatrixX& inMat, int channel1, int channel2);

    //linear T2 terms in T3
    void                            T1a_inverse(MatrixX& inMat, int channel1, int channel2);
    void                            T1b_inverse(MatrixX& inMat, int channel1, int channel2);

    //linear T3 terms in T3
    void                            T2c_inverse(MatrixX& inMat, int channel1, int channel2);
    void                            T2d_inverse(MatrixX& inMat, int channel1, int channel2);
    void                            T2e_inverse(MatrixX& inMat, int channel1, int channel2);

    //since the T3b-e diagrams are more finicky than the rest, we need a remapper for the first product
    spp::sparse_hash_map<unsigned long int, variable_type> T3D_remap;
    void                            T3b_Inverse_temp(MatrixX& inMat, int channel1, int channel2);
    void                            T3c_Inverse_temp(MatrixX inMat, int channel1, int channel2);
    void                            T3d_Inverse_temp(MatrixX& inMat, int channel1, int channel2);
    void                            T3e_Inverse_temp(MatrixX& inMat, int channel1, int channel2);

    //quadratic T2 terms in T3
    void                            T3b_inverse(MatrixX& inMat, int channel1, int channel2);
    void                            T3c_inverse(MatrixX& inMat, int channel1, int channel2);
    void                            T3d_inverse(MatrixX& inMat, int channel1, int channel2);
    void                            T3e_inverse(MatrixX& inMat, int channel1, int channel2);

    //T2*T3 terms in T3
    void                            T5a_inverse(MatrixX& inMat, int channel1, int channel2);
    void                            T5b_inverse(MatrixX& inMat, int channel1, int channel2);
    void                            T5c_inverse(MatrixX& inMat, int channel1, int channel2);
    void                            T5b_inverse_I(MatrixX& inMat, int channel);
    void                            T5c_inverse_I(MatrixX& inMat, int channel);
    void                            T5d_inverse(MatrixX& inMat, int channel1, int channel2);
    void                            T5e_inverse(MatrixX& inMat, int channel1, int channel2);
    void                            T5f_inverse(MatrixX inMat, int channel1, int channel2);
    void                            T5g_inverse(MatrixX& inMat, int channel1, int channel2);



    void                            make3x1Block_inverse(MatrixX inMat, int ku, int i1, int i2, int i3, int i4, spp::sparse_hash_map<unsigned long, variable_type> &T_list, bool add);
    void                            make3x1Block_inverse_D10b(MatrixX inMat, int ku, int i1, int i2, int i3, int i4, spp::sparse_hash_map<unsigned long, variable_type> &T_list, bool add);
    void                            make2x2Block_inverse(MatrixX inMat, int ku, int i1, int i2, int i3, int i4, spp::sparse_hash_map<unsigned long, variable_type> &T_list, bool add);

    void                            make3x3Block_inverse(MatrixX inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6, spp::sparse_hash_map<unsigned long int, variable_type>& T_list, bool add);
    void                            make3x3Block_inverse_I(MatrixX inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<variable_type>& T_vec, bool add);
    void                            make3x3Block_inverse_I_T1a(MatrixX inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<variable_type>& T_vec, bool add);
    void                            make3x3Block_inverse_I_T1b(MatrixX inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<variable_type>& T_vec, bool add);

    void                            addElementsT2(bool Pij, bool Pab);
    void                            addElementsT3_T1a();
    void                            addElementsT3_T1b();
    void                            addElementsT3_T2c();
    void                            addElementsT3_T2d();
    void                            addElementsT3_T2e();
    void                            addElementsT3_T3b();
    void                            addElementsT3_T3c();
    void                            addElementsT3_T3d();
    void                            addElementsT3_T3e();
    void                            addElementsT3_T5a();
    void                            addElementsT3_T5b();
    void                            addElementsT3_T5c();
    void                            addElementsT3_T5d();
    void                            addElementsT3_T5e();
    void                            addElementsT3_T5f();
    void                            addElementsT3_T5g();
    void                            addElementsT3(bool Pij, bool Pik, bool Pjk, bool Pab, bool Pac, bool Pbc);

    spp::sparse_hash_map<int, int> permuteT3(int index, spp::sparse_hash_map<int, int> indices);

    void setIntClass(class MakeIntMat* intClass);
    void setSystem(class System* system);
    void setElements_T2();
    void setElements_T3();
    MatrixX makeBlockMat(int index);
    void makeDenomMat();
    void makeDenomMat3();
};

#endif // MAKEAMPMAT_H
