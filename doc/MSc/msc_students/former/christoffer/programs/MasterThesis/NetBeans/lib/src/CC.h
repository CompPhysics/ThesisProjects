/* 
 * File:   CC.h
 * Author: chrishir
 *
 * Created on 6. desember 2011, 14:46
 */

#ifndef CC_H
#define	CC_H

#include <armadillo>
#include "System.h"
#include "GEMM.h"
#include <iostream>

namespace toffyrn
{
    namespace libUNK
    {

        class CC
        {
        public:
            CC(double precision = 0.001, int max_iter = 20);
            virtual ~CC();
            void set_system(System * sys);
            double solve_ground_state_energy(int debug_flags = DEBUG_NONE, std::ostream &debug_stream = std::cout);

            arma::field<arma::vec> get_timing_info() const
            {
                return timing_info;
            }

            GEMM const * get_mult() const
            {
                return mult;
            }

            void set_mult(GEMM * mult)
            {
                if (this->mult != NULL)
                    delete(this->mult);
                this->mult = mult;
            }

            void reset_timing_info();


            ////FLAGS FOR DEBUGGING
            static const int DEBUG_NONE = 0;
            static const int DEBUG_T1 = 1 << 1;
            static const int DEBUG_T2 = 1 << 2;
            static const int DEBUG_T_ALL = DEBUG_T1 | DEBUG_T2;
            static const int DEBUG_I1 = 1 << 3;
            static const int DEBUG_I2 = 1 << 4;
            static const int DEBUG_I3 = 1 << 5;
            static const int DEBUG_I4 = 1 << 6;
            static const int DEBUG_I5 = 1 << 7;
            static const int DEBUG_I6 = 1 << 8;
            static const int DEBUG_I7 = 1 << 9;
            static const int DEBUG_I8 = 1 << 10;
            static const int DEBUG_I9 = 1 << 11;
            static const int DEBUG_I10 = 1 << 12;
            static const int DEBUG_I11 = 1 << 13;
            static const int DEBUG_I_T1 = DEBUG_I1 | DEBUG_I2 | DEBUG_I3 | DEBUG_I4 | DEBUG_I5;
            static const int DEBUG_I_T2 = DEBUG_I6 | DEBUG_I7 | DEBUG_I8 | DEBUG_I9 | DEBUG_I10 | DEBUG_I11;
            static const int DEBUG_I_ALL = DEBUG_I_T1 | DEBUG_I_T2;
            static const int DEBUG_ENERGY = 1 << 14;
            static const int DEBUG_INTERACTIONS = 1 << 15;
            static const int DEBUG_D1 = 1 << 16;
            static const int DEBUG_D2 = 1 << 17;
            static const int DEBUG_ALL = DEBUG_T_ALL | DEBUG_I_ALL | DEBUG_ENERGY | DEBUG_INTERACTIONS | DEBUG_D1 | DEBUG_D2;
            //END OF DEBUG FLAGS



        private:
            double precision;
            int max_iter;
            System * sys;
            GEMM * mult;


            void init_additional_mappings();
            std::vector <arma::uvec> map_p_p1;
            arma::umat map_p_p1_inv;
            std::vector <arma::uvec> map_p_h1;
            arma::umat map_p_h1_inv;
            std::vector <arma::uvec> map_p1_h;
            arma::umat map_p1_h_inv;
            std::vector <arma::uvec> map_h;
            arma::umat map_h_inv;
            std::vector <arma::uvec> map_p;
            arma::umat map_p_inv;
            std::vector <arma::uvec> map_p_p_p1;
            arma::umat map_p_p_p1_inv;
            std::vector <arma::uvec> map_p_p_h1;
            arma::umat map_p_p_h1_inv;
            std::vector <arma::uvec> map_h_h_p1;
            arma::umat map_h_h_p1_inv;
            std::vector <arma::mat> v_d1l_em1;
            std::vector <arma::mat> v_abd1_e;
            std::vector <arma::mat> v_ad1_el1;
            std::vector <arma::mat> v_d_lme1;



            /**
             * Containing timing information in seconds.
             * 0: Time used calculating energy
             * 1-11: Time used calculating i1-i11
             * 12: Time used calculating T1
             * 13: Time used calculating T2
             * 14: Time used calculating D1
             * 15: Time used calculating D2
             **/
            arma::field<arma::vec> timing_info;


            /**
             * @brief Calculate intermediate 1.
             * 
             * i1_ad = f_ad + <de||al> t^e_l
             * 
             * @param t1_old T_1 amplitudes from last iteration.
             * @return Matrix i1
             */
            arma::mat c_i1(arma::mat const &t1_old);
            void c_i1_t2(arma::mat &i1, arma::mat const &t1_old);

            /**
             * @brief Calculate intermediate 2
             * 
             * i2_dl = f_dl + <de||lm> t^e_m
             * 
             * @param t1_old T_1 amplitudes from last iteration.
             * @return Matrix i2
             */
            arma::mat c_i2(arma::mat const &t1_old);

            /**
             * @brief Calculate intermediate 3
             * 
             * i3_li = f_li + <di||ml> t^d_m + 0.5 <de||lm> t^de_im + i2_dl t^d_i
             * 
             * @param i2 The second intermediate.
             * @param t1_old T_1 amplitudes from last iteration.
             * @param t2_old T_2 amplitudes from last iteration.
             * @return Matrix i3
             */
            arma::mat c_i3(
                    arma::mat const &i2,
                    arma::mat const &t1_old,
                    std::vector<arma::mat> const &t2_old);
            void c_i3_t3(arma::mat &i3, std::vector<arma::mat> const &t2_old);

            /**
             * @brief Calculate intermediate 4
             * 
             * i4^di_lm = i5^di_lm + 0.5 <ed||lm> t^e_i
             * 
             * It is stored as i4(lmd)^nu_mu.
             * 
             * @param i5 The fifth intermediate.
             * @param t1_old T_1 amplitudes from last iteration.
             * @return STL vector with one matrix per channel.
             */
            std::vector<arma::mat> c_i4(
                    std::vector<arma::mat> const &i5,
                    arma::mat const &t1_old);


            /**
             * @brief Calculate intermediate 5
             * 
             * i5^di_lm = - <di||lm> + 0.5 <ed||lm> t^e_i
             * 
             * It is stored as i5(lmd)^nu_mu
             * 
             * @param t1_old T_1 amplitudes from last iteration.
             * @return STL vector with one matrix per channel.
             */
            std::vector<arma::mat> c_i5(arma::mat const &t1_old);
            void c_i5_t2(std::vector<arma::mat> &i5, arma::mat const &t1_old);

            /**
             * @brief Calculate T_1 denominator
             * 
             * D_ai = - i1_aa + i3_ii
             * 
             * @param i1 The first intermediate
             * @param i3 The third intermediate
             * @return Matrix D1.
             */
            arma::mat c_d1(arma::mat const &i1, arma::mat const &i3);

            /**
             * @brief Calculate intermediate 6
             * 
             * i6^lm_ij = <lm||ij> + 0.5 <de||lm>t^de_ij + i5^di_lm t^d_j 
             *                                           - i5^dj_lm t^d_i
             * 
             * @param i5 The fifth intermediate
             * @param t1_old T_1 amplitudes from last iteration.
             * @param t2_old T_2 amplitudes from last iteration.
             * @return STL vector with one matrix per channel.
             */
            std::vector<arma::mat> c_i6(
                    std::vector<arma::mat> const &i5,
                    arma::mat const &t1_old,
                    std::vector<arma::mat> const &t2_old);
            void c_i6_t3(
                    std::vector<arma::mat> &i6,
                    std::vector<arma::mat> const &i5,
                    arma::mat const &t1_old);

            /**
             * @brief Calculate intermediate 7
             * 
             * i7_bd = i1_bd - i2_dl t^b_l + 0.5 <de||lm> t^eb_lm
             * 
             * @param i1 The first intermediate
             * @param i2 The second intermediate
             * @param t1_old T_1 amplitudes from last iteration
             * @param t2_old T_2 amplitudes from last iteration
             * @return Matrix i7.
             */
            arma::mat c_i7(
                    arma::mat const &i1,
                    arma::mat const &i2,
                    arma::mat const &t1_old,
                    std::vector<arma::mat> const &t2_old);
            void c_i7_t3(arma::mat &i7, std::vector<arma::mat> const &t2_old);


            /**
             * 
             * i9^bl_dj = - <bl||dj> + 0.5 <ed||bl> t^e_j
             * 
             * @param t2_old
             * @return 
             */
            std::vector<arma::mat> c_i9(arma::mat const &t2_old);
            void c_i9_t2(std::vector<arma::mat> &i9, arma::mat const &t1_old);

            /**
             * 
             * i8^bl_dj = i9^bl_dj + 0.5 <ed||bl> t^e_j + i4^dj_lm t^b_m
             *           + 0.5 <de||lm> t^eb_mj
             * 
             * @param i4
             * @param i9
             * @param t1_old
             * @param t2_old
             * @return 
             */
            std::vector<arma::mat> c_i8(
                    std::vector<arma::mat> const &i4,
                    std::vector<arma::mat> const &i9,
                    arma::mat const &t1_old,
                    std::vector<arma::mat> const &t2_old);
            void c_i8_t2(std::vector<arma::mat> &i8, arma::mat const &t1_old);
            void c_i8_t3(
                    std::vector<arma::mat> &i8,
                    std::vector<arma::mat> const &i4,
                    arma::mat const &t1_old);

            /**
             * @brief Calculate the fourth term of i8.
             * 
             * i8^bl_dj += <de||lm> t^eb_mj
             * 
             * @param i8 Intermediate i8.
             * @param t2_old Old T2 amplitude.
             */
            void c_i8_t4(
                    std::vector<arma::mat> &i8,
                    std::vector<arma::mat> const &t2_old);

            /**
             * 
             * i10^al_ij = <al||ij> + 0.5 <de||al> t^de_ij
             *          + i9^al_di t^d_j - i9^al_dj t^d_i + 0.5 i6^lm_ij t^a_m
             * 
             * @param i6
             * @param i9
             * @param t1_old
             * @param t2_old
             * @return 
             */
            std::vector<arma::mat> c_i10(
                    std::vector<arma::mat> const &i6,
                    std::vector<arma::mat> const &i9,
                    arma::mat const &t1_old,
                    std::vector<arma::mat> const &t2_old);

            /**
             * 
             * i11^ab_dj =  <ab||dj> + 0.5 <ab||de> t^e_j
             * 
             * @param t1_old
             * @return 
             */
            std::vector<arma::mat> c_i11(arma::mat const &t1_old);
            void c_i11_t2(std::vector<arma::mat> &i11, arma::mat const &t1_old);

            /**
             * 
             * d2^ab_ij = -0.5<ab||ab> + i3_ii + i3_jj 
             *            -0.5 i6^ij_ij - i7_bb - i7_aa
             * 
             * @param i3
             * @param i6
             * @param i7
             * @return 
             */
            std::vector<arma::mat> c_d2(
                    arma::mat const &i3,
                    std::vector<arma::mat> const &i6,
                    arma::mat const &i7);

            /**
             * 
             * t1^a_i D_ai = f_ai + <la||di> t^d_l + 0.5 <de||al> t^de_il
             *              + i1_ad t^d_i - i3_li t^a_l + 0.5 i4^di_lm t^da_lm
             *              + i2_dl t^ad_il + t1^a_i D_ai
             * 
             * @param i1
             * @param i2
             * @param i3
             * @param i4
             * @param t1_old
             * @param t2_old
             * @return 
             */
            arma::mat c_t1(
                    arma::mat const &i1,
                    arma::mat const &i2,
                    arma::mat const &i3,
                    std::vector<arma::mat> const &i4,
                    arma::mat const &d1,
                    arma::mat const &t1_old,
                    std::vector<arma::mat> const &t2_old);

            void c_t1_t2(arma::mat &t1_new, arma::mat const &t1_old);
            void c_t1_t3(arma::mat &t1_new, std::vector<arma::mat> const &t2_old);
            void c_t1_t6(
                    arma::mat &t1_new,
                    std::vector<arma::mat> const &i4,
                    std::vector<arma::mat> const &t2_old);
            void c_t1_t7(
                    arma::mat &t1_new,
                    arma::mat const &i2,
                    std::vector<arma::mat> const &t2_old);


            /**
             * 
             * t^ab_ij D^ab_ij = <ab||ij> + 0.5 <ab||de> t^de_ij
             *                  - i3_li t^ab_lj - i3_lj t^ab_il
             *                  + 0.5 i6^lm_ij t^ab_lm
             *                  + i7_bd t^ad_ij + i7_ad t^db_ij
             *                  + P_ij P_ab i8^bl_dj t^ad_il - P_ab i10^al_ij t^b_l
             *                  + P_ij i11^ab_dj t^d_i + t^ab_ij D^ab_ij
             * 
             * @param i3
             * @param i6
             * @param i7
             * @param i8
             * @param i10
             * @param i11
             * @param d2
             * @param t1_old
             * @param t2_old
             * @return 
             */
            std::vector<arma::mat> c_t2(
                    arma::mat const &i3,
                    std::vector<arma::mat> const &i6,
                    arma::mat const &i7,
                    std::vector<arma::mat> const &i8,
                    std::vector<arma::mat> const &i10,
                    std::vector<arma::mat> const &i11,
                    std::vector<arma::mat> const &d2,
                    arma::mat const &t1_old,
                    std::vector<arma::mat> const &t2_old);

            void c_t2_t3(
                    std::vector<arma::mat> &t2_new,
                    arma::mat const &i3,
                    std::vector<arma::mat> const &t2_old);
            void c_t2_t5(
                    std::vector<arma::mat> &t2_new,
                    arma::mat const &i7,
                    std::vector<arma::mat> const &t2_old);
            void c_t2_t6(
                    std::vector<arma::mat> &t2_new,
                    std::vector<arma::mat> const &i8,
                    std::vector<arma::mat> const &t2_old);
            void c_t2_t7(
                    std::vector<arma::mat> &t2_new,
                    std::vector<arma::mat> const &i10,
                    arma::mat const &t1_old);
            void c_t2_t8(
                    std::vector<arma::mat> &t2_new,
                    std::vector<arma::mat> const &i11,
                    arma::mat const &t1_old);

        };

    }
}

#endif	/* CC_H */

