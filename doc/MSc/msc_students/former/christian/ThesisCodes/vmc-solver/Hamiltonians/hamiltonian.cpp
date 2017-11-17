#include "hamiltonian.h"
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

Hamiltonian::Hamiltonian(System* system, bool analyticalKinetic, double alpha, double omega) {
    m_system = system;
    m_analyticalKinetic = analyticalKinetic;
    m_alpha = alpha;
    m_omega = omega;
    setUpHermitePolynomials();
}

double Hamiltonian::computeKineticEnergy(std::vector<Particle*> particles){
    // Compute the kinetic energy using numerical differentiation.

    double numberOfParticles = m_system->getNumberOfParticles();
    double numberOfDimensions = m_system->getNumberOfDimensions();
    double h = 1e-4;

    // Evaluate wave function at current step
    double waveFunctionCurrent = m_system->getWaveFunction()->evaluate(particles);
    double kineticEnergy = 0;

    for (int i=0; i < numberOfParticles; i++){
        for (int j=0; j < numberOfDimensions; j++){

            // Evaluate wave function at forward step
            particles[i]->adjustPosition(h, j);
            m_system->getWaveFunction()->updateDistances(i);
            m_system->getWaveFunction()->updateSPWFMat(i);
            m_system->getWaveFunction()->updateJastrow(i);
            double waveFunctionPlus = m_system->getWaveFunction()->evaluate(particles);

            // Evaluate wave function at backward step
            particles[i]->adjustPosition(-2*h, j);
            m_system->getWaveFunction()->updateDistances(i);
            m_system->getWaveFunction()->updateSPWFMat(i);
            m_system->getWaveFunction()->updateJastrow(i);
            double waveFunctionMinus = m_system->getWaveFunction()->evaluate(particles);

            // Part of numerical diff
            kineticEnergy -= (waveFunctionPlus - 2*waveFunctionCurrent + waveFunctionMinus);

            // Move particles back to original position
            particles[i]->adjustPosition(h, j);
            m_system->getWaveFunction()->updateDistances(i);
            m_system->getWaveFunction()->updateSPWFMat(i);
            m_system->getWaveFunction()->updateJastrow(i);
        }
    }
    // Other part of numerical diff. Also divide by evaluation of current wave function
    // and multiply by 0.5 to get the actual kinetic energy.
    kineticEnergy = 0.5*kineticEnergy / (waveFunctionCurrent*h*h);
    return kineticEnergy;
}

void Hamiltonian::setUpHermitePolynomials() {
    int maxImplemented = 50;

    m_hermitePolynomials = new HermitePolynomials*[maxImplemented];
    m_hermitePolynomialsDerivative = new HermitePolynomials*[maxImplemented];
    m_hermitePolynomialsDoubleDerivative = new HermitePolynomials*[maxImplemented];

    m_hermitePolynomials[0] = new HermitePolynomial_0(m_alpha, m_omega);
    m_hermitePolynomials[1] = new HermitePolynomial_1(m_alpha, m_omega);
    m_hermitePolynomials[2] = new HermitePolynomial_2(m_alpha, m_omega);
    m_hermitePolynomials[3] = new HermitePolynomial_3(m_alpha, m_omega);
    m_hermitePolynomials[4] = new HermitePolynomial_4(m_alpha, m_omega);
    m_hermitePolynomials[5] = new HermitePolynomial_5(m_alpha, m_omega);
    m_hermitePolynomials[6] = new HermitePolynomial_6(m_alpha, m_omega);
    m_hermitePolynomials[7] = new HermitePolynomial_7(m_alpha, m_omega);
    m_hermitePolynomials[8] = new HermitePolynomial_8(m_alpha, m_omega);
    m_hermitePolynomials[9] = new HermitePolynomial_9(m_alpha, m_omega);
    m_hermitePolynomials[10] = new HermitePolynomial_10(m_alpha, m_omega);
    m_hermitePolynomials[11] = new HermitePolynomial_11(m_alpha, m_omega);
    m_hermitePolynomials[12] = new HermitePolynomial_12(m_alpha, m_omega);
    m_hermitePolynomials[13] = new HermitePolynomial_13(m_alpha, m_omega);
    m_hermitePolynomials[14] = new HermitePolynomial_14(m_alpha, m_omega);
    m_hermitePolynomials[15] = new HermitePolynomial_15(m_alpha, m_omega);
    m_hermitePolynomials[16] = new HermitePolynomial_16(m_alpha, m_omega);
    m_hermitePolynomials[17] = new HermitePolynomial_17(m_alpha, m_omega);
    m_hermitePolynomials[18] = new HermitePolynomial_18(m_alpha, m_omega);
    m_hermitePolynomials[19] = new HermitePolynomial_19(m_alpha, m_omega);
    m_hermitePolynomials[20] = new HermitePolynomial_20(m_alpha, m_omega);
    m_hermitePolynomials[21] = new HermitePolynomial_21(m_alpha, m_omega);
    m_hermitePolynomials[22] = new HermitePolynomial_22(m_alpha, m_omega);
    m_hermitePolynomials[23] = new HermitePolynomial_23(m_alpha, m_omega);
    m_hermitePolynomials[24] = new HermitePolynomial_24(m_alpha, m_omega);
    m_hermitePolynomials[25] = new HermitePolynomial_25(m_alpha, m_omega);
    m_hermitePolynomials[26] = new HermitePolynomial_26(m_alpha, m_omega);
    m_hermitePolynomials[27] = new HermitePolynomial_27(m_alpha, m_omega);
    m_hermitePolynomials[28] = new HermitePolynomial_28(m_alpha, m_omega);
    m_hermitePolynomials[29] = new HermitePolynomial_29(m_alpha, m_omega);
    m_hermitePolynomials[30] = new HermitePolynomial_30(m_alpha, m_omega);
    m_hermitePolynomials[31] = new HermitePolynomial_31(m_alpha, m_omega);
    m_hermitePolynomials[32] = new HermitePolynomial_32(m_alpha, m_omega);
    m_hermitePolynomials[33] = new HermitePolynomial_33(m_alpha, m_omega);
    m_hermitePolynomials[34] = new HermitePolynomial_34(m_alpha, m_omega);
    m_hermitePolynomials[35] = new HermitePolynomial_35(m_alpha, m_omega);
    m_hermitePolynomials[36] = new HermitePolynomial_36(m_alpha, m_omega);
    m_hermitePolynomials[37] = new HermitePolynomial_37(m_alpha, m_omega);
    m_hermitePolynomials[38] = new HermitePolynomial_38(m_alpha, m_omega);
    m_hermitePolynomials[39] = new HermitePolynomial_39(m_alpha, m_omega);
    m_hermitePolynomials[40] = new HermitePolynomial_40(m_alpha, m_omega);
    m_hermitePolynomials[41] = new HermitePolynomial_41(m_alpha, m_omega);
    m_hermitePolynomials[42] = new HermitePolynomial_42(m_alpha, m_omega);
    m_hermitePolynomials[43] = new HermitePolynomial_43(m_alpha, m_omega);
    m_hermitePolynomials[44] = new HermitePolynomial_44(m_alpha, m_omega);
    m_hermitePolynomials[45] = new HermitePolynomial_45(m_alpha, m_omega);
    m_hermitePolynomials[46] = new HermitePolynomial_46(m_alpha, m_omega);
    m_hermitePolynomials[47] = new HermitePolynomial_47(m_alpha, m_omega);
    m_hermitePolynomials[48] = new HermitePolynomial_48(m_alpha, m_omega);
    m_hermitePolynomials[49] = new HermitePolynomial_49(m_alpha, m_omega);

    m_hermitePolynomialsDerivative[0] = new dell_HermitePolynomial_0(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[1] = new dell_HermitePolynomial_1(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[2] = new dell_HermitePolynomial_2(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[3] = new dell_HermitePolynomial_3(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[4] = new dell_HermitePolynomial_4(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[5] = new dell_HermitePolynomial_5(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[6] = new dell_HermitePolynomial_6(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[7] = new dell_HermitePolynomial_7(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[8] = new dell_HermitePolynomial_8(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[9] = new dell_HermitePolynomial_9(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[10] = new dell_HermitePolynomial_10(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[11] = new dell_HermitePolynomial_11(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[12] = new dell_HermitePolynomial_12(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[13] = new dell_HermitePolynomial_13(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[14] = new dell_HermitePolynomial_14(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[15] = new dell_HermitePolynomial_15(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[16] = new dell_HermitePolynomial_16(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[17] = new dell_HermitePolynomial_17(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[18] = new dell_HermitePolynomial_18(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[19] = new dell_HermitePolynomial_19(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[20] = new dell_HermitePolynomial_20(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[21] = new dell_HermitePolynomial_21(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[22] = new dell_HermitePolynomial_22(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[23] = new dell_HermitePolynomial_23(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[24] = new dell_HermitePolynomial_24(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[25] = new dell_HermitePolynomial_25(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[26] = new dell_HermitePolynomial_26(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[27] = new dell_HermitePolynomial_27(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[28] = new dell_HermitePolynomial_28(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[29] = new dell_HermitePolynomial_29(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[30] = new dell_HermitePolynomial_30(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[31] = new dell_HermitePolynomial_31(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[32] = new dell_HermitePolynomial_32(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[33] = new dell_HermitePolynomial_33(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[34] = new dell_HermitePolynomial_34(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[35] = new dell_HermitePolynomial_35(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[36] = new dell_HermitePolynomial_36(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[37] = new dell_HermitePolynomial_37(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[38] = new dell_HermitePolynomial_38(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[39] = new dell_HermitePolynomial_39(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[40] = new dell_HermitePolynomial_40(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[41] = new dell_HermitePolynomial_41(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[42] = new dell_HermitePolynomial_42(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[43] = new dell_HermitePolynomial_43(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[44] = new dell_HermitePolynomial_44(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[45] = new dell_HermitePolynomial_45(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[46] = new dell_HermitePolynomial_46(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[47] = new dell_HermitePolynomial_47(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[48] = new dell_HermitePolynomial_48(m_alpha, m_omega);
    m_hermitePolynomialsDerivative[49] = new dell_HermitePolynomial_49(m_alpha, m_omega);

    m_hermitePolynomialsDoubleDerivative[0] = new lapl_HermitePolynomial_0(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[1] = new lapl_HermitePolynomial_1(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[2] = new lapl_HermitePolynomial_2(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[3] = new lapl_HermitePolynomial_3(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[4] = new lapl_HermitePolynomial_4(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[5] = new lapl_HermitePolynomial_5(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[6] = new lapl_HermitePolynomial_6(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[7] = new lapl_HermitePolynomial_7(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[8] = new lapl_HermitePolynomial_8(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[9] = new lapl_HermitePolynomial_9(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[10] = new lapl_HermitePolynomial_10(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[11] = new lapl_HermitePolynomial_11(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[12] = new lapl_HermitePolynomial_12(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[13] = new lapl_HermitePolynomial_13(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[14] = new lapl_HermitePolynomial_14(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[15] = new lapl_HermitePolynomial_15(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[16] = new lapl_HermitePolynomial_16(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[17] = new lapl_HermitePolynomial_17(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[18] = new lapl_HermitePolynomial_18(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[19] = new lapl_HermitePolynomial_19(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[20] = new lapl_HermitePolynomial_20(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[21] = new lapl_HermitePolynomial_21(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[22] = new lapl_HermitePolynomial_22(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[23] = new lapl_HermitePolynomial_23(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[24] = new lapl_HermitePolynomial_24(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[25] = new lapl_HermitePolynomial_25(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[26] = new lapl_HermitePolynomial_26(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[27] = new lapl_HermitePolynomial_27(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[28] = new lapl_HermitePolynomial_28(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[29] = new lapl_HermitePolynomial_29(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[30] = new lapl_HermitePolynomial_30(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[31] = new lapl_HermitePolynomial_31(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[32] = new lapl_HermitePolynomial_32(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[33] = new lapl_HermitePolynomial_33(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[34] = new lapl_HermitePolynomial_34(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[35] = new lapl_HermitePolynomial_35(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[36] = new lapl_HermitePolynomial_36(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[37] = new lapl_HermitePolynomial_37(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[38] = new lapl_HermitePolynomial_38(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[39] = new lapl_HermitePolynomial_39(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[40] = new lapl_HermitePolynomial_40(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[41] = new lapl_HermitePolynomial_41(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[42] = new lapl_HermitePolynomial_42(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[43] = new lapl_HermitePolynomial_43(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[44] = new lapl_HermitePolynomial_44(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[45] = new lapl_HermitePolynomial_45(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[46] = new lapl_HermitePolynomial_46(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[47] = new lapl_HermitePolynomial_47(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[48] = new lapl_HermitePolynomial_48(m_alpha, m_omega);
    m_hermitePolynomialsDoubleDerivative[49] = new lapl_HermitePolynomial_49(m_alpha, m_omega);

}

