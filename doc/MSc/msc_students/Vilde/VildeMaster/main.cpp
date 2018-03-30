#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <random>
#include <vector>
#include <cmath>

using namespace std;
using namespace Eigen;

double interaction(VectorXd x, int nx);
void SGD(VectorXd &b, VectorXd &c, MatrixXd &w, VectorXd grad, double eta, int nx, int nh, ofstream &outfile);
void ASGD(VectorXd &b, VectorXd &c, MatrixXd &w, VectorXd grad, int cycles, int nx, int nh,
          double &asgd_X_prev, VectorXd &grad_prev, double &t_prev, double t, ofstream &outfile);

int main()
{
    // INITIALIZATIONS
    // Comment/uncomment line 169 to include/exclude interaction
    int nx = 4; // Number which represents particles*dimensions.
    int nh = 2; // Number of hidden units.
    int n_par = nx + nh + nx*nh;
    int n_cycles = 15000;  // 1000
    int n_samples = 10000;  // 100
    double sigma = 1.0; // Normal distribution visibles
    double omega = 1.0;
    double eta = 0.1; // SGD learning rate
    double x_mean;  // Normal distribution visibles
    double der1lnPsi;
    double der2lnPsi;
    // ASGD parameters
    double A = 20.0;
    double t_prev = A;
    double t = A;
    double asgd_X_prev;

    ofstream outfile;
    outfile.open("alphasNX" + to_string(nx) + "NH" + to_string(nh) + "CorrNO.txt");
    ofstream outfile2;
    outfile2.open("posAndHiddenNX" + to_string(nx) + "NH" + to_string(nh) + "CorrNO.txt");

    // DECLARATIONS
    // The Boltzmann machine / wave function
    VectorXd x;
    VectorXd h;
    VectorXd b;
    VectorXd c;
    MatrixXd w;
    x.resize(nx);
    h.resize(nh);
    b.resize(nx);
    c.resize(nh);
    w.resize(nx, nh);
    // Wf derived wrt rbm parameters, to be added up for each sampling
    VectorXd derPsi;
    VectorXd derPsi_temp;
    derPsi.resize(n_par);
    derPsi_temp.resize(n_par);
    // Local energy times wf derived wrt rbm parameters, to be added up for each sampling
    VectorXd EderPsi;
    EderPsi.resize(n_par);
    // Gradient, to be computed at each cycle, after all sampling iterations completed
    VectorXd grad;
    VectorXd grad_prev;
    grad.resize(n_par);
    grad_prev.resize(n_par);
    // Other
    VectorXd probHgivenX;
    VectorXd Q;
    probHgivenX.resize(nh);
    Q.resize(nh);

    // ASSIGNMENTS
    // The visible units/ positions
    std::uniform_real_distribution<double> distribution_initX(0,1);
    mt19937 rgen_Gibbs;
    std::random_device rd_Gibbs;
    rgen_Gibbs.seed(rd_Gibbs());
    for(int i=0; i<nx; i++){
        x(i)=distribution_initX(rgen_Gibbs);
    }
    // The rbm parameters
    int seed_initRBM=12345;
    float sigma_initRBM = 0.01;
    std::default_random_engine generator_initRBM(seed_initRBM);
    std::normal_distribution<double> distribution_initRBM(0,sigma_initRBM);
    for (int i=0; i<nx; i++){
        b(i) = distribution_initRBM(generator_initRBM);
        //outfile << b(i) << " ";
    }
    for (int i=0; i<nh; i++){
        c(i) = distribution_initRBM(generator_initRBM);
        //outfile << c(i) << " ";
    }
    for (int i=0; i<nx; i++){
        for (int j=0; j<nh; j++){
            w(i,j) = distribution_initRBM(generator_initRBM);
            //outfile << w(i,j) << " ";
        }
    }
    //outfile << '\n';

    random_device rd_h;
    default_random_engine generator_h(rd_h());
    std::uniform_real_distribution<double> distribution_setH(0,1);
    random_device rd_x;
    default_random_engine generator_x(rd_x());
    // START
    // Minimization loop
    for (int cycles=0; cycles<n_cycles; cycles++) {
        // Variables to store summations during sampling
        double Eloc = 0;
        double Eloc_temp = 0;
        double Eloc2 = 0;
        derPsi.setZero();
        EderPsi.setZero();


        // Sampling loop
        for (int samples=0; samples<n_samples; samples++) {
            // Compute prob of hidden given visible and set hidden values accordingly


            for (int j=0; j<nh; j++) {
                probHgivenX(j) = 1.0/(1 + exp(c(j) + x.transpose()*w.col(j)));
                h(j) = distribution_setH(generator_h) < probHgivenX(j);
                //if (cycles==59 && samples > 0.1*n_samples) {
                    //outfile2 << h(j) << " ";
                //}
                //cout << h(j) << endl;
            }
            // Set new positions (visibles) given hidden, according to normal distribution

            for (int i=0; i<nx; i++) {
                x_mean = b(i) + w.row(i)*h;
                //cout << b(i) << "  " << x_mean << endl;
                normal_distribution<double> distribution_x(x_mean, sigma);
                //cout << cycles << "   " << samples << "   " << i << "   " << x_mean << endl;
                x(i) = distribution_x(generator_x);
                //cout << cycles << "  " << samples << "  " << x(i) << "  " << x_mean << endl;
                //if (cycles==59 && samples > 0.1*n_samples) {
                    //outfile2 << x(i) << " ";
                //}
            }
            //if (cycles==59 && samples > 0.1*n_samples) {
            //    outfile2 << endl;
            //}
            if (samples > 0.1*n_samples) {
                // Compute the local energy
                Q = -c - (x.transpose()*w).transpose();
                Eloc_temp = 0;
                // Loop over the visibles (n_particles*n_coordinates) for the Laplacian
                for (int r=0; r<nx; r++) {
                    double sum1 = 0;
                    double sum2 = 0;
                    for (int j=0; j<nh; j++) {
                        sum1 += w(r,j)/(1.0+exp(-Q(j)));
                        sum2 += w(r,j)*w(r,j)*exp(Q(j))/((exp(Q(j))+1.0)*(exp(Q(j))+1.0));
                    }
                    der1lnPsi = -(x(r) - b(r)) - sum1;
                    der2lnPsi = -1.0 + sum2;
                    Eloc_temp += -der1lnPsi*der1lnPsi - der2lnPsi + omega*omega*x(r)*x(r);


                }
                Eloc_temp = 0.5*Eloc_temp;

                // With interaction:
                Eloc_temp += interaction(x, nx);


                Eloc += Eloc_temp;

                // Compute the 1/psi * dPsi/dalpha_i, that is Psi derived wrt each RBM parameter.
                for (int k=0; k<nx; k++) {
                    derPsi_temp(k) = (x(k) - b(k));
                }
                for (int k=nx; k<(nx+nh); k++) {
                    derPsi_temp(k) = -1.0/(1.0+exp(-Q(k-nx)));
                }
                int k=nx + nh;
                for (int i=0; i<nx; i++) {
                    for (int j=0; j<nh; j++) {
                        derPsi_temp(k) = x(i)/(1.0+exp(Q(j)));
                        k++;
                    }
                }
                // Add up values for expectation values
                derPsi += derPsi_temp;
                EderPsi += Eloc_temp*derPsi_temp;
                Eloc2 += Eloc_temp*Eloc_temp;
            }
        }
        // Compute expectation values
        double n_samp = n_samples - 0.1*n_samples;
        Eloc = Eloc/n_samp;
        Eloc2 = Eloc2/n_samp;
        derPsi = derPsi/n_samp;
        EderPsi = EderPsi/n_samp;

        double variance = Eloc2 - Eloc*Eloc;

        // Compute gradient
        grad = 2*(EderPsi - Eloc*derPsi);

        outfile << Eloc << " ";
        // Choose one of the methods:
        //SGD(b, c, w, grad, eta, nx, nh, outfile);
        ASGD(b, c, w, grad, cycles, nx, nh, asgd_X_prev, grad_prev, t_prev, t, outfile);

        double gradnorm = sqrt(grad.squaredNorm());
        cout << cycles << "   " << Eloc << "   " << variance << "   " << b(0) << //" "
             //<< b(1) << " "
             //<< c(0) << " " << c(1) << " " << c(2) << " " << c(3) << " "
             //<< w(0,0) << " " << w(0,1) << " " << w(0,2) << " " << w(0,3) << " "
             //<< w(1,0) << " " << w(1,1) << " " << w(1,2) << " " << w(1,3) <<
                endl;


    }

}






double interaction(VectorXd x, int nx) {
    double interaction_term = 0;
    double rx;
    double ry;
    double r_dist;
    for (int r=0; r<nx-2; r+=2) {
        for (int s=(r+2); s<nx; s+=2) {
            rx = (x(r) - x(s));
            ry = (x(r+1) - x(s+1));
            r_dist = sqrt(rx*rx + ry*ry);
            interaction_term += 1.0/r_dist;
        }

    }
    return interaction_term;
}





void SGD(VectorXd &b, VectorXd &c, MatrixXd &w, VectorXd grad, double eta, int nx, int nh, ofstream &outfile) {
    // Compute new parameters
    for (int i=0; i<nx; i++) {
        outfile << b(i) << " ";
        b(i) = b(i) - eta*grad(i);
    }
    for (int j=0; j<nh; j++) {
        outfile << c(j) << " ";
        c(j) = c(j) - eta*grad(nx + j);
    }
    int k = nx + nh;
    for (int i=0; i<nx; i++) {
        for (int j=0; j<nh; j++) {
            outfile << w(i,j) << " ";
            w(i,j) = w(i,j) - eta*grad(k);
            k++;
        }
    }
    outfile << '\n';

}

void ASGD(VectorXd &b, VectorXd &c, MatrixXd &w, VectorXd grad, int cycles, int nx, int nh,
          double &asgd_X_prev, VectorXd &grad_prev, double &t_prev, double t, ofstream &outfile) {
    // ASGD parameters
    double a = 0.1;
    double A = 20.0;
    double f_min = -0.5;
    double f_max = 2.0;
    double asgd_omega = 1.0;

    double f = f_min + (f_max - f_min)/(1 - (f_max/f_min)*exp(-asgd_X_prev/asgd_omega));
    t = 0;
    if (t < (t_prev + f)) t=t_prev+f;
    if (cycles==0 || cycles==1) t=A;
    double gamma = a/(t+A);

    //cout << f << " " << t << " " << asgd_X_prev << gamma << endl;
    // Compute new parameters
    for (int i=0; i<nx; i++) {
        outfile << b(i) << " ";
        b(i) = b(i) - gamma*grad(i);
    }
    for (int j=0; j<nh; j++) {
        outfile << c(j) << " ";
        c(j) = c(j) - gamma*grad(nx + j);
    }
    int k = nx + nh;
    for (int i=0; i<nx; i++) {
        for (int j=0; j<nh; j++) {
            outfile << w(i,j) << " ";
            w(i,j) = w(i,j) - gamma*grad(k);
            k++;
        }
    }
    outfile << '\n';
    asgd_X_prev = -grad.dot(grad_prev);
    grad_prev = grad;
    t_prev = t;


}

//void ASGD_test() {
    // Flipped Gaussian distribution f(x|mu, sigma**2) = -1/sqrt(2*pi*sigma**2)*exp(-((x-mu)**2)/(2*sigma**2))
    // Optimize with respect to x, which should go to zero.
//    double mu = 0.0;
//    double sigma = 1.0;

//}

