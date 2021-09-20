/**
 * @file   CLgemm.cpp
 * @author toffyrn
 * 
 * Created on 15. februar 2012, 12:49
 */

#include "CLgemm.h"

using namespace toffyrn::libUNK;
using namespace std;
using namespace arma;

CLgemm::CLgemm()
{
    context = createSomeContext();
    //Get selected device
    vector<cl::Device> device_vec = context.getInfo<CL_CONTEXT_DEVICES > ();
    cl::Device device = device_vec[0];
    //Create a CommandQueue
    queue = cl::CommandQueue(context, device);
    //Print som rather coooooool info!
    cout << "Initialised device " << device.getInfo<CL_DEVICE_NAME > ();
    int comp_units = device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS > ();
    cout << ".\nThis device has " << comp_units <<
            " compute units.\n";

    clAmdBlasSetup();
    tot_time = 0;
}

CLgemm::CLgemm(const CLgemm& orig)
{
}

CLgemm::~CLgemm()
{
    clAmdBlasTeardown();
}

void CLgemm::dgemm(
        arma::mat &res,
        arma::mat const &left,
        arma::mat const &right,
        double alpha,
        double beta,
        bool transL,
        bool transR)
{
    timer.tic();
    
    if (transL == false && transR == false)
    {
        //Check if res is correctly sized
        if (res.n_rows != left.n_rows || res.n_cols != right.n_cols)
            throw std::string("Badly shaped \"res\" in CLgemm::dgemm");

        int m = left.n_rows;
        int n = right.n_cols;
        int p = left.n_cols;
        int work = m * n * p;

        if (work < 35000000) //TODO: A more carefull analysis of this decision should be done.
        {
            res = alpha * left * right + beta * res;
        } else
        {
            //The memptr can not be const in cl, however it is expected not to change.
            cl_double *A_p = const_cast<double *> (left.memptr());
            cl_double *B_p = const_cast<double *> (right.memptr());
            cl_double *C_p = res.memptr();
            //Create CL buffers, pointing to Host memory.
            cl::Buffer A_cl(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof (*A_p) * left.n_elem, A_p);
            cl::Buffer B_cl(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof (*B_p) * right.n_elem, B_p);
            cl::Buffer C_cl(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof (*C_p) * res.n_elem, C_p);

            const size_t amd_M = left.n_rows;
            const size_t amd_N = right.n_cols;
            const size_t amd_P = right.n_rows;
            if (left.n_cols != right.n_rows)
                throw std::string("CLgemm: left and right matrix dimensions not compatible.");
            const clAmdBlasOrder amd_order = clAmdBlasColumnMajor;
            const clAmdBlasTranspose amd_transAll = clAmdBlasNoTrans;
            const cl_double amd_alpha = alpha;
            const cl_double amd_beta = beta;
            const size_t amd_lda = amd_M;
            const size_t amd_ldb = amd_P;
            const size_t amd_ldc = amd_M;

            cl::Event e;
            clAmdBlasDgemm(amd_order, amd_transAll, amd_transAll, amd_M, amd_N, amd_P,
                    amd_alpha, A_cl(), amd_lda, B_cl(), amd_ldb, amd_beta, C_cl(), amd_ldc,
                    1, &queue(), 0, NULL, &e());
            queue.enqueueReadBuffer(C_cl, true, 0, sizeof (*C_p) * res.n_elem, C_p); //BLOCKING

            //clReleaseMemObject();
        }
    } else if (transL == true && transR == false)
        res = alpha * trans(left) * right + beta * res;
    else if (transL == false && transR == true)
        res = alpha * left * trans(right) + beta * res;
    else if (transL == true && transR == true)
        res = alpha * trans(left) * trans(right) + beta * res;
    tot_time += timer.toc();

    return;
}

cl::Context CLgemm::createSomeContext()
{
    //Query CL for all available platforms.
    //A platform is supplied by a vendor, 
    //and can contain multiple devices.
    std::vector<cl::Platform> all_platforms;
    cl::Platform::get(&all_platforms);

    //Structure to store all devices and a corresponding ID
    int idx = 0;
    std::vector<cl::Device> all_devices;

    //Loop over all platforms.
    for (int i = 0; i < all_platforms.size(); i++)
    {
        //Print platform information
        cout << all_platforms[i].getInfo<CL_PLATFORM_NAME > () << "\t"
                << all_platforms[i].getInfo<CL_PLATFORM_VERSION > () << endl;

        //Retreive current platforms devices
        std::vector<cl::Device> devices;
        all_platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &devices);

        //Loop over all devices, and print their name
        for (int j = 0; j < devices.size(); j++)
        {
            cout << "  [[" << idx << "]]\t" << devices[j].getInfo<CL_DEVICE_NAME > () << endl;
            all_devices.push_back(devices[j]);
            idx++;
        }
    }

    //Let user choose device.
    cout << "Select device: ";
    cin >> idx;

    //Create a context with one selected device.
    std::vector<cl::Device> chosen_dev;
    chosen_dev.push_back(all_devices[idx]);
    cl::Context context(chosen_dev);

    return context;
}

