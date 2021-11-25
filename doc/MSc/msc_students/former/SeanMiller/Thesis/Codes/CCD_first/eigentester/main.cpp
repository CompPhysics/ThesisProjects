#include <iostream>
#include <math.h>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <map>
#include <unordered_map>
#include <chrono>

#include <iomanip>
#include <cstdlib>
#include <string>
#include <fstream>
#include <complex>

//#include <unsupported/Eigen/CXX11/Tensor>

//#include <mpi.h>
#include <omp.h>
#include <sparsepp/spp.h>
#include <sparsepp/spp_utils.h>
//#include <stdio.h>
#include "testclass.h"

typedef std::chrono::high_resolution_clock Clock;   //needed for timing

using namespace std;

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

/*bool myfunction (int i, int j) {
  return (i==j);
}*/
bool contractor(int i, int j){ return i==j; }

int* returnInt(int i){
    return &i;
}

//using Matrix = Eigen::MatrixXd;

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



int main(int argc, char** argv)
{


    int threads = 4;
    std::vector<spp::sparse_hash_map<int, int>> VP;
    VP.resize(threads);
    spp::sparse_hash_map<int, int> VS;
    std::unordered_map<int, int> VSS;

    int lim = 1e7; int in;
    for(int i=0;i<lim;i++){
        in = i;
        VS[in] = i;
        VSS[in] = i;
    }
    cout << VS.size() << " " << VSS.size() << endl;
    int something=0;
    auto T1 = Clock::now(); auto T2=Clock::now();
    T1 = Clock::now();

    for(int i=0;i<lim;i++){
        something=VS[i];
    }
    T2=Clock::now();

    std::cout << "time: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(T2 - T1).count()
              << " nanoseconds."
              << std::endl;
    T1 = Clock::now();

    for(int i=0;i<lim;i++){
        something=VSS[i];
    }
    T2=Clock::now();

    std::cout << "time: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(T2 - T1).count()
              << " nanoseconds."
              << std::endl;



    for (int j=0;j<threads;j++){
        VP[j] = VS;
    }

    int x=10; int y=10;
    Eigen::MatrixXi MatR;
    MatR.conservativeResize(x,y);
    MatR.setRandom();

    Eigen::MatrixXi MatW1;
    MatW1.conservativeResize(x,y);

    Eigen::MatrixXi MatW2;
    MatW2.conservativeResize(x,y);

    auto t1 = Clock::now();

    int thread;
    //spp::sparse_hash_map<unsigned long int, int>* map = nullptr;
#pragma omp parallel for num_threads(4) private(thread)
    for (int c=0;c<MatW1.cols();c++){
        thread = omp_get_thread_num();
        //map = &VP[thread];
        for(int r=0;r<MatW1.rows();r++){
            MatW1(r,c) = VP[thread][MatR(r,c)];
            //MatW1(r,c) = (int)map[(unsigned long int)MatR(r,c)];
        }
    }

    auto t2 = Clock::now();

    std::cout << "time: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()
              << " nanoseconds."
              << std::endl;

    t1 = Clock::now();

    for (int c=0;c<MatW2.cols();c++){
        for(int r=0;r<MatW2.rows();r++){
            MatW2(r,c) = VS[MatR(r,c)];
        }
    }

    t2 = Clock::now();

    std::cout << "time: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()
              << " nanoseconds."
              << std::endl;

#pragma omp parallel for num_threads(4)
    for (int c=0;c<MatW1.cols();c++){
        for(int r=0;r<MatW1.rows();r++){
            MatW1(r,c) = VSS[MatR(r,c)];
        }
    }
    cout << "success" << endl;


    long int A= 10;
    int counter = 0;
    while (A>0 && counter < 20){
        A *= 10;
        counter += 1;
        cout << A << endl;
    }

    std::cout << "------------" << std::endl;
    int Nh = 14;
    int Ns = 30;
    int count1 = 0;
    int count2 = 0;
    for (int a=Nh; a<Ns; a++){
        for (int b=a+1; b<Ns; b++){
            for (int c=b+1; c<Ns; c++){
                for (int i=0; i<Nh; i++){
                    for (int j=i+1; j<Nh; j++){
                        count1 ++;
                    }
                }
            }
        }
    }
    for (int a=Nh; a<Ns; a++){
        for (int b=Nh; b<Ns; b++){
            for (int c=Nh; c<Ns; c++){
                for (int i=0; i<Nh; i++){
                    for (int j=0; j<Nh; j++){
                        count2 ++;
                    }
                }
            }
        }
    }
    std::cout << "New way: " << count1 << std::endl;
    std::cout << "Old way: " << count2 << std::endl;

    //std::complex<double> hey;

    /*int max = 1e6;
    int threads = 4;
    int max4 = max/threads;
    //std::vector<std::unordered_map<int, double>> VP;
    std::unordered_map<int, double> VP[threads];
    //VP.resize(threads);
    for (int j=0; j<threads; j++){
        for (int i = j*max4; i<max4*(j+1); i++){
            double hey = i;
            VP[j][i] = hey;
        }
    }

    spp::sparse_hash_map<int, int> VSS;

    std::unordered_map<int, double> VS;
    for (int i = 0; i<max; i++){
        double hey = i;
        VS[i] = hey;
        VSS[i] = hey;
    }
    //cout << VS.size() << endl;

    int find = -100;
    double finds = 0;
    double findsP[threads];

    int times[threads];
    for (int j=0; j<threads; j++){
        times[j] = 0;
    }

    int x=1000; int y=100;
    Eigen::MatrixXi MatR;
    MatR.conservativeResize(x,y);
    MatR.setRandom();

    Eigen::MatrixXi MatW1;
    MatW1.conservativeResize(x,y);

    Eigen::MatrixXi MatW2;
    MatW2.conservativeResize(x,y);

    Eigen::MatrixXi MatB;
    MatB.conservativeResize(x,y);


    auto T1 = Clock::now();
    int step = MatR.cols()/threads;*/
   /* for (int k=0; k<threads; k++){
#pragma omp parallel num_threads(threads)
        {
            int thread = omp_get_thread_num();
            int val; int move = 0;
            int minCol = step*thread+step*k;
            int maxCol = step*(thread+1)+step*k;
            if (maxCol > MatR.cols()){move=-MatR.cols();}
            for (int c=minCol+move; c<maxCol+move; c++){
                for (int r=0; r<MatR.rows(); r++){
                    if (MatB(r,c)!=1){
                        val = VP[thread][MatR(r,c)];
                        if (val != 0){MatW(r,c) = val; MatB(r,c) = 1;}
                    }
                }
            }
        }
    }*/
    /*#pragma omp parallel for ordered num_threads(threads)
    for (int c=0; c<MatR.cols(); c++){
        for (int r=0; r<MatR.rows(); r++){
            MatW1(r,c) = VS[MatR(r,c)];
        }
    }*/

    /*auto T2 = Clock::now();
    std::cout << "time: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(T2 - T1).count()
              << " nanoseconds."
              << std::endl;

    auto T3 = Clock::now();
    //#pragma omp parallel for num_threads(threads)
    for (int c=0; c<MatR.cols(); c++){
        for (int r=0; r<MatR.rows(); r++){
            //cout << "hey" << endl;
            MatW2(r,c) = VSS[MatR(r,c)];
        }
    }
    auto T4 = Clock::now();
    std::cout << "time: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(T4 - T3).count()
              << " nanoseconds."
              << std::endl;*/

   /* for (int c=0; c<MatR.cols(); c++){
        for (int r=0; r<MatR.rows(); r++){
            if (MatW1(r,c)!=MatW2(r,c)){cout << "FUCK" << endl;}
        }
    }*/

    //cout << "start" << endl;
    //auto t5 = Clock::now();
    /*#pragma omp parallel for num_threads(threads)
    for (int j=0; j<threads; j++){
        auto t1 = Clock::now();
        findsP[j] = VP[j][find];
        auto t2 = Clock::now();
        times[j] = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        //std::cout << VP[j][find] << std::endl;
    }*/
    /*#pragma omp parallel num_threads(threads)
        findsP[omp_get_thread_num()] = VP[omp_get_thread_num()][find];*/


    /*auto t6 = Clock::now();
    cout << "end" << endl;

    int time=0;
    for (int j=0; j<threads; j++){
        if (time<times[j]){time = times[j];}
        cout << times[j] << endl;
    }

    std::cout << "time: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(t6 - t5).count()
              << " nanoseconds. Val: "
              << finds
              << std::endl;

    auto t3 = Clock::now(); cout << "start" << endl;
    finds = VS[find];
    auto t4 = Clock::now(); cout << "end" << endl;
    std::cout << "time: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count()
              << " nanoseconds. Val: "
              << finds
              << std::endl;*/


    //std::cout << 14 % 6 << std::endl;

    /*ofstream myfile;
    ostringstream os;
    os << "verytemp.txt";
    string s = os.str();

    for (int i=0; i<20; i++){
        myfile.open(s, ios::app);
        myfile <<  i << "\n" << std::flush;
        myfile.close();
    }*/

    /*int min = 14;
    int max = 2000;
    int index = 0;

    std::vector<int> checks;

    for (int i=0;i<30;i++){
        checks.push_back( i*10 );
    }

    int size = max*max*max*max;
    Eigen::MatrixXi mat;
    mat.conservativeResize(1, size);
    std::cout << "sup" << std::endl;


    auto t1 = Clock::now();

    for (int i1 = min; i1<max; i1++){
        for (int i2 = min; i2<max; i2++){
            for (int i3 = min; i3<max; i3++){
                int check = i1*1+i2*2+i3*3;
                auto it = std::find(checks.begin(), checks.end(), check);
                if (it != checks.end()){
                    mat(0,index) = check;
                    index ++;
                }
            }
        }
    }


    auto t2 = Clock::now();

    std::cout << "end. time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
              << " milliseconds" << std::endl;*/

    /*for (int m = 8; m<9; m++){
        Eigen::initParallel();

        omp_set_num_threads(m);
        Eigen::setNbThreads(m);

        int n = 4000;

        auto t1 = Clock::now();
        //Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
        //Eigen::MatrixXd B = Eigen::MatrixXd::Random(n,n);
        //Eigen::MatrixXd C;

        std::cout << "start" << std::endl;

        //C = A*B;
        //Eigen::MatrixXd C = A*B;
        Eigen::MatrixXd C = Eigen::MatrixXd::Random(n,n)*Eigen::MatrixXd::Random(n,n);
        auto t2 = Clock::now();


        std::cout << "end. time: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
                  << " milliseconds" << std::endl;
    }*/


    /*std::unordered_map<int, double> T;
    std::vector<int> D1;
    std::vector<int> D2;
    std::vector<int> D3;

    int limit = 1e5;

    for (int i = 0; i<limit; i++){
        T[i] = limit-1-i;
        D1.push_back(i);
        D2.push_back(limit-1-i);
        D3.push_back(i+1);
    }

    testClass* Tester = new testClass;

    Tester->D1 = D1;
    Tester->D2 = D2;

    auto t1 = Clock::now();

    for (int j=0; j<limit; j++){
        //std::cout << D3[T[j]] << std::endl;
        std::cout << D3[Tester->testy( j )] << std::endl;
    }
    //std::cout << "wut" << std::endl;
    auto t2 = Clock::now();

    std::cout << "end. time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
              << " milliseconds" << std::endl;*/

    //MPI_Init (&argc, &argv);	/* starts MPI */
    //int rank, size;
    //MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
    //MPI_Comm_size (MPI_COMM_WORLD, &size);	/* get number of processes */

    //cout << size << endl;
    //MPI_Finalize();
    /*Eigen::initParallel();
    Eigen::setNbThreads(4);

    int nthreads = Eigen::nbThreads( );
    std::cout << "THREADS = " << nthreads << std::endl; // returns '1'

    Eigen::Matrix<int, 3, 10> m1;
    Eigen::Matrix<int, 3, 4> m2;

    Eigen::Vector4i v1(1,2,3,4);
    Eigen::Vector4i v2(1,2,3,4);

    //Eigen::VectorXi v = Eigen::VectorXi(1,2,3);

    Eigen::Matrix<int, 5, 1> vec;
    vec << 1,2,3,4,5;


    m1 << 1,2,3,3,4,5,2,3,1,4,
                4,5,6,2,4,6,1,3,6,3,
                7,8,9,4,5,2,6,3,2,4;*/


    /*int xLim = 2e3;
    int yLim = 2e3;
    Eigen::MatrixXi M;
    std::vector<int> V;
    map<int, int> Map;
    M.conservativeResize(xLim,yLim);

    //#pragma omp parallel for
    for (int x=0; x<xLim; x++){
        for (int y=0; y<yLim; y++){
            //V.push_back( x+y );
            Map[x+y] = x+y;
        }
    }

    int a = 2147483647;
    std::cout << a << std::endl;

    auto t1 = Clock::now();

    //int omp_get_num_threads( );
    //cout << omp_get_num_threads() << endl;
    #pragma omp parallel for
    for (int x=0; x<xLim; x++){
        for (int y=0; y<yLim; y++){
            M(x,y) = Map[y*x];
            M(x,y) += Map[x+y];
        }
    }

    auto t2 = Clock::now();

    std::cout << "Total time used: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()
              << " " << std::endl;*/

    //MPI_Finalize();

    /*for (int i=0; i<1e3; i++){
        Vec.push_back(i+1);
    }
    int val = 9;
    int index;
    auto it = std::find(Vec.begin(), Vec.end(), val);
    if (it == Vec.end()){
      cout << "value not in map" << endl;// name not in vector
    }
    else{
      index = distance(Vec.begin(), it);
    }

    cout << Vec[index] << endl;

    Vec[index] = 12;

    cout << Vec[index] << endl;

*/

    /*map<int, double> places1;
    map<int, double> places2;

    for (int i=0; i<1e7; i++){
        places1[i] = (double) i;
        places2[i] = (double) 1e3-i;
    }

    cout << "done" << endl;
    cout << places1[1e6] << endl;*/
/*
    places1 = places2;
    //places1[9999] += 2;
    cout << places1[9999] << endl;


    std::vector<int> v3_temp;
    for (int i=0; i<1e3; i++){
        v3_temp.push_back(i);
    }
    cout << v3_temp[0] << endl;

    //VectorXf::Map(pointer, size)
      //  double* ptr = &v[0];
    //Eigen::VectorXi v3;

    int* ptr = &v3_temp[0];
    Eigen::Map<Eigen::VectorXi> v3(ptr, v3_temp.size());
    //cout << v3 << endl;


    std::vector<int> array2 = { 9, 7, 5, 3, 1 };*/



    /*Eigen::VectorXi vec1;
    vec1 << 1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,5,5,5,5,5,5,5,5,6,6,6,6;
    Eigen::VectorXi vec2;
    vec2 << -1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,1,1,1,2,2,2,2,3,3,3,3,3,4,4,4,5,5,5,5,5,5,5,5,6,6,6,6,7,7,7;*/

    //cout << vecDelta(v1+v1,v2+v2) << endl;

    /*for (int i= 0; i<1; i++){
        cout << i << endl;
    }*/

    /*m1 << 1,2,3,3,4,5,2,3,1,4,
            4,5,6,2,4,6,1,3,6,3,
            7,8,9,4,5,2,6,3,2,4;

    std::vector<int> v = {1,2,3,4,5,6,7,8};

    //Eigen::Vector4f v(1,2,3,4);
    v.erase( unique( v.begin(), v.end() ), v.end() );

    //cout << v.size() << endl;
    for (int i; i<v.size(); i++){
        cout << v[i] << endl;
    }
    int a = 5;
    cout << returnInt(a) << endl;

    m2 << 1,1,1,1,
            1,1,1,1,
            1,1,1,1;

    cout << ((m2.transpose())*m2).trace() << endl;*/
}
