//============================================================
//  PURPOSE: linear congruential random number generator
//           based on the URAND Fortran routine discussed in
//           Chapter 10 of "Computer Methods for Mathematical
//           Computations" by Forsythe, Malcolm, and Moler,
//           which is in turn based on the discussion given in
//           "The Art of Computer Programming: Volume 2 --
//           Seminumerical Algorithms" by Donald E. Knuth.
//  LANGUAGE: C++
//  COMPILER: xlC for IBM RS6000 PowerStations
//  AUTHOR:   Henry D. Shapiro, shapiro@cs.unm.edu
//  DATE LATEST REVISION: January 12, 1998
//  ADDRESS: Department of Computer Science,
//           The University of New Mexico
//           Albuquerque, NM 87131
//============================================================

#ifndef RANDOM_NUMBERS_H
#define RANDOM_NUMBERS_H

// A class, each instance of which provides a stream of floating point
//   random numbers.

class RandomStream {
public:
    typedef int int32; // the random number generator assumes 32 bit integers

    enum Randomization { RANDOMIZE };
    enum RandomErrors { negSeed };

    // Create a stream of random numbers with a reproducible seed.
    RandomStream(int32 init_seed = DEFAULT_SEED);
        // Seed should be positive.

    // Create a stream of random numbers with a random seed.
    RandomStream(Randomization);
        // Effectively, this routine seeds the stream with the current time.

    ~RandomStream();

    // Return the next pseudorandom real number in the half-open interval [0,1).
    double next();

private:
    // Define away copy constructor and operator= for class  RandomStream .
    RandomStream(const RandomStream&);
    RandomStream& operator=(const RandomStream&);

    // Record the previous seed given out for randomized RandomStreams,
    //   so that if the clock hasn't ticked, we can give out a different
    //   seed and get a different stream.
    static const int32 DEFAULT_SEED;
    static int32 last_cal_time;
    static int32 current_increment;

    int32 x_subk;
};
#endif
