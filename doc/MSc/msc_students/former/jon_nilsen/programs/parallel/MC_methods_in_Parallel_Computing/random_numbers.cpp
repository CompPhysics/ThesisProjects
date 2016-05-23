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
//  DATE LATEST REVISION: October 16, 2000
//  ADDRESS: Department of Computer Science,
//           The University of New Mexico
//           Albuquerque, NM 87131
//============================================================

// A class, each instance of which provides a stream of floating point
//   random numbers.
// The implementation assumes 32 bit two's complement arithmetic that is
//   reasonably well behaved.

#include <time.h>

#include "random_numbers.h"

const RandomStream::int32 RandomStream::DEFAULT_SEED = 1;
RandomStream::int32 RandomStream::current_increment = 1;

// Create a stream of random numbers with a reproducible seed.
RandomStream::RandomStream(int32 init_seed) : x_subk(init_seed) {
    // Seed should be positive.
    if (init_seed <= 0) throw negSeed;
} // RandomStream::RandomStream(int32 init_seed)

// Create a stream of random numbers with a random seed.
// Effectively, this routine seeds the stream with the current time
//   (in seconds).  If the time is unavailable use the default seed.
RandomStream::RandomStream(Randomization) {
    int32 cal_time = int32(time(0));
    // Do different things if time(0) works or fails
    if (cal_time >= 1)
        // time(0) works
        // Make sure that x_subk is not the same as it was in the previous
        //   instantiation, even if the clock hasn't ticked or if so many
        //   streams were requested that they have overwhelmed a clock tick.
        //   (The use of current_increment is essentially a bug fix for a
        //   bug discovered when a student asked for multiple streams in
        //   rapid succession.)
        x_subk = cal_time + current_increment;
    else
        // time(0) fails
        x_subk = DEFAULT_SEED + current_increment;
    current_increment++;

    // Protect against overflow.
    if (x_subk <= 0) throw negSeed;
} // RandomStream::RandomStream(Randomization);

RandomStream::~RandomStream() { }

// Return the next pseudorandom real number in the half-open interval [0,1).
double RandomStream::next() {
    /* A linear congruential random number generator, that is, one where
                       x[k+1] = a*x[k] + c    mod m
         Constants a, c, and m are based on the URAND Fortran routine discussed
         in Chapter 10 of "Computer Methods for Mathematical Computations"
         by Forsythe, Malcolm, and Moler, which is in turn based on the
         discussion given in "The Art of Computer Programming:
         Volume 2 -- Seminumerical Algorithms" by Donald E. Knuth.

         The routine presented here has been specialized to work on 32-bit
         two's-complement machines in which arithmetic is "well behaved".
         In particular, it is assumed that interrupts are ignored on
         integer overflow.  Other "abnormal" arithmetic behavior may
         invalidate the behavior of this code.

       m = 2^31, or one larger than the largest positive number
           representable on the machine
    */
    // halfm = m/2 = 2^30 = the largest positive power of 2 representable on
    //   the machine
    static const int32 halfm = int32(1) << 30;
    static const int32 a =  843314861;
    static const int32 c =  453816693;

    // Compute Next Random Number
    x_subk = a*x_subk + c;

    // Due to the way overflow is handled in the hardware, the modulo m
    //   has been done automatically, except the result can be negative.
    if (x_subk < 0)
      x_subk = (x_subk + halfm) + halfm; // x_subk + m, but done carefully

    // Return the answer as a real number in [0,1).
    return x_subk/(2.0*halfm);
} // double RandomStream::next()
