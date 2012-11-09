/*
***********************************************************************
	RNG.h -- A random number generator struct
    It creates and seeds a random number generator. 
    Defing it as "RNG rng;"  and then access the pointer as
    pointer as RNG.ptr (this is what you need to use GSL random number
    distributions).

	Author: Kousik Samanta // 2012-10-24
***********************************************************************
*/

#ifndef RNG_H_
#define RNG_H_

#include<iostream>
#include <unistd.h>    // for getpid to work
#include <gsl/gsl_rng.h>

using namespace std;

struct RNG
{
    // Define the pointer for the random number generator
    gsl_rng *ptr;
 
    // The constructor for the struct
    RNG()
    {
        ptr = gsl_rng_alloc( gsl_rng_mt19937 );  // allocation
        long seed = time(NULL) * getpid();     // the seed
        gsl_rng_set(ptr, seed);                   // seed is set
     }

     // The destructor
     ~RNG()
     {
        gsl_rng_free(ptr);
     }
};
       
#endif

