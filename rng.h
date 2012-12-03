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

#include <iostream>
#include <unistd.h>    // for getpid to work
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

struct RNG
{
    // Define the pointer for the random number generator
    gsl_rng *ptr;
 

    // The constructor for the struct
    RNG():
        ptr ( gsl_rng_alloc(gsl_rng_mt19937) ) // initialize pointer
    {
        unsigned long seed = time(NULL)*getpid(); // the seed
        gsl_rng_set(ptr, seed);                   // seed is set
    };


    // The destructor

    ~RNG()
    {
        gsl_rng_free(ptr);
    };


    // Get a random number sampled from a standard uniform distribution
    // (i.e., the one that lies in the interval [0,1] ).
    double Sample_uniform()
    {
        return gsl_rng_uniform(ptr);
    };
    

    // Get a random number sampled from a Gaussian distribution
    // defined by the mean and the standard deviation (sigma)
    double Sample_gaussian(double &sigma, double &mean)
    {
        return ( gsl_ran_gaussian(ptr, sigma) + mean );
     };


    // Overload Sample_gaussian
    double Sample_gaussian(const double &sigma, const double &mean)
    {
        return ( gsl_ran_gaussian(ptr, sigma) + mean );
    };

    // mean = 0 case
    double Sample_gaussian(const double &sigma)
    {
        return ( gsl_ran_gaussian(ptr, sigma) ) ; 
    };
};
       
#endif

