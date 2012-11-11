/*
***********************************************************************
	spinboson.h -- A class for the Spin Boson model
	Author: Kousik Samanta // 2012-10-04
***********************************************************************
*/

#ifndef SPINBOSON_H_
#define SPINBOSON_H_

#include<iostream>
#include<iomanip>
#include<fstream>
#include<complex>
#include <gsl/gsl_randist.h>
#include "rng.h"
#include "spinbosoninput.h"

using namespace std;

// Get a shorthand notation for double complex
typedef complex<double> dcomplex;	
typedef unsigned long ULONG;


// The SpinBoson class
//=====================================================================
class SpinBoson 
{
    private:

        // Constants for the spinboson problem
        double Er, kT, EPS0, V12, GAMMA, OMEGA0, OMEGA, M;
       
        // The input variables 
        //--------------------------------------------------------
        int surface;         // surface index 
        double x_val, p_val; // The dynamical variables x and p
        dcomplex c_val[2];   // coefs

        // Variables to store initial values
        int S_i;
        double x_i, p_i;
        dcomplex c_i[2];

        // Diabatic population
        double pop_d[2];

        // Other variables (computed by the methods in the class)
        //--------------------------------------------------------

        // V ==  PES, dVDx == dV/dx, dc == derivative coupling
        double V[2], dVdx[2], dc[2][2];  

        // Time derivatives
        double x_dot, p_dot;
        dcomplex c_dot[2];

        // Struct instances for the random number generators (rng.h).
        // All but the last are for Gaussian distributions.
        RNG rng_force, rng_uniform;

    public:

        // Constructor, copy constructor and destructor
        //--------------------------------------------------------
        SpinBoson(const SpinBosonInput&);
        SpinBoson(const SpinBoson&); 
        ~SpinBoson();						

       // The other methods
       //--------------------------------------------------------
       void Init_vars(double&, double&);
       void Init_vars(int&, double&, double&, dcomplex (&)[2]);
       void Get_PES();
       void Get_dVdx();
       void Get_derivative_coupling();
       void Get_time_derivatives(const double&);
       void Check_for_hopping(const double&);
       void Take_a_RK4_step(const double&, SpinBoson&);

       double Theta();
       double* Diabatic_pop();

       double Get_x(){return x_val; };
       void Print_terminal_points(ULONG&, ofstream& );
       void Print_PES(ofstream& );
}; 


#endif


