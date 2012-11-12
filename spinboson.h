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
        const double Er, kT, EPS0, V12, GAMMA, OMEGA0, OMEGA, M;
       
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

        // Rrandom number generator objects (see rng.h).
        RNG RNG_force, RNG_uniform;

    public:

        // The constructor method
        SpinBoson(const SpinBosonInput& SBI): 
            Er (SBI.Er),
            kT (SBI.kT),
            EPS0 (SBI.EPS0),
            V12 (SBI.V12),
            GAMMA (SBI.GAMMA),
            OMEGA0 (SBI.OMEGA0),
            OMEGA (SBI.OMEGA),
            M (SBI.OMEGA0*sqrt(SBI.Er/2.0))
         {};

        // The copy constructor method
        SpinBoson(const SpinBoson& SB): 
            Er (SB.Er),
            kT (SB.kT),
            EPS0 (SB.EPS0),
            V12 (SB.V12),
            GAMMA (SB.GAMMA),
            OMEGA0 (SB.OMEGA0),
            OMEGA (SB.OMEGA),
            M (SB.OMEGA0*sqrt(SB.Er/2.0))
         {};

        // The destructor method
        ~SpinBoson(){};						

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


