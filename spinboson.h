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


        // Other variables (computed by the methods in the class)
        //--------------------------------------------------------

        // V ==  PES, dVDx == dV/dx, dc == derivative coupling
        double V[2], dVdx[2], dc[2][2];  

        // Time derivatives
        double dxdt, dpdt;
        dcomplex dcdt[2];

       // Structs for the random number generators (see rng.h)
       RNG rng_x, rng_p, rng_surf, rng_hop, rng_force;  

    public:

        // First the constructor and then the destructor
        //--------------------------------------------------------
        SpinBoson(const SpinBosonInput&);
        ~SpinBoson();						

       // The other methods
       //--------------------------------------------------------
       void Set_random_xpc();
       void Set_specific_xpc(int, double, double, dcomplex *);
       void Get_PES();
       void Get_dVdx();
       void Get_derivative_coupling();
       void Get_time_derivatives(double);
       void Check_for_hopping(const double);
       void Take_a_Runge_Kutta_step(const double, SpinBoson&);
       double Cos_theta(double);
       double Diabatic_pop(char);

       void Print_xpc(ofstream& );
       void Print_PES(ofstream& );
}; 


#endif


