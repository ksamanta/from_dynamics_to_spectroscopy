/*
***********************************************************************
	spinboson.cpp -- Various "methods" for the SpinBoson class
	Author: Kousik Samanta // 2012-10-04
***********************************************************************
*/
#include "spinboson.h"


// The constructor and the destructor
//======================================================================

SpinBoson::SpinBoson(const SpinBosonInput & SBI)
{
    // Defin the SpinBoson constants
    Er 		= SBI.Er;
    kT 		= SBI.kT;
    EPS0 	= SBI.EPS0;
    V12		= SBI.V12;
    GAMMA  	= SBI.GAMMA;
    OMEGA0	= SBI.OMEGA0;
    OMEGA 	= SBI.OMEGA;
    M		= OMEGA0 * sqrt(Er/2.0);

	// Choose the random number generator as a pointer
  	rng_ptr_x 		= gsl_rng_alloc(gsl_rng_mt19937);  
  	rng_ptr_p 		= gsl_rng_alloc(gsl_rng_mt19937);  
  	rng_ptr_surf 	= gsl_rng_alloc(gsl_rng_mt19937);  
  	rng_ptr_hop 	= gsl_rng_alloc(gsl_rng_mt19937);  
  	rng_ptr_force 	= gsl_rng_alloc(gsl_rng_mt19937);  

	// Seed the generators
	size_t seed = time(NULL) * getpid();  // seed
  	gsl_rng_set(rng_ptr_x, seed);           
  	gsl_rng_set(rng_ptr_p, 2*seed);         
	gsl_rng_set(rng_ptr_surf, 3*seed);
	gsl_rng_set(rng_ptr_force, 4*seed);
	gsl_rng_set(rng_ptr_hop, 5*seed);
}


// The destructor
SpinBoson::~SpinBoson()
{
    gsl_rng_free(rng_ptr_x);
    gsl_rng_free(rng_ptr_p);
    gsl_rng_free(rng_ptr_surf);
    gsl_rng_free(rng_ptr_hop);
    gsl_rng_free(rng_ptr_force); 
}



// Other methods
//======================================================================


// Set_random_xpc --- set random values of x and p 
// (x is taken from a Boltzmann distribution of energy and p is from 
// a Maxwell-Boltzmann // distribution of p vectors)
//----------------------------------------------------------------------
void SpinBoson::Set_random_xpc()
{
    // Parameters for equilibrium position and moentum distributions
    const double SIGMA_X = sqrt(kT)/OMEGA0;
    const double SIGMA_P = sqrt(kT);
    const double MEAN_X  = -M / (OMEGA0*OMEGA0);

    // Set the random values of x and p
    x_val = gsl_ran_gaussian(rng_ptr_x, SIGMA_X) + MEAN_X;
    p_val = gsl_ran_gaussian(rng_ptr_p, SIGMA_P);   // Mean is 0

    // C_lower and C_upper from C_left=1, C_right=0
    double theta = Theta(x_val);
    c_val[0] = dcomplex( sin(theta/2.0), 0.0 );
    c_val[1] = dcomplex( cos(theta/2.0), 0.0 );

    // surface
    double probability_lower = sin(theta/2.0) * sin(theta/2.0);
	double rand_num = gsl_rng_uniform(rng_ptr_surf);
    if ( probability_lower >= rand_num ) 
        surface = 0;
    else 
        surface = 1;
}


// Set_specific_xpc --- set specific values of x, p and c
//----------------------------------------------------------------------
void SpinBoson::Set_specific_xpc(int Surface, double X, double P,
dcomplex* C)
{
    surface  = Surface;
    x_val 	 = X;
    p_val 	 = P;
    c_val[0] = C[0];
    c_val[1] = C[1];
}


// Theta --- where tan(theta) = 2.0*V12/ (2*M*x_val+EPS0)
//----------------------------------------------------------------------
double SpinBoson::Theta(double x)
{
    double height = 2.0 * V12;
    double base   = 2.0 * M * x + EPS0;
    return atan2(height, base);
}

		

// Get_PES -- Compute the adiabatic potential energy surfaces
//----------------------------------------------------------------------
void SpinBoson::Get_PES()
{
    double first_part = OMEGA0*OMEGA0*x_val*x_val - EPS0;	
    double second_part = sqrt(  (2.0*M*x_val + EPS0)
                        *(2.0*M*x_val + EPS0)
                        + 4.0*V12*V12	);


    V[0][0] = (first_part - second_part)/2.0;
    V[1][1] = (first_part + second_part)/2.0;
    V[0][1] = V[1][0] = 0.0;
}


// Get_dVdX -- Compute the spatial derivative of the adiabatic PESs
//----------------------------------------------------------------------
void SpinBoson::Get_dVdx()
{
    double first_part = OMEGA0*OMEGA0*x_val;
    double second_part = 2.0*M*(2.0*M*x_val+EPS0)
        / sqrt( (2.0*M*x_val+EPS0) * (2.0*M*x_val+EPS0) 
                + 4.0*V12*V12  );
	
    dVdx[0] = first_part - second_part;
    dVdx[1] = first_part + second_part;
}


// Get_derivative_coupling -- Compute the adiabatic derivative coupling: 
// just the d01 (AKA, d12) element
//----------------------------------------------------------------------
void SpinBoson::Get_derivative_coupling(){

	double numer = 2.0 * M * V12;   
	double denom = (2.0*M*x_val+EPS0) * (2.0*M*x_val+EPS0) 
					+ 4.0*V12*V12;
	double d01 = numer / denom;   // Check the sign

	dc[0][1] =  d01;
	dc[1][0] = -d01;
	dc[0][0] = dc[1][1] = 0.0;
}



// Get_time_derivatives ----Get dx/dt, dp/dt, dc/dt
//----------------------------------------------------------------------
void SpinBoson::Get_time_derivatives(double Random_force){

	// dxdt
	//..................................................................
	dxdt = p_val;   // mass is chosen to be 1.

	// dpdt  
	//..................................................................

	// Newtonian_force 
	Get_dVdx();
	double Newtonian_force = -dVdx[surface];

	// Frictional force
	double Frictional_force = -GAMMA * p_val;

	// Total force
	dpdt = Newtonian_force + Frictional_force + Random_force;
	

	// dc/dt
	//..................................................................

	// Define complex zero and i (imaginary unit)
	const dcomplex cZERO(0.0, 0.0);
	const dcomplex i(0.0, 1.0); // i

	// First, u = (U^dag) * dU/dT
	dcomplex iw_over_two = i * OMEGA / 2.0;
    double theta = Theta(x_val);
	dcomplex sw = iw_over_two * sin(theta);
	dcomplex cw = iw_over_two * cos(theta);

	dcomplex u[2][2] = { {cw, sw}, 
		                 {sw,-cw} 	};
	
	// Get the PES and the derivative couplings
	Get_PES();
	Get_derivative_coupling();
	
		
	// Now compute dc/dt
	dcomplex sum;
	for ( int k=0; k<2; k++ ){
		sum = cZERO;
		for (int j=0; j<2; j++){
			sum -= ( i*V[k][j] + dxdt*dc[k][j] + u[k][j] ) 
				 * c_val[j];
		}
		dcdt[k] = sum;
	}
}
		

	
// Check_for_hopping -- Check if a hopping (FSSH style) is possible,
// and if it is, then do hop and change the necessary variables
// (dt is the time-step)
//----------------------------------------------------------------------
void SpinBoson::Check_for_hopping(const double dt)
{	
	// The current and new (to jump on to) surfaces
	int k = surface;       // current surface
	int j = 1 - surface;   // the other surface
	int other_surface = j;

	// First compute the necessary metrics to judge if a hopping is
	// feasible

	// Get the densities
	dcomplex aa[2][2];
	
    for (int m=0; m<2; m++)
    {
        for (int n=0; n<2; n++)
            aa[m][n] = c_val[m] * conj(c_val[n]);
    }

	// b_jk (Eq. 14) and g_kj (Eq. 19) // Ref. Tully's paper
	Get_derivative_coupling();  // get dc(i,j)

	double b_jk = -2.0 * real( conj(aa[j][k]) * (p_val/M) * dc[j][k] );
	double g_kj = dt * b_jk / real(aa[k][k]);

    if (g_kj < 0.0) 
        g_kj = 0.0;
		

	// Generate a standard uniform random number 
	// and check for hopping against it
	double rand_num = gsl_rng_uniform(rng_ptr_hop);
	Get_PES();  // Get V(i,j)

	if (g_kj > rand_num) {

		double PE_gain 	= V[other_surface][other_surface] 
				- V[surface][surface];
		double KE_new = p_val * p_val / (2.0 * M) - PE_gain;

		// Update the momentum in case a hop is really feasible
		if ((surface == 0 && KE_new + OMEGA > 0.0) 
		|| 	(surface == 1 && KE_new > 0.0) ) {
			double pSign = p_val / abs(p_val);
			p_val = pSign * sqrt (2.0 * M * KE_new);
			surface = other_surface;
		}
	}
}	

	
// Take_a_Runge_Kutta_step -- Take a 4-th order Runge-Kutta step and 
// change the dyanmical variables
//----------------------------------------------------------------------
void SpinBoson::Take_a_Runge_Kutta_step(const double dt, 
SpinBoson & DummySB)
{		
	// Generate a random force for the Langevin-type dynamics
	double sigma_force = sqrt(2.0*GAMMA*kT/dt);
	double random_force = gsl_ran_gaussian(rng_ptr_force, sigma_force);

    // The initial dynamic variables for DummySB
	double Dummy_x = x_val;
    double Dummy_p = p_val;  
    dcomplex Dummy_c[2] = {c_val[0], c_val[1]};

	for (int RK_order = 1; RK_order <= 4; RK_order++){

		// Set/reset x, p and c in RK
		DummySB.Set_specific_xpc(surface, Dummy_x, Dummy_p, Dummy_c);

		// Get the time derivates
		DummySB.Get_time_derivatives(random_force);

       // Some necessary factors for RK integration (see below)
	   double k_fac, sum_fac;
       switch(RK_order) 
       {
           case 1: sum_fac = 1.0/6.0; k_fac = 0.5; break;
           case 2: sum_fac = 1.0/3.0; k_fac = 0.5; break;
           case 3: sum_fac = 1.0/3.0; k_fac = 1.0; break;
           case 4: sum_fac = 1.0/6.0; k_fac = 1.0; break; 
       }
				
	   // Compute Runge-Kutta parameters for x
       double kx = dt * DummySB.dxdt;  // RK parameter (k)
       x_val    += sum_fac * kx;   	   // Add to the integration sum
       Dummy_x   = x_val + k_fac * kx; // Get input for the next RK_order

       // Compute Runge-Kutta parameters for p
       double kp = dt * DummySB.dpdt;
       p_val    += sum_fac * kp;
       Dummy_p   = p_val + k_fac * kp;

		// Compute Runge-Kutta parameters for c
        for (int j=0; j<2; j++)
        {
           dcomplex kc = dt * DummySB.dcdt[j];
           c_val[j]   += sum_fac * kc;
           Dummy_c[j]  = c_val[j] + k_fac * kc;
       }
   }
}


// Get_diabatic_population --
//----------------------------------------------------------------------
double SpinBoson::Diabatic_pop(char well)
{
    // Theta is related to the eigenvectors
    const double cos_theta = cos( Theta(x_val) );

    // Define populations, and set the population of the current
    // adiabatic surface to be 1.
    double pop_d = 0.0, pop_a[2];
    pop_a[surface]   = 1.0;
    pop_a[1-surface] = 0.0;

    if (well == 'L' || well == 'l')  // left diabatic surface
    {
        pop_d =  pop_a[0] * ( 1 - cos_theta ) / 2.0 
               + pop_a[1] * ( 1 + cos_theta ) / 2.0;
    } 
    else if (well == 'R' || well == 'r') // right diabatic surface
    {
        pop_d =  pop_a[0] * ( 1 + cos_theta ) / 2.0 
               + pop_a[1] * ( 1 - cos_theta ) / 2.0;
    } 
    return pop_d;
}

// Print_xpc -- print out the dynamical variables, x, p and c
//----------------------------------------------------------------------
void SpinBoson::Print_xpc(ofstream & OutStream){

   OutStream << " sur: " << setw(2) << surface 
       << "   x: " << fixed << setw(10) << setprecision(3) << x_val 
       << "   p: " << fixed << setw(10) << setprecision(3) << p_val
       << "   C: " << scientific << setw(15) << setprecision(3) << c_val[0] 
       << "  " << scientific << setw(15) << setprecision(3) << c_val[1] 
       << "  tot pop: " << scientific << setw(15) << setprecision(3) 
       << sqrt( c_val[0]*conj(c_val[0]) + c_val[1]*conj(c_val[1]) ) 
       << endl;
}
		
// Print_PES --- print the adiabatic energies
//----------------------------------------------------------------------
void SpinBoson::Print_PES (ofstream & OutStream){
    Get_PES();
    Get_derivative_coupling();
    OutStream << x_val << " " << V[0][0] << "  " << V[1][1] << " "<< dc[0][1]<<endl;
}



