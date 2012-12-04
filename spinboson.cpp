/*
***********************************************************************
	spinboson.cpp -- Various "methods" for the SpinBoson class
	Author: Kousik Samanta // 2012-10-04.
***********************************************************************
*/

#include "spinboson.h"


//----------------------------------------------------------------------
// Init_vars --- Initialize the dynamic variables, x, p and c, and
// choose the adiabatic surface.
//----------------------------------------------------------------------

// Given initial x and p, compute c such that the starting point is on
// the left well. Depending on the c's, choose an adiabatic surface.

void SpinBoson::Init_vars(double &X, double &P)
{
    x_val = X;
    p_val = P;

    // C_lower and C_upper from C_left=1, C_right=0
    double theta = Theta();
    double pop_0 = (1.0 - cos(theta) ) / 2.0;
    double pop_1 = 1.0 - pop_0;
    c_val[0] = dcomplex( sqrt(pop_0), 0.0 );
    c_val[1] = dcomplex( sqrt(pop_1), 0.0);

    // surface
    double rand_num = RNG_uniform.Sample_uniform();
    if ( pop_0 >= rand_num ) 
        surface = 0;
    else 
        surface = 1;

   // Store these information as initial variables
   S_i    = surface;
   x_i    = x_val;
   p_i    = p_val;
   c_i[0] = c_val[0];
   c_i[1] = c_val[1];
}


// Overload it to set specific values of all the variables

void SpinBoson::Init_vars(int &S, double &X, double &P,
dcomplex (&C)[2] )
{
    surface  = S;
    x_val    = X;
    p_val    = P;
    c_val[0] = C[0];
    c_val[1] = C[1];
}


//----------------------------------------------------------------------
// Theta --- get theta = arctan( V12 / (M*x_val+EPS0/2.0) )
//----------------------------------------------------------------------

inline double SpinBoson::Theta()
{
    return atan2( V12, (M*x_val+EPS0/2.0) );
}

		
//----------------------------------------------------------------------
// Get_PES -- Compute the adiabatic potential energy surfaces
//----------------------------------------------------------------------

inline void SpinBoson::Get_PES()
{
    double first_part = (OMEGA0*OMEGA0*x_val*x_val - EPS0)/2.0;	
    double second_part = sqrt( (M*x_val+EPS0/2.0)*(M*x_val+EPS0/2.0)
          + V12*V12	);

    V[0] = first_part - second_part;
    V[1] = first_part + second_part;
}


//----------------------------------------------------------------------
// Get_grad_PES -- Compute the gradients of the adiabatic PESs
//----------------------------------------------------------------------

inline void SpinBoson::Get_grad_PES()
{
    double first_part = OMEGA0*OMEGA0*x_val;
    double second_part = M*(M*x_val+EPS0/2.0)
            / sqrt( (M*x_val+EPS0/2.0)*(M*x_val+EPS0/2.0) + V12*V12 );

    grad_PES[0] = first_part - second_part;
    grad_PES[1] = first_part + second_part;
}


//----------------------------------------------------------------------
// Get_derivative_coupling -- Compute the adiabatic derivative coupling: 
// We are actually computing 
//     U^H(dU/dt) = U^H [ (\partial U)/(\partial t) ]
//                + U^H [ (\partial U/(\partial x) ] * (p/m)
// where U diagonalizes the effective Hamiltonian in presence of the
// field. We are also assuming that the initial non-adiabatic coupling
// is not very strong, so it can be neglected. (In absense of the
// field, the first term is zero and the second is the same as the
// scalar product of regular derivative coupling and velocity.)
//----------------------------------------------------------------------

inline void SpinBoson::Get_derivative_coupling()
{
    // Constat x part
    double denom =( (M*x_val+EPS0/2.0)*(M*x_val+EPS0/2.0) + V12*V12 );
    double d01 = -0.5 * M * V12 / denom;

    dcomplex half_iw = 0.5 * dcomplex(0.0,1.0) * OMEGA;
    double cos_theta = cos( Theta() );
    double sin_theta = sin( Theta() );

    dc[0][0] = half_iw*cos_theta;
    dc[1][1] = half_iw*(-cos_theta);
    dc[0][1] = half_iw*sin_theta + p_val*d01; // m = 1
    dc[1][0] = half_iw*sin_theta - p_val*d01; // m = 1
}


//----------------------------------------------------------------------
// Get_time_derivatives ----Get dx/dt, dp/dt, dc/dt
//----------------------------------------------------------------------

void SpinBoson::Get_time_derivatives(const double &Random_force)
{
    //  dx/dt
    //..................................................................
    x_dot = p_val;   // mass is chosen to be 1.

    // dp/dt
    //..................................................................
    // Newtonian_force 
    Get_grad_PES();
    double Newtonian_force = -grad_PES[surface];

    //Frictional force
    double Frictional_force = -GAMMA * p_val;

    //Total force
    p_dot = Newtonian_force + Frictional_force + Random_force; //m=1


    // dc/dt
    //..................................................................

    //Get the PES and the derivative couplings
    Get_PES();
    Get_derivative_coupling();

    //Now compute dc/dt
    for ( int k=0; k<2; k++ )
    {
        dcomplex sum = dcomplex(0.0,1.0)*V[k]*c_val[k]; 
	    for (int j=0; j<2; j++)
            sum += dc[k][j] * c_val[j];
	    c_dot[k] = -sum;
	}
}
		

//----------------------------------------------------------------------
// Check_for_hopping -- Check if a hopping (FSSH style) is possible,
// and if it is, then do hop (i.e., change the surface)
//----------------------------------------------------------------------

void SpinBoson::Check_for_hopping(const double &dt)
{	
    // The current and new (to jump on to) surfaces
    int k = surface;      // current surface
    int j = 1 - surface;  // the other surface
    int other_surface = j;

    // Get the required elements of the density matrix, aa =  c x c^+
    dcomplex aa_kj = c_val[k] * conj( c_val[j] );  
    double real_aa_kk = real( c_val[k] * conj( c_val[k] ) ); 

    
    // dc here refers to the derivative coupling in presence of
    // the field, especially U^H(dU/dt) --- deriv w/ respect to t.
    // (See the actual function for details).
    Get_derivative_coupling();  // get dc(i,j)


    // b_jk (Eq. 14) and g_kj (Eq. 19) refer to Tully's FSSH paper:
    // "Molecular dynamics with electronic transitions"
    // J.C. Tully, J. Chem. Phys. 93, 1061 (1990).

    double b_jk = -2.0 * real( dc[j][k] * aa_kj );

    double g_kj; 
    if ( abs(real_aa_kk) < 1.0e-12 )
        g_kj = 0.0;
    else if ( b_jk/real_aa_kk < 0.0 )
        g_kj = 0.0;
    else
        g_kj = dt * b_jk / real_aa_kk;


    // Check for hopping against a random number
    double rand_num = RNG_uniform.Sample_uniform();

    if (g_kj > rand_num) 
    {
        Get_PES();  // Get V(i)
        double Delta_PE = V[other_surface] - V[surface]; //PE change
        double KE_new = (p_val*p_val/2.0) - Delta_PE; // for w=0 only

        // These are the energetic criteria for hoppping:
        // 1. surface=1: it's always ok to go down to surface=0.
        // 2. surface=0: a hop to surface=1 is feasible if mechanical 
        // energy is conserved (in absence of the laser field) OR the 
        // energy supplied by the laser field is greater than or equal
        // to the difference in the energy between the two surfaces.
        // If the excitation (1 <-- 0) happens due to the laser field
        // do not change the momentum.
        if (surface == 1 || 
           (surface == 0 && abs(OMEGA) < 1.0e-12 && KE_new > 1.0e-12))
        {
            double pSign = p_val / abs(p_val);
            p_val = pSign * sqrt(2.0*KE_new); //don't change direction
            surface = other_surface;
        }
        else if ( surface == 0 && OMEGA >= Delta_PE )
        {
            surface = other_surface;
        }
    }
}

	
//--------------------------------------------------------------------
// Take_a_RK4_step -- Take a 4-th order Runge-Kutta step and 
// change the dyanmical variables.
// NOTE: You need to initialize a dummy SpinBoson instance, DummySB
// (one of the arguments to be passed) BEFORE you can use this method)
//--------------------------------------------------------------------

void SpinBoson::Take_a_RK4_step(const double &dt, 
SpinBoson &DummySB)
{		
    // Generate a random force for the Langevin-type dynamics
    // (it's the same for all RK micro-steps corresponding to one 
    // single big RK step.)
    const double sigma_force  = sqrt(2.0*GAMMA*kT/dt);
    const double random_force =RNG_force.Sample_gaussian(sigma_force);


    // The initial dynamic variables for DummySB
    double Dummy_x = x_val;
    double Dummy_p = p_val;  
    dcomplex Dummy_c[2] = {c_val[0], c_val[1]};

    // Initialize the integration sums
    double sum_x = x_val;
    double sum_p = p_val;
    dcomplex sum_c[2] = { c_val[0], c_val[1] };

    // Now take the RK4 micro steps 
    for (int RK_order = 1; RK_order <= 4; RK_order++)
    {
        // Set/reset x, p and c in RK
        DummySB.Init_vars(surface, Dummy_x, Dummy_p, Dummy_c);

        // Get the time derivates
        DummySB.Get_time_derivatives(random_force);

        // Some necessary factors for RK4 integration
        double k_fac=0.5, sum_fac=1.0/6.0;//initilize to avoid warning
        switch(RK_order) 
        {
           case 1: sum_fac = 1.0/6.0; k_fac = 0.5; break;
           case 2: sum_fac = 1.0/3.0; k_fac = 0.5; break;
           case 3: sum_fac = 1.0/3.0; k_fac = 1.0; break;
           case 4: sum_fac = 1.0/6.0; k_fac = 1.0; break; 
        }
				
        // Compute Runge-Kutta parameters for x
        double kx = dt * DummySB.x_dot;  // RK parameter (k)
        sum_x    += sum_fac * kx;        // Update x
        Dummy_x   = x_val + k_fac * kx;  // inp for the next RK_order

        // Compute Runge-Kutta parameters for p
        double kp = dt * DummySB.p_dot;
        sum_p    += sum_fac * kp;
        Dummy_p   = p_val + k_fac * kp;

        // Compute Runge-Kutta parameters for c
        for (int j=0; j<2; j++)
        {
            dcomplex kc = dt * DummySB.c_dot[j];
            sum_c[j]   += sum_fac * kc; 
            Dummy_c[j]  = c_val[j] + k_fac * kc;
        }
    }

    // Return the integrations sums
    x_val = sum_x;
    p_val = sum_p;
    c_val[0] = sum_c[0];
    c_val[1] = sum_c[1];
}


//----------------------------------------------------------------------
// Get_diabatic_population -- returns pop_d, a private array 
//----------------------------------------------------------------------

double* SpinBoson::Diabatic_pop()
{
    const double sin2_half_theta = (1.0 - cos(Theta()) ) / 2.0 ;
    const double cos2_half_theta = (1.0 - sin2_half_theta);


    // Define populations, and set the population of the active
    // adiabatic surface to be 1 and that of the inactive to be 0.
    // (Suffix "d" means diabatic, and "a" adiabatic)

    double pop_a[2];
    pop_a[surface]   = 1.0;
    pop_a[1-surface] = 0.0;

    pop_d[0] =  pop_a[0]*sin2_half_theta + pop_a[1]*cos2_half_theta;
    pop_d[1] =  pop_a[0]*cos2_half_theta + pop_a[1]*sin2_half_theta;
    
    return pop_d;
}


//--------------------------------------------------------------------
// Print_terminal_points -- print out the variables at the current 
// time step
//--------------------------------------------------------------------

void SpinBoson::Print_terminal_points(unsigned long &index,
ofstream &OutStream)
{
    // Compute norms
    double norm_i = real(  sqrt( c_i[0]*conj(c_i[0]) 
                               + c_i[1]*conj(c_i[1]) )  );
    double norm_f = real(  sqrt( c_val[0]*conj(c_val[0]) 
                               + c_val[1]*conj(c_val[1]) )  );

    // Now print in the output stream
    OutStream << "Traj: "
        << setw(6) << index 
        << scientific << setprecision(3)
        << "  INITL  S: "  << setw(2) << S_i 
        << ",  x: " << setw(10) << x_i
        << ",  p: " << setw(10) << p_i 
        << ",  c0: (" << setw(10) << real(c_i[0])
            << ", " << setw(10) << imag(c_i[0]) << ")"
        << ",  c1: (" << setw(10) << real(c_i[1])
            <<", " << setw(10) << imag(c_i[1]) << ")"
        << ",  |1-Norm|: " << setw(10) << abs(1.0-norm_i)
        << endl 
    << setw(12) << " "
        << "  FINAL  S: "  << setw(2) << surface
        << ",  x: " << setw(10) << x_val
        << ",  p: " << setw(10) << p_val
        << ",  c0: (" << setw(10) << real(c_val[0])
            << ", " << setw(10) << imag(c_val[0]) << ")"
        << ",  c1: (" << setw(10) << real(c_val[1])
            << ", " << setw(10) << imag(c_val[1]) << ")"
        << ",  |1-Norm|: " << setw(10) << abs(1.0 - norm_f)
        << endl 
        << endl;
}
		

//----------------------------------------------------------------------
// Print_PES --- print the adiabatic energies
//----------------------------------------------------------------------

void SpinBoson::Print_PES (ofstream &OutStream)
{
    // Get the PES and derivative couplings first
    Get_PES();
    Get_derivative_coupling();

    // Now print
    OutStream << setw(16) << x_val << "  " 
    << scientific << setprecision(8)
    << setw(16) << V[0] << "  " 
    << setw(16) << V[1] << "  "
    << setw(16) << real(dc[0][1]) << " " 
    << setw(16) << imag(dc[0][1]) << endl;
}


