/*
**********************************************************************

    This is the driver program to run the Tully-style FSSH 
    trajectories for the Spin-Boson system.

    Author: Kousik Samanta 2012-10-25

**********************************************************************
*/

#include "spinboson.h"
#include <vector>
#include <omp.h>

using namespace std;

int main()
{
    // Obtain the constants first
    //------------------------------------------------------------

    // Gather the input parameters in the SpinBosonInput struct
    const SpinBosonInput  SBI("inp.txt");
    const double DT        = SBI.DT;
    const ULONG MAX_TRAJ   = SBI.MAX_TRAJ;
    const ULONG MAX_STEPS  = SBI.MAX_STEPS;
    const ULONG MAX_POINTS = SBI.MAX_POINTS;

    // Define the SpinBoson instance SB and a dummy instance, DummySB
    SpinBoson SB(SBI), DummySB(SBI);


    // Take a look at the PESs 
    //------------------------------------------------------------
    ofstream PESStream("pes.txt");  
    int surface = 0;
    dcomplex c[2] = { dcomplex(1.0,0.0), dcomplex(0.0,0.0)};
    double p = 20.0;
    for (long j=-6000; j< 6000; j++)
    {
        double x = double(j);
        SB.Init_vars(surface, x, p, c);
        SB.Print_PES(PESStream);
	}
    PESStream.close();


    // Now start the simulation
    //------------------------------------------------------------

    // First, generate MAX_TRAJ random intial x and p values 
    // (Generate them all at once so that they conform to the 
    // correct distributions).
    const double Er=SBI.Er, kT=SBI.kT, OMEGA0=SBI.OMEGA0;
    const double SIGMA_X = sqrt(kT)/OMEGA0, SIGMA_P = sqrt(kT);
    const double MEAN_X = -sqrt(Er/2.0)/OMEGA0; 
    vector<double> x_i(MAX_TRAJ), p_i(MAX_TRAJ);
    RNG rng_x, rng_p;
    
    for (ULONG traj = 1; traj <= MAX_TRAJ; traj++)
    {
        x_i[traj] = gsl_ran_gaussian(rng_x.ptr, SIGMA_X) + MEAN_X;
        p_i[traj] = gsl_ran_gaussian(rng_p.ptr, SIGMA_P); // MEAN_P=0
    }
    
    // Define registers to hold populations (initialized to zero)
    vector<double> partial_pop_left(MAX_STEPS);
    vector<double> partial_pop_right(MAX_STEPS);
    vector<double> total_pop_left(MAX_STEPS);
    vector<double> total_pop_right(MAX_STEPS);


    // Output stream
    ofstream OutStream("out.txt");  // Output file and the stream

    // Start OpenMP parallel block by explicitly specifying private 
    // and shared variables in order to minimize the risk of the 
    // so-called "race condition".
    #pragma omp parallel default(none) \
        firstprivate ( SB, DummySB ) \
        firstprivate ( partial_pop_left, partial_pop_right ) \
        shared (x_i, p_i, total_pop_left, total_pop_right, OutStream) 
    {
        // Get the number of threads and the id of the current thread
        ULONG num_threads = omp_get_num_threads();
        ULONG thread_id   = omp_get_thread_num();
        if ( thread_id == (num_threads-1) ) 
        {
           OutStream << "Details for the thread w/ highest ID (#" 
           << thread_id << ")\n\n";
        }

        // Give chunks of the trajectories to the running threads
        #pragma omp for 
        for (ULONG traj=1; traj <= MAX_TRAJ; traj++)
        {  
            // Feed in  the random initial values
            SB.Init_vars( x_i[traj], p_i[traj] );

            // Start the clock now
            for (ULONG t=0; t < MAX_STEPS; t++) 
            {
                // Add up the contributions to the total population due to
                // the trajectories handled by the current thread 
                double* pop_d = SB.Diabatic_pop();
                total_pop_left[t]  += pop_d[0];
                total_pop_right[t] += pop_d[1];

                // Check if hopping is feasible
                SB.Check_for_hopping(DT);

                // Take a Runge-Kutta step and update the dyanmical
                // variables
                SB.Take_a_RK4_step(DT, DummySB);
            }

            // In order to get some idea about the progress of the
            // calculation print the terminal points of the
            // trajectories handled by the last thread
            if ( thread_id == (num_threads-1) )
                SB.Print_terminal_points(traj, OutStream);
        }

        // Now add up the contributions from different threads
        // 'atomically' (in order to avoids the "race condition")
        for (ULONG t=0; t < MAX_STEPS; t++)
        {
            #pragma omp atomic
            total_pop_left[t]  += partial_pop_left[t];
            #pragma omp atomic
            total_pop_right[t] += partial_pop_right[t];
        }
        OutStream << "End of simulation.\n\n";
        OutStream.close();
    }


    // Print only  MAX_STEPS / (MAX_STEPS/MAX_POINTS) equally spaced
    // data points (Notice the integer division!).
    // Otherwise, the output file gets way too big!
    ULONG t_step = MAX_STEPS/MAX_POINTS;  // integer division
    if ( t_step == 0  ) t_step = 1;
    const int PRECISION = log10(MAX_STEPS);
 
    ofstream PopStream("pop.txt");  
    PopStream << "# Av pop for " << MAX_TRAJ << " trajectories:\n";

    for ( ULONG t=0; t < MAX_STEPS; t += t_step )
    {
        PopStream << scientific << setprecision(PRECISION) 
        << setw(20) << double(t)*DT << "  " << setprecision(8)
        << setw(20) << total_pop_left[t]/double(MAX_TRAJ)  << "  " 
        << setw(20) << total_pop_right[t]/double(MAX_TRAJ) << endl;
    }
    PopStream.close();
    return 0;
}

