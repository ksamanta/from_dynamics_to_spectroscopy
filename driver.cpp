/*
**********************************************************************

    This is the driver program to run the Tully-style FSSH 
    trajectories for the Spin-Boson system.

    Author: Kousik Samanta 2012-10-25

**********************************************************************
*/

#include "spinboson.h"
#include <vector>

#ifdef _OPENMP 
    #include <omp.h> 
#else 
    #define omp_get_thread_num() 0 
#endif

using namespace std;


int main()
{
    //------------------------------------------------------------
    // Obtain the constants first
    //------------------------------------------------------------

    // Gather the input parameters in the SpinBosonInput struct
    const SpinBosonInput  SBI("inp.txt");
    const double DT = SBI.DT;
    const unsigned long MAX_TRAJ   = SBI.MAX_TRAJ;
    const unsigned long MAX_STEPS  = SBI.MAX_STEPS;
    const unsigned long MAX_POINTS = SBI.MAX_POINTS;

    // Define the SpinBoson instance SB and a dummy instance, DummySB
    SpinBoson SB(SBI), DummySB(SBI);


    //------------------------------------------------------------
    // Take a look at the PESs 
    //------------------------------------------------------------
    ofstream PESStream("pes.txt");  
    int surface = 0;
    dcomplex c[2] = { dcomplex(1.0,0.0), dcomplex(0.0,0.0)}; 
    double p = 20.0; // some arbitrary value
    for (long j=-10000; j< 10000; j++)
    {
        double x = double(j);
        SB.Init_vars(surface, x, p, c);
        SB.Print_PES(PESStream);
	}
    PESStream.close();


    //------------------------------------------------------------
    // Now start the simulation
    //------------------------------------------------------------

    // First, generate MAX_TRAJ random intial x and p values 
    // for MAX_TRAJ trajectories.
    // (Generate them all at once so that they conform to the 
    // correct distributions).

    RNG RNG_x, RNG_p;
    vector<double> x_i(MAX_TRAJ), p_i(MAX_TRAJ);
    const double SIGMA_X = sqrt(SBI.kT)/SBI.OMEGA0;
    const double MEAN_X = -sqrt(SBI.Er/2.0)/SBI.OMEGA0; 
    const double SIGMA_P = sqrt(SBI.kT);
    
    for (unsigned long traj = 1; traj <= MAX_TRAJ; traj++)
    {
        x_i[traj] = RNG_x.Sample_gaussian(SIGMA_X, MEAN_X);
        p_i[traj] = RNG_p.Sample_gaussian(SIGMA_P);
    }
    

    // Define registers to hold populations (initialized to zero)
    vector<double> total_pop_left(MAX_STEPS);
    vector<double> total_pop_right(MAX_STEPS);


    // Output stream
    ofstream OutStream("out.txt");  // Output file and the stream


    // Start OpenMP parallel block by explicitly specifying private 
    // and shared variables in order to minimize the risk of the 
    // so-called "race condition".

    #pragma omp parallel default(none) \
    firstprivate( SB, DummySB ) \
    shared(x_i, p_i, total_pop_left, total_pop_right, OutStream)
    {
        // Get the number of threads and the id of the current thread

        unsigned long num_threads = omp_get_num_threads();
        unsigned long thread_id   = omp_get_thread_num();

        // I'm choosing to print the info for the tread w/ highest ID
        if ( thread_id == (num_threads-1) ) 
        {
           OutStream << "Details for the thread w/ highest ID (#" 
           << thread_id << ")\n\n";
        }


        // Registers to store populations for the current thread
        vector<double> partial_pop_left(MAX_STEPS);
        vector<double> partial_pop_right(MAX_STEPS);


        // Give chunks of the trajectories to the running threads

        #pragma omp for schedule(static)
        for (unsigned long traj=1; traj <= MAX_TRAJ; traj++)
        {  
            // Feed in  the random initial values
            SB.Init_vars( x_i[traj], p_i[traj] );


            // Start the clock now
            for (unsigned long t=0; t < MAX_STEPS; t++) 
            {
                // Add up the contributions to the total population
                // due to the trajectories handled by the current 
                // thread 

                double* pop_d = SB.Diabatic_pop();
                partial_pop_left[t]  += pop_d[0];
                partial_pop_right[t] += pop_d[1];


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

        for (unsigned long t=0; t < MAX_STEPS; t++)
        {
            #pragma omp atomic
            total_pop_left[t]  += partial_pop_left[t];

            #pragma omp atomic
            total_pop_right[t] += partial_pop_right[t];
        }
        OutStream << "End of simulation.\n\n";
        OutStream.close();
    }


    //---------------------------------------------------------------
    // Print approximately MAX_POINTS  equally spaced data points.
    // (The output file gets way too big if you print all of them!)
    //---------------------------------------------------------------

   // Make decision on the interval of time points to be printed
    unsigned long t_step = MAX_STEPS/MAX_POINTS;  // integer division
    if ( t_step == 0  ) t_step = 1;

    // Find out the number of digits the highest time value may have
    const unsigned int PRECISION = (unsigned int)  log10(MAX_STEPS) ;
 
    // The output stream
    ofstream PopStream("pop.txt");  
    PopStream << "# Av pop for " << MAX_TRAJ << " trajectories:\n";

    // Now print, and after that close the output stream
    for (unsigned long t=0; t < MAX_STEPS; t += t_step)
    {
        PopStream << scientific << setprecision(PRECISION) 
            << setw(20) << double(t)*DT << "  " << setprecision(8)
            << setw(20) << total_pop_left[t]/double(MAX_TRAJ)  << "  " 
            << setw(20) << total_pop_right[t]/double(MAX_TRAJ) <<endl;
    }
    PopStream.close();


    // All done
    return 0;
}

