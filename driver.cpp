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
    // Define the variables first
    //------------------------------------------------------------

    // Gather the input parameters in the SpinBosonInput struct
    const SpinBosonInput  SBInp("inp.txt");
    const ULONG MAX_STEPS = SBInp.MAX_STEPS;
    const ULONG MAX_TRAJ = SBInp.MAX_TRAJ;
    const double DT = SBInp.DT;


    // Define the SpinBoson instance SB and a dummy instance, DummySB
    SpinBoson SB(SBInp), DummySB(SBInp);

    // Define random number generators to generate initial x, p and c
    RNG rng_x, rng_p, rng_uniform;

    // Define registers to hold populations (initialized to zero)
    vector<double> partial_pop_left(MAX_STEPS);
    vector<double> partial_pop_right(MAX_STEPS);
    vector<double> total_pop_left(MAX_STEPS);
    vector<double> total_pop_right(MAX_STEPS);


    // Take a look at the PES first
    //------------------------------------------------------------
    ofstream PESStream("pes.txt");  
    dcomplex c_in[2] = { dcomplex(1.0,0.0), dcomplex(0.0,0.0)};
    for (ULONG j=-6000; j< 6000; j++)
    {
        SB.Init_vars(0, double(j), 20.0, c_in);
        SB.Print_PES(PESStream);
	}
    PESStream.close();


    // Now start the simulation
    //------------------------------------------------------------

    // Output stream
    ofstream OutStream("out.txt");  // Output file and the stream

    #pragma omp parallel firstprivate(partial_pop_left, partial_pop_right) 

    #pragma omp for firstprivate(SB, DummySB) 
    for (ULONG traj=1; traj <= MAX_TRAJ; traj++)
    {  
        int tid = omp_get_thread_num();
        ostringstream oss; 
        oss << " traj " << traj << " tid: " << tid << endl;
        OutStream << oss.str() << endl;

        // Set the random initial variables
        SB.Init_vars(rng_x, rng_p, rng_uniform);
 

        // Start the clock now
        for (ULONG t=0; t < MAX_STEPS; t++) 
        {
            // Add up the contributions to the total population due to
            // the trajectories handled by the current thread 
            total_pop_left[t]  += SB.Diabatic_pop('L');
            total_pop_right[t] += SB.Diabatic_pop('R');

            // Check if hopping is feasible
            SB.Check_for_hopping(DT);

            // Take a Runge-Kutta step and update the dyanmical
            // variables
            SB.Take_a_Runge_Kutta_step(DT, DummySB);
        }

        // Now add up the contributions from different threads
        // 'atomically' (in order to avoids the "race condiiton")
        for (ULONG t=0; t < MAX_STEPS; t++)
        {
            #pragma omp atomic
            total_pop_left[t]  += partial_pop_left[t];
            #pragma omp atomic
            total_pop_right[t] += partial_pop_right[t];
        }
    }

    ofstream PopStream("pop.txt");  // open/reopen for writing
    PopStream << "# Av pop for " << MAX_TRAJ << " trajectories:\n";

    for ( ULONG t=0; t < MAX_STEPS; t++ )
    {
        PopStream << scientific << setw(10) << t << "  " 
        << scientific << setw(20) << setprecision(8) 
            << total_pop_left[t]/double(MAX_TRAJ) << "   " 
        << scientific << setw(15) << setprecision(5) 
            << total_pop_right[t]/double(MAX_TRAJ) << endl;
    }
    PopStream.close();
    OutStream.close();
	return 0;
}

