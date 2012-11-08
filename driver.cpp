/*
**********************************************************************

    This is the driver program to run the Tully-style FSSH 
    trajectories for the Spin-Boson system.

    Author: Kousik Samanta 2012-10-25

**********************************************************************
*/

#include "spinboson.h"
#include<vector>
#include<omp.h>
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

    // Initialize the population registers with zero pop
    // Leave enough room to store all the time steps on each
    // trajectories
    vector<double> pop_left(MAX_TRAJ*MAX_STEPS);
    vector<double> pop_right(MAX_TRAJ*MAX_STEPS);


    // Define the SpinBoson instance SB and a dummy instance, DummySB
    SpinBoson SB(SBInp), DummySB(SBInp);

    // Take a look at the PES first
    //------------------------------------------------------------
    ofstream PESStream("pes.txt");  
    dcomplex c_in[2] = { dcomplex(1.0,0.0), dcomplex(0.0,0.0)};
    for (ULONG j=-6000; j< 6000; j++)
    {
        SB.Set_specific_xpc(0, double(j), 20.0, c_in);
        SB.Print_PES(PESStream);
	}
    PESStream.close();

    // Now start the simulation
    //------------------------------------------------------------

    // Output stream
    ofstream OutStream("out.txt");  // Output file and the stream

    #pragma omp parallel for firstprivate(SB, DummySB) 
    for (ULONG traj=1; traj <= MAX_TRAJ; traj++)
    {  
        int tid = omp_get_thread_num();
        ostringstream oss; 
        oss << " traj " << traj << " tid: " << tid << endl;
        OutStream << oss.str() << endl;

        // Generate starting position, x0 and momentum, p0
        SB.Set_random_xpc();
 
        //Print the details of the current time-step as a "footprint"
        //SB.Footprints("Start", traj, OutStream); 

        // Start the clock now
        for (ULONG t=0; t < MAX_STEPS; t++) 
        {
            // Update the population registers for this time step
            // pack the 2d array into a 1-d array
            ULONG traj_t = (traj-1) * MAX_TRAJ + t; 
            pop_left[traj_t]  = SB.Diabatic_pop('L');
            pop_right[traj_t] = SB.Diabatic_pop('R');

            // Check if hopping is feasible
            SB.Check_for_hopping(DT);

            // Take a Runge-Kutta step and update the dyanmical
            // variables
            SB.Take_a_Runge_Kutta_step(DT, DummySB);
        }
        // Print the endpoints
        //SB.Footprints("End", traj, OutStream);
    }

    ofstream PopStream("pop.txt");  // open/reopen for writing
    PopStream << "# Av pop for " << MAX_TRAJ << " trajectories:\n";

    for ( ULONG t=0; t < MAX_STEPS; t++ )
    {
        double tot_pop_left = 0.0;
        double tot_pop_right = 0.0;

        for (ULONG traj=1; traj <= MAX_TRAJ; traj++)
        {
            ULONG traj_t = (traj-1) * MAX_TRAJ + t; 
            tot_pop_left += pop_left[traj_t];
            tot_pop_right += pop_right[traj_t];
        }
 
        PopStream << scientific << setw(10) << t << "  " 
        << scientific << setw(20) << setprecision(8) 
        << tot_pop_left/double(MAX_TRAJ) << "   " 
        << scientific << setw(15) << setprecision(5) 
        << tot_pop_right/double(MAX_TRAJ) << endl;
    }
    PopStream.close();
    OutStream.close();
	return 0;
}

