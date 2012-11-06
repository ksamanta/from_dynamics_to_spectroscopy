/*
**********************************************************************

    This is the driver program to run the Tully-style FSSH 
    trajectories for the Spin-Boson system.

    Author: Kousik Samanta 2012-10-25

**********************************************************************
*/

#include "spinboson.h"
#include<vector>
using namespace std;

int main()
{
    // Define the variables first
    //------------------------------------------------------------

    // Gather the input parameters in the SpinBosonInput struct
    const SpinBosonInput  SBInp("inp.txt");
    const size_t MAX_STEPS = SBInp.MAX_STEPS;
    const size_t MAX_TRAJ = SBInp.MAX_TRAJ;
    const double DT = SBInp.DT;

    // Initialize the population registers with zero pop
    vector<double> total_pop_left(MAX_STEPS+1);
    vector<double> total_pop_right(MAX_STEPS+1);


    // Define the SpinBoson instance SB and a dummy instance, DummySB
    SpinBoson SB(SBInp), DummySB(SBInp);

    // Some dummy vecotrs
    complex_vec dummy_c(2), sum_c(2);

    // Take a look at the PES first
    //------------------------------------------------------------
    ofstream PESStream("pes.txt");  
    complex_vec c_in;
    c_in.push_back( dcomplex(1.0,0.0) );
    c_in.push_back( dcomplex(0.0,0.0) );
    for (long j=-6000; j< 6000; j++)
    {
        int surface=0;
        double x = double(j);
        double p = 20.0;
        SB.Set_specific_xpc(surface, x, p, c_in);
        SB.Print_PES(PESStream);
	}
    PESStream.close();

    // Now start the simulation
    //------------------------------------------------------------

    // Output stream
    ofstream OutStream("out.txt");  // Output file and the stream

    for (size_t traj=1; traj <= MAX_TRAJ; traj++)
    {
        // Generate starting position, x0 and momentum, p0
        SB.Set_random_xpc();
 
        //Print the details of the current time-step as a "footprint"
        SB.Footprints("Start", traj, OutStream); 

       // Start the clock now
       for (size_t t=0; t < MAX_STEPS; t++) 
       {

            // Update the population registers for this time step
            total_pop_left[t]  += SB.Diabatic_pop('L');
            total_pop_right[t] += SB.Diabatic_pop('R');

            // Check if hopping is feasible
            SB.Check_for_hopping(DT);

            // Take a Runge-Kutta step and update the dyanmical
            // variables
            SB.Take_a_Runge_Kutta_step(DT, DummySB, dummy_c, sum_c);
        }
        // Print the endpoints
        SB.Footprints("End", traj, OutStream);


        // print the average populations in increments of 100 trajectories
        // as well as after all the trajectories are done
        if ( traj%500 == 0 || traj == MAX_TRAJ)
        {
            ofstream PopStream("pop.txt");  // open/reopen for writing
            PopStream << "# Av pop for " << traj << " trajectories:\n";

            for ( size_t t=0; t < MAX_STEPS; t++ )
            {
                PopStream << scientific << setw(10) << t << "  " 
                << scientific << setw(20) << setprecision(8) 
                << total_pop_left[t]/double(traj) << "   " 
                << scientific << setw(15) << setprecision(5) 
                << total_pop_right[t]/double(traj) << endl;
            }
            PopStream.close();
        } 
    }
    OutStream.close();
	return 0;
}

