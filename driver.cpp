#include "spinboson.h"
#include<vector>

using namespace std;

int main()
{
	// Define parameters for simulation
	//------------------------------------------------------------

	// Output stream
	ofstream OutStream("out.txt");  // Output file and the stream

    // Gather the input parameters in the SpinBosonInput struct
    const SpinBosonInput  SBInp("inp.txt");
    const size_t MAX_STEPS = SBInp.MAX_STEPS;
    const size_t MAX_TRAJ = SBInp.MAX_TRAJ;
    const double DT = SBInp.DT;

	// Initiate the population registers to be zero
    vector<double> total_pop_left(MAX_STEPS+1);
    vector<double> total_pop_right(MAX_STEPS+1);


	// Define the SpinBoson instance SB and a temp instance DummySB
	SpinBoson SB(SBInp), DummySB(SBInp);

    // Take a look at the PES first
	//------------------------------------------------------------
	ofstream PESStream("pes.txt");  
	dcomplex c_in[2] = { dcomplex(1.0,0.0),
					     dcomplex(0.0,0.0)};
	for (long j=-6000; j< 6000; j++){
		SB.Set_specific_xpc(0, double(j), 20.0, c_in);
		SB.Print_PES(PESStream);
	}

	// Now start the simulation
	//------------------------------------------------------------
	for (size_t traj=1; traj <= MAX_TRAJ; traj++)
    {
		// Generate starting position, x0 and momentum, p0
		SB.Set_random_xpc();
        OutStream << "Start of traj# " << traj << endl;
        SB.Print_xpc(OutStream);

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
			SB.Take_a_Runge_Kutta_step(DT, DummySB);
		}

        // Print the endpoints of the trajectory
        OutStream << "End of traj# " << traj << endl;
        SB.Print_xpc(OutStream);

       // print the average populations in increments of 100 trajectories
       // as well as after all the trajectories are done
       if ( traj%100 == 0 || traj == MAX_TRAJ)
       {
           ofstream PopStream("pop.txt");  // open/reopen for writing
           PopStream << "# Av pop for " << traj << " trajectories:\n";

           for ( size_t t=0; t < MAX_STEPS; t++ )
           {
                 PopStream << scientific << setw(10) << t << "  " 
                     << scientific << setw(15) << setprecision(5) 
                         << total_pop_left[t]/double(traj) << "   " 
                     << scientific << setw(15) << setprecision(5) 
                         << total_pop_right[t]/double(traj) << endl;
           }
           PopStream.close();
       } 
	}

	return 0;
}

