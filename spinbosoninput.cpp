/*
================================================================================

 This contains the constructors and destructor of the SpinBosonInput
 struct. The main purpose of this is to initialize the necessary
 variables or read them from a file.

 The default values for the input are given as sample below:

------------  Start sample input -------------------------------
    Er		  = 2.39e-2;  // Reorganization energy
    kT		  = 9.5e-4;   // Boltzmann const * Temp
    EPS0	  = 1.8e-2;   // Energy diff between the well minima
    V12		  = 1.25e-5;  // diabatic coupling
    GAMMA	  = 2.34e-6;  // frictional constant
    OMEGA0	  = 2.734e-6; // Freq. of the oscillator
    OMEGA	  = 1.0e0;	  // Freq. of incident radiation
    DT		  = 1.25; 	  // time step
    MAX_STEPS = 10000;	  // Maximum time to simulate
    MAX_TRAJ  = 10000;	  // Max number of trajectories to run
------------  End sample input -------------------------------
 
 NOTE:
 1. Key and the corresponding value must be on the same line
 2. Make sure to put an assignment operator ("=" or ":" in between the
    key and the value. 
 3. Don't put a non-whitespace character before the key.
 2. Blank lines and whitespaces are OK. 
 3. There are three ways to define a SpinBosonInput object: 

    A. Calling it witout an argument:
   ------------------------------------------------------------------
      If you just want to initialize the default values then 
    define a SpinBosonInput struct object without an argument as:

        SpinBosonInput SBI();  // SBI is the object

    B. Calling it by using an ifstream object as an argument:
   ------------------------------------------------------------------
    If you want to read the values from an input file (say, "inp.txt"),
    then, first include the fstream header as "#include<fstream>" on
    top, and then

      ifstream InpStream("inp.txt") // InpStream is an ifstream object
      SpinBosonInput SBI(InpStream) // SBI is again the SpinBosonInput
                                    // object

    C. Calling it by using a filename (a string) as an input:
   ------------------------------------------------------------------
    Alternatively it's also possible to initialize the values from an
    input file by directly passing the file name:

      SpinBosonInput SBI("inp.txt") // don't forget the double quotes!
      

 Author: Kousik Samanta (2012-10-19)
================================================================================
*/

#include<string>
#include<sstream>
#include "spinbosoninput.h"

//===========================================================================
// Inline scripts to trim a string
// Ref: http://www.cplusplus.com/faq/sequences/strings/trim/
//===========================================================================

inline string trim_right_copy(
    const string& s,
    const string& delimiters = " \f\n\r\t\v" )
{
    return s.substr( 0, s.find_last_not_of( delimiters ) + 1 );
}

inline string trim_left_copy(
    const string& s,
    const string& delimiters = " \f\n\r\t\v" )
{
    return s.substr( s.find_first_not_of( delimiters ) );
}

inline string trim_copy(
    const string& s,
    const string& delimiters = " \f\n\r\t\v" )
{
    return trim_left_copy( trim_right_copy( s, delimiters ), delimiters );
}

//===========================================================================




// Default constructor for the SpinBosonInput struct
//===========================================================================

SpinBosonInput::SpinBosonInput()
{
    // Default values
    Er        = 2.39e-2;    // Reorganization energy
    kT        = 9.5e-4;     // Boltzmann const * Temp
    EPS0      = 1.8e-2;     // Energy diff between the well minima
    V12       = 1.25e-5;    // diabatic coupling
    GAMMA     = 2.34e-6;    // frictional constant
    OMEGA0    = 2.734e-6;   // Freq. of the oscillator
    OMEGA     = 1.0e0;	    // Freq. of incident radiation
    DT        = 1.25;       // time step
    MAX_STEPS = 10000;	    // Maximum time to simulate
    MAX_TRAJ  = 10000;	    // Max number of trajectories to run
}


// Standard constructor with an ifstream object (input stream) as an argument
//===========================================================================
SpinBosonInput::SpinBosonInput(ifstream &InpStream)
{
    // Load the default values first
    *this = SpinBosonInput();

    // Define a string to hold the lines
    string line;

    // Read key in a loop
    while ( getline(InpStream, line)  ) 
    {

        // Check if an assignment op ('=' or ':' ) is in the line
        unsigned long loc_of_assignment_op = line.find_first_of("=:");
        if ( loc_of_assignment_op != string::npos ) 
        {
            string key  = line.substr(0, loc_of_assignment_op);
            string val = line.substr(loc_of_assignment_op+1); 
           
            // Check wheter the key is non-blank and val is numeric
            if (key.find_first_not_of(" \f\n\r\t\v") != string::npos
            && val.find_first_of("0123456789") != string::npos) 
            {
                // Strip the key off of the whitespaces around it
                string trimmed_key = trim_copy(key);

                // Define a stringstream object for the value, val
                istringstream ss_val(val);

                // Now assign values
                if (		trimmed_key == "Er"		)  ss_val >> Er;
                else if (trimmed_key == "kT"		)  ss_val >> kT;
                else if (trimmed_key == "EPS0"	)  ss_val >> EPS0;
                else if (trimmed_key == "V12"	)  ss_val >> V12;
                else if (trimmed_key == "GAMMA"	)  ss_val >> GAMMA;
                else if (trimmed_key == "OMEGA0"	)  ss_val >> OMEGA0;
                else if (trimmed_key == "OMEGA"	)  ss_val >> OMEGA;
                else if (trimmed_key == "DT"		)  ss_val >> DT;
                else if (trimmed_key == "MAX_STEPS") ss_val>>MAX_STEPS;
                else if (trimmed_key == "MAX_TRAJ")  ss_val >> MAX_TRAJ;
            }
        }
    }
}


// Standard constructor with a string (input file name) as an argument
//===========================================================================
SpinBosonInput::SpinBosonInput(string FileName)
{
    ifstream InpStream(FileName.c_str()); // pre-C++11 takes C-style str 
    *this = SpinBosonInput(InpStream);
}



// destructor
//===========================================================================
SpinBosonInput::~SpinBosonInput(){}

			
		

