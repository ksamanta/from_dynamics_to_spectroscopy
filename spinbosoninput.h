/*
====================================================================
 This is the header for SpinBosonInput class. 
 This class processes the input for the Spin-Boson problem.
 See the usage guide at the top of SpinBosonInput.cpp (method file)

 Author: Kousik Samanta (2012-10-19)
====================================================================
*/

#ifndef SpinBosonInput_H_
#define SpinBosonInput_H_

#include<iostream>
#include<fstream>

using namespace std;

struct SpinBosonInput {

    double Er, kT, EPS0, V12, GAMMA, OMEGA0, OMEGA, DT;
    unsigned long MAX_TRAJ, MAX_STEPS, MAX_POINTS;

    SpinBosonInput();	     	// set the default values stored
    SpinBosonInput(ifstream& );	// call via ifstream object
    SpinBosonInput(string);		// call via a string (filename)
    ~SpinBosonInput();
};

#endif

