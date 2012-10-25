/*
***********************************************************************
	util.cpp -- A bunch of utility routines
	Author: Kousik Samanta // 2012-10-04
***********************************************************************
*/
#include<iostream>
#include<fstream>
#include<string>
#include<iomanip>
#include<complex>
using namespace std;



// Print a matrix
//---------------------------------------------------------------------

// real case
void print_mat (ofstream  & OutStream, string header, double * array, 
  int num_row, int num_col) {
	
	// Print the header first
	OutStream << header << endl;

	// Header row only if more than one column
	if (num_col > 1) {
		OutStream << setw(10) << "          ";
		for (int icol=0; icol<num_col; icol++)
			OutStream << setw(10) << icol;
		OutStream << endl;
	}
	
	// Now print the array
	for (int irow=0; irow<num_row; irow++) {
		// The Row label
		OutStream << "Row: " << setw(5) << irow << "|  ";

			
		// Now print each column entry for the current row
		for (int icol=0; icol<num_col; icol++) {
			OutStream << scientific << setw(10) << setprecision(3) 
				<< array[irow*num_col + icol];
		}
		OutStream << std::endl;
	}
	OutStream << std::endl;
}	


// complex case
void print_mat (ofstream  & OutStream, string header, complex<double> * array, 
  int num_row, int num_col) {
	
	// Print the header first
	OutStream << header << endl;

	// Header row only if more than one column
	if (num_col > 1) {
		OutStream << setw(10) << "          ";
		for (int icol=0; icol<num_col; icol++)
			OutStream << setw(25) << icol;
		OutStream << endl;
	}
	
	// Now print the array
	for (int irow=0; irow<num_row; irow++) {
		// The Row label
		OutStream << "Row: " << setw(5) << irow << "|  ";

			
		// Now print each column entry for the current row
		for (int icol=0; icol<num_col; icol++) {
			OutStream << scientific << setw(25) << setprecision(3) 
				<< array[icol*num_row + irow];
		}
		OutStream << std::endl;
	}
	OutStream << std::endl;
}	


