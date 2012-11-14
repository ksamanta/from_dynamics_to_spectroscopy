From dynamics to spectroscopy
=============================

Spectroscopy using non-adiabatic dynamics

1. driver.cpp
-------------
The main driver routine. It parses the input gathering the constants for
the spin-boson problem. Then it initializes each trajectory with random
x and p on the left diabatic well. Then it propagates these
quantum-classical trajectories using 4-th order Runge-Kutta numerical
integration technique at each time-step.


2. spinbosoninput.h and spinbosoninput.cpp
------------------------------------------
Thsese files constitute the structure SpinBosonInput. Basically, the
constructor of the structure parses the input file provided to store the
constants for the spin-boson problem.


3. spinboson.h and spinboson.cpp
--------------------------------
These constitute the main SpinBoson class. This class is initialized by
a SpinBosonInput object. The methods of this class are capable of doing
the following
  -Set specific values of x, p and c (as well as surface).
  -Compute the PESs, gradient of PESs and derivative coupling.
  -Compute the time derivatives, dx/dt, dp/dt and dc/dt.
  -Using the method to compute time derivative, take Runge-Kutta step.


4. rng.h
--------
It constitutes the structure, RNG that creates a pointer for a random
number generator as well as seeds the generator. It also has two
different methods to get a random a number from a standard uniform and
gaussian distributions.

