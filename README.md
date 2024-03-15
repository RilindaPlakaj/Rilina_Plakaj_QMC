# Rilina_Plakaj_QMC
Quantum Monte Carlo (QMC) Program
This program implements a Quantum Monte Carlo (QMC) simulation for various atomic and molecular systems. QMC methods are powerful computational techniques used to approximate solutions to the Schrödinger equation for many-body quantum systems.

Usage
Compile the program:  ifort -o main_qmc.x main_qmc.f90 subroutines.f90 functions.f90
Run the executable:  ./main_qmc.x

Follow the prompts printed on the screen to select the molecule of interest:
For Hydrogen atom, press 1
For Helium atom, press 2
For H2+ ion, press 3
For H2 molecule, press 4
For H3+ ion, press 5
Input Files
The program reads data from input files based on the selected molecule:
H_atom.dat for Hydrogen atom.
He_atom.dat for Helium atom.
H2+_ion.dat for H2+ ion.
H2_molecule.dat for H2 molecule.
H3+_ion.dat for H3+ ion.

Subroutines
pdmc: Performs the main Quantum Monte Carlo simulation.
drift: Evaluates the drift vector.
random_gauss: Generates random numbers following a Gaussian distribution.
ave_error: Computes the average and statistical error of an array.

Functions
potential: Computes the potential energy.
psi: Computes the wave function.
kinetic: Computes the local kinetic energy.
e_loc: Computes the local energy.

Output
Energy: Final estimated energy of the system.
Acceptance Rate (A): calculates the average and the variance of a given vector.

Notes
Ensure correct input files are provided for desired molecules.
Adjust simulation parameters (nmax, nruns, dt, tau, E_ref) as needed for accurate results.
