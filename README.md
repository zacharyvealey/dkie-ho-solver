# dkie-ho-solver
A program to solve for the deuterium kinetic isotope effects (DKIEs) for a given potential.

The program contained within dkie_potential_random_gen.py will randomly generate a potential that describes a molecular system from a range of realistic parameters. The functional form of the potential is given by a harmonic potential perturbed by a Gaussian contribution: 

V(x, k, b, A) = 1/2 * k * x** 2 + A * exp(-(x/b)** 2)

where k is the force constant, A & b are molecular parameters, and x describes the coordinate of interest. The program solves the potential using reformulated analytical matrix elements for a Gaussian perturbation (including correction for phase difference and redefinition of b) as first proposed in "Gaussian Perturbtaions to Harmonic Oscillators, S. Bell, J. Phys B 2, 1001 (1969).

The program performs this simulation using reduced masses of 1 and 2amu in order to compare the effects of an interaction dominated by proton or deuteron motion. For the symmetric double-minimum wells that are produced, the tunneling splitting is calculated from the bifurcation of the energy levels, and is directly proportional to the rate of the reaction. To determine the DKIE, the ratio of the tunneling splittings for the two species is taken.

The program outputs the relevant parameters describing the electronic potential as well as the vibrational frequencies and tunneling splittings deteremined for the isotopologs.

* Note that this program illustrates a continuum from high-barrier/deep-tunneling dynamics to the limit of no-barrier/pure-harmonic motion. As such, the meaning of the tunneling splitting becomes fluid as the zero-point energy surpasses the barrier height. Inspection of the DKIE values reflects this where phenomena dominated by tunneling interactions have DKIEs >> 7 whereas phenomena that approach pure harmonic behavior proceed toward Sqrt(2).
