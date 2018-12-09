# A program to solve a symmetric double-minimum well for tunneling splittings.

import random
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import dkie_solver_function as solver

# Create a potential of the form:
# V(x, k, b, A) = 1/2 * k_force_constant * x**2 + A * exp(-(x/b)**2);

# Specify the number of harmonic oscillator basis functions to use.
nh_order = 16

# Physical Constants
h_planck = 6.626070040 * 10**-34
h_bar_planck = h_planck / (2 * math.pi)
cvel = 299792458

# Randomly pick constants (from a range of realistic constraints)
# to describe the potential.
k_force_constant = random.uniform(1, 100)
bAng = random.uniform(0.01, 0.5)    # in units of Angstroms
b = bAng * 10**-10                  # in meters
Acm = random.uniform(100, 10000)    # in wavenumbers
A = Acm * (h_planck * cvel * 100)   # in Hz


# Solve the potential and ouput tunneling splittings.
results = solver.solve_potential(nh_order, k_force_constant, b, Acm)

# Display the results.
print("RESULTS FROM SOLVING A SYM. DOUBLE-MINIMUM POTENTIAL IN AN HO BASIS")
print("--------------------------------------------------------------------\n")
print("Potential Form: " + "V = 1/2 * k * x**2 + A * exp(-(x/b)**2)")
print("Number of Harmonic Oscillator Basis Functions: " + str(results[0])+"\n")
print("Force Constant, k (N/m): " + "\t" + str(round(results[1],5)))
print("Parameter b (m): " + "\t\t" + str(round(results[2],15)))
print("Parameter A (N*m): " + "\t\t" + str(round(results[3],24)))
print("Minima Position (m): " + "\t\t" + str(round(results[4],15)))
print("Barrier Height (cm-1): " + "\t\t" + str(round(results[6],5)) + "\n")

print("PROTON RESULTS:")
print("---------------------")
print("Linear Freq: " + "\t" + str(round(results[8],5)))
#print("Number of Im. Freq: " + "\t" + str(results[11]))
print("Vibrational Frequencies (cm-1): " + "\n\t" + 
        str([round(x,2) for x in results[9]]))
print("Tunneling Splittings (cm-1): " + "\n\t" + 
        str([round(x,3) for x in results[17][1]]))

print("\nDEUTERON RESULTS:")
print("---------------------")
print("Linear Freq: " + "\t" + str(round(results[13],5)))
#print("Number of Im. Freq: " + "\t" + str(results[16]))
print("Vibrational Frequencies (cm-1): " + "\n\t" + 
        str([round(x,2) for x in results[14]]))
print("Tunneling Splittings (cm-1): " + "\n\t" + 
        str([round(x,3) for x in results[18][1]]))
        
print("\nDEUTERIUM KINETIC ISOTOPE EFFECTS:")
print("--------------------------------------")
print("DKIEs: " + "\n\t" + 
        str([round(x,3) for x in results[19]]))
                
# Plot the potential and wavefunctions.

# Create a list of evenly-spaced numbers over the range
if results[4] > 1 * 10**-11:
        if results[1] < 50:
                x_range = 3 * results[4]  
        else:
                x_range = 2 * results[4] 
else:
        x_range = 5 * 10**-11

x = np.linspace(-1 * x_range, x_range, 1000)

potential = 0.5 * k_force_constant * x**2 + A * np.exp(-(x/b)**2) - results[5]
potential /= (100 * h_planck * cvel)

# Add the locations of the energy levels to the plot.
for i in range(len(results[9])):
    if (i % 2) == 0:
        plt.hlines(results[9][i], -1 * x_range, x_range, colors = 'b')
    else:
        plt.hlines(results[9][i], -1 * x_range, x_range, colors = 'r')
plt.plot(x, potential, 'k-')

# Set the y scale.
if results[6] > 300:
        plt.ylim(-100, 2 * results[6])
else:
        plt.ylim(-100, 2000)

# Display the plot
plt.show()                   


