# A program to solve a symmetric double-minimum well for tunneling splittings.

import math
import numpy as np
import dkie_functions as df

def solve_potential(nh_order, k_force_constant, b, Acm):

    # Define physical constants needed.
    amu = 1.660539040 * 10**-27
    h_planck = 6.626070040 * 10**-34
    h_bar_planck = h_planck / (2 * math.pi)
    cvel = 299792458

    # Specify the parameters used to construct double-minimum potential well.
    A = Acm * (h_planck * cvel * 100)

    # Determine the location of the minima and the zero-energy offset.
    e_offset_loc = df.energy_offset_location(k_force_constant, A, b)
    e_offset = df.energy_offset(k_force_constant, A, b)
    e_offset_cm = e_offset / (h_planck * cvel * 100)

    # Calculate the barrier height of the potential.
    eff_barr_height = A - e_offset
    eff_barr_height_cm = Acm - e_offset_cm

    proton_rates = []
    deuteron_rates = []

    for m in range(1,3):
        mass = m * amu

        # Compute the characteristics describing the potential.
        angular_freq = math.sqrt(k_force_constant / mass)
        linear_freq = angular_freq / (2 * math.pi * cvel * 100)
        alpha = (mass * k_force_constant / h_bar_planck**2)**(1/4)

        # Create Hamiltonian.
        Hm = df.create_hamiltonian(nh_order, linear_freq, alpha, b, Acm)

        # Obtain the eigenvalues and eigenvectors of the Hamiltonian.
        vals, vecs, nrot, ovals, tun_splits = df.get_eigenvalues(Hm, e_offset_cm)

        if m == 1:
            proton_rates.append(m)
            proton_rates.append(tun_splits[:5])
            proton_ovals = ovals
            proton_vecs = vecs
            proton_nrot = nrot
            proton_ang_freq = angular_freq
            proton_lin_freq = linear_freq
        elif m == 2:
            deuteron_rates.append(m)
            deuteron_rates.append(tun_splits[:5])
            deuteron_ovals = ovals
            deuteron_vecs = vecs
            deuteron_nrot = nrot
            deuteron_ang_freq = angular_freq
            deuteron_lin_freq = linear_freq

    dkie = []
    if len(proton_rates[1]) == len(deuteron_rates[1]):
        for i in range(len(proton_rates[1])):
            dkie.append(proton_rates[1][i] / deuteron_rates[1][i])
    else:
        print("You have a different number of tunneling splittings for the " +
                "proton and deuteron!")

    results = [nh_order, k_force_constant, b, A, e_offset_loc, e_offset,
                eff_barr_height_cm, proton_ang_freq, proton_lin_freq,
                proton_ovals, proton_vecs, proton_nrot, deuteron_ang_freq,
                deuteron_lin_freq, deuteron_ovals, deuteron_vecs, 
                deuteron_nrot, proton_rates, deuteron_rates, dkie]

    return results





