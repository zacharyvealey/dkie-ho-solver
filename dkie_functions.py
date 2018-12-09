# A repository of functions useful for calculating DKIEs.

import math
import numpy as np

def energy_offset_location(k_force_constant, A, b):
    """A function to calculate the x-coordinates of the potential minima."""
    if k_force_constant <= (2 * A) / b**2:
        dom = np.log((2 * A) / (k_force_constant * b**2))
        e_offset_loc = b * math.sqrt(dom)
    else:
        e_offset_loc = 0
        
    return e_offset_loc

def energy_offset(k_force_constant, A, b):
    """A function to determine the minimum energy of the potential."""
    if k_force_constant <= (2 * A) / b**2:
        dom = np.log((2 * A) / (k_force_constant * b**2))
        e_offset = 0.5 * b**2 * k_force_constant * (1 + dom)
    else:
        e_offset = A
        
    return e_offset
    
def exp_ho_integral(vp, v, b):
    """
    The analytical solution for integrals involved in using a symmetric
    double-minimum potential of the form V = 0.5kx^2 + A*exp(-(x/b)^2)
    """
    if (vp % 2) == (v % 2):
        pL = v % 2
        
        # Preparing terms to input for calculating integral value.
        exp_1 = (vp - v) / 2
        exp_2 = (vp + v + 1) / 2
        
        sq_num = math.factorial((vp - pL) / 2) * math.factorial((v - pL) / 2)
        sq_denom = math.gamma((vp + pL + 1) / 2) * math.gamma((v + pL + 1) / 2)
        b_num = b**(2 * pL + 1)
        b_denom = (1 + b**2)**(exp_2)
        
        # Perform a sum from 0 to sum limit.
        sum_int = 0
        sum_limit = math.floor(min((vp - pL) / 2,(v - pL) / 2) + 1)
        
        for s in range(sum_limit):
            sum_num = math.gamma(exp_2 - s)
            sum_denom = math.factorial(s) * math.factorial((vp - pL) / 2 - s) 
            sum_denom *= math.factorial((v - pL) / 2 - s)
            sum_int += (b**4 - 1)**s * sum_num / sum_denom
        
        # Construct the matrix element for index vp,v of the Hamiltonian.        
        exp_ho_integral = (-1)**exp_1 * math.sqrt(sq_num / sq_denom) 
        exp_ho_integral *= (b_num / b_denom) * sum_int
        
    else:
        exp_ho_integral = 0
        
    return exp_ho_integral

def create_hamiltonian(nh_order, linear_freq, alpha, b, Acm):
    """
    A function to intialize and create a Hamiltonian describing the potential.
    """
    # Initialize the Hamiltonian.
    array_freq = []
    for diag in range(nh_order):
        array_freq.append(diag + 1 / 2)

    Hm = linear_freq * np.diag(array_freq, k=0)

    # Create Hamiltonian.
    for i in range(nh_order):
        for j in range(i, nh_order):
            Hm[i][j] += Acm * exp_ho_integral(i, j, alpha * b)
            if i != j:
                Hm[j][i] = Hm[i][j]
                
    return Hm
 
def get_eigenvalues(Hm, e_offset_cm):
    """ 
    A function to diagnolize the Hamiltonian and extract energies and
    eigenvectors.
    """
    # Obtain the eigenvalues and eigenvectors of the Hamiltonian.
    vals, vecs = np.linalg.eigh(Hm, UPLO='L')
    vals.sort()
    nrot = 0

    # Permute (eigen) vals and vecs in the event of imaginary frequencies.
    while vals[0] < 0:
        vals.append(vals.pop(0))
        vecs.append(vec.pop(0))
        nrot += 1
    
    # Correct the eigenvalues with the energy offset and report the tunneling
    # splittings.
    ovals = vals - e_offset_cm
    tun_splits = []

    for level in range(0, int(len(ovals)), 2):
        tun_splits.append(ovals[level+1] - ovals[level])
        
    return vals, vecs, nrot, ovals, tun_splits
 
 
