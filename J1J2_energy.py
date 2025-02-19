from Overlap import bcs_overlap, Sz,overlap_J1J2
import numpy as np
from scipy.optimize import minimize
from parameter import Nx,Ny,J
 

def energy_bcs_J1J2(eta,J):
    return (overlap_J1J2(J,eta,eta)/(bcs_overlap(eta,eta)))

# initial guess for eta

eta_initial = np.random.uniform(-1,1,Nx*Ny)   # intializes random eta
 
result = minimize(energy_bcs_J1J2, eta_initial, args=(J,), method='BFGS')

eta_optimized = result.x
final_energy = result.fun
print(f"sBCS J1J2 Energy: {final_energy}")
#print(eta_optimized)

Sz_sum = 0
for i in range(1,Nx*Ny+1):  # This calculates Sz for each site 

    #print(Sz(eta_optimized,eta_optimized,i)/bcs_overlap(eta_optimized,eta_optimized))
    Sz_sum += Sz(eta_optimized,eta_optimized,i)/bcs_overlap(eta_optimized,eta_optimized)

print("The global Sz is:",Sz_sum)


with open("energy_J1J2.txt","a") as file: # open a file and write the results to it
    file.write(f"{J} {final_energy}\n")