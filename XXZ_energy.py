from Overlap import *
import numpy as np
from scipy.optimize import minimize
from parameter import Delta
 

def energy_bcs_XXZ(eta,Delta):
    return (overlap_XXZ(Delta,eta,eta))/(bcs_overlap(eta,eta))

# initial guess for eta

eta_initial = np.random.uniform(-1,1,Nx*Ny)   # intializes random eta
 
result = minimize(energy_bcs_XXZ, eta_initial, args=(Delta,), method='BFGS')

eta_optimized = result.x
final_energy = result.fun
print(f"sBCS XXZ Energy: {final_energy}")
#print(eta_optimized)

Sz_sum = 0
for i in range(1,Nx*Ny+1):  # This calculates Sz for each site 

    #print(Sz(eta_optimized,eta_optimized,i)/bcs_overlap(eta_optimized,eta_optimized))
    Sz_sum += Sz(eta_optimized,eta_optimized,i)/bcs_overlap(eta_optimized,eta_optimized)

print("The global Sz is:",Sz_sum)
 



