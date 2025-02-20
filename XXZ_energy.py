from Overlap import *
import numpy as np
from scipy.optimize import minimize
import os
#from parameter import Delta
#import sys
Delta = -2.0 #float(sys.argv[1])  # Read Delta from command-line argument


def energy_bcs_XXZ(eta,Delta):
    return (overlap_XXZ(Delta,eta,eta))/(bcs_overlap(eta,eta))

eta_file = "eta_optimized.txt"

if os.path.exists(eta_file): # checks if the file exists or not
    eta_initial = np.loadtxt(eta_file) # load eta_optimized from the file

else:   # otherwise do an intial guess
    eta_initial = np.random.uniform(-1,1,Nx*Ny) 


# initialize minimization 
result = minimize(energy_bcs_XXZ, eta_initial, args=(Delta,), method='BFGS')

eta_optimized = result.x
final_energy = result.fun
print(f"sBCS XXZ Energy:{Delta}        {final_energy}")
#print(eta_optimized)

Sz_sum = 0
for i in range(1,Nx*Ny+1):  # This calculates Sz for each site 

    #print(Sz(eta_optimized,eta_optimized,i)/bcs_overlap(eta_optimized,eta_optimized))
    Sz_sum += Sz(eta_optimized,eta_optimized,i)/bcs_overlap(eta_optimized,eta_optimized)

#print("The global Sz is:",Sz_sum)

#np.savetxt(eta_file,eta_optimized) # save the eta_optimized in the .txt

with open(eta_file,"w") as file: # save the eta_optimized in the .txt
    for eta in eta_optimized:
        file.write(f"{eta}\n")

with open("energy_XXZ.txt","a") as file:  # saves the result to .txt file
    file.write(f"{Delta}    {final_energy}\n")
 


