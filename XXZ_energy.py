from Overlap import *
import numpy as np
from scipy.optimize import minimize
import os
from parameter import *
#import sys
#Delta = -1.0 #float(sys.argv[1])  # Read Delta from command-line argument
#Nsite = Nx*Ny
#print(f"The total number of sites are {Nsite}")

def energy_bcs_XXZ(eta,Delta):
    return (overlap_XXZ(Delta,eta,eta))/(bcs_overlap(eta,eta))

eta_file = "eta_optimized.txt"

if os.path.exists(eta_file): # checks if the file exists or not
    eta_initial = np.loadtxt(eta_file) # load eta_optimized from the file

else:   # otherwise do an intial guess
    eta_initial = np.random.uniform(-1,1,Nx*Ny) 

def Sz_sum(eta):
    Sz_global= 0
    for i in range(1,Nsites+1):  # This calculates Sz for each site 
    #print(Sz(eta_optimized,eta_optimized,i)/bcs_overlap(eta_optimized,eta_optimized))  
        Sz_global += Sz(eta,eta,i)/bcs_overlap(eta,eta)
    return Sz_global

# Define constraint dictionary
constraint = {'type': 'eq', 'fun': Sz_sum}

# Perform minimization with Sz = 0 constraint
result = minimize(energy_bcs_XXZ, eta_initial, args=(Delta,), method='SLSQP', constraints=constraint)

# initialize minimization 
#result = minimize(energy_bcs_XXZ, eta_initial, args=(Delta,), method='BFGS')

eta_optimized = result.x
final_energy = result.fun

print(f"The Energy for {Delta}:  {final_energy}")
#print(f"The global Sz is: {Sz_sum(eta_optimized)}")

#np.savetxt(eta_file,eta_optimized) # save the eta_optimized in the .txt

with open(eta_file,"w") as file: # save the eta_optimized in the .txt
    for eta in eta_optimized:
        file.write(f"{eta}\n")

with open("energy_XXZ_1D_12.txt","a") as file:  # saves the result to .txt file
    file.write(f"{Delta}    {final_energy:.12f}\n")

    
  