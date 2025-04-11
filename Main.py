from overlap_uv import *
from parameter import *
import numpy as np
from Gradient import Total_Gradient_XXZ_1D
from scipy.optimize import minimize
import os
import sys




def Energy(var,Nsites,Delta):  # function to calculate XXZ_Energy in 1D
    ham.theta = theta
    ham.Nsites = Nsites
    theta = var[:Nsites]
    phi = var[Nsites:]
    return (ham.J1J2_2D_overlap(Delta))  #Change the Hamiltonian accordingly
    #return (ham.XXZ_overlap(Delta))
    
def Sz_sum(theta,Nsites):
    ham.theta = theta
    Sz_global= 0
    for i in range(Nsites):  # This calculates Sz for each site 
    #print(Sz(eta_optimized,eta_optimized,i)/bcs_overlap(eta_optimized,eta_optimized))  
        Sz_global += ham.Sz(i)
    return Sz_global/ham.bcs_overlap()


#---------------------------------------------------------------------------------------------------------#
'''This section loads the etas of the previous calculation as an initial guess for the current calculation.
If the file doesn't exist, it creates a folder'''


'''theta_file = 'theta_opt.txt'

if os.path.exists(theta_file):
    theta_initial = np.loadtxt(theta_file)

else:
     theta_initial = np.random.uniform(-np.pi,np.pi,Nsites)'''

best_obj = np.inf
best_theta = None

for i in range(100):
    theta0 = np.random.uniform(-np.pi,np.pi, Nsites)

    #theta0 = np.array([np.pi/4,-np.pi/4,np.pi/4,-np.pi/4,np.pi/4,-np.pi/4,
                      # np.pi/4,-np.pi/4,np.pi/4,-np.pi/4,np.pi/4,-np.pi/4,
                      # np.pi/4,-np.pi/4,np.pi/4,-np.pi/4])
    ham = BCSHamiltonian(theta0, Nx, Ny)  # The Hamiltonian
# ----------------------------------------------------------------------------------------------------------------#
# Prepare the result data to write to file
    '''log_file = 'log.txt'
    sys.stdout = open(log_file,'w')'''
# ----------------------------------------------------------------------------------------------------------------#
    E_before = Energy(theta0,Nsites,Delta)


#-----------------------------------------------------------------------------------------------------------------#
 
# Define constraint dictionary
    constraint = {'type': 'eq', 'fun': lambda theta :Sz_sum(theta,Nsites)}  #Using Lambda function since Sz_sum has two arguments

# Perform minimization with Sz = 0 constraint
    #result = minimize(XXZ_1D_energy, theta0, args=(Nsites,Delta,), method='trust-constr',jac=Total_Gradient_XXZ_1D, constraints=constraint,\
    #              options={'xtol':1e-14,'maxiter':2000})
    result = minimize(Energy, theta0, args=(Nsites,Delta,), constraints=constraint,\
                  options={'xtol':1e-14,'maxiter':2000})
    #theta_optimized = result.x
    final_energy = result.fun
    print(f"run {i+1}: Energy_before : {E_before} Energy : {final_energy}")

    if final_energy < best_obj:
        best_obj = final_energy
        best_theta = result.x

#--------------------------------------------------------------------------------------------------------------------#

 
'''with open(theta_file,'w') as file:     #overwriting the etas in the .txt file
    for theta in theta_optimized:
        file.write(f"{theta}\n")'''

#sys.stdout.close()
'''for theta in theta_optimized:
    print(f"{np.cos(theta)}\t{np.sin(theta)}\n")'''

with open('energy_4x4_XXZ_new.txt',"a") as file:      # writing the energy to a text file
    file.write(f"{Delta}   {best_obj:.12f}\n")

sys.stdout = sys.__stdout__

if result.success:
    print('optimization successful!')
else:
    print("Warning! optimization failed: ", result.message)
ham.theta = best_theta
print(ham.bcs_overlap())
print(f"The Energy for {Delta}:  {best_obj:.12f}")
print(f"The global Sz is: {Sz_sum(best_theta,Nsites)}")