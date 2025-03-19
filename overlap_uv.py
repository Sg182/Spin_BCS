from neighbor import *
import numpy as np
# Only for Real u,v

def bcs_overlap(theta,Nsites):   # computes overlap between two BCS states
    overlap = 1
    
    for k in range(Nsites):
        overlap *= ((np.cos(theta[k]))**(2) + (np.sin(theta[k]))**(2))
    return overlap


def S_x(theta,N,p):
    cos_part_p = (np.cos(theta[p]))
    sin_part_p = (np.sin(theta[p]))
    prefactor =  (1/2)*(2*cos_part_p*sin_part_p)
    prefactor_Sx = prefactor*bcs_overlap(theta,N)
    return prefactor_Sx

def Sz(theta,Nsites,i):     # calculates Sz for site i
    prefactor_Sz = (1/2)*((-1*(np.cos(theta[i]))**(2) + (np.sin(theta[i]))**(2)))/\
        ((np.cos(theta[i]))**(2) + (np.sin(theta[i]))**(2))
    return bcs_overlap(theta,Nsites)*prefactor_Sz

def Splus_Sminus(theta,p,q):
    cos_part_p = (np.cos(theta[p]))
    sin_part_p = (np.sin(theta[p]))
    cos_part_q =  (np.cos(theta[q]))
    sin_part_q = (np.sin(theta[q]))
    prefactor_splus_sminus = cos_part_p*sin_part_p*cos_part_q*sin_part_q
    return prefactor_splus_sminus



def S_zS_z(theta,p,q):
    cos_part_p = (np.cos(theta[p]))**2
    sin_part_p = (np.sin(theta[p]))**2
    cos_part_q =  (np.cos(theta[q]))**2
    sin_part_q = (np.sin(theta[q]))**2  
    numerator = ((-1*cos_part_p + sin_part_p)*(-1*cos_part_q + sin_part_q))
    denominator = ((cos_part_p + sin_part_p)*(cos_part_q+sin_part_q))
    prefactor = (1/4)*((numerator)/(denominator))
    prefactor_SzS_z = prefactor
    return prefactor_SzS_z



def XXZ_1D_overlap(theta,Nsites,Delta):  # function to calculate XXZ_Energy in 1D
    sum_energy = 0
    for i in range(1,Nsites+1):
        x,y = inverse_mapping(i,Ny)  # tuple unpacking

        neighbour_sites = neighbor_square(x,y,Nx,Ny)  # returns neighboring sites for i
        for j in neighbour_sites:
            if i>= j:
                continue

            sum_splus_sminus = Splus_Sminus(theta,i-1,j-1)
            sum_SzSz = S_zS_z(theta,i-1,j-1)
            sum_energy += sum_splus_sminus + Delta*sum_SzSz

    return sum_energy
    

'''cos_part_p = (np.cos(np.pi/2))
sin_part_p = (np.sin(np.pi/2))
cos_part_q =  (np.cos(np.pi/2))
sin_part_q = (np.sin(np.pi/2))
print(cos_part_p)

prefactor_splus_sminus = cos_part_p*sin_part_p*cos_part_q*sin_part_q
print(prefactor_splus_sminus)

cos_part_p = (np.cos(np.pi/4))**2
sin_part_p = (np.sin(np.pi/4))**2
cos_part_q =  (np.cos(np.pi/4))**2
sin_part_q = (np.sin(np.pi/4))**2  
numerator = ((-cos_part_p + sin_part_p)*(-cos_part_q + sin_part_q))
denominator = ((cos_part_p + sin_part_p)*(cos_part_q+sin_part_q))
prefactor = (1/4)*((numerator)/denominator)
print(denominator)'''
    