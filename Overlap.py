from neighbor import * # imports Nx,Ny


def overlap_XXZ(Delta,epsilon,eta):
    total_overlap = 0

   #nearest_neighbours = neighbor_square(x,y,Nx,Ny)

    for i in range(1,Nx*Ny+1):

        x,y = inverse_mapping(i,Ny)  # tuple unpacking

        neighbour_sites = neighbor_square(x,y,Nx,Ny)

        for j in neighbour_sites:

            if i>= j:
                continue
                 # ensuring i < j to avoid double counting
            prefactor_SzSz = (1/4)*((-1+epsilon[i-1].conjugate()*eta[i-1])*(-1+epsilon[j-1].conjugate()*eta[j-1]) )/((1 + epsilon[i-1].conjugate()*eta[i-1])*(1 + epsilon[j-1].conjugate()*eta[j-1])) # For Delta*S^z_i*S^z_j

            prefactor_SxSx = (1/4)*((epsilon[i-1].conjugate()+eta[i-1])*(epsilon[j-1].conjugate()+eta[j-1]))/((1 + epsilon[i-1].conjugate()*eta[i-1])*(1 + epsilon[j-1].conjugate()*eta[j-1]))

            prefactor_SySy = (-1/4)*((-1*epsilon[i-1].conjugate()+eta[i-1])*(-1*epsilon[j-1].conjugate()+eta[j-1]))/((1 + epsilon[i-1].conjugate()*eta[i-1])*(1 + epsilon[j-1].conjugate()*eta[j-1]))

            total_overlap += prefactor_SxSx + prefactor_SySy + Delta*prefactor_SzSz # adding all the contributions due to summation in the Hamiltonian
            #print(total_overlap)
    # BCS overlap

    bcs_overlap = 1
    
    for k in range(Nx*Ny):
        bcs_overlap *= ( 1+ epsilon[k].conjugate()*eta[k])

    return total_overlap*bcs_overlap

def bcs_overlap(epsilon,eta):   # computes overlap between two BCS states
    overlap = 1
    for k in range(Nx*Ny):
        overlap *= ( 1+ (epsilon[k].conjugate())*eta[k])
    return overlap

def Sz(epsilon,eta,i):     # calculates Sz for site i
    prefactor_Sz = (1/2)*(-1+epsilon[i-1].conjugate()*eta[i-1])/(1 + epsilon[i-1].conjugate()*eta[i-1]) 
    bcs_overlap = 1
    for k in range(Nx*Ny):
        bcs_overlap *= ( 1+ epsilon[k].conjugate()*eta[k])
    return bcs_overlap*prefactor_Sz



                             