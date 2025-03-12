def bcs_overlap(epsilon,eta):
    overlap = 1
    for k in range(Nx*Ny):
        overlap *= (1+ (epsilon[k].conjugate())*eta[k])
    return overlap

def Sz_p(epsilon,eta,i):
    overlap_Sz = (1/2)*(-1+epsilon[i-1].conjugate()*eta[i-1])/(1 + epsilon[i-1].conjugate()*eta[i-1])

    return overlap_Sz*bcs_overlap(epsilon,eta)

def Sz_pSz_q(epsilon,eta,i,j):

    overlap_Sz_pSz_q = (1/4)*((-1+epsilon[i-1].conjugate()*eta[i-1])*(-1+epsilon[j-1].conjugate()*eta[j-1]))\
        /((1 + epsilon[i-1].conjugate()*eta[i-1])*(1 + epsilon[j-1].conjugate()*eta[j-1]))
    return overlap_Sz_pSz_q*bcs_overlap(epsilon,eta)

def Splus_pSminus_q(epsilon,eta,i,j):   #S+_pS-_q
    overlap_SpSq =  ((epsilon[i-1].conjugate())*(eta[j-1]))\
        /((1 + epsilon[i-1].conjugate()*eta[i-1])*(1 + epsilon[j-1].conjugate()*eta[j-1]))  
    return overlap_SpSq  