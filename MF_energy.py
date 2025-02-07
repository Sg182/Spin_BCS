from Overlap import overlap_XXZ, bcs_overlap
from scipy.optimize import minimize
from parameter import *

def energy_bcs(Delta,eta):
    return (overlap_XXZ(Delta,eta,eta))/(bcs_overlap(eta,eta))

# initial guess for eta

eta_initial = [complex(1,0.5) for i in range(Nx*Ny)]

result = minimize(energy_bcs, eta_initial, args=(Delta), method='BFGS')

eta_optimized = result.x
final_energy = result.fun

print(f"BCS Energy: {final_energy}")


