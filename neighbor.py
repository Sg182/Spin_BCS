from parameter import Nx,Ny

def inverse_mapping(i, Ny):
     
    x = (i - 1) // Ny + 1  # Row index (1-based)
    y = (i - 1) % Ny + 1   # Column index (1-based)
    return (x, y)



def neighbor_square(x,y,Nx,Ny): ## Periodic Boundary cond.
    sites = list()

    def map_square(x,y,Nx,Ny):
        return ((x-1)*Ny +y)     # if x represents row, y represents


    if (x +1 >= Nx):
        sites.append(map_square(Nx,y,Nx,Ny))

    else:
        sites.append(map_square(x+1,y,Nx,Ny))

    if (x -1 <=0):
        sites.append(map_square(Nx,y,Nx,Ny))

    else:
        sites.append(map_square(x-1,y,Nx,Ny))

    if (y +1 >= Ny):
        sites.append(map_square(x,Ny,Nx,Ny))

    else:
        sites.append(map_square(x,y+1,Nx,Ny))

    if (y -1 <=0):
        sites.append(map_square(x,Ny,Nx,Ny))

    else:
        sites.append(map_square(x,y-1,Nx,Ny))

    return sites

'''def neighbor_square(x, y, Nx, Ny):
     

    # Nested site_index function
    def site_index(x, y, Nx, Ny):
        
        return (x - 1) * Ny + y # x represents rows, y represents column
    
    sites = []

     
    new_x = (x % Nx) + 1  # Wrap around
    sites.append(site_index(new_x, y, Nx, Ny))

    # Left neighbor (x-1, y), wrapping around if x-1 < 1
    new_x = (x - 2) % Nx + 1  # Wrap around
    sites.append(site_index(new_x, y, Nx, Ny))

     
    new_y = (y % Ny) + 1  # Wrap around
    sites.append(site_index(x, new_y, Nx, Ny))

    # Down neighbor (x, y-1), wrapping around if y-1 < 1
    new_y = (y - 2) % Ny + 1  # Wrap around
    sites.append(site_index(x, new_y, Nx, Ny))

    return list(sites)'''

def second_neighbor(x,y,Nx,Ny):  ## PBC # x rep rows, y columns    ## ((x+-1)-1)%Nx + 1
    sites = []

    def map_square(x,y,Nx,Ny):
        return ((x-1)*Ny +y)
    if Nx != 1:

    # Diagonal neighbors (x ± 1, y ± 1)
        sites.append(map_square((x + 1 - 1) % Nx + 1, (y + 1 - 1) % Ny + 1, Nx, Ny))
        sites.append(map_square((x + 1 - 1) % Nx + 1, (y - 1 - 1) % Ny + 1, Nx, Ny))
        sites.append(map_square((x - 1 - 1) % Nx + 1, (y + 1 - 1) % Ny + 1, Nx, Ny))
        sites.append(map_square((x - 1 - 1) % Nx + 1, (y - 1 - 1) % Ny + 1, Nx, Ny))

    else:
    # Next-nearest along the same row or column (x ± 2, y ± 2)
    ##sites.add(map_square((x + 2 - 1) % Nx + 1, y, Nx, Ny))
    ##sites.add(map_square((x - 2 - 1) % Nx + 1, y, Nx, Ny))
        sites.append(map_square(x, (y + 2 - 1) % Ny + 1, Nx, Ny))
        sites.append(map_square(x, (y - 2 - 1) % Ny + 1, Nx, Ny))

    return list(sites)

 
 
def neighbor_triangular(x, y, Nx, Ny):
    """Find nearest neighbors for a triangular lattice (1-based indexing) with PBC."""

    def map_square(x, y, Nx, Ny):
        return (x - 1) * Ny + y

    sites = []

    if (x-1) % 2 == 0:  # even row
        
        sites.append(map_square(x, (y + 1 - 1) % Ny + 1, Nx, Ny))# right
        
        sites.append(map_square((x - 1 - 1) % Nx + 1, (y + 1 - 1) % Ny + 1, Nx, Ny))# top-right
    
        sites.append(map_square((x - 1 - 1) % Nx + 1, y, Nx, Ny))# top-left
        
        sites.append(map_square(x, (y - 1 - 1) % Ny + 1, Nx, Ny))# left
        
        sites.append(map_square((x + 1 - 1) % Nx + 1, y, Nx, Ny))# bottom-left
        
        sites.append(map_square((x + 1 - 1) % Nx + 1, (y + 1 - 1) % Ny + 1, Nx, Ny))# bottom-right
    else:  # odd row
        
        sites.append(map_square(x, (y + 1 - 1) % Ny + 1, Nx, Ny))# right
        
        sites.append(map_square((x - 1 - 1) % Nx + 1, y, Nx, Ny))# top-right
        
        sites.append(map_square((x - 1 - 1) % Nx + 1, (y - 1 - 1) % Ny + 1, Nx, Ny))# top-left
        
        sites.append(map_square(x, (y - 1 - 1) % Ny + 1, Nx, Ny))# left
        
        sites.append(map_square((x + 1 - 1) % Nx + 1, (y - 1 - 1) % Ny + 1, Nx, Ny)) # bottom-left
 
        sites.append(map_square((x + 1 - 1) % Nx + 1, y, Nx, Ny)) # bottom-right

    return sites

#neig  = neighbor_triangular(2,4,4,4)
#print(neig)

neighsq = neighbor_square(1,2,2,4)
print(neighsq)