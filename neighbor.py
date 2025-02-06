from parameter import Nx,Ny

def inverse_mapping(i, Ny):
     
    x = (i - 1) // Ny + 1  # Row index (1-based)
    y = (i - 1) % Ny + 1   # Column index (1-based)
    return (x, y)



def neighbor_square(x,y,Nx,Ny):
    sites = list()

    def map_square(x,y,Nx,Ny):
        return ((x-1)*Nx +y)     # if x represents row, y represents


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
        """Maps (x, y) coordinates to a single index."""
        return (x - 1) * Nx + y
    
    sites = []

    # Right neighbor (x+1, y), wrapping around if x+1 > Nx
    new_x = (x % Nx) + 1  # Wrap around
    sites.append(site_index(new_x, y, Nx, Ny))

    # Left neighbor (x-1, y), wrapping around if x-1 < 1
    new_x = (x - 2) % Nx + 1  # Wrap around
    sites.append(site_index(new_x, y, Nx, Ny))

    # Up neighbor (x, y+1), wrapping around if y+1 > Ny
    new_y = (y % Ny) + 1  # Wrap around
    sites.append(site_index(x, new_y, Nx, Ny))

    # Down neighbor (x, y-1), wrapping around if y-1 < 1
    new_y = (y - 2) % Ny + 1  # Wrap around
    sites.append(site_index(x, new_y, Nx, Ny))

    return sites'''

Nx = 4
Ny = 4
N = Nx*Ny

nearest_neighbour = neighbor_square(3,1,Nx,Ny)
print(nearest_neighbour)