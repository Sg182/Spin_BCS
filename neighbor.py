from parameter import Nx,Ny

def inverse_mapping(i, Ny):
     
    x = (i - 1) // Ny + 1  # Row index (1-based)
    y = (i - 1) % Ny + 1   # Column index (1-based)
    return (x, y)



def neighbor_square(x, y, Nx, Ny, periodic=True):
    sites = []

    def map_square(x, y, Nx, Ny):
        return (x - 1) * Ny + y  # Correct: Ny for column-major indexing

    if periodic:
        # Right neighbor (x+1, y)
        new_x = (x % Nx) + 1
        sites.append(map_square(new_x, y, Nx, Ny))

        # Left neighbor (x-1, y)
        new_x = (x - 2) % Nx + 1
        sites.append(map_square(new_x, y, Nx, Ny))

        # Up neighbor (x, y+1)
        new_y = (y % Ny) + 1
        sites.append(map_square(x, new_y, Nx, Ny))

        # Down neighbor (x, y-1)
        new_y = (y - 2) % Ny + 1
        sites.append(map_square(x, new_y, Nx, Ny))

    else:
        # Right neighbor (x+1, y)
        if x + 1 <= Nx:
            sites.append(map_square(x + 1, y, Nx, Ny))

        # Left neighbor (x-1, y)
        if x - 1 >= 1:
            sites.append(map_square(x - 1, y, Nx, Ny))

        # Up neighbor (x, y+1)
        if y + 1 <= Ny:
            sites.append(map_square(x, y + 1, Nx, Ny))

        # Down neighbor (x, y-1)
        if y - 1 >= 1:
            sites.append(map_square(x, y - 1, Nx, Ny))

    return sites

    

     

'''def neighbor_square(x, y, Nx, Ny):
     

    # Nested site_index function
    def site_index(x, y, Nx, Ny):
        
        return (x - 1) * Nx + y
    
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

    return sites'''

def second_neighbor(x, y, Nx, Ny, periodic=True):  
    sites = []

    def map_square(x, y, Nx, Ny):
        return (x - 1) * Ny + y  # row-major indexing

    if periodic:
        if Nx != 1:
            # Diagonal neighbors (x ± 1, y ± 1)
            sites.append(map_square((x + 1 - 1) % Nx + 1, (y + 1 - 1) % Ny + 1, Nx, Ny))
            sites.append(map_square((x + 1 - 1) % Nx + 1, (y - 1 - 1) % Ny + 1, Nx, Ny))
            sites.append(map_square((x - 1 - 1) % Nx + 1, (y + 1 - 1) % Ny + 1, Nx, Ny))
            sites.append(map_square((x - 1 - 1) % Nx + 1, (y - 1 - 1) % Ny + 1, Nx, Ny))
        else:
            sites.append(map_square(x, (y + 2 - 1) % Ny + 1, Nx, Ny))
            sites.append(map_square(x, (y - 2 - 1) % Ny + 1, Nx, Ny))

    else:
        if Nx != 1:
            # Diagonal neighbors (x ± 1, y ± 1), open boundaries
            if x + 1 <= Nx and y + 1 <= Ny:
                sites.append(map_square(x + 1, y + 1, Nx, Ny))
            if x + 1 <= Nx and y - 1 >= 1:
                sites.append(map_square(x + 1, y - 1, Nx, Ny))
            if x - 1 >= 1 and y + 1 <= Ny:
                sites.append(map_square(x - 1, y + 1, Nx, Ny))
            if x - 1 >= 1 and y - 1 >= 1:
                sites.append(map_square(x - 1, y - 1, Nx, Ny))
        else:
            # Next-nearest along column (since Nx = 1)
            if y + 2 <= Ny:
                sites.append(map_square(x, y + 2, Nx, Ny))
            if y - 2 >= 1:
                sites.append(map_square(x, y - 2, Nx, Ny))

    return sites

 
 
def neighbor_triangular(x, y, Nx, Ny,periodic=True):
    """Find nearest neighbors for a triangular lattice (1-based indexing) with PBC."""

    def map_square(x, y, Nx, Ny):
        return (x - 1) * Ny + y

    sites = []

    if periodic:  # PBC

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

    else:   #OBC
        if (x - 1) % 2 == 0:  # even row
        # Right neighbor (x, y+1)
            if y + 1 <= Ny:
                sites.append(map_square(x, y + 1, Nx, Ny))

        # Top-right neighbor (x-1, y+1)
            if x - 1 >= 1 and y + 1 <= Ny:
                sites.append(map_square(x - 1, y + 1, Nx, Ny))

        # Top-left neighbor (x-1, y)
            if x - 1 >= 1:
                sites.append(map_square(x - 1, y, Nx, Ny))

        # Left neighbor (x, y-1)
            if y - 1 >= 1:
                sites.append(map_square(x, y - 1, Nx, Ny))

        # Bottom-left neighbor (x+1, y)
            if x + 1 <= Nx:
                sites.append(map_square(x + 1, y, Nx, Ny))

        # Bottom-right neighbor (x+1, y+1)
            if x + 1 <= Nx and y + 1 <= Ny:
                sites.append(map_square(x + 1, y + 1, Nx, Ny))
        else:  # odd row
        # Right neighbor (x, y+1)
            if y + 1 <= Ny:
                sites.append(map_square(x, y + 1, Nx, Ny))

        # Top-right neighbor (x-1, y)
            if x - 1 >= 1:
                sites.append(map_square(x - 1, y, Nx, Ny))

        # Top-left neighbor (x-1, y-1)
            if x - 1 >= 1 and y - 1 >= 1:
                sites.append(map_square(x - 1, y - 1, Nx, Ny))

        # Left neighbor (x, y-1)
            if y - 1 >= 1:
                sites.append(map_square(x, y - 1, Nx, Ny))

        # Bottom-left neighbor (x+1, y-1)
            if x + 1 <= Nx and y - 1 >= 1:
                sites.append(map_square(x + 1, y - 1, Nx, Ny))

        # Bottom-right neighbor (x+1, y)
            if x + 1 <= Nx:
                sites.append(map_square(x + 1, y, Nx, Ny))

    return sites

def second_neighbor_triangular(x, y, Nx, Ny, periodic=True):
    """Second nearest neighbors for triangular lattice (1-based indexing)."""

    def map_square(x, y, Nx, Ny):
        return (x - 1) * Ny + y

    sites = []

    if periodic:
        if (x - 1) % 2 == 0:  # even row
            neighbor_offsets = [(-2, 0), (-1, -1), (-1, 2), 
                                
                                (1, -1), (1, 2), (2, 0)]
        else:  # odd row
            neighbor_offsets = [(-2, 0), (-1, -2), (-1, 1),
                                 
                                (1, -2), (1, 1), (2, 0)]

        for dx, dy in neighbor_offsets:
            new_x = (x - 1 + dx) % Nx + 1
            new_y = (y - 1 + dy) % Ny + 1
            sites.append(map_square(new_x, new_y, Nx, Ny))
            print(map_square(new_x, new_y, Nx, Ny))
        

    else:  # open boundary conditions
        if (x - 1) % 2 == 0:  # even row
            neighbor_offsets = [(-2, 0), (-1, -1), (-1, 2),
                                 
                                (1, -1), (1, 2), (2, 0)]
        else:  # odd row
            neighbor_offsets = [(-2, 0), (-1, -2), (-1, 1),
                                 
                                (1, -2), (1, 1), (2, 0)]

        for dx, dy in neighbor_offsets:
            new_x = x + dx
            new_y = y + dy
            if 1 <= new_x <= Nx and 1 <= new_y <= Ny:
                sites.append(map_square(new_x, new_y, Nx, Ny))

    return sites



     

neig_pbc  = neighbor_triangular(3,3,4,4,periodic=True)
neig_obc  = second_neighbor_triangular(3,3,4,4,periodic=True)
print(neig_pbc)
print(neig_obc)

 