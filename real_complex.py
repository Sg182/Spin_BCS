import numpy as np

def real_to_complex(arr):

    '''from 2n real numbers, return n complex number'''
    
    arr = np.asarray(arr) # turns a list into numpy array
    

    if len(arr) % 2 != 0:
        raise ValueError("The number of elements in the array must be even")

    return arr[::2] + 1j*arr[1::2]

def complex_to_real(arr):

    ''' intakes n complex numbers and returns 2n real numbers'''
    arr = np.asarray(arr)
     
    arr1 = np.column_stack((arr.real,arr.imag)) # see np.column_stack((arr1,arr2)) method
    arr_real = arr1.ravel()   # array.ravel() flattens the array,

    return arr_real
    

x = np.array([1.0, 2.0, 3.0, 4.0,1.0,6])
z = real_to_complex(x)         # [1.+2.j 3.+4.j]
x_back = complex_to_real(z)    # [1. 2. 3. 4.]

print("z:", z)
print("x_back:", x_back)


