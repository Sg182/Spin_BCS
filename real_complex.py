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
    if len(arr) % 2 != 0:
        raise ValueError("The number of elemets in the array must be even")
    arr1 = np.column_stack((arr.real,arr.imag)) # see np.column_stack((arr1,arr2)) method
    arr_real = arr1.ravel()   # array.ravel() flattens the array,

    return arr_real

#arr = real_to_complex([1,2,3,4,8,9])
#print(complex_to_real(arr))

