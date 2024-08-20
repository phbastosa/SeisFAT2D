import numpy as np

def read_binary_array(n1,filename):
    return np.fromfile(filename, dtype = np.float32, count = n1)    

def read_binary_matrix(n1,n2,filename):
    data = np.fromfile(filename, dtype = np.float32, count = n1*n2)   
    return np.reshape(data, [n1, n2], order='F')

def get_analytical_refractions(v, z, x):

    refracted_waves = np.zeros((len(z), len(x)))

    for n in range(len(z)):
        refracted_waves[n,:] = x / v[n+1]
        for i in range(n+1):
            angle = np.arcsin(v[i] / v[n+1])
            refracted_waves[n,:] += 2.0*z[i]*np.cos(angle) / v[i]
    
    return refracted_waves