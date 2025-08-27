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

def compute_stiffness(vp, vs, ro, ep, dl, tht):
    
    SI = 1e9

    C = np.zeros((3,3))
    M = np.zeros((3,3))

    c11 = 0; c13 = 0; c15 = 0
    c33 = 0; c35 = 0; c55 = 0

    SI = 1e9

    c33 = ro*vp**2 / SI
    c55 = ro*vs**2 / SI

    c11 = c33*(1.0 + 2.0*ep)

    c13 = np.sqrt((c33 - c55)**2 + 2.0*dl*c33*(c33 - c55)) - c55

    C[0,0] = c11; C[0,1] = c13; C[0,2] = c15  
    C[1,0] = c13; C[1,1] = c33; C[1,2] = c35  
    C[2,0] = c15; C[2,1] = c35; C[2,2] = c55     

    tht = np.radians(tht)

    c = np.cos(tht)
    s = np.sin(tht)

    sin2 = np.sin(2.0*tht)
    cos2 = np.cos(2.0*tht)

    M = np.array([[     c**2,     s**2, sin2],
                  [     s**2,     c**2,-sin2],
                  [-0.5*sin2, 0.5*sin2, cos2]])
    
    Cr = (M @ C @ M.T) * SI

    return Cr
