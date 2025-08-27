import numpy as np

def catch_parameter(parameters, target):
    file = open(parameters, "r")
    for line in file.readlines():
        if line[0] != "#":
            splitted = line.split()
            if len(splitted) != 0:
                if splitted[0] == target: 
                    return splitted[2]     

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

def get_analytical_time_reflections(v, z, x):

    Tint = 2.0 * z * v[:-1]
    Vrms = np.zeros(len(z))

    reflections = np.zeros((len(z), len(x)))
    for i in range(len(z)):
        Vrms[i] = np.sqrt(np.sum(v[:i+1]**2 * Tint[:i+1]) / np.sum(Tint[:i+1]))
        reflections[i] = np.sqrt(x**2 + 4.0*np.sum(z[:i+1])**2) / Vrms[i]

    return reflections 

def get_analytical_amps_reflections(vp, ro, z):
    
    reflectivity = np.zeros(len(z))    
    reflectivity = (vp[1:]*ro[1:] - vp[:-1]*ro[:-1]) / (vp[1:]*ro[1:] + vp[:-1]*ro[:-1])

    reflections = np.zeros_like(reflectivity)
    transmission = np.zeros_like(reflectivity)

    reflections[0] = reflectivity[0]
    transmission[0] = 1.0 - reflectivity[0]

    for i in range(1, len(reflectivity)):
        reflections[i] = transmission[i-1]*reflectivity[i]
        transmission[i] = transmission[i-1]*(1.0 - reflectivity[i])

        for j in range(i,0,-1):
            reflections[i] *= 1.0 - reflectivity[i - j]

    return reflections

def get_ricker_wavelet(nt, dt, fmax):
    fc = fmax / (3.0*np.sqrt(np.pi))
    arg = np.pi*((np.arange(nt) - 0.5*nt)*dt*fc*np.pi)**2
    return (1.0 - 2.0*arg)*np.exp(-arg)

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
