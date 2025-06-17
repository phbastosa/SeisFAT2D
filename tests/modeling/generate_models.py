import numpy as np

nx = 201
nz = 201

dx = 50.0
dz = 50.0

vp = np.array([1800])
vs = np.array([1060])
ro = np.array([2250])
z = np.array([])

E = np.array([0.1])
D = np.array([0.0])

tht = np.array([0.0]) * np.pi/180.0

V = np.zeros((nz, nx))

C11 = np.zeros_like(V)
C13 = np.zeros_like(V)
C15 = np.zeros_like(V)
C33 = np.zeros_like(V)
C35 = np.zeros_like(V)
C55 = np.zeros_like(V)

C = np.zeros((3,3))
M = np.zeros((3,3))

c11 = c13 = c15 = 0
c33 = c35 = 0
c55 = 0

SI = 1e9

for i in range(len(vp)):
    
    layer = int(np.sum(z[:i])/dz)

    c33 = ro[i]*vp[i]**2 / SI
    c55 = ro[i]*vs[i]**2 / SI

    c11 = c33*(1.0 + 2.0*E[i])

    c13 = np.sqrt((c33 - c55)**2 + 2.0*D[i]*c33*(c33 - c55)) - c55

    C[0,0] = c11; C[0,1] = c13; C[0,2] = c15;  
    C[1,0] = c13; C[1,1] = c33; C[1,2] = c35  
    C[2,0] = c15; C[2,1] = c35; C[2,2] = c55; 

    c = np.cos(tht[i])
    s = np.sin(tht[i])

    sin2 = np.sin(2.0*tht[i])
    cos2 = np.cos(2.0*tht[i])

    M = np.array([[     c**2,     s**2, sin2],
                  [     s**2,     c**2,-sin2],
                  [-0.5*sin2, 0.5*sin2, cos2]])

    Cr = (M @ C @ M.T) * SI

    V[layer:] = vp[i]

    C11[layer:] = Cr[0,0]; C13[layer:] = Cr[0,1]; C15[layer:] = Cr[0,2] 
    C33[layer:] = Cr[1,1]; C35[layer:] = Cr[1,2]; C55[layer:] = Cr[2,2]

V.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/modeling_test_vp.bin")

C11.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/modeling_test_C11.bin")
C13.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/modeling_test_C13.bin")
C15.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/modeling_test_C15.bin")
C33.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/modeling_test_C33.bin")
C35.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/modeling_test_C35.bin")
C55.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/modeling_test_C55.bin")
