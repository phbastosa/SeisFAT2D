import sys; sys.path.append("../src/")

import numpy as np

import functions as pyf

x_max = 2e4
z_max = 5e3

dh = 10.0

nx = int((x_max / dh) + 1)
nz = int((z_max / dh) + 1)

vp = np.array([1500, 1700, 1900, 2300, 3000, 3500])
vs = np.array([   0, 1100, 1250, 1400, 1800, 2100])
ro = np.array([1000, 2250, 2280, 2300, 2370, 2450])

z = np.array([200, 500, 1000, 1500, 1500])

epsilon = np.array([0, 0, 0, 0, 0, 0, 0, 0])
delta = np.array([0, 0, 0, 0, 0, 0, 0, 0])
tilt = np.array([0, 0, 0, 0, 0, 0, 0, 0])

c11 = np.zeros_like(vp)
c13 = np.zeros_like(vp)
c15 = np.zeros_like(vp)
c33 = np.zeros_like(vp)
c35 = np.zeros_like(vp)
c55 = np.zeros_like(vp)

for k in range(len(vp)):
    
    C = pyf.compute_stiffness(vp[k], vs[k], ro[k], epsilon[k], delta[k], tilt[k])
    
    c11[k] = C[0,0] 
    c13[k] = C[0,1] 
    c15[k] = C[0,2] 
    c33[k] = C[1,1]
    c35[k] = C[1,2]
    c55[k] = C[2,2]

Vp = np.zeros((nz,nx)) + 1500

C11 = np.zeros((nz,nx))
C13 = np.zeros((nz,nx))
C15 = np.zeros((nz,nx))
C33 = np.zeros((nz,nx))
C35 = np.zeros((nz,nx))
C55 = np.zeros((nz,nx))

for i in range(len(z)): 

    Vp[int(np.sum(z[:i+1]/dh)):] = vp[i+1]

    C11[int(np.sum(z[:i+1]/dh)):] = c11[i+1]
    C13[int(np.sum(z[:i+1]/dh)):] = c13[i+1]
    C15[int(np.sum(z[:i+1]/dh)):] = c15[i+1]
    C33[int(np.sum(z[:i+1]/dh)):] = c33[i+1]
    C35[int(np.sum(z[:i+1]/dh)):] = c35[i+1]
    C55[int(np.sum(z[:i+1]/dh)):] = c55[i+1]

Vp.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/modeling_test_vp.bin")

C11.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/modeling_test_C11.bin")
C13.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/modeling_test_C13.bin")
C15.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/modeling_test_C15.bin")
C33.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/modeling_test_C33.bin")
C35.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/modeling_test_C35.bin")
C55.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/modeling_test_C55.bin")
