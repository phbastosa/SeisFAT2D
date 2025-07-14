import numpy as np
import matplotlib.pyplot as plt

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

    return C

def createGaussian(nx,dx,A,xc,sigx):
    x = np.arange(nx)*dx
    surface = A*np.exp(-0.5*(x - xc)**2/sigx**2)
    return surface

x_max = 2e4
z_max = 5e3

dh = 100.0

nx = int((x_max / dh) + 1)
nz = int((z_max / dh) + 1)

ns = 181
nr = 401

SPS = np.zeros((ns, 2))
RPS = np.zeros((nr, 2))
XPS = np.zeros((ns, 3))

SPS[:, 0] = np.linspace(1000, 19000, ns) 
SPS[:, 1] = 1000.0 

RPS[:, 0] = np.linspace(0, 20000, nr)
RPS[:, 1] = 10.0 

XPS[:, 0] = np.arange(ns)
XPS[:, 1] = np.zeros(ns) 
XPS[:, 2] = np.zeros(ns) + nr

np.savetxt("../inputs/geometry/anisoTomo_SPS.txt", SPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/anisoTomo_RPS.txt", RPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/anisoTomo_XPS.txt", XPS, fmt = "%.0f", delimiter = ",")

vp = 1650
vs = 1050
ro = 2250

Vp = np.zeros((nz,nx)) + 1500
Vs = np.zeros((nz,nx)) 
Ro = np.zeros((nz,nx)) + 1000

Ep = np.zeros((nz,nx))
Dl = np.zeros((nz,nx))
Tht = np.zeros((nz,nx))

dgvp = 50.0
dgvs = 50.0
dgro = 10.0
wb = 1000.0

for i in range(nz):
    if i >= wb/dh:    
        Vp[i] = vp + (i*dh - wb)*dgvp/dh 
        Vs[i] = vs + (i*dh - wb)*dgvs/dh
        Ro[i] = ro + (i*dh - wb)*dgro/dh

C11 = np.zeros_like(Vp) 
C13 = np.zeros_like(Vp) 
C15 = np.zeros_like(Vp) 
C33 = np.zeros_like(Vp)
C35 = np.zeros_like(Vp)
C55 = np.zeros_like(Vp)

for i in range(nz):
    for j in range(nx):

        C = compute_stiffness(Vp[i,j], Vs[i,j], Ro[i,j], Ep[i,j], Dl[i,j], Tht[i,j])
        
        C11[i,j] = C[0,0] 
        C13[i,j] = C[0,1] 
        C15[i,j] = C[0,2] 
        C33[i,j] = C[1,1]
        C35[i,j] = C[1,2]
        C55[i,j] = C[2,2]

Vp.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/init_vp.bin")

C11.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/init_C11.bin")
C13.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/init_C13.bin")
C15.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/init_C15.bin")
C33.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/init_C33.bin")
C35.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/init_C35.bin")
C55.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/init_C55.bin")

plt.figure(figsize = (10,8))
plt.subplot(211)
plt.imshow(Vp, cmap = "jet", extent = [0, (nx-1)*dh, (nz-1)*dh, 0])
plt.plot(RPS[:,0], RPS[:,1], "o", color = "green")
plt.plot(SPS[:,0], SPS[:,1], "o", color = "gray")

plt.subplot(212)
plt.imshow(Ep, cmap = "jet", extent = [0, (nx-1)*dh, (nz-1)*dh, 0])
plt.plot(RPS[:,0], RPS[:,1], "o", color = "green")
plt.plot(SPS[:,0], SPS[:,1], "o", color = "gray")

plt.tight_layout()
plt.show()

radius = 1000

dvp = np.array([500, -500])
dvs = np.array([300, -300])
dro = np.array([ 50, -50])

dep = np.array([ 0.0, 0.0])
ddl = np.array([ 0.0, 0.0])

circle_centers = np.array([[3000, 8000],
                           [3000, 12000]])

x, z = np.meshgrid(np.arange(nx)*dh, np.arange(nz)*dh)

for k in range(len(circle_centers)):
    
    distance = np.sqrt((x - circle_centers[k,1])**2 + (z - circle_centers[k,0])**2)

    Vp[distance <= radius] += dvp[k]
    Vs[distance <= radius] += dvs[k]
    Ro[distance <= radius] += dro[k]

    Ep[distance <= radius] += dep[k]
    Dl[distance <= radius] += ddl[k]

for i in range(nz):
    for j in range(nx):

        C = compute_stiffness(Vp[i,j], Vs[i,j], Ro[i,j], Ep[i,j], Dl[i,j], Tht[i,j])
        
        C11[i,j] = C[0,0] 
        C13[i,j] = C[0,1] 
        C15[i,j] = C[0,2] 
        C33[i,j] = C[1,1]
        C35[i,j] = C[1,2]
        C55[i,j] = C[2,2]

Vp.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/true_vp.bin")

C11.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/true_C11.bin")
C13.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/true_C13.bin")
C15.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/true_C15.bin")
C33.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/true_C33.bin")
C35.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/true_C35.bin")
C55.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/true_C55.bin")

plt.figure(figsize = (10,8))
plt.subplot(211)
plt.imshow(Vp, cmap = "jet", extent = [0, (nx-1)*dh, (nz-1)*dh, 0])
plt.plot(RPS[:,0], RPS[:,1], "o", color = "green")
plt.plot(SPS[:,0], SPS[:,1], "o", color = "gray")

plt.subplot(212)
plt.imshow(Ep, cmap = "jet", extent = [0, (nx-1)*dh, (nz-1)*dh, 0])
plt.plot(RPS[:,0], RPS[:,1], "o", color = "green")
plt.plot(SPS[:,0], SPS[:,1], "o", color = "gray")

plt.tight_layout()
plt.show()
