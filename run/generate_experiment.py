import numpy as np

depth = np.array([1000, 1500, 1900, 2600, 3100, 3700])
vp = np.array([1650, 1800, 2000, 2200, 2500, 3000])

vs = 0.7*vp
rho = 310*vp**0.25

dh = 12.5

x = 2e4
z = 4e3

nx = int(x / dh) + 1
nz = int(z / dh) + 1

model_vp = np.zeros((nz,nx)) + 1500.0
model_vs = np.zeros((nz,nx))
model_rho = np.zeros((nz,nx)) + 1000.0

for i, k in enumerate(depth):
    model_vp[int(k/dh):] = vp[i]
    model_vs[int(k/dh):] = vs[i]
    model_rho[int(k/dh):] = rho[i]

model_vp.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/model_vp_{nz}x{nx}_{dh:.1f}m.bin")
model_vs.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/model_vs_{nz}x{nx}_{dh:.1f}m.bin")
model_rho.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/model_rho_{nz}x{nx}_{dh:.1f}m.bin")

ns = 477
nr = 797

spread = 321

SPS = np.zeros((ns, 2))
SPS[:,0] = np.linspace(0, 11900, ns)
SPS[:,1] = np.zeros(ns) + 25.0

RPS = np.zeros((nr, 2))
RPS[:,0] = np.linspace(100, 20000, nr)
RPS[:,1] = np.zeros(nr) + 25.0

XPS = np.zeros((ns, 3))
XPS[:,0] = np.linspace(ns-1, 0, ns)
XPS[:,1] = (nr - spread) - np.arange(ns) 
XPS[:,2] = nr - np.arange(ns)

np.savetxt("../inputs/geometry/user_SPS.txt", SPS, delimiter = ",", fmt = "%.2f")
np.savetxt("../inputs/geometry/user_RPS.txt", RPS, delimiter = ",", fmt = "%.2f")
np.savetxt("../inputs/geometry/user_XPS.txt", XPS, delimiter = ",", fmt = "%.0f")

