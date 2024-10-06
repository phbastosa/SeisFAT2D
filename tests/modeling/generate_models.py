import numpy as np

x_max = 2e4
z_max = 5e3

dh = 10.0

nx = int((x_max / dh) + 1)
nz = int((z_max / dh) + 1)

nInterfaces = 8
depth_beg = 1000.0
depth_end = 4500.0

velp_beg = 1650.0
velp_end = 3000.0

z = np.linspace(depth_beg, depth_end, nInterfaces)
vp = np.linspace(velp_beg, velp_end, nInterfaces)
vs = 0.7 * vp 
rho = 310*vp**0.25

model_vp = np.zeros((nz, nx)) + 1500.0
model_vs = np.zeros((nz, nx))
model_rho = np.zeros((nz, nx)) + 1000.0

for i, depth in enumerate(z):
    model_vp[int(depth/dh):] = vp[i] 
    model_vs[int(depth/dh):] = vs[i] 
    model_rho[int(depth/dh):] = rho[i] 

model_vp.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/modeling_test_vp_{nz}x{nx}_{dh:.0f}m.bin")
model_vs.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/modeling_test_vs_{nz}x{nx}_{dh:.0f}m.bin")
model_rho.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/modeling_test_rho_{nz}x{nx}_{dh:.0f}m.bin")
