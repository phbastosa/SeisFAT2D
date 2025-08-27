import numpy as np

x_max = 2e4
z_max = 5e3

dh = 10.0

nx = int((x_max / dh) + 1)
nz = int((z_max / dh) + 1)

vp_model = np.zeros((nz, nx)) + 1500
vs_model = np.zeros((nz, nx)) 
rho_model = np.zeros((nz, nx)) + 1000

v = np.array([1500, 1700, 1900, 2300, 3000, 3500])
z = np.array([200, 500, 1000, 1500, 1500])

for i in range(len(z)):
    vp_model[int(np.sum(z[:i+1]/dh)):] = v[i+1]
    vs_model[int(np.sum(z[:i+1]/dh)):] = 0.7*v[i+1]
    rho_model[int(np.sum(z[:i+1]/dh)):] = 310*v[i+1]**0.25

vp_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/modeling_test_vp_model_{nz}x{nx}_{dh:.0f}m.bin")
vs_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/modeling_test_vs_model_{nz}x{nx}_{dh:.0f}m.bin")
rho_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/modeling_test_rho_model_{nz}x{nx}_{dh:.0f}m.bin")