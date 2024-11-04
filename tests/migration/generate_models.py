import numpy as np

x_max = 1e4
z_max = 1e3

dh = 5.0

nx = int((x_max / dh) + 1)
nz = int((z_max / dh) + 1)

amps = np.array([2.0, 3.0, 2.0])
center_x = np.array([0.25*x_max, 0.5*x_max, 0.75*x_max])
sigma_x = np.array([0.02*x_max, 0.05*x_max, 0.02*x_max])

x = np.arange(nx)*dh
z = np.zeros_like(x)
for i in range(len(sigma_x)):
    z += -0.4*z_max*(0.5*amps[i]*np.exp(-0.5*((x - center_x[i])/sigma_x[i])**2) - 0.75)

vp_model = np.zeros((nz, nx)) + 2000

for i in range(nx):
    vp_model[int(z[i]/dh):, i] = 2500

vs_model = 0.7*vp_model
rho_model = 310*vp_model**0.25

vp_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/migration_test_vp_model_{nz}x{nx}_{dh:.0f}m.bin")
vs_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/migration_test_vs_model_{nz}x{nx}_{dh:.0f}m.bin")
rho_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/migration_test_rho_model_{nz}x{nx}_{dh:.0f}m.bin")
