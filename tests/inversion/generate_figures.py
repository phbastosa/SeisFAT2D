import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

parameters = str(sys.argv[1])

SPS = np.loadtxt(pyf.catch_parameter(parameters, "SPS"), delimiter = ",", comments = "#", dtype = float) 
RPS = np.loadtxt(pyf.catch_parameter(parameters, "RPS"), delimiter = ",", comments = "#", dtype = float) 
XPS = np.loadtxt(pyf.catch_parameter(parameters, "XPS"), delimiter = ",", comments = "#", dtype = int) 

nx = int(pyf.catch_parameter(parameters, "x_samples"))
nz = int(pyf.catch_parameter(parameters, "z_samples"))

dx = float(pyf.catch_parameter(parameters, "x_spacing"))
dz = float(pyf.catch_parameter(parameters, "z_spacing"))

init_vp = pyf.read_binary_matrix(nz, nx, f"../inputs/models/inversion_test_init_vp.bin")
true_vp = pyf.read_binary_matrix(nz, nx, f"../inputs/models/inversion_test_true_vp.bin")
tomo_vp = pyf.read_binary_matrix(nz, nx, f"../outputs/models/inversion_test_tomography_iso_final_model_vp_{nz}x{nx}.bin")

m2km = 1e-3

xloc = np.linspace(0, nx-1, 6)
xlab = np.around(xloc*dx*m2km, decimals = 1)

zloc = np.linspace(0, nz-1, 6)
zlab = np.around(zloc*dz*m2km, decimals = 1)

fig, ax = plt.subplots(figsize = (15, 4), ncols = 3)

im = ax[0].imshow(true_vp, aspect = "auto", cmap = "jet", vmin = 1500, vmax = 2000)

cbar = plt.colorbar(im)
cbar.set_label("Velocity P [m/s]")

ax[0].plot(RPS[:, 0]/dx, RPS[:, 1]/dz, "ob")
ax[0].plot(SPS[:, 0]/dx, SPS[:, 1]/dz, "or")

ax[0].set_xticks(xloc)
ax[0].set_yticks(zloc)
ax[0].set_xticklabels(xlab)    
ax[0].set_yticklabels(zlab)    
ax[0].set_ylabel("Depth [km]", fontsize = 15)
ax[0].set_xlabel("Distance [km]", fontsize = 15)

im = ax[1].imshow(init_vp, aspect = "auto", cmap = "jet", vmin = 1500, vmax = 2000)

cbar = plt.colorbar(im)
cbar.set_label("Velocity P [m/s]")

ax[1].plot(RPS[:, 0]/dx, RPS[:, 1]/dz, "ob")
ax[1].plot(SPS[:, 0]/dx, SPS[:, 1]/dz, "or")

ax[1].set_xticks(xloc)
ax[1].set_yticks(zloc)
ax[1].set_xticklabels(xlab)    
ax[1].set_yticklabels(zlab)    
ax[1].set_ylabel("Depth [km]", fontsize = 15)
ax[1].set_xlabel("Distance [km]", fontsize = 15)

im = ax[2].imshow(tomo_vp, aspect = "auto", cmap = "jet", vmin = 1500, vmax = 2000)

cbar = plt.colorbar(im)
cbar.set_label("Velocity P [m/s]")

ax[2].plot(RPS[:, 0]/dx, RPS[:, 1]/dz, "ob")
ax[2].plot(SPS[:, 0]/dx, SPS[:, 1]/dz, "or")

ax[2].set_xticks(xloc)
ax[2].set_yticks(zloc)
ax[2].set_xticklabels(xlab)    
ax[2].set_yticklabels(zlab)    
ax[2].set_ylabel("Depth [km]", fontsize = 15)
ax[2].set_xlabel("Distance [km]", fontsize = 15)

fig.tight_layout()
plt.savefig("inversion_test_results.png", dpi = 200)
plt.show()