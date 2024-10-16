import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import pyFunctions as pyf

SPS = np.loadtxt("../inputs/geometry/modeling_test_SPS.txt", delimiter = ",", comments = "#", dtype = float) 
RPS = np.loadtxt("../inputs/geometry/modeling_test_RPS.txt", delimiter = ",", comments = "#", dtype = float) 
XPS = np.loadtxt("../inputs/geometry/modeling_test_XPS.txt", delimiter = ",", comments = "#", dtype = int) 

m2km = 1e-3

nx = 2001
nz = 501 

dh = 10.0

model_vp = pyf.read_binary_matrix(nz, nx, f"../inputs/models/modeling_test_vp_model_{nz}x{nx}_{dh:.0f}m.bin")
model_vs = pyf.read_binary_matrix(nz, nx, f"../inputs/models/modeling_test_vs_model_{nz}x{nx}_{dh:.0f}m.bin")
model_rho = pyf.read_binary_matrix(nz, nx, f"../inputs/models/modeling_test_rho_model_{nz}x{nx}_{dh:.0f}m.bin")

xloc = np.linspace(0, nx-1, 11)
xlab = np.array(xloc*dh*m2km, dtype = int)

zloc = np.linspace(0, nz-1, 5)
zlab = np.array(zloc*dh*m2km, dtype = int)

fig, ax = plt.subplots(figsize = (15, 7), ncols = 1, nrows = 3)

im = ax[0].imshow(model_vp, aspect = "auto", cmap = "Greys")

cbar = plt.colorbar(im)
cbar.set_label("Velocity P [m/s]")

ax[0].plot(RPS[:, 0]/dh, RPS[:, 1]/dh, "ob")
ax[0].plot(SPS[:, 0]/dh, SPS[:, 1]/dh, "or")

ax[0].set_xticks(xloc)
ax[0].set_yticks(zloc)
ax[0].set_xticklabels(xlab)    
ax[0].set_yticklabels(zlab)    
ax[0].set_ylabel("Depth [km]", fontsize = 15)
ax[0].set_xlabel("Distance [km]", fontsize = 15)

im = ax[1].imshow(model_vs, aspect = "auto", cmap = "Greys")

cbar = plt.colorbar(im)
cbar.set_label("Velocity S [m/s]")

ax[1].plot(RPS[:, 0]/dh, RPS[:, 1]/dh, "ob")
ax[1].plot(SPS[:, 0]/dh, SPS[:, 1]/dh, "or")

ax[1].set_xticks(xloc)
ax[1].set_yticks(zloc)
ax[1].set_xticklabels(xlab)    
ax[1].set_yticklabels(zlab)    
ax[1].set_ylabel("Depth [km]", fontsize = 15)
ax[1].set_xlabel("Distance [km]", fontsize = 15)

im = ax[2].imshow(model_rho, aspect = "auto", cmap = "Greys")

cbar = plt.colorbar(im)
cbar.set_label("Density [kg/m³]")

ax[2].plot(RPS[:, 0]/dh, RPS[:, 1]/dh, "ob")
ax[2].plot(SPS[:, 0]/dh, SPS[:, 1]/dh, "or")

ax[2].set_xticks(xloc)
ax[2].set_yticks(zloc)
ax[2].set_xticklabels(xlab)    
ax[2].set_yticklabels(zlab)    
ax[2].set_ylabel("Depth [km]", fontsize = 15)
ax[2].set_xlabel("Distance [km]", fontsize = 15)

fig.tight_layout()
plt.savefig("modeling_test_models.png", dpi = 200)
plt.show()

nt = 10001
dt = 1e-3

nShots = 3
nStations = 401

offset = np.arange(nStations)

fig, ax = plt.subplots(figsize = (15, 7), ncols = 3, nrows = 1)

xloc = np.linspace(0, nStations-1, 5)
xlab = np.linspace(50, 19950, 5, dtype = int)

tloc = np.linspace(0, nt-1, 11)
tlab = np.linspace(0, (nt-1)*dt, 11, dtype = int)

for i in range(nShots):
    
    eikonal = pyf.read_binary_array(nStations, f"../outputs/syntheticData/eikonal_iso_nStations401_shot_{i+1}.bin")
    elastic = pyf.read_binary_matrix(nt, nStations, f"../outputs/syntheticData/elastic_iso_data_nStations401_nSamples10001_shot_{i+1}.bin")

    scale = np.std(elastic)

    im = ax[i].imshow(elastic, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale)

    ax[i].plot(offset, eikonal / dt, "--r")

    ax[i].set_xticks(xloc)
    ax[i].set_yticks(tloc)
    ax[i].set_xticklabels(xlab)    
    ax[i].set_yticklabels(tlab)    
    ax[i].set_ylabel("Time [s]", fontsize = 15)
    ax[i].set_xlabel("Distance [m]", fontsize = 15)

fig.tight_layout()
plt.savefig("modeling_test_data.png", dpi = 200)
plt.show()




