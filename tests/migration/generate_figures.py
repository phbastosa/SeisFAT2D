import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

parameters = str(sys.argv[1])

SPS = np.loadtxt(pyf.catch_parameter(parameters,"SPS"), delimiter = ",", dtype = np.float32) 
RPS = np.loadtxt(pyf.catch_parameter(parameters,"RPS"), delimiter = ",", dtype = np.float32) 
XPS = np.loadtxt(pyf.catch_parameter(parameters,"XPS"), delimiter = ",", dtype = np.int32) 

nx = 501
nz = 101 

dh = 10.0

ns = 61
nr = 501

model_vp = pyf.read_binary_matrix(nz, nx, f"../inputs/models/migration_test_vp.bin")
seismic = pyf.read_binary_matrix(nz, nx, f"../outputs/seismic/migration_test_kirchhoff_result_101x501.bin")

xloc = np.linspace(0, nx-1, 11)
xlab = np.array(xloc*dh, dtype = int)

zloc = np.linspace(0, nz-1, 6)
zlab = np.array(zloc*dh, dtype = int)

fig, ax = plt.subplots(figsize = (15, 5))

im = ax.imshow(model_vp, aspect = "auto", cmap = "jet")

ax.imshow(seismic, aspect = "auto", cmap = "Greys", alpha = 0.5)

cbar = plt.colorbar(im)
cbar.set_label("Velocity P [m/s]")

ax.plot(RPS[:, 0]/dh, RPS[:, 1]/dh, "ob")
ax.plot(SPS[:, 0]/dh, SPS[:, 1]/dh, "or")

ax.set_xticks(xloc)
ax.set_yticks(zloc)
ax.set_xticklabels(xlab)    
ax.set_yticklabels(zlab)    
ax.set_ylabel("Depth [km]", fontsize = 15)
ax.set_xlabel("Distance [km]", fontsize = 15)

fig.tight_layout()
plt.savefig("migration_test_setup.png", dpi = 200)
plt.show()

