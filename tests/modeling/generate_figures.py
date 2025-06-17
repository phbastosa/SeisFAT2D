import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

SPS = np.loadtxt("../inputs/geometry/modeling_test_SPS.txt", delimiter = ",", comments = "#", dtype = float) 
RPS = np.loadtxt("../inputs/geometry/modeling_test_RPS.txt", delimiter = ",", comments = "#", dtype = float) 
XPS = np.loadtxt("../inputs/geometry/modeling_test_XPS.txt", delimiter = ",", comments = "#", dtype = int) 

m2km = 1e-3

nx = 201
nz = 201 

dh = 50.0

nr = len(RPS)

model_vp = pyf.read_binary_matrix(nz, nx, f"../inputs/models/modeling_test_vp.bin")

xloc = np.linspace(0, nx-1, 11)
xlab = np.array(xloc*dh*m2km, dtype = int)

zloc = np.linspace(0, nz-1, 6)
zlab = np.array(zloc*dh*m2km, dtype = int)

fig, ax = plt.subplots(figsize = (15, 7))

im = ax.imshow(model_vp, aspect = "auto", cmap = "Greys")

cbar = plt.colorbar(im)
cbar.set_label("Velocity P [m/s]")

ax.plot(RPS[:, 0]/dh, RPS[:, 1]/dh, "ob")
ax.plot(SPS[0]/dh, SPS[1]/dh, "or")

ax.set_xticks(xloc)
ax.set_yticks(zloc)
ax.set_xticklabels(xlab)    
ax.set_yticklabels(zlab)    
ax.set_ylabel("Depth [km]", fontsize = 15)
ax.set_xlabel("Distance [km]", fontsize = 15)

fig.tight_layout()
plt.savefig("modeling_test_models.png", dpi = 200)

eikonal_iso = pyf.read_binary_array(nr, f"../outputs/syntheticData/eikonal_iso_nStations{nr}_shot_1.bin")
eikonal_ani = pyf.read_binary_array(nr, f"../outputs/syntheticData/eikonal_ani_nStations{nr}_shot_1.bin")

offset = np.sqrt((SPS[0] - RPS[:,0])**2 + (SPS[1] - RPS[:,1])**2)

analyticalT = offset / model_vp[0,0] 

fig, ax = plt.subplots(figsize = (16,6))
  
ax.plot(eikonal_iso, ".", label = "Eikonal Isotropic")
ax.plot(eikonal_ani, ".", label = r"Eikonal Anisotropic $\epsilon = 0.1$")
ax.plot(analyticalT, ".", label = "Analytical isotropic")

ax.set_ylabel("Time [s]", fontsize = 15)
ax.set_xlabel("Trace index", fontsize = 15)
ax.legend(loc = "upper right", fontsize = 15)

ax.invert_yaxis()
fig.tight_layout()
plt.savefig("modeling_test_data.png", dpi = 300)
