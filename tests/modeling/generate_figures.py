import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

SPS = np.loadtxt("../inputs/geometry/modeling_test_SPS.txt", delimiter = ",", comments = "#", dtype = float) 
RPS = np.loadtxt("../inputs/geometry/modeling_test_RPS.txt", delimiter = ",", comments = "#", dtype = float) 
XPS = np.loadtxt("../inputs/geometry/modeling_test_XPS.txt", delimiter = ",", comments = "#", dtype = int) 

m2km = 1e-3

nx = 2001
nz = 501 

dh = 10.0

nt = 10001
dt = 1e-3

ns = 3
nr = 401

model_vp = pyf.read_binary_matrix(nz, nx, f"../inputs/models/modeling_test_vp.bin")

xloc = np.linspace(0, nx-1, 11)
xlab = np.array(xloc*dh*m2km, dtype = int)

zloc = np.linspace(0, nz-1, 6)
zlab = np.array(zloc*dh*m2km, dtype = int)

fig, ax = plt.subplots(figsize = (15, 5))

im = ax.imshow(model_vp, aspect = "auto", cmap = "Greys")

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
plt.savefig("modeling_test_setup.png", dpi = 200)
plt.show()

v = np.array([1500, 1700, 1900, 2300, 3000, 3500])
z = np.array([200, 500, 1000, 1500, 1500])

fig, ax = plt.subplots(figsize = (15, 7), nrows = 3)

eikonal_an = np.zeros(nr)

for i in range(ns):

    x = np.sqrt((SPS[i,0] - RPS[:,0])**2 + (SPS[i,1] - RPS[:,1])**2)

    refractions = pyf.get_analytical_refractions(v,z,x)

    for k in range(nr):
        eikonal_an[k] = min(x[k]/v[0], np.min(refractions[:,k]))

    eikonal_nu = pyf.read_binary_array(nr, f"../outputs/data/modeling_test_eikonal_iso_nStations{nr}_shot_{i+1}.bin")

    ax[i].plot(eikonal_an - eikonal_nu, "k")

    ax[i].set_ylabel("(Ta - Tn) [ms]", fontsize = 15)
    ax[i].set_xlabel("Channel index", fontsize = 15)
    
    ax[i].set_yticks(np.linspace(-0.005, 0.005, 5))
    ax[i].set_yticklabels(np.linspace(-5, 5, 5, dtype = float))

    ax[i].set_xticks(np.linspace(0, nr, 11))
    ax[i].set_xticklabels(np.linspace(0, nr, 11, dtype = int))

    ax[i].set_xlim([0, nr])

    ax[i].invert_yaxis()

fig.tight_layout()
plt.savefig("modeling_test_accuracy.png", dpi = 200)
plt.show()