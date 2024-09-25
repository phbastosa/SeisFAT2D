import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import pyFunctions as pyf

SPS = np.loadtxt("../inputs/geometry/modeling_test_SPS.txt", delimiter = ",", comments = "#", dtype = float) 
RPS = np.loadtxt("../inputs/geometry/modeling_test_RPS.txt", delimiter = ",", comments = "#", dtype = float) 
XPS = np.loadtxt("../inputs/geometry/modeling_test_XPS.txt", delimiter = ",", comments = "#", dtype = int) 

m2km = 1e-3

nx = 4001
nz = 601 

dh = 5.0

model = pyf.read_binary_matrix(nz, nx, f"../inputs/models/modeling_test_model_{nz}x{nx}_{dh:.0f}m.bin")

xloc = np.linspace(0, nx-1, 11)
xlab = np.array(xloc*dh*m2km, dtype = int)

zloc = np.linspace(0, nz-1, 5)
zlab = np.array(zloc*dh*m2km, dtype = int)

fig, ax = plt.subplots(figsize = (15, 5), ncols = 1, nrows = 2)

im = ax[0].imshow(model, aspect = "auto", cmap = "Greys")

ax[0].plot(RPS[:, 0]/dh, RPS[:, 1]/dh, "ob")
ax[0].plot(SPS[:, 0]/dh, SPS[:, 1]/dh, "ok")

ax[0].set_xticks(xloc)
ax[0].set_xticklabels(xlab)    
ax[0].set_yticks(zloc)
ax[0].set_yticklabels(zlab)    
ax[0].set_ylabel("Depth [km]", fontsize = 15)
ax[0].set_xlabel("Distance [km]", fontsize = 15)

for i in range(len(SPS)):
    
    data = pyf.read_binary_array(len(RPS), f"../outputs/syntheticData/travel_time_shot_{i}.bin")
    
    ax[1].plot(RPS[XPS[i,1]:XPS[i,2],0], data)

ax[1].set_xticks(xloc*dh)
ax[1].set_xticklabels(xlab)    
ax[1].set_ylim([0, 3.2])
ax[1].set_xlim([0, (nx-1)*dh])
ax[1].set_ylabel("Time [s]", fontsize = 15)
ax[1].set_xlabel("Distance [km]", fontsize = 15)
ax[1].invert_yaxis()

fig.tight_layout()
plt.savefig("modeling_test.png", dpi = 200)
plt.show()
    






