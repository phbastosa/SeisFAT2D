import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import pyFunctions as pyf

SPS = np.loadtxt("../inputs/geometry/SeisFAT2D_SPS.txt", delimiter = ",", comments = "#", dtype = float) 
RPS = np.loadtxt("../inputs/geometry/SeisFAT2D_RPS.txt", delimiter = ",", comments = "#", dtype = float) 
XPS = np.loadtxt("../inputs/geometry/SeisFAT2D_XPS.txt", delimiter = ",", comments = "#", dtype = int) 

nd = len(RPS)

nx = 4001
nz = 1001 

dh = 5.0

velocities = np.array([2000, 2500, 3000, 4200])
thickness = np.array([475, 1000, 2500])
distance = np.linspace(50, 19950, nd)

model = pyf.read_binary_matrix(nz, nx, f"../inputs/model/modeling_test_model_{nz}x{nx}_{dh:.0f}m.bin")


fig, ax = plt.subplots(figsize = (15, 8), ncols = 1, nrows = 3)

xloc = np.linspace(0, nx-1, 11)
xlab = np.array(xloc*dh/1000, dtype = int)

zloc = np.linspace(0, nz-1, 5)
zlab = np.array(zloc*dh/1000, dtype = int)

im = ax[0].imshow(model, aspect = "auto", cmap = "Greys")

ax[0].plot(RPS[XPS[0,1]:XPS[0,2],0]/dh, RPS[XPS[0,1]:XPS[0,2],1]/dh, "kv")

ax[0].set_xticks(xloc)
ax[0].set_xticklabels(xlab)    

ax[0].set_yticks(zloc)
ax[0].set_yticklabels(zlab)    

ax[0].set_ylabel("Depth [km]", fontsize = 15)
ax[0].set_xlabel("Distance [km]", fontsize = 15)

for i in range(len(XPS)):

    ax[0].plot(SPS[XPS[i,0],0]/dh, SPS[XPS[i,0],1]/dh, "o", markersize = 10)
    
    x = SPS[XPS[i,0],0] - RPS[XPS[i,1]:XPS[i,2],0]
    z = SPS[XPS[i,0],1] - RPS[XPS[i,1]:XPS[i,2],1]    

    offset = np.sqrt(x**2 + z**2)    

    direct_wave = offset / velocities[0]
    refractions = pyf.get_analytical_refractions(velocities, thickness, offset)

    dcalc = pyf.read_binary_array(nd, f"../outputs/syntheticData/travel_time_{nd}_stations_shot_{i+1}.bin")

    dtrue = np.zeros_like(dcalc)

    for i in range(nd):
        dtrue[i] = min(direct_wave[i], np.min(refractions[:,i]))

    ax[1].plot(distance, dtrue, "k")
    ax[1].plot(distance, dcalc)

    ax[2].plot(distance, 1000*(dtrue - dcalc))

ax[1].set_xticks(xloc*dh)
ax[1].set_xticklabels(xlab)    
ax[1].set_ylim([0, 7])
ax[1].set_xlim([0, (nx-1)*dh])
ax[1].set_ylabel("Time [s]", fontsize = 15)
ax[1].set_xlabel("Distance [km]", fontsize = 15)
ax[1].invert_yaxis()

ax[2].set_xticks(xloc*dh)
ax[2].set_xticklabels(xlab)    
ax[2].set_ylim([-0.5, 3.0])
ax[2].set_xlim([0, (nx-1)*dh])
ax[2].set_ylabel("Error (Ta - Tn) [ms]", fontsize = 15)
ax[2].set_xlabel("Distance [km]", fontsize = 15)
ax[2].invert_yaxis()

fig.tight_layout()
plt.show()

# direct_wave = x / v[0]
# refractions = pyf.get_analytical_refractions(v,z,x)

# dtrue = np.zeros_like(dcalc)

# for i in range(nd):
#     dtrue[i] = min(direct_wave[i], np.min(refractions[:,i]))
    






