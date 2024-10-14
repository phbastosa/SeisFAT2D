import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import pyFunctions as pyf

SPS = np.loadtxt("../inputs/geometry/inversion_test_SPS.txt", delimiter = ",", comments = "#", dtype = float) 
RPS = np.loadtxt("../inputs/geometry/inversion_test_RPS.txt", delimiter = ",", comments = "#", dtype = float) 
XPS = np.loadtxt("../inputs/geometry/inversion_test_XPS.txt", delimiter = ",", comments = "#", dtype = int) 

m2km = 1e-3

nx = 2001
nz = 501 

dh = 10.0

trueModel = pyf.read_binary_matrix(nz, nx, "../inputs/models/inversion_test_true_model_501x2001_10m.bin")
initModel = pyf.read_binary_matrix(nz, nx, "../inputs/models/inversion_test_init_model_501x2001_10m.bin")



fig, ax = plt.subplots(nrows = 2, figsize = (15,6))

im = ax[0].imshow(trueModel, aspect = "auto", cmap = "jet")

ax[0].set_xlabel("Distance [m]", fontsize = 15)
ax[0].set_ylabel("Depth [m]", fontsize = 15)

ax[0].plot(RPS[:,0]/dh, RPS[:,1]/dh, "ok")
ax[0].plot(SPS[:,0]/dh, SPS[:,1]/dh, "og")

cbar = plt.colorbar(im)

im = ax[1].imshow(initModel, aspect = "auto", cmap = "jet")

ax[1].plot(RPS[:,0]/dh, RPS[:,1]/dh, "ok")
ax[1].plot(SPS[:,0]/dh, SPS[:,1]/dh, "og")

ax[1].set_xlabel("Distance [m]", fontsize = 15)
ax[1].set_ylabel("Depth [m]", fontsize = 15)

cbar = plt.colorbar(im)

fig.tight_layout()
plt.show()