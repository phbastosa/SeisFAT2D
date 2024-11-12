import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

SPS = np.loadtxt("../inputs/geometry/migration_test_SPS.txt", delimiter = ",", comments = "#", dtype = float) 
RPS = np.loadtxt("../inputs/geometry/migration_test_RPS.txt", delimiter = ",", comments = "#", dtype = float) 
XPS = np.loadtxt("../inputs/geometry/migration_test_XPS.txt", delimiter = ",", comments = "#", dtype = int) 

m2km = 1e-3

nx = 2001
nz = 201 

dh = 5.0

model = pyf.read_binary_matrix(nz, nx, "../inputs/models/migration_test_vp_model_201x2001_5m.bin")
image = pyf.read_binary_matrix(nz, nx, "../outputs/migratedImages/kirchhoff_result_201x2001.bin")

xloc = np.linspace(0, nx-1, 11)
xlab = np.array(xloc*dh*m2km, dtype = int)

zloc = np.linspace(0, nz-1, 5)
zlab = np.array(zloc*dh*m2km)

fig, ax = plt.subplots(nrows = 2, figsize = (15, 5))

im = ax[0].imshow(model, aspect = "auto", cmap = "jet", vmin = 2000, vmax = 2500)

ax[0].set_xticks(xloc)
ax[0].set_yticks(zloc)
ax[0].set_xticklabels(xlab)    
ax[0].set_yticklabels(zlab)    
ax[0].set_ylabel("Depth [km]", fontsize = 15)
ax[0].set_xlabel("Distance [km]", fontsize = 15)

ax[0].plot(RPS[:,0]/dh, RPS[:,1]/dh, "ok")
ax[0].plot(SPS[:,0]/dh, SPS[:,1]/dh, "og")

ax[0].set_title("Model", fontsize = 15)

cbar = plt.colorbar(im)
cbar.set_label("Velocity [m/s]")

im = ax[1].imshow(image, aspect = "auto", cmap = "Greys")

ax[1].plot(RPS[:,0]/dh, RPS[:,1]/dh, "ok")
ax[1].plot(SPS[:,0]/dh, SPS[:,1]/dh, "og")

ax[1].set_xticks(xloc)
ax[1].set_yticks(zloc)
ax[1].set_xticklabels(xlab)    
ax[1].set_yticklabels(zlab)    
ax[1].set_ylabel("Depth [km]", fontsize = 15)
ax[1].set_xlabel("Distance [km]", fontsize = 15)

ax[1].set_title("Image", fontsize = 15)

cbar = plt.colorbar(im)
cbar.set_label("Amplitude")

fig.tight_layout()
plt.savefig("migration_test_results.png", dpi = 300)
plt.show()
