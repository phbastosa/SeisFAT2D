import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

parameters = str(sys.argv[1])

SPS = np.loadtxt(pyf.catch_parameter(parameters,"SPS"), delimiter = ",", dtype = np.float32) 
RPS = np.loadtxt(pyf.catch_parameter(parameters,"RPS"), delimiter = ",", dtype = np.float32) 
XPS = np.loadtxt(pyf.catch_parameter(parameters,"XPS"), delimiter = ",", dtype = np.int32) 

nx = int(pyf.catch_parameter(parameters, "x_samples"))
nz = int(pyf.catch_parameter(parameters, "z_samples"))

dx = float(pyf.catch_parameter(parameters, "x_spacing"))
dz = float(pyf.catch_parameter(parameters, "z_spacing"))

image_folder = pyf.catch_parameter(parameters, "output_image_folder")

model_vp = pyf.read_binary_matrix(nz, nx, pyf.catch_parameter(parameters, "vp_model_file"))
seismic = pyf.read_binary_matrix(nz, nx, image_folder + f"kirchhoff_section_{nz}x{nx}.bin")

xloc = np.linspace(0, nx-1, 11)
xlab = np.array(xloc*dx, dtype = int)

zloc = np.linspace(0, nz-1, 6)
zlab = np.array(zloc*dz, dtype = int)

scale = np.std(seismic)

fig, ax = plt.subplots(figsize = (15, 5))

im = ax.imshow(model_vp, aspect = "auto", cmap = "jet")

ax.imshow(seismic, aspect = "auto", cmap = "Greys", vmin = -scale, vmax = scale, alpha = 0.8)

cbar = plt.colorbar(im)
cbar.set_label("Velocity P [m/s]")

ax.plot(RPS[:, 0]/dx, RPS[:, 1]/dz, "ob")
ax.plot(SPS[:, 0]/dx, SPS[:, 1]/dz, "or")

ax.set_xticks(xloc)
ax.set_yticks(zloc)
ax.set_xticklabels(xlab)    
ax.set_yticklabels(zlab)    
ax.set_ylabel("Depth [km]", fontsize = 15)
ax.set_xlabel("Distance [km]", fontsize = 15)

fig.tight_layout()
plt.savefig("migration_test_result.png", dpi = 200)
plt.show()

