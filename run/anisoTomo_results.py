import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable

nx = 201
nz = 51

dh = 100

model_true = np.fromfile(f"../inputs/models/true_vp.bin", count = nx*nz, dtype = np.float32).reshape([nz,nx], order = "F")
model_init = np.fromfile(f"../inputs/models/init_vp.bin", count = nx*nz, dtype = np.float32).reshape([nz,nx], order = "F")

model_iso = np.fromfile(f"../outputs/recoveredModels/tomography_iso_final_model_vp_51x201.bin", count = nx*nz, dtype = np.float32).reshape([nz,nx], order = "F")
# model_ani = np.fromfile(f"../outputs/recoveredModels/ani_tomography_iso_final_model_101x101.bin", count = nx*nz, dtype = np.float32).reshape([nz,nx], order = "F")
model_ani = np.zeros_like(model_true)

SPS = np.loadtxt("../inputs/geometry/anisoTomo_SPS.txt", delimiter = ",", dtype = float)
RPS = np.loadtxt("../inputs/geometry/anisoTomo_RPS.txt", delimiter = ",", dtype = float)

m2km = 1e-3

vmin = np.min(model_true)
vmax = np.max(model_true)

dvmin =-500
dvmax = 500

cmap = mpl.colormaps["jet"]
norm = mpl.colors.Normalize(vmin*m2km, vmax*m2km)

dcmap = mpl.colormaps["bwr"]
dnorm = mpl.colors.Normalize(dvmin*m2km, dvmax*m2km)

xloc = np.linspace(0, nx-1, 6)
xlab = np.linspace(0, m2km*(nx-1)*dh, 6, dtype = int)

zloc = np.linspace(0, nz-1, 6)
zlab = np.linspace(0, m2km*(nz-1)*dh, 6, dtype = int)

fig, ax = plt.subplots(ncols = 2, nrows = 3, figsize = (19, 7))

im = ax[0,0].imshow(model_true, aspect = "auto", cmap = cmap, vmin = vmin, vmax = vmax)
ax[0,0].plot(np.zeros(nz) + 0.5*(nx-1), np.arange(nz), "--k")
ax[0,0].plot(SPS[:,0]/dh, SPS[:,1]/dh, "o", color = "gray")
ax[0,0].plot(RPS[:,0]/dh, RPS[:,1]/dh, "o", color = "green")
ax[0,0].set_xlabel("x [km]", fontsize = 15)
ax[0,0].set_ylabel("z [km]", fontsize = 15)
ax[0,0].set_xticks(xloc)
ax[0,0].set_yticks(zloc)
ax[0,0].set_xticklabels(xlab)
ax[0,0].set_yticklabels(zlab)
divider = make_axes_locatable(ax[0,0])
cax = divider.append_axes("right", size = "2.5%", pad = 0.1)
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm = norm, cmap = cmap), cax = cax, ticks = np.linspace(vmin*m2km, vmax*m2km, 5), orientation = "vertical")
cbar.ax.set_yticklabels(np.around(np.linspace(vmin*m2km, vmax*m2km, 5), decimals = 1))
cbar.set_label("Velocity [km/s]", fontsize = 15)

im = ax[0,1].imshow(model_true - model_init, aspect = "auto", cmap = dcmap, vmin = dvmin, vmax = dvmax)
ax[0,1].plot(np.zeros(nz) + 0.5*(nx-1), np.arange(nz), "--k")
ax[0,1].plot(SPS[:,0]/dh, SPS[:,1]/dh, "o", color = "gray")
ax[0,1].plot(RPS[:,0]/dh, RPS[:,1]/dh, "o", color = "green")
ax[0,1].set_xlabel("x [km]", fontsize = 15)
ax[0,1].set_ylabel("z [km]", fontsize = 15)
ax[0,1].set_xticks(xloc)
ax[0,1].set_yticks(zloc)
ax[0,1].set_xticklabels(xlab)
ax[0,1].set_yticklabels(zlab)
divider = make_axes_locatable(ax[0,1])
cax = divider.append_axes("right", size = "2.5%", pad = 0.1)
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm = dnorm, cmap = dcmap), cax = cax, ticks = np.linspace(dvmin*m2km, dvmax*m2km, 5), orientation = "vertical")
cbar.ax.set_yticklabels(np.around(np.linspace(dvmin*m2km, dvmax*m2km, 5), decimals = 1))
cbar.set_label("Velocity [km/s]", fontsize = 15)

im = ax[1,0].imshow(model_iso, aspect = "auto", cmap = cmap, vmin = vmin, vmax = vmax)
ax[1,0].plot(np.zeros(nz) + 0.5*(nx-1), np.arange(nz), "--k")
ax[1,0].plot(SPS[:,0]/dh, SPS[:,1]/dh, "o", color = "gray")
ax[1,0].plot(RPS[:,0]/dh, RPS[:,1]/dh, "o", color = "green")
ax[1,0].set_xlabel("x [km]", fontsize = 15)
ax[1,0].set_ylabel("z [km]", fontsize = 15)
ax[1,0].set_xticks(xloc)
ax[1,0].set_yticks(zloc)
ax[1,0].set_xticklabels(xlab)
ax[1,0].set_yticklabels(zlab)
divider = make_axes_locatable(ax[1,0])
cax = divider.append_axes("right", size = "2.5%", pad = 0.1)
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm = norm, cmap = cmap), cax = cax, ticks = np.linspace(vmin*m2km, vmax*m2km, 5), orientation = "vertical")
cbar.ax.set_yticklabels(np.around(np.linspace(vmin*m2km, vmax*m2km, 5), decimals = 1))
cbar.set_label("P wave velocity [km/s]", fontsize = 15)

im = ax[1,1].imshow(model_iso - model_init, aspect = "auto", cmap = dcmap, vmin = dvmin, vmax = dvmax)
ax[1,1].plot(np.zeros(nz) + 0.5*(nx-1), np.arange(nz), "--k")
ax[1,1].plot(SPS[:,0]/dh, SPS[:,1]/dh, "o", color = "gray")
ax[1,1].plot(RPS[:,0]/dh, RPS[:,1]/dh, "o", color = "green")
ax[1,1].set_xlabel("x [km]", fontsize = 15)
ax[1,1].set_ylabel("z [km]", fontsize = 15)
ax[1,1].set_xticks(xloc)
ax[1,1].set_yticks(zloc)
ax[1,1].set_xticklabels(xlab)
ax[1,1].set_yticklabels(zlab)
divider = make_axes_locatable(ax[1,1])
cax = divider.append_axes("right", size = "2.5%", pad = 0.1)
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm = dnorm, cmap = dcmap), cax = cax, ticks = np.linspace(dvmin*m2km, dvmax*m2km, 5), orientation = "vertical")
cbar.ax.set_yticklabels(np.around(np.linspace(dvmin*m2km, dvmax*m2km, 5), decimals = 1))
cbar.set_label("Velocity [km/s]", fontsize = 15)


# ax[1,1].imshow(model_ani, cmap = "jet", vmin = vmin, vmax = vmax)
# ax[1,1].plot(np.zeros(nz) + 0.5*(nx-1), np.arange(nz), "--k")
# ax[1,1].plot(SPS[:,0]/dh, SPS[:,1]/dh, "o", color = "gray")
# ax[1,1].plot(RPS[:,0]/dh, RPS[:,1]/dh, "o", color = "green")
# ax[1,1].set_xlabel("x [km]", fontsize = 15)
# ax[1,1].set_ylabel("z [km]", fontsize = 15)
# ax[1,1].set_xticks(xloc)
# ax[1,1].set_yticks(zloc)
# ax[1,1].set_xticklabels(xlab)
# ax[1,1].set_yticklabels(zlab)
# divider = make_axes_locatable(ax[1,1])
# cax = divider.append_axes("right", size = "2.5%", pad = 0.1)
# cbar = fig.colorbar(mpl.cm.ScalarMappable(norm = norm, cmap = cmap), cax = cax, ticks = np.linspace(vmin*m2km, vmax*m2km, 5), orientation = "vertical")
# cbar.ax.set_yticklabels(np.around(np.linspace(vmin*m2km, vmax*m2km, 5), decimals = 1))
# cbar.set_label("P wave velocity [km/s]", fontsize = 15)

fig.tight_layout()
# plt.savefig("motivation2D.png", dpi = 300)
plt.show()

# depth = np.arange(nz)*dh*m2km

# fig, ax = plt.subplots(figsize = (4,8))

# ax.plot(model_true[:,int(0.5*nx)] - model_init[:,int(0.5*nx)], depth, color = "black", label = "Reference")
# ax.plot(model_iso[:,int(0.5*nx)] - model_init[:,int(0.5*nx)], depth, color = "green", label = "Isotropic")
# ax.plot(model_ani[:,int(0.5*nx)] - model_init[:,int(0.5*nx)], depth, color = "red", label = "Anisotropic")

# ax.set_xlim([-1000,1000])
# ax.set_ylim([0, m2km*(nz-1)*dh])
# ax.set_xlabel("P wave velocity [m/s]", fontsize = 15)
# ax.set_ylabel("Depth [km]", fontsize = 15)

# ax.legend(loc = "upper right", fontsize = 12)

# ax.invert_yaxis()
# fig.tight_layout()
# plt.savefig("welllogs2D.png", dpi = 300)
# plt.show()

# residuo_iso = np.loadtxt("../outputs/convergence/tomography_iso_convergence_5_iterations.txt", dtype = float)
# residuo_ani = np.loadtxt("../outputs/convergence/ani_tomography_iso_convergence_5_iterations.txt", dtype = float)

# residuo_iso *= 100.0/np.max(residuo_iso)
# residuo_ani *= 100.0/np.max(residuo_ani)

# fig, ax = plt.subplots(figsize = (10,4))

# ax.plot(residuo_iso, "--or", label = "Isotropic")
# ax.plot(residuo_ani, "--og", label = "Anisotropic")

# ax.set_ylabel("Residuals [%]", fontsize = 15)
# ax.set_xlabel("Iterations", fontsize = 15)

# ax.legend(loc = "upper right", fontsize = 12)
# fig.tight_layout()
# plt.show()