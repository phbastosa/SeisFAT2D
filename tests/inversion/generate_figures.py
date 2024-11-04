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

diffModel = trueModel - initModel

xloc = np.linspace(0, nx-1, 11)
xlab = np.array(xloc*dh*m2km, dtype = int)

zloc = np.linspace(0, nz-1, 5)
zlab = np.array(zloc*dh*m2km, dtype = int)

fig, ax = plt.subplots(nrows = 3, figsize = (15, 8))

im = ax[0].imshow(trueModel, aspect = "auto", cmap = "jet", vmin = 1500, vmax = 4000)

ax[0].set_xticks(xloc)
ax[0].set_yticks(zloc)
ax[0].set_xticklabels(xlab)    
ax[0].set_yticklabels(zlab)    
ax[0].set_ylabel("Depth [km]", fontsize = 15)
ax[0].set_xlabel("Distance [km]", fontsize = 15)

ax[0].plot(RPS[:,0]/dh, RPS[:,1]/dh, "ok")
ax[0].plot(SPS[:,0]/dh, SPS[:,1]/dh, "og")

ax[0].set_title("Reference model", fontsize = 15)

cbar = plt.colorbar(im)
cbar.set_label("Velocity [m/s]")

im = ax[1].imshow(initModel, aspect = "auto", cmap = "jet", vmin = 1500, vmax = 4000)

ax[1].plot(RPS[:,0]/dh, RPS[:,1]/dh, "ok")
ax[1].plot(SPS[:,0]/dh, SPS[:,1]/dh, "og")

ax[1].set_xticks(xloc)
ax[1].set_yticks(zloc)
ax[1].set_xticklabels(xlab)    
ax[1].set_yticklabels(zlab)    
ax[1].set_ylabel("Depth [km]", fontsize = 15)
ax[1].set_xlabel("Distance [km]", fontsize = 15)

ax[1].set_title("Initial model", fontsize = 15)

cbar = plt.colorbar(im)
cbar.set_label("Velocity [m/s]")

im = ax[2].imshow(diffModel, aspect = "auto", cmap = "bwr", vmin = -600, vmax = 600)

ax[2].plot(RPS[:,0]/dh, RPS[:,1]/dh, "ok")
ax[2].plot(SPS[:,0]/dh, SPS[:,1]/dh, "og")

ax[2].set_xticks(xloc)
ax[2].set_yticks(zloc)
ax[2].set_xticklabels(xlab)    
ax[2].set_yticklabels(zlab)    
ax[2].set_ylabel("Depth [km]", fontsize = 15)
ax[2].set_xlabel("Distance [km]", fontsize = 15)

ax[2].set_title("Difference model", fontsize = 15)

cbar = plt.colorbar(im)
cbar.set_label("Velocity [m/s]")

fig.tight_layout()
plt.savefig("configuration.png", dpi = 300)
plt.show()


leastSquaresModel = pyf.read_binary_matrix(nz, nx, "../outputs/recoveredModels/least_squares_final_model_501x2001.bin")
adjointStateModel = pyf.read_binary_matrix(nz, nx, "../outputs/recoveredModels/adjoint_state_final_model_501x2001.bin")

leastSquaresDiff = leastSquaresModel - initModel
adjointStateDiff = adjointStateModel - initModel

fig, ax = plt.subplots(nrows = 2, figsize = (15, 8))

im = ax[0].imshow(leastSquaresDiff, aspect = "auto", cmap = "bwr", vmin = -600, vmax = 600)

ax[0].set_xticks(xloc)
ax[0].set_yticks(zloc)
ax[0].set_xticklabels(xlab)    
ax[0].set_yticklabels(zlab)    
ax[0].set_ylabel("Depth [km]", fontsize = 15)
ax[0].set_xlabel("Distance [km]", fontsize = 15)

ax[0].plot(RPS[:,0]/dh, RPS[:,1]/dh, "ok")
ax[0].plot(SPS[:,0]/dh, SPS[:,1]/dh, "og")

ax[0].set_title("Least-Squares Tomography", fontsize = 15)

cbar = plt.colorbar(im)
cbar.set_label("Velocity [m/s]")

im = ax[1].imshow(adjointStateDiff, aspect = "auto", cmap = "bwr", vmin = -600, vmax = 600)

ax[1].plot(RPS[:,0]/dh, RPS[:,1]/dh, "ok")
ax[1].plot(SPS[:,0]/dh, SPS[:,1]/dh, "og")

ax[1].set_xticks(xloc)
ax[1].set_yticks(zloc)
ax[1].set_xticklabels(xlab)    
ax[1].set_yticklabels(zlab)    
ax[1].set_ylabel("Depth [km]", fontsize = 15)
ax[1].set_xlabel("Distance [km]", fontsize = 15)

ax[1].set_title("Adjoint-State Tomography", fontsize = 15)

cbar = plt.colorbar(im)
cbar.set_label("Velocity [m/s]")

fig.tight_layout()
plt.savefig("model_results.png", dpi = 300)
plt.show()

leastSquaresCurve = np.loadtxt("../outputs/convergence/least_squares_convergence_5_iterations.txt", dtype = np.float32)
adjointStateCurve = np.loadtxt("../outputs/convergence/adjoint_state_convergence_5_iterations.txt", dtype = np.float32)

plt.figure(1, figsize = (8,5))

plt.plot(leastSquaresCurve, "--ob")
plt.plot(adjointStateCurve, "--og")

plt.xlabel("Iterations", fontsize = 15)
plt.ylabel(r"$||d^{obs} - d^{cal}||^2_2$", fontsize = 15)

plt.tight_layout()
plt.savefig("convergence.png", dpi = 300)
plt.show()
