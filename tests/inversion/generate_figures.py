import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import functions as pyf

SPS = np.loadtxt("../inputs/geometry/inversion_test_SPS.txt", delimiter = ",", comments = "#", dtype = float) 
RPS = np.loadtxt("../inputs/geometry/inversion_test_RPS.txt", delimiter = ",", comments = "#", dtype = float) 
XPS = np.loadtxt("../inputs/geometry/inversion_test_XPS.txt", delimiter = ",", comments = "#", dtype = int) 

m2km = 1e-3

nx = 201
nz = 51 

dh = 100.0

trueModel = pyf.read_binary_matrix(nz, nx, f"../inputs/models/inversion_test_true_model_{nz}x{nx}_{dh:.0f}m.bin")
initModel = pyf.read_binary_matrix(nz, nx, f"../inputs/models/inversion_test_init_model_{nz}x{nx}_{dh:.0f}m.bin")

diffModel = trueModel - initModel

xloc = np.linspace(0, nx-1, 11)
xlab = np.array(xloc*dh*m2km, dtype = int)

zloc = np.linspace(0, nz-1, 6)
zlab = np.array(zloc*dh*m2km, dtype = int)

fig, ax = plt.subplots(nrows = 3, figsize = (14, 8))

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

im = ax[2].imshow(diffModel, aspect = "auto", cmap = "bwr", vmin = -500, vmax = 500)

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
plt.savefig("inversion_test_configuration.png", dpi = 300)
plt.show()

leastSquaresModel = pyf.read_binary_matrix(nz, nx, "../outputs/recoveredModels/inversion_test_least_squares_final_model_51x201.bin")
adjointStateModel = pyf.read_binary_matrix(nz, nx, "../outputs/recoveredModels/inversion_test_adjoint_state_final_model_51x201.bin")

leastSquaresDiff = leastSquaresModel - initModel
adjointStateDiff = adjointStateModel - initModel

fig, ax = plt.subplots(nrows = 2, figsize = (14, 8))

im = ax[0].imshow(leastSquaresDiff, aspect = "auto", cmap = "bwr", vmin = -500, vmax = 500)

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

im = ax[1].imshow(adjointStateDiff, aspect = "auto", cmap = "bwr", vmin = -0.5, vmax = 0.5)

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
plt.savefig("inversion_test_results.png", dpi = 300)
plt.show()

leastSquaresCurve = np.loadtxt("../outputs/convergence/inversion_test_least_squares_convergence_5_iterations.txt", dtype = np.float32)
adjointStateCurve = np.loadtxt("../outputs/convergence/inversion_test_adjoint_state_convergence_1_iterations.txt", dtype = np.float32)

plt.figure(1, figsize = (15,5))

plt.plot(leastSquaresCurve / np.max(leastSquaresCurve) * 100, "--ob", label = "least-squares tomography")
plt.plot(adjointStateCurve / np.max(adjointStateCurve) * 100, "--og", label = "adjoint-state tomography")

plt.xlabel("Iterations", fontsize = 15)
plt.ylabel(r"$||d^{obs} - d^{cal}||^2_2$ [%]", fontsize = 15)

plt.ylim([-2, 102])
plt.xlim([-0.05, len(leastSquaresCurve)-0.95])

plt.legend(loc = "upper right")
plt.tight_layout()
plt.savefig("inversion_test_convergence.png", dpi = 300)
plt.show()
