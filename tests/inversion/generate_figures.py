import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import pyFunctions as pyf

SPS = np.loadtxt("../inputs/geometry/SeisFAT2D_SPS.txt", delimiter = ",", comments = "#", dtype = float) 
RPS = np.loadtxt("../inputs/geometry/SeisFAT2D_RPS.txt", delimiter = ",", comments = "#", dtype = float) 
XPS = np.loadtxt("../inputs/geometry/SeisFAT2D_XPS.txt", delimiter = ",", comments = "#", dtype = int) 

nx = 2001
nz = 501 
dh = 10.0

true_model = pyf.read_binary_matrix(nz, nx, f"../inputs/model/inversion_trueModelTest_{nz}x{nx}_{dh:.0f}m.bin")

radius = 1250

circle_centers = np.array([[3000, 7500],
                           [3000, 12500]])

circle1 = plt.Circle((circle_centers[0,1], circle_centers[0,0]), radius, color = "k", ls = "--", fill = False)
circle2 = plt.Circle((circle_centers[1,1], circle_centers[1,0]), radius, color = "k", ls = "--", fill = False)

fig, ax = plt.subplots(figsize = (12, 3))

ax.imshow(true_model, aspect = "auto", cmap = "jet", extent = [0, 20000, 5000, 0])

ax.plot(SPS[:,0], SPS[:,1], "o", color = "black", markersize = 5)
ax.plot(RPS[:,0], RPS[:,1], "o", color = "green", markersize = 5)

ax.add_patch(circle1)
ax.add_patch(circle2)
ax.set_xlim([0,20000])
ax.set_ylim([0,5000])
ax.invert_yaxis()
fig.tight_layout()
plt.show()