import numpy as np
import matplotlib.pyplot as plt

NNZ = 4031901

nx = 101
nz = 101

n_model = nx*nz

SPS = np.loadtxt("../inputs/geometry/anisoTomo_SPS.txt", delimiter = ",", dtype = float)
RPS = np.loadtxt("../inputs/geometry/anisoTomo_RPS.txt", delimiter = ",", dtype = float)
XPS = np.loadtxt("../inputs/geometry/anisoTomo_XPS.txt", delimiter = ",", dtype = int)

spread = XPS[0,2] - XPS[0,1]
n_data = len(SPS) * spread

N = n_data
M = 3*n_model

G = np.zeros((N,M))
m = np.zeros(M)

iA = np.fromfile(f"iA_{NNZ}.bin", dtype = np.int32, count = NNZ)
jA = np.fromfile(f"jA_{NNZ}.bin", dtype = np.int32, count = NNZ)
vA = np.fromfile(f"vA_{NNZ}.bin", dtype = np.float32, count = NNZ)

for p in range(NNZ):
    if iA[p] < n_data:
        G[iA[p],jA[p]] = vA[p]

I = np.sum(G, axis = 0)

Is = np.reshape(I[:n_model], [nz,nx], order = "F")
Ie = np.reshape(I[n_model:2*n_model], [nz,nx], order = "F")
Id = np.reshape(I[2*n_model:], [nz,nx], order = "F")

folder = "../outputs/recoveredModels/tomography_vti_final_model_"

S = 1.0/np.fromfile(f"{folder}vp_101x101.bin", dtype = np.float32, count = nx*nz).reshape([nz,nx], order = "F")
E = np.fromfile(f"{folder}ep_101x101.bin", dtype = np.float32, count = nx*nz).reshape([nz,nx], order = "F")
D = np.fromfile(f"{folder}dl_101x101.bin", dtype = np.float32, count = nx*nz).reshape([nz,nx], order = "F")

plt.subplot(231)
plt.imshow(Is)

plt.subplot(232)
plt.imshow(Ie)

plt.subplot(233)
plt.imshow(Id)

plt.subplot(234)
plt.imshow(S)

plt.subplot(235)
plt.imshow(E)

plt.subplot(236)
plt.imshow(D)

plt.show()