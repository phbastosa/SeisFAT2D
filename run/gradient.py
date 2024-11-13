import numpy as np
import matplotlib.pyplot as plt

nx = 203
nz = 53

adjoint_comp = np.fromfile("adjoint_comp.bin", dtype = np.float32, count = nx*nz)
adjoint_comp = np.reshape(adjoint_comp, [nz,nx], order = "F")

adjoint_grad = np.fromfile("adjoint_grad.bin", dtype = np.float32, count = nx*nz)
adjoint_grad = np.reshape(adjoint_grad, [nz,nx], order = "F")

adjoint_full = np.zeros((nz,nx))

adj_max = np.max(adjoint_grad)
adj_min = np.min(adjoint_grad)

comp_max = np.max(adjoint_comp)
comp_min = np.min(adjoint_comp)

for i in range(nz):
    for j in range(nx):
       
        adjoint_full[i,j] = (adjoint_grad[i,j] - adjoint_comp[i,j]) / (adjoint_comp[i,j] + 1e-6)


gradient = np.fromfile("gradient.bin", dtype = np.float32, count = (nx-2)*(nz-2))
gradient = np.reshape(gradient, [nz-2,nx-2], order = "F")

fft_gradient = np.fft.fftshift(np.fft.fftn(gradient))

fft_gradient[25,100] = 0.0

gradient_filtered = np.real(np.fft.ifftn(np.fft.fftshift(fft_gradient)))*gradient

plt.subplot(311)
plt.imshow(gradient)

plt.subplot(312)
plt.imshow(np.abs(fft_gradient))

plt.subplot(313)
plt.imshow(gradient_filtered)

plt.show()


# plt.subplot(311)
# plt.imshow(adjoint_grad)

# plt.subplot(312)
# plt.imshow(adjoint_comp)

# plt.subplot(313)
# plt.imshow(adjoint_full)

# plt.show()