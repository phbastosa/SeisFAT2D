import numpy as np

nx = 4001
nz = 1001
dh = 5.0

z = np.array([500, 1000, 2500])
v = np.array([2000, 2500, 3000, 4200])

model = np.ones((nz, nx)) * v[0]

for i in range(len(z)):
    model[int(np.sum(z[:i+1])/dh):] = v[i+1]

model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/model/modeling_test_model_{nz}x{nx}_{dh:.0f}m.bin")


