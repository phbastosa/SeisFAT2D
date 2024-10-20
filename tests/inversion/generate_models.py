import numpy as np

x_max = 2e4
z_max = 5e3

dh = 10.0

nx = int((x_max / dh) + 1)
nz = int((z_max / dh) + 1)

true_model = np.zeros((nz, nx)) + 1500
init_model = np.zeros((nz, nx)) + 1500

dv = 5.0
vi = 1650.0
wb = 1000.0

for i in range(nz):
    if i > wb/dh:
        true_model[i] = vi + (i*dh - wb)*dv/dh 
        init_model[i] = vi + (i*dh - wb)*dv/dh

radius = 1000

velocity_variation = np.array([500, 500])

circle_centers = np.array([[3000, 10000],
                           [3000, 10000]])

x, z = np.meshgrid(np.arange(nx)*dh, np.arange(nz)*dh)

for k, dv in enumerate(velocity_variation):
    
    distance = np.sqrt((x - circle_centers[k,1])**2 + (z - circle_centers[k,0])**2)

    true_model[distance <= radius] += dv

true_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/inversion_test_true_model_{nz}x{nx}_{dh:.0f}m.bin")
init_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/inversion_test_init_model_{nz}x{nx}_{dh:.0f}m.bin")
