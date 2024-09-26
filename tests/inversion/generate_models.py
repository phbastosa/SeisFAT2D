import numpy as np

x_max = 2e4
z_max = 3e3

dh = 5.0

nx = int((x_max / dh) + 1)
nz = int((z_max / dh) + 1)

true_model = np.zeros((nz, nx))
init_model = np.zeros((nz, nx))

dv = 20
vi = 2000

for i in range(nz):
    true_model[i] = vi + i*dv/dh
    init_model[i] = vi + i*dv/dh

radius = 1000

velocity_variation = np.array([-500, 500])

circle_centers = np.array([[1500, 8000],
                           [1500, 12000]])

x, z = np.meshgrid(np.arange(nx)*dh, np.arange(nz)*dh)

for k, dv in enumerate(velocity_variation):
    
    distance = np.sqrt((x - circle_centers[k,1])**2 + (z - circle_centers[k,0])**2)

    true_model[distance <= radius] += dv

true_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/inversion_test_true_model_{nz}x{nx}_{dh:.0f}m.bin")
init_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/inversion_test_init_model_{nz}x{nx}_{dh:.0f}m.bin")


