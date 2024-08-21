import numpy as np

nx = 2001
nz = 501

dh = 10.0

init_model = np.zeros((nz, nx))
true_model = np.zeros((nz, nx))

dv = 50
vi = 2000
wb = 1000

for i in range(nz):
    init_model[i] = vi + i*dv/dh if i*dh > wb else 1500.0
    true_model[i] = vi + i*dv/dh if i*dh > wb else 1500.0

radius = 1250

velocity_variation = np.array([-500, 500])

circle_centers = np.array([[3000, 7500],
                           [3000, 12500]])

x, z = np.meshgrid(np.arange(nx)*dh, np.arange(nz)*dh)

for k, dv in enumerate(velocity_variation):
    
    distance = np.sqrt((x - circle_centers[k,1])**2 + (z - circle_centers[k,0])**2)

    true_model[distance <= radius] += dv

true_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/model/inversion_trueModelTest_{nz}x{nx}_{dh:.0f}m.bin")
init_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/model/inversion_initModelTest_{nz}x{nx}_{dh:.0f}m.bin")
