import numpy as np

x_max = 1e4
z_max = 1e3

dh = 10.0

nx = int((x_max / dh) + 1)
nz = int((z_max / dh) + 1)

vp = np.zeros((nz, nx))

dvx = 10.0
dvz = 20.0
vi = 2000.0

for i in range(nz):
    for j in range(nx):
        vp[i,j] = vi + i*dvz/dh + j*dvx/dh

radius = 250

velocity_variation = np.array([-500, 500])

circle_centers = np.array([[500, 4000],
                           [500, 6000]])

x, z = np.meshgrid(np.arange(nx)*dh, np.arange(nz)*dh)

for k, dv in enumerate(velocity_variation):
    
    distance = np.sqrt((x - circle_centers[k,1])**2 + (z - circle_centers[k,0])**2)

    vp[distance <= radius] += dv

vs = vp / 1.7
rho = 310.0*vp**0.25

vp.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/modeling_test_vp_{nz}x{nx}_{dh:.0f}m.bin")
vs.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/modeling_test_vs_{nz}x{nx}_{dh:.0f}m.bin")
rho.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/modeling_test_rho_{nz}x{nx}_{dh:.0f}m.bin")


