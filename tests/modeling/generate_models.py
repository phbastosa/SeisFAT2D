import numpy as np

x_max = 2e4
z_max = 5e3

dh = 10.0

nx = int((x_max / dh) + 1)
nz = int((z_max / dh) + 1)

vp_model = np.zeros((nz, nx)) + 1500
vs_model = np.zeros((nz, nx)) 
rho_model = np.zeros((nz, nx)) + 1000

dv = 5.0
vi = 1650.0
wb = 1000.0

for i in range(nz):
    if i > wb/dh:
        vp_model[i] = vi + (i*dh - wb)*dv/dh 

radius = 1000

velocity_variation = np.array([-500, 500])

circle_centers = np.array([[3000, 8000],
                           [3000, 12000]])

x, z = np.meshgrid(np.arange(nx)*dh, np.arange(nz)*dh)

for k, dv in enumerate(velocity_variation):
    
    distance = np.sqrt((x - circle_centers[k,1])**2 + (z - circle_centers[k,0])**2)

    vp_model[distance <= radius] += dv


vs_model[int(wb/dh):] = 0.7*vp_model[int(wb/dh):]
rho_model[int(wb/dh):] = 310*vp_model[int(wb/dh):]**0.25

vp_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/modeling_test_vp_model_{nz}x{nx}_{dh:.0f}m.bin")
vs_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/modeling_test_vs_model_{nz}x{nx}_{dh:.0f}m.bin")
rho_model.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/models/modeling_test_rho_model_{nz}x{nx}_{dh:.0f}m.bin")
