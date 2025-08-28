import sys; sys.path.append("../src/")

import numpy as np
import functions as pyf
import scipy.ndimage as sp

parameters = str(sys.argv[1])

nx = int(pyf.catch_parameter(parameters, "x_samples"))
nz = int(pyf.catch_parameter(parameters, "z_samples"))

dx = float(pyf.catch_parameter(parameters, "x_spacing"))
dz = float(pyf.catch_parameter(parameters, "z_spacing"))

init_vp = np.zeros((nz,nx)) + 1500
true_vp = np.zeros((nz,nx)) + 1500

r = 250
xc = 1000
zc = 1000

x,z = np.meshgrid(np.arange(nx)*dx, np.arange(nz)*dz)

distance = np.sqrt((x - xc)**2 + (z - zc)**2)

true_vp[distance < r] += 500

true_vp = 1.0 / sp.gaussian_filter(1.0 / true_vp, 10.0)

init_vp.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/inversion_test_init_vp.bin")
true_vp.flatten("F").astype(np.float32, order = "F").tofile("../inputs/models/inversion_test_true_vp.bin")
