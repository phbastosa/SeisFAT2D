import numpy as np

ns = 11
nr = 400

SPS = np.zeros((ns, 2))
RPS = np.zeros((nr, 2))
XPS = np.zeros((ns, 3))

SPS[:, 0] = np.linspace(2500, 17500, ns) 
SPS[:, 1] = 10.0 

RPS[:, 0] = np.linspace(25, 19975, nr)
RPS[:, 1] = 10.0 

XPS[:, 0] = np.arange(ns)
XPS[:, 1] = np.arange(ns) * 20
XPS[:, 2] = np.arange(ns) * 20 + 200

np.savetxt("../inputs/geometry/modeling_test_SPS.txt", SPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/modeling_test_RPS.txt", RPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/modeling_test_XPS.txt", XPS, fmt = "%.0f", delimiter = ",")