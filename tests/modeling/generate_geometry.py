import numpy as np

ns = 3
nr = 401

SPS = np.zeros((ns, 2))
RPS = np.zeros((nr, 2))
XPS = np.zeros((ns, 3))

SPS[:, 0] = np.linspace(100, 19900, ns) 
SPS[:, 1] = 1000.0 

RPS[:, 0] = np.linspace(0, 20000, nr)
RPS[:, 1] = 10.0 

XPS[:, 0] = np.arange(ns)
XPS[:, 1] = np.zeros(ns) 
XPS[:, 2] = np.zeros(ns) + nr 

np.savetxt("../inputs/geometry/modeling_test_SPS.txt", SPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/modeling_test_RPS.txt", RPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/modeling_test_XPS.txt", XPS, fmt = "%.0f", delimiter = ",")