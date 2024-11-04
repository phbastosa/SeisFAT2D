import numpy as np

ns = 180
nr = 201

spread = 22

SPS = np.zeros((ns, 2), dtype = float)
RPS = np.zeros((nr, 2), dtype = float)
XPS = np.zeros((ns, 3), dtype = int)

SPS[:,0] = np.linspace(525, 9475, ns) 
SPS[:,1] = 0.0 

RPS[:,0] = np.linspace(0, 10000, nr)
RPS[:,1] = 0.0 

XPS[:, 0] = np.arange(ns)
XPS[:, 1] = np.arange(ns)  
XPS[:, 2] = np.arange(ns) + spread

np.savetxt("../inputs/geometry/migration_test_SPS.txt", SPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/migration_test_RPS.txt", RPS, fmt = "%.2f", delimiter = ",")
np.savetxt("../inputs/geometry/migration_test_XPS.txt", XPS, fmt = "%.0f", delimiter = ",")

