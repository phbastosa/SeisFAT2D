import sys; sys.path.append("../src/")

import numpy as np
import functions as pyf

parameters = str(sys.argv[1])

ns = 61
nr = 501

SPS = np.zeros((ns, 2))
RPS = np.zeros((nr, 2))
XPS = np.zeros((ns, 3))

SPS[:, 0] = np.linspace(1000, 4000, ns) 
SPS[:, 1] = 10.0 

RPS[:, 0] = np.linspace(0, 5000, nr)
RPS[:, 1] = 10.0 

spread = 2000
ds = 50
dr = 10

XPS[:, 0] = np.arange(ns)
XPS[:, 1] = np.arange(ns)*ds/dr 
XPS[:, 2] = np.arange(ns)*ds/dr + spread/dr + 1 

np.savetxt(pyf.catch_parameter(parameters, "SPS"), SPS, fmt = "%.2f", delimiter = ",")
np.savetxt(pyf.catch_parameter(parameters, "RPS"), RPS, fmt = "%.2f", delimiter = ",")
np.savetxt(pyf.catch_parameter(parameters, "XPS"), XPS, fmt = "%.0f", delimiter = ",")
