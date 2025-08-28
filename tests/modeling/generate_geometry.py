import sys; sys.path.append("../src/")

import numpy as np
import functions as pyf

parameters = str(sys.argv[1])

ns = 3
nr = 401

SPS = np.zeros((ns, 2))
RPS = np.zeros((nr, 2))
XPS = np.zeros((ns, 3))

SPS[:, 0] = np.linspace(100, 19900, ns) 
SPS[:, 1] = 0.0 

RPS[:, 0] = np.linspace(0, 20000, nr)
RPS[:, 1] = 0.0 

XPS[:, 0] = np.arange(ns)
XPS[:, 1] = np.zeros(ns) 
XPS[:, 2] = np.zeros(ns) + nr 

np.savetxt(pyf.catch_parameter(parameters, "SPS"), SPS, fmt = "%.2f", delimiter = ",")
np.savetxt(pyf.catch_parameter(parameters, "RPS"), RPS, fmt = "%.2f", delimiter = ",")
np.savetxt(pyf.catch_parameter(parameters, "XPS"), XPS, fmt = "%.0f", delimiter = ",")
