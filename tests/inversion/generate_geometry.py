import sys; sys.path.append("../src/")

import numpy as np
import functions as pyf

parameters = str(sys.argv[1])

ns = 200
nr = 200

SPS = np.zeros((ns, 2))
RPS = np.zeros((nr, 2))
XPS = np.zeros((ns, 3))

SPS[:, 0] = 5.0 
SPS[:, 1] = np.linspace(5, 1995, ns) 

RPS[:, 0] = 1995.0
RPS[:, 1] = np.linspace(5, 1995, nr) 

XPS[:, 0] = np.arange(ns)
XPS[:, 1] = np.zeros(ns) 
XPS[:, 2] = np.zeros(ns) + nr 

np.savetxt(pyf.catch_parameter(parameters, "SPS"), SPS, fmt = "%.2f", delimiter = ",")
np.savetxt(pyf.catch_parameter(parameters, "RPS"), RPS, fmt = "%.2f", delimiter = ",")
np.savetxt(pyf.catch_parameter(parameters, "XPS"), XPS, fmt = "%.0f", delimiter = ",")
