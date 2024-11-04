import sys; sys.path.append("../src/")

import numpy as np
import pyFunctions as pyf

nt = 4001
dt = 5e-4

ns = 180
nr = 22

RPS = np.loadtxt("../inputs/geometry/migration_test_RPS.txt", delimiter = ",", dtype = float)
SPS = np.loadtxt("../inputs/geometry/migration_test_SPS.txt", delimiter = ",", dtype = float)
XPS = np.loadtxt("../inputs/geometry/migration_test_XPS.txt", delimiter = ",", dtype = int)

for s in range(ns):

    sx = SPS[s,0]
    sz = SPS[s,1]

    rx = RPS[XPS[s,1]:XPS[s,2], 0]
    rz = RPS[XPS[s,1]:XPS[s,2], 1]

    distance = np.sqrt((sx - rx)**2 + (sz - rz)**2)

    time = distance / 2000.0

    ts = np.array(time/dt + 0.15/dt, dtype = int)

    elastic = pyf.read_binary_matrix(nt, nr, f"../inputs/data/migration_test_elastic_iso_nStations{nr}_nSamples{nt}_shot_{s+1}.bin")
    
    for r in range(nr):
        elastic[:ts[r], r] = 0.0

    elastic.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/data/migration_test_elastic_iso_nStations{nr}_nSamples{nt}_shot_{s+1}.bin")
