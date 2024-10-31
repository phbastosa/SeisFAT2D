import sys; sys.path.append("../src/")

import numpy as np
import pyFunctions as pyf

import matplotlib.pyplot as plt

nt = 5001
dt = 1e-3

ns = 181
nr = 401

XPS = np.loadtxt("../inputs/geometry/migration_test_XPS.txt", delimiter = ",", dtype = int)
RPS = np.loadtxt("../inputs/geometry/migration_test_RPS.txt", delimiter = ",", dtype = float)
SPS = np.loadtxt("../inputs/geometry/migration_test_SPS.txt", delimiter = ",", dtype = float)

for s in range(ns):

    sx = SPS[s, 0]
    sz = SPS[s, 1]

    rx = RPS[XPS[s,1]:XPS[s,2], 0]
    rz = RPS[XPS[s,1]:XPS[s,2], 1]

    distance = np.sqrt((sx - rx)**2 + (sz - rz)**2)

    time = distance / 1500.0

    ts = np.array(time/dt + 0.5/dt, dtype = int)

    elastic = pyf.read_binary_matrix(nt, nr, f"../outputs/syntheticData/elastic_iso_nStations{nr}_nSamples{nt}_shot_{s+1}.bin")
    
    for r in range(nr):
        elastic[:ts[r], r] = 0.0

    elastic.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/data/migration_test_data_nStations{nr}_nSamples{nt}_shot_{s+1}.bin")
