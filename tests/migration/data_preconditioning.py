import sys; sys.path.append("../src/")

import numpy as np
import pyFunctions as pyf

nt = 5001
dt = 1e-3

ns = 181
nr = 401

for s in range(ns):

    eikonal = pyf.read_binary_array(nr, f"../outputs/syntheticData/eikonal_iso_nStations{nr}_shot_{s+1}.bin")
    elastic = pyf.read_binary_matrix(nt, nr, f"../outputs/syntheticData/elastic_iso_nStations{nr}_nSamples{nt}_shot_{s+1}.bin")

    ts = np.array(eikonal/dt + 0.8/dt, dtype = int)

    for r in range(nr):
        elastic[:ts[r], r] = 0.0

    elastic.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/data/migration_test_data_nStations{nr}_nSamples{nt}_shot_{s+1}.bin")
