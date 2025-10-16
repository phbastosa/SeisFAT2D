import sys; sys.path.append("../src/")

import numpy as np
import functions as pyf

SPS = np.loadtxt("../../anisoModels/BP2D/geometry/BP2D_streamer_SPS.txt", delimiter = ",", dtype = np.float32) 
RPS = np.loadtxt("../../anisoModels/BP2D/geometry/BP2D_streamer_RPS.txt", delimiter = ",", dtype = np.float32) 
XPS = np.loadtxt("../../anisoModels/BP2D/geometry/BP2D_streamer_XPS.txt", delimiter = ",", dtype = np.int32) 

nt = 6001
dt = 5e-4

ns = len(SPS)
nr = len(RPS)

spread = XPS[0,2]

for sId in range(ns):

    data_file = f"../../anisoModels/BP2D/acoustic/seismogram_nt6001_nr116_500us_shot_{sId+1}.bin"

    data = pyf.read_binary_matrix(nt, spread, data_file)

    data *= 1.0 / np.max(np.abs(data))

    travel_time = np.sqrt((SPS[sId,0] - RPS[XPS[sId,1]:XPS[sId,2],0])**2 + 
                          (SPS[sId,1] - RPS[XPS[sId,1]:XPS[sId,2],1])**2) / 1200

    tId = np.array((travel_time + 0.2) / dt, dtype = int)

    sigma = 50
    for rId in range(spread):
        if tId[rId] >= nt:
            tId[rId] = nt-1 

        data[:tId[rId], rId] *= np.exp(-((np.arange(tId[rId])-tId[rId]) / sigma)**2)        

    output_path = f"../inputs/data/seismogram_BP2D_shot_{sId+1}.bin"
    
    data.flatten("F").astype(np.float32, order = "F").tofile(output_path)

import matplotlib.pyplot as plt

fig, ax = plt.subplots()

ax.imshow(data, aspect = "auto", cmap = "Greys")

plt.show()