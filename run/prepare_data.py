import numpy as np

nt = 10001
dt = 5e-4
nr = 340

SPS = np.loadtxt("../inputs/geometry/WASMEM2D_SPS.txt", dtype = float, delimiter = ",")
RPS = np.loadtxt("../inputs/geometry/WASMEM2D_RPS.txt", dtype = float, delimiter = ",")
XPS = np.loadtxt("../inputs/geometry/WASMEM2D_XPS.txt", dtype = int, delimiter = ",")

vwb = 1400.0

for sId in range(len(SPS)):

    data_path = f"../inputs/data/elastic_iso_nStations340_nSamples10001_shot_{sId+1}.bin"

    seismic = np.fromfile(data_path, dtype = np.float32, count = nt*nr).reshape([nt,nr], order = "F")

    seismic *= 1.0 / np.max(np.abs(seismic))

    offset = np.sqrt((SPS[sId,0] - RPS[XPS[sId,1]:XPS[sId,2],0])**2+
                     (SPS[sId,1] - RPS[XPS[sId,1]:XPS[sId,2],1])**2)  

    twb = offset / vwb + 0.05

    for i in range(nr):
        seismic[:int(twb[i]/dt),i] = 0.0

    seismic.flatten("F").astype(np.float32, order = "F").tofile(f"../inputs/data/SeisFAT2D_kdm_ISO_input_shot_{sId+1}.bin")
