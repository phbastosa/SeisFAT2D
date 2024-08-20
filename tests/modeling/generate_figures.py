import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt
import pyFunctions as pyf

nd = 401

dcalc = pyf.read_binary_array(nd, "../outputs/syntheticData/travel_time_401_stations_shot_1.bin")

z = np.array([500, 1000, 2500])
v = np.array([2000, 2500, 3000, 4200])
x = np.linspace(50, 19950, nd)

direct_wave = x / v[0]
refractions = pyf.get_analytical_refractions(v,z,x)

dtrue = np.zeros_like(dcalc)

for i in range(nd):
    dtrue[i] = min(direct_wave[i], np.min(refractions[:,i]))
    

plt.plot(dtrue - dcalc)
plt.ylim([-0.003, 0.003])
plt.gca().invert_yaxis()
plt.show()