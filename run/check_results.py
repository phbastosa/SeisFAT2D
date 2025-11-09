import sys; sys.path.append("../src/")

import numpy as np
import matplotlib.pyplot as plt

import functions as pyf

parameters = str(sys.argv[1])

nx = int(pyf.catch_parameter(parameters, "x_samples"))
nz = int(pyf.catch_parameter(parameters, "z_samples"))

dx = float(pyf.catch_parameter(parameters, "x_spacing"))
dz = float(pyf.catch_parameter(parameters, "z_spacing")) 

nt = int(pyf.catch_parameter(parameters, "time_samples")) 
dt = float(pyf.catch_parameter(parameters, "time_spacing"))

max_angle = float(pyf.catch_parameter(parameters, "mig_max_angle"))

da = float(pyf.catch_parameter(parameters, "mig_angle_spacing"))

nang = int(max_angle / da) + 1

sps_path = pyf.catch_parameter(parameters, "SPS")
rps_path = pyf.catch_parameter(parameters, "RPS")
xps_path = pyf.catch_parameter(parameters, "XPS")

SPS = np.loadtxt(sps_path, delimiter = ",", dtype = np.float32) 
RPS = np.loadtxt(rps_path, delimiter = ",", dtype = np.float32) 
XPS = np.loadtxt(xps_path, delimiter = ",", dtype = np.int32) 

ds = abs(SPS[1,0] - SPS[0,0])
dr = abs(RPS[1,0] - RPS[0,0])

nr = XPS[0,2]
ns = len(XPS)

offset = SPS[0,0] - RPS[XPS[0,1]:XPS[0,2],0]

offset_min = np.min(offset)
offset_max = np.max(offset)

#----------------------------------------------------------------------------------

vp_path = pyf.catch_parameter(parameters, "vp_model_file")

vp = pyf.read_binary_matrix(nz, nx, vp_path)

fig, ax = plt.subplots(figsize = (17,5))

im = ax.imshow(vp, aspect = "auto", cmap = "jet", extent = [0,(nx-1)*dx, (nz-1)*dz, 0], vmin = 1400, vmax = 2600)

cbar = plt.colorbar(im, ax = ax, pad = 0.01)
cbar.set_label("P wave velocity [m/s]", fontsize = 15)

ax.set_xlabel("Distance [m]", fontsize = 15)
ax.set_ylabel("Depth [m]", fontsize = 15)

ax.plot(RPS[XPS[0,1]:XPS[0,2],0]+20, RPS[XPS[0,1]:XPS[0,2],1]+50, "vr", markersize = 5)
ax.plot(SPS[0,0]+20, SPS[0,1]+50, "ok", markersize = 5)

ax.plot(RPS[XPS[-1,1]:XPS[-1,2],0]-20, RPS[XPS[-1,1]:XPS[-1,2],1]+50, "vr", markersize = 5)
ax.plot(SPS[-1,0]-20, SPS[-1,1]+50, "ok", markersize = 5)

ax.set_xlim([0,(nx-1)*dx])

fig.tight_layout()
plt.show()

dCMP = 0.5*min(ds,dr)

cmp_beg = 0.5*(SPS[0,0] + RPS[0,0])
cmp_end = 0.5*(SPS[-1,0] + RPS[-1,0])

nCMP = int((cmp_end - cmp_beg) / dCMP) + 1

fold = np.zeros(nCMP, dtype = int)
xCMP = np.arange(cmp_beg, cmp_end + dCMP, dCMP)

for sId in range(len(XPS)):

    spread = RPS[int(XPS[sId,1]):int(XPS[sId,2]),0]

    for rId, xrec in enumerate(spread):

        midpoint = int(rId + 2.0*(ds/dr)*sId)

        fold[midpoint] += 1

fig, ax = plt.subplots(figsize = (10,5))

ax.plot(xCMP, fold)
ax.set_xlim([0,(nx-1)*dx])
ax.set_ylim([0,np.max(fold)])

ax.set_xlabel("Distance [m]", fontsize = 15)
ax.set_ylabel("Receivers per CMP position", fontsize = 15)
ax.set_title("Geometry Coverage", fontsize = 18)

fig.tight_layout()
plt.show()

#----------------------------------------------------------------------------------

sId = 100

raw_path = f"../../FWI2D/outputs/data/seismogram_nt{nt}_nr{nr}_{dt*1e6:.0f}us_shot_{sId+1}.bin"
mig_path = f"../inputs/data/seismogram_input_shot_{sId+1}.bin"

gather_raw = pyf.read_binary_matrix(nt, nr, raw_path)
gather_mig = pyf.read_binary_matrix(nt, nr, mig_path)

gather_raw *= 100.0 / np.max(np.abs(gather_raw))
gather_mig *= 100.0

scale = 1.0

fig, ax = plt.subplots(ncols = 2, figsize = (10,8))

im = ax[0].imshow(np.flip(gather_raw, axis = 1), aspect = "auto", cmap = "Greys", extent = [offset_min, offset_max, (nt-1)*dt, 0], vmin = -scale, vmax = scale)
cbar = plt.colorbar(im, ax = ax[0], pad = 0.02, shrink = 0.5)
cbar.set_label("Raw data amplitude", fontsize = 15)

ax[0].set_title("Raw data", fontsize = 18)
ax[0].set_xlabel("Offset [m]", fontsize = 15)
ax[0].set_ylabel("Two Way Time [s]", fontsize = 15)

im = ax[1].imshow(np.flip(gather_mig, axis = 1), aspect = "auto", cmap = "Greys", extent = [offset_min, offset_max, (nt-1)*dt, 0], vmin = -scale, vmax = scale)
cbar = plt.colorbar(im, ax = ax[1], pad = 0.02, shrink = 0.5)
cbar.set_label("Mig data amplitude", fontsize = 15)

ax[1].set_title("Mig data", fontsize = 18)
ax[1].set_xlabel("Offset [m]", fontsize = 15)
ax[1].set_ylabel("Two Way Time [s]", fontsize = 15)

fig.tight_layout()
plt.show()

#----------------------------------------------------------------------------------

IDKDM_path = f"../outputs/seismic/IDKDM_result_{nz}x{nx}.bin"
ADKDM_path = f"../outputs/seismic/ADKDM_result_{nz}x{nCMP}x{nang}.bin"

IDKDM = pyf.read_binary_matrix(nz, nx, IDKDM_path)
ADKDM = pyf.read_binary_volume(nz, nang, nCMP, ADKDM_path)

cmpId = 360

fig, ax = plt.subplots(figsize = (17,5))

im = ax.imshow(IDKDM, aspect = "auto", cmap = "Greys", extent = [0,(nx-1)*dx, (nz-1)*dz, 0])
cbar = plt.colorbar(im, ax = ax, pad = 0.01)
cbar.set_label("Amplitude", fontsize = 15)
ax.set_xlabel("Distance [m]", fontsize = 15)
ax.set_ylabel("Depth [m]", fontsize = 15)
ax.set_title("IDKDM", fontsize = 18)

fig.tight_layout()
plt.show()

fig, ax = plt.subplots(ncols = 5, figsize = (17,5))

for i in range(len(ax)):

    index = cmpId + i - int(0.5*len(ax))

    im = ax[i].imshow(ADKDM[:,:,index], aspect = "auto", cmap = "Greys", extent = [0, max_angle, (nz-1)*dz, 0])
    cbar = plt.colorbar(im, ax = ax[i], pad = 0.02, shrink = 0.5)
    cbar.set_label("Amplitude", fontsize = 15)
    ax[i].set_xlabel("Angle [°]", fontsize = 15)
    ax[i].set_ylabel("Depth [m]", fontsize = 15)
    ax[i].set_title(f"xCMP = {xCMP[index]} m", fontsize = 15)

fig.tight_layout()
plt.show()

#----------------------------------------------------------------------------------

max_it = int(pyf.catch_parameter(parameters, "mig_max_iteration"))

IDLSKDM_path = f"../outputs/seismic/IDLSKDM_result_{nz}x{nx}_iteration_{max_it}.bin"
ADLSKDM_path = f"../outputs/seismic/ADLSKDM_result_{nz}x{nCMP}x{nang}_iteration_{max_it}.bin"

IDLSKDM = pyf.read_binary_matrix(nz, nx, IDLSKDM_path)
ADLSKDM = pyf.read_binary_volume(nz, nang, nCMP, ADLSKDM_path)

fig, ax = plt.subplots(figsize = (17,5))

im = ax.imshow(IDLSKDM, aspect = "auto", cmap = "Greys", extent = [0,(nx-1)*dx, (nz-1)*dz, 0])
cbar = plt.colorbar(im, ax = ax, pad = 0.01)
cbar.set_label("Amplitude", fontsize = 15)
ax.set_xlabel("Distance [m]", fontsize = 15)
ax.set_ylabel("Depth [m]", fontsize = 15)
ax.set_title("IDLSKDM", fontsize = 18)

fig.tight_layout()
plt.show()

fig, ax = plt.subplots(ncols = 5, figsize = (17,5))

for i in range(len(ax)):

    index = cmpId + i - int(0.5*len(ax))

    im = ax[i].imshow(ADLSKDM[:,:,index], aspect = "auto", cmap = "Greys", extent = [0, max_angle, (nz-1)*dz, 0])
    cbar = plt.colorbar(im, ax = ax[i], pad = 0.02, shrink = 0.5)
    cbar.set_label("Amplitude", fontsize = 15)
    ax[i].set_xlabel("Angle [°]", fontsize = 15)
    ax[i].set_ylabel("Depth [m]", fontsize = 15)
    ax[i].set_title(f"xCMP = {xCMP[index]} m", fontsize = 15)

fig.tight_layout()
plt.show()

#----------------------------------------------------------------------------------

IDLSKDM_convergence = np.loadtxt(f"../outputs/residuo/IDLSKDM_convergence_{max_it}_iterations.txt")
ADLSKDM_convergence = np.loadtxt(f"../outputs/residuo/ADLSKDM_convergence_{max_it}_iterations.txt")

IDLSKDM_convergence *= 1.0 / np.max(ADLSKDM_convergence)
ADLSKDM_convergence *= 1.0 / np.max(ADLSKDM_convergence)

fig, ax = plt.subplots(figsize = (10,5))

ax.plot(IDLSKDM_convergence, "o--", label = "IDLSKDM")
ax.plot(ADLSKDM_convergence, "o--", label = "ADLSKDM")

ax.legend(loc = "upper right", fontsize = 15)

ax.set_title("Convergence curve", fontsize = 18)
ax.set_xlabel("Iteration number", fontsize = 15)
ax.set_ylabel(r"$|Lm - d|^2_2$ normalized")

fig.tight_layout()
plt.show()

