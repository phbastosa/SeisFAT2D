import numpy as np
import matplotlib.pyplot as plt

i = complex(0,1)

nt = 1001
dt = 1e-3
fm = 30.0
t0 = 0.15

fc = fm / (3.0 * np.sqrt(np.pi)) 

td = np.arange(nt) * dt - t0
arg = np.pi**3 * fc**2 * td**2

signal_original = (1.0 - 2.0*arg)*np.exp(-arg)

freq = np.fft.fftfreq(nt, dt)

fft_signal = np.fft.fft(signal_original)

w = 2.0 * np.pi * freq

fft_signal = np.sqrt(i*w) * fft_signal

signal_modified = (1/10) * np.real(np.fft.ifft(fft_signal))


plt.plot(td, signal_original)
plt.plot(td, signal_modified)

plt.show()
