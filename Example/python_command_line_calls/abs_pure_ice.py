"""
Absorption coefficient of pure ice
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

l = flick.run("iop absorption_length 350e-9 700e-9 900 pure_ice")
wl = l[:,0]*1e9
a = 1/l[:,1]

fig, ax = plt.subplots() 
ax.semilogy(wl,a, linewidth=1)
ax.grid()
ax.set_xlabel("Wavelength [nm]")
ax.set_ylabel("Pure ice absorption coefficient [m$^{-1}$]")
plt.show();

