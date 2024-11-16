"""
Single scattering albedo of snow
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

command = "iop absorption_length 300e-9 700e-9 900"
l = flick.run("iop absorption_length 300e-9 800e-9 900 pure_ice")
wl = l[:,0]*1e9
r = 1e-3
area = np.pi*r**2
vol = 4/3*np.pi*r**3
a = 1/l[:,1]*vol
b = 2*area
ssalb = b/(a+b)

fig, ax = plt.subplots() 
ax.semilogy(wl,a, linewidth=1)
ax.grid()
ax.set_xlabel("Wavelength [nm]")
ax.set_ylabel("Single scattering albedo of snow")
plt.show();
