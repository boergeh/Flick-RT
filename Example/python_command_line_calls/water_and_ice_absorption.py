"""
Shows the pure water and ice absorption coefficient
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
from matplotlib import cm
import flick

sub_command = "iop absorption_length 300e-9 700e-9 900"
l_w = flick.run(sub_command+" pure_water 273.15 35")
l_i = flick.run(sub_command+" pure_ice")
wl_w = l_w[:,0]*1e9
a_w = 1/l_w[:,1]
wl_i = l_i[:,0]*1e9
a_i = 1/l_i[:,1]

fig, ax = plt.subplots() 
ax.semilogy(wl_w,a_w, linewidth=2,label="pure water")
ax.semilogy(wl_i,a_i, linewidth=2,label="pure ice")
ax.legend()
ax.grid()
ax.set_xlabel("Wavelength [nm]")
ax.set_ylabel("Absorption coefficient [m$^{-1}$]")
ax.set_ylim([4e-4, 1])
plt.show();
fig.savefig("water_and_ice_absorption.png", dpi=300)
