"""
Plots absorption optical tickness of the atmosphere
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick


ot = flick.absorption_optical_thickness("ao_config",750e-9,800e-9,25).atmosphere()
fig, ax = plt.subplots(1,1)
fig.set_size_inches(6,5)
ax.semilogy(ot[:,0]*1e6,ot[:,1])
ax.grid()
ax.set_xlabel(r"Wavelength [$\mu$m]")
ax.set_ylabel("Absorption optical thickness")
plt.show()
fig.savefig("optical_thickness.png",dpi=300)








