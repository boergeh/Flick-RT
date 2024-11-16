"""
Wavelength-dependent phase function of spheres
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

s = "450e-9 650e-9 3 pure_water 294 0"
scat_l = flick.run("iop scattering_length "+ s)
ref_idx = flick.run("iop refractive_index "+ s)

ds = []
for i in range(3):
    wl = ref_idx[i][0]
    idx = ref_idx[i][1]
    p = flick.run("mie "+str(wl)+" 1.0 "+str(idx)+" 800e-9 0.03 1 scattering_matrix_element 0 0 100")
    ds.append(p)

plt.style.use('dark_background')
fig, ax = plt.subplots()   
c = [[0,0,1],[0,1,0],[1,0,0]] 
for i in range(3):
    x = ds[i][:,0]*180/np.pi
    y = ds[i][:,1]
    ax.semilogy(x, y, color=c[i][:],linewidth=2)
ax.grid()
ax.set_xlabel("angle [degrees]")
ax.set_ylabel(r"differential scattering cross section [m$^2$$\,$sr$^{-1}$]")
plt.show();
