"""
Shows wavelength-dependent phase function
"""
import numpy as np
import matplotlib.pyplot as plt
import os

s = "450e-9 650e-9 3 pure_water 294 0 > tmp.txt"
os.system("flick iop scattering_length "+s)
scat_l = np.loadtxt("tmp.txt")

os.system("flick iop refractive_index "+s)
ref_idx = np.loadtxt("tmp.txt")

ds = []
for i in range(3):
    wl = ref_idx[i][0]
    idx = ref_idx[i][1]
    os.system("flick mie "+str(wl)+" 1.0 "+str(idx)+" 800e-9 0.03 1 scattering_matrix_element 0 0 100 > tmp.txt")
    p = np.loadtxt("tmp.txt")
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
ax.set_ylabel("differential scattering cross section [m$^2$sr$^{-1}$]")
plt.show();
#fig.savefig("material_phase_function.pdf", bbox_inches='tight')
