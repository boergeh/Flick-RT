import numpy as np
import matplotlib.pyplot as plt
import os

fig, ax = plt.subplots(2,2,sharex=True)
fig.tight_layout(pad=3)
fig.set_size_inches(7.5, 5)

os.system('flick text xy ./absorption_coefficient.txt > xy_tmp.txt')
xy = np.loadtxt("xy_tmp.txt")
ax[0,0].loglog(xy[:,0], xy[:,1],'b-')
ax[0,0].grid()
ax[0,0].set_ylabel("Absorption coefficient [m$^{-1}$]")

os.system('flick text xy ./refractive_index.txt > xy_tmp.txt')
xy = np.loadtxt("xy_tmp.txt")
ax[0,1].loglog(xy[:,0], xy[:,1],'b-')
ax[0,1].grid()
ax[0,1].set_ylabel("Refractive index")
ax[0,1].set_yticks([0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

os.system('flick text xy ./temperature_correction.txt > xy_tmp.txt')
xy = np.loadtxt("xy_tmp.txt")
ax[1,0].semilogx(xy[:,0], xy[:,1],'b-')
ax[1,0].grid()
ax[1,0].set_xlabel("Wavelength [m]")
ax[1,0].set_ylabel("Temperature correction [m$^{-1}$K$^{-1}$]")

os.system('flick text xy ./salinity_correction.txt > xy_tmp.txt')
xy = np.loadtxt("xy_tmp.txt")
ax[1,1].semilogx(xy[:,0], xy[:,1],'b-')
ax[1,1].grid()
ax[1,1].set_xlabel("Wavelength [m]")
ax[1,1].set_ylabel("Salinity correction [m$^{-1}$Lg$^{-1}$]")
ax[0,1].set_xticks([1e-8, 1e-6, 1e-4, 1e-2, 1e0, 1e2])

os.system('rm -f ./xy_tmp.txt')
plt.show()
#fig.savefig("plot_data.pdf", bbox_inches='tight')
