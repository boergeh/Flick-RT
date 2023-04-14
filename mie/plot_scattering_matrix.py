import numpy as np
import matplotlib.pyplot as plt
import os

fig, ax = plt.subplots(2,2,sharex=False)
fig.tight_layout(pad=3)
fig.set_size_inches(7.5, 5)

s1 = "flick mie 550e-9 1.0 1.33+0i 0.25e-6 0.0 0.1 scattering_matrix_element "
s2 = " 200 > xy_tmp.txt"

os.system(s1+"0 0"+s2)
xy = np.loadtxt("xy_tmp.txt")
ang = xy[:,0]*180/np.pi
F11 = xy[:,1];
ax[0,0].semilogy(ang, F11,'b-')
ax[0,0].grid()
ax[0,0].set_ylabel("F$_{11}$ [m$^{-2}$ sr$^{-1}$]")

os.system(s1+"0 1"+s2)
xy = np.loadtxt("xy_tmp.txt")
ax[0,1].plot(ang, xy[:,1]/F11,'b-')
ax[0,1].grid()
ax[0,1].set_ylabel("F$_{12}$/F$_{11}$")

os.system(s1+"2 2"+s2)
xy = np.loadtxt("xy_tmp.txt")
ax[1,0].plot(ang, xy[:,1]/F11,'b-')
ax[1,0].grid()
ax[1,0].set_ylabel("F$_{33}$/F$_{11}$")

os.system(s1+"3 2"+s2)
xy = np.loadtxt("xy_tmp.txt")
ax[1,1].plot(ang, xy[:,1]/F11,'b-')
ax[1,1].grid()
ax[1,1].set_ylabel("F$_{43}$/F$_{11}$")

for i in range(2):
    for j in range(2):
        ax[i,j].set_xticks([0, 45, 90, 135, 180])
        ax[i,j].set_xlabel("Scattering angle")

plt.show()
#fig.savefig("plot_data.pdf", bbox_inches='tight')
