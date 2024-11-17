import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

fig, ax = plt.subplots(2,2,sharex=False)
fig.tight_layout(pad=3)
fig.set_size_inches(7.5, 5)

s1 = "mie 550e-9 1.0 1.33+0i 0.4e-6 0.0 0.1 scattering_matrix_element "
n_angles = 128

xy = flick.run(s1+"0 0 "+str(n_angles))
ang = xy[:,0]*180/np.pi
F11 = xy[:,1];
ax[0,0].semilogy(ang, F11,'b-')
ax[0,0].grid()
ax[0,0].set_ylabel("F$_{11}$ [m$^{-2}$ sr$^{-1}$]")

xy = flick.run(s1+"0 1 "+str(n_angles))
ax[0,1].plot(ang, xy[:,1]/F11,'b-')
ax[0,1].grid()
ax[0,1].set_ylabel("F$_{12}$/F$_{11}$")

xy = flick.run(s1+"2 2 "+str(n_angles))
ax[1,0].plot(ang, xy[:,1]/F11,'b-')
ax[1,0].grid()
ax[1,0].set_ylabel("F$_{33}$/F$_{11}$")

xy = flick.run(s1+"3 2 "+str(n_angles))
ax[1,1].plot(ang, xy[:,1]/F11,'b-')
ax[1,1].grid()
ax[1,1].set_ylabel("F$_{43}$/F$_{11}$")

for i in range(2):
    for j in range(2):
        ax[i,j].set_xticks([0, 45, 90, 135, 180])
        ax[i,j].set_xlabel("Scattering angle")

if __name__ == "__main__":
    plt.show()

