import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

colors = ['r','g','b']
lms = ["L","M","S"]
data1 = [] 
for i in range(len(lms)):
    xy = flick.run("filter cone_"+lms[i])
    data1.append(xy)
    
xyz = ["x","y","z"]
data2 = [] 
for i in range(len(lms)):
    xy = flick.run("filter "+xyz[i]+"_bar")
    data2.append(xy)

data3 = flick.run("radiator cie_d65")
data4 = flick.run("radiator cie_a")

fig, ax = plt.subplots(2,2,constrained_layout=True)
fig.set_size_inches(7,7)
for i in range(len(data1)):
    x = data1[i][:,0]*1e9
    y = data1[i][:,1]
    ax[0][0].plot(x, y,'-', linewidth=1.5, label=lms[i],color=colors[i])
ax[0][0].set_title("Spectral sensitivity of cone cells",fontsize=9)
            
xyz = [r"$\bar{x}$",r"$\bar{y}$",r"$\bar{z}$"]
for i in range(len(data2)):
    x = data2[i][:,0]*1e9
    y = data2[i][:,1]
    ax[0][1].plot(x, y,'-', linewidth=1.5, label=xyz[i],color=colors[i])
ax[0][1].set_title(r"CIE color matching functions, 1931 2$\degree$",fontsize=9)

ax[1][0].plot(data3[:,0]*1e9,data3[:,1],label="CIE standard illuminant D65")
xyz = flick.chromaticity(data3)    
rgb = flick.rgb(data3) 
ax[1][0].set_facecolor(rgb)
ax[1][0].set_title('xy = '+str(np.round(xyz[0:2],3))+ \
                   ',  RGB = '+str(np.round(rgb,3)),fontsize=9)
           
ax[1][1].plot(data4[:,0]*1e9,data4[:,1],'-',
              label="CIE standard illuminant A")
xyz = flick.chromaticity(data4)
rgb = flick.rgb(data4) 
ax[1][1].set_facecolor(rgb)
ax[1][1].set_title('xy = '+str(np.round(xyz[0:2],3)) + \
                    ',  RGB = '+str(np.round(rgb,3)),fontsize=9)
for i in range(2):
    for j in range(2):
        ax[i][j].legend()
        ax[i][j].grid()
        ax[i][j].set_xlabel("Wavelength [nm]")
        ax[i][j].set_xlim([375,730])

if __name__ == "__main__":
    plt.show()
