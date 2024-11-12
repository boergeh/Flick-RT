"""

"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
from matplotlib import cm
import flick

colors = ['r','g','b']

lms = ["L","M","S"]
c = [] 
for i in range(len(lms)):
    xy = flick.run("filter cone_"+lms[i])
    c.append(xy)
    
xyz = ["x","y","z"]
t = [] 
for i in range(len(lms)):
    xy = flick.run("filter tristimulus_"+xyz[i])
    t.append(xy)

fig, ax = plt.subplots(2,constrained_layout=True)
fig.set_size_inches(5,7)
for i in range(len(c)):
    x = c[i][:,0]*1e9
    y = c[i][:,1]*1e-9
    ax[0].plot(x, y,'-', linewidth=1.5, label=lms[i],color=colors[i])
ax[0].set_title("Spectral sensitivity of cone cells")


xyz = [r"$\bar{x}$",r"$\bar{y}$",r"$\bar{z}$"]
for i in range(len(c)):
    x = t[i][:,0]*1e9
    y = t[i][:,1]*1e-9
    ax[1].plot(x, y,'-', linewidth=1.5, label=xyz[i],color=colors[i])
ax[1].set_title("Tristimulus functions")

for i in range(2):
    ax[i].legend()
    ax[i].grid()
    ax[i].set_xlabel("Wavelength [nm]")
    ax[i].set_xlim([375,710])
plt.savefig('photometry.png',dpi=300)
plt.show();

