"""
Make the Flick logo
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import flick
from matplotlib import cm
from matplotlib.colors import LightSource

recalculate = True
black_background = False 

f = "config"
if not os.path.exists(f):
    flick.run("accurt -g toa_reflectance "+f)
        
n_angles = 100
n_polar = 50;
n_azimuth = 150
flick.config(f,"stream_upper_slab_size", flick.to_streams(n_angles))    
flick.config(f,"detector_wavelengths", 300e-9)
flick.config(f,"n_angles", n_angles)
flick.config(f,"n_heights", 8)
flick.config(f,"detector_height", -1)
flick.config(f,"detector_orientation", "up")
flick.config(f,"source_zenith_angle", 0)
flick.config(f,"aerosol_od", 0)
flick.config(f,"cloud_liquid", 0)
flick.config(f,"nap_concentration", 0)
flick.config(f,"chl_concentration", 0)
flick.config(f,"cdom_440", 0)
flick.config(f,"bottom_depth", 10000)
flick.config(f,"detector_radiance_distribution_override",
             np.array([n_polar, n_azimuth]))

if (recalculate):
    print("Using Flick to find UVR distributon at 1 m depth ...")
    r = flick.run("accurt "+f)
    np.save("logo",r)
else:
    r = np.load("logo.npy")

scaling_factor = 1.0    
theta, phi = np.linspace(0, np.pi, n_polar),np.linspace(-np.pi, np.pi, n_azimuth)
THETA, PHI = np.meshgrid(theta, phi)
R = r.transpose()/r.max()*scaling_factor
X = R * np.sin(THETA) * np.cos(PHI)
Y = R * np.sin(THETA) * np.sin(PHI)
Z = R * np.cos(THETA)

if black_background:
    plt.style.use('dark_background')
ls = LightSource(azdeg=0, altdeg=55)
rgb = ls.shade(R, cmap=cm.hot, vert_exag=0, blend_mode='soft')
fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=rgb,
                       linewidth=0, antialiased=False, shade=False,alpha=1)

font = {'fontname':'Arial'}
ax.text(0,0.77,-0.35,"flick",fontsize=180,fontweight='heavy',**font)

ax.plot(0,1.94,0.67,'o',color='black',markerfacecolor=[0.81,0,0],markeredgewidth=0,
        markersize=35,zorder=10)
ax.view_init(elev=0, azim=0, roll=0)


ax.set_aspect('equal')
fig.set_size_inches(9,3.5,forward=True)
plt.subplots_adjust(left=-0.3,
                    bottom=-0.3,
                    right=0.9,
                    top=1.28)
plt.axis('off')
plt.show()

if black_background:
    fig.savefig("logo_black.png",dpi=500)
    fig.savefig("logo_black.jpg",dpi=500)
else:
    fig.savefig("logo_white.png",dpi=500)
    fig.savefig("logo_white.jpg",dpi=500)











