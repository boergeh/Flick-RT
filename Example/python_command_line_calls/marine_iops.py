"""
Plot inherent optical properties of ocean materials. SI-units and
degrees are used unless otherwise specified.

"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

station = "ECOSENS_HF22_D1"
from_wl = 300e-9
to_wl = 1050e-9
wl_grid = np.linspace(from_wl,to_wl,300);
wl_width = 10e-9

meta_path = "../../material/marine_particles/iop_tables"
meta = flick.ocean_meta(meta_path+"/"+station+"_meta.txt")
iops = flick.marine_iops(station,meta.spm,from_wl,to_wl,200)
iops.set_b_scaling_factor(1)
a_w = iops.a_water()
a_spm = iops.a_spm()
a_cdom = iops.a_cdom()

fig, ax = plt.subplots(2,2)
fig.set_size_inches(7.5,6)
ax[0,0].plot(a_w[:,0]*1e9,a_w[:,1],'-',label="water")
ax[0,0].plot(a_spm[:,0]*1e9,a_spm[:,1],'-',label="SPM")
ax[0,0].plot(a_cdom[:,0]*1e9,a_cdom[:,1],'-',label="CDOM")
ax[0,0].legend()
ax[0,0].grid()
ax[0,0].set_xlabel('Wavelength [nm]')
ax[0,0].set_ylabel('Absorption coefficient [m$^{-1}$]')
ymax = np.max(a_cdom[:,1])*0.3
ax[0,0].set_ylim([-0.05*ymax,1.05*ymax]) 

b_w = iops.b_water()
b_spm = iops.b_spm()
ax[1,0].plot(b_w[:,0]*1e9,b_w[:,1],'-',label="water")
ax[1,0].plot(b_spm[:,0]*1e9,b_spm[:,1],'-',label="SPM")
ax[1,0].legend()
ax[1,0].grid()
ax[1,0].set_xlabel('Wavelength [nm]')
ax[1,0].set_ylabel('Scattering coefficient [m$^{-1}$]')

p = iops.volume_scattering_function(0)
ax[0,1].semilogy(p[:,0],p[:,1],label="measured")    
terms = [4,16,64]
for i in range(len(terms)):
    n = terms[i]
    g = iops.asymmetry_factor(n)
    f = iops.volume_scattering_scaling_factor(n)
    bb = iops.back_scattering_coefficient(n)
    p = iops.volume_scattering_function(n)
    l = "n="+str(n)+", g="+"{:.2f}".format(g)+", f="+"{:.2f}".format(f)+", b$_b$="+"{:.2f}".format(bb)
    ax[0,1].semilogy(p[:,0],p[:,1],label=l)    
ax[0,1].grid()
ax[0,1].legend()
ax[0,1].set_xlabel(r'Scattering angle [$\degree$]')
ax[0,1].set_ylabel('VSF at 515 nm [m$^{-1}$sr$^{-1}$]')
ax[0,1].set_xticks([0,45,90,135,180])

ax[1,1].axis("off")
ax[1,1].text(0.0,0.0,station+" - "+meta.to_string())

plt.subplots_adjust(left=0.09, bottom=0.08, right=0.98, top=0.98,
                    wspace=0.3, hspace=0.2)
plt.show()
fig.savefig("marine_iops_"+station+".png",dpi=300)








