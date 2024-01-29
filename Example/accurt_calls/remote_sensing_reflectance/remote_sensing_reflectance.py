"""
Remote sensing reflectance spectrum with a given spectral
band width. See the flick_tmp/config file, which will be generated
after the first run, for documentation on all variables that may be
set with the set function used in this example. SI-units and degrees
are used unless otherwise specified.
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

wl_grid = np.linspace(300e-9,950e-9,20);
f = flick.remote_sensing_reflectance()

f.set_n_angles(90)
f.set("aerosol_od", 0)
f.set("cloud_liquid", 0)
f.set("nap_concentration", 3e-3)
f.set("chl_concentration", 2e-6)
f.set("cdom_440", 0.1)
    
Rsr = f.spectrum(wl_grid)
Rsr[:,0] = Rsr[:,0]*1e9 # to nm

fig, ax = plt.subplots()
fig.set_size_inches(5,4)
ax.plot(Rsr[:,0],Rsr[:,1])
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel('Remote sensing reflectance [sr$^{-1}$]')
ax.grid()
plt.subplots_adjust(left=0.17, bottom=0.12, right=0.96, top=0.96,
                    wspace=0, hspace=0)
plt.show()
fig.savefig("remote_sensing_reflectance.pdf", bbox_inches='tight')








