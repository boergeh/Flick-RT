"""
Sub-surface nadir radiance spectrum with a given spectral
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

time_point_utc = "2022 5 20 11 24 0"
latitude = "59.878675"
longitude = "5.655558"
wl_grid = np.linspace(300e-9,950e-9,20);
wl_width = 10e-9
f = flick.ocean_nadir_radiance()
f.set_n_angles(100)
f.set("aerosol_od", 0)
f.set("cloud_liquid", 0)
f.set("nap_concentration", 3e-3)
f.set("chl_concentration", 2e-6)
f.set("cdom_440", 0.1)

L = f.spectrum(wl_grid, wl_width, time_point_utc, latitude, longitude)
L = f.to_mW_per_m2_nm_sr(L)      

fig, ax = plt.subplots()
fig.set_size_inches(5,4)
ax.plot(L[:,0],L[:,1])
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel('Ocean nadir radiance [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
ax.grid()
plt.show()
fig.savefig("ocean_radiance.pdf", bbox_inches='tight')




