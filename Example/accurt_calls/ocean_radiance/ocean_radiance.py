"""
Modeled and measured sub-surface nadir radiance spectrum with a
given spectral band width. See the flick_tmp/config file that will be
generated after the first run for documentation on all variables that
may be set with the 'set' function used in this example. SI-units and
degrees are used unless otherwise specified.
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

meta = flick.ocean_meta('HF22_D001_ocean_meta.txt')
Lm = flick.table("HF22_D001_ocean_radiance.txt")

wl_grid = np.linspace(300e-9,800e-9,20);
wl_width = 10e-9
f = flick.ocean_nadir_radiance()
f.set_n_angles(100)
f.set("aerosol_od", 0)
f.set("cloud_liquid", 0)
f.set("mp_names", "SD16_VF18") #To be replaced with Hardangerfjord sample
f.set("mp_concentrations", 20e-3)
f.set("detector_height", -0.5)

L = f.spectrum(wl_grid, wl_width, meta.time_point_utc, meta.latitude, meta.longitude)
L = f.to_mW_per_m2_nm_sr(L)

fig, ax = plt.subplots()
fig.set_size_inches(5,4)
line1 = ax.plot(L[:,0],L[:,1],label='Modeled with SD16 marine particles')
line2 = ax.plot(Lm[:,0],Lm[:,1],label='Measured in Hardangerfjorden')
ax.legend(loc='upper right')
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel('Ocean nadir radiance [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
ax.grid()
plt.show()
fig.savefig("ocean_radiance.pdf", bbox_inches='tight')




