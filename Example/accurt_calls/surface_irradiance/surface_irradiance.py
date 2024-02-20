"""
Surface downward irradiance spectrum with a given spectral
band width. See the flick_tmp/config file, which will be generated
after the first run, for documentation on all variables that may be
set with the set function used in this example. SI-units and degrees
are used unless otherwise specified.
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

time_point_utc = "2022 5 20 11 24 0"
latitude = "59.878675"
longitude = "5.655558"
wl_grid = np.linspace(280e-9,950e-9,30);
wl_width = 10e-9
f = flick.surface_irradiance()

f.set_n_angles(50)
f.set("aerosol_od", 0)
f.set("cloud_liquid", 0)
f.set("ozone", 0.003)
f.set("water_vapor", 10)
    
E = f.spectrum(wl_grid, wl_width, time_point_utc, latitude, longitude)
E = f.to_W_per_m2_nm(E)      

fig, ax = plt.subplots()
fig.set_size_inches(5,4)
ax.plot(E[:,0],E[:,1])
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel('Surface irradiance [W m$^{-2}$ nm$^{-1}$]')
ax.grid()
plt.show()
fig.savefig("surface_irradiance.pdf", bbox_inches='tight')
