"""
Top-of-atmosphere upwelling radiance divided by downwelling
irradiance. Refer to the flick_tmp/config file (generated after the
first run) for documentation on all variables that can be configured
using the set function, as shown in this example. SI units and degrees
are used throughout unless stated otherwise.

"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.environ['FLICK_PATH'])
from python_script import flick

wl_grid = np.linspace(400e-9,700e-9,30);
solar_zenith_angle = 50
f = flick.toa_reflectance()
f.set_n_angles(50)
f.set("aerosol_od", 0)
f.set("cloud_liquid", 0)
f.set("ozone", 0.003)
    
R = f.spectrum(wl_grid, solar_zenith_angle)

plt.plot(R[:,0]*1e9,R[:,1])
plt.xlabel('Wavelength [nm]')
plt.ylabel('Top-of-atmosphere reflectance [sr$^{-1}$]')
plt.grid()

if __name__ == "__main__":
    plt.show()
