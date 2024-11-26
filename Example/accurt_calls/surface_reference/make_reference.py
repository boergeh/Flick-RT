"""
Flick solar radiation reference spectrum just above ocean surface in W/m2/nm
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+'/python_script')
import flick

# When Sun-Earth distance is one astronomical unit
time_point_1AU = 20240404000000

# Wavelengths used to compute transmittance of atmosphere
wl_low = 290e-9
wl_high = 1040e-9
n_wls = int((wl_high-wl_low)/5e-9)
wls = np.linspace(wl_low, wl_high, n_wls)

# Spectral resolution
wl_width = 10e-9

# Run Flick with AccuRT
f = flick.ocean_downward_plane_irradiance()
f.set('aerosol_od', 0.08)
f.set('cloud_liquid', 0)
f.set('detector_height', 0.1)
f.set('cdom_440', 0.1)
f.set('chl_concentration', 1e-6)
f.set('nap_concentration', 1e-3)
f.set_override_sun_zenith_angle(60)
E = f.spectrum(wls, wl_width, time_point_1AU)
E = f.to_W_per_m2_nm(E)

# Save to ascii file
file_name = 'surface_reference.txt' 
np.savetxt(file_name, E, fmt=['%6.2f ','%8.3e'])

# Insert header
with open(file_name, 'r') as file:
    original_content = file.read()
with open('header.txt', 'r') as file:
    header = file.read()
with open(file_name, 'w') as file:
    file.write(header + '\n' + original_content)

# Show spectrum
plt.plot(E[:,0],E[:,1])
plt.xlabel('Wavelength [nm]')
plt.ylabel(r'Irradiance [W$\,$m$^{-2}\,$nm$^{-1}$]')
plt.grid()
if __name__ == "__main__":
    plt.show()
