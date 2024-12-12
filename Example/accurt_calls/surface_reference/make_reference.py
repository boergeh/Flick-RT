"""
Flick solar radiation reference spectrum just above ocean surface in
W/m2/nm.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
import sys
sys.path.append(os.environ['FLICK_PATH']+'/python_script')
import flick

path = os.environ['FLICK_PATH']+"/Example/accurt_calls/surface_reference"
os.chdir(path)

full_spectrum = False
save_data = False
overwrite_current_standard = False

# When Sun-Earth distance is one astronomical unit
time_point_1AU = 20240404000000

# Wavelengths used to compute transmittance of atmosphere
if full_spectrum:
    wl_low = 290e-9
    wl_high = 1040e-9
    n_wls = 230
else:
    wl_low = 400e-9
    wl_high = 700e-9
    n_wls = 30

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
wls = flick.atmosphere_wavelengths(wl_low, wl_high, n_wls)
E = f.spectrum(wls, wl_width, time_point_1AU)
E = f.to_W_per_m2_nm(E)

# Current flick reference
Er = flick.run('radiator surface_reference')
Er = f.to_W_per_m2_nm(Er)

if save_data:
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

if overwrite_current_standard:
    shutil.copy('surface_reference.txt','../../../radiator/.')

# Show spectrum
plt.plot(E[:,0],E[:,1],label='Suggested reference')
plt.plot(Er[:,0],Er[:,1],label='Current reference')
plt.xlabel('Wavelength [nm]')
plt.ylabel(r'Irradiance [W$\,$m$^{-2}\,$nm$^{-1}$]')
plt.grid()
plt.legend()
if __name__ == "__main__":
    plt.show()
