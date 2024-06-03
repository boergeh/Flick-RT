"""
Compute upward nadir radiance and downward irradiance spectra in the
atmosphere-ocean system. See the flick_tmp/config file that will be
generated after the first run for documentation on all variables that
may be set with the 'set' function used in this script. SI-units (mks)
and degrees are used unless otherwise specified. Two text files
containing computed radiance and irradiance spectra will be created
and stored in the output directory.

User generated data files may be added to the input directory. 

Default settings gives output radiation spectra inside the water
column with negative irradiance_height, but spectra from any height in
the atmosphere and water column may be generated. For example can
top-of-atmosphere radiation be computed by setting irradiance_height =
120e3.

Generated radiance and irradiance spectra can be viewed afterwards by
running plot_radiation.py.

Remote sensing reflectance can be obtained by setting the
irradiance_height just above the water surface, and set it higher than
the detector separation to ensure atmospheric radiance, and to remove
surface specular reflections, set
f.set('subtract_specular_radiance','true'). The remote sensing
reflectance can the be viewed afterwards with plot_radiation_ratio.py

The n_wavelengths variable gives number of wavelength used to generate
the transmittance from top of the atmosphere to the position of the
radiation output. A high resolution solar spectrum is then multiplied
in and the resulting spectrum is smoothed with a triangular function
with a given band_width.

"""
import numpy as np
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+'/python_script')
import flick

# - begin user input -
station = 'ECOSENS_HF22_D1'
irradiance_height = -0.3
detector_separation = 0.37
first_wavelength = 320e-9
last_wavelength = 850e-9
n_wavelengths = 60
band_width = 10e-9
compute_at_satellite_wavelengths = False
# - end user input -

meta = flick.ocean_meta('input/'+station+'_meta.txt')

def radiation(f, detector_height, wavelengths):
    # - begin user input -
    f.set('aerosol_od', 0.28)
    f.set('aerosol_ratio', 1)
    f.set('relative_humidity', 0.5)
    f.set('cloud_liquid', 0)
    f.set('ozone', 0.004)
    f.set('pressure', 1000e2)
    f.set('temperature', 273+15)
    f.set('water_vapor',30)
    f.set('mp_names', './input/'+station) 
    f.set('mp_concentrations', meta.spm)
    f.set('mcdom_names', './input/'+station)
    f.set('mcdom_scaling_factors', 1.0)
    f.set('subtract_specular_radiance','false')
    # - end user input -
    f.set('detector_height', detector_height)
    return f.spectrum(wavelengths, band_width, meta.time_point_utc,
                      meta.latitude, meta.longitude)

wl = np.linspace(first_wavelength, last_wavelength, n_wavelengths)


if compute_at_satellite_wavelengths:
    instrument_wl = flick.table('input/'+station+'_toa_radiance.txt')[:,0]
else:
    instrument_wl = flick.table('input/'+station+'_ocean_irradiance.txt')[:,0]
wl = flick.move_closest_values_in_a_to_those_in_b(wl,instrument_wl)

height = irradiance_height
file_name = 'output/'+station+'_computed_irradiance_W_per_m2_nm.txt'
f = flick.ocean_downward_plane_irradiance()
E = radiation(f, height, wl)
E = f.to_W_per_m2_nm(E)
np.savetxt(file_name, E, fmt=['%6.2f ','%8.3e'])


if compute_at_satellite_wavelengths:
    instrument_wl = flick.table('input/'+station+'_toa_radiance.txt')[:,0]
else:
    instrument_wl = flick.table('input/'+station+'_ocean_radiance.txt')[:,0]

wl = flick.move_closest_values_in_a_to_those_in_b(wl,instrument_wl)
height = irradiance_height-detector_separation
file_name = 'output/'+station+'_computed_radiance_mW_per_m2_nm_sr.txt'
f = flick.ocean_nadir_radiance()
L = radiation(f, height, wl)
L = f.to_mW_per_m2_nm_sr(L)
np.savetxt(file_name, L, fmt=['%6.2f ','%8.3e'])
