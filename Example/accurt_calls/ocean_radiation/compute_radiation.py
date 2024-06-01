"""
Computes upward nadir radiance and downward irradiance
spectra. See the flick_tmp/config file that will be generated after
the first run for documentation on all variables that may be set with
the 'set' function used in this script. SI-units and degrees are used
unless otherwise specified. Two text files containing the computed
spectra will be created.

"""
import numpy as np
#import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+'/python_script')
import flick

station = 'ECOSENS_HF22_D1'
irradiance_height = -0.3
detector_separation = 0.37
first_wavelength = 320e-9
last_wavelength = 850e-9
n_wavelengths = 60
band_width = 10e-9
meta = flick.ocean_meta(station+"_meta.txt")

def radiation(f, detector_height, wavelengths):
    f.set('detector_height', detector_height)
    f.set('aerosol_od', 0.28)
    f.set('aerosol_ratio', 1)
    f.set('relative_humidity', 0.5)
    f.set('cloud_liquid', 0)
    f.set('ozone', 0.004)
    f.set('pressure', 1000e2)
    f.set('temperature', 273+15)
    f.set('water_vapor',30)
    f.set('mp_names', station) 
    f.set('mp_concentrations', meta.spm)
    f.set('mcdom_names', station)
    f.set('mcdom_scaling_factors', 1.0)
    f.set('subtract_specular_radiance','false')
    return f.spectrum(wavelengths, band_width, meta.time_point_utc,
                      meta.latitude, meta.longitude)

wl = np.linspace(first_wavelength, last_wavelength, n_wavelengths)

ramses_wl = flick.table(station+"_ocean_irradiance.txt")[:,0]
wl = flick.move_closest_values_in_a_to_those_in_b(wl,ramses_wl)
height = irradiance_height
file_name = station+'_computed_irradiance_W_per_m2_nm.txt'
f = flick.ocean_downward_plane_irradiance()
E = radiation(f, height, wl)
E = f.to_W_per_m2_nm(E)
np.savetxt(file_name, E, fmt=['%6.2f ','%8.3e'])

ramses_wl = flick.table(station+"_ocean_radiance.txt")[:,0]
wl = flick.move_closest_values_in_a_to_those_in_b(wl,ramses_wl)
height = irradiance_height-detector_separation
file_name = station+'_computed_radiance_mW_per_m2_nm_sr.txt'
f = flick.ocean_nadir_radiance()
L = radiation(f, height, wl)
L = f.to_mW_per_m2_nm_sr(L)
np.savetxt(file_name, L, fmt=['%6.2f ','%8.3e'])

#radiation_type = [flick.ocean_downward_plane_irradiance(),
#     flick.ocean_nadir_radiance()]
#detector_height = [irradiance_height,
#                   irradiance_height - detector_separation]
#wavelengths = [[300e-9, 400e-9], [300e-9, 400e-9]]
#
#def change_units(radiation_type, spectrum): 
#    if i==0:
#        return radiation_type.to_W_per_m2_nm(spectrum)
#    if i==1:
#        return radiation_type.to_mW_per_m2_nm_sr(spectrum)
   
#for i in range(2):
#    r = radiation(radiation_type[i], detector_height[i], wavelengths[i])
#    if i==0:
#        r = radiation_type[i].to_W_per_m2_nm(r)
#        file_name = station+'_computed_irradiance.txt'
#    elif i==1:
#        r = radiation_type[i].to_mW_per_m2_nm_sr(r)
#        file_name = station+'_computed_radiance.txt'
#    np.savetxt(file_name, r, fmt=['%6.2f ','%8.3e'])

    

#wavelengths = wavelengths().linspace_with([])

#class wavelengths:
#    first = 320e-9
#    last = 850e-9
#    n_values = 10
#    band_width = 10e-9

#    def linspace_with(self,values):
#        a = np.linspace(self.first,self.last,self.n_values)
#        b = values
#        return move_closest_values_in_a_to_those_in_b(a,b)

    

#radiance_height = irradiance_height - 0.37


#radiation = ['radiance','irradiance']
#wl_grid = np.linspace(from_wl,to_wl, n_wls);
#wl_grid_L = wl_grid 
#wl_grid_E = wl_grid 
#meta = flick.ocean_meta(station+"_meta.txt")

#wl_grid_L = wl_grid_L[np.where((wl_grid_L > from_wl) & (wl_grid_L < to_wl))]
#wl_grid_E = wl_grid_E[np.where((wl_grid_E > from_wl) & (wl_grid_E < to_wl))]

#f = [flick.ocean_downward_plane_irradiance(),
#      flick.ocean_nadir_radiance()]

#        E = f.to_W_per_m2_nm(E)
#    else:
#        L = f.spectrum(wl_grid_L, wl_width, tp, meta.latitude, meta.longitude)
#        L = f.to_mW_per_m2_nm_sr(L)



#if kind=='irradiance':
#        f = flick.ocean_downward_plane_irradiance()
#        wavelengths = 
#        f.set("detector_height", irradiance_height)
#    elif kind=='radiance':
#        f = flick.ocean_nadir_radiance()
#        f.set("detector_height", irradiance_height-detector_separation)
