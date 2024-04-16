"""
Top-of-atmosphere radiance spectrum with a given spectral band
width. See the flick_tmp/config file, which will be generated after
the first run, for documentation on all variables that may be set with
the set function used in this example. SI-units and degrees are used
unless otherwise specified. Measurements by OLCI on Sentinel 3,
https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/olci/
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

station = "ECOSENS_HF22_D1"
from_wl = 380e-9
#to_wl = 710e-9
to_wl = 1030e-9
meta = flick.ocean_meta(station+"_meta.txt")
toa_meta = flick.toa_meta(station+"_toa_meta.txt")
Lm = flick.table(station+"_toa_radiance.txt")
wl_grid = np.linspace(from_wl,to_wl,65);
wl_width = 5e-9
polar_viewing_angle = 180-toa_meta.observation_polar_angle # check this
azimuth_viewing_angle = toa_meta.observation_azimuth_angle # check this
f = flick.toa_radiance(polar_viewing_angle, azimuth_viewing_angle)

f.set_n_angles(200)
f.set("n_heights",8)
f.set("gas_spectral_region","uv_vis_toa")
#f.set("gases","o2")
f.set("aerosol_od", 0.28)
f.set("aerosol_ratio", 1)
f.set("cloud_liquid", 0e-6)
f.set("ozone", 0.004)
f.set("pressure", 1010e2)
f.set("temperature", 273+15)
f.set("water_vapor",30)
f.set("mp_names", station)
abs_scale = 1
f.set("mp_concentrations", meta.spm*abs_scale)
f.set("mp_scattering_scaling_factors", 1.0/abs_scale)
f.set("mcdom_names", station)
f.set("mcdom_scaling_factors", 1)
f.set("concentration_relative_depths",[0,0.7,0.71,1])
f.set("concentration_scaling_factors",[1,1,1,1])

tp = toa_meta.time_point_utc
L = f.spectrum(wl_grid, wl_width, tp, meta.latitude, meta.longitude)
L = f.to_mW_per_m2_nm_sr(L)      

fig, ax = plt.subplots()
fig.set_size_inches(5,4)
ax.plot(L[:,0],L[:,1],label='modeled')
ax.plot(Lm[:,0],Lm[:,1],'o',fillstyle='none',label='measured')
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel('TOA radiance [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
ax.grid()
ax.legend()
plt.show()
fig.savefig("toa_radiance.png", dpi=300)








