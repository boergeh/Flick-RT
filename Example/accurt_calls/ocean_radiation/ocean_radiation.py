"""
Modeled and measured sub-surface nadir radiance and irradiance
spectrum with a given spectral band width. See the flick_tmp/config
file that will be generated after the first run for documentation on
all variables that may be set with the 'set' function used in this
example. SI-units and degrees are used unless otherwise specified.

"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

station = "ECOSENS_HF22_D1"
meta = flick.ocean_meta(station+"_meta.txt")
Lm = flick.table(station+"_ocean_radiance.txt")
Em = flick.table(station+"_ocean_irradiance.txt")

wl_grid = np.linspace(300e-9,800e-9,67);
wl_width = 10e-9
f = [flick.ocean_downward_plane_irradiance(),
      flick.ocean_nadir_radiance()]

top_sensor_depth = 0.15
for i in range(2):
    if i==0:
        f = flick.ocean_downward_plane_irradiance()
        f.set("detector_height", -top_sensor_depth)
    else:
        f = flick.ocean_nadir_radiance()
        f.set("detector_height", -top_sensor_depth-0.4)

    f.set_n_angles(100)
    f.set("aerosol_od", 0.28)
    f.set("aerosol_ratio", 1)
    f.set("cloud_liquid", 2e-5)
    f.set("ozone", 0.003)
    f.set("pressure", 1010e2)
    f.set("temperature", 273+15)
    f.set("water_vapor",30)
    f.set("bottom_depth",200)
    f.set("mp_names", station) 
    abs_scale = 1
    f.set("mp_concentrations", meta.spm*abs_scale)
    f.set("mp_scattering_scaling_factors", 1/abs_scale)
    f.set("mcdom_names", station)
    f.set("mcdom_scaling_factors", 1)
    f.set("concentration_relative_depths",[0,0.07,0.071,1])
    f.set("concentration_scaling_factors",[1,1,1,1])
    tp = meta.time_point_utc
    if i==0:
        E = f.spectrum(wl_grid, wl_width, tp, meta.latitude, meta.longitude)
        E = f.to_W_per_m2_nm(E)
    else:
        L = f.spectrum(wl_grid, wl_width, tp, meta.latitude, meta.longitude)
        L = f.to_mW_per_m2_nm_sr(L)

fig, ax = plt.subplots(2,1,sharex=True)
fig.set_size_inches(5,7)
fig.tight_layout(pad=1.0)
line1 = ax[0].plot(L[:,0],L[:,1],label='Modeled')
line2 = ax[0].plot(Lm[:,0],Lm[:,1],label='Measured')
ax[0].legend(loc='upper right')
ax[0].set_ylabel('Nadir radiance [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
ax[0].grid()
line1 = ax[1].plot(E[:,0],E[:,1],label='Modeled')
line2 = ax[1].plot(Em[:,0],Em[:,1]*1e-3,label='Measured')
ax[1].legend(loc='upper right')
ax[1].set_xlabel('Wavelength [nm]')
ax[1].set_ylabel('Downward irradiance [W m$^{-2}$ nm$^{-1}$]')
ax[1].grid()
plt.subplots_adjust(left=0.12, bottom=0.07, right=0.98, top=0.98,
                    wspace=0, hspace=0.04)
plt.show()
fig.savefig("ocean_radiation.pdf", bbox_inches='tight')


