"""
Plot measured and computed downward irradiance and nadir radiance
"""
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

station = "ECOSENS_HF22_D1"
E = flick.table(station+"_computed_irradiance_W_per_m2_nm.txt")
L = flick.table(station+"_computed_radiance_mW_per_m2_nm_sr.txt")
L[:,1] = interp1d(L[:,0], L[:,1])(E[:,0])
to_watts = 1e-3;
fig, ax = plt.subplots()
fig.set_size_inches(5,4)
ax.semilogy(E[:,0],to_watts*L[:,1]/E[:,1])
ax.set_ylabel('Radiance irradiance ratio [sr$^{-1}$]')
ax.grid()
ax.set_xlabel('Wavelength [nm]')
plt.show()
fig.savefig(station+'_plotted_radiation_ratio.pdf')



