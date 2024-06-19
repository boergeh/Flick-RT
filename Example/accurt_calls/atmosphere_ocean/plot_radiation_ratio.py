"""
Plot ratio of nadir radiance to downward irradiance.
Run compute_radiation.py first.
"""
#from scipy.interpolate import interp1d
import numpy as np
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

station = "ECOSENS_HF22_D1"

E = flick.table('output/'+station+"_computed_irradiance_W_per_m2_nm.txt")
L = flick.table('output/'+station+"_computed_radiance_mW_per_m2_nm_sr.txt")
wls = L[:,0]
E_new = np.interp(wls,E[:,0],E[:,1])
to_watts = 1e-3;
fig, ax = plt.subplots()
fig.set_size_inches(4,5)
ax.plot(wls,to_watts*L[:,1]/E_new*1000)
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel(r'Remote sensing reflectance $\times$ 1000 [sr$^{-1}$]')

ax.set_xscale('log')
ax.set_yscale('log')
ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
ax.tick_params(axis='both', which='major', labelsize=9)
ax.tick_params(axis='both', which='minor', labelsize=9)
ax.yaxis.set_minor_formatter(mticker.FormatStrFormatter("%3.3g"))
ax.minorticks_on()
ax.grid(True, which='both')
for label in ax.yaxis.get_ticklabels('minor')[1::2]:
        label.set_visible(False)

plt.subplots_adjust(left=0.15, bottom=0.09, right=0.98, top=0.98)
plt.show()
fig.savefig('output/'+station+'_plotted_radiation_ratio.pdf')



