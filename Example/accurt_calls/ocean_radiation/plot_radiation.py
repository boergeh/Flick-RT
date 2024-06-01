"""
Plot measured and computed downward irradiance and nadir radiance
"""
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

station = "ECOSENS_HF22_D1"
Em = flick.table(station+"_ocean_irradiance.txt")
Lm = flick.table(station+"_ocean_radiance.txt")
Ec = flick.table(station+"_computed_irradiance_W_per_m2_nm.txt")
Lc = flick.table(station+"_computed_radiance_mW_per_m2_nm_sr.txt")

fig, ax = plt.subplots(2,1,sharex=True)
fig.set_size_inches(5,7)
fig.tight_layout(pad=1.0)
line1 = ax[0].semilogy(Lc[:,0],Lc[:,1],label='Computed')
line2 = ax[0].semilogy(Lm[:,0],Lm[:,1],label='Measured')
ax[0].legend(loc='upper right')
ax[0].set_ylabel('Nadir radiance [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
ax[0].grid()
line1 = ax[1].semilogy(Ec[:,0],Ec[:,1],label='Computed')
line2 = ax[1].semilogy(Em[:,0],Em[:,1]*1e-3,label='Measured')
ax[1].legend(loc='upper right')
ax[1].set_xlabel('Wavelength [nm]')
ax[1].set_ylabel('Downward irradiance [W m$^{-2}$ nm$^{-1}$]')
ax[1].grid()
plt.subplots_adjust(left=0.15, bottom=0.07, right=0.98, top=0.98,
                    wspace=0, hspace=0.04)
plt.show()
fig.savefig(station+'_plotted_radiation.pdf')



