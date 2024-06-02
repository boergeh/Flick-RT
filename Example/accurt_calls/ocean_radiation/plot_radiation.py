"""
Plot measured and computed downward irradiance and nadir radiance.

"""
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

compute_spectra = False

if compute_spectra:
    exec(open("compute_radiation.py").read())
else:
    station = "ECOSENS_HF22_D1"

Em = flick.table(station+"_ocean_irradiance.txt")
Lm = flick.table(station+"_ocean_radiance.txt")
Ec = flick.table(station+"_computed_irradiance_W_per_m2_nm.txt")
Lc = flick.table(station+"_computed_radiance_mW_per_m2_nm_sr.txt")

fig, ax = plt.subplots(2,1,sharex=True)
fig.set_size_inches(5,7)
line1 = ax[0].plot(Lc[:,0],Lc[:,1],label='Computed')
line2 = ax[0].plot(Lm[:,0],Lm[:,1],label='Measured')
ax[0].set_ylabel('Nadir radiance [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
ax[0].grid()
line1 = ax[1].plot(Ec[:,0],Ec[:,1],label='Computed')
line2 = ax[1].plot(Em[:,0],Em[:,1]*1e-3,label='Measured')
ax[1].set_xlabel('Wavelength [nm]')
ax[1].set_ylabel('Downward irradiance [W m$^{-2}$ nm$^{-1}$]')
for i in range(2):
    ax[i].set_xscale('log')
    ax[i].set_yscale('log')
    ax[i].xaxis.set_minor_formatter(mticker.ScalarFormatter())
    ax[i].yaxis.set_minor_formatter(mticker.ScalarFormatter())
    ax[i].tick_params(axis='both', which='major', labelsize=10)
    ax[i].tick_params(axis='both', which='minor', labelsize=9)
    ax[i].yaxis.set_minor_formatter(mticker.FormatStrFormatter("%3.3g"))
    ax[i].minorticks_on()
    ax[i].grid(True, which='both')
    for label in ax[i].yaxis.get_ticklabels('minor')[1::2]:
        label.set_visible(False)
    
plt.subplots_adjust(left=0.15, bottom=0.07, right=0.98, top=0.98,
                    wspace=0, hspace=0.04)
plt.show()
fig.savefig(station+'_plotted_radiation.pdf')



