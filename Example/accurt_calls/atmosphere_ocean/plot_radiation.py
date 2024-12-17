"""
Plot computed and measured radiation.
Run compute_radiation.py first.

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

path = os.environ['FLICK_PATH']+"/Example/accurt_calls/atmosphere_ocean"
os.chdir(path)

station = "ECOSENS_HF22_D1"
show_computed_data = True
show_satellite_data = False
show_ramses_data = True

fig, ax = plt.subplots(2,1,sharex=True)
fig.set_size_inches(4.4,7.5)
if show_computed_data:
    Lc = flick.table('output/'+station+"_computed_radiance_mW_per_m2_nm_sr.txt")    
    line1 = ax[0].plot(Lc[:,0],Lc[:,1],label='Computed')
if show_ramses_data:
    Lm = flick.table('input/'+station+"_ocean_radiance.txt")
    line2 = ax[0].plot(Lm[:,0],Lm[:,1],label='Ramses')
if show_satellite_data:
    Ls = flick.table('input/'+station+"_toa_radiance.txt")
    line3 = ax[0].plot(Ls[:,0],Ls[:,1],'o',fillstyle='none',label='Satellite')
ax[0].set_ylabel('Nadir radiance [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
ax[0].grid()
ax[0].legend()
ax[0].set_ylim([0.5,200])

def PAR(E):
    x = np.linspace(400,700,300)
    y = np.interp(x, E[:,0], E[:,1])
    return np.trapezoid(y,x)
    
if show_computed_data:
    Ec = flick.table('output/'+station+"_computed_irradiance_W_per_m2_nm.txt")
    line1 = ax[1].plot(Ec[:,0],Ec[:,1],label='Computed (PAR: '+str(round(PAR(Ec)))+' Wm$^{-2}$)')
if show_ramses_data:
    Em = flick.table('input/'+station+"_ocean_irradiance.txt")
    Em[:,1] *= 1e-3
    line2 = ax[1].plot(Em[:,0],Em[:,1],label='Ramses (PAR: '+str(round(PAR(Em)))+' Wm$^{-2}$)')
ax[1].set_xlabel('Wavelength [nm]')
ax[1].set_ylabel('Downward irradiance [W m$^{-2}$ nm$^{-1}$]')
ax[1].legend()
ax[1].set_ylim([5e-2,3])
for i in range(2):
    ax[i].set_xscale('log')
    ax[i].set_yscale('log')
    ax[i].xaxis.set_minor_formatter(mticker.ScalarFormatter())
    ax[i].yaxis.set_minor_formatter(mticker.ScalarFormatter())
    ax[i].tick_params(axis='both', which='major', labelsize=9)
    ax[i].tick_params(axis='both', which='minor', labelsize=9)
    ax[i].yaxis.set_minor_formatter(mticker.FormatStrFormatter("%3.3g"))
    ax[i].minorticks_on()
    ax[i].grid(True, which='both')
    for label in ax[i].yaxis.get_ticklabels('minor')[1::2]:
        label.set_visible(False)
    
plt.subplots_adjust(left=0.15, bottom=0.07, right=0.97, top=0.98,
                    wspace=0, hspace=0.04)

fig.savefig('output/'+station+'_plotted_radiation.pdf')

if __name__ == "__main__":
    plt.show()




