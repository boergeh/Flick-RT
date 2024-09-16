import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.ticker as mticker
import netCDF4

nc = netCDF4.Dataset('./S3A_OL_SRF_20160713_mean_rsr.nc4')
#nc = netCDF4.Dataset('./S3B_OL_SRF_0_20180109_mean_rsr.nc4')

n_bands = 21;
center_wl = np.zeros((n_bands,2))
for i in range(n_bands):
    wl = nc['mean_spectral_response_function_wavelength'][i,:]
    srf = nc['mean_spectral_response_function'][i,:]
    center_wl[i,:] = [(wl[-1]+wl[0])/2, i]
    data = np.column_stack([wl, srf])
    np.savetxt('srf/band_'+str(i)+'.txt',data,fmt='%.6g  ')
np.savetxt('srf/center_wavelength.txt',center_wl,fmt='%.0f  %.6g')

fig, ax = plt.subplots()
fig.set_size_inches(6,1.0)
from_wl = 380
to_wl = 1050
for i in range(21):
    wl = nc['mean_spectral_response_function_wavelength'][i,:]
    srf = nc['mean_spectral_response_function'][i,:]
    color = cm.rainbow(i*30)
    if wl[-1] > 720:
        color = cm.rainbow(i*30+30)
    plt.semilogx(wl,srf,color=color)
ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
ax.tick_params(axis='both', which='major', labelsize=9)
ax.tick_params(axis='both', which='minor', labelsize=9)
ax.set_xlim([from_wl,to_wl])
ax.minorticks_on()
ax.grid(True, which='both')
plt.show()

