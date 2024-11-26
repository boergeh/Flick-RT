import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

def plot(spectrum, label):
    to_nm = 1e9;
    to_per_nm = 1e-9
    n_photons = flick.save_and_run('filter',spectrum,'n_photons 400e-9 700e-9')
    n = int(n_photons/6.023e23*1e6)
    plt.plot(spectrum[:,0]*to_nm, spectrum[:,1]*to_per_nm,label=label
             +', PAR: '+str(n)+r' $\mu$mol$\,$m$^{-2}\,$s$^{-1}$')
    
toa = flick.run('radiator toa_solar')
surf = flick.run('radiator surface_reference')

plot(toa,r'TOA (high resolution), $\theta=0\degree$')
plot(surf,r'Surface (low resolution), $\theta=60\degree$')
plt.legend()
plt.xlim([290, 1050])
plt.xlabel('Wavelength [nm]')
plt.ylabel(r'Irradiance [W$\,$m$^{-2}\,$nm$^{-1}$]')
plt.grid()

if __name__ == "__main__":
    plt.show()
