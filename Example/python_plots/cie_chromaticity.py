"""
Compute the CIE 1931 two-degree xy chromaticity values and
corresponding RGB values from a tabulated input spectrum.

"""

import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
import sys
sys.path.append(os.environ['FLICK_PATH']+'/python_script')
import flick

path = os.environ['FLICK_PATH']+"/Example/python_plots"
os.chdir(path)

def nm_to_m(E):
    E[:,0] *= 1e-9
    E[:,1] *= 1e9
    return E

def m_to_nm(E):
    E[:,0] *= 1e9
    E[:,1] *= 1e-9
    return E

# Read two-column ascii file
E = flick.read('flick_standard_spectrum.txt')
E = nm_to_m(E)
xyz = flick.chromaticity(E)
xy = xyz[0:2]
rgb = flick.rgb(E)
E = m_to_nm(E)

# Show spectrum with xy coordinate and RGB color
fig, ax = plt.subplots()
ax.plot(E[:,0],E[:,1])
ax.set_title('Flick standard irradiance spectrum, CIE 1931 xy = ' + \
          str(np.round(xy,3))+ \
          ',  RGB = '+str(np.round(rgb,3)),fontsize=9)
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel(r'Irradiance [W$\,$m$^{-2}\,$nm$^{-1}$]')
ax.grid()
ax.set_facecolor(rgb)
if __name__ == "__main__":
    plt.show()
