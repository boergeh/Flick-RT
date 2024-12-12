"""
Wavelength sampling based on second derivative of atmospheric
transmittance. Ensures grid points in absorption bands.
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

path = os.environ['FLICK_PATH']+"/Example/python_plots"
os.chdir(path)

full_spectrum = False

def transmittance(optical_thickness):
    T = optical_thickness
    T[:,1] = np.exp(-optical_thickness[:,1])
    return T

if full_spectrum:
    wl_low = 290e-9
    wl_high = 1040e-9
    n_wls = 230
else:
    wl_low = 600e-9
    wl_high = 800e-9
    n_wls = 20
    
ot = flick.atmosphere_optical_thickness("optical_thickness_config",
                                          wl_low,wl_high,n_wls*10).attenuation();
T = transmittance(ot)
Tc = flick.save_and_run('filter',T,'curvature_sampled '+str(n_wls))

fig, ax = plt.subplots(1,1)
fig.set_size_inches(6,5)
ax.plot(T[:,0]*1e6,T[:,1],'-',label='High resolution transmittance')
ax.plot(Tc[:,0]*1e6,Tc[:,1],'.-',label='Curvature sampled')
ax.grid()
ax.legend() 
ax.set_xlabel(r"Wavelength [$\mu$m]")
ax.set_ylabel("Atmospheric transmittance")
if __name__ == "__main__":
    plt.show()








