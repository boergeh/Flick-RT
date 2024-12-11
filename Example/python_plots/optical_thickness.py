"""
Curvature sampled absorption optical thickness of the atmosphere
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

path = os.environ['FLICK_PATH']+"/Example/python_plots"
os.chdir(path)

n_wls = 300
T = flick.absorption_optical_thickness("optical_thickness_config",500e-9,900e-9,n_wls).atmosphere()
Tc = flick.save_and_run('filter',T,'curvature_sampled '+str(int(n_wls/10)))

fig, ax = plt.subplots(1,1)
fig.set_size_inches(6,5)
ax.plot(T[:,0]*1e6,T[:,1],'-')
ax.plot(Tc[:,0]*1e6,Tc[:,1],'.-')
ax.grid()
ax.set_xlabel(r"Wavelength [$\mu$m]")
ax.set_ylabel("Atmosphere absorption optical thickness")
if __name__ == "__main__":
    plt.show()








