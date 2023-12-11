"""
Nadir sub-surface radiance spectrum with a given spectral resolution
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess

def run_os(command):
    subprocess.run(command, shell=True, check=True, universal_newlines=True)
    subprocess.CalledProcessError    

def configure():
    wavelengths = np.linspace(300e-9,500e-9,30);    
    wls_str = "\""+" ".join(str(x) for x in wavelengths)+"\""
    command = "flick text config_template set DETECTOR_WAVELENGTHS "+wls_str+" > updated_config"
    run_os(command)

def calculate():
    print("calculating ...")
    command = "flick accurt updated_config > output_spectrum" 
    run_os(command)    

def plot():
    r = np.loadtxt("output_spectrum")
    x = r[:,0]*1e9
    y = r[:,1]
    fig, ax = plt.subplots()
    ax.plot(x, y,linewidth=1)
    ax.grid()
    ax.set_ylabel("$L_u(ocean) / E_d(toa)$  [sr$^{-1}$]")
    ax.set_xlabel("Wavelength [nm]")
    plt.show();

    
if not os.path.exists("./config_template"):
    run_os("flick accurt -g ocean_radiance config_template")
configure()
calculate()
plot()




