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

def to_string(s):
    if isinstance(s, np.ndarray):
        return "\""+" ".join(str(x) for x in s)+"\""
    return str(s)

def set_config(name,value):
    command = "flick text config set "+name+" "+to_string(value)+" > config_tmp"  
    run_os(command)
    run_os("mv -f config_tmp config")
    
def configure():
    wavelengths = np.linspace(300e-9,800e-9,9);
    set_config("DETECTOR_WAVELENGTHS", wavelengths)
    set_config("cdom_440", 1)
    set_config("chl_concentration", 0)
    set_config("detector_height", -0.5)
    set_config("mp_concentrations", 0)
    
def calculate():
    print("calculating ...")
    command = "flick accurt config > calculated_spectrum" 
    run_os(command)    

def plot():
    r = np.loadtxt("calculated_spectrum")
    x = r[:,0]*1e9
    y = r[:,1]
    fig, ax = plt.subplots()
    ax.plot(x, y,linewidth=1)
    ax.grid()
    ax.set_ylabel("$L_u(ocean) / E_d(toa)$  [sr$^{-1}$]")
    ax.set_xlabel("Wavelength [nm]")
    plt.show();

    
if not os.path.exists("./config"):
    run_os("flick accurt -g ocean_radiance config")
configure()
calculate()
plot()




