"""
Nadir sub-surface radiance spectrum with a given spectral resolution
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess

class py_accurt:
    def __run_os(command):
        subprocess.run(command, shell=True, check=True, universal_newlines=True)
        subprocess.CalledProcessError

    
    
def run_os(command):
    subprocess.run(command, shell=True, check=True, universal_newlines=True)
    subprocess.CalledProcessError    

def to_string(numbers):
    if isinstance(numbers, np.ndarray):
        return "\""+" ".join(str(x) for x in numbers)+"\""
    return str(numbers)

def set_config(name, value):
    command = "flick text config set "+name+" "+to_string(value)+" > config_tmp"
    run_os(command)
    run_os("mv -f config_tmp config")

def to_streams(n_angles):
    n_streams = np.floor(n_angles**(1/1.6))    
    if (n_streams % 2) != 0:
        n_streams += 1
    return str(n_streams).rstrip('0').rstrip('.');

def configure():
    wavelengths = np.linspace(300e-9,800e-9,30);
    n_angles = 52
    set_config("stream_upper_slab_size", to_streams(n_angles))
    set_config("detector_wavelengths", wavelengths)
    set_config("n_angles", n_angles)
    set_config("n_heights", 3)
    set_config("cdom_440", 0.01)
    set_config("nap_concentration", 0.0)
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




