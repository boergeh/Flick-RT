"""
Nadir sub-surface radiance spectrum with a given spectral resolution
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import flick

def configure(file_name,zenith_angle):
    f = file_name
    if not os.path.exists(f):
        flick.run("accurt -g ocean_radiance "+f)
        
    wavelengths = np.linspace(310e-9,810e-9,30);
    n_angles = 52
    flick.config(f,"stream_upper_slab_size", flick.to_streams(n_angles))
    flick.config(f,"detector_wavelengths", wavelengths)
    flick.config(f,"n_angles", n_angles)
    flick.config(f,"n_heights", 3)
    flick.config(f,"source_zenith_angle", zenith_angle)
    flick.config(f,"cdom_440", 0.15)
    flick.config(f,"nap_concentration", 0.15e-3)
    flick.config(f,"chl_concentration", 0.9e-6)
    flick.config(f,"detector_height", -0.4)
    flick.config(f,"mp_concentrations", 0)
    
def relative_radiance(file_name):
    print("calculating ...")
    return flick.run("accurt "+file_name)

def radiance(relative_radiance, solar_zenith_angle):
    toa_irradiance = flick.run("radiator toa-solar")
    wl = toa_irradiance[:,0];
    r = np.empty([len(wl),2])
    r[:,0] = wl 
    r[:,1] = np.interp(wl,relative_radiance[:,0],relative_radiance[:,1]) * toa_irradiance[:,1]*np.cos(solar_zenith_angle*np.pi/180)
  
    return r

def smooth(radiance):
    np.savetxt("tmp_radiance",radiance)
    wl = np.linspace(300e-9,800e-9,100);
    radiance = np.empty([len(wl),2])
    radiance[:,0] = wl
    for i in range(len(wl)):
        radiance[i,1] = flick.run("filter tmp_radiance gaussian_mean "+str(wl[i])+" 10e-9")
    return radiance

def sun_earth_distance_correction(radiance,time):
    distance_au = flick.run("sun_position distance "+time)
    radiance[:,1] = radiance[:,1] * (1/distance_au)**2
    return radiance

def to_mW_per_m2_per_nm(radiance):
    radiance[:,0] = radiance[:,0]*1e9
    radiance[:,1] = radiance[:,1]*1e-9*1e3
    return radiance

def measured_spectrum():
    r = np.loadtxt("measured_spectrum")
    return r[:,[0,1]]

def plot(radiance):
    m = measured_spectrum()
    x = radiance[:,0]
    y = radiance[:,1]
    fig, ax = plt.subplots()
    line1=ax.plot(m[:,0],m[:,1],linewidth=1,label="measurements")
    line2=ax.plot(x, y,linewidth=1,label="calculations")
    ax.grid()
    ax.set_ylabel("Nadir radiance [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]")
    ax.set_xlabel("Wavelength [nm]")
    ax.legend()
    plt.show()
    fig.savefig("ocean_radiance.pdf", bbox_inches='tight')
    

time_utc = "2023 6 21 13 0 0"
latitude = "60.39"
longitude = "5.32"
config_file_name = "config"
solar_zenith_angle = flick.run("sun_position zenith_angle "+time_utc+" "+latitude+" "+longitude)
configure(config_file_name,solar_zenith_angle[0][0])
r = relative_radiance(config_file_name)
r = radiance(r, solar_zenith_angle)
r = smooth(r)
r = sun_earth_distance_correction(r,time_utc)
r = to_mW_per_m2_per_nm(r)
plot(r)




