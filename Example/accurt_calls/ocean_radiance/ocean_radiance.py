"""
Nadir sub-surface radiance spectrum with a given spectral resolution
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import flick

def configure(file_name, zenith_angle):
    f = file_name
    flick.run("accurt -g ocean_radiance "+f)
        
    wavelengths = np.linspace(300e-9,950e-9,20);
    n_angles = 30
    flick.config(f,"stream_upper_slab_size", flick.to_streams(n_angles))    
    flick.config(f,"detector_wavelengths", wavelengths)
    flick.config(f,"n_angles", n_angles)
    flick.config(f,"n_heights", 3)
    flick.config(f,"pressure", 1000e2)
    flick.config(f,"source_zenith_angle", zenith_angle)
    flick.config(f,"aerosol_od", 0.07)
    flick.config(f,"cloud_liquid", 0)
    flick.config(f,"cdom_440", 0.09)
    flick.config(f,"nap_concentration", 2.5e-3)
    flick.config(f,"chl_concentration", 2e-6)
    
def relative_radiance(file_name):
    print("calculating ...")
    return flick.run("accurt "+file_name)

def radiance(relative_radiance, solar_zenith_angle):
    toa_irradiance = flick.run("radiator toa-solar")
    wl = toa_irradiance[:,0];
    r = np.empty([len(wl),2])
    r[:,0] = wl 
    r[:,1] = np.interp(wl,relative_radiance[:,0],relative_radiance[:,1]) \
        * toa_irradiance[:,1]*np.cos(solar_zenith_angle*np.pi/180)
    return r

def smooth(radiance):
    np.savetxt(".tmp_radiance",radiance)
    wl = np.linspace(300e-9,900e-9,100);
    radiance = np.empty([len(wl),2])
    radiance[:,0] = wl
    for i in range(len(wl)):
        radiance[i,1] = flick.run("filter .tmp_radiance gaussian_mean "+
                                  str(wl[i])+" 10e-9")
    return radiance

def sun_earth_distance_correction(radiance,time):
    distance_au = flick.run("sun_position distance "+time)
    radiance[:,1] = radiance[:,1] * (1/distance_au)**2
    return radiance

def to_mW_per_m2_per_nm(radiance):
    radiance[:,0] = radiance[:,0]*1e9
    radiance[:,1] = radiance[:,1]*1e-9*1e3
    return radiance

def measured_radiance():
    r = np.loadtxt("measured_radiance")
    return r[:,[0,1]]

def measured_toa_reflectance():
    r = np.loadtxt("measured_toa_reflectance")
    return r

def ocean_radiance():
    flick.config(config_file_name, "subtract_specular_radiance", "false")
    flick.config(config_file_name, "reference_detector_height", 120e3)
    flick.config(config_file_name, "detector_height", -0.3)
    r = relative_radiance(config_file_name)
    r = radiance(r, solar_zenith_angle)
    r = smooth(r)
    r = sun_earth_distance_correction(r, time_utc)
    r = to_mW_per_m2_per_nm(r)
    return r

def toa_reflectance():
    flick.config(config_file_name, "reference_detector_height", 120e3)
    flick.config(config_file_name, "detector_height", 120e3)
    flick.config(config_file_name, "subtract_specular_radiance", "false")
    r = relative_radiance(config_file_name)
    r[:,0] = r[:,0]*1e9
    return r

def remote_sensing_reflectance():
    flick.config(config_file_name, "reference_detector_height", 0.01)
    flick.config(config_file_name, "detector_height", 0.01)
    flick.config(config_file_name, "subtract_specular_radiance", "true")
    r = relative_radiance(config_file_name)
    r[:,0] = r[:,0]*1e9
    return r

def plot_ocean_radiance(ax, radiance):
    m = measured_radiance()
    x = radiance[:,0]
    y = radiance[:,1]
    line1=ax.plot(m[:,0],m[:,1],linewidth=1,label="measurements")
    line2=ax.plot(x, y,linewidth=1,label="calculations")
    ax.grid()
    ax.set_ylabel("Nadir radiance [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]")
    ax.legend()

def plot_reflectances(ax, toa, rrs, toam):
    line1=ax.plot(toa[:,0],toa[:,1],linewidth=1,label= \
                  "top of atmosphere, nadir")
    line2=ax.plot(rrs[:,0],rrs[:,1],linewidth=1,label="ocean remote sensing")
    line3=ax.plot(toam[:,0],toam[:,1],marker='o',markersize=2,linewidth=1,
                  label="Sentinel 3, TOA, slant viewing")
    ax.grid()
    ax.set_ylabel("Nadir reflectance [sr$^{-1}$]")
    ax.set_xlabel("Wavelength [nm]")
    ax.legend()
    
    
time_utc = "2022 5 20 11 24 0"
latitude = "59.878675"
longitude = "5.655558"
config_file_name = ".tmp_flick_config"
solar_zenith_angle = flick.run("sun_position zenith_angle "+time_utc+" "+
                               latitude+" "+longitude)
configure(config_file_name, solar_zenith_angle[0][0])

a = ocean_radiance()
b = toa_reflectance()
c = remote_sensing_reflectance()
d = measured_toa_reflectance()

fig, ax = plt.subplots(2,1,sharex=True)
fig.set_size_inches(5,6)
plot_ocean_radiance(ax[0],a)
plot_reflectances(ax[1],b,c,d)
plt.subplots_adjust(left=0.15,
                    bottom=0.08,
                    right=0.98,
                    top=0.98,
                    wspace=0.0,
                    hspace=0.04)
plt.show()
fig.savefig("ocean_radiance.pdf", bbox_inches='tight')








