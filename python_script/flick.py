"""
flick wrapper

example:
  import flick
  m = flick.run("radiator planck 6000 9")
  print(m)
"""
import numpy as np
import subprocess
import os

def run(arguments):
    flick_output = _run_os("flick "+arguments)
    if flick_output.stdout:
        return _to_matrix(flick_output)

def _run_os(command):
    c = subprocess.run(command, stdout=subprocess.PIPE, shell=True, check=True)
    return c

def _to_matrix(flick_output):
    try:
        l = flick_output.stdout.decode('utf-8');
        l = l.splitlines()
        l = list(filter(None, l))
        if not l:
            return np.ndarray(shape=(0,0),dtype=float) 
        numbers_first_line = [float(x) for x in l[0].split()]
        rows = len(l);
        cols = len(numbers_first_line)
        matrix = np.ndarray(shape=(rows,cols),dtype=float)
        for i in range(rows):
            floats = [float(x) for x in l[i].split()]
            matrix[i] = floats
        return matrix
    except Exception as e:
        raise Exception(flick_output.stdout.decode('utf-8')) 

def _to_string(numbers):
    if isinstance(numbers, np.ndarray) or isinstance(numbers, list):
        return "\""+" ".join(str(x) for x in numbers)+"\""
    return str(numbers)

def config(file_name, parameter_name, value):
    command = "flick text "+file_name+" set "+parameter_name+" "+ \
        _to_string(value)+" > "+file_name+"_tmp"
    _run_os(command)
    _run_os("mv -f "+file_name+"_tmp "+ file_name)

def to_streams(n_angles):
    n_streams = np.floor(n_angles**(1/1.6))    
    if (n_streams % 2) != 0:
        n_streams += 1
    return str(n_streams).rstrip('0').rstrip('.');

class basic_radiation:
    _tmpdir = "flick_tmp"
    _config_name = _tmpdir+"/config"
    def _ensure_config_exists(self):
        if not os.path.isfile(self._config_name):
            _run_os("mkdir -p "+self._tmpdir)
            run("accurt -g toa_reflectance "+self._config_name)
            
    def set(self, config_parameter, value):
        self._ensure_config_exists()
        config(self._config_name, config_parameter, value)
        
    def _relative_spectrum(self):
        self._ensure_config_exists()
        return run("accurt "+self._config_name)

    def toa_zenith_irradiance(self, time_point_utc):
        distance_au = run("sun_position distance "+time_point_utc)
        r = run("radiator toa-solar")
        r[:,1] = r[:,1] * (1/distance_au)**2
        return r

    def sun_zenith_angle(self,time_point_utc, latitude, longitude):
        a = run("sun_position zenith_angle "+time_point_utc+" "+ \
                               latitude+" "+longitude)
        return a[0][0]
        
    def _absolute_spectrum(self, wl_grid, wl_width, time_point_utc,
                           latitude, longitude):
        self.set("wavelengths", wl_grid)
        L_r = self._relative_spectrum()
        F_0 = self.toa_zenith_irradiance(time_point_utc);
        wl = F_0[:,0];
        L_r = self._interpolate(L_r, wl)
        a = self.sun_zenith_angle(time_point_utc, latitude, longitude)
        s = L_r[:,1] * F_0[:,1] * np.cos(a*np.pi/180);
        points = np.vstack((wl,s)).T
        return self.smooth(points, wl_grid[0], wl_grid[-1], wl_width)

    def _interpolate(self, points, x):       
         x0 = points[:,0] 
         y0 = points[:,1]
         new_points = np.empty([len(x),2])
         new_points[:,0] = x
         new_points[:,1] = np.interp(x,x0,y0)
         return new_points

    def smooth(self, spectrum, from_wl, to_wl, wl_width):
        np.savetxt(self._tmpdir+"/spectrum", spectrum)
        n_points = round(2*(to_wl - from_wl)/wl_width)
        wl = np.linspace(from_wl, to_wl, n_points)
        spectrum = np.empty([len(wl),2])
        spectrum[:,0] = wl
        for i in range(len(wl)):
            spectrum[i,1] = run("filter "+self._tmpdir+"/spectrum gaussian_mean "+ \
                      str(wl[i])+" "+str(wl_width))
        return spectrum

    def set_n_angles(self,n_angles):
         self.set("n_angles", round(n_angles)) 
         self.set("stream_upper_slab_size", to_streams(round(n_angles))) 


class radiance_distribution(basic_radiation):
    def __init__(self, wl, n_polar, n_azimuth):
        self.wl = wl
        self.n_polar = n_polar
        self.n_azimuth = n_azimuth
        self.set("detector_wavelengths",self.wl)
        self.set("detector_radiance_distribution_override",
                 [self.n_polar, self.n_azimuth])
        self.set_n_angles(16**1.6)
        
    def values(self):
        return self._relative_spectrum().transpose()

    def polar_angles(self):
        return np.linspace(0, np.pi, self.n_polar)

    def azimuth_angles(self):
        return np.linspace(-np.pi, np.pi, self.n_azimuth)

    def radiance_surface(self,r):
        theta = self.polar_angles()
        phi = self.azimuth_angles()
        t, p = np.meshgrid(theta, phi)
        x = r * np.sin(t) * np.cos(p)
        y = r * np.sin(t) * np.sin(p)
        z = r * np.cos(t)
        return x,y,z
    
    
class relative_radiation(basic_radiation):
    def spectrum(self):
        return self._relative_spectrum()


class absolute_radiation(basic_radiation):
    def spectrum(self,wl_grid, wl_width, time_point_utc,
                                       latitude, longitude):
        return self._absolute_spectrum(wl_grid, wl_width, time_point_utc,
                                       latitude, longitude)

    
class radiance(absolute_radiation):
    def __init__(self, polar_viewing_angle, azimuth_viewing_angle):
        self.set_n_angles(16**1.6)
        self.set("detector_orientation_override",[polar_viewing_angle,
                                                  azimuth_viewing_angle])

    def to_mW_per_m2_nm_sr(self,spectrum):
        spectrum[:,0] = spectrum[:,0]*1e9
        spectrum[:,1] = spectrum[:,1]*1e-6
        return spectrum

    
class toa_radiance(radiance):
    pass


class ocean_nadir_radiance(radiance):
    def __init__(self):
        self.set("detector_height", -0.5)
        self.set("detector_orientation_override",[0,0])

        
class remote_sensing_reflectance(relative_radiation):
    def __init__(self):
        self.set("detector_height", 0.1)
        self.set("subtract_specular_radiance", "true")
        self.set_n_angles(16**1.6)

        
class surface_irradiance(absolute_radiation):
    def __init__(self):
        self.set("detector_height", 0.1)
        self.set("detector_orientaion", "up")
        self.set("detector_type", "irradiance")
        self.set_n_angles(8**1.6)

    def to_W_per_m2_nm(self,spectrum):
        spectrum[:,0] = spectrum[:,0]*1e9
        spectrum[:,1] = spectrum[:,1]*1e-3
        return spectrum
    




