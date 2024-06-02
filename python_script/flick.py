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
import scipy.special

def run(arguments):
    flick_output = _run_os("flick "+arguments)
    if flick_output.stdout:
        return _to_matrix(flick_output)

def table(file_name):
    return run("text "+file_name+" matrix")

def move_closest_values_in_a_to_those_in_b(a,b):
    i = 0
    j = 0
    while i < len(a)-1 and j < len(b):
        if (a[i] < b[j] and a[i+1] > b[j]):
            if (abs(a[i]-b[j]) < abs(a[i+1]-b[j])):
                a[i] = b[j]
            else:
                a[i+1] = b[j] 
            j += 1
        i += 1
    return a

class absorption_optical_thickness:
    def __init__(self, ao_config, from_wl, to_wl, n_wls):
        self._ao_config = ao_config
        self._wls = str(from_wl)+" "+str(to_wl)+" "+str(n_wls)

    def atmosphere(self):
        command = "iop absorption_optical_thickness_120000 "+self._wls+ \
            " atmosphere_ocean "+self._ao_config+" 0"
        return run(command)
    
        
class marine_iops:
    def __init__(self, name, spm, from_wl, to_wl, n_wls):
        self._name = name
        self._spm = spm
        self._wls = str(from_wl)+" "+str(to_wl)+" "+str(n_wls)
        
    def a_water(self):
        return self._coefficient("absorption","pure_water")
    
    def a_spm(self):
        return self._coefficient("absorption","marine_particles")
    
    def a_cdom(self):
        return self._coefficient("absorption","marine_cdom")

    def b_water(self):
        return self._coefficient("scattering","pure_water")

    def b_spm(self):
        return self._coefficient("scattering","marine_particles")
    
    def set_b_scaling_factor(self,f):
        self._b_scaling = f

    def volume_scattering_function(self,n_terms,use_cutoff_angle=0):
        p = []
        if n_terms == 0:
            n_points = 1000
            p = run("iop scattering_ab_"+str(n_points)+self._expansion_command())
        else:
            n_points = 200+n_terms**1.6
            p = run("iop scattering_ab_fitted_"+str(n_terms)+"_"+ \
                    str(n_points)+"_"+str(use_cutoff_angle)+ \
                    self._expansion_command())
        p = np.flipud(p)
        p[:,0] = np.arccos(p[:,0])*180/np.pi
        return p

    def back_scattering_coefficient(self, n_terms):
        k = self._expansion_factors(n_terms)
        bb = 0;
        for i in range(len(k)):
            bb += 2*np.pi*k[i,0]*self._legendre_first_half_integral(i)
        return bb

    def _legendre_both_halfs_integral(self,term_number):
        i = term_number
        return self._legendre_first_half_integral(i) + \
            self._legendre_second_half_integral(i)

    def _legendre_second_half_integral(self,term_number):
        # https://math.stackexchange.com/questions/35804/
        # integrating-legendre-polynomials-over-half-range
        l = term_number
        if l == 0:
            return 1
        if l % 2 == 0:
            return 0
        else: # odd terms
            n = (l-1)/2
            return (-1)**n/(2**(2*n+1)*(n+1))*scipy.special.comb(2*n, n)
        
    def _legendre_first_half_integral(self,term_number):
        return (-1)**term_number * self._legendre_second_half_integral(term_number)
        
    def asymmetry_factor(self,n_terms):
        k = self._expansion_factors(n_terms)
        return k[1,0]/k[0,0]/3
     
    def volume_scattering_scaling_factor(self,n_terms):
        k = self._expansion_factors(n_terms)
        l = run("iop scattering_length "+self._expansion_command())
        b = 1/l[0,1]
        f = k[0,0]/b*4*np.pi
        return f

    def _expansion_factors(self,n_terms):
        return run("iop wigner_alpha_beta_"+str(n_terms)+ \
                   self._expansion_command())
        
    def _expansion_command(self):
        return " 515e-9 515e-9 1 marine_particles "+str(self._name)+" "+ \
            str(self._spm)+" "+str(self._b_scaling)
     
    def _coefficient(self,length_type, material):
        command = ""
        if material=="marine_particles":
            command = self._name +" "+str(self._spm)+" "+str(self._b_scaling)
        elif material=="marine_cdom":
            command = self._name +" 1"
        elif material=="pure_water":
            command = "280 30"
        c = run("iop "+length_type+"_length "+self._wls+" "+material+" "+command)
        c[:,1] = 1/c[:,1]
        return c
        
        
class ocean_meta:
    def __init__(self, file_name):
        m = run("text "+file_name+" xy 14")
        self.latitude = m[0,1]
        self.longitude = m[1,1]
        self.time_point_utc = m[2,1]
        self.n_days = m[3,1]
        self.depth = m[4,1]
        self.temperature = m[5,1]+273.15
        self.salinity = m[6,1]
        self.spm = m[7,1]*1e-3
        self.chl = m[8,1]*1e-6
        self.poc = m[9,1]*1e-3

    def to_string(self):
        s = (
        r"$\it{metadata}$" + "\n" +    
        "latitude: " + str(self.latitude) + " $\degree$"+"\n" +
        "longitude: " + str(self.longitude) + " $\degree$"+"\n" +
        "UTC time point: " + str(round(self.time_point_utc)) + "\n" +
        "days since year 2000: " + str(self.n_days) + "\n" +
        "depth: " + str(self.depth) + " m\n" +
        "temperature: " + str(round(self.temperature-273.15,1)) + " $\degree$C\n" +
        "salinity: " + str(self.salinity) + " psu\n" +
        "SPM: " + str(self.spm*1e3) + " g m$^{-3}$\n" +
        "CHL: " + str(self.chl*1e6) + " mg m$^{-3}$\n" +
        "POC: " + str(self.poc*1e3) + " g m$^{-3}$\n" +
            ""
        )
        return s
            

class toa_meta:
    def __init__(self, file_name):
        m = run("text "+file_name+" xy 14")
        self.observation_polar_angle = m[0,1]
        self.observation_azimuth_angle = m[1,1]
        self.solar_zenith_angle = m[2,1]
        self.solar_azimuth_angle = m[3,1]
        self.time_point_utc = m[4,1]
    
def _run_os(command):
    c = subprocess.run("nice "+command, stdout=subprocess.PIPE, shell=True, check=True)
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
    def _generate_config(self, config_type):
        _run_os("mkdir -p "+self._tmpdir)
        run("accurt -g "+config_type+" "+self._config_name)
            
    def set(self, config_parameter, value):
        config(self._config_name, config_parameter, value)
        
    def _relative_spectrum(self, wl_grid, source_zenith_angle):
        self.set("source_zenith_angle", source_zenith_angle)
        self.set("detector_wavelengths", wl_grid)
        return run("accurt "+self._config_name)

    def _to_spaced_string(self,time_point):
        if isinstance(time_point, str):
            return time_point
        s = str(time_point)
        sec = "0"
        if len(s) > 13:
            sec = s[12:14]
        return s[0:4]+" "+s[4:6]+" "+s[6:8]+" "+s[8:10]+" "+s[10:12]+" "+sec     

    def toa_zenith_irradiance(self, time_point_utc):
        distance_au = run("sun_position distance "+ \
                          self._to_spaced_string(time_point_utc))
        r = run("radiator toa-solar")
        r[:,1] = r[:,1] * (1/distance_au)**2
        return r

    def sun_zenith_angle(self,time_point_utc, latitude, longitude):
        command = "sun_position zenith_angle "+ \
        self._to_spaced_string(time_point_utc)+" "+ \
        str(latitude)+" "+str(longitude)
        a = run(command)
        return a[0][0]
    
    def _absolute_spectrum(self, wl_grid, wl_width, time_point_utc,
                           latitude, longitude):
        a = self.sun_zenith_angle(time_point_utc, latitude, longitude)
        L_r = self._relative_spectrum(wl_grid,a)
        F_0 = self.toa_zenith_irradiance(time_point_utc);
        wl = F_0[:,0];
        L_r = self._interpolate(L_r, wl)
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
        n_points = round(3*(to_wl - from_wl)/wl_width)
        wl = np.linspace(from_wl, to_wl, n_points)
        spectrum = np.empty([len(wl),2])
        spectrum[:,0] = wl
        for i in range(len(wl)):
            spectrum[i,1] = run("filter "+self._tmpdir+"/spectrum triangular "+ \
                                str(wl[i])+" "+str(wl_width))
            #spectrum[i,1] = run("filter "+self._tmpdir+"/spectrum gaussian_mean "+ \
             #         str(wl[i])+" "+str(wl_width))
        return spectrum

    def set_n_angles(self,n_angles):
         self.set("n_angles", round(n_angles)) 
         self.set("stream_upper_slab_size", to_streams(round(n_angles))) 


class radiance_distribution(basic_radiation):
    def __init__(self, n_polar, n_azimuth):
        self._generate_config("toa_reflectance")
        self.n_polar = n_polar
        self.n_azimuth = n_azimuth
        self.set("detector_radiance_distribution_override",
                 [n_polar, n_azimuth])
        self.set_n_angles(16**1.6)
        
    def values(self, wl_grid, source_zenith_angle):
        return self._relative_spectrum(wl_grid,source_zenith_angle).transpose()

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
    def spectrum(self, wl_grid, source_zenith_angle):
        return self._relative_spectrum(wl_grid, source_zenith_angle)


class absolute_radiation(basic_radiation):
    def spectrum(self,wl_grid, wl_width, time_point_utc,
                                       latitude, longitude):
        return self._absolute_spectrum(wl_grid, wl_width, time_point_utc,
                                       latitude, longitude)

    
class radiance(absolute_radiation):
    def __init__(self, polar_viewing_angle, azimuth_viewing_angle):
        self._generate_config("toa_reflectance")
        self.set_n_angles(16**1.6)
        self.set("detector_orientation_override",[polar_viewing_angle,
                                                  azimuth_viewing_angle])

    def to_mW_per_m2_nm_sr(self,spectrum):
        spectrum[:,0] = spectrum[:,0]*1e9
        spectrum[:,1] = spectrum[:,1]*1e-6
        return spectrum

class plane_irradiance(absolute_radiation):
    def __init__(self):
        self._generate_config("ocean_radiance")
        self.set("detector_type","plane_irradiance")
        self.set("detector_orientation","up")
        self.set_n_angles(16**1.6)
        
    def to_W_per_m2_nm(self,spectrum):
        spectrum[:,0] = spectrum[:,0]*1e9
        spectrum[:,1] = spectrum[:,1]*1e-9
        return spectrum
    
    
class toa_radiance(radiance):
    pass


class ocean_nadir_radiance(radiance):
    def __init__(self):
        self._generate_config("ocean_radiance")

        
class ocean_downward_plane_irradiance(plane_irradiance):
    pass


class ocean_upward_plane_irradiance(plane_irradiance):
    def __init__(self):
        self.set("detector_orientation","down")

        
class remote_sensing_reflectance(relative_radiation):
    def __init__(self):
        self._generate_config("rs_reflectance")
        self.set_n_angles(16**1.6)

        
class surface_irradiance(absolute_radiation):
    def __init__(self):
        self._generate_config("boa_transmittance")
        self.set_n_angles(8**1.6)

    def to_W_per_m2_nm(self,spectrum):
        spectrum[:,0] = spectrum[:,0]*1e9
        spectrum[:,1] = spectrum[:,1]*1e-9
        return spectrum

    
class snow_transmittance(relative_radiation):
    def __init__(self):
        self._generate_config("toa_reflectance")
        self.set("detector_height",0)
        self.set("detector_orientation","up")
        self.set("reference_detector_height", 1.01)
        self.set("detector_type","plane_irradiance")
        
    




