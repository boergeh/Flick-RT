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

def _run(pre_arg, data, post_arg):
    file_name = 'flick_python_run_tmp.txt'
    np.savetxt(file_name, data)
    matrix = run(pre_arg+' '+file_name+' '+post_arg)
    os.remove(file_name)
    return matrix

def chromaticity(spectrum):
    return _run('filter',spectrum,'chromaticity')

def rgb(spectrum):
    return _run('filter',spectrum,'rgb')
    
def table(file_name):
    return run("text "+file_name+" matrix")

def include_given_values(lower_limit,upper_limit,n,v):
    v = [x for x in v if lower_limit < x < upper_limit]
    v = np.insert(v,0,lower_limit)
    v = np.insert(v,len(v),upper_limit)
    if len(v) >= n:
        step = round(len(v)/n)
        return v[::step]   
    i = 0
    while len(v) < n:
        if i > len(v)-2:
            i = 0
        v = np.insert(v,i+1, v[i]+(v[i+1]-v[i])/2)
        i += 2
    return v

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
        self._b_scaling = 1
        self._salinity = 30
        self._temperature = 290
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

    def set_salinity(self, s):
        self._salinity = s

    def set_temperature(self, t):
        self._temperature = t

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

    def back_scattering_fraction(self,n_terms):
        k = self._expansion_factors(n_terms)
        bb = 0;
        for i in range(len(k)):
            bb += 2*np.pi*k[i,0]*self._legendre_first_half_integral(i)
        l = run("iop scattering_length "+self._expansion_command())
        b = 1/l[1]
        return bb/b

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
        b = 1/l[1]
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
            command = str(self._temperature) +" "+str(self._salinity)
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
        "latitude: " + str(self.latitude) + r" $\degree$"+"\n" +
        "longitude: " + str(self.longitude) + r" $\degree$"+"\n" +
        "UTC time point: " + str(round(self.time_point_utc)) + "\n" +
        "days since year 2000: " + str(self.n_days) + "\n" +
        "depth: " + str(self.depth) + " m\n" +
        "temperature: " + str(round(self.temperature-273.15,1)) + r" $\degree$C"+"\n" +
        "salinity: " + str(self.salinity) + " psu\n" +
        "SPM: " + str(self.spm*1e3) + r" g m$^{-3}$"+"\n" +
        "CHL: " + str(self.chl*1e6) + r" mg m$^{-3}$"+"\n" +
        "POC: " + str(self.poc*1e3) + r" g m$^{-3}$"+"\n" +
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
        num_rows, num_cols = matrix.shape
        if num_rows == 1:
            return matrix[0]
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
    _use_sentinel3_srf = False
    _override_sun_zenith_angle = np.nan
        
    def _generate_config(self, config_type):
        _run_os("mkdir -p "+self._tmpdir)
        run("accurt -g "+config_type+" "+self._config_name)
            
    def set(self, config_parameter, value):
        config(self._config_name, config_parameter, value)

    def use_sentinel3_srf(self,tf):
        self._use_sentinel3_srf = tf
        
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

    def sun_azimuth_angle(self, time_point_utc, latitude, longitude):
        command = "sun_position azimuth_angle "+ \
            self._to_spaced_string(time_point_utc)+" "+ \
            str(latitude)+" "+str(longitude)
        a = run(command)
        return a[0][0]

    def set_override_sun_zenith_angle(self, angle):
        self._override_sun_zenith_angle = angle
        
    def _absolute_spectrum(self, wl_grid, wl_width, time_point_utc,
                           latitude, longitude):
        if np.isnan(self._override_sun_zenith_angle):
            a = self.sun_zenith_angle(time_point_utc, latitude, longitude)
        else:
            a = self._override_sun_zenith_angle

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
        n_points = round(5*(to_wl - from_wl)/wl_width)
        wl = np.linspace(from_wl, to_wl, n_points)
        spectrum = np.empty([len(wl),2])
        spectrum[:,0] = wl
        for i in range(len(wl)):
            if self._use_sentinel3_srf == True:
                spectrum[i,1] = run("filter "+self._tmpdir+"/spectrum sentinel3 "+ \
                                str(wl[i]))
            else:
                spectrum[i,1] = run("filter "+self._tmpdir+"/spectrum triangular "+ \
                                    str(wl[i])+" "+str(wl_width))
                
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
        return self._relative_spectrum(wl_grid,source_zenith_angle)

    def polar_angles(self):
        return np.linspace(0, np.pi, self.n_polar)

    def azimuth_angles(self):
        return np.linspace(-np.pi, np.pi, self.n_azimuth)

    def radiance_surface(self,r):
        theta = self.polar_angles()
        phi = self.azimuth_angles()
        rt = r.transpose() 
        t, p = np.meshgrid(theta, phi)
        x = rt * np.sin(t) * np.cos(p)
        y = rt * np.sin(t) * np.sin(p)
        z = rt * np.cos(t)
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
        
    




