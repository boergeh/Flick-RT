"""

"""

import numpy as np

class iops:
    def __init__(self, absorption_coefficient, scattering_coefficient, asymmetry_factor):
        self._k_abs = absorption_coefficient
        self._k_sca = scattering_coefficient
        self._g = asymmetry_factor

    def asymmetry_factor(self):
        return self._g
    
    def single_scattering_albedo(self):
        return self._k_sca / (self._k_abs + self._k_sca)

    def extinction_coefficient(self):
        return self._k_sca + self._k_abs
    
    def optical_depth(self, geometrical_depth):
        return self.extinction_coefficient()*geometrical_depth

    def scaled_optical_depth(self, geometrical_depth):
        return (1-self._g)*self.optical_depth(geometrical_depth)

    
class iops_generator:  
    def __init__(self, volume_fraction, grain_radius, refractive_index, wavelength):
        self._volume_fraction = volume_fraction
        self._grain_radius = grain_radius
        self._refractive_index = refractive_index
        self._wavelength = wavelength

    def material_absorption_coefficient(self):
        return 4*np.pi*np.imag(self._refractive_index)/self._wavelength
    
    #k_ext = 3*volume_fraction
    #k_abs = 0
    #k_sca = k_ext - k_abs
    #g = 0.75
    # eq 2.39
    #n = 1.31
    #rho = 0.0123+0.1622*(n-1)

    
class low_absorption_spheres(iops_generator):
    """
    E.g. Bohren, C.F. and Clothiaux, E.E., 2006. Fundamentals of
    atmospheric radiation: an introduction with 400 problems. John
    Wiley & Sons.

    """
    def get_iops(self) :
        k_abs = self._volume_fraction * self.material_absorption_coefficient()
        k_sca = 3/2*self._volume_fraction/self._grain_radius
        g = 0.89
        return iops(k_abs, k_sca, g)    

       
class two_stream_white_on_black:
    """
    Two-stream approximation for layer of white particles on black bottom
    E.g. Bohren, C.F. and Clothiaux, E.E., 2006. Fundamentals of
    atmospheric radiation: an introduction with 400 problems. John
    Wiley & Sons.
    """
    def __init__(self, iops, thickness):
        self._iops = iops
        self._thickness = thickness
        
    def reflectivity(self):
        tau_s = self._iops.scaled_optical_depth(self._thickness)  
        return tau_s / (2+tau_s)
    
    def transmissivity(self):
        tau_s = self._iops.scaled_optical_depth(self._thickness)  
        return 2 / (2+tau_s)
    
class grey_on_grey:
    # Implementation based on
    # Kokhanovsky, A.A., 2021. Snow optics. Springer.

    def __init__(self, iops, thickness, bottom_albedo):
        self._iops = iops
        self._L = thickness
        self._A = bottom_albedo

    def reflection_function(self, ksi, eta):
        # Reflection function
        # R, eq 3.199
        return 0

    def transmission_function(self, ksi, eta):
        # Transmission function
        # ksi is cosine of source zenith angle
        # eta is absolute value of cosine of observation angle
        # T, eq 3.201
        A = self._A
        r = self._diffuse_reflectance_black()
        r_p = r
        t = self._diffuse_transmittance_black()
        return self._diffuse_transmittance_black() + A*r_p*t/(1-A*r)
        #return self._transmission_function_black(ksi, eta) + A*r_p*t/(1-A*r)
        
    def diffuse_reflectance(self):
        # Diffuse reflectance approximation of, eq 3.201
        A = self._A
        r = self._diffuse_reflectance_black()
        t = self._diffuse_transmittance_black()
        return r + A*t**2/(1-A*r)

    def diffuse_transmittance(self):
        # Diffuse transmittance approximation of, eq 3.201
        A = self._A
        r = self._diffuse_reflectance_black()
        t = self._diffuse_transmittance_black()
        return t/(1-A*r)

    def _transmission_function_black(self, ksi, eta):
        # Transmission function for zero bottom albedo and a weakly absorbing layer
        # T_black, eq 3.191
        x = self._x()
        y = self._y()
        return np.sinh(y)/np.sinh(x+self._b()*y)*self._escape(ksi)*self._escape(eta)
    
    def _diffuse_reflectance_black(self):
        # Diffuse reflectance (spherical albedo) for zero bottom albedo. Isotropic incidence.
        # r, eq 3.169 
        x = self._x()
        y = self._y()
        z = x+y
        if z <= 0:
            return 0
        return np.sinh(x)/np.sinh(x+y)
    
    def _diffuse_transmittance_black(self):
        # Diffuse transmittance for zero bottom albedo. Isotropic incidence.
        # t, eq 3.170
        x = self._x()
        y = self._y()
        z = x + y;
        if z<=0:
            return 1
        return np.sinh(y)/np.sinh(x+y)
    
    def _escape(self, mu_0):
        # Escape function for Henyey-Greenstein with g=0.75
        # mu_0 is cosine of source zenith angle
        # u_0, eq 3.106
        return 3/5*mu_0+(1+np.sqrt(mu_0))/3
    
    def _r_0(self):
        # below eq 3.167
        return np.exp(-self._y())
    
    def _diffusion_exponent(self):
        # k, eq 3.172
        omega_0 = self._iops.single_scattering_albedo()
        g = self._iops.asymmetry_factor()
        
        k = self._iops.extinction_coefficient()
        #print("h = ",self._L)
        #print("g = ",g)
        #print("sscoalb = ",1-omega_0)
        #print("b = ",omega_0*k)
        return np.sqrt(3*(1-g)*(1-omega_0))
    
    def _similarity(self):
        # s, eq 3.69
        beta = self._absorption_probability()
        return np.sqrt(beta/(3*(1-self._iops.asymmetry_factor())))
    
    def _absorption_probability(self):
        # beta, eq 3.44
        return 1-self._iops.single_scattering_albedo()

    def _b(self):
        # eq 3.186
        return 1.072
        
    def _gamma(self):
        # below eq 3.167
        k = self._diffusion_exponent()
        k_ext = self._iops.extinction_coefficient()
        return k * k_ext
        
    def _x(self):
        # eq 3.175
        return self._gamma()*self._L
        
    def _y(self):
        # below eq 3.172
        return 4*self._similarity()
    
        
"""
refractive_index = 1.31+1e-9j
#print(np.imag(refractive_index))
iops = low_absorption_spheres(volume_fraction=0.2,
                              grain_radius=300e-6,
                              refractive_index=1.31+1e-9j,
                              wavelength=500e-9).get_iops()

thickness = 0.1

wob = two_stream_white_on_black(iops, thickness)
gog = grey_on_grey(iops, thickness, bottom_albedo=0.001)


#print(iops.scaled_optical_depth(2))
print(wob.reflectivity())
print(wob.transmissivity())
print(gog._diffuse_reflectance_black())
print(gog._diffuse_transmittance_black())
print(" ")
print(wob.reflectivity()+wob.transmissivity())
print(gog._diffuse_transmittance_black()+gog._diffuse_reflectance_black())
print(" ")

print(gog.diffuse_reflectance()/wob.reflectivity())
print(gog.diffuse_transmittance()/wob.transmissivity())
print(" ")

print(gog.diffuse_reflectance())
print(gog.diffuse_transmittance())
print(gog.diffuse_reflectance()+gog.diffuse_transmittance())

"""
