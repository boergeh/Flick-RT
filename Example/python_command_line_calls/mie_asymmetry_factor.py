"""
Shows that the asymmetry factor is growing very slowly with radius for large spheres
"""
import numpy as np
from scipy.special import roots_legendre
import matplotlib.pyplot as plt
import os

n_points = 2*2048
roots, weights = roots_legendre(n_points)
radius = np.logspace(np.log10(1e-6),np.log10(300e-6),100)
print(radius)
g = []
for r in radius:
    s = "flick mie 500e-9 1.0 1.33 "+str(r)+" 0 0.1 scattering_matrix_element 0 0 "+str(n_points)+" > xy_tmp.txt"
    os.system(s)
    xy = np.loadtxt("xy_tmp.txt")
    mu = np.cos(xy[:,0])
    p = xy[:,1]
    g.append(np.sum(weights*p*mu)/np.sum(weights*p))
    print(r)

fig, ax = plt.subplots()
ax.semilogx(radius*1e6, g,'b-')
ax.grid()
ax.set_title("full mie, $\lambda=500$ nm, refractive index = 1.33")
ax.set_xlabel("radius [$\mu m$]")
ax.set_ylabel("asymmetry factor")
ax.set_ylim([0.65, 1])
ax.set_xticks([1,3,10,30,100,300])
plt.show();
#fig.savefig("mie_asymmetry_factor.pdf", bbox_inches='tight')


