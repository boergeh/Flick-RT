"""
Shows that the asymmetry factor is growing slowly with radius
for large spheres
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import sys
import os
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

n_points = 1e4
radius = np.logspace(np.log10(1e-6), np.log10(50e-6), 30)
print("radius [m]:")
g = []
for r in radius:
    s = "mie 500e-9 1.0 1.33 "+str(r)+ \
        " 0 1 scattering_matrix_element 0 0 "+str(n_points)
    xy = flick.run(s)
    mu = np.cos(xy[:,0])
    s00 = xy[:,1]
    x = mu[::-1]
    f = s00[::-1]
    g.append(integrate.simpson(x*f,x=x)/integrate.simpson(f,x=x))
    print(r)

fig, ax = plt.subplots()
ax.semilogx(radius*1e6, g,'b-')
ax.grid()
ax.set_title(r"full mie, $\lambda=500$ nm, refractive index = 1.33")
ax.set_xlabel(r"radius [$\mu m$]")
ax.set_ylabel("asymmetry factor")
ax.set_xlim([0.8, 120])
ax.set_ylim([0.65, 1])
ax.set_xticks([1,3,10,30,100])
plt.show();



