"""
Shows how the asymmetry factor of snow varies with number of volume scattering
function expansion terms
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
from matplotlib import cm
import flick

n_terms = [10, 30]
sub_command = " 500e-9 500e-9 1 snow full_mie 0.2 5e-6 0"
m = [] 
for i in range(len(n_terms)):
    xy = flick.run("iop wigner_alpha_beta_" + str(n_terms[i])+sub_command)
    m.append(xy)

g = []
for i in range(len(m)):
    alpha = m[i][:,0]
    g.append(alpha[1]/alpha[0]/3)

m = [] 
for i in range(len(n_terms)):
    xy = flick.run("iop scattering_ab_fitted_" + str(n_terms[i])+sub_command)
    m.append(xy)

exact = flick.run("iop scattering_ab_500_" + str(n_terms[i])+sub_command)

fig, ax = plt.subplots() 
line = []
ax.semilogy(exact[:,0],exact[:,1], linewidth=2,label="exact")
for i in range(len(m)):
    x = m[i][:,0]
    vsf = m[i][:,1]
    ax.semilogy(x, vsf, linewidth=1.5, label="N = "+ \
                str(n_terms[i])+", g="+f"{g[i]:.3f}")

ax.legend()
ax.grid()
ax.set_title("Small snow grains with radius 0.05 mm and volume fraction 0.2")
ax.set_xlabel("Cosine of scattering angle")
ax.set_ylabel("Volume scattering function [m$^{-1}$sr$^{-1}$]")
plt.show();
fig.savefig("snow_asymmetry_factor.pdf", bbox_inches='tight')
