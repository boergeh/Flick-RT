"""
Shows how the asymmetry factor of snow varies with number of volume scattering
function expansion terms
"""
import numpy as np
import matplotlib.pyplot as plt
import os

n_terms = [64, 32, 16, 8]
sub_command = " 500e-9 500e-9 0.2 snow fast_mie 0.2 1000e-6 0 > tmp.txt"
m = [] 
for i in range(len(n_terms)):
    os.system("flick iop wigner_alpha_beta_" + str(n_terms[i])+sub_command)
    m.append(np.loadtxt("tmp.txt"))

g = []
for i in range(len(m)):
    alpha = m[i][:,0]
    g.append(alpha[1]/alpha[0]/3)

m = [] 
for i in range(len(n_terms)):
    os.system("flick iop scattering_ab_500_" + str(n_terms[i])+sub_command)
    m.append(np.loadtxt("tmp.txt"))

fig, ax = plt.subplots() 
line = []     
for i in range(len(m)):
    x = m[i][:,0]
    vsf = m[i][:,1]
    ax.semilogy(x, vsf, linewidth=2,label="g="+f"{g[i]:.3f}"+", N="+str(n_terms[i]))

ax.legend()
ax.grid()
ax.set_title("Snow grains with radius 1 mm and volume fraction 0.2  ")
ax.set_xlabel("Cosine of scattering angle")
ax.set_ylabel("Volume scattering function [m$^{-1}$sr$^{-1}$]")
plt.show();
#fig.savefig("snow_asymmetry_factor.pdf", bbox_inches='tight')
