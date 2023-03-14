import numpy as np
import matplotlib.pyplot as plt
import os

fig, ax = plt.subplots(2,2,sharex=True)
fig.tight_layout(pad=2.3)
fig.set_size_inches(7.5, 5)


n_terms = [[8, 16],
           [32, 64]]

for i in range(2):
    for j in range(2):
        os.system('flick text xy ./petzold_phase_function.txt > tmp.txt')
        xy = np.loadtxt("tmp.txt")
        ax[i,j].semilogy(xy[:,0], xy[:,1],'b-',linewidth=1,label="Petzold")
        
        os.system("flick delta_fit "+str(n_terms[i][j])+
                  " function_values petzold_phase_function.txt > tmp.txt")
        xy = np.loadtxt("tmp.txt")
       
        ax[i,j].semilogy(xy[:,0], xy[:,1],'r-',linewidth=1,label="delta-fit")
        ax[i,j].legend(loc='upper left')
        ax[i,j].set_xlim(-1.02, 1.02)
        ax[i,j].set_ylim(4e-4, 3e3) 
        ax[i,j].grid()
        ax[i,j].set_title(str(n_terms[i][j])+" terms")
        if i == 1:
            ax[i,j].set_xlabel("Cosine of scattering angle")
        if j == 0:
            ax[i,j].set_ylabel("Phase function [sr$^{-1}$]")
        
os.system('rm -f ./tmp.txt')
plt.show()
fig.savefig("delta_fit_plot.pdf", bbox_inches='tight')
