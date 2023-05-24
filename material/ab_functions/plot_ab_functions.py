import numpy as np
import matplotlib.pyplot as plt
import os

fig, ax = plt.subplots(3,2,sharex=True)
fig.tight_layout(pad=2.0)
fig.set_size_inches(7.5, 7.5)
fig.subplots_adjust(hspace=0.1)

scaling_factor = np.loadtxt("scaling_factor.txt")
n_terms = int(np.loadtxt("n_terms.txt"))
p = np.loadtxt("a0.txt")
pf = np.loadtxt("a_fitted0.txt")
p = p[:,1];
pf = pf[:,1];

n = 0
for i in range(2):
    for j in range(2):
        if n==0:
            xy = np.loadtxt("a"+str(n)+".txt")
            ax[i,j].semilogy(xy[:,0], xy[:,1],'b-')
            xy = np.loadtxt("a_fitted"+str(n)+".txt")
            ax[i,j].semilogy(xy[:,0], xy[:,1],'r-')
            ax[i,j].set_ylabel("$a_"+str(n+1)+"$",math_fontfamily='dejavuserif')
            ax[i,j].legend(['original','fitted Wigner d-functions'])
            ax[i,j].text(-0.95,0.5,"fitted terms: "+str(n_terms) + "\nscaling factor: " +
                         str(round(scaling_factor*100)/100),
                         bbox={'facecolor': 'white', 'alpha': 0.6, 'pad': 0.5,
                               'edgecolor': 'grey','boxstyle':'round'})
        else:
            xy = np.loadtxt("a"+str(n)+".txt")
            ax[i,j].plot(xy[:,0], xy[:,1]/p,'b-')
            xy = np.loadtxt("a_fitted"+str(n)+".txt")
            ax[i,j].plot(xy[:,0], xy[:,1]/pf,'r-')
            ax[i,j].set_ylabel("$a_"+str(n+1)+"\, /\, a_1$",math_fontfamily='dejavuserif')    
            
        ax[i,j].grid()    
        n += 1

n=0        
for i in [2]:
    for j in range(2):
        xy = np.loadtxt("b"+str(n)+".txt")
        ax[i,j].plot(xy[:,0], xy[:,1]/p,'b-')
        xy = np.loadtxt("b_fitted"+str(n)+".txt")
        ax[i,j].plot(xy[:,0], xy[:,1]/pf,'r-')
        
        ax[i,j].grid()
        ax[i,j].set_ylabel("$b_"+str(n+1)+"\, /\, a_1$",math_fontfamily='dejavuserif')
        ax[i,j].set_xlabel("$\cos(\Theta)$")    
        n += 1

plt.show()
fig.savefig("plot_ab_functions.pdf", bbox_inches='tight')
