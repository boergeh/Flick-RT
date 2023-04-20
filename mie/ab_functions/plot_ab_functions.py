import numpy as np
import matplotlib.pyplot as plt
import os

fig, ax = plt.subplots(3,2,sharex=True)
fig.tight_layout(pad=0.7)
fig.set_size_inches(7.5, 7.5)

n = 0
for i in range(2):
    for j in range(2):
        xy = np.loadtxt("a"+str(n)+".txt")
        ax[i,j].plot(xy[:,0], xy[:,1],'b-')
        xy = np.loadtxt("a_scaled"+str(n)+".txt")
        ax[i,j].plot(xy[:,0], xy[:,1],'r-')
        xy = np.loadtxt("a_scaled_fitted"+str(n)+".txt")
        ax[i,j].plot(xy[:,0], xy[:,1],'g-')
        #ax[i,j].grid()
        ax[i,j].set_title("a$_"+str(n+1)+"$")
        n += 1

n=0        
for i in [2]:
    for j in range(2):
        xy = np.loadtxt("b"+str(n)+".txt")
        ax[i,j].plot(xy[:,0], xy[:,1],'b-')
        xy = np.loadtxt("b_scaled"+str(n)+".txt")
        ax[i,j].plot(xy[:,0], xy[:,1],'r-')
        xy = np.loadtxt("b_scaled_fitted"+str(n)+".txt")
        ax[i,j].plot(xy[:,0], xy[:,1],'g-')
        
        #ax[i,j].grid()
        ax[i,j].set_title("b$_"+str(n+1)+"$")
        n += 1

plt.show()
fig.savefig("plot_ab_functions.pdf", bbox_inches='tight')
