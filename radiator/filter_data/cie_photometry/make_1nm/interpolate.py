import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy import integrate

m = np.loadtxt('cie_lms_cf_2deg_no_header_5nm.txt')

wl1 = 370
wl2 = 830
n = (wl2-wl1)+1
xnew = np.linspace(wl1,wl2,n)
d = []
for i in range(3):
    x = m[:,0]
    y = m[:,i+1]
    idx = np.where(y>0)
    x = x[idx]
    y = y[idx]
    spl = CubicSpline(x,np.log(y))
    ynew = np.exp(spl(xnew))
    if i==2:
        idx = np.where(xnew > 700)
        print(idx)
        ynew[idx] = ynew[idx[0][0]]
        print(ynew[idx[0][0]])
    plt.semilogy(m[:,0],m[:,1+i],'o')
    plt.semilogy(xnew,ynew,'-')
    d.append(ynew)
plt.show()

data = np.stack([xnew, d[0], d[1], d[2]], axis=1)
np.savetxt('cie_lms_cf_2deg_no_header_1nm.txt', data, fmt='%.0f  %.4e  %.4e  %.4e')

