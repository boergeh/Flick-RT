import numpy as np
import matplotlib.pyplot as plt
import os

#os.system('make_and_test')
# run ./make_and_test before running this python script

h = np.loadtxt("histogram_view.txt")
x_center = h[1:,0]
y_center = h[0,1:]
z = h[1:,1:]

dx = (x_center[1]-x_center[0])
dy = (y_center[1]-y_center[0])
x = x_center-dx/2;
y = y_center-dy/2;
x = np.append(x,x[-1]+dx)
y = np.append(y,y[-1]+dy)

fig, ax = plt.subplots()
im = ax.pcolormesh(x,y,np.transpose(z))
if np.size(x) < 5:
    ax.set_xticks(x)
    ax.set_yticks(y)
ax.grid()
ax.set_xlabel("x-bins")
ax.set_ylabel("y-bins")
ax.set_title("histogram_view")
fig.colorbar(im)
plt.show()

