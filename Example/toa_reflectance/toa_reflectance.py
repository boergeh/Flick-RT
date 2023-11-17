"""
Shows top-of-atmosphere nadir reflectance [1/sr] above the ocean.
"""
import numpy as np
import matplotlib.pyplot as plt
import os

wavelengths = np.linspace(400e-9,700e-9,5);
if not os.path.exists("./config"):
    os.system("flick accurt -g toa_reflectance config")
    
wls_str = "\""+" ".join(str(x) for x in wavelengths)+"\""
command = "flick text config set DETECTOR_WAVELENGTHS "+wls_str+" > new_config"
os.system(command)
print("calculating ...")
os.system("flick accurt new_config > tmp.txt")
r = np.loadtxt("tmp.txt")
x = r[:,0]*1e9
y = r[:,1]

fig, ax = plt.subplots()
ax.plot(x, y,linewidth=1)
ax.grid()
ax.set_title("Top-of-atmosphere nadir reflectance")
ax.set_ylabel("$L_u / E_d$  [sr$^{-1}$]")
ax.set_xlabel("Wavelength [nm]")
plt.show();

