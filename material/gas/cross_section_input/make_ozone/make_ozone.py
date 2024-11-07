import numpy as np
import matplotlib.pyplot as plt
import math

# Compressed ozone molecular absorption cross section spectra based on
# extensive measurements by:

# Gorshelev, V., Serdyuchenko, A., Weber, M., Chehade, W. and Burrows,
# J.P., 2014. High spectral resolution ozone absorption
# cross-sections–Part 1: Measurements, data analysis and comparison
# with previous measurements around 293 K. Atmospheric Measurement
# Techniques, 7(2), pp.609-624.

# Data from:
# https://doi.org/10.5281/zenodo.5793207

T = np.array([193,203,213,223,233,243,253,263,273,283,293])
C_min = 1e-29
T_0 = 293

def read_o3_absorption_file(filename):
    header = []
    data = []
    header_lines = 45
    with open(filename, 'r') as file:
        for _ in range(header_lines):
            line = file.readline().strip()
            header.append(line)
        for line in file:
            data_row = list(map(float, line.strip().split()))
            data.append(data_row)
    return header, data

def insert_header(file_name, header_file_name):
    with open(file_name, 'r') as file:
        original_content = file.read()
    with open(header_file_name, 'r') as file:
        header = file.read()
    with open(file_name, 'w') as file:
        file.write(header + '\n' + original_content)
        
def temperature_dependence(C):
    if any(x<=C_min for x in C):
        lC_0 = np.mean(C)
        dlC_dT = 0
        if lC_0 < C_min:
            lC_0 = C_min
        lC_0 = np.log(lC_0)
    else:    
        p = np.polyfit(T-T_0, np.log(C), 1)
        lC_0 = p[1]
        dlC_dT = p[0]
    return lC_0, dlC_dT

def compress(x,y,width):
    n_points = int(5*(x[-1]-x[0])/width)+2
    print(n_points)
    x_new = np.linspace(x[0],x[-1],n_points)
    y_new = np.zeros(len(x_new))
    for i in range(len(x_new)):
        if i == 0:
            x_low = x_new[0]
            x_high = x_new[2]
        elif i == len(x_new)-1:
            x_low = x_new[-3]
            x_high = x_new[-1]
        else:
            x_low = x_new[i-1]
            x_high = x_new[i+1]
        i_low = np.argmin(np.abs(x - x_low))
        i_high = np.argmin(np.abs(x - x_high))
        pf = np.polyfit(x[i_low:i_high], y[i_low:i_high], 3)
        p = np.poly1d(pf)
        y_new[i] = p(x_new[i])
        if math.isnan(y_new[i]):
            y_new[i] = np.log(C_min)
    return x_new, y_new

header, data = read_o3_absorption_file("SerdyuchenkoGorshelev5digits_latest.dat")
data = np.array(data)
n_rows = 88615
wls = np.zeros(n_rows)
C_T0 = np.zeros(n_rows)
slope = np.zeros(n_rows)
to_m2 = 1e-4
to_m = 1e-9
for row in range(n_rows):
    wls[row] = data[row,0]*to_m
    C = np.flipud(data[row,1:])*to_m2
    lC_0, dlC_dT = temperature_dependence(C)
    C_T0[row] =  np.exp(lC_0)
    slope[row] = dlC_dT


spectral_width = 0.001
x1, y1 = compress(np.log(wls), np.log(C_T0), spectral_width)
x2, y2 = compress(np.log(wls), slope, spectral_width)

wls_c = np.exp(x1)
C_T0_c = np.exp(y1)
slope_c = y2

d1 = np.stack([wls_c, C_T0_c], axis=1)
d2 = np.stack([wls_c, slope_c], axis=1)
file1 = "ozone_cross_section_293K.txt"
file2 = "ozone_temperature_dependence.txt"
np.savetxt(file1, d1, fmt='%.4e   %.3e')
np.savetxt(file2, d2, fmt='%.4e   %.2e')

insert_header(file1, "header1.txt")
insert_header(file2, "header2.txt")

T0 = 293
T1 = 293
T2 = 193
plt.loglog(wls*1e9,C_T0*np.exp(slope*(T1-T0))*1e4,label='Original, T = 293 K')
plt.loglog(wls_c*1e9,C_T0_c*np.exp(slope_c*(T1-T0))*1e4,label='Compressed, T = 293 K',linewidth=0.5)
plt.loglog(wls_c*1e9,C_T0_c*np.exp(slope_c*(T2-T0))*1e4,label='Compressed, T = 193 K',linewidth=0.5)
plt.xlabel('Wavelength [nm]')
plt.ylabel(r'Molecular ozone cross section [cm$^2$]')
#plt.grid()
plt.grid(which='both')
plt.legend()
plt.show()
