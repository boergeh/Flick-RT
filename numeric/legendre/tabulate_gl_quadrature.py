import numpy as np
from scipy.special import roots_legendre

for i in range(14):
    n = 2**i
    roots, weights = roots_legendre(n)
    filename = "./gl_quadrature/q_"+str(n)+".txt"
    data = np.column_stack([roots, weights])
    np.savetxt(filename, data, fmt=['%0.15e','%0.15e'])
