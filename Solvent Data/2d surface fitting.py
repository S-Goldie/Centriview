# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 22:35:31 2024

@author: Goldie
"""

#!/usr/bin/evn python

import numpy as np
import scipy.linalg
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# some 3-dim points
density_data = np.genfromtxt('IPA Density.csv', delimiter=',', skip_header=1)
viscosity_data = np.genfromtxt('IPA Viscosity.csv', delimiter=',', skip_header=1)

# regular grid covering the domain of the data
X,Y = np.meshgrid(np.arange(0, 1.01, 0.01), np.arange(0, 51, 1))
XX = X.flatten()
YY = Y.flatten()

def surface_plot(data, order, label):    
    if order == 1:
        # best-fit linear plane
        A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
        C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])    # coefficients
        
        # evaluate it on grid
        Z = C[0]*X + C[1]*Y + C[2]
        
        # or expressed using matrix/vector product
        #Z = np.dot(np.c_[XX, YY, np.ones(XX.shape)], C).reshape(X.shape)
    
    elif order == 2:
        # best-fit quadratic curve
        A = np.c_[np.ones(data.shape[0]), data[:,:2], np.prod(data[:,:2], axis=1), data[:,:2]**2]
        C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])
        
        # evaluate it on a grid
        Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY**2], C).reshape(X.shape)
    
    elif order == 3:
        # M = [ones(size(x)), x, y, x.^2, x.*y, y.^2, x.^3, x.^2.*y, x.*y.^2, y.^3]
        A = np.c_[np.ones(data.shape[0]), data[:,:2], data[:,0]**2, np.prod(data[:,:2], axis=1), \
                  data[:,1]**2, data[:,0]**3, np.prod(np.c_[data[:,0]**2,data[:,1]],axis=1), \
                  np.prod(np.c_[data[:,0],data[:,1]**2],axis=1), data[:,1]**3]
        C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])
        
        Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX**2, XX*YY, YY**2, XX**3, XX**2*YY, XX*YY**2, YY**3], C).reshape(X.shape)
    
    # plot points and fitted surface
    fig = plt.figure(label)
    ax = fig.add_subplot(projection = '3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
    ax.scatter(data[:,0], data[:,1], data[:,2], c='r', s=50)
    plt.xlabel('IPA wt. Content')
    plt.ylabel('Temperature / $^oC$')
    ax.set_zlabel(label)
    ax.axis('equal')
    ax.axis('tight')
    plt.show()
    
    return(C)

def property_function(w, T, c):
    return(c[0] + w*c[1] + T*c[2] + (w**2)*c[3] + w*T*c[4] + (T**2)*c[5] + (w**3)*c[6] + (w**2)*T*c[7] + w*(T**2)*c[8] + (T**3)*c[9])

density_c = surface_plot(density_data, 3, 'Density')
viscosity_c = surface_plot(viscosity_data, 3, 'Viscosity')

IPAWaterDensity = [0.999581, -0.0867712, -0.000109338, -0.171483, -0.00165038, -1.24842e-06, 0.0592866, 0.00105765, -2.63718e-07, -1.86795e-08]
IPAWaterViscosity = [1.78798, 14.9948, -0.0677034, -13.7734, -0.283823, 0.00133142, 1.61501, 0.194239, 0.000875951, -9.46676e-06]

weight_frac = 0.6
temp = 10

density_test = property_function(weight_frac, temp, density_c)
#literature values = 0.901
viscosity_test = property_function(weight_frac, temp, viscosity_c)
#literature values = 3.1

