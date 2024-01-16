# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 00:41:12 2023

@author: Goldie
"""

import streamlit as st
import streamlit.components.v1 as components
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.style as mplstyle
import numpy as np
import mpld3
import time

#Constants for model
#Current material parameters for NMP centrifugation


pS = 1595           #surfactant layer density, from Hirsham using anhydrous cholate on surface
d0 = 4.25E-10       #surfactant layer thickness, from Hirsham

st.title('Modelling Centrifugation of 2D Nanomaterials')
st.markdown('*citation goes here*')
st.subheader('IN DEVELOPEMENT')
st.write('This interactive tool allows you to visualise the changes to a nanomaterial during centrifugation.\
             By changing the experiment parameters like rpm, centrifuge time and rotor dimensions the live plots \
             will display the relative populations of different sizes flakes within different fractions. \
             Currently this tool does not account for starting dispersions. This is particularly significant \
             because the form of equation extends symmetrically out to extreme aspect ratios; however, \
             such flake dimensions are unrealistic in real nanoflake populations')
st.write('The complete theoretical considerations will be published in a manuscript currently in preparation. \
              For a basic overview however, the rate of sedimentation within a centrifugal field can be defined according \
              to nanoflake size. The relative movement of each nanoflake within a dispersions can be then compared to \
              the tube dimensions. This produces a function defining the relative population of each nanoflake remaining \
              in the supernatant.')
st.latex(r'''F_n = \begin{cases} \frac{r_2 e^{-a} - r_1}{r_2 - r_1}, & \text{if } a > 0 \\ 0, & \text{if } a \leq 0 \end{cases}''')

st.write('where')
st.latex(r'''a = \frac{t \omega^2 (N d_1 \Delta\rho_{NS} + 2 d_0 \Delta\rho_S) \sqrt{\frac{Lw}{2}}}{12 \eta \sqrt[3]{\frac{3}{4\pi}}}''')

st.markdown(
'''Defining all parameters:
- $r_1 $ = axial distance from center of rotation to top of the liquid surface $/ m$
- $r_2 $ = axial distance from center of rotation to bottom of the centrifuge tube $/ m$
- $t $ = time $/ s$
- $\omega $ = angular velocity $/ rad \space s^{-1}$
- $N $ = number of layers
- $d_1 $ = monolayer thickness $m$
- $ \Delta p_{NS} $ = difference in density between nanosheet and liquid $/ kg \space m^{-3}$
- $d_0 $ = surfactant layer thickness $/ m$
- $\Delta p_S $ = difference in density between surfactant layer and liquid $/ kg \space m^{-3}$
- $L $ = nanosheet length $/ m$
- $w $ = nanosheet width $/ m$
- $\eta $ = liquid viscosity $/ kg \space m^{-1} \space s^{-1}$)
'''
    )

st.caption('Caption text to go underneath the figures being plotted')
#Imput user parameters
material_input = st.selectbox('Pick a material:',['WS2','MoS2','Graphene','NiPS3','MnPS3','GaS','BN','Ni(OH)2','Cu(OH)2','RuCl3','CrTe3'])
solvent_input = st.selectbox('Pick a solvent:',['Water','NMP','CHP','IPA'])
rpm_lower = st.slider('Select rpm of supernatant retention (use left/right arrow keys for fine adjustment)', min_value=100, max_value=30000)
rpm_higher = st.slider('Select rpm of sedimentation retention (use left/right arrow keys for fine adjustment)', min_value=100, max_value=30000)
time_hour = st.slider('Select time in minutes (use left/right arrow keys for fine adjustment)', min_value=10, max_value=360, value=120) / 60


material_properties = {
    'Graphene':[2260, 3.35E-10, 2.5, 280],
    'WS2': [7500, 6.18E-10, 1.75, 46],
    'MoS2': [5060, 6.15E-10, 1.8, 45],
    'NiPS3': [3180, 6.6E-10, 2.0, 45],
    'MnPS3': [2900, 6.8E-10, 2.0, 45],
    'GaS': [3860, 7.5E-10, 2.0, 45],
    'BN': [2100, 3.4E-10, 2.0, 45],
    'Ni(OH)2': [4100, 4.6E-10, 2.0, 45],
    'Cu(OH)2': [3370, 5E-10, 2.0, 45],
    'RuCl3': [3260, 6E-10, 2.0, 45],
    'CrTe3': [4700, 1.1E-9, 2.0, 45]
    }

solvent_properties = {
    'Water': [1.3E-3, 1000],
    'NMP': [2.3E-3, 1030],
    'CHP': [2.1E-3, 1000],
    'IPA': [3.31E-3, 785]
    }

pNS = material_properties[material_input][0]
d1 = material_properties[material_input][1]

n = solvent_properties[solvent_input][0]
pL = solvent_properties[solvent_input][1]

def supernatant_fraction_arrays(layer_number, lateral_size, time_hour, rpm):
    "Produces an array of fraction remaining in the supernatant from arrays of thickness and lateral size"
    r1, r2 = 0.071, 0.1
    time_s = 3600 * time_hour
    w = (rpm * 2 * np.pi) / 60
    experiment_constant = time_s*w*w/(12*n*np.cbrt(3/(4*np.pi)))
    surfactant_density_term = 2*d0*(pS-pL)
    nanosheet_density_term = pNS-pL
    x_array = [experiment_constant*(surfactant_density_term + nanosheet_density_term * x * d1) for x in layer_number]
    velocity_matrix = np.outer(lateral_size*1e-9/2, x_array)
    fraction_matrix = np.maximum((r2*np.exp(-velocity_matrix)-r1)/(r2-r1), 0)
    return(fraction_matrix)

def cutoff_area(layer_number, time_hour, rpm):
    r1, r2 = 0.071, 0.1
    time_s = 3600 * time_hour
    h = layer_number * d1
    w = (rpm * 2 * np.pi) / 60
    return (((np.log(r2 / r1) * (12 * n * np.cbrt(3 / (4 * np.pi))) / (time_s * w * w * (h * (pNS - pL) + 2 * d0 * (pS - pL)))) ** 2) * 1e18)

time_4 = time.time()

#Define the plot parameters and apply the model to the plot
dummy_layer_numbers = np.arange(start=0.5, stop=20.1, step=0.1)
dummy_lateral_size = np.arange(start=10, stop=1501, step=1)

X, Y = np.meshgrid(dummy_layer_numbers, dummy_lateral_size)

f1 = supernatant_fraction_arrays(dummy_layer_numbers, dummy_lateral_size, time_hour, rpm_lower)
f2 = supernatant_fraction_arrays(dummy_layer_numbers, dummy_lateral_size, time_hour, rpm_higher)

time_5 = time.time()
print(f"Fractions calculation: {time_5 - time_4}")

Z = (1-f2)*f1

time_6 = time.time()
print(f"Final calculation: {time_6 - time_5}")

#cutoff_values = [cutoff_area(m, time_hour, rpm_higher) for m in dummy_layer_numbers]

time_7 = time.time()
print(f"Cutoff Calculation: {time_7 - time_6}")
    
mplstyle.use('fast')
fig, ax = plt.subplots()
CS = ax.contour(X, Y, Z, 6, cmap=cm.coolwarm)
ax.set(xlim=(0, 20), ylim=(0, 1500),  xlabel="Layer Number", ylabel="Lateral Size / nm")
plt.tight_layout()
st.subheader('Flake Size Probability Distribution')
st.write('Discounting the starting distribution, which will bias all separations towards \
         smaller flakes, since small flakes dominate most LPE samples, we can produce a \
         probability plot showing the relative population of each flake size expected for given centrifuge conditions. \
         The contour plot shows a clear band of preferential separation with cooler, blue colours showing the decay \
         in probability of flakes outside the preferred size range.')
st.pyplot(fig)

time_8 = time.time()
print(f"Plotting Time: {time_8 - time_7}")