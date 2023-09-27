# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 00:41:12 2023

@author: Goldie
"""

import streamlit as st
import streamlit.components.v1 as components
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import mpld3

#Constants for model
#Current material parameters for NMP centrifugation
pS = 1595           #surfactant layer density, from Hirsham using anhydrous cholate on surface
d0 = 4.25E-10       #surfactant layer thickness, from Hirsham

def supernatant_frac(areas, layer_number, time_hour, rpm):
    "Return the fraction of a flake population remaning in the supernatant"
    r1, r2 = 0.071, 0.1
    time_s = 3600 * time_hour
    h = layer_number * d1
    w = (rpm * 2 * np.pi) / 60
    a = time_s * w * w * (h * (pNS - pL) + 2 * d0 * (pS - pL )) * np.sqrt(areas * 1e-18) / (12 * n * np.cbrt(3 / (4 * np.pi)))
    return max((r2*np.exp(-a)-r1)/(r2-r1), 0)

vector_supernat = np.vectorize(supernatant_frac)

def cutoff_area(layer_number, time_hour, rpm):
    r1, r2 = 0.071, 0.1
    time_s = 3600 * time_hour
    h = layer_number * d1
    w = (rpm * 2 * np.pi) / 60
    return (((np.log(r2 / r1) * (12 * n * np.cbrt(3 / (4 * np.pi))) / (time_s * w * w * (h * (pNS - pL) + 2 * d0 * (pS - pL)))) ** 2) * 1e18)

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
material = st.selectbox('Pick a material:',['WS2','MoS2','Graphene','NiPS3','MnPS3','GaS','BN','Ni(OH)2','Cu(OH)2','RuCl3','CrTe3'])
solvent = st.selectbox('Pick a solvent:',['Water','NMP','CHP','IPA'])
rpm_lower = st.slider('Select rpm of supernatant retention (use left/right arrow keys for fine adjustment)', min_value=100, max_value=30000)
rpm_higher = st.slider('Select rpm of sedimentation retention (use left/right arrow keys for fine adjustment)', min_value=100, max_value=30000)
time_hour = st.slider('Select time in minutes (use left/right arrow keys for fine adjustment)', min_value=10, max_value=360, value=120) / 60

if material == 'Graphene':
    d1 = 3.35E-10   #monolayer thickness
    k1 = 2.5        #L/W
    k2 = 280        #L/h
    pNS = 2260      #structure density
if material == 'WS2':
    d1 = 6.18E-10    #monolayer of WS2, J. A. Wilson and A. D. Yoffe, Advances in Physics 18, 193 (1969).
    k1 = 1.75       #L/W
    k2 = 46         #L/h
    pNS = 7500      #structure density
if material == 'MoS2':
    d1 = 6.15E-10   #monolayer thickness, J. A. Wilson and A. D. Yoffe, Advances in Physics 18, 193 (1969).
    k1 = 1.8       
    k2 = 45         
    pNS = 5060      #structure density
if material == 'NiPS3':
    d1 = 6.6E-10
    pNS = 3180
if material == 'MnPS3':
    d1 = 6.8E-10
    pNS = 2900
if material == 'GaS':
    d1 = 7.5E-10
    pNS = 3860
if material == 'BN':
    d1 = 3.4E-10
    pNS = 2100
if material == 'NiOH2':
    d1 = 4.6E-10
    pNS = 4100
if material == 'CuOH2':
    d1 = 5E-10
    pNS = 3370
if material == 'RuCl3':
    d1 = 6E-10
    pNS = 3260
if material == 'CrTe3':
    d1 = 1.1E-9
    pNS = 4700    

if solvent == 'Water':
    n = 1.3E-3
    pL = 1000
if solvent == 'NMP':
    n = 2.3E-3
    pL = 1030
if solvent == 'CHP':
    n = 2.1E-3
    pL = 1000
if solvent == 'IPA':
    n = 3.31E-3
    pL = 785

#Define the plot parameters and apply the model to the plot
dummy_layer_numbers = np.arange(start=0.5, stop=50, step=0.1)
dummy_areas = np.arange(0, 500000, 1000)

X, Y = np.meshgrid(dummy_layer_numbers, dummy_areas)

f1 = vector_supernat(Y,X,time_hour,rpm_lower)
f2 = vector_supernat(Y,X,time_hour,rpm_higher)

Z = (1-f2)*f1

cutoff_values = [cutoff_area(m, time_hour, rpm_higher) for m in dummy_layer_numbers]
    
#fig = plt.figure()
#plt.title('Predicted Cut-Off Size', loc='center')
#plt.plot(dummy_layer_numbers, cutoff_values, 'k--', label='Cut-off')
#plt.ylabel("Area",fontsize=14)
#plt.xlabel("Layer Number",fontsize=14)
#plt.ylim(0, 5e+5)
#plt.tight_layout()
#plt.legend()

fig, ax = plt.subplots()
CS = ax.contour(X, Y, Z, 6, cmap=cm.coolwarm)
ax.set(xlim=(0, 20), ylim=(0, 100000),  xlabel="Layer Numbers", ylabel="Area")
plt.tight_layout()
st.subheader('Flake Size Probability Distribution')
st.write('Discounting the starting distribution, which will bias all separations towards \
         smaller flakes, since small flakes dominate most LPE samples, we can produce a \
         probability plot showing the relative population of each flake size expected for given centrifuge conditions. \
         The contour plot shows a clear band of preferential separation with cooler, blue colours showing the decay \
         in probability of flakes outside the preferred size range.')
st.pyplot(fig)

           
