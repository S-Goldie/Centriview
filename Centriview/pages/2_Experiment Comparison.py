# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 14:24:38 2024

@author: Stuart Goldie
"""

import streamlit as st
import streamlit.components.v1 as components
import matplotlib.pyplot as plt
import numpy as np
import mpld3
from mpld3 import plugins

###BACK END###
#Storage of parameters for calculations
#"Keys" : ['Density / kg m-3', 'Layer Thickness / m', 'k1', 'k2']
MATERIALS = {
    "Graphene" :[2260, 3.35E-10, 2.5, 280],
    "WS2" :     [7500, 6.18E-10, 1.75, 46],
    "MoS2" :    [5060, 6.15E-10, 1.8, 45],
    "NiPS3" :   [3180, 6.6E-10, 2.0, 45],
    "MnPS3" :   [2900, 6.8E-10, 2.0, 45],
    "GaS" :     [3860, 7.5E-10, 2.0, 45],
    "BN" :      [2100, 3.4E-10, 2.0, 45],
    "NiOH2" :   [4100, 4.6E-10, 2.0, 45],
    "CuOH2" :   [3370, 5E-10, 2.0, 45],
    "RuCl3" :   [3260, 6E-10, 2.0, 45],
    "CrTe3" :   [4700, 1.1E-9, 2.0, 45]
    }

#"Keys" : ['density_m', 'density_c', 'viscosity_A', 'viscosity_Ea']
#Constants produce values in common units, not SI. Conversion done in the function
SOLVENTS = {
    "IPA" :     [-8.97149E-4, 0.80402, 3.05818E-4, 21846.037],
    "THF" :     [-0.00108, 0.90959, 0.0445145, 5845.749],
    "NMP" :     [-8.81678E-4, 1.04999, 0.015821, 11553.57],
    "DLG (Cyrene™)" : [-9.73846E-4, 1.27159, 9.089E-04, 23410],
    "Toluene" : [-9.47951E-4, 0.88605, 0.0185636, 8438.939],
    "DMF" :     [-9.62251E-4, 0.96862, 0.021271, 9047.422],
    "Ethanol" : [-9.11871E-4, 0.80885, 0.003494, 14223.19],
    }

solvent_list = ['Water'] + list(SOLVENTS.keys()) + ['IPA:Water Mixture', 'Other']
materials_list = list(MATERIALS.keys()) + ['Other']

#Constants for two phase mixtures, fitted using a crude two-dimensional polynomial
IPAWaterDensity = [0.999581, -0.0867712, -0.000109338, -0.171483, -0.00165038, -1.24842e-06, 0.0592866, 0.00105765, -2.63718e-07, -1.86795e-08]
IPAWaterViscosity = [1.78798, 14.9948, -0.0677034, -13.7734, -0.283823, 0.00133142, 1.61501, 0.194239, 0.000875951, -9.46676e-06]

def solvent_viscosity(temp, solvent):
    """Returns the solvent viscosity in SI units given correct constants"""
    R = 8.314
    A = SOLVENTS[solvent][2]
    E_a = SOLVENTS[solvent][3]
    Kelvin = 273.15 + temp
    centipoise = A*np.exp(E_a/(R*Kelvin))
    return(centipoise/1000)

def solvent_density(temp, solvent):
    """Returns the solvent density in g cm-3"""
    m = SOLVENTS[solvent][0]
    c = SOLVENTS[solvent][1]
    return((m * temp + c)*1000)

#water follows different trends
def water_viscosity(temp):
    """Returns viscosity of water using third order exponential"""
    kelvin = temp+273.15
    A = 2.3168E-12
    B = 4729.66128
    C = 0.04265
    D = -2.14282E-5
    centipoise = A*np.exp(B/kelvin+C*kelvin+D*kelvin**2)
    return(centipoise / 1000)

def water_density(temp):
    """Returns the density of water using third order polynomial"""
    C = 0.99997
    B1 = 3.1087E-5
    B2 = -6.45401E-6
    B3 = 2.09654E-8
    return((C+B1*temp+B2*temp**2+B3*temp**3)*1000)

def mixed_property_function(w, T, c):
    return(c[0] + w*c[1] + T*c[2] + (w**2)*c[3] + w*T*c[4] + (T**2)*c[5] + (w**3)*c[6] + (w**2)*T*c[7] + w*(T**2)*c[8] + (T**3)*c[9])

def speed_equivalence(t_1, w_1, n_1, t_2, n_2, p_NS1, p_NS2, p_l1, p_l2):
    """Returns the angular rotation to match experiments with constant geometry"""
    k = (p_NS1-p_l1)/(p_NS2-p_l2)
    w2 = k * (t_1*n_2*w_1**2)/(t_2*n_1)
    return np.sqrt(w2)
    
#define functions for calculation of experiment coorespondance assuming constant or variable rotor geometry
#variable geometry is complex and requires the user specifies some variables to make it solvable

#once user inputs have been selected, define correct parameters and calculate experimental parameters
#plot time vs rpm for equivalent seperation (interactive plot) and highlight the minimum time at hgihest speed possible for the rotor selected


###FRONT END###
#Could do with deciding on a standard page banner with title, citation etc.
#Text outlining the basic idea behind the calculation and brief summary of the key parameters and operation of the calculation

#User inputs:
    #from experiment to match material, solvent, temperature, (mixture composition), material, change rotor?
    #if rotor swap then rotor and fill height before and after
    #intended experiment needs material, solvent, temperature, (mixture composition), rotor and fill height if changed

#output of interactive plot and minimum time at high speed

#discussion of equations with derivation to clarify operation of calculator and what all the terms mean

##_FRONT END CODE_##
st.title('Consistent Centrifugation Calculator')
st.subheader('CONFIDENTIAL')
st.markdown('_No information contained herein to be shared outside EU Project: 2D-PRINTABLE – GA No: 101135196_')
st.write('This calculation tool is designed to compare the separation expected during a centrifugation process, \
        and predict the experimental conditions required to replicate such separation when using different solvents, \
        materials or even different temperatures. Currently a working progress, these calculations \
        are underpinned by the theoretical work undertaken in the publication listed on the home page.')

st.subheader('Experiment Calculator')
st.write('To match experiment conditions, select the conditions from the previous, old experiment to be matched and \
        the conditions intended for the new experiment.\
        Check the required units carefully, this is a common error when comparing information from different sources.')


col1, col2 = st.columns(2)
col1.write("Original Experiment to be Replicated")
col2.write("New Experimental Conditions")

material1 = col1.selectbox('Material:',materials_list, key='material1selection')
material2 = col2.selectbox('Material:',materials_list, key='material2selection')

if material1 == 'Other':
    density_material1 = col1.number_input('Density of Custom Material in $g cm^{-3}$', key='matdensity1')*1000
else:
    density_material1 = MATERIALS[material1][0]
if material2 == 'Other':
    density_material2 = col2.number_input('Density of Custom Material in $g cm^{-3}$', key='matdensity2')*1000
else:
    density_material2 = MATERIALS[material2][0]

temp1 = col1.slider('Temperature (use left/right arrow keys for fine adjustment', min_value=5, max_value=50, key='temperature1')
temp2 = col2.slider('Temperature (use left/right arrow keys for fine adjustment', min_value=5, max_value=50, key='temperature2')
solvent1 = col1.selectbox('Solvent', solvent_list, key='oldsolvent')
solvent2 = col2.selectbox('Solvent', solvent_list, key='newsolvent')

if solvent1 == 'Water':
    density1 = water_density(temp1)
    viscosity1 = water_viscosity(temp1)
elif solvent1 == 'Other':
    density1 = col1.number_input('Density of Custom Solvent in $g cm^{-3}$', key='soldensity1')*1000
    viscosity1 = col1.number_input('Viscosity of Custom Solvent in $cP$', key='solvisc1')/1000
elif solvent1 == 'IPA:Water Mixture':
    composition1 = col1.slider('Weight Content of Alcohol', min_value=0.0, max_value=1.0, value=0.5, key='composition1')
    density1 = mixed_property_function(composition1, temp1, IPAWaterDensity)*1000
    viscosity1 = mixed_property_function(composition1, temp1, IPAWaterViscosity)/1000
else:
    density1 = solvent_density(temp1, solvent1)
    viscosity1 = solvent_viscosity(temp1, solvent1)

if solvent2 == 'Water':
    density2 = water_density(temp2)
    viscosity2 = water_viscosity(temp2)
elif solvent2 == 'Other':
    density2 = col2.number_input('Density of Custom Solvent in $g cm^{-3}$', key='soldensity2')*1000
    viscosity2 = col2.number_input('Viscosity of Custom Solvent in $cP$', key='solvisc2')/1000
elif solvent2 == 'IPA:Water Mixture':
    composition2 = col2.slider('Weight Content of Alcohol', min_value=0.0, max_value=1.0, value=0.5, key='composition2')
    density2 = mixed_property_function(composition2, temp2, IPAWaterDensity)*1000
    viscosity2 = mixed_property_function(composition2, temp2, IPAWaterViscosity)/1000
else:
    density2 = solvent_density(temp2, solvent2)
    viscosity2 = solvent_viscosity(temp2, solvent2)

  
time1 = col1.slider('Previous experiment time in minutes', min_value=10, max_value=540, value=120)
col1.markdown(f'Hours: {time1/60:.1f}')
rpm1 = col1.number_input('Previous rpm')

##_BACK END CALCULATIONS_##

w1 = (rpm1 * 2 * np.pi) / 60


time2 = np.arange(10, 540, 1)
w2 = speed_equivalence(time1, w1, viscosity1, time2, viscosity2, density_material1, density_material2, density1, density2)
rpm2 = w2*60/(2*np.pi)

speed_match_rads = speed_equivalence(time1, w1, viscosity1, time1, viscosity2, density_material1, density_material2, density1, density2)
speed_match = speed_match_rads*60/(2*np.pi)
st.markdown(f'__Speed needed to match same time : {speed_match:.0f} rpm__')

fig, ax = plt.subplots()
line, = ax.plot(time2, rpm2)
ax.set_xlabel('Time / min')
ax.set_ylabel('Angular velocity / rpm')

fig_html = mpld3.fig_to_html(fig)
components.html(fig_html, height=600)

st.markdown("""
            ### References
            *Solvent density and viscosity at varying temperatures are calculated 
            from linear and Arrhenius fits respectively from large data sets 
            downloaded from Reaxys. Full data sets and constants found on GitHub.*\n
            *Material data taken from cif files available on ICSD; thermal expansion assumed negligible.*
            """)