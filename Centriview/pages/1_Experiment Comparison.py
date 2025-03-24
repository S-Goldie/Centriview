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

#Import Backend Data from Properties File#
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)

sys.path.append(parent_dir)
import properties as prop
solvent_list = ['Water'] + list(prop.SOLVENTS.keys()) + ['IPA:Water Mixture', 'Other']
materials_list = list(prop.MATERIALS.keys()) + ['Other']

def speed_equivalence(t_1, w_1, n_1, t_2, n_2, p_NS1, p_NS2, p_l1, p_l2):
    """Returns the angular rotation to match experiments with constant geometry"""
    k = (p_NS1-p_l1)/(p_NS2-p_l2)
    w2 = k * (t_1*n_2*w_1**2)/(t_2*n_1)
    return np.sqrt(w2)

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
    thickness_material1 = col1.number_input('Layer Thickness of Custom Material in $nm$', key='mat1_d')*1E-9
else:
    density_material1 = prop.MATERIALS[material1][0]
    thickness_material1 = prop.MATERIALS[material1][1]
if material2 == 'Other':
    density_material2 = col2.number_input('Density of Custom Material in $g cm^{-3}$', key='matdensity2')*1000
    thickness_material2 = col2.number_input('Layer Thickness of Custom Material in $nm$', key='mat2_d')*1E-9
else:
    density_material2 = prop.MATERIALS[material2][0]
    thickness_material2 = prop.MATERIALS[material2][1]

temp1 = col1.slider('Temperature (use left/right arrow keys for fine adjustment', min_value=5, max_value=50, key='temperature1')
temp2 = col2.slider('Temperature (use left/right arrow keys for fine adjustment', min_value=5, max_value=50, key='temperature2')
solvent1 = col1.selectbox('Solvent', solvent_list, key='oldsolvent')
solvent2 = col2.selectbox('Solvent', solvent_list, key='newsolvent')

if solvent1 == 'Water':
    density1 = prop.water_density(temp1)
    viscosity1 = prop.water_viscosity(temp1)
elif solvent1 == 'Other':
    density1 = col1.number_input('Density of Custom Solvent in $g cm^{-3}$', key='soldensity1')*1000
    viscosity1 = col1.number_input('Viscosity of Custom Solvent in $cP$', key='solvisc1')/1000
elif solvent1 == 'IPA:Water Mixture':
    composition1 = col1.slider('Weight Content of Alcohol', min_value=0.0, max_value=1.0, value=0.5, key='composition1')
    density1 = prop.mixed_property_function(composition1, temp1, IPAWaterDensity)*1000
    viscosity1 = prop.mixed_property_function(composition1, temp1, IPAWaterViscosity)/1000
else:
    density1 = prop.solvent_density(temp1, solvent1)
    viscosity1 = prop.solvent_viscosity(temp1, solvent1)

if solvent2 == 'Water':
    density2 = prop.water_density(temp2)
    viscosity2 = prop.water_viscosity(temp2)
elif solvent2 == 'Other':
    density2 = col2.number_input('Density of Custom Solvent in $g cm^{-3}$', key='soldensity2')*1000
    viscosity2 = col2.number_input('Viscosity of Custom Solvent in $cP$', key='solvisc2')/1000
elif solvent2 == 'IPA:Water Mixture':
    composition2 = col2.slider('Weight Content of Alcohol', min_value=0.0, max_value=1.0, value=0.5, key='composition2')
    density2 = prop.mixed_property_function(composition2, temp2, IPAWaterDensity)*1000
    viscosity2 = prop.mixed_property_function(composition2, temp2, IPAWaterViscosity)/1000
else:
    density2 = prop.solvent_density(temp2, solvent2)
    viscosity2 = prop.solvent_viscosity(temp2, solvent2)

  
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
