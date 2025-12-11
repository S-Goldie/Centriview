# -*- coding: utf-8 -*-
"""
Created on Sat Jul 12 14:51:42 2025

@author: Stuart Goldie
"""

import streamlit as st
import streamlit.components.v1 as components
import numpy as np
import plotly.graph_objects as go

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

def central_speed(time, target_N, target_L, ratio, p_l, n, p_NS, h_NS, r1, r2):
    """Returns the central speed, in rads per second, required for a target flake size given experimental parameters for a range of times"""
    w2t = (12*n*np.sqrt(2*ratio)*np.cbrt(3/(4*np.pi))) * np.log(2*r2/(r2+r1)) / (target_L*target_N*h_NS*(p_NS-p_l))
    return np.sqrt(w2t/time)
    

##_FRONT END CODE_##
st.title('Centrifuge Speed Calculator')
st.markdown("""
This calculator allows you to estimate the central speed, around which a two-step centrifuge process can be planned to optimise the yield of a target flake size. 
         
To do a 2-step cascade, you need to use 2 rpms and extract the supernatant from the slower, and sediment from the higher. These speeds should be equidistant 
            above and below the central value calculated here. When deciding the separation of the upper and lower rpms, you need to consider that making them 
            closer will give a narrower size distribution but less mass. Making them further apart will do the opposite.

Centrifugation can only select flakes already present within a sample, and the true selection is a product of $N×L$. For a more complete visualisation of 
         the size separation see:
""")
st.page_link("pages/4_Fraction Population Prediction.py", label="Size Selection Visualisation", icon="📊")
st.markdown('For a more comprehensive discussion see the Theoretical Discussion page or the accompanying publication available through the links below.')
st.page_link("pages/5_Theoretical Discussion.py", label="Theoretical Discussion", icon="📃")
st.page_link("https://arxiv.org/abs/2503.05111", label="Publication", icon="📜")

st.subheader('Experiment Calculator')
st.write('To estimate the central rpm and time, enter the following details:')
st.write("""
         - **Material and Solvent:** Select from common options or enter custom values. Density and viscosity are automatically calculated for standard choices.
         - **Rotor Dimensions:** Input the radii from the centrifuge axis to the top of the liquid and sediment for accurate rpm estimation.
         - **Target Nanosheet Size:** Specify the desired layer number and lateral size. The calculator uses the aspect ratio from the selected material.
         - **Time and Speed Relationship:** The graph shows how centrifuge time and speed are linked; any point on the curve will yield the target flake size.
         
         *When entering all values check the required units carefully, this is a common error when comparing information from different sources.*""")


material1 = st.selectbox('Material:',materials_list, key='material1selection')

if material1 == 'Other':
    density_material = st.number_input('Density of Custom Material in $g cm^{-3}$', key='matdensity1')*1000
    thickness_material = st.number_input('Layer Thickness of Custom Material in $nm$', key='mat1_d')*1E-9
    ratio = st.slider('Estimated length/width aspect ratio', min_value=1.0, max_value=20.0, value=2.0)
    st.markdown('*If aspect ratios are unknown, typical values for many 2D materials are around 2*')
else:
    density_material = prop.MATERIALS[material1][0]
    thickness_material = prop.MATERIALS[material1][1]
    ratio = prop.MATERIALS[material1][2]


solvent1 = st.selectbox('Solvent', solvent_list, key='solvent')
temp1 = st.slider('Temperature / $^o C$', min_value=5, max_value=50, key='temperature1')


if solvent1 == 'Water':
    density1 = prop.water_density(temp1)
    viscosity1 = prop.water_viscosity(temp1)
elif solvent1 == 'Other':
    density1 = st.number_input('Density of Custom Solvent in $g cm^{-3}$', key='soldensity1')*1000
    viscosity1 = st.number_input('Viscosity of Custom Solvent in $cP$', key='solvisc1')/1000
elif solvent1 == 'IPA:Water Mixture':
    composition1 = st.slider('Weight Content of Alcohol', min_value=0.0, max_value=1.0, value=0.5, key='composition1')
    density1 = prop.mixed_property_function(composition1, temp1, 'Density')*1000
    viscosity1 = prop.mixed_property_function(composition1, temp1, 'Viscosity')/1000
else:
    density1 = prop.solvent_density(temp1, solvent1)
    viscosity1 = prop.solvent_viscosity(temp1, solvent1)

  
st.markdown('Information about the rotor dimensions are required to accurately identify the rpm. For more information on measuring these parameters, see the schematic on Theoretical Discussion.')
col1, col2 = st.columns(2)
r1 = col1.number_input('$R_1$ - Radius from axis to top of the liquid in mm', value=71.0)/1000
r2 = col2.number_input('$R_2$ - Radius from axis to top of the sediment in mm', value=100.0)/1000


st.markdown('Enter the desired nanosheet size. Remember the estimation will only work for nanosheets that are present in a significant population in the starting dispersion.')
col3, col4 = st.columns(2)
target_N = col3.number_input('Desired nanosheet layer number', key='target_thickness', step=1)
target_L = col4.number_input('Desired nanosheet length in nm', key='target_length')*1E-9
st.markdown('*Note: length is entered and lateral size is estimated using the length/width aspect ratio from the material selected.*')


time1 = st.slider('Ideal experiment time in minutes', min_value=10, max_value=540, value=120)
st.markdown(f'Hours: {time1/60:.1f}')

##_BACK END CALCULATIONS_##
if (target_L*target_N*thickness_material*(density_material-density1)) != 0:
    time2 = np.arange(10, 540, 1)
    w2 = central_speed(time2*60, target_N, target_L, ratio, density1, viscosity1, density_material, thickness_material, r1, r2)
    rpm2 = w2*60/(2*np.pi)

    speed_match_rads = central_speed(time1*60, target_N, target_L, ratio, density1, viscosity1, density_material, thickness_material, r1, r2)
    speed_match = speed_match_rads*60/(2*np.pi)
    st.markdown(f'__Speed needed to match previous time : {speed_match:.0f} rpm__')


    fig = go.Figure()
    fig.add_trace(go.Scatter(x=time2, y=rpm2, mode='lines', name='Angular Velocity'))
    fig.update_layout(
        xaxis_title='Time / min',
        yaxis_title='Angular velocity / rpm',
        title='Angular Velocity vs Time',
        height=600
    )

    st.plotly_chart(fig)

    st.caption('Since centrifuge time and angular velocity are linked, any combination of time and speed along the plotted line can be used.')

st.markdown("""
            ### References
            *Solvent density and viscosity at varying temperatures are calculated 
            from linear and Arrhenius fits respectively from large data sets 
            downloaded from Reaxys. Full data sets and constants found on [GitHub](https://github.com/S-Goldie/Centriview).*
            
            *Material data taken from cif files available on ICSD; thermal expansion assumed negligible.*
            """)

st.divider()
st.image('https://2d-printable.eu/storage/sites/18/2023/02/Funding_statement-1-768x161.png', caption=None, width=254)
