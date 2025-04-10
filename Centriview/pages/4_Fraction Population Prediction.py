# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 00:41:12 2023

@author: Stuart Goldie
"""

import streamlit as st
import streamlit.components.v1 as components
import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.style as mplstyle
import mpld3


##_BACKEND DATA_##
#Import backend data from properties file#
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)

sys.path.append(parent_dir)
import properties as prop

pS = 1595           #surfactant layer density, from Hirsham using anhydrous cholate on surface
d0 = 4.25E-10       #surfactant layer thickness, from Hirsham

solvent_list = ['Water'] + list(prop.SOLVENTS.keys()) + ['IPA:Water Mixture','Other']
materials_list = list(prop.MATERIALS.keys()) + ['Other']

solvent_list = ['Water'] + list(prop.SOLVENTS.keys()) + ['IPA:Water Mixture','Other']
materials_list = list(prop.MATERIALS.keys()) + ['Other']

##_Define key functions for calculations_##

def supernatant_fraction_arrays(layer_number, lateral_size, time_hour, rpm, r1, r2):
    "Produces an array of fraction remaining in the supernatant from arrays of thickness and lateral size"
    time_s = 3600 * time_hour
    w = (rpm * 2 * np.pi) / 60
    experiment_constant = time_s*w*w/(12*n*np.cbrt(3/(4*np.pi)))
    surfactant_density_term = 2*d0*(pS-pL)
    nanosheet_density_term = pNS-pL
    x_array = [experiment_constant*(surfactant_density_term + nanosheet_density_term * x * d1) for x in layer_number]
    velocity_matrix = np.outer(lateral_size*1e-9/np.sqrt(2), x_array)
    fraction_matrix = np.maximum((r2*np.exp(-velocity_matrix)-r1)/(r2-r1), 0)
    return(fraction_matrix)

def fraction_linear(layer_number, time_hour, rpm):
    time_s = 3600 * time_hour
    w = (rpm * 2 * np.pi) / 60
    experiment_constant = time_s*w*w*k2/(12*n*np.sqrt(2*k1)*np.cbrt(3/(4*np.pi)))
    surfactant_density_term = 2*d0*(pS-pL)
    nanosheet_density_term = pNS-pL
    height = d1 * layer_number
    array = -1*experiment_constant*(height**2 * nanosheet_density_term + height*surfactant_density_term)
    fraction = np.maximum((r2*np.exp(array)-r1)/(r2-r1), 0)
    return(fraction)

def cutoff_area(layer_number, time_hour, rpm, r1, r2):
    time_s = 3600 * time_hour
    h = layer_number * d1
    w = (rpm * 2 * np.pi) / 60
    return (((np.log(r2 / r1) * (12 * n * np.cbrt(3 / (4 * np.pi))) / (time_s * w * w * (h * (pNS - pL) + 2 * d0 * (pS - pL)))) ** 2) * 1e18)

##_FRONT END_##
st.title('Fraction Population Prediction')

st.write("""
         This interactive tool allows you to visualize the changes to the size distribution of a 2D
         nanomaterial during centrifugation. By changing the experiment parameters 
         like rpm, centrifuge time and rotor dimensions the live plots will display 
         the relative population function according to flake size for a sample 
         collected between the two speeds. For detailed information on this see the 
         Theoretical Discussion page or the accompanying publication available through the links below.
         """)
st.page_link("pages/3_Theoretical Discussion.py", label="Theoretical Discussion", icon="📃")
st.page_link("https://arxiv.org/abs/2503.05111", label="Publication", icon="📜")

st.subheader("Experimental Parameters")
    
#Input user parameters
material = st.selectbox('Pick a material:', materials_list)
if material == 'Other':
    pNS = st.number_input('Density of Custom Material in $g$ $cm^{-3}$', key='matdensity1', value=2.0) * 1000
    d1 = st.number_input('Monolayer thickness from unit cell in Å', key='matthickness', value=3.0) * 1E-10
    aspect_flag = False
elif len(prop.MATERIALS[material]) > 2:
    pNS = prop.MATERIALS[material][0]
    d1 = prop.MATERIALS[material][1]
    k1 = prop.MATERIALS[material][2]
    k2 = prop.MATERIALS[material][3]
    aspect_flag = True
else:
    pNS = prop.MATERIALS[material][0]
    d1 = prop.MATERIALS[material][1]
    aspect_flag = False

solvent = st.selectbox('Pick a solvent:',solvent_list)

st.markdown("""
         _Other materials and solvents can be entered manually if the key constants are known. Note the 
         importance of entering values according to the correct units._
         """)

temp = st.slider('Temperature (use left/right arrow keys for fine adjustment', min_value=5, max_value=50, key='temperature1')

if solvent == 'Water':
    pL = prop.water_density(temp)
    n = prop.water_viscosity(temp)
elif solvent == 'Other':
    pL = st.number_input('Density of Custom Solvent in $g$ $cm^{-3}$', key='soldensity1')*1000
    n = st.number_input('Viscosity of Custom Solvent in $cP$', key='solvisc1')/1000
elif solvent == 'IPA:Water Mixture':
    composition = st.slider('Weight Content of Alcohol', min_value=0.0, max_value=1.0, value=0.5, key='composition2')
    pL = prop.mixed_property_function(composition, temp, 'Density')*1000
    n = prop.mixed_property_function(composition, temp, 'Viscosity')/1000
else:
    pL = prop.solvent_density(temp, solvent)
    n = prop.solvent_viscosity(temp, solvent)

st.markdown("""Illustrations of rotor geometry are included in the Theoretical Discussion. 
            These values are typical for a full tube in a fixed angle rotor, but for best 
            results the exact geometry should be measured.""")
r1 = st.number_input('$R_1$ - Radius from axis to top of the liquid in cm', value=7.1)/100
r2 = st.number_input('$R_2$ - Radius from axis to top of the sediment in cm', value=10)/100

time_hour = st.slider('Select time in minutes (use left/right arrow keys for fine adjustment)', min_value=10, max_value=540, value=120) / 60
st.markdown(f'Hours: {time_hour:.1f} h')

st.markdown('''
            This is intended to be a two-step process. The first, lower speed step is to 
            remove the sediment and then the second, higher speed step is to sediment the 
            desired material and discard the remaining supernatant.
            ''')
col1, col2 = st.columns(2)
rpm_lower = col1.slider('Select rpm of supernatant retention', min_value=100, max_value=40000)
rpm_higher = col2.slider('Select rpm of sedimentation retention', min_value=100, max_value=40000)

#Define the plot parameters and apply the model to the plot
dummy_layer_numbers = np.arange(start=0.5, stop=20.1, step=0.1)
dummy_lateral_size = np.arange(start=10, stop=1501, step=1)
if aspect_flag == True:
    aspect_ratio_lengths = 1e09*dummy_layer_numbers*d1*k2/np.sqrt(k1)

X, Y = np.meshgrid(dummy_layer_numbers, dummy_lateral_size)

f1 = supernatant_fraction_arrays(dummy_layer_numbers, dummy_lateral_size, time_hour, rpm_lower, r1, r2)
f2 = supernatant_fraction_arrays(dummy_layer_numbers, dummy_lateral_size, time_hour, rpm_higher, r1, r2)

Z = (1-f2)*f1

#cutoff_values = [cutoff_area(m, time_hour, rpm_higher, r1, r2) for m in dummy_layer_numbers]
    
st.subheader('Flake Size Probability Distribution')
st.markdown('''
            Using the experimental parameters entered above, the plots below show the _change_ in population 
         distribution of different sized nanosheets. For more information on interpreting these plots see the
         Theoretical Discussion page or click `See more information` below.
            ''')

if aspect_flag == False:
    st.write('''__For a new material with unknown aspect-ratios the 1D trend cannot be estimated.__''')
elif aspect_flag == True:
    with st.expander("See more information on the population function and linear plot"):
        st.markdown("""
            The function calculated from the experimental parameters above describes the _change_ in 
            nanosheet size distribution after centrifugation. A high percentage remaining indicates 
            retention of nanosheets of that size, while a low percentage indicates they are likely 
            to be discarded. By changing the experimental parameters, it is possible to predict the 
            effect of different centrifugation conditions. In a rigorous treatment, this function of 
            centrifuge parameters and nanosheet sizes should be multiplied by the starting nanosheet size 
            distribution to predict the final distribution. However, since this is very rarely known, 
            we simply plot the change that a given centrifuge process will cause.
                
            To more easily visualize this change to the nanosheet size distribution, we can use the
            most common aspect ratios known for liquid-phase exfoliated nanosheets to reduce the 
            dimensionality into a 1D line plot.$^{[1]}$
            
            In this plot, the population remaining in the final sediment after the two-step process
            described above is shown as a function only of layer number, while the flake lateral 
            area is a fixed aspect ratio of thickness.
            """)
    
    f1_linear = fraction_linear(dummy_layer_numbers, time_hour, rpm_lower)
    f2_linear = fraction_linear(dummy_layer_numbers, time_hour, rpm_higher)
    f2_sed = (1-f2_linear)*f1_linear
    
    fig3 = go.Figure()
    fig3.add_trace(go.Scatter(x=dummy_layer_numbers, y=f2_sed*100, mode='lines', line=dict(color='grey', dash='dash')))
    fig3.update_layout(
        xaxis_title="Layer Number",
        yaxis_title="% Population Remaining",
        xaxis=dict(range=[0, 20]),
        template="plotly_white"
    )
    st.plotly_chart(fig3)

    st.caption('''
               A line plot showing the population change function using common aspect-ratio's
               known for liquid-phase exfoliated nanosheets to reduce the dimensionality.$^{[1]}$
               This fixed aspect ratio is plotted as a dashed grey line on the contour plot below. 
               ''')
    
with st.expander("See more information and to view 3D plot"):
    st.markdown(r"""
        For nanosheets prepared by other techniques, the flake area and thickness may not be related 
        by the same aspect ratios as those known for liquid-phase exfoliation samples. In this case, the full 
        3D surface describing the percentage of flakes remaining as a function of lateral size and thickness, 
        as separate variables, is a more accurate description.   
                    
        To more easily render such 3D surfaces as the conditions are changed, they are shown here as 2D contour plots.
        As described in the caption, to read these 2D plots, the regions beyond the dark blue lines correspond 
        to the dark purple regions of the 3D plot where the flake population is zero. No flakes of such size will 
        remain after centrifugation. Inside the two orange contour lines, the highest population will be retained.
        (If only one orange line is shown, the speeds selected are too low to remove the largest flakes.)
                
        To more easily understand this plot, it can be useful to refer to the line plot above. This follows 
        the 3D surface over a path described by fixed aspect ratios, shown as the grey dashed
        line through the contour plot. The peak, corresponding to the highest population remaining, should fall 
        within the boundaries of the orange lines, while the function decays to zero beyond the blue boundaries.
                
        Fully interactive example 3D surface and contour plots are also shown on the Theoretical Discussion page.
        """)
    
    if st.button("Generate 3D Surface Plot"):
        fig3d = go.Figure(data=[go.Surface(
            z=Z*100,
            x=X,
            y=Y,
            colorscale='Viridis',
            showscale=False
        )])
        fig3d.update_layout(
            scene=dict(
                xaxis_title="Layer Number",
                yaxis_title="Lateral Size / nm",
                zaxis_title="% Remaining",
            ),
            title="3D Surface Plot"
        )
        st.plotly_chart(fig3d)

        st.caption(r"""
                   3D surface plot showing the fraction of nanosheets remaining in the supernatant after the two 
                   centrifuge steps described above, as defined by both lateral size, $\sqrt{LW}$, and layer 
                   number $N$. The dark purple region indicates flakes so large they are completely removed by 
                   the first centrifuge process, while the decay at small flake sizes shows flakes so small they 
                   are not appreciably concentrated by even the higher speed.

                   """)

mplstyle.use('fast')
fig1, ax1 = plt.subplots()
CS = ax1.contour(X, Y, Z, 6, cmap=cm.coolwarm)
ax1.set(xlim=(0, 20), ylim=(0, 1500),  xlabel="Layer Number", ylabel="Lateral Size / nm")
if aspect_flag == True:
    ax1.plot(dummy_layer_numbers, aspect_ratio_lengths, linestyle='--', linewidth=1, color='grey')
plt.tight_layout()
st.pyplot(fig1)

st.caption(r'''
           Contour plot of the 3D surface that describes the nanosheet population remaining as defined by 
           both lateral size, $\sqrt{LW}$, and layer number $N$. Regions outside the dark blue lines represent 
           flakes so large they are completely removed by the first centrifuge process. The contour lines of 
           changing color represent the gradient of the surface, where changing flake size has a large 
           influence on the population remaining. The region containing the flakes that are more retained 
           by the selected centrifuge experiment is enclosed by the darkest orange lines. 
           This represents the combination of flake thickness and lateral size that is most concentrated.

           _A 3D surface plot is available within the information panel above, but this may reduce performance._ 
           ''')

st.markdown("""
            ### References
            [1] = https://pubs.acs.org/doi/full/10.1021/acsnano.9b02234 \n
            *Solvent density and viscosity at varying temperatures are calculated 
            from linear and Arrhenius fits respectively from large data sets 
            downloaded from Reaxys. Full data sets and constants found on GitHub.*\n
            *Material data taken from cif files available on ICSD; thermal expansion assumed negligible.*
            """)
