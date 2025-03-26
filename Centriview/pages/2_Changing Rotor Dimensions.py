# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 16:30:53 2024

@author: Stuart Goldie
"""

import streamlit as st
import streamlit.components.v1 as components
#import matplotlib.pyplot as plt
#from matplotlib import cm
#import matplotlib.style as mplstyle
import numpy as np
from scipy.optimize import minimize
#import mpld3
import plotly.graph_objects as go

###BACK END###

#Import Backend Data from Properties File#
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

#parameters1 = [mat1_pNS, mat1_d1, density1, viscosity1, r1, r2]

def supernatant_fraction_arrays(dummy_layer_numbers, dummy_lateral_size, time_min, rpm, params):
    "Produces an array of fraction remaining in the supernatant from arrays of thickness and lateral size"
    pNS, d1, pL, n, r1, r2 = params[0], params[1], params[2], params[3], params[4], params[5]
    time_s = 60 * time_min
    w = (rpm * 2 * np.pi) / 60
    experiment_constant = time_s*w*w/(12*n*np.cbrt(3/(4*np.pi)))
    surfactant_density_term = 2*d0*(pS-pL)
    nanosheet_density_term = pNS-pL
    x_array = [experiment_constant*(surfactant_density_term + nanosheet_density_term * x * d1) for x in dummy_layer_numbers]
    velocity_matrix = np.outer(dummy_lateral_size*1e-9/np.sqrt(2), x_array)
    fraction_matrix = np.maximum((r2*np.exp(-velocity_matrix)-r1)/(r2-r1), 0)
    return(fraction_matrix)
    
#Currently no weighting on this function
def difference_optimisation(target_rpm):
    "Compares the difference between two different sedimentation processes for optimisation"
    fraction_matrix1 = supernatant_fraction_arrays(dummy_layer_numbers, dummy_lateral_size, time_min, previous_rpm, parameters1)
    fraction_matrix2 = supernatant_fraction_arrays(dummy_layer_numbers, dummy_lateral_size, time_min, target_rpm, parameters2)
    weight = 0.5 * (fraction_matrix1+fraction_matrix2)
    #need to try and include the weighting function
    return(np.sum(weight*abs(fraction_matrix1-fraction_matrix2)))

###FRONT END###

st.title('Consistent Separation with Different Rotors')
st.markdown('## Beta - This page is still in Development')

st.markdown('This calculation tool is designed to compare the separation expected during a centrifugation process, \
        and find the most suitable conditions to replicate that using different parameters including different rotors and filling heights. \
        This final point is significant because an exact experimental equivalence can only be defined when the ratio of $R_1$ to $R_2$ \
        are equal between experiments. When different rotor geometery and fill heights are used only an approximate solution can be found.')

st.subheader('Experiment Calculator')
st.write("""
        To match experiments, select the conditions from the previous and new experiment. Changing too many parameters simultaneously 
        will cause greater error and instability in the suggested solution.

         _Note: because of the nature of the approximation, the solver will attempt to minimise the difference between the two experiments. 
         This may require longer than other calculations within this WebApp._
        """)


col1, col2 = st.columns(2)
col1.write("Original Experiment to be Replicated")
col2.write("New Experimental Conditions")

material1 = col1.selectbox('Material:',materials_list, key='material1selection')
material2 = col2.selectbox('Material:',materials_list, key='material2selection')

if material1 == 'Other':
    mat1_pNS = col1.number_input('Density of Custom Material in $g cm^{-3}$', key='mat1_pNS')*1000
    mat1_d1 = col1.number_input('Layer Thickness of Custom Material in $nm$', key='mat1_d1')*1E-9
else:
    mat1_pNS = prop.MATERIALS[material1][0]
    mat1_d1 = prop.MATERIALS[material1][1]

if material2 == 'Other':
    mat2_pNS = col2.number_input('Density of Custom Material in $g cm^{-3}$', key='mat2_pNS')*1000
    mat2_d1 = col2.number_input('Layer Thickness of Custom Material in $nm$', key='mat2_d1')*1E-9
else:
    mat2_pNS = prop.MATERIALS[material2][0]
    mat2_d1 = prop.MATERIALS[material2][1]

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
    density1 = prop.mixed_property_function(composition1, temp1, 'Density')*1000
    viscosity1 = prop.mixed_property_function(composition1, temp1, 'Viscosity')/1000
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
    density2 = prop.mixed_property_function(composition2, temp2, 'Density')*1000
    viscosity2 = prop.mixed_property_function(composition2, temp2, 'Viscosity')/1000
else:
    density2 = prop.solvent_density(temp2, solvent2)
    viscosity2 = prop.solvent_viscosity(temp2, solvent2)

r1 = col1.number_input('$R_1$ - Radius from axis to top of the liquid in mm', value=71)/1000
r2 = col1.number_input('$R_2$ - Radius from axis to top of the sediment in mm', value=100)/1000

q1 = col2.number_input('$Q_1$ - Radius from axis to top of the liquid in mm', value=71)/1000
q2 = col2.number_input('$Q_2$ - Radius from axis to top of the sediment in mm', value=100)/1000

parameters1 = [mat1_pNS, mat1_d1, density1, viscosity1, r1, r2]
parameters2 = [mat2_pNS, mat2_d1, density2, viscosity2, q1, q2]
  
time_min = col1.slider('Previous experiment time in minutes', min_value=10, max_value=540, value=120)
col1.markdown(f'Hours: {time_min/60:.1f}')
previous_rpm = col1.number_input('Previous rpm', value=1000)

#Define the plot parameters and optimise the difference to find the optimum rotation speed
dummy_layer_numbers = np.arange(start=1, stop=21, step=0.1)
dummy_lateral_size = np.arange(start=100, stop=1500, step=10)

X, Y = np.meshgrid(dummy_layer_numbers, dummy_lateral_size)

res = minimize(difference_optimisation, previous_rpm, method='nelder-mead')
optimal_rpm = res.x[0]

st.markdown(f"__Optimal rpm to match these experimental conditions (including time): {int(optimal_rpm)} rpm__")
st.write("The plots below show the difference between the two experimental conditions for comparison, and a plot of time/speed to allow consideration of the desired balance between achievable speed and time.")
    
fraction_matrix1 = supernatant_fraction_arrays(dummy_layer_numbers, dummy_lateral_size, time_min, previous_rpm, parameters1)
fraction_matrix2 = supernatant_fraction_arrays(dummy_layer_numbers, dummy_lateral_size, time_min, optimal_rpm, parameters2)


max_rpm = st.number_input('Maximum rpm:', value = previous_rpm)

time2 = np.arange(10, 540, 1)
rpm2 = np.sqrt(optimal_rpm**2 * time_min / time2)
optimal_time = optimal_rpm**2 * time_min / max_rpm**2

st.markdown(f"__Optimal time at maximum rpm to match the original experiment as entered above: {int(optimal_time)} min__")

st.caption('''
    The relationship between rpm and time is such that both can be changed together while still producing the same experimental result.
    In this case, the user can manipulate the plot to find the desired balance between rotation time and speed. Given a maximum rpm above, 
    the required experimental time is also shown.
    ''')

fig2 = go.Figure()

fig2.add_trace(go.Scatter(x=time2, y=rpm2, mode='lines', name='RPM vs Time', line=dict(color='black')))
fig2.add_trace(go.Scatter(x=[0, 540], y=[max_rpm, max_rpm], mode='lines', name='Max RPM', line=dict(color='grey', dash='dash')))
fig2.add_trace(go.Scatter(x=[optimal_time, optimal_time], y=[max_rpm * 0.8, max_rpm * 1.2], mode='lines', name='Optimal Time', line=dict(color='grey', dash='dash')))

fig2.update_layout(
    xaxis_title="Time / min",
    yaxis_title="Angular velocity / rpm",
    xaxis=dict(range=[0, 540]),
    yaxis=dict(range=[0, max_rpm * 1.2]),
    showlegend=True,
    template="plotly_white"
)

st.plotly_chart(fig2)


with st.expander("Click for more information on the experiment matching"):
    st.markdown("""
        ### Overview of how the experiment matching works
         To match different centrifuge experiments with so many changes, an analytical solution is not possible. 
        Instead, we model the change in population of different nanosheet sizes caused by the original centrifuge experiment. 
        The same calculation is then performed to identify the population change of different sized nanosheets following the 
        new centrifuge conditions proposed. This population change is mathematically defined by the function introduced on the Theoretical Discussion page.
        """)
    st.page_link("pages/3_Theoretical Discussion.py", label="Theoretical Discussion", icon="📃")
    st.markdown(r"""
        By iterating over different experimental conditions, a solution can be found that minimises the difference between these two 
        functions over a wide range of relevant nanosheet sizes. In other words, this solution describes the experimental condition that 
        should produce the same thickness and flake size as the original centrifuge process, accounting for differences in material 
        and solvent properties, and centrifuge hardware available.
                
        However, in most cases it is not possible to exactly match the nanosheet movement between 
        experiments for all nanosheet sizes. This is especially notable for materials with very 
        different layer thicknesses, because we are trying to optimise for separatoin according 
        to layer number and flake area.
                
        To visualise the quality of the solution, the percentage difference between the two functions 
        is plotted below for all nanosheet sizes. This heatmap identifies combinations of flake size 
        and thickness that have greater or lesser errors between the experiments.
        """)
    
    difference_matrix = np.abs(fraction_matrix1 - fraction_matrix2)
    
    fig1 = go.Figure(data=go.Heatmap(
        z=difference_matrix * 100,
        x=dummy_layer_numbers,
        y=dummy_lateral_size,
        colorscale='RdPu',
        colorbar=dict(title='% Error')
    ))

    fig1.update_layout(
        xaxis_title="Layer Number",
        yaxis_title="Lateral Size / nm",
        xaxis=dict(range=[1, 15]),
        yaxis=dict(range=[100, 1500]),
        template="plotly_white"
    )
    fig1.data[0].update(zmin=0.1, zmax=10)

    st.plotly_chart(fig1)
    
    st.caption('''
        A heat map of the percentage difference between the two experiments. Regions of greater 
        error are shown in darker colour.
        ''')