# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 00:41:12 2023

@author: Stuart Goldie
"""

import streamlit as st
import streamlit.components.v1 as components
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.style as mplstyle
import numpy as np
import mpld3


##_BACKEND DATA AND FUNCTIONS_##
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

solvent_list = ['Water'] + list(SOLVENTS.keys()) + ['IPA:Water Mixture','Other']
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

def speed_equivalence(t_1, w_1, n_1, t_2, n_2, p_NS1, p_NS2, p_l1, p_l2):
    """Returns the angular rotation to match experiments with constant geometry"""
    k = np.sqrt((p_NS1-p_l1)/(p_NS2-p_l2))
    w2 = k * (t_1*n_2*w_1**2)/(t_2*n_1)
    return np.sqrt(w2)

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

def mixed_property_function(w, T, c):
    return(c[0] + w*c[1] + T*c[2] + (w**2)*c[3] + w*T*c[4] + (T**2)*c[5] + (w**3)*c[6] + (w**2)*T*c[7] + w*(T**2)*c[8] + (T**3)*c[9])


pS = 1595           #surfactant layer density, from Hirsham using anhydrous cholate on surface
d0 = 4.25E-10       #surfactant layer thickness, from Hirsham


##_FRONT END_##
st.title('Modelling Centrifugation of 2D Nanomaterials')
st.subheader('CONFIDENTIAL')
st.markdown('_No information contained herein to be shared outside EU Project: 2D-PRINTABLE – GA No: 101135196_')
st.write("""
         This interactive tool allows you to visualise the changes to a 
         nanomaterial during centrifugation. By changing the experiment parameters 
         like rpm, centrifuge time and rotor dimensions the live plots will display 
         the relative population function according to flake size for a sample 
         collected between the two speeds, assuming all other conditions are kept 
         constant. For more information on this see the Theoretical Discussion page.\n
         Other materials are possible if the key constants are known. Note the 
         importance of entering values according to the correct units.
         """)
    
#Imput user parameters
material = st.selectbox('Pick a material:',materials_list)
if material == 'Other':
    pNS = st.number_input('Density of Custom Material in $g$ $cm^{-3}$', key='matdensity1', value=2.0)*1000
    d1 = st.number_input('Monolayer thickness from unit cell in Å', key='matthickness', value=3.0)*1E-10
else:
    pNS = MATERIALS[material][0]
    d1 = MATERIALS[material][1]
    k1 = MATERIALS[material][2]
    k2 = MATERIALS[material][3]

solvent = st.selectbox('Pick a solvent:',solvent_list)
temp = st.slider('Temperature (use left/right arrow keys for fine adjustment', min_value=5, max_value=50, key='temperature1')
if solvent == 'Water':
    pL = water_density(temp)
    n = water_viscosity(temp)
elif solvent == 'Other':
    pL = st.number_input('Density of Custom Solvent in $g$ $cm^{-3}$', key='soldensity1')*1000
    n = st.number_input('Viscosity of Custom Solvent in $cP$', key='solvisc1')/1000
elif solvent == 'IPA:Water Mixture':
    composition = st.slider('Weight Content of Alcohol', min_value=0.0, max_value=1.0, value=0.5, key='composition2')
    pL = mixed_property_function(composition, temp, IPAWaterDensity)*1000
    n = mixed_property_function(composition, temp, IPAWaterViscosity)/1000
else:
    pL = solvent_density(temp, solvent)
    n = solvent_viscosity(temp, solvent)

st.markdown('''*Illustrations of rotor geometry included in Theoretical Discussion*''')
r1 = st.number_input('$R_1$ - Radius from axis to top of the liquid in cm', value=7.1)/100
r2 = st.number_input('$R_2$ - Radius from axis to top of the sediment in cm', value=10)/100

time_hour = st.slider('Select time in minutes (use left/right arrow keys for fine adjustment)', min_value=10, max_value=540, value=120) / 60
st.markdown(f'Hours: {time_hour:.1f} *(adjusted in minutes)*')

st.markdown('''*Supernatant should be lower than sediment. See Discussion for scheme of intended cascade*''')
col1, col2 = st.columns(2)
rpm_lower = col1.slider('Select rpm of supernatant retention (use left/right arrow keys for fine adjustment)', min_value=100, max_value=30000)
rpm_higher = col2.slider('Select rpm of sedimentation retention (use left/right arrow keys for fine adjustment)', min_value=100, max_value=30000)

#Define the plot parameters and apply the model to the plot
dummy_layer_numbers = np.arange(start=0.5, stop=20.1, step=0.1)
dummy_lateral_size = np.arange(start=10, stop=1501, step=1)
if material != 'Other':
    aspect_ratio_lengths = 1e09*dummy_layer_numbers*d1*k2/np.sqrt(k1)

X, Y = np.meshgrid(dummy_layer_numbers, dummy_lateral_size)

f1 = supernatant_fraction_arrays(dummy_layer_numbers, dummy_lateral_size, time_hour, rpm_lower, r1, r2)
f2 = supernatant_fraction_arrays(dummy_layer_numbers, dummy_lateral_size, time_hour, rpm_higher, r1, r2)

Z = (1-f2)*f1

#cutoff_values = [cutoff_area(m, time_hour, rpm_higher, r1, r2) for m in dummy_layer_numbers]
    
mplstyle.use('fast')
fig1, ax1 = plt.subplots()
CS = ax1.contour(X, Y, Z, 6, cmap=cm.coolwarm)
ax1.set(xlim=(0, 20), ylim=(0, 1500),  xlabel="Layer Number", ylabel="Lateral Size / nm")
if material != 'Other':
    ax1.plot(dummy_layer_numbers, aspect_ratio_lengths, linestyle='--', linewidth=1, color='grey')
plt.tight_layout()
st.subheader('Flake Size Probability Distribution')
with st.expander("See more information"):
    st.markdown("""
         The change in nanosheet distribution following the modelled centrifuge 
         process is shown below as a 2D contour plot (showing thickness and lateral
         sizes independently) and a 1D line plot (where lateral size and thickness 
         are correlated according to previously published trends).$^{[1]}$\n
         This function should be considered the _change_ in nanosheet population
         following centrifugation. If a completely uniform dispersion with equal
         population of all possible flake sizes were centrifuged this function 
         would perfectly describe the resulting nanosheet size distribution. However, 
         since all experimentally prepared samples have an initial size distribution 
         this guide merely illustrates the change expected.\n
         The contour plot shows regions with no population remaining outside the 
         dark blue lines, with the population retained increasing towards orange at 
         the centre. The distribution along the grey dashed line (most common 
         aspect-ratio) is shown as a simple 1D plot underneath.
         """)

st.pyplot(fig1)
st.caption('''
           Colour contour plot showing the population distribution of flake sizes given 
           experimental details selected above. Lateral flake size refers to $\sqrt{LW}$. 
           The darker, orange represents a high population decaying down through 
           blue. The grey dashed line shows the mean aspect ratio of liquid phase 
           exfoliated nanosheets.''')

if material == 'Other':
    st.write('''__For a new material with unknown aspect-ratio's the 1D trend cannot be estimated.__''')
elif material != 'Other':
    with st.expander("See more information"):
        st.markdown("""
                    To better illustrate this 2D contour plot, consider the mean nanosheet 
                    aspect ratio shown along the grey dashed line. The population change 
                    may be uniform between large-and-thin sheets and small-and-thick; 
                    however experimentally these are very uncommon. Moving further from 
                    the grey line, it is unlikely that any such flakes would exist in the 
                    starting dispersion and therefore would remain absent from the final product.\n
                    Plotting the population change distribution along this grey line, see below, 
                    shows the expected change in nanosheet size distribution for the most 
                    common nanosheet aspect ratios.
                    """)
    
    f1_linear = fraction_linear(dummy_layer_numbers, time_hour, rpm_lower)
    f2_linear = fraction_linear(dummy_layer_numbers, time_hour, rpm_higher)
    f2_sed = (1-f2_linear)*f1_linear
    
    fig2, ax2 = plt.subplots()
    ax2.set(xlim=(0,20), xlabel="Layer Number", ylabel="% Population Remaining")
    ax2.plot(dummy_layer_numbers, f2_sed*100, color='grey', linestyle='--')
    st.pyplot(fig2)
    st.caption('''
               A line plot showing the population change function along the grey-dotted 
               line of the contour plot above. The aspect-ratios are taken from 
               literature results of known materials. $^{[1]}$
               ''')

st.markdown("""
            ### References
            [1] = https://pubs.acs.org/doi/full/10.1021/acsnano.9b02234 \n
            *Solvent density and viscosity at varying temperatures are calculated 
            from linear and Arrhenius fits respectively from large data sets 
            downloaded from Reaxys. Full data sets and constants found on GitHub.*\n
            *Material data taken from cif files available on ICSD; thermal expansion assumed negligible.*
            """)