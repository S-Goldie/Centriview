# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 14:24:38 2024

@author: Goldie
"""

import streamlit as st
import streamlit.components.v1 as components
import matplotlib.pyplot as plt
import numpy as np
import mpld3

###BACK END###
#database import solid material constants, assuming negligable changes with temperature
#database import rotor geometry (calculations to resolve angle and fill height required?)
#solvent parameter calculation - select solvent of interest to define the relationship/model to use.
#For pure solvents use accurate temperature dependence, for solvent mixtures use less accurate 3d plotted model

#define functions for calculation of experiment coorespondance assuming constant or variable rotor geometry
#constant geometry is fairly straightforward and just needs material parameters
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
st.markdown('*citation does here*')
st.write('This calculation tool is designed to compare the separation expected during a centrifugation process, \
        and predict the experimental conditions required to replicate such separation when using different solvents, \
        materials or even different centrifuge rotors and geometries. Currently a working progress, these calculations \
        are underpinned by the theoretical work undertaken in the publication listed above.')

st.subheader('Experiment Calculator')
st.write('To match experiment conditions, select the conditions from the previous, old experiment to be matched and \
        the conditions intended for the new experiment.\
        At time of writing, many experimental parameters remain a working progress and may result in errors')


col1, col2 = st.columns(2)
col1.write("Original Experiment to be Replicated")
col2.write("New Experimental Conditions")

material1 = col1.selectbox('Material:',['WS2','MoS2','Graphene','NiPS3','MnPS3','GaS','BN','Ni(OH)2','Cu(OH)2','RuCl3','CrTe3'], key='material1selection')
material2 = col2.selectbox('Material:',['WS2','MoS2','Graphene','NiPS3','MnPS3','GaS','BN','Ni(OH)2','Cu(OH)2','RuCl3','CrTe3'], key='material2selection')
temp1 = col1.slider('Temperature (use left/right arrow keys for fine adjustment', min_value=5, max_value=50, key='temperature1')
temp2 = col2.slider('Temperature (use left/right arrow keys for fine adjustment', min_value=5, max_value=50, key='temperature2')
solvent1 = col1.selectbox('Solvent', ['Water', 'NMP', 'IPA', 'IPA/Water Mixture'], key='oldsolvent')
solvent2 = col2.selectbox('Solvent', ['Water', 'NMP', 'IPA', 'IPA/Water Mixture'], key='newsolvent')
if 'Mixture' in solvent1:
    composition1 = col1.slider('Weight Content of First Phase', min_value=0.0, max_value=1.0, value=0.5, key='composition1')
if 'Mixture' in solvent2:
    composition2 = col2.slider('Weight Content of First Phase', min_value=0.0, max_value=1.0, value=0.5, key='composition2')
    
time1 = col1.slider('Previous experiment time in minutes', min_value=10, max_value=360, value=120)
rpm1 = col1.number_input('Previous rpm')

st.markdown('*under construction*')

#Theoretical Discussion underpinning this calculation
st.subheader('Theoretical Background')

st.write('The fundamental derivation and proof of the equations that describe the motion and separation of particles \
            within the centrifuge are detailed in our paper [citation] and will not be repeated here. The key expression \
            on which we build is:')

st.latex(r'''\begin{aligned} F_{\text {liq }}(A, h, \omega, t)=\left\{\begin{array}{rr}
\frac{R_2 e^{-S(h, A) \omega^2 t}-R_1}{R_2-R_1}, & S(h, A) \omega^2 t \leq \ln \left(\frac{R_2}{R_1}\right) \\
0, & S(h, A) \omega^2 t \geq \ln \left(\frac{R_2}{R_1}\right) \end{array}\right. \\ & \end{aligned}''')
st.write('where')
st.latex(r'''S(h, A)=\frac{1}{12 \eta \sqrt[3]{\frac{3}{4 \pi}}} \sqrt[2]{A}(h(\rho_{N S}-\rho_l)+2d(\rho_S-\rho_l))''')

st.markdown('__Consistent Rotor Geometry__')
st.markdown('Assuming a consistent geometry between experiments, that is the rotor and fill height remain unchanged, \
        then $R_1$ and $R_2$ become identical across experiments and only the sedimentation, $S \omega^2 t$, matters. For \
        identical separation this must also be equal. Using subscripts to denote the separation to be replicated (1) and the new\
        experiment (2) and cancelling all identical parameters ($\sqrt{A}$ and $12 \sqrt[3]{3/{4 \pi}}$) this becomes:')

st.latex(r'''\frac{t_1 \omega_1^2}{\eta_1} (h(\rho_{NS_1}-\rho_{l_1})+2d(\rho_S-\rho_{l_1})) = 
         \frac{t_2 \omega_2^2}{\eta_2} (h(\rho_{NS_2}-\rho_{l_2})+2d(\rho_S-\rho_{l_2}))''')

st.markdown('$d$ and $\\rho_s$ refer to surfactant layer thickness and hydration density. We have shown this is only significant \
            for graphene in aqueous surfactant solutions. Materials with thicker monolayers and greater density are relatively \
            unaffected by surfactant coatings so these terms will be neglected for simplicity.')

st.latex(r'''\frac{t_1 \omega_1^2}{\eta_1} h(\rho_{NS_1}-\rho_{l_1}) = 
         \frac{t_2 \omega_2^2}{\eta_2} h(\rho_{NS_2}-\rho_{l_2})''')

st.markdown('Dividing through by $h$ (we are looking for identical sedimentation behaviour with respect to flake size between experiments)\
            and rearranging produces an expression for centrifuge conditions that will produce identical sedimentation.')

st.latex(r'''\frac{t_1 \omega_1^2 \eta_2}{t_2 \omega_2^2 \eta_1}=\sqrt{\frac{\rho_{NS_2}-\rho_{L_2}}{\rho_{NS_1}-\rho_{L_1}}}''')

st.markdown('''where:
- $t$ = time
- $\omega$ = angular velocity (rpm)
- $\eta$ = viscosity
- $\\rho_{NS}$ = density of nanosheets
- $\\rho_l$ = density of solvent
            ''')

st.markdown('__Variable Rotor Geometry__')

st.write()