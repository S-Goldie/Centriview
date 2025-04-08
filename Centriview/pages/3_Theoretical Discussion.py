# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 16:30:53 2024

@author: Stuart Goldie
"""

import streamlit as st
import numpy as np
from streamlit import columns
import plotly.graph_objects as go

##Back End##
#Import constants and define functions to make the example 3D plot
layer_numbers = np.arange(0.5,20,0.04)
areas = np.arange(5000, 500000, 1000)
X, Y = np.meshgrid(layer_numbers, areas)

#Current example parameters for WS2 in NMP centrifugation
pL = 1028           #liquid density, assuming linear interpolation of light and heavy water at 36 oC
pS = 1595           #surfactant layer density, from Hirsham using anhydrous cholate on surface
d0 = 4.25E-10       #surfactant layer thickness, from Hirsham
n = 2.3E-3        #viscosity, approximate value at temperature and mixture [J. Kestin, Viscosity of light and heavy water and their mixtures,Physica A: Statistical Mechanics and its Applications,Volume 134, Issue 1,1985]
d1 = 6.18E-10       #monolayer of WS2, J. A. Wilson and A. D. Yoffe, Advances in Physics 18, 193 (1969).
pNS = 7500          #structure density

def supernatant_frac(areas, layer_number, time_hour, rpm):
    "Return the fraction of a flake population remaning in the supernatant"
    r1, r2 = 0.071, 0.1
    time_s = 3600 * time_hour
    h = layer_number * d1
    w = (rpm * 2 * np.pi) / 60
    a = time_s * w * w * (h * (pNS - pL) + 2 * d0 * (pS - pL )) * np.sqrt(areas * 1e-18) / (12 * n * np.cbrt(3 / (4 * np.pi)))
    return max((r2*np.exp(-a)-r1)/(r2-r1), 0)

vector_supernat = np.vectorize(supernatant_frac)

f1 = vector_supernat(Y,X,2,1800)
f2 = vector_supernat(Y,X,2,3600)

Z = (1-f2)*f1

## Front End ##
#Display text and discussion

st.title('Theoretical Discussion')

st.write("""
         The complete theoretical considerations are published in a pre-print available [here](https://arxiv.org/abs/2503.05111).
         Here is given a basic overview to faciliate an understanding of the calculations being run.

         The speed at which a nanosheet sinks in a centrifugal field can be defined from its size. 
         The relative movement of each size of nanosheet can be then compared to the tube dimensions. For a uniform starting distribution, 
         any nanosheets that move further than the length of the tube are considered sedimented. Flakes that sink less than this distance will be
         reduced in population according to a ratio of their movement to the total centrifuge tube length. 
         
         For example, a nanosheet that sinks half the length of the centrifuge tube will have its population halved. 
         This follows since half the population of nanosheets will start lower than halfway 
         down the tube, and therefore terminate in the sediment at the bottom by the end of the experiment. Whilst the other half, 
         starting further up the tube, will sink, but will not sink all the way to the sediment by the end of the experiment, 
         thereby remaining in the supernatant. This idea is used to produce a function defining the relative population of each 
         nanoflake remaining in the supernatant.
         """)

st.latex(r'''\begin{aligned} F_{\text {liq }}(A, h, \omega, t)=\left\{\begin{array}{rr}
\frac{R_2 e^{-S(h, A) \omega^2 t}-R_1}{R_2-R_1}, & S(h, A) \omega^2 t \leq \ln \left(\frac{R_2}{R_1}\right) \\
0, & S(h, A) \omega^2 t \geq \ln \left(\frac{R_2}{R_1}\right) \end{array}\right. \\ & \end{aligned}''')
st.write('where')
st.latex(r'''S(h, A)=\frac{1}{12 \eta \sqrt[3]{\frac{3}{4 \pi}}} \sqrt[2]{A}(h(\rho_{N S}-\rho_l)+2d(\rho_S-\rho_l))''')


st.markdown(
'''Defining all parameters, rotor geometry illustrated below:
- $R_1 $ = axial distance from center of rotation to top of the liquid surface $/ m$
- $R_2 $ = axial distance from center of rotation to bottom of the centrifuge tube $/ m$
- $t $ = time $/ s$
- $\omega $ = angular velocity $/ rad \space s^{-1}$
- $N $ = number of layers
- $d_1 $ = monolayer thickness $m$
- $ \Delta \\rho_{NS} $ = difference in density between nanosheet and liquid $/ kg \space m^{-3}$
- $d_0 $ = surfactant layer thickness $/ m$
- $\Delta \\rho_S $ = difference in density between surfactant layer and liquid $/ kg \space m^{-3}$
- $L $ = nanosheet length $/ m$
- $w $ = nanosheet width $/ m$
- $\eta $ = liquid viscosity $/ kg \space m^{-1} \space s^{-1}$)
'''
    )

st.image('https://github.com/S-Goldie/Centriview/blob/main/Centriview/radius_scheme.png?raw=true',caption='Schematic of radii depending on filling height and rotor geometry.')

#Theoretical Discussion underpinning the experiment comparison
st.subheader('Experiment Comparison Background')

st.write("""Building on the expression introduced above, we can consider the 
         effective separation of nanosheets in two different experiments. Equating $F_{liq}$
         from the two experiments, conditions can be identified that produce 
         consistent nanosheet size distributions from both.""")

st.markdown('__Consistent Rotor Geometry__')
st.markdown(r'Assuming a consistent geometry between experiments, that is the rotor and fill height remain unchanged, \
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

st.latex(r'''\frac{t_1 \omega_1^2 \eta_2}{t_2 \omega_2^2 \eta_1}=\frac{\rho_{NS_2}-\rho_{L_2}}{\rho_{NS_1}-\rho_{L_1}}''')

st.markdown('''where:
- $t$ = time
- $\omega$ = angular velocity (rpm)
- $\eta$ = viscosity
- $\\rho_{NS}$ = density of nanosheets
- $\\rho_l$ = density of solvent
            
To compare time (${t_2}$) and angular velocity ($\omega_2$), we plot the rearranged form of this expression:
''')

st.latex(r'''\boldsymbol{\omega_2} = \sqrt{\frac{t_1 \omega_1^2 \eta_2}{\boldsymbol{t_2} \eta_1} \frac{\rho_{NS_1}-\rho_{L_1}}{\rho_{NS_2}-\rho_{L_2}}}''')

#Theoretical Discussion explaining fraction populations
st.subheader('Sediment Fraction Populations')
st.markdown("""
            The expression above is a relative population function of nanosheets 
            remaining in the supernatant liquid after a given centrifuge experiment.
            That is for every nanosheet size, $F_L(\omega)$ returns a number 
            between 0 and 1 that is the fraction of the initial population of 
            that size of nanosheets remaining. Since all flakes are assumed retained
            across the two fractions, sediment and supernatant, the sediment 
            function is therefore by definition:
            """)
st.latex(r'''F_{sed}=1-F_{liq}''')
st.markdown("""
            If multiple centrifuge steps are combined, greater control and separation 
            can be achieved. Most importantly by removing both very large and very 
            small nanosheets from a sample. The effect of subsequent steps can be 
            described as the relative fraction after each experiment.
            """) 
st.image('https://github.com/S-Goldie/Centriview/blob/main/Centriview/cascade_scheme.png?raw=true', caption='Scheme of a three-step centrifuge process \
         highlighting the different fractions that are extracted by separation \
         of sediment and supernatant after each centrifuge step.')
st.markdown("""
            In the Fraction Population Prediction tool it is assumed that the 
            lower speed is used to remove the sediment, leaving only the supernatant 
            ($F_L(\omega_1)$) which is then centrifuged at higher speed ($\omega_2$) 
            and the sediment from this step is retained. In this way, as illustrated, the resulting relative population 
            function is the product ($F_L(\omega_1)F_s(\omega_2)$) of these sequential steps.
            """)

st.subheader('Example 3D Plot')
st.markdown("""
            The following plots are examples of the 3D surface described by the equations above, linking the flake area 
            and layer number to the fraction retained following a centrifuge cascade. If the initial, starting distribution 
            is known, this surface should be applied to the starting distribution to accurately describe the resulting flake 
            size distribution after the centrifuge processing.

            Fortunately, for many experimentally prepared dispersions, the starting distribution is so wide that the change 
            caused by the centrifuge is a dominant factor in the final distribution.
            """)

col1, col2 = st.columns(2)

# 3D Surface Plot in the first column
with col1:
    fig3d = go.Figure(data=[go.Surface(
        z=Z*100,
        x=layer_numbers,
        y=areas * 1e-6,
        colorscale='Viridis',
        showscale=False
    )])
    fig3d.update_layout(
        scene=dict(
            xaxis_title="Layer Numbers",
            yaxis_title="Area / sq. microns",
            zaxis_title="% Remaining",
            zaxis=dict(range=[0, 100]),
        ),
        title="3D Surface Plot"
    )
    st.plotly_chart(fig3d)
    st.caption("""Example 3D surface plot showing the fraction of nanosheets remaining in the supernatant after two centrifuge steps.
            This example data is modelling WS$_2$, centrifuged in NMP for two hours at 1800 rpm from which the supernatant is retained, and then at 3600 rpm from which the sediment is retained.""")

# 2D Contour Plot in the second column
with col2:
    fig = go.Figure(data=go.Contour(
        z=Z,
        x=layer_numbers,
        y=areas * 1e-6,
        colorscale='Viridis',
        contours=dict(showlabels=True),
        showscale=True,
        colorbar=dict(title="% Remaining")
    ))
    fig.update_layout(
        xaxis_title="Layer Numbers",
        yaxis_title="Area / sq. microns",
        title="2D Contour Plot"
    )
    st.plotly_chart(fig)
    st.caption("""Example 2D contour plot showing the same data as the 3D surface plot. Due to simplicity of plotting and reading, the 2D contour plot is preferred for data visualisation.""")