# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 16:30:53 2024

@author: Stuart Goldie
"""

import streamlit as st

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
