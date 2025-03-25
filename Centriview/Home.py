# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 15:33:08 2024

@author: Stuart Goldie
"""

import streamlit as st

st.title('Centriview')

st.markdown("""
            Centriview is an interactive tool that uses equations of motion 
            for 2D nanosheets under centrifugation to predict the outcome of experiments.
            This tool can predict the population function in terms of flake size 
            for a given experiment, and identify suitable centrifuge conditions
            to replicate different experiments.

            ### How to use
            Use the sidebar to navigate to the different tools available.

            1. Experiment comparison - this tool allows you to identify the centrifuge conditions required to replicate the flake size seperation achieved with a different material or solvent system.
            2. Changing Rotor Dimensions - as above, with additional flexibility to change the rotor dimensions. _Note this can only be an approximate solution_
            3. Theoretical Discussion - a brief overview of the equations and theory behind the tool. For a full discussion, see the publication linked below.
            4. Fraction Population Prediction - this tool allows you to visualise the change in flake size distribution for a given experiment.
            """)

st.image('https://github.com/S-Goldie/Centriview/blob/main/Centriview/flakes_falling.png?raw=true', width=360)
st.markdown("""
            ### Want to know more?
            - The equations and full discussion are published at [arXve](https://arxiv.org/abs/2503.05111). If
            you make use of this tool, please reference the above publication in
            your own works.
            - The code and database for this webapp is made freely available and open-source 
            on [GitHub](https://github.com/S-Goldie/Centriview).

            ### Disclaimer
            _This tool is made freely available under a creative 
            commons cc-by licence.  The licensor offers the licensed material 
            as-is and as-available, and makes no representations or warranties 
            of any kind concerning the licensed material, whether express, implied, 
            statutory, or other. This includes, without limitation, warranties 
            of title, merchantability, fitness for a particular purpose, 
            non-infringement, absence of latent or other defects, accuracy, or 
            the presence or absence of errors, whether or not known or discoverable._
            
            _The tools and data files made available within this repository were 
            supported by funding from the European Union (2D-PRINTABLE HE: 101135196)._
            """)

st.image('https://upload.wikimedia.org/wikipedia/commons/e/e0/Flag_of_Europe_color_with_border.svg', caption=None, width=254)
