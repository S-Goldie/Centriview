# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 15:33:08 2024

@author: Stuart Goldie
"""

import streamlit as st

st.title('Centriview')

st.markdown(r"""
            ## Understanding and modelling the centrifugation of 2D nanosheets

            *Stuart Goldie,$^{1}$ Steffen Ott,$^{2}$ Anthony Dawson,$^{3}$ Tamara Starke,$^{2}$ 
            Cian Gabbett,$^{3}$ Victor Vega Mayoral,$^{4}$ Kevin Synnatschke,$^{2,3,5}$ Marilia Horn,$^{1,6}$ 
            Jonathan N. Coleman,$^{3,*}$ Claudia Backes,$^{1,2,*}$*

            Centriview is an interactive tool that uses equations of motion
            for 2D nanosheets under centrifugation to predict the outcome of experiments.
            This tool can predict the population function in terms of flake size 
            for a given experiment, and identify suitable centrifuge conditions
            to replicate different experiments.
            """)

st.image('https://github.com/S-Goldie/Centriview/blob/main/Centriview/flakes_falling.png?raw=true', width=360)

st.markdown("""
            ### How to use
            Use the sidebar to navigate to the different tools available.

            1. Centrifuge Speed Calculation - calculate the central rpm of a 2 step cascade to optimise for a target flake size
            2. Experiment comparison - identify the centrifuge conditions required to replicate seperation achieved with a different material or solvent system.
            3. Changing Rotor Dimensions - as above, with additional flexibility to change the rotor dimensions. _Note this can only be an approximate solution_
            4. Theoretical Discussion - a brief overview of the equations and theory behind the tool. For a full discussion, see the publication linked below.
            5. Fraction Population Prediction - visualise the change in flake size distribution for a given experiment.
            """)

st.markdown(r"""
            ### Want to know more?
            - The equations and full discussion are published at [arXve](https://arxiv.org/abs/2503.05111). If
            you make use of this tool, please reference the above publication in
            your own works.
            - The code and database for this webapp is made freely available and open-source 
            on [GitHub](https://github.com/S-Goldie/Centriview).

            ---
            ### Author Affiliations
            _\*Corresponding authors: [Jonathan N. Coleman](https://orcid.org/0000-0001-9659-9721), [Claudia Backes](https://orcid.org/0000-0002-4154-0439)_
            1)	Physical Chemistry of Nanomaterials and CINSaT, Kassel University, Heinrich-Plett Str. 40, 34132 Kassel, Germany
            2)	Applied Physical Chemistry, Heidelberg University, Im Neuenheimer Feld 253, 69120 Heidelberg, Germany
            3)	School of Physics and CRANN, Trinity College, Dublin 2, Ireland
            4)	Instituto Madrileño de Estudios Avanzados en Nanociencia (IMDEA), C/ Faraday 9, 28049 Madrid, Spain
            5)	Chair for Molecular Functional Materials, Dresden University of Technology, Stadtgutstr. 59, 01217 Dresden, Germany
            6)	University of Münster, Corrensstr. 3, 48149 Münster, Germany

            
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
            supported by funding from the European Union (2D-PRINTABLE HE: 101135196). 
            Views and opinions expressed are however those of the author(s) only and do 
            not necessarily reflect those of the European Union or the European Commission. 
            Neither the European Union nor the granting authority can be held responsible for them._
            """)

st.image('https://2d-printable.eu/storage/sites/18/2023/02/Funding_statement-1-768x161.png', caption=None, width=254)
