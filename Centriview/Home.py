# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 15:33:08 2024

@author: Stuart Goldie
"""

import streamlit as st

st.title('Centriview Home')

st.markdown("""
            ## 2D-PRINTABLE – GA No: 101135196 – CONFIDENTIAL
            *Manuscript in preparation*
            
            Centriview is an interactive tool that uses equations of motion 
            for 2D nanosheets under centrifugation to predict the outcome of experiments.
            This tool can predict the population function in terms of flake size 
            for a given experiment, and identify suitable centrifuge conditions
            to replicate previous preparations with new materials.
            """)
st.image('flakes_falling.png')
st.markdown("""
            ### Want to know more?
            - The equations and full discussion are published at [doi here]. If
            you make use of this tool, please reference the above publication in
            your own works.
            - The code and database for this webapp is made freely available as open-source 
            on GitHub: [address here when public]
            """)


#FUTURE LIABILITY DISCLAIMER
#            *Disclaimer: THIS TOOL IS MADE FREELY AVAILABLE UNDER A CREATIVE 
#            COMMONS CC-BY LICENCE.  THE LICENSOR OFFERS THE LICENSED MATERIAL 
#            AS-IS AND AS-AVAILABLE, AND MAKES NO REPRESENTATIONS OR WARRANTIES 
#            OF ANY KIND CONCERNING THE LICENSED MATERIAL, WHETHER EXPRESS, IMPLIED, 
#            STATUTORY, OR OTHER. THIS INCLUDES, WITHOUT LIMITATION, WARRANTIES 
#            OF TITLE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, 
#            NON-INFRINGEMENT, ABSENCE OF LATENT OR OTHER DEFECTS, ACCURACY, OR 
#            THE PRESENCE OR ABSENCE OF ERRORS, WHETHER OR NOT KNOWN OR DISCOVERABLE.*