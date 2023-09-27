# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 00:41:12 2023

@author: Goldie
"""

import streamlit as st
import streamlit.components.v1 as components
import matplotlib.pyplot as plt
import mpld3

st.title('First Trial Page Title')
st.header('Trial of Heading')
st.subheader('Smaller Heading Test')
st.write('This is where we would start to explain how the web app should work and how the equations work. Exactly how much of the theory we try and explain is probably limited since ht is what the paper will be for')
st.markdown('_Markdown_') # see #*


st.caption('Caption text to go underneath the figures being plotted')
st.latex(r''' e^{i\pi} + 1 = 0 ''')

material = st.radio('Pick a material:',['WS2','MoS2','Gra'])
rpm = st.slider('Select rpm (use left/right arrow keys for fine adjustment)', min_value=0, max_value=60000)

fig = plt.figure()
plt.plot([1,2,3,4,5])

st.subheader('Static Pyplot')
st.pyplot(fig)

st.subheader('Interactive HTML Plot')
fig_html = mpld3.fig_to_html(fig)
components.html(fig_html, height=600)