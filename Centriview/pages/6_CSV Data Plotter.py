# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 17:54:00 2024

@author: Stuart Goldie
"""

import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

##_FRONT END CODE_##
st.title('CSV Data Plotter')
st.markdown("""
This tool allows you to import two CSV arrays and visualize the data in multiple ways:

1. **Scatter Plot**: Plots data from the first CSV as x-axis and the second CSV as y-axis (index-matched)
2. **Ratio Trend**: Calculates the ratio of the two values (y/x) and plots the trend as a function of column number

Upload two CSV files below to get started.
""")

st.subheader('Upload CSV Files')
st.write('Upload two CSV files. The data will be matched by index (row number).')

col1, col2 = st.columns(2)

with col1:
    file1 = st.file_uploader('Upload First CSV (X-axis data)', type=['csv'], key='file1')
    
with col2:
    file2 = st.file_uploader('Upload Second CSV (Y-axis data)', type=['csv'], key='file2')

if file1 is not None and file2 is not None:
    try:
        # Read the CSV files
        df1 = pd.read_csv(file1)
        df2 = pd.read_csv(file2)
        
        st.success('Files uploaded successfully!')
        
        # Display file information
        col3, col4 = st.columns(2)
        with col3:
            st.write(f'**First CSV Shape:** {df1.shape[0]} rows × {df1.shape[1]} columns')
            with st.expander('Preview First CSV'):
                st.dataframe(df1.head())
        
        with col4:
            st.write(f'**Second CSV Shape:** {df2.shape[0]} rows × {df2.shape[1]} columns')
            with st.expander('Preview Second CSV'):
                st.dataframe(df2.head())
        
        # Let user select columns
        st.subheader('Select Data Columns')
        col5, col6 = st.columns(2)
        
        with col5:
            x_column = st.selectbox('Select column from first CSV (X-axis):', df1.columns, key='x_col')
        
        with col6:
            y_column = st.selectbox('Select column from second CSV (Y-axis):', df2.columns, key='y_col')
        
        # Extract the selected columns
        x_data = df1[x_column].values
        y_data = df2[y_column].values
        
        # Ensure data lengths match
        min_length = min(len(x_data), len(y_data))
        if len(x_data) != len(y_data):
            st.warning(f'Data lengths differ. Using first {min_length} elements from each array.')
            x_data = x_data[:min_length]
            y_data = y_data[:min_length]
        
        # Remove any NaN values
        valid_indices = ~(np.isnan(x_data) | np.isnan(y_data))
        x_data_clean = x_data[valid_indices]
        y_data_clean = y_data[valid_indices]
        
        if len(x_data_clean) == 0:
            st.error('No valid data points found. Please check your CSV files.')
        else:
            st.success(f'Using {len(x_data_clean)} valid data points.')
            
            # Calculate ratios (avoiding division by zero)
            with np.errstate(divide='ignore', invalid='ignore'):
                ratios = y_data_clean / x_data_clean
                # Replace inf and -inf with NaN
                ratios = np.where(np.isinf(ratios), np.nan, ratios)
            
            # Remove NaN ratios for plotting
            valid_ratio_indices = ~np.isnan(ratios)
            ratios_clean = ratios[valid_ratio_indices]
            indices_clean = np.arange(len(x_data_clean))[valid_ratio_indices]
            
            # Create plots
            st.subheader('Data Visualization')
            
            # Scatter plot
            st.markdown('#### Scatter Plot: Y vs X')
            fig1 = go.Figure()
            fig1.add_trace(go.Scatter(
                x=x_data_clean, 
                y=y_data_clean, 
                mode='markers',
                marker=dict(size=8, color='blue', opacity=0.6),
                name='Data Points'
            ))
            fig1.update_layout(
                xaxis_title=f'{x_column} (from first CSV)',
                yaxis_title=f'{y_column} (from second CSV)',
                title='Scatter Plot: Second CSV vs First CSV',
                height=500,
                hovermode='closest'
            )
            st.plotly_chart(fig1, use_container_width=True)
            
            # Ratio trend plot
            st.markdown('#### Ratio Trend: Y/X vs Index')
            fig2 = go.Figure()
            fig2.add_trace(go.Scatter(
                x=indices_clean, 
                y=ratios_clean, 
                mode='lines+markers',
                marker=dict(size=6, color='red'),
                line=dict(color='red', width=2),
                name='Ratio Trend'
            ))
            fig2.update_layout(
                xaxis_title='Column Number (Index)',
                yaxis_title=f'Ratio ({y_column} / {x_column})',
                title='Ratio Trend as Function of Column Number',
                height=500,
                hovermode='closest'
            )
            st.plotly_chart(fig2, use_container_width=True)
            
            # Statistics
            st.subheader('Statistics')
            col7, col8, col9 = st.columns(3)
            
            with col7:
                st.metric('Mean Ratio', f'{np.mean(ratios_clean):.4f}')
                st.metric('Median Ratio', f'{np.median(ratios_clean):.4f}')
            
            with col8:
                st.metric('Std Dev Ratio', f'{np.std(ratios_clean):.4f}')
                st.metric('Min Ratio', f'{np.min(ratios_clean):.4f}')
            
            with col9:
                st.metric('Max Ratio', f'{np.max(ratios_clean):.4f}')
                st.metric('Valid Points', f'{len(ratios_clean)}')
            
            # Download processed data
            st.subheader('Download Processed Data')
            
            # Create dataframe with results
            results_df = pd.DataFrame({
                'Index': indices_clean,
                f'{x_column}_X': x_data_clean[valid_ratio_indices],
                f'{y_column}_Y': y_data_clean[valid_ratio_indices],
                'Ratio_Y_X': ratios_clean
            })
            
            csv = results_df.to_csv(index=False)
            st.download_button(
                label="Download Results as CSV",
                data=csv,
                file_name="processed_data_with_ratios.csv",
                mime="text/csv"
            )
            
    except Exception as e:
        st.error(f'Error processing files: {str(e)}')
        st.write('Please ensure your CSV files are properly formatted.')

else:
    st.info('Please upload both CSV files to begin analysis.')

st.divider()
st.image('https://2d-printable.eu/storage/sites/18/2023/02/Funding_statement-1-768x161.png', caption=None, width=254)
