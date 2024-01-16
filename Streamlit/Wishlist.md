# Simple wishlist of plots and features for Streamlit App

## Features
* selectable solvent, material and centrifuge conditions with live plots (conditions being: time, speed, filling height, possibility of changing rotor) PROGRESS: real time plotting is faster for contour plot.
* temperature and solvent mixtures for the contour plot.
* Q: how best to include the rotor distances R1 and R2 for the general contour plot. Use ratio or actual fill height and rotor angle, or even import the rotor angle from properties and use fill height as an input?
* conversion of centrifuge condtions between different materials and solvents. Currently this is working progress for idential rotor geometry, but solving for different rotor geometry is proving very difficult.

## Plots
* contour plot of fraction population
* contour plots as above, except using a 2D log-normal function to model the changes to an approximate distribution
* plot of most likely flake size, LN, N, L etc. for different centrifuge speeds
