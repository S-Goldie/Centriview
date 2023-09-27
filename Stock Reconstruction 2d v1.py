# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 16:03:16 2023

@author: Stuart Goldie
"""

#Stock reconstruction based on sediment fractions

#import data sets as numpy arrays
#create 3d histograms from the data as arrays of counts with identical bin ranges for each data set
#correct each histogram with the weighting yield from each sediment
#sum corrected arrays together
#use these corrected arrays as counts for production of data points using random number generator to populate a data cloud

import numpy as np
import matplotlib.pyplot as plt

#Specify the prefixes used in the file names. These will remain as data labels:
fraction_list = ('400', '1000', '5000', '10000', '30000')

#Specify the mass yield from each sediment for weighting of fractions:
mass_list = (175.9, 9.76, 7.71, 2.25, 2.0)

def load_data(filename):
    """Import data as 2d array. a is flake index, b is measurement: 0=length, 1=width, 2=layer number"""
    return np.genfromtxt(filename, skip_header=1, delimiter=",")

experimental_data = {}
for speed in fraction_list:
    experimental_data[str(speed)] = load_data(f'{speed}_AFM.csv')

def massvol_weighting(data_set, mass):
    flake_volume = np.mean(data_set[:,0] * data_set[:,1] * data_set[:,2])
    return mass / flake_volume

raw_weighting_list = []
for i, speed in enumerate(fraction_list):
    raw_weighting_list.append(massvol_weighting(experimental_data[str(speed)], mass_list[i]))

weighting_list = []
for j in raw_weighting_list:
    weighting_list.append(int(j / min(raw_weighting_list)))

histogram_bins = [25, 50]
histogram_range = [[0,1000],[0,50]]

#%%Note on Histogram bin definition from numpy docs.
#https://numpy.org/doc/stable/reference/generated/numpy.histogram.html
#All but the last (righthand-most) bin is half-open. In other words, if bins is:
#[1, 2, 3, 4] then the first bin is [1, 2) (including 1, but excluding 2) and the second [2, 3).
#The last bin, however, is [3, 4], which includes 4.
#This boundary definition is consistent with the random number generator used to produce the data cloud
#%%

def weight_histogram(data_set, weighting):
    """Take data_set and return 3d histogram array multiplied by effective density weighting"""
    histogram = np.histogram2d(data_set[:,0], data_set[:,2], histogram_bins, histogram_range)
    return histogram[0] * weighting

def histogram_reconstruction(*args):
    sum_list = []
    for hist in args:
        sum_list.append(hist)
    return sum(sum_list)

weighted_histograms = {}

#Calculate the first histogram explicitly so the bin lists are callable later.
histogram_first = np.histogram2d(experimental_data[str(fraction_list[0])][:,0], experimental_data[str(fraction_list[0])][:,2], histogram_bins, histogram_range)
flake_volume = np.mean(experimental_data[str(fraction_list[0])][:,0] * experimental_data[str(fraction_list[0])][:,1] * experimental_data[str(fraction_list[0])][:,2])
weighted_histograms[str(fraction_list[0])] = (histogram_first[0] * weighting_list[0])

for index, speed in enumerate(fraction_list[1:]):
    weighted_histograms[str(speed)] = weight_histogram(experimental_data[str(speed)], weighting_list[index])

complete_histogram = weighted_histograms[str(fraction_list[0])]
for i in fraction_list[1:]:
    complete_histogram = complete_histogram + weighted_histograms[str(i)]

#Populate a new data set using random numbers constrained by the bin limits of the histogram.
#The histogram counts determine how many random numbers to generate for each bin range
reconstructed_lengths = []
reconstructed_thickness = []
for index,value in np.ndenumerate(complete_histogram):
    for _ in range(int(value)):
        reconstructed_lengths.append(np.random.default_rng().uniform(histogram_first[1][index[0]], histogram_first[1][index[0] + 1]))
        reconstructed_thickness.append(np.random.default_rng().uniform(histogram_first[2][index[1]], histogram_first[2][index[1] + 1]))

output_data = np.asarray([reconstructed_lengths, reconstructed_thickness])
processed_headings = ('Length, Layers')
processed_output = [list(n) for n in zip(*output_data)]
np.savetxt('Reconstructed 2D Data.csv', (processed_output), delimiter=',', header=processed_headings)

plt.figure('Reconstructed Data Cloud', figsize=(4, 3), dpi=300)
plt.scatter(reconstructed_lengths, reconstructed_thickness, s=1, color='k')
plt.xlabel('Length / nm')
plt.ylabel('Number of Layers')
plt.tight_layout()
plt.savefig('Reconstructed Data.png')

plt.figure('Original Data Clouds', figsize=(4, 3), dpi=300)
for speed in fraction_list:
    plt.scatter(experimental_data[str(speed)][:,0], experimental_data[str(speed)][:,2], s=1, label=speed)
plt.xlabel('Length / nm')
plt.ylabel('Number of Layers')
plt.legend()
plt.tight_layout()
plt.savefig('Original Data.png')
