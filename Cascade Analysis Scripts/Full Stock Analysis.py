# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 01:12:34 2023

Calculation of fraction remaining for each measured nanosheet within a dispersion, plotting 2D histograms for each fraction following a cascade taking into account the evolving fractions of each material.

@author: Goldie
"""

import numpy as np
import matplotlib.pyplot as plt
import csv

# Load data from CSV files. The data must be stored with each flake on a 
#separate row and columns in the order: length , width, layer number. (NOT THICKNESS)
#Data stored as [a,b] where a specifies an individual flake, b specifies the 
#measurement: 0=length, 1=width, 2=layer number
def load_data(filename):
    return np.genfromtxt(filename, skip_header=1, delimiter=",")

#The speed list is in rpm; file names for the experimental data of each fraction must match to this list
#The first speed is assumed to be the initial removal of unexfoliated flakes, and therefore the sediment is not analysed
speed_list = (900, 1800, 2900, 6500, 9000, 15000)
material = 'WS2'
solvent = 'NMP'

stock_data = load_data('Stock_AFM.csv')
experimental_data = {}
for speed in speed_list[1:]:
    experimental_data[str(speed)] = load_data(f'{speed}_AFM.csv')

#Surfactant parameters are for sodium cholate in water. For graphene, different
#surfactant and solvent systems will cause signficant error. Recommend setting the thickness to zero in NMP.
pS = 1595           #surfactant layer density, from Hirsham using anhydrous cholate on surface
d0 = 4.25E-10       #surfactant layer thickness, from Hirsham

if material == 'Graphite':
    d1 = 3.35E-10   #monolayer thickness
    k1 = 2.5        #L/W
    k2 = 280        #L/h
    pNS = 2260      #nanosheet density
if material == 'WS2':
    d1 = 6.18E-10    #monolayer of WS2, J. A. Wilson and A. D. Yoffe, Advances in Physics 18, 193 (1969).
    k1 = 1.75       #L/W
    k2 = 46         #L/h
    pNS = 7500      #structure density
if material == 'MoS2':
    d1 = 6.15E-10   #monolayer thickness, J. A. Wilson and A. D. Yoffe, Advances in Physics 18, 193 (1969).
    k1 = 1.8       
    k2 = 45         
    pNS = 5060      #structure density
if material == 'NiPS3':
    d1 = 6.6E-10
    pNS = 3180
if material == 'MnPS3':
    d1 = 6.8E-10
    pNS = 2900
if material == 'GaS':
    d1 = 7.5E-10
    pNS = 3860
if material == 'BN':
    d1 = 3.4E-10
    pNS = 2100
if material == 'NiOH2':
    d1 = 4.6E-10
    pNS = 4100
if material == 'CuOH2':
    d1 = 5E-10
    pNS = 3370
if material == 'RuCl3':
    d1 = 6E-10
    pNS = 3260
if material == 'CrTe3':
    d1 = 1.1E-9
    pNS = 4700    

if solvent == 'Water':
    n = 1.3E-3
    pL = 1000
if solvent == 'NMP':
    n = 2.3E-3
    pL = 1030
if solvent == 'CHP':
    n = 2.1E-3
    pL = 1000
if solvent == 'IPA':
    n = 3.31E-3
    pL = 785

def supernatant_frac(length,width,layer_number, time_hour, rpm):
    """"Calculate the proportion of flakes of a given size remaining in the supernatant"""
    r1, r2 = 0.071, 0.1
    #r2 = 0.085
    time_s = 3600 * time_hour
    h = layer_number * d1
    w = (rpm * 2 * np.pi) / 60
    a = time_s * w * w * (h * (pNS - pL) + 2 * d0 * (pS - pL )) * np.sqrt( length * width / 2) / (12 * n * np.cbrt(3 / (4 * np.pi)))
    return max((r2*np.exp(-a)-r1)/(r2-r1), 0)

def cutoff_area(layer_number, time_hour, rpm):
    """Calculate the absolute maximum flake size, above which 100% sedimentation is guaranteed"""
    r1, r2 = 0.071, 0.1
    time_s = 3600 * time_hour
    h = layer_number * d1
    w = (rpm * 2 * np.pi) / 60
    return (((np.log(r2 / r1) * (12 * n * np.cbrt(3 / (4 * np.pi))) / (time_s * w * w * (h * (pNS - pL) + 2 * d0 * (pS - pL)))) ** 2) * 1e18)

fractions_calculated = np.zeros([len(speed_list), len(stock_data)])
supernatants_compounded = np.zeros([len(speed_list), len(stock_data)])
sediments_compounded = np.zeros([len(speed_list)-1, len(stock_data)])

sedimented_flakes = np.zeros([(len(speed_list)-1)*2, len(stock_data)])

#summary statistics saved in a dictionary. keys: central rcf, experimental mean thickness, 
#experimental mean length, predicted mean thickness, predicted mean length
summary_statistics = []

#Bins for plotting the bivariate histograms. No affect on numerical analysis.
bin_range = [[0,800],[0,40]]

for i, speed in enumerate(speed_list):
    for j in np.arange(len(stock_data)):
        fractions_calculated[i, j] = supernatant_frac(stock_data[j,0] * 1e-9, stock_data[j,1] * 1e-9, stock_data[j,2], 2, speed)
    if i == 0:
        supernatants_compounded[0] = fractions_calculated[0]
    elif i > 0:
        supernatants_compounded[i] = fractions_calculated[i] * supernatants_compounded[i-1]
        sediments_compounded[i-1] = supernatants_compounded[i-1] - supernatants_compounded[i]
        
        fraction_statistics = {'central rcf': (speed + speed_list[i-1]) / 2,
                               'exp <L>': np.average(experimental_data[str(speed)][:,0]),
                               'exp <N>': np.average(experimental_data[str(speed)][:,2]),
                               'Comp <L>': np.average(stock_data[:,0], weights=sediments_compounded[i-1]),
                               'Comp <N>': np.average(stock_data[:,2], weights=sediments_compounded[i-1])}
        summary_statistics.append(fraction_statistics)
        
        plt.figure(f'Predicted Sediment {speed_list[i-1]}-{speed_list[i]}', figsize=(3.3, 3))
        plt.hist2d(stock_data[:,0], stock_data[:,2], bins=[50,50], range=bin_range, cmin=0.1, weights=sediments_compounded[i-1])
        plt.xlabel("L / nm",fontsize=14)
        plt.ylabel("Layer Number",fontsize=14)
        plt.tight_layout()
        plt.colorbar()
        #plt.savefig(f'Predicted Sediment {speed_list[i-1]}-{speed_list[i]}.png', format='png', dpi=600)
        
        layer_numbers = []
        areas = []
        for num, weight in np.ndenumerate(sediments_compounded[i-1]):
            if weight != 0:
                layer_numbers.append(stock_data[num, 2][0])
                areas.append(stock_data[num, 0][0]*stock_data[num,1][0]/2)
            
        try:
            sedimented_flakes[(i-1)*2] = np.array(layer_numbers)
            sedimented_flakes[((i-1)*2)+1] = np.array(areas)
        except ValueError:
            invalid_points = len(stock_data) - len(layer_numbers)
            layer_array = np.append(np.array(layer_numbers), np.zeros([invalid_points]))
            areas_array = np.append(np.array(areas), np.zeros([invalid_points]))
            sedimented_flakes[(i-1)*2] = layer_array
            sedimented_flakes[((i-1)*2)+1] = areas_array            
                
        dummy_layer_numbers = np.arange(start=1, stop=max(layer_numbers), step=0.1)
        cutoff_values = [cutoff_area(m, 2, speed_list[i - 1]) for m in dummy_layer_numbers]
           
        plt.figure(f'Predicted Data Cloud {speed_list[i-1]}-{speed_list[i]}')
        plt.scatter(layer_numbers, areas, s=4, marker='x', color='blue', label='Predicted')
        plt.scatter(experimental_data[str(speed)][:,2], experimental_data[str(speed)][:,0]*experimental_data[str(speed)][:,1]/2, 
                    s=4, marker='x', color='red', label='Experimental')
        plt.plot(dummy_layer_numbers, cutoff_values, 'r--')
        plt.ylabel("Area $nm^2$",fontsize=14)
        plt.xlabel("Layer Number",fontsize=14)
        plt.ylim(0,max(areas)+1000)
        plt.tight_layout()
        #plt.savefig(f'Predicted Data Cloud {speed_list[i-1]}-{speed_list[i]}.png', format='png', dpi=600)

     
#plt.figure('Stock Data Cloud')
#plt.scatter(stock_data[:, 2], stock_data[:, 0]*stock_data[:,1]/2, s=4, marker='x', color='blue')
#plt.ylabel("Area $nm^2$",fontsize=14)
#plt.xlabel("Layer Number",fontsize=14)
#plt.tight_layout()
#plt.savefig('Stock Data Cloud.png', format='png', dpi=600)
        
with open('Summary Stats.csv', 'wt') as fout:        #create file to save histogram data into using the original data file name
    cout = csv.DictWriter(fout, ['central rcf', 'exp <N>', 'Comp <N>', 'exp <L>', 'Comp <L>'])
    cout.writeheader()              #Using headers consistent with dictionary keys
    cout.writerows(summary_statistics)    #each dictionary in the results list is converted into one line in the results table

output_fractions = np.concatenate((np.transpose(stock_data), fractions_calculated, sediments_compounded), axis=0)
fraction_headings = 'Length, Width, Layer Number, ' + ', '.join([f'{speed}rpm Fraction Weighting' for speed in speed_list]) + ', ' + ', '.join(
    [f'{speed1}-{speed2}rpm Sediment Weighting' for speed1, speed2 in zip(speed_list[:-1], speed_list[1:])])
processed_output = [list(n) for n in zip(*output_fractions)]
np.savetxt('Fractions Analysed.csv', (processed_output), delimiter=',', header=fraction_headings)

sedimented_headings = ', '.join([f'{speed1}-{speed2}rpm N, {speed1}-{speed2}rpm A' for speed1, speed2 in zip(speed_list[:-1], speed_list[1:])])
flakes_output = [list(m) for m in zip(*sedimented_flakes)]
np.savetxt('Flakes in Fractions.csv', (flakes_output), delimiter=',', header=sedimented_headings)


def plot_histogram_heatmap(data, title, filename):
    """Plot a bivariant histograms as a light/dark heatmap of length and layer number"""
    plt.figure(title, figsize=(3.3, 3))
    plt.hist2d(data[:, 0], data[:, 2], bins=[50, 50], range=bin_range, cmin=1)
    plt.xlabel("L / nm", fontsize=14)
    plt.ylabel("Layer Number", fontsize=14)
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(filename, dpi=600)

#plot_histogram_heatmap(stock_data, 'Stock Dispersion', 'Stock Dispersion.png')
#for speed in speed_list[1:]:
#    plot_histogram_heatmap(experimental_data[str(speed)], f'Measured Sediment: {speed}', f'Measured Sediment: {speed}.png')
