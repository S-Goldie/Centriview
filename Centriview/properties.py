# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 22:18:40 2024

@author: Stuart Goldie
"""
euler = 2.718281828459045235360287471352

##_BACKEND DATA AND FUNCTIONS_##
#Storage of parameters for calculations
#"Keys" : ['Density / kg m-3', 'Layer Thickness / m', 'k1', 'k2']
MATERIALS = {
    "Graphene" :[2260, 3.35E-10, 2.5, 280],
    "WS2" :     [7500, 6.18E-10, 1.75, 46],
    "MoS2" :    [5060, 6.15E-10, 1.8, 45],
    "NiPS3" :   [3180, 6.6E-10, 2.0, 45],
    "MnPS3" :   [2900, 6.8E-10, 2.0, 45],
    "GaS" :     [3860, 7.5E-10, 2.0, 45],
    "BN" :      [2100, 3.4E-10, 2.0, 45],
    "NiOH2" :   [4100, 4.6E-10, 2.0, 45],
    "CuOH2" :   [3370, 5E-10, 2.0, 45],
    "RuCl3" :   [3260, 6E-10, 2.0, 45],
    "CrTe3" :   [4700, 1.1E-9, 2.0, 45],
    "Black Phosphorus" :    [2690, 5E-10],
    "Antimony" :[6694, 5.63E-10]
    }

#"Keys" : ['density_m', 'density_c', 'viscosity_A', 'viscosity_Ea']
#Constants produce values in common units, not SI. Conversion done in the function
SOLVENTS = {
    "IPA" :     [-8.97149E-4, 0.80402, 3.05818E-4, 21846.037],
    "THF" :     [-0.00108, 0.90959, 0.0445145, 5845.749],
    "NMP" :     [-8.81678E-4, 1.04999, 0.015821, 11553.57],
    "DLG (Cyrene™)" : [-9.73846E-4, 1.27159, 9.089E-04, 23410],
    "Toluene" : [-9.47951E-4, 0.88605, 0.0185636, 8438.939],
    "DMF" :     [-9.62251E-4, 0.96862, 0.021271, 9047.422],
    "Ethanol" : [-9.11871E-4, 0.80885, 0.003494, 14223.19],
    "Methanol" : [-9.38703E-4, 0.81022, 0.0105775, 9821.7161],
    "DCM" : [-0.0021, 1.36941, 0.0265109, 6821.8779],
    "Butan-2-ol" : [-8.684E-4, 0.82422, 0.000112, 25247.25]
    }

solvent_list = ['Water'] + list(SOLVENTS.keys()) + ['IPA:Water Mixture','Other']
materials_list = list(MATERIALS.keys()) + ['Other']

def mixed_property_function(w, T, parameter):
    """Returns the density or viscosity depending on parameter flag, of a water-IPA binary mixture using crude 2d polynomial"""
    if parameter == 'Density':
        c = [0.999581, -0.0867712, -0.000109338, -0.171483, -0.00165038, -1.24842e-06, 0.0592866, 0.00105765, -2.63718e-07, -1.86795e-08]
    elif parameter == 'Viscosity':
        c = [1.78798, 14.9948, -0.0677034, -13.7734, -0.283823, 0.00133142, 1.61501, 0.194239, 0.000875951, -9.46676e-06]
    return(c[0] + w*c[1] + T*c[2] + (w**2)*c[3] + w*T*c[4] + (T**2)*c[5] + (w**3)*c[6] + (w**2)*T*c[7] + w*(T**2)*c[8] + (T**3)*c[9])

def solvent_viscosity(temp, solvent):
    """Returns the solvent viscosity in SI units given correct constants"""
    R = 8.314
    A = SOLVENTS[solvent][2]
    E_a = SOLVENTS[solvent][3]
    Kelvin = 273.15 + temp
    centipoise = A*euler**(E_a/(R*Kelvin))
    return(centipoise/1000)

def solvent_density(temp, solvent):
    """Returns the solvent density in g cm-3"""
    m = SOLVENTS[solvent][0]
    c = SOLVENTS[solvent][1]
    return((m * temp + c)*1000)

#water follows different trends
def water_viscosity(temp):
    """Returns viscosity of water using third order exponential"""
    kelvin = temp+273.15
    A = 2.3168E-12
    B = 4729.66128
    C = 0.04265
    D = -2.14282E-5
    centipoise = A*euler**(B/kelvin+C*kelvin+D*kelvin**2)
    return(centipoise / 1000)

def water_density(temp):
    """Returns the density of water using third order polynomial"""
    C = 0.99997
    B1 = 3.1087E-5
    B2 = -6.45401E-6
    B3 = 2.09654E-8
    return((C+B1*temp+B2*temp**2+B3*temp**3)*1000)
