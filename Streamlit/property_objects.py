# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 19:58:23 2023

__author__ = "Stuart Goldie"
__copyright__ = "Copyright 2023, University of Kassel"
__credits__ = ["Stuart Goldie","Claudia Backes"]

__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Stuart Goldie"
__email__ = "stuart.goldie@uni-kassel.de"
__status__ = "Production"

Database of material and solvent properties
"""

class material:
    def __init__(self, name, density, thickness, k1, k2):
        self.name = name
        self.density = density
        self.thickness = thickness
        self.k1 = k1
        self.k2 = k2
        
    def __str__(self):
        return f"{self.name} has a density of {self.density} kg m^-3, layer thickness of {self.thickness} m, and aspect ratios l/w {self.k1}, l/h {self.k2}"

class solvent:
    def __init__(self, name, viscosity, density):
        self.name = name
        self.viscosity = viscosity
        self.density = density
        
    def __str__(self):
        return f"{self.name} has a viscosity of {self.viscosity} Pa s and a density of {self.density} kg m^-3"

class rotor:
    inner_radius = None
    def __init__(self, name, angle, radius, diameter, max_speed, manufacturer):
        self.name = name
        self.angle = angle
        self.radius = radius
        self.diameter = diameter
        self.max_speed = max_speed
        self.manufacturer = manufacturer


#Materials#
GRAPHITE = material('Graphite',2260, 3.35E-10, 2.5, 280)
WS2 = material('WS2', 7500, 6.18E-10, 1.75, 46)
MoS2 = material('MoS2', 5060, 6.15E-10, 1.8, 45)
#All materials below need k1 and k2 updating
NiPS3 = material('NiPS3', 3180, 6.6E-10, 2.0, 45)
MnPS3 = material('MnPS3', 2900, 6.8E-10, 2.0, 45)
GaS = material('GaS', 3860, 7.5E-10, 2.0, 45)
BN = material('BN', 2100, 3.4E-10, 2.0, 45)
NiOH2 = material('NiOH2', 4100, 4.6E-10, 2.0, 45)
CuOH2 = material('CuOH2', 3370, 5E-10, 2.0, 45)
RuCl3 = material('RuCl3', 3260, 6E-10, 2.0, 45)
CrTe3 = material('CrTe3', 4700, 1.1E-9, 2.0, 45)


#Solvent Parameters#

WATER = solvent('Water', 1.3E-3, 1000)
NMP = solvent('NMP', 2.3E-3, 1030)
CHP = solvent('CHP', 2.1E-3, 1000)
IPA = solvent('IPA', 3.31E-3, 785)

'''
Database of rotor dimensions
'''

#all tube dimensions given in degrees and millimeters.
#some rotors have outer and inner rings, only those with identical angles are included
JA30_50 = rotor('JA-30.50 Ti', 34, 108, 25, 30000, 'Beckman-Coulter')
JA25_50 = rotor('JA-25.50', 34, 108, 29, 25000, 'Beckman-Coulter')
JA25_15 = rotor('JA-25.15', 25, 106, 18, 25000, 'Beckman-Coulter')
JA25_15.inner_radius = 52
JA21 = rotor('JA-21', 40, 102, 16, 21000, 'Beckman-Coulter')
JA20_1 = rotor('JA-20.1', 23, 115, 18, 20000, 'Beckman-Coulter')
JA20_1.inner_radius = 98
JA20 = rotor('JA-20', 34, 108, 29, 20000, 'Beckman-Coulter')
JA18 = rotor('JA-18', 23, 132, 38, 18000, 'Beckman-Coulter')
JA17 = rotor('JA-17', 25, 123, 29, 17000, 'Beckman-Coulter')
JA14_50 = rotor('JA-14.50', 35, 160, 30, 14000, 'Beckman-Coulter')
JA12 = rotor('JA-12', 35, 144, 30, 12000, 'Beckman-Coulter')
JA10 = rotor('JA-10', 25, 158, 69, 10000, 'Beckman-Coulter')

VFC8_50 = rotor('VFC 8.50', 25, 104, 29, 11360, 'Beckman-Coulter')
VF6_94 = rotor('VF 6.94', 25, 106, 38, 10000, 'Beckman-Coulter')
VF48_2 = rotor('VF 48.2', 53, 100, 11, 13500, 'Beckman-Coulter')
VFC24_15 = rotor('VFC 24.15', 40, 126, 17, 9000, 'Beckman-Coulter')
VF100_2 = rotor('VF 100.2', 45, 163, 11, 6500, 'Beckman-Coulter')
VF100_2.inner_radius=151
VF6_250 = rotor('VF 6.250', 30, 145, 61.5, 5450, 'Beckman-Coulter')

TA15_15 = rotor('TA-15-1.5', 45, 100, 11, 15000, 'Beckman-Coulter')
TA14_50 = rotor('TA-14-50', 25, 96, 29, 14000, 'Beckman-Coulter')

FX6100 = rotor('FX6100', 25, 98, 38, 10200, 'Beckman-Coulter')
FX301_5 = rotor('FX301.5', 45, 100, 11, 16000, 'Beckman-Coulter')
FX241_5P = rotor('FX241.5P', 53, 66, 11, 14800, 'Beckman-Coulter')
FX121_5P = rotor('FX121.5P', 45, 62, 11, 14800, 'Beckman-Coulter')
 

#'VF 48.2' and 'VFC 24.15' has two different angles depending on row; parameters set according to highest angle
#'VFC 24.15' has two different angles due to inner and outer rings
#'JA-18.1' lists two different angles,unclear why.


'''Swinging Bucket rotors and other large capacity unsuitable for model:
'JLA-16.250', 'VS 4.750', 'VS 4.750-96', 'VS 2.5-96', 'TA-10-250'

All TLA range and Type## rotors are used for ultracentrifugation.
MLA used for large volumes in ultracentrifuge
Not included, for now, or routine sample preparation
'''

#TLA_100 = rotor('TLA-100', 30, 38.9, 7, 'Beckman-Coulter') [missing max rpm]

'''Considering Hettich, datasheet and rotor dimensions are complicated by the selection of inserts and adapters for different rotor sizes
I will briefly attempt to summarise the rotor dimensions for un-adapted rotors but do not have the interest in making a full database of all possible inserts'''
#Universal320 Centrifuge
HET_1556 = rotor('1556, 6-place', 35, 115, 35, 9000, 'Hettich')
HET_1613 = rotor('1613, 12-place', 35, 103, 17, 6000, 'Hettich')
HET_1615 = rotor('1615, 12-place', 35, 103, 17, 12000, 'Hettich')
HET_1627 = rotor('1027, 18-place', 45, 18, 17, 14150, 'Hettich')  #clearly the radius is incorrect in Hettich documentation
HET_1553 = rotor('1553, 30-place', 45, 97, 11, 14150, 'Hettich')
HET_1552 = rotor('1552, 24-place', 50, 87, 11, 16000, 'Hettich')

#Rotina380 Centrifuge
HET_1720 = rotor('1720, 6-place', 45, 121, 38, 11000, 'Hettich')
HET_1792 = rotor('1792, 6-place', 45, 122, 38, 11000, 'Hettich')
HET_1789A = rotor('1789-A, 30-place', 45, 97, 11, 15000, 'Hettich')

#Rotina420 Centrifuge
ROT420_4795 = rotor('4795, 4-place', 25, 102, 38, 9500, 'Hettich')
ROT420_4794 = rotor('4794, 6-place', 45, 122, 38, 11000, 'Hettich')
ROT420_4790A = rotor('4790-A, 30-place', 45, 97, 11, 15000, 'Hettich')

#Rotata460 Centrifuge
ROT460_4489A = rotor('4489-A, 30-place', 45, 97, 11, 15000, 'Hettich')
ROT460_5645 = rotor('5645, 6-place', 25, 122, 38, 8500, 'Hettich')
ROT460_5615 = rotor('5615, 6-place', 45, 122, 38, 11500, 'Hettich')

'''ThermoFisher rotors'''
HIGHCON2 = rotor('HIGHConic™ II', 45, 126, 38, 10350, 'ThermoFisher')
HIGHCON3 = rotor('HIGHConic™ III', 45, 120, 30, 9500, 'ThermoFisher')
MICROCLICK24 = rotor('MicroClick24', 45, 85, 11, 16000, 'ThermoFisher')
MICROCLICK18 = rotor('MicroClick18', 45, 102, 17, 15000, 'ThermoFisher')
MICROCLICK30 = rotor('MicroClick30', 45, 99, 11, 14000, 'ThermoFisher')
MICROL30 = rotor('Microliter30', 45, 100, 11, 15200, 'ThermoFisher')
#Fiberlite rotors have no clear radius stated on ThermoFisher website
FIBERLITE15_6 = rotor('Fiberlite™ F15-6', 25, 98, 38, 15000, 'ThermoFisher')
FIBERLITE13_14 = rotor('Fiberlite™ F13-14', 34, 153, 30, 10000, 'ThermoFisher')
FIBERLITE15_8 = rotor('Fiberlite™ F15-8', 25, 104, 30, 14500, 'ThermoFisher')
FIBERLITE21_48 = rotor('Fiberlite™ F21-48', 45, 97, 11, 15200, 'ThermoFisher')

FIBERLITE10 = rotor('Fiberlite™ F10-6', 45, 122, 38, 10500, 'ThermoFisher')

