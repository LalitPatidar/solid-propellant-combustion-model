import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import ScalarFormatter

Tinit = 298 # Kelvin

# Melt layer thickness
PPrasad = [1.0, 5.0, 20.0, 50.0, 90.0]
RPrasad = [65, 25, 14, 7, 4]

PZenin = [1.0, 5.0, 20.0, 70.0]
RZenin = [70, 30, 20, 15]

version = 'Present work'

Pversion = []
Xmversion = []

for P in ['1', '5', '20', '50', '90']:
    file_name = './'+ P +'atm/mass_fractions_liquid.txt'
    with open(file_name,'r') as File:
        lines = File.readlines()
        
        Xm = (float(lines[-1].split()[0]))*10000.00 # micrometer
        
        Pversion.append(float(P))
        Xmversion.append(Xm)

fig, ax = plt.subplots()
ax.xaxis.set_minor_locator(AutoMinorLocator())

ax.plot(Pversion,Xmversion, label = version, marker='*')
ax.plot(PPrasad,RPrasad, linestyle='--', label = 'Prasad et al.$^a$', marker='s')
ax.plot(PZenin,RZenin,linestyle='None', label = 'Zenin$^b$', marker='d')

ax.set_xlabel('Pressure (atm)')
ax.set_ylabel('Melt layer thickness ('r'$\mu$m)')
ax.ticklabel_format(style='plain')

ax.legend(loc='upper center',bbox_to_anchor=(0.5,1.10),ncol=4)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(direction='in',which='both')


plt.savefig('Xmelt-vs-P.pdf')

# Delta layer thickness
PPrasad = [1.0, 5.0, 20.0, 50.0, 90.0]
RPrasad = [195, 53, 23, 16, 6]

PZenin = [1.0, 5.0, 20.0, 70.0]
RZenin = [250, 70, 35, 28]

version = 'Present work'

Pversion = []
Xmversion = []

for P in ['1', '5', '20', '50', '90']:
    file_name = './'+ P +'atm/delta_layer_solution.txt'
    with open(file_name,'r') as File:
        lines = File.readlines()
        
        Xm = (float(lines[-1].split()[0]))*10000.00 # micrometer
        
        Pversion.append(float(P))
        Xmversion.append(Xm)

fig, ax = plt.subplots()
ax.xaxis.set_minor_locator(AutoMinorLocator())

ax.plot(Pversion,Xmversion, label = version, marker='*')
ax.plot(PPrasad,RPrasad, linestyle='--', label = 'Prasad', marker='s')
ax.plot(PZenin,RZenin,linestyle='None', label = 'Zenin', marker='d')

ax.set_xlabel('Pressure (atm)')
ax.set_ylabel('Delta layer thickness ('r'$\mu$m)')
ax.ticklabel_format(style='plain')

ax.legend(loc='upper center',bbox_to_anchor=(0.5,1.10),ncol=4)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(direction='in',which='both')


plt.savefig('Xdelta-vs-P.pdf')

# Decomposition % in liquid
PPrasad = [1.0, 5.0, 20.0, 50.0, 90.0]
RPrasad = [59, 49, 11, 0.2, 0]

#PZenin = [1.0, 5.0, 20.0, 70.0]
#RZenin = [250, 70, 35, 28]

version = 'Present work'

Pversion = []
Xmversion = []

for P in ['1', '5', '20', '50', '90']:
    file_name = './'+ P +'atm/mass_fractions_liquid.txt'
    with open(file_name,'r') as File:
        lines = File.readlines()
        
        Xm = 100.0 - (float(lines[-1].split()[3]))*100.00 # %Decomposition in liquid
        
        Pversion.append(float(P))
        Xmversion.append(Xm)

fig, ax = plt.subplots()
ax.xaxis.set_minor_locator(AutoMinorLocator())

ax.plot(Pversion,Xmversion, label = version, marker='*')
ax.plot(PPrasad,RPrasad, linestyle='--', label = 'Prasad', marker='s')
#ax.plot(PZenin,RZenin,linestyle='None', label = 'Zenin', marker='d')

ax.set_xlabel('Pressure (atm)')
ax.set_ylabel('Mass fraction of HMX at liquid-gas interface (%)')
ax.ticklabel_format(style='plain')

ax.legend(loc='upper center',bbox_to_anchor=(0.5,1.10),ncol=4)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(direction='in',which='both')

plt.savefig('Y_HMX_liq-vs-P.pdf')
