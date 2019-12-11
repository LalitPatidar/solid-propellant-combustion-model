import csv
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import AutoMinorLocator

mrk = [">", "<", "d", "o", "s", "^"]

fig, ax = plt.subplots()
ax.xaxis.set_minor_locator(AutoMinorLocator())

Pressures = ['./1atm','./5atm','./20atm','./50atm','./90atm']
for P in Pressures: 
    x = np.asarray([])
    temperature = np.asarray([])

    with open(P+'/gas_phase_solution.csv','r') as File:
        reader = csv.reader(File)
        rownum = 0
        
        for row in reader:
            rownum = rownum + 1
            if rownum == 1:
                header = row
            else:
                x = np.append(x,float(row[0]))
                temperature = np.append(temperature,float(row[3]))

    ax.plot(1000*x,temperature,label=P[2:],marker = mrk[Pressures.index(P)],markevery=10)

ax.ticklabel_format(style='sci',axis='x',scilimits=(-1,1))
ax.set_xlim(0.0,1.0)
ax.set_xlabel('x(mm)')
ax.set_ylabel('Temperature (K)')
ax.legend(loc='upper center',bbox_to_anchor=(0.5,1.10),ncol=5)

ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(direction='in',which='both')

#textstr = '\n'.join(['','Computational'])
#ax.text(0.75,0.92,textstr,transform=ax.transAxes)

plt.savefig('Temperature-profiles-gas.pdf')


# Liquid temperature profiles
fig, ax = plt.subplots()
ax.xaxis.set_minor_locator(AutoMinorLocator())

Pressures = ['./1atm','./5atm','./20atm','./50atm','./90atm']
for P in Pressures: 
    x = np.asarray([])
    temperature = np.asarray([])

    with open(P+'/mass_fractions_liquid.txt','r') as File:
        lines = File.readlines()
        
        for line in lines[1:]:
            x = np.append(x,float(line.split()[0]))
            temperature = np.append(temperature,float(line.split()[1]))

    ax.plot(10000*x,temperature,label=P[2:],marker = mrk[Pressures.index(P)],markevery=20)
    plt.axvline(10000*x[-1],linestyle='--',color='black',linewidth=0.5)

#ax.ticklabel_format(style='sci',axis='x',scilimits=(-1,1))
ax.set_xlim(0.0,70)
ax.set_xlabel('x ('r'$\mu$m)')
ax.set_ylabel('Temperature (K)')
ax.legend(loc='upper center',bbox_to_anchor=(0.5,1.10),ncol=5)

ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(direction='in',which='both')

textstr = '\n'.join(['','90atm'])
ax.text(0.062,0.95,textstr,transform=ax.transAxes)

textstr = '\n'.join(['','50atm'])
ax.text(0.1,0.85,textstr,transform=ax.transAxes)

textstr = '\n'.join(['','20atm'])
ax.text(0.18,0.75,textstr,transform=ax.transAxes)

textstr = '\n'.join(['','5atm'])
ax.text(0.41,0.55,textstr,transform=ax.transAxes)

textstr = '\n'.join(['','1atm'])
ax.text(0.90,0.40,textstr,transform=ax.transAxes)

plt.savefig('Temperature-profiles-liquid.pdf')
