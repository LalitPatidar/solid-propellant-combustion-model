import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import ScalarFormatter
import scipy.optimize as optimization

Tinit = 298 # Kelvin

# Burn-rate
PPrasad = [1.0, 5.0, 20.0, 50.0, 90.0]
RPrasad = [0.035, 0.153, 0.401, 0.896, 1.509]

PZenin = [1.0, 5.0, 20.0, 70.0]
RZenin = [0.035, 0.150, 0.40, 1.0]

PAtwood = [2.37,   3.4,    6.80,  10.2,   10.2,   13.61, 13.61, 13.61,    20.41,  27.22,  27.22,  34.02, 47.63,  61.24, 61.24, 102.0, 102.0  ]
RAtwood = [0.0955, 0.1011, 0.167, 0.2921, 0.2459, 0.2616, 0.3404, 0.2921, 0.4445, 0.5994, 0.5461, 0.655, 0.8357, 1.143, 1.034, 1.4895, 1.6759 ]

PCaltech = [1.0, 5.0, 20.0, 50.0, 90.0]
RCaltech = [0.03065, 0.1546, 0.5979, 1.4364, 2.59]

PYetter = [1.0, 5.0, 20.0, 50.0, 90.0]
RYetter = [0.0387, 0.1760, 0.5820, 1.2459, 2.005]


version = 'Present work (T$_0$ = 298K)'

Pversion = []
Rversion = []
for P in ['1', '5', '20', '50', '90']:
    file_name = './'+ P +'atm/run.log'
    with open(file_name,'r') as File:
        lines = File.readlines()
        
        r = (float(lines[-1].split()[1]) + float(lines[-1].split()[2]))/2.0/1.9
        
        Pversion.append(float(P))
        Rversion.append(r)

def logRb(logP, n, lna):
    return  n*logP + lna

coeff = optimization.curve_fit(logRb,np.log(Pversion),np.log(Rversion))

n = coeff[0][0]
a = np.exp(coeff[0][1])

print(a)
print(n)

rb = a*Pversion**(n)

fig, ax = plt.subplots()
ax.xaxis.set_minor_locator(AutoMinorLocator())

ax.plot(Pversion,Rversion, label = version, marker='*')
#ax.plot(Pversion,rb, label = 'log-Rb', marker='>')
ax.plot(PCaltech,RCaltech, linestyle='--', label = 'Caltech-mechanism$^a$', marker='^')
ax.plot(PYetter,RYetter, linestyle=':',label = 'Yetter-mechanism$^b$', marker='v')
ax.plot(PPrasad,RPrasad, linestyle='None', label = 'Prasad et al.$^c$', marker='s',markersize=8)
ax.plot(PZenin,RZenin,linestyle='None', label = 'Zenin$^d$', marker='d',markersize=8)
ax.plot(PAtwood,RAtwood,linestyle='None', label = 'Atwood et al.$^e$',marker='o')

ax.set_xlabel('Pressure (atm)')
ax.set_ylabel('Burn-rate (cm/s)')
ax.ticklabel_format(style='plain')

ax.legend(loc='upper center',bbox_to_anchor=(0.5,1.22),ncol=2)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(direction='in',which='both')

plt.xscale('log')
plt.yscale('log')

ax.xaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.grid(which='major',linestyle='-')
ax.grid(which='minor',linestyle=':')

textstr = '\n'.join(['r$_b$ = ap$^n$\n','a = %6.3f\nn = %6.3f' %(a,n)])
ax.text(0.7,0.4,textstr,transform=ax.transAxes)

plt.savefig('Rb-vs-aPn.pdf', bbox_inches='tight')


