#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import FormatStrFormatter 

x = np.asarray([])
temperature = np.asarray([])
CH2O = np.asarray([])
CO = np.asarray([])
HCN = np.asarray([])
N2O = np.asarray([])
NO = np.asarray([])
NO2 = np.asarray([])
H2O = np.asarray([])
CO2 = np.asarray([])
HMX = np.asarray([])
c_HONO = np.asarray([])
t_HONO = np.asarray([])
N2 = np.asarray([])
ONNO2 = np.asarray([])

with open('./1atm/mass_fractions_liquid.txt','r') as File:
    lines = File.readlines()
    
    header = lines[0].split()
    
    iCH2O = header.index('CH2O')
    iCO = header.index('CO')
    iHCN = header.index('HCN')
    iN2O = header.index('N2O')
    iNO = header.index('NO')
    iNO2 = header.index('NO2')
    iH2O = header.index('H2O')
    iCO2 = header.index('CO2')
    iHMX = header.index('HMX')
    ic_HONO = header.index('c_HONO')
    it_HONO = header.index('t_HONO')
    iN2 = header.index('N2')
    iONNO2 = header.index('ONNO2')
    
    new_lines = lines[1::1]

    for line in new_lines:
        x = np.append(x,10.0*float(line.split()[0])) # x in mm
        temperature = np.append(temperature,float(line.split()[1]))
        CH2O = np.append(CH2O,float(line.split()[iCH2O]))
        CO = np.append(CO,float(line.split()[iCO]))
        HCN = np.append(HCN,float(line.split()[iHCN]))
        N2O = np.append(N2O,float(line.split()[iN2O]))
        NO = np.append(NO,float(line.split()[iNO]))
        NO2 = np.append(NO2,float(line.split()[iNO2]))
        H2O = np.append(H2O,float(line.split()[iH2O]))
        CO2 = np.append(CO2,float(line.split()[iCO2]))
        HMX = np.append(HMX,float(line.split()[iHMX]))
        c_HONO = np.append(c_HONO,float(line.split()[ic_HONO]))
        t_HONO = np.append(t_HONO,float(line.split()[it_HONO]))
        N2 = np.append(N2,float(line.split()[iN2]))
        ONNO2 = np.append(ONNO2,float(line.split()[iONNO2]))

markpts = 21
fig, ax1 = plt.subplots()

ax1.plot(x,HMX,label='HMX',marker='*',color='black',markevery=markpts)
ax1.set_ylabel('Mass fraction of HMX in liquid')
ax1.legend(loc='center left',ncol=1)
#ax1.ticklabel_format(style='sci',axis='y',scilimits=(-1,1))
ax1.set_ylim(0.72,1.0)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax1.set_xlabel('x(mm)')

ax2 = ax1.twinx()

ax2.plot(x,CH2O,label='CH2O',marker=">",markevery=markpts)
ax2.plot(x,CO,label='CO',marker="p",markevery=markpts)
ax2.plot(x,HCN,label='HCN',marker="<",markevery=markpts)
ax2.plot(x,N2O,label='N2O',marker="D",markevery=markpts)
ax2.plot(x,NO,label='NO',marker="^",markevery=markpts)
ax2.plot(x,NO2,label='NO2',marker="o",markevery=markpts)
ax2.plot(x,H2O,label='H2O',marker="s",markevery=markpts)
ax2.plot(x,CO2,label='CO2',marker="d",markevery=markpts)
#ax2.plot(x,c_HONO,label='cis-HONO',marker="s",mfc='None')
#ax2.plot(x,t_HONO,label='trans-HONO',marker="d",mfc='None')
ax2.plot(x,N2,label='N2',marker="<",mfc='None',markevery=markpts)
ax2.plot(x,ONNO2,label='ONNO2',marker="D",mfc='None',markevery=markpts)


#ax2.ticklabel_format(style='sci',axis='y',scilimits=(-1,1))
ax2.legend(loc='upper center',bbox_to_anchor=(0.5,1.2),ncol=5)
ax2.set_ylabel('Mass fraction of species in liquid')

ax1.xaxis.set_ticks_position('both')
#ax1.yaxis.set_ticks_position('both')
ax1.tick_params(direction='in',which='both')
ax2.tick_params(direction='in',which='both')

plt.savefig('Liquids.pdf',bbox_inches='tight')

