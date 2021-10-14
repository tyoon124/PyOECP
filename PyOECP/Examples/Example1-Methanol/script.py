#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 23:41:29 2021
"""

from PyOECP import References
from PyOECP import Transform
import matplotlib.pyplot as plt
import numpy as np

''' Example 1 Methanol
This script tries to convert the reflection coefficients from VNAs.
The data files and VNAs are as follows.
VNA 
Short: LS11Short.csv
Open: LS11Open.csv
Acetone: LS11Acetone.csv
Water: LS11Water.csv
Methanol: LS11Methanol.csv

'''

''' 1.1 Low Frequency Data '''
T = 25

address = 'data/low/'
S11r0 = References.Parser(address + 'S11Short.csv')
S11r1 = References.Parser(address + 'S11Open.csv')
S11r2 = References.Parser(address + 'S11Water.csv')
S11r3 = References.Parser(address + 'S11Acetone.csv')
S11m = References.Parser(address + 'S11Methanol.csv')

frequency = S11r1[:,0]

TransformModel = Transform.Antenna(frequency,S11m,S11r0,S11r1,S11r2,S11r3,
                                   m2='Open',m3='Water_Kaatze',m4='Acetone_Onimisi',temperature=T,
                                   Window=101,concentrations=[None,None,None,None])

AntennaE = TransformModel.Calculate()

spacing = 10

TransformModel1 = Transform.Capacitance(frequency,S11m,S11r0,S11r1,S11r2,
                                        m1='Short',m2='Open',m3='Water_Kaatze',Window=91)

CapacitanceE = TransformModel1.Calculate()


plt.figure(figsize=(5,4),dpi=300)
plt.gcf().subplots_adjust(bottom=0.15,left=0.15)

font = {'size':15}
plt.rc('font', **font)
plt.rcParams['font.family'] = 'serif'

"""Let's visualize the data."""
plt.semilogx(frequency[::spacing],np.real(AntennaE)[::spacing],'o',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.5,markersize=8,label="$\epsilon'$ (Antenna)")
plt.semilogx(frequency[::spacing],-np.imag(AntennaE)[::spacing],'o',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.5,markersize=8,label="$\epsilon''$ (Antenna)")

plt.semilogx(frequency[::spacing],np.real(CapacitanceE)[::spacing],'s',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.5,markersize=8,label="$\epsilon'$ (Capacitance)")
plt.semilogx(frequency[::spacing],-np.imag(CapacitanceE)[::spacing],'s',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.5,markersize=8,label="$\epsilon'$ (Capacitance)")

Theoretical = References.Methanol_Barthel(frequency,temperature=T)['epsilon']
plt.semilogx(frequency,np.real(Theoretical),color='red',linewidth=2,label="$\epsilon'$ (Literature)")
plt.semilogx(frequency,-np.imag(Theoretical),'--',color='blue',linewidth=2,label="$\epsilon''$ (Literature)")

plt.xlabel("frequency [Hz]")
plt.ylabel("$\epsilon$")
plt.ylim([0,45])
plt.legend(loc='upper right', ncol=2, fontsize='xx-small',edgecolor='k')
plt.show()

''' 1.2 High Frequency Data '''
address = 'data/high/'
S11r0 = References.Parser(address + 'S11Short.csv')
S11r1 = References.Parser(address + 'S11Open.csv')
S11r2 = References.Parser(address + 'S11Water.csv')
S11r3 = References.Parser(address + 'S11Acetone.csv')
S11m = References.Parser(address + 'S11Methanol.csv')

frequency = S11r1[:,0]

TransformModel = Transform.Antenna(frequency,S11m,S11r0,S11r1,S11r2,S11r3,
                                   m2='Open',m3='Water_Kaatze',m4='Acetone_Onimisi',temperature=T,
                                   Window=101,concentrations=[None,None,None,None])

AntennaE = TransformModel.Calculate()

spacing = 10

TransformModel1 = Transform.Capacitance(frequency,S11m,S11r0,S11r1,S11r2,
                                        m1='Short',m2='Open',m3='Water_Kaatze',Window=91)

CapacitanceE = TransformModel1.Calculate()


plt.figure(figsize=(5,4),dpi=300)
plt.gcf().subplots_adjust(bottom=0.15,left=0.15)

font = {'size':15}
plt.rc('font', **font)
plt.rcParams['font.family'] = 'serif'

"""Let's visualize the data."""
plt.semilogx(frequency[::spacing],np.real(AntennaE)[::spacing],'o',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.5,markersize=8,label="$\epsilon'$ (Antenna)")
plt.semilogx(frequency[::spacing],-np.imag(AntennaE)[::spacing],'o',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.5,markersize=8,label="$\epsilon''$ (Antenna)")

plt.semilogx(frequency[::spacing],np.real(CapacitanceE)[::spacing],'s',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.5,markersize=8,label="$\epsilon'$ (Capacitance)")
plt.semilogx(frequency[::spacing],-np.imag(CapacitanceE)[::spacing],'s',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.5,markersize=8,label="$\epsilon'$ (Capacitance)")

Theoretical = References.Methanol_Barthel(frequency,temperature=T)['epsilon']
plt.semilogx(frequency,np.real(Theoretical),color='red',linewidth=2,label="$\epsilon'$ (Literature)")
plt.semilogx(frequency,-np.imag(Theoretical),'--',color='blue',linewidth=2,label="$\epsilon''$ (Literature)")

plt.xlabel("frequency [Hz]")
plt.ylabel("$\epsilon$")
plt.ylim([0,45])
plt.legend(loc='upper right', ncol=2, fontsize='xx-small',edgecolor='k')
plt.show()