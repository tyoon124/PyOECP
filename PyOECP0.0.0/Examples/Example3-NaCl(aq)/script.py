#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 23:41:29 2021
"""

from PyOECP import Models
from PyOECP import References
from PyOECP import Transform

import matplotlib.pyplot as plt
import numpy as np

np.random.seed(426082)

''' Example 3 Aqueous sodium chloride solution
The goals of this script is as follows.

(1) We will see if the selection of reference liquids can have a significant impact on the accuracy.
(2) We will check if the implemention of the EP correction algorithm is okay.

For goal (1), we will do the following things.
(a) Convert the reflection coefficients in NaCl (aq) 0.18 M solution from a VNA. 
    We will use (open, short, water, and acetone) in the first set.
(b) Convert the reflection coefficients in NaCl (aq) 0.18 M solution from a VNA.
    We will use (open, short, water, and NaCl (aq) 0.09 M) in the second set.
(c) Compare the results from (a) and (b).

For goal (2), we will do the following things.
(a) Use NaCl (aq) 0.09 M as a reference liquid. Peyman model is used.
(3) Subtract the specific conductance contribution and EP contribution from NaCl (aq) 0.18 M.

The data files and VNAs are as follows.
Short: S11Short.csv
Open: S11Open.csv
Acetone: S11Acetone.csv
Water: S11Water.csv
NaCl (aq) 0.09 M: S11NaClL.csv
NaCl (aq) 0.18 M: S11NaClH.csv
'''

T = 25

''' 3.1. Let's examine the influence of references. '''
address = 'data/'
S11r0 = References.Parser(address + 'S11Short.csv')
S11r1 = References.Parser(address + 'S11Open.csv')
S11r21 = References.Parser(address + 'S11NaClL.csv') 
S11r22 = References.Parser(address + 'S11Water.csv') 
S11r3 = References.Parser(address + 'S11Acetone.csv')

S11m = References.Parser(address + 'S11NaClH.csv')

frequency = S11r1[:,0]

TransformModel = Transform.Antenna(frequency,S11m,S11r0,S11r1,S11r21,S11r3,
                                   m2='Open',m3='NaClAqueous_Peyman',m4='Acetone_Wei',temperature=T,
                                   Window=101,concentrations=[None,None,0.09,None])
AntennaData1 = TransformModel.Calculate()

plt.figure(figsize=(5,4),dpi=300)
plt.gcf().subplots_adjust(bottom=0.15,left=0.15)

font = {'size':15}
plt.rc('font', **font)
plt.rcParams['font.family'] = 'serif'

''' Let's visualize the data. '''
spacing = 10
plt.semilogx(frequency[::spacing],np.real(AntennaData1)[::spacing],'o',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.5,markersize=8,label="$\epsilon'$ (NaCl (aq) ref.)")
plt.semilogx(frequency[::spacing],-np.imag(AntennaData1)[::spacing],'o',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.5,markersize=8,label="$\epsilon''$ (NaCl (aq) ref.)")

TransformModel = Transform.Antenna(frequency,S11m,S11r0,S11r1,S11r22,S11r3,
                                   m2='Open',m3='Water_Kaatze',m4='Acetone_Wei',temperature=T,
                                   Window=101,concentrations=[None,None,None,None])
AntennaData2 = TransformModel.Calculate()

''' Let's visualize the data. '''
spacing = 10
plt.semilogx(frequency[::spacing],np.real(AntennaData2)[::spacing],'s',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.5,markersize=8,label="$\epsilon'$ (Water ref.)")
plt.semilogx(frequency[::spacing],-np.imag(AntennaData2)[::spacing],'s',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.5,markersize=8,label="$\epsilon''$ (Water ref.)")

Peyman = References.NaClAqueous_Peyman(frequency,c=0.18)['epsilon']
plt.semilogx(frequency,np.real(Peyman),'r',linewidth=2,label="$\epsilon'$ (Peyman, c=0.18 M)")
plt.semilogx(frequency,-np.imag(Peyman),'b',linewidth=2, label="$\epsilon''$ (Peyman, c=0.18 M)")

plt.xlabel("frequency [Hz]")
plt.ylabel("$\epsilon$")
plt.legend(loc='upper right', ncol=1, fontsize='xx-small',edgecolor='k')

plt.show()

''' 3.2 EP Correction '''
''' 3.2.1. Load data '''
address = 'data/'
S11r0 = References.Parser(address + 'S11Short.csv')
S11r1 = References.Parser(address + 'S11Open.csv')
S11r2 = References.Parser(address + 'S11Water.csv') 
S11r3 = References.Parser(address + 'S11Acetone.csv')

S11m = References.Parser(address + 'S11NaClH.csv')


frequency = S11r1[:,0]

TransformModel = Transform.Antenna(frequency,S11m,S11r0,S11r1,S11r2,S11r3,
                                   m2='Open',m3='Water_Kaatze',m4='Acetone_Wei',temperature=T,
                                   Window=101,concentrations=[None,None,None,None])
AntennaData = TransformModel.Calculate()

''' 3.2.1. Fit a model. 
We will fit a Cole-Cole relaxation model in conjunction with the electrode polarization model.'''

par = Models.Parameters()
''''3.2.2. Initial estimation of the conductance.'''
relaxation = -np.imag(AntennaData)
relaxation = relaxation[frequency<5e8]
frequencyL = frequency[frequency<5e8]

conductance = lambda x, a: x/(2*np.pi*a*8.8541878128-12)
from scipy.optimize import curve_fit
popt = curve_fit(conductance,frequencyL,relaxation)
conductance = popt[0]

''' 3.2.3. Set the initial model parameters. We will use the Cole-Cole model. '''
par.Set('ei',1.0)
par.Set('conductance',conductance)
par.Set('As',[0.2])
par.Set('magnitudes',np.max(np.real(AntennaData))-1.0)
par.Set('times',1/(2*np.pi*frequency[np.imag(AntennaData)==np.min(np.imag(AntennaData))]))
par.Set('CPEs',[1.0,1.0])
par = par.Parameters()

production = 30000
Trial = Models.MCMC(frequency,AntennaData,par,production,
                    lb=None,ub=None,control=None,burnin=None,Rate=0.010)
chi2s, Chain, par2 = Trial.Run()

FittedData = Models.Discrete(frequency,par2)

''' 3.2.4. Compare the result by excluding the electrode polarization effect and conductance contribution. '''

par3 = par2
par3['CPEs'] = np.array([],ndmin=1)
par3['conductance'] = np.array([],ndmin=1)
Data = Models.Discrete(frequency,par3)

plt.figure(figsize=(5,4),dpi=300)
plt.gcf().subplots_adjust(bottom=0.15,left=0.15)

spacing = 10
plt.semilogx(frequency[::spacing],np.real(AntennaData)[::spacing],'o',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.5,markersize=8,label="$\epsilon'$ (Antenna)")
plt.semilogx(frequency[::spacing],-np.imag(AntennaData)[::spacing],'o',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.5,markersize=8,label="$\epsilon''$ (Antenna)")
plt.semilogx(frequency[::spacing],np.real(Data)[::spacing],'s',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.5,markersize=8,label="$\epsilon'$ (EP & $\kappa$ correction)")
plt.semilogx(frequency[::spacing],-np.imag(Data)[::spacing],'s',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.5,markersize=8,label="$\epsilon''$ (EP & $\kappa$ correction)")

plt.xlabel("frequency [Hz]")
plt.ylabel("$\epsilon$")
plt.legend(loc='upper right', ncol=1, fontsize='xx-small',edgecolor='k')

r1 = np.real(AntennaData)
j1 = -np.imag(AntennaData)
r2 = np.real(Data)
j2 = -np.imag(Data)