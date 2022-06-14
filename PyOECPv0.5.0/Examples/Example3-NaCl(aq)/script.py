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
import copy

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
(a) Use NaCl (aq) 1.44 M as a reference liquid. Peyman model is used.
(3) Subtract the specific conductance contribution and EP contribution from NaCl (aq) 0.18 M.

The data files and VNAs are as follows.
Short: S11Short.csv
Open: S11Open.csv
Acetone: S11Acetone.csv
Water: S11Water.csv
NaCl (aq) 0.09 M: S11NaClL1.csv
NaCl (aq) 0.18 M: S11NaClL2.csv
NaCl (aq) 1.44 M: S11NaClH.csv
'''

np.random.seed(881346)

T = 25
epsilon0 = 8.8541878128e-12

''' 3.1. Let's examine the influence of references. '''
address = 'data/'
S11r0 = References.Parser(address + 'S11ShortL.csv')
S11r1 = References.Parser(address + 'S11OpenL.csv')
S11r21 = References.Parser(address + 'S11NaClL1.csv') 
S11r22 = References.Parser(address + 'S11WaterL.csv') 
S11r3 = References.Parser(address + 'S11AcetoneL.csv')

S11m = References.Parser(address + 'S11NaClL2.csv')

frequency = S11r1[:,0]

TransformModel = Transform.Marsland(frequency,S11m,S11r0,S11r1,S11r21,S11r3,
                                   m2='Open',m3='NaClAqueous_Peyman',m4='Acetone_Wei',temperature=T,
                                   Window=101,concentrations=[None,None,0.09,None])
MarslandData1 = TransformModel.Calculate()

TransformModel = Transform.Marsland(frequency,S11m,S11r0,S11r1,S11r22,S11r3,
                                   m2='Open',m3='Water_Kaatze',m4='Acetone_Wei',temperature=T,
                                   Window=101,concentrations=[None,None,None,None])
MarslandData2 = TransformModel.Calculate()

''' Check if the Komarov model yields the similar result. '''
TransformModel = Transform.Komarov(frequency, S11m, S11r1, S11r3, S11r21,
                                   'Open','Acetone_Wei','NaClAqueous_Peyman',
                                   0.3,0.8,2.1+0*1j,M=50,Window=51,
                                   concentrations=[None,None,0.09])

KomarovData1 = TransformModel.epsilon

TransformModel = Transform.Komarov(frequency, S11m, S11r1, S11r3, S11r22,
                                   'Open','Acetone_Wei','Water_Kaatze',
                                   0.3,0.8,2.1+0*1j,M=50,Window=51,
                                   concentrations=[None,None,None])

KomarovData2 = TransformModel.epsilon

# Reference data generation.
Peyman = References.NaClAqueous_Peyman(frequency,c=0.18)['epsilon']

''' Let's visualize the data. '''
fig, (ax1,ax2) = plt.subplots(2,1)
fig.set_size_inches((5,8))
fig.set_dpi(300)
font = {'size':15}
plt.rc('font', **font)
plt.rcParams['font.family'] = 'serif'

spacing = 10

ax1.semilogx(frequency[::spacing],np.real(MarslandData1)[::spacing],'o',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.0,markersize=7,label="$\epsilon'$ (Marsland, NaCl ref.)")
ax1.semilogx(frequency[::spacing],-np.imag(MarslandData1)[::spacing],'o',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.0,markersize=7,label="$\epsilon''$ (Marsland, NaCl ref.)")
ax1.semilogx(frequency[::spacing],np.real(MarslandData2)[::spacing],'s',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.0,markersize=7,label="$\epsilon'$ (Marsland, Water ref.)")
ax1.semilogx(frequency[::spacing],-np.imag(MarslandData2)[::spacing],'s',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.0,markersize=7,label="$\epsilon''$ (Marsland, Water ref.)")
ax1.semilogx(frequency,np.real(Peyman),'r',linewidth=2,label="$\epsilon'$ (Peyman, c=0.18 M)")

ax1.semilogx(frequency[::spacing],np.real(KomarovData1)[::spacing],'>',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.0,markersize=7,label="$\epsilon'$ (Komarov, NaCl ref.)")
ax1.semilogx(frequency[::spacing],-np.imag(KomarovData1)[::spacing],'>',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.0,markersize=7,label="$\epsilon''$ (Komarov, NaCl ref.)")
ax1.semilogx(frequency[::spacing],np.real(KomarovData2)[::spacing],'p',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.0,markersize=7,label="$\epsilon'$ (Komarov, Water ref.)")
ax1.semilogx(frequency[::spacing],-np.imag(KomarovData2)[::spacing],'p',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.0,markersize=7,label="$\epsilon''$ (Komarov, Water ref.)")
ax1.semilogx(frequency,-np.imag(Peyman),'b',linewidth=2, label="$\epsilon''$ (Peyman, c=0.18 M)")

#ax1.xlabel("frequency [Hz]")
ax1.set_ylabel("$\epsilon$")
ax1.set_ylim([0,160])
ax1.legend(loc='upper right', ncol=2, fontsize=7,edgecolor='k')
ax1.text(-0.25,1,'(a)',transform=ax1.transAxes)
#ax1.show()
#plt.savefig('Figure6.pdf',dpi=300)


''' 3.2. High-concentration NaCl (aq)
Let's examine the influence of references. '''
address = 'data/'
S11r0 = References.Parser(address + 'S11ShortH.csv')
S11r1 = References.Parser(address + 'S11OpenH.csv')
S11r2 = References.Parser(address + 'S11WaterH.csv') 
S11r3 = References.Parser(address + 'S11AcetoneH.csv')

S11m = References.Parser(address + 'S11NaClH.csv')

frequency = S11r1[:,0]

TransformModel = Transform.Komarov(frequency, S11m, S11r0, S11r1, S11r2, 
                                   'Short','Open','Water_Kaatze',
                                   0.3,0.8,2.1+0*1j,M=50,Window=81,
                                   concentrations=[None,None,None],temperature=T)

KomarovData = TransformModel.epsilon

Peyman = References.NaClAqueous_Peyman(frequency,c=1.44,temperature=T)['epsilon']

''' Let's visualize the data in terms of the permittivity. '''

spacing = 8

par = Models.Parameters()
''''3.2.2. Initial estimation of the conductance.'''
relaxation = -np.imag(KomarovData)
relaxation = relaxation[frequency<5e8]
frequencyL = frequency[frequency<5e8]

conductance = lambda x, a: x/(2*np.pi*a*epsilon0)
from scipy.optimize import curve_fit
popt = curve_fit(conductance,frequencyL,relaxation)
conductance = popt[0]

''' 3.2.3. Set the initial model parameters. We will use the Cole-Cole model. '''
par.Set('ei',1.0)
par.Set('conductance',conductance)
par.Set('As',[0.6])
par.Set('magnitudes',70.0)
par.Set('times',1/(2*np.pi*2e10))
par.Set('CPEs',[400,.9])
lb = copy.deepcopy(par)
ub = copy.deepcopy(par)
par = par.Parameters()
lb.Set('conductance',par['conductance']*0.5)
ub.Set('conductance',par['conductance']*1.5)
lb.Set('magnitudes',30.0)
ub.Set('magnitudes',90.0)
lb.Set('times',1/(2*np.pi*1e10))
ub.Set('times',1/(2*np.pi*3e10))
lb.Set('As',0.0)
ub.Set('As',1.0)
lb.Set('CPEs',[0.0,0.0])
ub.Set('CPEs',[np.inf,1.0])

production = 30000

lb = lb.Parameters()
ub = ub.Parameters()

Names = lb.keys()
for Name in Names:
    if lb[Name] is None:
        continue
    else:
        for ele in range(len(lb[Name])):
            lb[Name][ele] = 0.0
            ub[Name][ele] = np.inf

Trial = Models.MCMC(frequency,KomarovData,par,production,
                    lb=lb,ub=ub,control=None,burnin=None,Rate=0.05)
chi2s, Chain, par2 = Trial.Run()

FittedData = Models.Discrete(frequency,par2)

''' 3.2.4. Compare the result by excluding the electrode polarization effect and conductance contribution. '''
par3 = par2
par3['CPEs'] = np.array([],ndmin=1)
par3['conductance'] = np.array([],ndmin=1)
Data = Models.Discrete(frequency,par3)

plt.figure(figsize=(5,4),dpi=300)
plt.gcf().subplots_adjust(bottom=0.15,left=0.15)

spacing1 = 8
spacing2 = 10

ax2.semilogx(frequency[::spacing1],np.real(KomarovData)[::spacing1],'o',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.0,markersize=7,label="$\epsilon'$ (Komarov)")
ax2.semilogx(frequency[::spacing1],-np.imag(KomarovData)[::spacing1],'o',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.0,markersize=7,label="$\epsilon''$ (Komarov)")

ax2.semilogx(frequency[::spacing2],np.real(Data)[::spacing2],'s',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.0,markersize=7,label="$\epsilon'$ (EP & $\kappa$ correction)")
ax2.semilogx(frequency[::spacing2],-np.imag(Data)[::spacing2],'s',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.0,markersize=7,label="$\epsilon''$ (EP & $\kappa$ correction)")

ax2.semilogx(frequency,np.real(Peyman),'r',linewidth=2,label="$\epsilon'$ (Peyman, c=1.44 M)")
ax2.semilogx(frequency,-np.imag(Peyman),'b',linewidth=2, label="$\epsilon''$ (Peyman, c=1.44 M)")

ax2.set_xlabel("frequency [Hz]")
ax2.set_ylabel("$\epsilon$")
ax2.set_ylim([0,300])
ax2.legend(loc='upper right', ncol=1, fontsize='xx-small',edgecolor='k')
ax2.text(-0.25,1,'(b)',transform=ax2.transAxes)
fig.savefig('Figure6.pdf',format='pdf',dpi=300,bbox_inches='tight')