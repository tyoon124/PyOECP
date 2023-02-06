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

TransformModel = Transform.Marsland(frequency,S11m,S11r0,S11r1,S11r2,S11r3,
                                   m2='Open',m3='Water_Kaatze',m4='Acetone_Onimisi',temperature=T,
                                   Window=81,concentrations=[None,None,None,None])

MarslandE = TransformModel.Calculate()

spacing = 10

TransformModel1 = Transform.Stuchly(frequency,S11m,S11r0,S11r1,S11r2,
                                        m1='Short',m2='Open',m3='Water_Kaatze',Window=51)

StuchlyE = TransformModel1.Calculate()

Komarov = Transform.Komarov(frequency, S11m, S11r1, S11r2, S11r3,
                            'Open','Water_Kaatze','Acetone_Onimisi',
                            1,3.8,2.1+0*1j,M=50,Window=51)

KomarovE = Komarov.epsilon

fig, (ax1,ax2) = plt.subplots(2,1)
fig.set_size_inches((5,8))
fig.set_dpi(300)
font = {'size':15}
plt.rc('font', **font)
plt.rcParams['font.family'] = 'serif'

'''Let's visualize the data.'''
ax1.semilogx(frequency[::spacing],np.real(MarslandE)[::spacing],'o',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.0,markersize=7,label="$\epsilon'$ (Marsland)")
ax1.semilogx(frequency[::spacing],-np.imag(MarslandE)[::spacing],'o',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.0,markersize=7,label="$\epsilon''$ (Marsland)")

ax1.semilogx(frequency[::spacing],np.real(StuchlyE)[::spacing],'s',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.0,markersize=7,label="$\epsilon'$ (Stuchly)")
ax1.semilogx(frequency[::spacing],-np.imag(StuchlyE)[::spacing],'s',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.0,markersize=7,label="$\epsilon'$ (Stuchly)")

ax1.semilogx(frequency[::spacing],np.real(KomarovE)[::spacing],'^',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.0,markersize=7,label="$\epsilon'$ (Komarov)")
ax1.semilogx(frequency[::spacing],-np.imag(KomarovE)[::spacing],'^',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.0,markersize=7,label="$\epsilon'$ (Komarov)")

Theoretical = References.Methanol_Barthel(frequency,temperature=T)['epsilon']
ax1.semilogx(frequency,np.real(Theoretical),color='red',linewidth=1.0,label="$\epsilon'$ (Literature)")
ax1.semilogx(frequency,-np.imag(Theoretical),'--',color='blue',linewidth=1.0,label="$\epsilon''$ (Literature)")

ax1.set_ylabel("$\epsilon$")
ax1.set_ylim([0,50])
ax1.legend(loc='upper right', ncol=2, fontsize='xx-small',edgecolor='k')
ax1.text(-0.25,1,'(a)',transform=ax1.transAxes)

''' 1.2 High Frequency Data '''
address = 'data/high/'
S11r0 = References.Parser(address + 'S11Short.csv')
S11r1 = References.Parser(address + 'S11Open.csv')
S11r2 = References.Parser(address + 'S11Water.csv')
S11r3 = References.Parser(address + 'S11Acetone.csv')
S11m = References.Parser(address + 'S11Methanol.csv')

frequency = S11r1[:,0]

TransformModel = Transform.Marsland(frequency,S11m,S11r0,S11r1,S11r2,S11r3,
                                   m2='Open',m3='Water_Kaatze',m4='Acetone_Onimisi',temperature=T,
                                   Window=101,concentrations=[None,None,None,None])

MarslandE = TransformModel.Calculate()

spacing = 10

TransformModel1 = Transform.Stuchly(frequency,S11m,S11r0,S11r1,S11r2,
                                        m1='Short',m2='Open',m3='Water_Kaatze',Window=51)

StuchlyE = TransformModel1.Calculate()

Komarov = Transform.Komarov(frequency, S11m, S11r1, S11r3, S11r2,
                            'Open','Acetone_Onimisi','Water_Kaatze',
                            0.3,0.8,2.1+0*1j,M=50,Window=51)

KomarovE = Komarov.epsilon

Theoretical = References.Methanol_Barthel(frequency,temperature=T)['epsilon']

'''Let's visualize the data.'''
ax2.semilogx(frequency[::spacing],np.real(MarslandE)[::spacing],'o',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.0,markersize=7,label="$\epsilon'$ (Marsland)")
ax2.semilogx(frequency[::spacing],-np.imag(MarslandE)[::spacing],'o',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.0,markersize=7,label="$\epsilon''$ (Marsland)")

ax2.semilogx(frequency[::spacing],np.real(StuchlyE)[::spacing],'s',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.0,markersize=7,label="$\epsilon'$ (Stuchly)")
ax2.semilogx(frequency[::spacing],-np.imag(StuchlyE)[::spacing],'s',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.0,markersize=7,label="$\epsilon'$ (Stuchly)")

ax2.semilogx(frequency[::spacing],np.real(KomarovE)[::spacing],'^',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.0,markersize=7,label="$\epsilon'$ (Komarov)")
ax2.semilogx(frequency[::spacing],-np.imag(KomarovE)[::spacing],'^',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.0,markersize=7,label="$\epsilon'$ (Komarov)")

ax2.semilogx(frequency,np.real(Theoretical),color='red',linewidth=1.0,label="$\epsilon'$ (Literature)")
ax2.semilogx(frequency,-np.imag(Theoretical),'--',color='blue',linewidth=1.0,label="$\epsilon''$ (Literature)")

ax2.set_xlabel("frequency [Hz]")
ax2.set_ylabel("$\epsilon$")
ax2.set_ylim([0,50])
ax2.legend(loc='upper right', ncol=2, fontsize='xx-small',edgecolor='k')
ax2.text(-0.25,1,'(b)',transform=ax2.transAxes)

plt.savefig('Figure3.pdf',dpi=300)

