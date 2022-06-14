# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 22:56:38 2021

@author: tae-jun_yoon
"""

from PyOECP import Models

import matplotlib.pyplot as plt
import numpy as np

''' Example 4: Continuous relaxation model. 
We are going to use a synthetic data prepared in Example 2.
'''

font = {'size':15}
plt.rc('font', **font)
plt.rcParams['font.family'] = 'serif'

''' 4.1 Synthetic Data: Two-Debye model '''

frequency = np.logspace(np.log10(1e7),np.log10(10e10),201)

Synthetic = Models.Parameters()
Synthetic.Set('ei',1.3)
Synthetic.Set('magnitudes',[50.0,10.0])

relaxation = [1/(2*np.pi*4e10),1/(2*np.pi*2e10)]

Synthetic.Set('times',relaxation)

par = Synthetic.Parameters()
epsilon = Models.Discrete(frequency,par)

span = np.logspace(-12,-9,801)
Probability, yfit, res = Models.Continuous(frequency,epsilon - par['ei'],span,alpha=0.05)

spacing = 10

plt.figure(figsize=(5,4),dpi=300)
plt.gcf().subplots_adjust(bottom=0.15,left=0.15)

plt.semilogx(frequency,np.real(yfit)+par['ei'],'r',label="$\epsilon'$ (Continuous)")
plt.semilogx(frequency[::spacing],np.real(epsilon)[::spacing],'o',
             markerfacecolor='None',markeredgecolor='red',
             markeredgewidth=1.5,markersize=8,label="$\epsilon'$ (Synthetic)")
plt.semilogx(frequency,-np.imag(yfit),'b',label="$\epsilon''$ (Continuous)")
plt.semilogx(frequency[::spacing],-np.imag(epsilon)[::spacing],'o',
             markerfacecolor='None',markeredgecolor='blue',
             markeredgewidth=1.5,markersize=8,label="$\epsilon''$ (Synthetic)")

plt.xlabel("frequency [Hz]")
plt.ylabel("$\epsilon$")
plt.ylim([0,80])
plt.legend(loc='upper center', ncol=2, fontsize='xx-small',edgecolor='k')
plt.show()

plt.figure(figsize=(5,4),dpi=300)
plt.gcf().subplots_adjust(bottom=0.15,left=0.15)

plt.plot([relaxation[0]*1e12,relaxation[0]*1e12],[-1,10],'r--',label="$\\tau_{1}$ (Synthetic)")
plt.plot([relaxation[1]*1e12,relaxation[1]*1e12],[-1,10],'b--',label="$\\tau_{2}$ (Synthetic)")
plt.semilogx(span/1e-12,Probability,'k',label="Continuous model")
plt.ylabel(r"$g(\tau)$")
plt.xlabel(r"$\tau$ [ps]")
plt.xlim([2,1e2])
plt.ylim([-0.25,8])
plt.legend(loc='upper right', ncol=1, fontsize='xx-small',edgecolor='k')
plt.show()