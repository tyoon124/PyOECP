#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 23:41:29 2021
"""

from PyOECP import Models

import matplotlib.pyplot as plt
import numpy as np

''' Example 2 Model fit
This script tries to fit a model using the MCMC method.
'''

np.random.seed(426082)

''' 2.1 Havriliak-Negami model '''

''' Step 1: Make a synthetic data. '''

frequency = np.logspace(np.log10(1e8),np.log10(10e10),201)

Synthetic = Models.Parameters()
Synthetic.Set('ei',2.0)
Synthetic.Set('magnitudes',50.0)

relaxation = 1/(2*np.pi*5e9)

Synthetic.Set('times',relaxation)
Synthetic.Set('As',0.15)
Synthetic.Set('Bs',0.08)

par = Synthetic.Parameters()
epsilon = Models.Discrete(frequency,par)

spacing = 10

''' Step 2: Fit a Discrete Model '''

Estimation = Models.Parameters()
ei = 2.5
Estimation.Set('ei',ei)
Estimation.Set('magnitudes',np.real(epsilon)[0]-ei)

MaximumFrequency = frequency[np.imag(epsilon)==np.min(np.imag(epsilon))]
relaxation = 1/(2*np.pi*MaximumFrequency)

Estimation.Set('times',relaxation)
Estimation.Set('As',0.3)
Estimation.Set('Bs',0.12)

InitialPar = Estimation.Parameters()

production = 100000
burnin = 5000

Trial = Models.MCMC(frequency,epsilon,InitialPar,production,
                       lb=None,ub=None,control=None,burnin=burnin,
                       Rate=0.01)
chi2s, Chains, par2 = Trial.Run()

FitResults = Models.Discrete(frequency,par2)

plt.figure(figsize=(5,4),dpi=300)
plt.gcf().subplots_adjust(bottom=0.15,left=0.15)

plt.semilogx(frequency,np.real(FitResults),color='r',label="$\epsilon'$ (Discrete Model)")
plt.semilogx(frequency,-np.imag(FitResults),color='b',label="$\epsilon''$ (Discrete Model)")
plt.plot(frequency[::spacing],np.real(epsilon)[::spacing],'o',markeredgecolor='r',
         markerfacecolor='None',markersize=8,label="$\epsilon'$ (Synthetic)")
plt.plot(frequency[::spacing],-np.imag(epsilon)[::spacing],'o',markeredgecolor='b',
         markerfacecolor='None',markersize=8,label="$\epsilon''$ (Synthetic)")
plt.xlabel('frequency [Hz]')
plt.ylabel(r'$\epsilon$')

plt.legend(loc='upper right', ncol=1, fontsize='xx-small',edgecolor='k')

plt.show()

''' Quantify the parameter uncertainty. '''
Trial.Examine(Chains,par)

print('------------------- Parameters and Errors -------------------')
for ele in Chains.keys():
    if len(Chains[ele]) != 0 or not all(np.isnan(Chains[ele])):
        if Chains[ele].shape[1] == 1:
            print(ele,':',np.mean(Chains[ele]),np.std(Chains[ele]))
        else:
            for Sub in range(Chains[ele].shape[1]):
                print(ele,':',np.mean(Chains[ele][:,Sub]),
                      np.std(Chains[ele][:,Sub]))


''' 2.2 Two-Debye relaxation model. '''

"""" Step 1: Make a synthetic data. """

frequency = np.logspace(np.log10(1e7),np.log10(10e10),201)

Synthetic = Models.Parameters()
Synthetic.Set('ei',1.3)
Synthetic.Set('magnitudes',[50.0,10.0])

relaxation = [1/(2*np.pi*4e10),1/(2*np.pi*2e10)]

Synthetic.Set('times',relaxation)

par = Synthetic.Parameters()
epsilon = Models.Discrete(frequency,par)


spacing = 10

""" Step 2: Fit a Discrete Model """

Estimation = Models.Parameters()
ei = 1.5
Estimation.Set('ei',ei)
Estimation.Set('magnitudes',[55.0,5.0])

relaxation = [1/(2*np.pi*8e10),1/(2*np.pi*7e9)]

Estimation.Set('times',relaxation)

InitialPar = Estimation.Parameters()


production = 100000

Trial = Models.MCMC(frequency,epsilon,InitialPar,production,
                       lb=None,ub=None,control=None,burnin=None,
                       Rate=0.01)
chi2s, Chains, par2 = Trial.Run()   

FitResults = Models.Discrete(frequency,par2)

plt.figure(figsize=(5,4),dpi=300)
plt.gcf().subplots_adjust(bottom=0.15,left=0.15)

plt.semilogx(frequency,np.real(FitResults),color='r',label="$\epsilon'$ (Discrete Model)")
plt.semilogx(frequency,-np.imag(FitResults),color='b',label="$\epsilon''$ (Discrete Model)")
plt.plot(frequency[::spacing],np.real(epsilon)[::spacing],'o',markeredgecolor='r',
         markerfacecolor='None',markersize=8,label="$\epsilon'$ (Synthetic)")
plt.plot(frequency[::spacing],-np.imag(epsilon)[::spacing],'o',markeredgecolor='b',
         markerfacecolor='None',markersize=8,label="$\epsilon''$ (Synthetic)")
plt.xlabel('frequency [Hz]')
plt.ylabel(r'$\epsilon$')

plt.legend(loc='center left', ncol=1, fontsize='xx-small',edgecolor='k')

Trial.Examine(Chains,par)

print('------------------- Parameters and Errors -------------------')
for ele in Chains.keys():
    try:
        if len(Chains[ele]) != 0 or not all(np.isnan(Chains[ele])):
            if Chains[ele].shape[1] == 1:
                print(ele,':',np.mean(Chains[ele]),np.std(Chains[ele]))
            else:
                for Sub in range(Chains[ele].shape[1]):
                    print(ele,':',np.mean(Chains[ele][:,Sub]),
                          np.std(Chains[ele][:,Sub]))
    except TypeError:
        continue
