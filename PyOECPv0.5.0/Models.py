#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 29 14:16:07 2021

@author: tae-jun_yoon
"""

'''
This script contains relaxation models.
'''

import numpy as np
import sys
import matplotlib.pyplot as plt
import copy
from tqdm import tqdm

class Parameters:
    
    def __init__(self):                
        self.par = {}
        self.par['ei'] = np.array([])
        self.par['magnitudes'] = np.array([])
        self.par['times'] = np.array([])
        self.par['As'] = np.array([])
        self.par['Bs'] = np.array([])
        self.par['CPEs'] = np.array([])
        self.par['conductance'] = np.array([])
    
    def Set(self,name,val):
        ''' Set the parameter values. '''
        try:
            self.par[name] = np.array(val,ndmin=1)
            return None
        except KeyError:
            print('The key variables are ei, magnitudes, times, As, Bs, CPEs, and conductance.')
    
    def Parameters(self):
        ''' Check if "par" is correctly set, and return "par" '''
        Names = ['As','Bs']
        for Name in Names:
            if len(self.par[Name]) < len(self.par['magnitudes']):                
                for i in range(len(self.par['magnitudes'])-len(self.par[Name])):
                    self.par[Name] = np.append(self.par[Name],None)       
            
        return self.par
        

def Discrete(frequency, par):
    ''' Multiple discrete relaxation model '''    
    
    Omega = 2*np.pi*frequency
    
    ei = par['ei']
    magnitudes = par['magnitudes']
    times = par['times']
    As = par['As']
    Bs = par['Bs']
    CPEs = par['CPEs']
    conductance = par['conductance']
    
    epsilon = ei*np.ones((len(frequency),),dtype=complex)
    
    e0 = 8.8541878128e-12
        
    for ele in range(len(magnitudes)):
        A = As[ele]
        B = Bs[ele]
        M = magnitudes[ele]
        T = times[ele]
        if A is not None and B is not None:
            # Havriliak-Negami Model
            epsilon += M/((1+1j*T*Omega)**(1-A))**(1-B)
        elif A is None and B is not None:
            # Cole-Davidson Model
            epsilon += M/((1+1j*T*Omega))**(1-B)
        elif A is not None and B is None: 
            # Cole-Cole Model
            epsilon += M/((1+1j*T*Omega)**(1-A))
        else:
            # Debye Model
            epsilon += M/((1+1j*T*Omega))
    
    if len(conductance) != 0:
        epsilon = epsilon + conductance/(1j*Omega*e0)
    
    if len(CPEs) != 0:        
        try:
            if len(CPEs) != 2:
                sys.exit('Wrong parameter set for CPE representation. (2 Parameters Required.)')
            else:                
                CPEs = CPEs[0]*(1j*Omega)**(1-CPEs[1])
                epsilon = epsilon/(1+epsilon/CPEs)
        except TypeError:
            print('Wrong parameter set for CPE representation. (2 Parameters Required.)')
    
    return epsilon

def Continuous(frequency,epsilon,parameters,alpha = 0.05):    
    ''' Construct continuous relaxation time distribution based on the method by Zasetsky-Buchner.
    Citation information 
    Zasetsky, A. Y., & Buchner, R. (2010). 
    Quasi-linear least squares and computer code for numerical evaluation of relaxation time distribution from broadband dielectric spectra. 
    Journal of Physics: Condensed Matter, 23(2), 025903.
    Note: The EP correction should be performed a priori to use this function for conductive (lossy) samples.
    '''
    from scipy.optimize import nnls
    
    b = np.concatenate((np.real(epsilon),np.imag(epsilon),np.zeros(len(parameters),)))
    X = np.zeros((len(b),len(parameters)))
    htau = np.log10(parameters[1]) - np.log10(parameters[0])
    for i in range(len(frequency)):
        y = 1/(1+1j*2*np.pi*frequency[i]*parameters)         
        X[i,:] = np.real(y)
        X[len(frequency)+i,:] = np.imag(y)
    for i in range(1,len(parameters)-1):
        X[len(frequency)*2+i,i] = -2*alpha
        X[len(frequency)*2+i,i-1] = 1*alpha
        X[len(frequency)*2+i,i+1] = 1*alpha    
    X[len(frequency)*2,0] = htau
    X[-1,-1] = htau
    p1, res = nnls(X,b)        
    
    yfit1 = np.matmul(X[:2*len(frequency),:],p1)
    yfit1 = yfit1[:len(frequency)] + 1j*yfit1[len(frequency):] 
    
    return p1, yfit1, res

class MCMC:
    ''' Perform Markov chain Monte Carlo simulation for estimating the parameters. '''
    def __init__(self,frequency,data,par,production,
                 lb=None,ub=None,control=None,burnin=None,Rate=0.10,
                 weight=None):
        ''' Input arguments
        data: Experimental data to fit.
        par: Model parameters. This variable should be generated from Parameters.Parameters().
        (Dictionary variable)
        production: The length of the production run.
        lb: lower bound. If not given, -infinity is applied.
        ub: upper bound. If not given, infinity is applied.
        control: Determine variables to be changed. If one of the variables is set to be true, the corresponding variable changes.        
        burnin: The length of the burnin run. If not given, 1/4 of the production run.        
        Accept: acceptance rate when rejected.
        Rate: the magnitude of change (Gaussian width)
        Weight: the uncertainty of the measured data can be used as a weight.
        '''
        
        self.frequency = frequency
        self.data = data        
        self.par = copy.deepcopy(par)
        self.production = production
        if burnin is None:
            burnin = int(production/4)
        self.burnin = burnin
        self.weight = weight
        self.Names = list(par.keys())
        
        if lb is None:            
            lb = copy.deepcopy(par)
            for Name in self.Names:
                if lb[Name] is None:
                    continue
                else:
                    for ele in range(len(lb[Name])):
                        lb[Name][ele] = 0.0
        
        lb['ei'] = np.array([1.0])
        
        if ub is None:            
            ub = copy.deepcopy(par)
            for Name in self.Names:
                if ub[Name] is None:
                    continue
                else:                    
                    for ele in range(len(ub[Name])):
                        ub[Name][ele] = np.inf
        
        for ind in range(len(ub['magnitudes'])):
            ub['magnitudes'][ind] = np.max(np.real(data)) - 1
            
        if control is None:
            control = copy.deepcopy(par)
            for Name in self.Names:
                if len(control[Name]) == 0:
                    control[Name] = False
                else:                    
                    for ele in range(len(control[Name])):
                        if self.par[Name][ele] is None:
                            control[Name][ele] = False
                        else:
                            control[Name][ele] = True
                    
        self.lb = lb
        self.ub = ub
        self.control = control        
        
        self.InitialRate = copy.deepcopy(control)
        for Name in self.Names:
            if self.InitialRate[Name] is False:
                self.InitialRate[Name] = 0.0
            else:                
                for ele in range(len(self.InitialRate[Name])):
                    if self.par[Name][ele] is False:
                        self.InitialRate[Name][ele] = 0.0
                    else:
                        self.InitialRate[Name][ele] = Rate
        
        self.AcceptanceCounter = 0
                
    def Select(self):
        '''         
        (1) randomly select a parameter
        (2) If the parameter is controllable, return its indices.
        (3) Otherwise, select other controllable variable.        
        '''
        
        while True:
            Name = np.random.choice(self.Names)
            if self.control[Name] is False:
                continue
            else:
                ind = np.random.choice(range(len(self.control[Name])))
                if self.control[Name][ind] is False:
                    continue
                else:
                    break
                                    
        return Name, ind
        
    def Change(self,Name,ind):
        ''' Change a parameter
            (a) Direction of the change is random.        
            (b) The magnitude of the change is restricted to be (0,rate)
        '''
        par1 = copy.deepcopy(self.par)        
                
        change = np.random.normal(0.,self.InitialRate[Name][ind])                    
        
        Par = par1[Name][ind]*(1+change)
        if Par > self.lb[Name][ind] and Par < self.ub[Name][ind]:
            par1[Name][ind] = Par
        
        return par1
    
    def chi2(self,data1):
        ''' Calculate the chi2. We use minimum chi-squared method. '''
        if self.weight is None:
            chi2 = np.sum(np.absolute(self.data - data1)**2)
        else:
            chi2 = np.sum(np.absolute(self.data - data1)**2/np.absolute(self.weight)**2)
        return chi2        
    
    def Run(self):
        ''' Perform Markov chain Monte Carlo. 
        If the likelihood ratio is larger than unity, accept it! 
        If not, generate a random number. 
        If the random number is smaller than the acceptance ratio, 
        accept it, although the chi2 increases. 
        Otherwise, reject the result.'''
        
        ''' During the burnin run, parameters are not saved. '''
        data0 = Discrete(self.frequency,self.par)        
        chi2b = self.chi2(data0)             
        for run in tqdm(range(self.burnin)):
            Name, ind = self.Select()
            par1 = self.Change(Name,ind)
            data1 = Discrete(self.frequency,par1)
            chi2a = self.chi2(data1)
            if chi2b > chi2a:
                self.par = copy.deepcopy(par1)
                chi2b = chi2a                
            else:
                criterion = np.random.uniform()
                ratio = np.exp(-(chi2a+chi2b)/4)
                if criterion < ratio:
                    self.par = copy.deepcopy(par1)
                    chi2b = chi2a                    
                else:
                    continue
        
        ''' During the production run, parameters and chi2s are saved. '''
        chi2s = np.zeros((int(self.production),))
        parameters = []        
        for run in tqdm(range(self.production)):
            Name, ind = self.Select()
            par1 = self.Change(Name,ind)
            data1 = Discrete(self.frequency,par1)
            chi2a = self.chi2(data1)
            if chi2b > chi2a:                
                self.par = copy.deepcopy(par1)
                chi2s[run] = chi2a
                chi2b = chi2a
                parameters.append(self.par)
                self.AcceptanceCounter += 1
            else:
                criterion = np.random.uniform()
                ratio = np.exp(-(chi2a+chi2b)/4)
                if criterion < ratio:
                    self.par = copy.deepcopy(par1)
                    chi2s[run] = chi2a
                    chi2b = chi2a
                    parameters.append(self.par)
                    self.AcceptanceCounter += 1
                else:
                    parameters.append(self.par)
                    chi2s[run] = chi2b                    
                    continue         
        print("\nAcceptance Rate: "+str(self.AcceptanceCounter/(self.production)))
        minID = np.argmin(chi2s)
        par2 = parameters[minID]
        
        ''' Align chains according to their names. The chain can be used for error analysis.'''
        Chain = {}                
        for Name in list(self.Names):
            Chain[Name] = list(map(lambda x: x[Name], parameters))
            Chain[Name] = np.asarray(Chain[Name])
            
        return chi2s, Chain, par2
    
    def Examine(self,Chain,TruePar=None,filename=None):
        ''' Examine the chain. '''
        plt.rcParams['font.family'] = 'serif'
        ChainNum = 0
        Controllables = []
        for Name in list(Chain.keys()):
            if Chain[Name].shape[1] == 0:
                continue
            else:                                         
                if Chain[Name].shape[1] == 1 or Chain[Name][0][0] is not None:
                    Controllables.append(Name)
                    ChainNum += 1
                else:        
                    for Sub in range(Chain[Name].shape[1]):
                        if Chain[Name][0,Sub] is not None:                            
                            Controllables.append(Name)
                            ChainNum += 1        
        figure, axis = plt.subplots(ChainNum, 
                                    figsize=(8,int(2*ChainNum)),
                                    sharex=True)
        for ele in range(ChainNum):
            ax = axis[ele]
            
            name = Controllables[ele]                                    
            if Chain[name].shape[1] == 1:
                ax.plot(Chain[name],"k",alpha=0.8)
            else:
                for Sub in range(Chain[name].shape[1]):
                    ax.plot(Chain[name][:,Sub],"k",alpha=0.8)
            if TruePar is not None:
                if TruePar[name] is not None:                     
                    for Sub in range(len(TruePar[name])):
                        par2 = np.repeat(TruePar[name][Sub],len(Chain[name]))
                        ax.plot(par2,"r")
            ax.set_xlim(0, len(Chain[name]))
            ax.set_ylabel(name)
            ax.yaxis.set_label_coords(-0.1, 0.5)        
        axis[-1].set_xlabel("step number")
        figure.tight_layout()
        if filename is not None:
            figure.savefig(filename,dpi=500)
        else:
            plt.show()
        
        
        
        