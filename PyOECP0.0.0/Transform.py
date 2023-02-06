#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 29 18:13:24 2021

@author: tae-jun_yoon
"""

import numpy as np

from scipy.signal import savgol_filter
from scipy.optimize import newton
from PyOECP import References

def ListReferences():
    AvailableReferences = dir(References)
    for EachReference in AvailableReferences:
        if '_' in EachReference and '__' not in EachReference:
            print(EachReference)
            
def ReferenceDetail(ReferenceName):
    from pprint import pprint
    ''' Print out detailed information about the reference spectra. 
    The Reference_name should be "string".
    '''    
    import matplotlib.pyplot as plt
    plt.figure(figsize=(5,5),dpi=250)
    frequency = np.array([1e9])
    data = eval('References.'+ReferenceName)(frequency)
    minimum_frequency = data['minFREQ']
    maximum_frequency = data['maxFREQ']
    frequency = np.logspace(np.log10(minimum_frequency),np.log10(maximum_frequency),100)
    data = eval('References.'+ReferenceName)(frequency)
    epsilon = data['epsilon']   
    plt.semilogx(frequency,np.real(epsilon),'r')
    plt.semilogx(frequency,-np.imag(epsilon),'b')
    plt.title(ReferenceName)
    plt.xlabel('frequency [Hz]')
    plt.ylabel('Complex permittivity')
    data.pop('epsilon')
    data.pop('frequency')
    pprint(data)
    plt.show()
    return None
    
class Capacitance:
    ''' Capacitance model based on Marsland & Evans's work. It is also denoted as M&E simple in some literature. '''
    
    def __init__(self,frequency,S11m,S11r1,S11r2,S11r3,
                 m1='Short',m2=None,m3=None,temperature=25,
                 Window=None,concentrations=None):
        self.frequency = frequency
        self.S11m = S11m[:,1] + 1j*S11m[:,2]
        self.S11r1 = S11r1[:,1] + 1j*S11r1[:,2]
        self.S11r2 = S11r2[:,1] + 1j*S11r2[:,2]
        self.S11r3 = S11r3[:,1] + 1j*S11r3[:,2]
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.temperature = temperature
        self.Window = Window
        self.concentrations = concentrations                

    def Calculate(self):
        if self.concentrations is None:
            self.concentrations = -1 * np.ones((4,))
                
        func = getattr(References,self.m2)    
        eps2 = func(self.frequency,self.temperature,self.concentrations[1])['epsilon']
        func = getattr(References,self.m3)    
        eps3 = func(self.frequency,self.temperature,self.concentrations[2])['epsilon']
                
        # Calculate Gn from four references
        d13 = self.S11r1-self.S11r3
        d21 = self.S11r2-self.S11r1 
        d32 = self.S11r3-self.S11r2        
        
        # Initial guess (Capacitance model)
        dm1 = self.S11m - self.S11r1
        dm2 = self.S11m - self.S11r2
        dm3 = self.S11m - self.S11r3
        
        epsilon = -(dm2*d13)/(dm1*d32)*eps3 - (dm3*d21)/(dm1*d32)*eps2
                
        if self.Window is not None:        
            e1 = savgol_filter(np.real(epsilon),self.Window,3)
            e2 = savgol_filter(np.imag(epsilon),self.Window,3)
            epsilon = e1 + 1j*e2
            
        return epsilon
        

class Antenna:
    ''' Antenna model 
    Since this model is more robust than capacitance model and faster than Komarov model, 
    this model is recommended to be a default model.
    Citation Information
    Marsland, T. P., & Evans, S. (1987, August). Dielectric measurements with an open-ended coaxial probe. In IEE Proceedings H (Microwaves, Antennas and Propagation) (Vol. 134, No. 4, pp. 341-349). IET Digital Library.
    To use the reference solutions (either binary or pure) conveniently, "concentrations" are always used as basic inputs.
    If a function called does not require concentration (e.g., pure liquids), assign -1 for them. This is not required when all reference liquids are pure.    
    '''
    
    def __init__(self,frequency,S11m,S11r1,S11r2,S11r3,S11r4,
                 m1='Short',m2=None,m3=None,m4=None,temperature=25,
                 Window=None,concentrations=None):
        self.frequency = frequency
        self.S11m = S11m[:,1] + 1j*S11m[:,2]
        self.S11r1 = S11r1[:,1] + 1j*S11r1[:,2]
        self.S11r2 = S11r2[:,1] + 1j*S11r2[:,2]
        self.S11r3 = S11r3[:,1] + 1j*S11r3[:,2]
        self.S11r4 = S11r4[:,1] + 1j*S11r4[:,2]
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.m4 = m4
        self.temperature = temperature
        self.Window = Window
        self.concentrations = concentrations
                
    def Mother(self,x,a,b):
        return a*x**(5/2) + x + b
    
    def Son(self,x,a,b):
        return a*(5/2)*x**(3/2) + 1
    
    def Calculate(self):
        if self.concentrations is None:
            self.concentrations = -1 * np.ones((4,))
        
        ''' We don't need eps1 if we stick to Marsland-Evans model. (Short)
        func = getattr(References,self.m1)
        eps1 = func(self.temperature,self.concentrations[0])['epsilon']
        '''
        func = getattr(References,self.m2)    
        e2 = func(self.frequency,self.temperature,self.concentrations[1])['epsilon']
        func = getattr(References,self.m3)    
        e3 = func(self.frequency,self.temperature,self.concentrations[2])['epsilon']
        func = getattr(References,self.m4)    
        e4 = func(self.frequency,self.temperature,self.concentrations[3])['epsilon']
        
        # Calculate Gn from four references
        d13 = -self.S11r1+self.S11r3
        d21 = -self.S11r2+self.S11r1 
        d32 = -self.S11r3+self.S11r2
        d41 = -self.S11r4+self.S11r1
        d42 = -self.S11r4+self.S11r2
        d43 = -self.S11r4+self.S11r3
        
        
        
        Gn = -(d41*d32*e4+d42*d13*e3+d43*d21*e2)/(d41*d32*e4**(5/2)+d42*d13*e3**(5/2)+d43*d21*e2**(5/2))
        
        yy2 = e2 + Gn*e2**(5/2)
        yy3 = e3 + Gn*e3**(5/2)
        
        # Initial guess (Capacitance model)
        dm1 = self.S11m - self.S11r1
        dm2 = self.S11m - self.S11r2
        dm3 = self.S11m - self.S11r3
        em = -(dm2*d13)/(dm1*d32)*e3 - (dm3*d21)/(dm1*d32)*e2
        
        b = dm2*d13/(dm1*d32)*yy3 + dm3*d21/(dm1*d32)*yy2
        
        epsilon = np.copy(em)
        for i in range(len(em)):
            epsilon[i] = newton(self.Mother, em[i], args=(Gn[i],b[i]), maxiter=10000)
            
        if self.Window is not None:        
            e1 = savgol_filter(np.real(epsilon),self.Window,3)
            e2 = savgol_filter(np.imag(epsilon),self.Window,3)
            epsilon = e1 + 1j*e2
            
        return epsilon
