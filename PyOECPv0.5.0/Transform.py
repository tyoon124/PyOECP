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
from scipy.special import factorial
from scipy.special import factorial2
from scipy.linalg import solve
from scipy.interpolate import splprep, splev

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
    
class Stuchly:
    ''' Capacitance model based on Marsland & Evans's work. 
    It is also denoted as M&E simple in some literature. '''
    
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
        
class Marsland:
    ''' Antenna model 
    The Lumped capacitance model proposed by Marsland and Evans.
    It is recommended as one of two default models.
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
        """
        if self.Window is not None:
            tck, u = splprep([self.frequency, np.real(epsilon)])
            _, e1 = splev(u, tck)
            tck, u = splprep([self.frequency, np.imag(epsilon)])
            _, e2 = splev(u, tck)
            epsilon = e1 + 1j*e2
        """  
        return epsilon

class Komarov:
    ''' This class implements the methodology proposed 
    by Komarov.
    
    This model is recommended as one of two default models.
    
    If you use this package, please cite the following article.
    [1] Komarov, S. A., Komarov, A. S., Barber, D. G., 
        Lemes, M. J., & Rysgaard, S. (2016). 
        Open-ended coaxial probe technique for 
        dielectric spectroscopy of 
        artificially grown sea ice. 
        IEEE Transactions on Geoscience and Remote Sensing, 
        54(8), 4941-4951.
    [2] Yoon, T. J., Maerzke, K. A., Currier, R. P., 
        & Findikoglu, A. T. (2021). 
        PyOECP: A flexible open-source software library 
        for estimating and modeling the complex permittivity 
        based on the open-ended coaxial probe (OECP) technique. 
        arXiv preprint arXiv:2109.14889.
        
    -----------------------
    Parameters
    (1) frequency: numpy array (vector)
    (2) S11m: MUT reflection coefficient
    (3) S11X: reflection coefficients of reference materials
    (4) epsilonX: complex permittivity of reference materials     
    (5) Ri and Ro: inner and outer radius of the OECP. 
                   Units should be in mm.
    (6) epsilonC: Dielectric constant of the insulating material. The default value is for Teflon.
    (7) M: The number of terms. The default value is 50.
    ------------
    Constants
    (1) c: speed of light (m/s)
    
    ------------
    '''
        
    def __init__(self,frequency,S11m,
                 S11r1,S11r2,S11r3,
                 m1,m2,m3,
                 Ri, Ro, epsilonC,
                 temperature=25,Window=None,
                 concentrations=None,
                 M=25):
        self.frequency = frequency
        self.omega = 2*np.pi*self.frequency
        self.S11m = S11m[:,1] + 1j*S11m[:,2]
        self.S11r1 = S11r1[:,1] + 1j*S11r1[:,2]
        self.S11r2 = S11r2[:,1] + 1j*S11r2[:,2]
        self.S11r3 = S11r3[:,1] + 1j*S11r3[:,2]
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.epsilonC = epsilonC
        self.temperature = temperature
        self.Window = Window
        if concentrations is None:
            self.concentrations = -1 * np.ones((4,))
        else:
            self.concentrations = concentrations
            
        ''' Unit conversion of the inner and outer radii '''
        self.Ri = Ri * 1e-3
        self.Ro = Ro * 1e-3
        self.M = M
                
        ''' Constants in SI units '''
        self.epsilon0 = 8.8541878128e-12 
        self.c = 299792458
        
        self.k0 = self.omega / self.c
        self.kc = self.k0 * np.sqrt(self.epsilonC)
        
        func = getattr(References,self.m1)    
        self.epsilon1 = func(self.frequency,self.temperature,self.concentrations[0])['epsilon']
        func = getattr(References,self.m2)    
        self.epsilon2 = func(self.frequency,self.temperature,self.concentrations[1])['epsilon']
        func = getattr(References,self.m3)    
        self.epsilon3 = func(self.frequency,self.temperature,self.concentrations[2])['epsilon']
        
        self.epsilon1 *= self.epsilon0
        self.epsilon2 *= self.epsilon0
        self.epsilon3 *= self.epsilon0
        self.epsilonC = epsilonC*self.epsilon0
        
        ''' Calculate G and I coefficients '''
        self.G = self.Gs()
        self.I = self.Is()
        
        ''' Calculate the transform coefficients '''
        self.y1 = self.mother(self.k0*np.sqrt(self.epsilon1))
        self.y2 = self.mother(self.k0*np.sqrt(self.epsilon2))
        self.y3 = self.mother(self.k0*np.sqrt(self.epsilon3))
        
        self.As, self.Bs, self.Cs = self.Bilinear()
        
        ''' Calculate the measured (normalized) admittance '''
        self.y = (self.As*self.S11m - self.Bs) / (self.S11m - self.Cs)
        
        ''' Initialize the wavenumber '''
        self.ko = np.copy(self.As)
        for i in range(len(self.ko)):
            initials = np.roots([self.I[1],
                                 self.I[0],
                                 (-1j*np.pi*self.kc[i]*
                                  np.log(self.Ro/self.Ri)*
                                  self.y[i])])
            initials = np.concatenate((np.sqrt(initials),
                                       -np.sqrt(initials)))
            self.ko[i] = initials[np.logical_and(np.real(initials)>0,
                                                   np.imag(initials)<0)]
        
        ''' Now apply the Newton-Raphson method '''
        etol = 1e-10
        while True:
            self.kn = np.copy(self.ko)
            for i in range(len(self.ko)):
                self.kn = (self.ko - 
                           (self.mother(self.ko)-self.y)/
                           (self.son(self.ko)))
            criterion = np.sum(np.absolute(self.kn - self.ko))
            if criterion < etol:
                break
            else:
                self.ko = self.kn
        
        self.epsilon = (self.kn / self.k0)**2
        self.epsilon = self.epsilon/self.epsilon0
        
        if self.Window is not None:        
            e1 = savgol_filter(np.real(self.epsilon),self.Window,3)
            e2 = savgol_filter(np.imag(self.epsilon),self.Window,3)
            self.epsilon = e1 + 1j*e2
        """
        if self.Window is not None:
            tck, u = splprep([self.frequency, np.real(self.epsilon)])
            _, e1 = splev(u, tck)
            tck, u = splprep([self.frequency, np.imag(self.epsilon)])
            _, e2 = splev(u, tck)
            self.epsilon = e1 + 1j*e2
        """
        
    def Gs(self):
        ''' Calculate the coefficients G '''
        G = np.zeros((len(self.frequency),),dtype=complex)
        for m in range(1,self.M+1):
            term1 = (-1)**(m+1)/(2**(m+2))*factorial(m)/factorial2(2*m+1)
            term2 = 0
            for n in range(1,m+1):
                term21 = self.Ri**(2*n) - self.Ro**(2*n)
                term22 = self.Ri**(2*(m-n+1)) - self.Ro**(2*(m-n+1))
                term23 = (factorial(n)*factorial(m-n+1))**2
                term2 = term2 + term21*term22/term23
            G[m-1] = term1 * term2
            
        return G

    def Is(self):
        ''' Calculate the coefficients I '''
        I = np.zeros((len(self.frequency),),dtype=complex)
        theta = np.linspace(0,np.pi,10000)
        intgr = self.Ri**2+self.Ro**2-2*self.Ri*self.Ro*np.cos(theta)
        for m in range(1,self.M+1):
            term1 = (-1)**(m+1)/((2*m-1)*factorial(2*m-1))
            intgr = np.sqrt(intgr)**(2*m-1)
            term21 = 2 * np.trapz(y=intgr,x=theta)
            term22 = 2**(3*m-1)*factorial(m-1)*(self.Ri**(2*m-1)+self.Ro**(2*m-1))/factorial2(2*m-1)
            term2 = term21 - term22
            I[m-1] = term1*term2
        
        return I
    
    def mother(self, k):
        ''' noramlized input admittance in terms of wavenumber '''
        term1 = k**2/(self.kc*np.log(self.Ro/self.Ri))
        term2 = np.zeros((len(self.frequency),),
                         dtype=complex)
        for m in range(1,self.M+1):
            term2 = (term2 + 
                     self.G[m-1]*k**(2*m+1) - 
                     1j/np.pi*self.I[m-1]*k**(2*m-2))
        y = term1 * term2
        
        return y
    
    def son(self,k):
        ''' Differentiation of the yb equation '''
        term1 = 1/(self.kc*np.log(self.Ro/self.Ri))
        term2 = np.zeros((len(self.frequency),),
                         dtype=complex)
        for m in range(1,self.M+1):
            term2 = (term2 + 
                     (2*m+3)*self.G[m-1]*k**(2*m+2) - 
                     (2*m)*1j/np.pi*self.I[m-1]*k**(2*m-1))
        dydk = term1 * term2
        return dydk
    
    def Bilinear(self):
        ''' Obtain the coefficients for bilinear transform '''
        As = np.zeros((len(self.frequency),),dtype=complex)
        Bs = np.copy(As)
        Cs = np.copy(As)
            
        for i in range(len(self.frequency)):
            if np.isinf(self.y1[0]):
                Cs[i] = self.S11r1[i]
                A = np.array([[self.S11r2[i], -1],
                              [self.S11r3[i], -1]],
                             dtype = complex)
                b = np.array([self.y2[i]*self.S11r2[i]-self.y2[i]*Cs[i],
                              self.y3[i]*self.S11r3[i]-self.y3[i]*Cs[i]],
                              dtype=complex).T
                As[i],Bs[i] = solve(A, b)
            elif np.isinf(self.y2[0]):
                Cs[i] = self.S11r2[i]
                A = np.array([[self.S11r1[i], -1],
                              [self.S11r3[i], -1]],
                             dtype = complex)
                b = np.array([self.y1[i]*self.S11r1[i]-self.y1[i]*Cs[i],
                              self.y3[i]*self.S11r3[i]-self.y3[i]*Cs[i]],
                              dtype=complex).T
                As[i],Bs[i] = solve(A, b)
            elif np.isinf(self.y2[0]):
                Cs[i] = self.S11r3[i]
                A = np.array([[self.S11r1[i], -1],
                              [self.S11r2[i], -1]],
                             dtype = complex)
                b = np.array([self.y1[i]*self.S11r1[i]-self.y1[i]*Cs[i],
                              self.y2[i]*self.S11r2[i]-self.y2[i]*Cs[i]],
                              dtype=complex).T
                As[i],Bs[i] = solve(A, b)
            else:
                A = np.array([[self.S11r1[i], -1, self.y1[i]],
                              [self.S11r2[i], -1, self.y2[i]],
                              [self.S11r3[i], -1, self.y3[i]]],
                             dtype = complex)
                b = np.array([self.y1[i]*self.S11r1[i],
                              self.y2[i]*self.S11r2[i],
                              self.y3[i]*self.S11r3[i]],dtype=complex).T
                As[i],Bs[i],Cs[i] = solve(A, b)
        
        return As, Bs, Cs