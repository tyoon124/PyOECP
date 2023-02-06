#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 29 14:15:06 2021

@author: tae-jun_yoon
"""

'''
This script contains reference liquid data.
The functions are in the form of "Material_author".
Citation information is also included.
Please cite the original articles when you use the data.
Temperature is in celsius (C) and frequency is in Hz.
To print out what substances are available, use "Transform.List_References()" function.
All interpolations are performed in "linear" mode if not remarked elsewhere.
'''

import numpy as np
from PyOECP import Models
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d

def Parser(filename,header=8,footer=1):
    ''' This function loads the reflection coefficient data '''
    infile = np.genfromtxt(filename,skip_header = header,
                           skip_footer = footer,
                           delimiter=',')
    frequency = infile[:,0]
    real = infile[:,1] 
    imaginary = infile[:,2]
    
    reflection = np.column_stack((frequency,real,imaginary))
    return reflection        

def Open(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = -np.inf
    data['maxFREQ'] = np.inf
    data['Citation'] = 'None'
    data['minTemperature (C)'] = -np.inf
    data['maxTemperature (C)'] = np.inf
    data['frequency'] = frequency
    data['Temperature (C)'] = 'Measured'
    data['epsilon'] = np.ones((len(frequency)),dtype=complex)
        
    return data
    
def Short(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = -np.inf
    data['maxFREQ'] = np.inf
    data['Citation'] = 'None'
    data['minTemperature (C)'] = -np.inf
    data['maxTemperature (C)'] = np.inf
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
            
    return data

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Pure solvents
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
def Acetone_Onimisi(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = .1e9
    data['maxFREQ'] = 10e9
    data['Citation'] = 'Onimisi, M. Y., Ikyumbur, J. T., Abdu, S. G., & Hemba, E. C. (2016). Frequency and temperature effect on dielectric properties of acetone and dimethylformamide. Physical Science International Journal, 1-8.'
    data['minTemperature (C)'] = 10
    data['maxTemperature (C)'] = 50
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
        
    Temperature = np.arange(10,55,10)
    epsilons = np.array([22.25,21.13,20.20,18.83,17.63])
    epsinf = np.array([8.69,4.55,3.34,2.70,1.32])
    tau = np.array([9.22,4.05,3.12,2.07,1.43])*1e-12  # Misprint correction
    
    f1 = interp1d(Temperature,epsilons)
    f2 = interp1d(Temperature,epsinf)
    f3 = interp1d(Temperature,tau)        
    
    par = Models.Parameters()
     
    eS, ei, times = f1(temperature), f2(temperature), f3(temperature)        
    par.Set('ei',ei)    
    par.Set('times',times)
    par.Set('magnitudes',eS - ei)
    
    par = par.Parameters()        
        
    data['epsilon'] = Models.Discrete(frequency,par)
        
    return data

def NNdimethylFormamide_Onimisi(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = .1e9
    data['maxFREQ'] = 10e9
    data['Citation'] = 'Onimisi, M. Y., Ikyumbur, J. T., Abdu, S. G., & Hemba, E. C. (2016). Frequency and temperature effect on dielectric properties of acetone and dimethylformamide. Physical Science International Journal, 1-8.'
    data['minTemperature (C)'] = 10
    data['maxTemperature (C)'] = 50
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
        
    Temperature = np.arange(10,55,10)
    epsilons = np.array([40.34,38.45,36.69,35.35,33.36])
    epsinf = np.array([3.51,2.91,3.10,2.98,3.00])
    tau = np.array([1.117,1.075,0.93,0.85,0.76])*1e-11
    
    f1 = interp1d(Temperature,epsilons)
    f2 = interp1d(Temperature,epsinf)
    f3 = interp1d(Temperature,tau)        
    
    par = Models.Parameters()
     
    eS, ei, times = f1(temperature), f2(temperature), f3(temperature)        
    par.Set('ei',ei)    
    par.Set('times',times)
    par.Set('magnitudes',eS - ei)
    
    par = par.Parameters()        
        
    data['epsilon'] = Models.Discrete(frequency,par)
        
    return data

def Butanol1_Gregory(frequency,temperature = 25,c=None):
    data = {}
    data['minFREQ'] = .1e9
    data['maxFREQ'] = 4e9
    data['Citation'] = 'Gregory, A. P., & Clarke, R. N. (2012). Tables of the complex permittivity of dielectric reference liquids at frequencies up to 5 GHz.'
    data['minTemperature (C)'] = 10
    data['maxTemperature (C)'] = 50
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remarks'] = 'See Note 26 in the original report.'
    
    from scipy.interpolate import interp1d
    
    if temperature >= 45 and temperature <= 50:    
        Temperature1 = np.arange(45,55,5)
        epsilons = np.array([15.03,14.44])
        fr = np.array([0.728,0.883])*1e9
        epsH = np.array([3.418,3.416])
        Gamma = np.array([0.051,0.047])
        f1 = interp1d(Temperature1,epsilons)
        f2 = interp1d(Temperature1,fr)
        f3 = interp1d(Temperature1,epsH)
        f4 = interp1d(Temperature1,Gamma)        
        epsilons,fr,epsH, Gamma = f1(temperature), f2(temperature), f3(temperature), f4(temperature)    
        
        data['epsilon'] =  epsH + (epsilons - epsH) / (1 + 1j*frequency/fr) - 1j*frequency*Gamma/1e9
        
    elif temperature < 45 and temperature >= 10:
        Temperature2 = np.arange(10,45,5)
        epsilons = np.array([19.54, 18.86, 18.19, 17.53, 16.89, 16.26, 15.65])
        fr1 = np.array([0.167, 0.207, 0.257, 0.318, 0.393, 0.483, 0.592])*1e9
        epsilonh = np.array([3.528, 3.527, 3.505, 3.506, 3.501, 3.493, 3.484])
        fr2 = np.array([4.878, 5.259, 5.865, 5.979, 6.433, 7.234, 8.139])*1e9
        epsinf = np.array([2.914, 2.914, 2.9, 2.921, 2.923, 2.902, 2.883])
        
        f1 = interp1d(Temperature2,epsilons)
        f2 = interp1d(Temperature2,fr1)
        f3 = interp1d(Temperature2,epsilonh)
        f4 = interp1d(Temperature2,fr2)        
        f5 = interp1d(Temperature2,epsinf)
        
        epsilons,fr1,epsilonh,fr2,epsinf = f1(temperature),f2(temperature),f3(temperature),f4(temperature),f5(temperature)
        data['epsilon'] = epsinf + (epsilons-epsilonh)/(1+1j*frequency/fr1) + (epsilonh-epsinf)/(1+1j*frequency/fr2)
        
    else:
        print('out of scope.')
        return None

    return data    

def DimethylSulfoxide_Gregory(frequency,temperature = 25,c=None):
    data = {}
    data['minFREQ'] = .1e9
    data['maxFREQ'] = 5e9
    data['Citation'] = 'Gregory, A. P., & Clarke, R. N. (2012). Tables of the complex permittivity of dielectric reference liquids at frequencies up to 5 GHz.'
    data['minTemperature (C)'] = 20
    data['maxTemperature (C)'] = 50
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remarks'] = 'See Note 4 and 33 in the original report.'
        
    Temperature = np.arange(20,55,5)
    epsilons = np.array([47.13,46.49,45.86,45.19,44.53,43.86,43.19])
    epsinf = np.array([6.802,6.501,6.357,5.984,5.828,5.637,5.410])
    fr = np.array([7.555,8.323,9.077,9.924,10.733,11.588,12.477])*1e9
    
    f1 = interp1d(Temperature,epsilons)
    f2 = interp1d(Temperature,epsinf)
    f3 = interp1d(Temperature,fr)
    
    epsilons = f1(temperature)
    epsinf = f2(temperature)
    fr = f3(temperature)
    
    data['epsilon'] =  epsinf + (epsilons - epsinf) / (1 + 1j*frequency/fr)
        
    return data

def Ethanediol_Gregory(frequency,temperature = 25,c=None):
    data = {}
    data['minFREQ'] = .1e9
    data['maxFREQ'] = 5e9
    data['Citation'] = 'Gregory, A. P., & Clarke, R. N. (2012). Tables of the complex permittivity of dielectric reference liquids at frequencies up to 5 GHz.'
    data['minTemperature (C)'] = 10
    data['maxTemperature (C)'] = 50
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remarks'] = 'See Note 4, 29, and 33 and Section 10 in the original report.'
        
    Temperature = np.arange(10,55,5)
    epsilons = np.array([44.31,43.08,41.89,40.75,39.66,38.60,37.63,36.63,35.71])
    epsinf = np.array([4.543,4.568,4.745,4.700,4.712,4.739,4.840,4.782,4.856])
    fr = np.array([0.599,0.764,0.962,1.190,1.451,1.756,2.101,2.484,2.912])*1e9
    beta = np.array([0.833,0.841,0.856,0.859,0.864,0.871,0.879,0.880,0.885])
    
    f1 = interp1d(Temperature,epsilons)
    f2 = interp1d(Temperature,epsinf)
    f3 = interp1d(Temperature,fr)
    f4 = interp1d(Temperature,beta)
    
    epsilons = f1(temperature)
    epsinf = f2(temperature)
    fr = f3(temperature)
    beta = f4(temperature)
    
    data['epsilon'] =  epsinf + (epsilons - epsinf) / (1 + 1j*frequency/fr)**beta
        
    return data

def Ethanol_Gregory(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = .1e9
    data['maxFREQ'] = 4e9
    data['Citation'] = 'Gregory, A. P., & Clarke, R. N. (2012). Tables of the complex permittivity of dielectric reference liquids at frequencies up to 5 GHz.'
    data['minTemperature (C)'] = 10
    data['maxTemperature (C)'] = 50
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remarks'] = 'See Note 25 in the original report.'
    
    from scipy.interpolate import interp1d
    
    Temperature = np.arange(10,55,5)
    epsilons = np.array([26.79, 25.95, 25.16, 24.43, 23.65, 22.88, 22.16, 21.45, 20.78])
    fr = np.array([0.596, 0.7, 0.829, 0.964, 1.124, 1.303, 1.511, 1.745, 2.01])*1e9
    epsH = np.array([4.624, 4.59, 4.531, 4.505, 4.471, 4.439, 4.41, 4.394, 4.378])
    Gamma = np.array([0.075, 0.071, 0.059, 0.056, 0.054, 0.053, 0.05, 0.049, 0.044])
    f1 = interp1d(Temperature,epsilons)
    f2 = interp1d(Temperature,fr)
    f3 = interp1d(Temperature,epsH)
    f4 = interp1d(Temperature,Gamma)
    
    epsilons,fr,epsH, Gamma = f1(temperature), f2(temperature), f3(temperature), f4(temperature)    
    
    data['epsilon'] =  epsH + (epsilons - epsH) / (1 + 1j*frequency/fr) - 1j*frequency*Gamma/1e9
    
    return data

def Methanol_Gregory(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = .1e9
    data['maxFREQ'] = 5e9
    data['Citation'] = 'Gregory, A. P., & Clarke, R. N. (2012). Tables of the complex permittivity of dielectric reference liquids at frequencies up to 5 GHz.'
    data['minTemperature (C)'] = 10
    data['maxTemperature (C)'] = 50
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remarks'] = 'See Note 33 in the original report.'
        
    Temperature = np.arange(10,55,5)
    epsilons = np.array([35.74, 34.68, 33.64, 32.66, 31.69, 30.78, 29.85, 28.95, 28.19])
    epsinf = np.array([5.818, 5.698, 5.654, 5.563, 5.45, 5.388, 5.251, 5.107, 5.224])
    fr = np.array([2.262, 2.532, 2.822, 3.141, 3.49, 3.862, 4.283, 4.738, 5.175])*1e9
    
    f1 = interp1d(Temperature,epsilons)
    f2 = interp1d(Temperature,epsinf)
    f3 = interp1d(Temperature,fr)
    
    epsilons = f1(temperature)
    epsinf = f2(temperature)
    fr = f3(temperature)
    
    data['epsilon'] =  epsinf + (epsilons - epsinf) / (1 + 1j*frequency/fr)
    
    return data

def Propanol1_Gregory(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = .1e9
    data['maxFREQ'] = 4e9
    data['Citation'] = 'Gregory, A. P., & Clarke, R. N. (2012). Tables of the complex permittivity of dielectric reference liquids at frequencies up to 5 GHz.'
    data['minTemperature (C)'] = 10
    data['maxTemperature (C)'] = 65
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remarks'] = 'See Note 26 in the original report.'
    
    from scipy.interpolate import interp1d
    
    if temperature >= 35 and temperature <= 50:    
        Temperature1 = np.arange(35,55,5)
        epsilons = np.array([19.07, 18.4, 17.76, 17.11])
        fr = np.array([0.709, 0.852, 1.019, 1.218])*1e9
        epsH = np.array([3.712, 3.697, 3.692, 3.684])
        Gamma = np.array([0.055, 0.051, 0.049, 0.047])
        f1 = interp1d(Temperature1,epsilons)
        f2 = interp1d(Temperature1,fr)
        f3 = interp1d(Temperature1,epsH)
        f4 = interp1d(Temperature1,Gamma)        
        epsilons,fr,epsH, Gamma = f1(temperature), f2(temperature), f3(temperature), f4(temperature)    
        
        data['epsilon'] =  epsH + (epsilons - epsH) / (1 + 1j*frequency/fr) - 1j*frequency*Gamma/1e9
        
    elif temperature < 35 and temperature >= 10:
        Temperature2 = np.arange(10,40,5)
        epsilons = np.array([22.61, 21.88, 21.15, 20.42, 19.75, 19.07])
        fr1 = np.array([0.268, 0.328, 0.398, 0.485, 0.586, 0.707])*1e9
        epsilonh = np.array([3.862, 3.844, 3.821, 3.804, 3.791, 3.774])
        fr2 = np.array([6.871, 7.302, 7.629, 8.29, 8.86, 9.638])*1e9
        epsinf = np.array([3.126, 3.121, 3.142, 3.129, 3.105, 3.076])
        
        f1 = interp1d(Temperature2,epsilons)
        f2 = interp1d(Temperature2,fr1)
        f3 = interp1d(Temperature2,epsilonh)
        f4 = interp1d(Temperature2,fr2)        
        f5 = interp1d(Temperature2,epsinf)
        
        epsilons,fr1,epsilonh,fr2,epsinf = f1(temperature),f2(temperature),f3(temperature),f4(temperature),f5(temperature)
        data['epsilon'] = epsinf + (epsilons-epsilonh)/(1+1j*frequency/fr1) + (epsilonh-epsinf)/(1+1j*frequency/fr2)
        
    else:
        print('out of scope.')
        return None

    return data    
    
def Propanol2_Gregory(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = .1e9
    data['maxFREQ'] = 4e9
    data['Citation'] = 'Gregory, A. P., & Clarke, R. N. (2012). Tables of the complex permittivity of dielectric reference liquids at frequencies up to 5 GHz.'
    data['minTemperature (C)'] = 10
    data['maxTemperature (C)'] = 50
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remarks'] = 'See Note 26 and 30 in the original report.'
    
    from scipy.interpolate import interp1d
    
    if temperature >= 30 and temperature <= 50:    
        Temperature1 = np.arange(30,55,5)
        epsilons = np.array([18.37, 17.65, 16.93, 16.21, 15.5])
        fr = np.array([0.565, 0.702, 0.87, 1.072, 1.315])*1e9
        epsH = np.array([3.466, 3.462, 3.458, 3.454, 3.451])
        Gamma = np.array([0.052, 0.047, 0.042, 0.038, 0.035])
        f1 = interp1d(Temperature1,epsilons)
        f2 = interp1d(Temperature1,fr)
        f3 = interp1d(Temperature1,epsH)
        f4 = interp1d(Temperature1,Gamma)        
        epsilons,fr,epsH, Gamma = f1(temperature), f2(temperature), f3(temperature), f4(temperature)    
        
        data['epsilon'] =  epsH + (epsilons - epsH) / (1 + 1j*frequency/fr) - 1j*frequency*Gamma/1e9
        
    elif temperature <= 25 and temperature >= 10:
        Temperature2 = np.arange(10,30,5)
        epsilons = np.array([21.73, 20.89, 20.11, 19.3])
        fr1 = np.array([0.217, 0.277, 0.351, 0.443])*1e9
        epsilonh = np.array([3.573, 3.566, 3.557, 3.551])
        fr2 = np.array([5.037, 5.545, 5.721, 5.999])*1e9
        epsinf = np.array([3.045, 3.035, 3.057, 3.065])
        
        f1 = interp1d(Temperature2,epsilons)
        f2 = interp1d(Temperature2,fr1)
        f3 = interp1d(Temperature2,epsilonh)
        f4 = interp1d(Temperature2,fr2)        
        f5 = interp1d(Temperature2,epsinf)
        
        epsilons,fr1,epsilonh,fr2,epsinf = f1(temperature),f2(temperature),f3(temperature),f4(temperature),f5(temperature)
        data['epsilon'] = epsinf + (epsilons-epsilonh)/(1+1j*frequency/fr1) + (epsilonh-epsinf)/(1+1j*frequency/fr2)
    else:
        print('out of scope.')
        return None

    return data    

def Acetonitrile_Stoppa(frequency,temperature = 25,c=None):
    data = {}
    data['minFREQ'] = .1e9
    data['maxFREQ'] = 80e9
    data['Citation'] = 'Stoppa, A., Nazet, A., Buchner, R., Thoman, A., & Walther, M. (2015). Dielectric response and collective dynamics of acetonitrile. Journal of Molecular Liquids, 212, 963-968.'
    data['minTemperature (C)'] = -5
    data['maxTemperature (C)'] = 65
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
        
    T = temperature + 273.15    
    eSC = np.array([278029.3114,6030.657247,32.807284,-0.067743])
    eS = eSC[0]/T**2 + eSC[1]/T + eSC[2] + eSC[3]*T
    mSC = [0,9074.259906,14.619628,-0.041952]
    magnitudes = mSC[0]/T**2 + mSC[1]/T + mSC[2] + mSC[3]*T
    tSC = [0,6798.434016,-32.648416,0.044182]
    times = tSC[0]/T**2 + tSC[1]/T + tSC[2] + tSC[3]*T
    times = times*1E-12    
    ei = eS - magnitudes
    
    par = Models.Parameters()    
    par.Set('magnitudes',magnitudes)
    par.Set('times',times)
    par.Set('ei',ei)
    par = par.Parameters()
    
    data['epsilon'] = Models.Discrete(frequency, par)
    
    return data

def Methanol_Barthel(frequency,temperature=25,c=None):
    
    data = {}
    data['minFREQ'] = 1e8
    data['maxFREQ'] = 293e9
    data['Citation'] = '(1) Barthel, J., Bachhuber, K., Buchner, R., & Hetzenauer, H. (1990). Dielectric spectra of some common solvents in the microwave region. Water and lower alcohols. Chemical physics letters, 165(4), 369-373., (2) Topics in current chemistry, Vol. I1 1, ed. F.L. Boschke (Springer, Berlin, 1983) pp. 35-144.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 32.50
    tau1 = 51.5E-12
    eps2 = 5.91
    tau2 = 7.09E-12
    eps3 = 4.90
    tau3 = 1.12E-12
    epsinf = 2.79
    
    g1 = (eps1-eps2)/(eps1-epsinf)
    g2 = (eps2-eps3)/(eps1-epsinf)
    g3 = (eps3-epsinf)/(eps1-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + (eps1-epsinf)*((g1/(1+1j*omega*tau1)) + 
                                  (g2/(1+1j*omega*tau2)) + 
                                  (g3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data

def Ethanol_Barthel(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = 1e8
    data['maxFREQ'] = 89e9
    data['Citation'] = '(1) Barthel, J., Bachhuber, K., Buchner, R., & Hetzenauer, H. (1990). Dielectric spectra of some common solvents in the microwave region. Water and lower alcohols. Chemical physics letters, 165(4), 369-373., (2) Topics in current chemistry, Vol. I1 1, ed. F.L. Boschke (Springer, Berlin, 1983) pp. 35-144.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 24.32
    tau1 = 163e-12
    eps2 = 4.49
    tau2 = 8.97e-12
    eps3 = 3.82
    tau3 = 1.81e-12
    epsinf = 2.69
    
    g1 = (eps1-eps2)/(eps1-epsinf)
    g2 = (eps2-eps3)/(eps1-epsinf)
    g3 = (eps3-epsinf)/(eps1-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + (eps1-epsinf)*((g1/(1+1j*omega*tau1)) + 
                                  (g2/(1+1j*omega*tau2)) + 
                                  (g3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data

def Propanol1_Barthel(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = 1e8
    data['maxFREQ'] = 89e9
    data['Citation'] = '(1) Barthel, J., Bachhuber, K., Buchner, R., & Hetzenauer, H. (1990). Dielectric spectra of some common solvents in the microwave region. Water and lower alcohols. Chemical physics letters, 165(4), 369-373., (2) K. Riederer, Ph.D. Thesis, University of Regensburg (1981).'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 20.43
    tau1 = 329e-12
    eps2 = 3.74
    tau2 = 15.1e-12
    eps3 = 3.20
    tau3 = 2.40e-12
    epsinf = 2.44
    
    de1 = (eps1-eps2)
    de2 = (eps2-eps3)
    de3 = (eps3-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + ((de1/(1+1j*omega*tau1)) + 
                    (de2/(1+1j*omega*tau2)) + 
                    (de3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data

def Propanol2_Barthel(frequency,temperature=25,c=None):
    data = {}    
    data['minFREQ'] = 1e8
    data['maxFREQ'] = 89e9
    data['Citation'] = '(1) Barthel, J., Bachhuber, K., Buchner, R., & Hetzenauer, H. (1990). Dielectric spectra of some common solvents in the microwave region. Water and lower alcohols. Chemical physics letters, 165(4), 369-373., (2) K. Riederer, Ph.D. Thesis, University of Regensburg (1981).'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 19.40
    tau1 = 359e-12
    eps2 = 3.47
    tau2 = 14.5e-12
    eps3 = 3.04
    tau3 = 1.96e-12
    epsinf = 2.42
    
    g1 = (eps1-eps2)/(eps1-epsinf)
    g2 = (eps2-eps3)/(eps1-epsinf)
    g3 = (eps3-epsinf)/(eps1-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + (eps1-epsinf)*((g1/(1+1j*omega*tau1)) + 
                                  (g2/(1+1j*omega*tau2)) + 
                                  (g3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data
    
def Pentanol1_DAprano(frequency,temperature=25,c=None):
    data = {}    
    data['minFREQ'] = 0.01e9
    data['maxFREQ'] = 10e9
    data['Citation'] = "D'Aprano, A., Donato, I. D., D'Arrigo, G., Bertolini, D., Cassettari, M., & Salvetti, G. (1985). Molecular association and dynamics in n-pentanol and 2-methyl-2-butanol: Ultrasonic, dielectric and viscosity studies at various temperatures. Molecular Physics, 55(2), 475-488."
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 14.70
    tau1 = 754e-12
    eps2 = 3.18
    tau2 = 26e-12
    eps3 = 2.61
    tau3 = 2.4e-12
    epsinf = 2.14
    
    g1 = (eps1-eps2)/(eps1-epsinf)
    g2 = (eps2-eps3)/(eps1-epsinf)
    g3 = (eps3-epsinf)/(eps1-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + (eps1-epsinf)*((g1/(1+1j*omega*tau1)) + 
                                  (g2/(1+1j*omega*tau2)) + 
                                  (g3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data
    
def Ethanediol_Barthel(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = 0.01e9
    data['maxFREQ'] = 70e9
    data['Citation'] = 'J. Barthel, H.J. Gores, G. Schmeer, and R. Wachter, in: F.L. Boschke (ed.), Topics in Current Chemistry 111, 35-144 (1983).'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 42.80
    tau1 = 145e-12
    eps2 = 7.30
    tau2 = 10e-12
    eps3 = 3.8
    tau3 = 0.0 # None
    epsinf = 2.14
    
    g1 = (eps1-eps2)/(eps1-epsinf)
    g2 = (eps2-eps3)/(eps1-epsinf)
    g3 = (eps3-epsinf)/(eps1-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + (eps1-epsinf)*((g1/(1+1j*omega*tau1)) + 
                                  (g2/(1+1j*omega*tau2)) + 
                                  (g3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data

def Formamide_Barthel(frequency,temperature=25,c=None):
    data = {}    
    data['minFREQ'] = 1e8
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Barthel, J., Bachhuber, K., Buchner, R., Gill, J. B., & Kleebauer, M. (1990). Dielectric spectra of some common solvents in the microwave region. Dipolar aprotic solvents and amides. Chemical physics letters, 167(1-2), 62-66.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 108.80
    tau1 = 37.3e-12
    alpha1 = 0.0057
    beta1 = 1.0
    eps2 = 7.08
    tau2 = 1.16e-12
    eps3 = 4.48
    tau3 = 0.0 # None
    epsinf = 4.48
    
    g1 = (eps1-eps2)/(eps1-epsinf)
    g2 = (eps2-eps3)/(eps1-epsinf)
    g3 = (eps3-epsinf)/(eps1-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + (eps1-epsinf)*((g1/(1+(1j*omega*tau1)**(1-alpha1))**beta1) + 
                                  (g2/(1+1j*omega*tau2)) + 
                                  (g3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data

def NmethylFormamide_Barthel(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = 1e8
    data['maxFREQ'] = 293e9
    data['Citation'] = 'Barthel, J., Bachhuber, K., Buchner, R., Gill, J. B., & Kleebauer, M. (1990). Dielectric spectra of some common solvents in the microwave region. Dipolar aprotic solvents and amides. Chemical physics letters, 167(1-2), 62-66.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 183.3
    tau1 = 128e-12
    eps2 = 6.13
    tau2 = 7.93e-12
    eps3 = 4.60
    tau3 = 0.78e-12
    epsinf = 3.20
    
    g1 = (eps1-eps2)/(eps1-epsinf)
    g2 = (eps2-eps3)/(eps1-epsinf)
    g3 = (eps3-epsinf)/(eps1-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + (eps1-epsinf)*((g1/(1+(1j*omega*tau1))) + 
                                  (g2/(1+1j*omega*tau2)) + 
                                  (g3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data

def Acetone_Wei(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = 1e8
    data['maxFREQ'] = 20e9
    data['Citation'] = 'Wei, Y. Z., & Sridhar, S. (1989). Technique for measuring the frequency‐dependent complex dielectric constants of liquids up to 20 GHz. Review of scientific instruments, 60(9), 3041-3046.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    
    eps1 = 21.2
    tau1 = 3.3e-12
    einf = 1.9
    omega = 2*np.pi*frequency
    
    eps = einf + (eps1-einf)/(1+(1j*omega*tau1))
    data['epsilon'] = eps
    
    return data

def Acetone_CalderWood(frequency,temperature=20,c=None):
    data = {}
    data['minFREQ'] = 2.9e9
    data['maxFREQ'] = 24e9
    data['Citation'] = 'Calderwood, J. H., & Smyth, C. P. (1956). Microwave Absorption and Molecular Structure in Liquids. XV. The Critical Wave Lengths of Some Short-chain Aliphatic and Cyclic Ketones and of Phenyl Ether1. Journal of the American Chemical Society, 78(7), 1295-1297.'
    data['minTemperature (C)'] = 20
    data['maxTemperature (C)'] = 20
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 21
    tau1 = 3.3e-12
    eps2 = 2
    tau2 = 0 # None
    eps3 = 2
    tau3 = 0 # None
    epsinf = 2
    
    g1 = (eps1-eps2)/(eps1-epsinf)
    g2 = (eps2-eps3)/(eps1-epsinf)
    g3 = (eps3-epsinf)/(eps1-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + (eps1-epsinf)*((g1/(1+(1j*omega*tau1))) + 
                                  (g2/(1+1j*omega*tau2)) + 
                                  (g3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data

def NNdimethylFormamide_Barthel(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = 1e8
    data['maxFREQ'] = 293e9
    data['Citation'] = 'Barthel, J., Bachhuber, K., Buchner, R., Gill, J. B., & Kleebauer, M. (1990). Dielectric spectra of some common solvents in the microwave region. Dipolar aprotic solvents and amides. Chemical physics letters, 167(1-2), 62-66.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 37.24
    tau1 = 10.4e-12
    eps2 = 4.38
    tau2 = 0.76e-12
    eps3 = 2.94
    tau3 = 0.0 # None
    epsinf = 2.94
    
    g1 = (eps1-eps2)/(eps1-epsinf)
    g2 = (eps2-eps3)/(eps1-epsinf)
    g3 = (eps3-epsinf)/(eps1-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + (eps1-epsinf)*((g1/(1+(1j*omega*tau1))) + 
                                  (g2/(1+1j*omega*tau2)) + 
                                  (g3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data

def NNdimethylAcetamide_Barthel(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = 1e8
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Barthel, J., Bachhuber, K., Buchner, R., Gill, J. B., & Kleebauer, M. (1990). Dielectric spectra of some common solvents in the microwave region. Dipolar aprotic solvents and amides. Chemical physics letters, 167(1-2), 62-66.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 38.43
    tau1 = 16e-12
    eps2 = 4.10
    tau2 = 1.33e-12
    eps3 = 3.04
    tau3 = 0.0 # None
    epsinf = 3.04
    
    g1 = (eps1-eps2)/(eps1-epsinf)
    g2 = (eps2-eps3)/(eps1-epsinf)
    g3 = (eps3-epsinf)/(eps1-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + (eps1-epsinf)*((g1/(1+(1j*omega*tau1))) + 
                                  (g2/(1+1j*omega*tau2)) + 
                                  (g3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data

def Benzonitrile_Poley(frequency,temperature=20,c=None):
    data = {}
    data['minFREQ'] = 0.8e9
    data['maxFREQ'] = 38e9
    data['Citation'] = 'Poley, J. P. (1955). Microwave dispersion of some polar liquids. Applied Scientific Research, Section A, 4(1), 337-387.'
    data['minTemperature (C)'] = 20
    data['maxTemperature (C)'] = 20
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 25.6
    tau1 = 38e-12
    eps2 = 3.9
    tau2 = 0.0 # None
    eps3 = 3.9
    tau3 = 0.0 # None
    epsinf = 3.9
    
    g1 = (eps1-eps2)/(eps1-epsinf)
    g2 = (eps2-eps3)/(eps1-epsinf)
    g3 = (eps3-epsinf)/(eps1-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + (eps1-epsinf)*((g1/(1+(1j*omega*tau1))) + 
                                  (g2/(1+1j*omega*tau2)) + 
                                  (g3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data

def Nitromethane_Chandra(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = 1.9e9
    data['maxFREQ'] = 32e9
    data['Citation'] = 'Chandra, S., & Nath, D. (1969). Dielectric Relaxation in Nitroalkanes. The Journal of Chemical Physics, 51(12), 5299-5304.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 36
    tau1 = 3.9e-12
    alpha1 = 0.07
    beta1 = 1.0
    eps2 = 2.0
    tau2 = 0.0 # None
    eps3 = 2.0
    tau3 = 0.0 # None
    epsinf = 2.0
    
    g1 = (eps1-eps2)/(eps1-epsinf)
    g2 = (eps2-eps3)/(eps1-epsinf)
    g3 = (eps3-epsinf)/(eps1-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + (eps1-epsinf)*((g1/(1+(1j*omega*tau1)**(1-alpha1))**beta1) + 
                                  (g2/(1+1j*omega*tau2)) + 
                                  (g3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data

def Pyrrolidone_Dachwitz(frequency,temperature=20,c=None):
    data = {}
    data['minFREQ'] = 1e8
    data['maxFREQ'] = 36e9
    data['Citation'] = 'Dachwitz, E., & Stockhausen, M. (1985). Dielectric Relaxation of 2‐Pyrrolidinone and Some N‐Substituted, Related Compounds in the Liquid State. Berichte der Bunsengesellschaft für physikalische Chemie, 89(9), 959-961.'
    data['minTemperature (C)'] = 20
    data['maxTemperature (C)'] = 20
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 28.8
    tau1 = 253e-12
    alpha1 = 0.0
    beta1 = 1.0
    eps2 = 11.7
    tau2 = 19e-12
    eps3 = 5.8
    tau3 = 0.0 # None
    epsinf = 5.8
    
    g1 = (eps1-eps2)/(eps1-epsinf)
    g2 = (eps2-eps3)/(eps1-epsinf)
    g3 = (eps3-epsinf)/(eps1-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + (eps1-epsinf)*((g1/(1+(1j*omega*tau1)**(1-alpha1))**beta1) + 
                                  (g2/(1+1j*omega*tau2)) + 
                                  (g3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data
    
def NmethylPyrrolidone_Dachwitz(frequency,temperature=20,c=None):
    data = {}
    data['minFREQ'] = 1e8
    data['maxFREQ'] = 36e9
    data['Citation'] = 'Dachwitz, E., & Stockhausen, M. (1985). Dielectric Relaxation of 2‐Pyrrolidinone and Some N‐Substituted, Related Compounds in the Liquid State. Berichte der Bunsengesellschaft für physikalische Chemie, 89(9), 959-961.'
    data['minTemperature (C)'] = 20
    data['maxTemperature (C)'] = 20
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 34.1
    tau1 = 21e-12
    alpha1 = 0.0
    beta1 = 1.0
    eps2 = 4.1
    tau2 = 0.0 # None
    eps3 = 4.1
    tau3 = 0.0 # None
    epsinf = 4.1
    
    g1 = (eps1-eps2)/(eps1-epsinf)
    g2 = (eps2-eps3)/(eps1-epsinf)
    g3 = (eps3-epsinf)/(eps1-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + (eps1-epsinf)*((g1/(1+(1j*omega*tau1)**(1-alpha1))**beta1) + 
                                  (g2/(1+1j*omega*tau2)) + 
                                  (g3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data

def PropyleneCarbonate_Barthel(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = 1e8
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Barthel, J., Bachhuber, K., Buchner, R., Gill, J. B., & Kleebauer, M. (1990). Dielectric spectra of some common solvents in the microwave region. Dipolar aprotic solvents and amides. Chemical physics letters, 167(1-2), 62-66.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 64.88
    tau1 = 43.1e-12
    alpha1 = 0.00
    beta1 = 0.904
    eps2 = 4.14
    tau2 = 0.0 # None
    eps3 = 4.14
    tau3 = 0.0 # None
    epsinf = 4.14
    
    g1 = (eps1-eps2)/(eps1-epsinf)
    g2 = (eps2-eps3)/(eps1-epsinf)
    g3 = (eps3-epsinf)/(eps1-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + (eps1-epsinf)*((g1/(1+(1j*omega*tau1)**(1-alpha1))**beta1) + 
                                  (g2/(1+1j*omega*tau2)) + 
                                  (g3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data

def DimethylSulfoxide_Barthel(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = 1e8
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Barthel, J., Bachhuber, K., Buchner, R., Gill, J. B., & Kleebauer, M. (1990). Dielectric spectra of some common solvents in the microwave region. Dipolar aprotic solvents and amides. Chemical physics letters, 167(1-2), 62-66.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remark'] = 'Data collected in Barthel, J. M. G., & Buchner, R. (1991). High frequency permittivity and its use in the investigation of solution properties. Pure and Applied Chemistry, 63(10), 1473-1482.'
    
    eps1 = 46.40
    tau1 = 20.5e-12
    alpha1 = 0.00
    beta1 = 0.888
    eps2 = 4.16
    tau2 = 0.0 # None
    eps3 = 4.16
    tau3 = 0.0 # None
    epsinf = 4.16
    
    g1 = (eps1-eps2)/(eps1-epsinf)
    g2 = (eps2-eps3)/(eps1-epsinf)
    g3 = (eps3-epsinf)/(eps1-epsinf)
    
    omega = 2*np.pi*frequency
    
    eps = epsinf + (eps1-epsinf)*((g1/(1+(1j*omega*tau1)**(1-alpha1))**beta1) + 
                                  (g2/(1+1j*omega*tau2)) + 
                                  (g3/(1+1j*omega*tau3)))
    data['epsilon'] = eps
    
    return data

def EthylAcetate_Birajdar(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = 1e7
    data['maxFREQ'] = 30e9
    data['Citation'] = 'Birajdar, S. S., Suryawanshi, D. B., Deshmukh, A. R., Shinde, R. V., Ingole, S. A., & Kumbharkhane, A. C. (2021). Dielectric relaxation behaviour of ethyl acetate-xylene mixtures using time domain reflectometry. Physics and Chemistry of Liquids, 59(4), 503-511.'
    data['minTemperature (C)'] = 15
    data['maxTemperature (C)'] = 25
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    
    Temperature = np.array([15,20,25])
    epsilons = np.array([7.11,6.49,6.04])
    epsinf = np.array([2.89,2.20,2.20])
    tau = np.array([4.25,4.22,4.19])*1e-12
    
    f1 = interp1d(Temperature,epsilons)
    f2 = interp1d(Temperature,epsinf)
    f3 = interp1d(Temperature,tau)        
    
    par = Models.Parameters()
     
    eS, ei, times = f1(temperature), f2(temperature), f3(temperature)        
    par.Set('ei',ei)    
    par.Set('times',times)
    par.Set('magnitudes',eS - ei)
    
    par = par.Parameters()        
        
    data['epsilon'] = Models.Discrete(frequency,par)
        
    return data

def Xylene_Birajdar(frequency,temperature=25,c=None):
    data = {}
    data['minFREQ'] = 1e7
    data['maxFREQ'] = 30e9
    data['Citation'] = 'Birajdar, S. S., Suryawanshi, D. B., Deshmukh, A. R., Shinde, R. V., Ingole, S. A., & Kumbharkhane, A. C. (2021). Dielectric relaxation behaviour of ethyl acetate-xylene mixtures using time domain reflectometry. Physics and Chemistry of Liquids, 59(4), 503-511.'
    data['minTemperature (C)'] = 15
    data['maxTemperature (C)'] = 25
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    
    Temperature = np.array([15,20,25])
    epsilons = np.array([2.81,2.40,2.23])
    epsinf = np.array([2.57,2.02,1.65])
    tau = np.array([4.20,3.58,2.75])*1e-12
    
    f1 = interp1d(Temperature,epsilons)
    f2 = interp1d(Temperature,epsinf)
    f3 = interp1d(Temperature,tau)        
    
    par = Models.Parameters()
     
    eS, ei, times = f1(temperature), f2(temperature), f3(temperature)        
    par.Set('ei',ei)    
    par.Set('times',times)
    par.Set('magnitudes',eS - ei)
    
    par = par.Parameters()        
        
    data['epsilon'] = Models.Discrete(frequency,par)
        
    return data

def Water_Kaatze(frequency,temperature=25,c=None):    
    data = {}
    data['minFREQ'] = 1.1
    data['maxFREQ'] = 57e9
    data['Citation'] = 'Kaatze, U. (1989). Complex permittivity of water as a function of frequency and temperature. Journal of Chemical and Engineering Data, 34(4), 371-374.'
    data['minTemperature (C)'] = -4.1
    data['maxTemperature (C)'] = 60    
    ei = 5.77 - 2.74E-2*temperature    
    eS = 10**(1.94404-1.991E-3*temperature)
    magnitudes = eS - ei
    times = (3.745E-15)*(1+(7E-5)*(temperature+273.15-300.65)**2)*np.exp(2.2957E3/(temperature+273.15))         
    
    par = Models.Parameters()    
    par.Set('magnitudes',magnitudes)
    par.Set('times',times)
    par.Set('ei',ei)    
    par = par.Parameters()
    
    epsilon = Models.Discrete(frequency, par)    
    data['epsilon'] = epsilon
    data['frequency'] = frequency    
    data['Temperature (C)'] = temperature
            
    return data

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Aqueous Electrolytes
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
def FormicAcidAqueous_Kaatze(frequency,temperature=25,c=1.0):
    data = {}
    data['minFREQ'] = 3e6
    data['maxFREQ'] = 70e9
    data['Citation'] = 'Kaatze, U., Menzel, K., & Pottel, R. (1991). Broad-band dielectric spectroscopy on carboxylic acid/water mixtures. Dependence upon composition. The Journal of Physical Chemistry, 95(1), 324-331.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.095
    data['maxConcentration (x)'] = 1.0
    data['Remarks'] = 'c in mole fraction of the substance in water. For conversion to mol/L, see the original article.'

    Concentration = np.array([0.095, 0.249, 0.432, 0.666, 0.837, 1])
    einf = np.array([5.1, 4.6, 5.2, 5.2, 4.7, 4.4])
    de1 = np.array([34.6, 12.7, 7.3, 0, 0, 0])
    tau1 = np.array([5.9, 4.9, 4.2, 0, 0, 0])*1e-12
    de2 = np.array([39.8, 63.4, 65.2, 66, 60, 52.6])
    tau2 = np.array([20.7, 31.1, 37.8, 54.4, 66.1, 76.7])*1e-12
    beta2 = np.array([0, 0.17, 0.09, 0.18, 0.16, 0.14])
    
    f1 = interp1d(Concentration,einf)
    f2 = interp1d(Concentration,de1)
    f3 = interp1d(Concentration,tau1)        
    f4 = interp1d(Concentration,de2)        
    f5 = interp1d(Concentration,tau2)        
    f6 = interp1d(Concentration,beta2)
        
    einf,de1,tau1,de2,tau2,beta2 = f1(c), f2(c), f3(c), f4(c), f5(c), f6(c)
    
    omega = 2*np.pi*frequency
    data['epsilon'] = einf + de1/(1+1j*omega*tau1) + de2/(1+1j*omega*tau2)**(1-beta2)
    
    return data

def AceticAcidAqueous_Kaatze(frequency,temperature=25,c=1.0):
    data = {}
    data['minFREQ'] = 3e6
    data['maxFREQ'] = 70e9
    data['Citation'] = 'Kaatze, U., Menzel, K., & Pottel, R. (1991). Broad-band dielectric spectroscopy on carboxylic acid/water mixtures. Dependence upon composition. The Journal of Physical Chemistry, 95(1), 324-331.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.02
    data['maxConcentration (x)'] = 1.0
    data['Remarks'] = 'c in mole fraction of the substance in water. For conversion to mol/L, see the original article.'

    Concentration = np.array([0.02, 0.045, 0.082, 0.136, 0.178, 0.275, 0.45, 0.613, 0.794, 0.899, 1])
    einf = np.array([5.1, 4.9, 4.7, 4.4, 4.1, 3.8, 3.3, 2.9, 3, 2.8, 2.6])
    de1 = np.array([62.9, 51.8, 38.8, 27.5, 18.6, 10.3, 2.3, 0, 0, 0, 0])
    tau1 = np.array([8.5, 9, 9.3, 9.9, 10.8, 10.7, 9.1, 0, 0, 0, 0])*1e-12
    de2 = np.array([8, 16.5, 25.8, 31.8, 35.9, 36.6, 31.5, 23.6, 13.6, 8.4, 4])
    tau2 = np.array([36.2, 39.8, 45.6, 58.1, 76.8, 93.9, 132.7, 179.9, 203.5, 199.7, 177.4])*1e-12
    beta2 = np.array([0.32, 0.33, 0.33, 0.34, 0.38, 0.38, 0.4, 0.45, 0.45, 0.48, 0.54])
    
    f1 = interp1d(Concentration,einf)
    f2 = interp1d(Concentration,de1)
    f3 = interp1d(Concentration,tau1)        
    f4 = interp1d(Concentration,de2)        
    f5 = interp1d(Concentration,tau2)        
    f6 = interp1d(Concentration,beta2)
        
    einf,de1,tau1,de2,tau2,beta2 = f1(c), f2(c), f3(c), f4(c), f5(c), f6(c)
    
    omega = 2*np.pi*frequency
    data['epsilon'] = einf + de1/(1+1j*omega*tau1) + de2/(1+1j*omega*tau2)**(1-beta2)
    
    return data

def PropionicAcidAqueous_Kaatze(frequency,temperature=25,c=0.899):
    data = {}
    data['minFREQ'] = 3e6
    data['maxFREQ'] = 70e9
    data['Citation'] = 'Kaatze, U., Menzel, K., & Pottel, R. (1991). Broad-band dielectric spectroscopy on carboxylic acid/water mixtures. Dependence upon composition. The Journal of Physical Chemistry, 95(1), 324-331.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.028
    data['maxConcentration (x)'] = 0.899
    data['Remarks'] = 'c in mole fraction of the substance in water. For conversion to mol/L, see the original article.'

    Concentration = np.array([0.028, 0.055, 0.108, 0.177, 0.263, 0.428, 0.599, 0.745, 0.899])
    einf = np.array([4.9, 4.6, 4.1, 3.6, 3.2, 2.8, 2.3, 2.4, 2.1])
    de1 = np.array([57.6, 44.7, 27.5, 16.2, 7.4, 2.1, 0, 0, 0])
    tau1 = np.array([9.7, 10.5, 11.1, 12.3, 12.3, 13, 0, 0, 0])*1e-12
    de2 = np.array([10.7, 18.4, 26.2, 27.5, 26.1, 20.1, 13.1, 7, 2.9])
    tau2 = np.array([74.6, 71, 81.5, 113.5, 136.2, 187, 229, 230, 204])*1e-12
    beta2 = np.array([0.56, 0.47, 0.43, 0.47, 0.47, 0.51, 0.58, 0.6, 0.72])
    
    f1 = interp1d(Concentration,einf)
    f2 = interp1d(Concentration,de1)
    f3 = interp1d(Concentration,tau1)        
    f4 = interp1d(Concentration,de2)        
    f5 = interp1d(Concentration,tau2)        
    f6 = interp1d(Concentration,beta2)
        
    einf,de1,tau1,de2,tau2,beta2 = f1(c), f2(c), f3(c), f4(c), f5(c), f6(c)
    
    omega = 2*np.pi*frequency
    data['epsilon'] = einf + de1/(1+1j*omega*tau1) + de2/(1+1j*omega*tau2)**(1-beta2)
    
    return data

def ButyricAcidAqueous_Kaatze(frequency,temperature=25,c=0.566):
    data = {}
    data['minFREQ'] = 3e6
    data['maxFREQ'] = 70e9
    data['Citation'] = 'Kaatze, U., Menzel, K., & Pottel, R. (1991). Broad-band dielectric spectroscopy on carboxylic acid/water mixtures. Dependence upon composition. The Journal of Physical Chemistry, 95(1), 324-331.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.088
    data['maxConcentration (x)'] = 0.566
    data['Remarks'] = 'c in mole fraction of the substance in water. For conversion to mol/L, see the original article.'

    Concentration = np.array([0.088, 0.211, 0.342, 0.436, 0.566])
    einf = np.array([4, 3, 2.5, 2.2, 1.9])
    de1 = np.array([31.4, 10.5, 3.4, 1.4, 0])
    tau1 = np.array([11.1, 12.2, 14.1, 16.1, 0])*1e-12
    de2 = np.array([18.5, 20.4, 16.7, 13.7, 10])
    tau2 = np.array([140, 187, 246, 258, 312])*1e-12
    beta2 = np.array([0.58, 0.57, 0.41, 0.38, 0.35])
    
    f1 = interp1d(Concentration,einf)
    f2 = interp1d(Concentration,de1)
    f3 = interp1d(Concentration,tau1)        
    f4 = interp1d(Concentration,de2)        
    f5 = interp1d(Concentration,tau2)        
    f6 = interp1d(Concentration,beta2)
        
    einf,de1,tau1,de2,tau2,beta2 = f1(c), f2(c), f3(c), f4(c), f5(c), f6(c)
    
    omega = 2*np.pi*frequency
    data['epsilon'] = einf + de1/(1+1j*omega*tau1) + de2/(1+1j*omega*tau2)**(1-beta2)
    
    return data

def NaClAqueous_Peyman(frequency,temperature=25,c=0):
    data = {}
    data['minFREQ'] = 0.13e9
    data['maxFREQ'] = 20e9
    data['Citation'] = 'Peyman, A., Gabriel, C., & Grant, E. H. (2007). Complex permittivity of sodium chloride solutions at microwave frequencies. Bioelectromagnetics: Journal of the Bioelectromagnetics Society, The Society for Physical Regulation in Biology and Medicine, The European Bioelectromagnetics Association, 28(4), 264-274.'
    data['minTemperature (C)'] = 5
    data['maxTemperature (C)'] = 35
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['minX [wt%]'] = 0
    data['maxX [wt%]'] = 5
    data['Remark'] = 'Concentration is in the molarity of NaCl [mol/L].'
        
    ei = 5.77 - 0.0274*temperature    
    if c < 1:        
        epssw = 10**(1.94404 - 1.991E-3*temperature)
        
        eS = epssw*(1.0 - 3.742E-4*temperature*c + 0.034*c**2 - 0.178*c + 1.515E-4*temperature - 4.929E-6*temperature**2)
        
        tauw = 3.745E-15*(1+7E-5*(temperature-27.5)**2)*np.exp(2.2957E3/(temperature+273.15))        
        times = tauw*(1.012 - 5.282E-3*temperature*c + 0.032*c**2 - 0.01*c - 1.724E-3*temperature + 3.766E-5*temperature**2)
        
        conductance = 0.174*temperature*c - 1.582*c**2 + 5.923*c        
        As = -6.348E-4*temperature*c - 5.1E-2*c**2 + 9E-2*c         
    else:
        eS = 84.328 + 0.117*temperature*c + 0.77*c**2 - 13.257*c - 0.207*temperature - 4.859E-3*temperature**2
        times = (17.76 + 0.022*temperature*c + 0.09*c**2 - 1.222*c - 0.525*temperature + 5.361E-3*temperature**2)*1E-12
        conductance = 0.061*temperature*c - 0.667*c**2 + 6.485*c - 0.439 - 1.2E-2*temperature + 3.374E-3*temperature**2
        As = 0.011 + 4.326E-4*temperature*c + 4.431E-3*c**2 + 4.754E-3*c + 1.82E-3*temperature - 6.154E-5*temperature**2
    
    magnitudes = eS - ei
    
    par = Models.Parameters()    
    par.Set('magnitudes',magnitudes)
    par.Set('times',times)
    par.Set('ei',ei)
    par.Set('As',As)
    par.Set('conductance',conductance)
    par = par.Parameters()
    
    data['epsilon'] = Models.Discrete(frequency, par)
    
    return data

def KClAqueous_Chen(frequency,temperature=25,c=1.423):
    data = {}
    data['minFREQ'] = 0.2e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Chen, T., Hefter, G., & Buchner, R. (2003). Dielectric spectroscopy of aqueous solutions of KCl and CsCl. The Journal of Physical Chemistry A, 107(20), 4025-4031.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.0519
    data['maxConcentration (x)'] = 1.423
    data['Remarks'] = 'c in mol/L of the substance in water. No interpolation is allowed at this stage (see the original article).'
    
    if c >=0.0519 and c<=1.423:
        Concentration = np.array([0.0519, 0.103, 0.1504, 0.2497, 0.3153, 0.4408, 0.6008, 0.7931, 1.104, 1.423])
        einf = np.array([6.35, 5.8, 5.21, 5.04, 4.85, 4.47, 4.31, 4.79, 5.08, 5.66])
        es = np.array([77.73, 77.41, 76.94, 76.15, 75.63, 74.52, 72.97, 71.23, 68.2, 65.93])
        tau1 = np.array([8.24, 8.17, 8.13, 7.99, 7.92, 7.8, 7.71, 7.57, 7.46, 7.41])*1e-12
        a1 = np.array([0.002, 0.011, 0.013, 0.025, 0.028, 0.034, 0.034, 0.036, 0.036, 0.034])
        ke = np.array([0.678, 1.3, 1.85, 2.99, 3.7, 5.06, 6.78, 8.79, 11.8, 15])
    else:
        print('out of range. See the remarks.')
        return None
    
    f1 = interp1d(Concentration,einf)
    f2 = interp1d(Concentration,es)
    f3 = interp1d(Concentration,tau1)        
    f4 = interp1d(Concentration,a1)
    f5 = interp1d(Concentration,ke)        
        
    einf,es,tau1,a1,ke = f1(c), f2(c), f3(c), f4(c), f5(c)
    
    omega = 2*np.pi*frequency
    data['epsilon'] = einf + (es-einf)/(1+(1j*omega*tau1)**(1-a1))
    
    return data

def CsClAqueous_Chen(frequency,temperature=25,c=1.423):
    data = {}
    data['minFREQ'] = 0.2e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Chen, T., Hefter, G., & Buchner, R. (2003). Dielectric spectroscopy of aqueous solutions of KCl and CsCl. The Journal of Physical Chemistry A, 107(20), 4025-4031.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.0504
    data['maxConcentration (x)'] = 1.760
    data['Remarks'] = 'c in mol/L of the substance in water.'
    
    if c >=0.0519 and c<=1.423:
        Concentration = np.array([0.0504, 0.0994, 0.15, 0.2499, 0.4016, 0.6009, 0.7964, 1.106, 1.506, 1.76])
        einf = np.array([6.12, 5.88, 5.74, 4.73, 5.33, 5.34, 4.32, 5.65, 6.33, 5.77])
        es = np.array([78.05, 77.67, 77.3, 76.54, 75.48, 73.87, 72.33, 70.23, 67.2, 65.49])
        tau1 = np.array([8.27, 8.22, 8.17, 7.99, 7.97, 7.85, 7.7, 7.68, 7.67, 7.66])*1e-12
        a1 = np.array([0.008, 0.017, 0.017, 0.025, 0.03, 0.032, 0.03, 0.038, 0.028, 0.009])
        ke = np.array([0.673, 1.29, 1.88, 3.03, 4.71, 6.86, 8.9, 11.9, 16.1, 18.6])
    else:
        print('out of range. See the remarks.')
        return None
    
    f1 = interp1d(Concentration,einf)
    f2 = interp1d(Concentration,es)
    f3 = interp1d(Concentration,tau1)        
    f4 = interp1d(Concentration,a1)
    f5 = interp1d(Concentration,ke)        
        
    einf,es,tau1,a1,ke = f1(c), f2(c), f3(c), f4(c), f5(c)
    
    omega = 2*np.pi*frequency
    data['epsilon'] = einf + (es-einf)/(1+(1j*omega*tau1)**(1-a1))
    
    return data

def MgSO4Aqueous_Buchner(frequency,temperature=25,c=0.01691):
    data = {}
    data['minFREQ'] = 0.2e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Buchner, R., Chen, T., & Hefter, G. (2004). Complexity in “simple” electrolyte solutions: ion pairing in MgSO4 (aq). The Journal of Physical Chemistry B, 108(7), 2365-2375.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.01691
    data['maxConcentration (x)'] = 2.236
    data['Remarks'] = 'c in mol/L of the substance in water. No interpolation is allowed, since the functional form is very complex. (See the original article)'
    
    Concentration = np.array([0.01691, 0.0295, 0.05007, 0.07296, 0.07744, 0.0998, 0.1498, 0.2001, 0.2706, 0.3634, 0.5002, 0.6472, 0.8206, 1.12, 1.351, 1.599, 1.897, 2.236])
    if c not in Concentration:
        print('No interpolation is allowed in this substance at this stage. Available concentrations are:'+str(Concentration)+'. See the cited article.')
        return None
    else:
        einf = np.array([4.6, 5.03, 6.36, 5.43, 5.44, 6.34, 6.71, 5.1, 6.75, 5.08, 6.37, 4.44, 4.71, 4.04, 5.77, 6.43, 5.78, 6.96])
        e1 = np.array([79.64, 81.24, 83.09, 84.69, 84.77, 85.92, 87.27, 88.31, 88.36, 87.05, 86.08, 83.77, 80.91, 76.75, 76.09, 73.8, 72.05, 70.92])
        tau1 = np.array([257, 276, 307, 315, 324, 412, 408, 500, 421, 400, 295, 300, 286, 288, 311, 426, 459, 499])*1e-12
        e2 = np.array([77.86, 77.59, 78.03, 78.55, 78.62, 81.19, 82.45, 83.74, 83.27, 83.01, 81.9, 80.34, 78.52, 73.18, 69.25, 67.44, 62.43, 57.18])
        tau2 = np.array([0, 0, 130, 129, 113, 156, 151, 147, 135, 126, 131, 117, 119, 115, 117, 124, 132, 122])*1e-12
        e3 = np.array([77.86, 77.59, 76.82, 76.3, 75.82, 75.58, 74.63, 73.61, 72.42, 70.7, 69.65, 67, 64.9, 60.47, 59.34, 56.32, 53.15, 48.87])
        tau3 = np.array([0, 22.3, 22, 23.4, 22, 16.2, 18.4, 19.6, 25.8, 20.8, 29.1, 21.7, 22.4, 21.9, 26.6, 24.6, 26.3, 23.2])*1e-12
        e4 = np.array([77.86, 76.9, 76.79, 74.04, 75.37, 74.8, 72.62, 71.35, 69.93, 66.93, 64.42, 59.22, 55.1, 48.32, 47.25, 42.67, 38.74, 32.37])
        tau4 = np.array([8.28, 8.29, 8.29, 8.03, 8.21, 8.25, 8.21, 8.13, 8.14, 8.16, 8.01, 7.94, 7.89, 7.75, 7.74, 7.74, 7.68, 7.66])*1e-12
        e5 = np.array([5.74, 6.52, 6.36, 6.12, 6.16, 6.34, 6.71, 6.12, 6.75, 6.71, 6.37, 6.81, 7.13, 7.03, 5.77, 6.43, 5.78, 6.96])
        tau5 = np.array([1, 1, 0, 1, 1, 0, 0, 1, 0, 1.5, 0, 0.7, 1, 1, 0, 0, 0, 0])*1e-12
                
        
        einf = einf[Concentration==c]
        e1 = e1[Concentration==c]
        tau1 = tau1[Concentration==c]
        e2 = e2[Concentration==c]
        tau2 = tau2[Concentration==c]
        e3 = e3[Concentration==c]
        tau3 = tau3[Concentration==c]
        e4 = e4[Concentration==c]
        tau4 = tau4[Concentration==c]
        e5 = e5[Concentration==c]
        tau5 = tau5[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + (e1-e2)/(1+1j*omega*tau1)
                           + (e2-e3)/(1+1j*omega*tau2)
                           + (e3-e4)/(1+1j*omega*tau3)
                           + (e4-e5)/(1+1j*omega*tau4)
                           + (e5-einf)/(1+1j*omega*tau5))
        return data
    
def LiClAqueous_Wachter(frequency,temperature=25,c=0.0498):
    data = {}
    data['minFREQ'] = 0.2e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Wachter, W., Fernandez, Š., Buchner, R., & Hefter, G. (2007). Ion association and hydration in aqueous solutions of LiCl and Li2SO4 by dielectric spectroscopy. The Journal of Physical Chemistry B, 111(30), 9010-9017.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.0498
    data['maxConcentration (x)'] = 0.985
    data['Remarks'] = 'c in mol/L of the substance in water. No interpolation is allowed (See the original article).'
    
    Concentration = np.array([0.0498, 0.0879, 0.0995, 0.149, 0.217, 0.298, 0.396, 0.487, 0.669, 0.786, 0.985])
    if c not in Concentration:
        print('No interpolation is allowed in this substance at this stage. Available concentrations are:'+str(Concentration)+'. See the cited article.')
        return None
    else:
        einf = np.array([6.83, 4.35, 6.59, 6.17, 4.61, 5.42, 5.3, 3.78, 2.88, 5.45, 4.81])
        e1 = np.array([77.86, 77.3, 77.26, 76.64, 75.61, 74.91, 73.25, 72, 69.86, 68.1, 65.59])
        tau1 = np.array([327, 220, 252, 206, 156, 238, 158, 140, 240, 276, 17.5])*1e-12
        e2 = np.array([77.54, 76.79, 76.66, 75.81, 74.6, 73.79, 72.55, 71.22, 69.16, 67.9, 58.03])
        tau2 = np.array([8.52, 8.41, 8.44, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 7.32])*1e-12
        e3 = np.array([6.83, 6.04, 6.59, 76.14, 74.36, 71.18, 68.81, 68.37, 64.05, 60.48, 5.58])
        tau3 = np.array([0, 0.5, 0, 8.37, 8.33, 7.97, 7.83, 7.95, 7.75, 7.37, 0.5])*1e-12
        e4 = np.array([6.83, 4.35, 6.59, 6.17, 6.02, 5.42, 5.3, 5.55, 6.3, 5.45, 4.81])
        tau4 = np.array([0, 0, 0, 0, 0.5, 0, 0, 0.5, 0.5, 0, 0])*1e-12
        Ke = np.array([0.515, 0.871, 0.977, 1.42, 1.99, 2.65, 3.41, 4.09, 5.36, 6.13, 7.36])
        
        einf = einf[Concentration==c]
        e1 = e1[Concentration==c]
        tau1 = tau1[Concentration==c]
        e2 = e2[Concentration==c]
        tau2 = tau2[Concentration==c]
        e3 = e3[Concentration==c]
        tau3 = tau3[Concentration==c]
        e4 = e4[Concentration==c]
        tau4 = tau4[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + (e1-e2)/(1+1j*omega*tau1)
                           + (e2-e3)/(1+1j*omega*tau2)
                           + (e3-e4)/(1+1j*omega*tau3)
                           + (e4-einf)/(1+1j*omega*tau4))
        return data

def Li2SO4Aqueous_Wachter(frequency,temperature=25,c=0.0747):
    data = {}
    data['minFREQ'] = 0.2e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Wachter, W., Fernandez, Š., Buchner, R., & Hefter, G. (2007). Ion association and hydration in aqueous solutions of LiCl and Li2SO4 by dielectric spectroscopy. The Journal of Physical Chemistry B, 111(30), 9010-9017.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.0747
    data['maxConcentration (x)'] = 1.9
    data['Remarks'] = 'c in mol/L of the substance in water. No interpolation is allowed (See the original article).'
    
    Concentration = np.array([0.0747, 0.149, 0.298, 0.494, 0.785, 1.17, 1.9])
    if c not in Concentration:
        print('No interpolation is allowed in this substance at this stage. Available concentrations are:'+str(Concentration)+'. See the cited article.')
        return None
    else:
        einf = np.array([3.4, 1.6, 5.83, 4.84, 3.42, 4.42, 3.93])
        e1 = np.array([78.95, 78.26, 76.52, 73.87, 69.86, 65.34, 57.96])
        tau1 = np.array([222, 163, 140, 130, 125, 130, 142])*1e-12
        e2 = np.array([76.01, 74.24, 71.7, 69.15, 65.8, 61.3, 53.43])
        tau2 = np.array([8.35, 45.2, 23.7, 28.8, 30, 27.4, 26])*1e-12
        e3 = np.array([6.05, 74.19, 67.34, 63.51, 57.37, 49.59, 37.75])
        tau3 = np.array([0.32, 8.4, 7.94, 7.98, 7.93, 7.71, 7.36])*1e-12
        e4 = np.array([3.4, 6.72, 5.64, 6.34, 6.97, 7.17, 7.66])
        tau4 = np.array([0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])*1e-12
        Ke = np.array([1.12, 1.95, 3.26, 4.71, 6.22, 7.49, 8.39]) # S/m
        
        einf = einf[Concentration==c]
        e1 = e1[Concentration==c]
        tau1 = tau1[Concentration==c]
        e2 = e2[Concentration==c]
        tau2 = tau2[Concentration==c]
        e3 = e3[Concentration==c]
        tau3 = tau3[Concentration==c]
        e4 = e4[Concentration==c]
        tau4 = tau4[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + (e1-e2)/(1+1j*omega*tau1)
                           + (e2-e3)/(1+1j*omega*tau2)
                           + (e3-e4)/(1+1j*omega*tau3)
                           + (e4-einf)/(1+1j*omega*tau4))
        return data
    
def Al2SO43Aqueous_Schrodle(frequency,temperature=25,c=0.0115):
    data = {}
    data['minFREQ'] = 0.2e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Schr\"odle, S., Rudolph, W. W., Hefter, G., & Buchner, R. (2007). Ion association and hydration in 3: 2 electrolyte solutions by dielectric spectroscopy: Aluminum sulfate. Geochimica et Cosmochimica Acta, 71(22), 5287-5300.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.0115
    data['maxConcentration (x)'] = 0.6506
    data['Remarks'] = 'c in mol/L of the substance in water. No interpolation is allowed (See the original article).'
    
    Concentration = np.array([0.0115, 0.0425, 0.0855, 0.1054, 0.1601, 0.2343, 0.4302, 0.6506])
    if c not in Concentration:
        print('No interpolation is allowed in this substance at this stage. Available concentrations are:'+str(Concentration)+'. See the cited article.')
        return None
    else:
        einf = np.array([3.53, 4.1, 3.01, 2.61, 2.64, 3.33, 3.48, 6.44])
        e1 = np.array([82.73, 89.48, 92.39, 92.79, 92.89, 91.7, 83.84, 73.91])
        tau1 = np.array([352, 357, 326, 288, 291, 270, 298, 346])*1e-12
        e2 = np.array([78.34, 81.05, 83.02, 82.21, 83.09, 81.32, 74.35, 67.64])
        tau2 = np.array([179, 174, 170, 163, 161, 157, 154, 172])*1e-12
        e3 = np.array([77.32, 75.4, 73.12, 72.32, 69.76, 67.08, 59.63, 52.84])
        tau3 = np.array([0, 17.5, 16.2, 16.4, 16.8, 20.4, 21, 22.1])*1e-12
        e4 = np.array([77.32, 72.95, 68.76, 66.32, 63.06, 57.05, 46.74, 36.83])
        tau4 = np.array([8.3, 8.09, 8.13, 8.05, 8.08, 7.81, 8.02, 8.28])*1e-12
        e5 = np.array([6.21, 6.18, 6.44, 6.47, 6.67, 6.65, 7.31, 7.87])
        tau5 = np.array([0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4])*1e-12
                
        
        einf = einf[Concentration==c]
        e1 = e1[Concentration==c]
        tau1 = tau1[Concentration==c]
        e2 = e2[Concentration==c]
        tau2 = tau2[Concentration==c]
        e3 = e3[Concentration==c]
        tau3 = tau3[Concentration==c]
        e4 = e4[Concentration==c]
        tau4 = tau4[Concentration==c]
        e5 = e5[Concentration==c]
        tau5 = tau5[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + (e1-e2)/(1+1j*omega*tau1)
                           + (e2-e3)/(1+1j*omega*tau2)
                           + (e3-e4)/(1+1j*omega*tau3)
                           + (e4-e5)/(1+1j*omega*tau4)
                           + (e5-einf)/(1+1j*omega*tau5))
        return data

def NaBrAqueous_Wachter(frequency,temperature=25,c=0.052):
    data = {}
    data['minFREQ'] = 0.2e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Wachter, W., Kunz, W., Buchner, R., & Hefter, G. (2005). Is there an anionic Hofmeister effect on water dynamics? Dielectric spectroscopy of aqueous solutions of NaBr, NaI, NaNO3, NaClO4, and NaSCN. The Journal of Physical Chemistry A, 109(39), 8675-8683.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.052
    data['maxConcentration (x)'] = 1.404
    data['Remarks'] = 'c in mol/L of the substance in water. No interpolation is allowed (See the original article). Even though the three-Debye model yields bettera greement with the experimental data, only the Cole-Cole model is used here due to its simpleness.'
    
    Concentration = np.array([0.052, 0.1004, 0.1501, 0.2538, 0.3529, 0.5037, 0.6521, 0.7996, 0.9989, 1.2, 1.404])
    if c not in Concentration:
        print('No interpolation is allowed in this substance at this stage. Available concentrations are:'+str(Concentration)+'. See the cited article.')
        return None
    else:
        einf = np.array([5.6, 5.29, 4.93, 5.19, 5.04, 5.2, 5.11, 5.41, 5.64, 5.84, 5.43])
        e1 = np.array([77.8, 77.2, 76.6, 74.8, 74.4, 72.3, 70.4, 68.9, 66.6, 64.6, 62.4])
        tau1 = np.array([8.23, 8.19, 8.09, 8.07, 7.96, 7.9, 7.75, 7.71, 7.55, 7.46, 7.19])*1e-12
        a1 = np.array([0.006, 0.011, 0.018, 0.019, 0.03, 0.028, 0.026, 0.034, 0.032, 0.031, 0.042])
        Ke = np.array([0.5605, 1.038, 1.534, 2.481, 3.348, 4.704, 5.978, 7.199, 8.748, 10.24, 11.62]) # S/m
        
        einf = einf[Concentration==c]
        e1 = e1[Concentration==c]
        tau1 = tau1[Concentration==c]
        a1 = a1[Concentration==c]
        Ke = Ke[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + (e1-einf)/(1+(1j*omega*tau1)**(1-a1)))
        return data
    
def NaNO3Aqueous_Wachter(frequency,temperature=25,c=0.0502):
    data = {}
    data['minFREQ'] = 0.2e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Wachter, W., Kunz, W., Buchner, R., & Hefter, G. (2005). Is there an anionic Hofmeister effect on water dynamics? Dielectric spectroscopy of aqueous solutions of NaBr, NaI, NaNO3, NaClO4, and NaSCN. The Journal of Physical Chemistry A, 109(39), 8675-8683.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.052
    data['maxConcentration (x)'] = 1.502
    data['Remarks'] = 'c in mol/L of the substance in water. No interpolation is allowed (See the original article). Even though the three-Debye model yields bettera greement with the experimental data, only the Cole-Cole model is used here due to its simpleness.'
    
    Concentration = np.array([0.0502, 0.1508, 0.3508, 0.6505, 1, 1.201, 1.502])
    if c not in Concentration:
        print('No interpolation is allowed in this substance at this stage. Available concentrations are:'+str(Concentration)+'. See the cited article.')
        return None
    else:
        einf = np.array([5.52, 5.47, 4.88, 4.99, 5.3, 4.94, 5.21])
        e1 = np.array([77.9, 76.6, 74.7, 71.7, 68.8, 66.9, 64.8])
        tau1 = np.array([8.21, 8.12, 7.93, 7.74, 7.61, 7.47, 7.41])*1e-12
        a1 = np.array([0.012, 0.017, 0.028, 0.037, 0.044, 0.054, 0.06])
        Ke = np.array([0.524, 1.44, 3.075, 5.231, 7.427, 8.534, 10.08]) # S/m
        
        einf = einf[Concentration==c]
        e1 = e1[Concentration==c]
        tau1 = tau1[Concentration==c]
        a1 = a1[Concentration==c]
        Ke = Ke[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + (e1-einf)/(1+(1j*omega*tau1)**(1-a1)))
        return data
    
def NaIAqueous_Wachter(frequency,temperature=25,c=0.0495):
    data = {}
    data['minFREQ'] = 0.2e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Wachter, W., Kunz, W., Buchner, R., & Hefter, G. (2005). Is there an anionic Hofmeister effect on water dynamics? Dielectric spectroscopy of aqueous solutions of NaBr, NaI, NaNO3, NaClO4, and NaSCN. The Journal of Physical Chemistry A, 109(39), 8675-8683.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.0495
    data['maxConcentration (x)'] = 1.511
    data['Remarks'] = 'c in mol/L of the substance in water. No interpolation is allowed (See the original article). Even though the three-Debye model yields bettera greement with the experimental data, only the Cole-Cole model is used here due to its simpleness.'
    
    Concentration = np.array([0.0495, 0.148, 0.3478, 0.6528, 1.002, 1.511])
    if c not in Concentration:
        print('No interpolation is allowed in this substance at this stage. Available concentrations are:'+str(Concentration)+'. See the cited article.')
        return None
    else:
        einf = np.array([5.28, 4.81, 4.79, 5.34, 5.05, 6.47])
        e1 = np.array([77.7, 76.6, 74.2, 70, 65.5, 60.9])
        tau1 = np.array([8.13, 8, 7.86, 7.59, 7.19, 7])*1e-12
        a1 = np.array([0.01, 0.021, 0.029, 0.025, 0.039, 0.041])
        Ke = np.array([0.638, 1.521, 3.387, 6.025, 8.765, 12.58]) # S/m
        
        einf = einf[Concentration==c]
        e1 = e1[Concentration==c]
        tau1 = tau1[Concentration==c]
        a1 = a1[Concentration==c]
        Ke = Ke[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + (e1-einf)/(1+(1j*omega*tau1)**(1-a1)))
        return data

def NaClO4Aqueous_Wachter(frequency,temperature=25,c=0.045):
    data = {}
    data['minFREQ'] = 0.2e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Wachter, W., Kunz, W., Buchner, R., & Hefter, G. (2005). Is there an anionic Hofmeister effect on water dynamics? Dielectric spectroscopy of aqueous solutions of NaBr, NaI, NaNO3, NaClO4, and NaSCN. The Journal of Physical Chemistry A, 109(39), 8675-8683.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.045
    data['maxConcentration (x)'] = 1.5
    data['Remarks'] = 'c in mol/L of the substance in water. No interpolation is allowed (See the original article). Even though the three-Debye model yields bettera greement with the experimental data, only the Cole-Cole model is used here due to its simpleness.'
    
    Concentration = np.array([0.045, 0.135, 0.3, 0.6, 0.9, 1.5])
    if c not in Concentration:
        print('No interpolation is allowed in this substance at this stage. Available concentrations are:'+str(Concentration)+'. See the cited article.')
        return None
    else:
        einf = np.array([5.51, 4.92, 4.89, 4.76, 4.82, 5.21])
        e1 = np.array([77.8, 76.6, 74.8, 71.1, 67.9, 62])
        tau1 = np.array([8.15, 7.97, 7.79, 7.45, 7.26, 6.88])*1e-12
        a1 = np.array([0.005, 0.019, 0.029, 0.038, 0.045, 0.054])
        Ke = np.array([0.461, 1.29, 2.58, 4.92, 6.94, 10.4]) # S/m
        
        einf = einf[Concentration==c]
        e1 = e1[Concentration==c]
        tau1 = tau1[Concentration==c]
        a1 = a1[Concentration==c]
        Ke = Ke[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + (e1-einf)/(1+(1j*omega*tau1)**(1-a1)))
        return data

def NaSCNAqueous_Wachter(frequency,temperature=25,c=0.0497):
    data = {}
    data['minFREQ'] = 0.2e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Wachter, W., Kunz, W., Buchner, R., & Hefter, G. (2005). Is there an anionic Hofmeister effect on water dynamics? Dielectric spectroscopy of aqueous solutions of NaBr, NaI, NaNO3, NaClO4, and NaSCN. The Journal of Physical Chemistry A, 109(39), 8675-8683.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.0497
    data['maxConcentration (x)'] = 1.7
    data['Remarks'] = 'c in mol/L of the substance in water. No interpolation is allowed (See the original article). Even though the three-Debye model yields bettera greement with the experimental data, only the Cole-Cole model is used here due to its simpleness.'
    
    Concentration = np.array([0.0497, 0.1493, 0.3972, 0.8073, 1.231, 1.7])
    if c not in Concentration:
        print('No interpolation is allowed in this substance at this stage. Available concentrations are:'+str(Concentration)+'. See the cited article.')
        return None
    else:
        einf = np.array([5.54, 5.14, 5.14, 4.93, 5.32, 5.51])
        e1 = np.array([77.6, 76.4, 73.6, 68.6, 64, 59.6])
        tau1 = np.array([8.12, 8.03, 7.82, 7.42, 7.26, 6.93])*1e-12
        a1 = np.array([0.004, 0.014, 0.025, 0.034, 0.037, 0.052])
        Ke = np.array([0.509, 1.39, 3.4, 6.54, 9.13, 11.7]) # S/m
        
        einf = einf[Concentration==c]
        e1 = e1[Concentration==c]
        tau1 = tau1[Concentration==c]
        a1 = a1[Concentration==c]
        Ke = Ke[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + (e1-einf)/(1+(1j*omega*tau1)**(1-a1)))
        return data

def NiSO4Aqueous_Chen(frequency,temperature=25,c=0.02472):
    data = {}
    data['minFREQ'] = 0.2e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Chen, T., Hefter, G., & Buchner, R. (2005). Ion association and hydration in aqueous solutions of nickel (II) and cobalt (II) sulfate. Journal of Solution Chemistry, 34(9), 1045-1066.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.02472
    data['maxConcentration (x)'] = 1.402
    data['Remarks'] = 'c in mol/L of the substance in water. No interpolation is allowed (See the original article).'
    
    Concentration = np.array([0.02472, 0.04993, 0.1, 0.1502, 0.1998, 0.2503, 0.3702, 0.499, 0.5003, 0.6465, 0.6703, 0.8014, 1, 1.011, 1.2, 1.402])
    if c not in Concentration:
        print('No interpolation is allowed in this substance at this stage. Available concentrations are:'+str(Concentration)+'. See the cited article.')
        return None
    else:
        einf = np.array([5.86, 6.17, 5.71, 4.44, 6.8, 6.21, 5.06, 6.16, 5.84, 4.49, 7.07, 4.28, 5.26, 5.4, 6.31, 4.16])
        e1 = np.array([81.39, 83.88, 86.68, 87.94, 89.23, 89.24, 89.45, 88.47, 87.96, 87.12, 87.55, 85.19, 81.74, 81.86, 80.22, 78.19])
        tau1 = np.array([566, 550, 455, 504, 489, 525, 665, 518, 535, 500, 488, 500, 436, 508, 439, 521])*1e-12
        e2 = np.array([80.09, 80.89, 82.78, 84.63, 84.77, 85.85, 86.79, 85.73, 86.08, 83.95, 83.5, 82.85, 80.45, 80.49, 77.61, 74.96])
        tau2 = np.array([261, 185, 176, 166, 161, 150, 145, 136, 139, 130, 128, 127, 132, 133, 129, 129])*1e-12
        e3 = np.array([77.79, 77.05, 76.07, 75.06, 74.21, 73.55, 72, 70.35, 70.47, 68.71, 68.08, 66.89, 65.01, 65.04, 63.1, 60.94])
        tau3 = np.array([0, 24.5, 22, 20.6, 19.9, 22.5, 19.7, 22.4, 23.5, 22, 22, 21.5, 24.3, 24.8, 24.3, 22.9])*1e-12
        e4 = np.array([77.79, 76.63, 75.4, 73.39, 70.44, 69.62, 64.96, 63.66, 63.71, 59.78, 58.78, 56.34, 54.05, 54.12, 50.4, 46.25])
        tau4 = np.array([8.27, 8.23, 8.15, 8.14, 8.03, 7.98, 7.83, 7.81, 7.81, 7.76, 7.74, 7.7, 7.68, 7.68, 7.66, 7.67])*1e-12
        e5 = np.array([5.86, 6.17, 5.71, 6.45, 6.8, 6.21, 6.29, 6.16, 5.84, 6.26, 7.07, 6.64, 5.26, 5.4, 6.31, 7.02])
        tau5 = np.array([0, 0, 0, 0.4, 0, 0, 0.4, 0, 0, 0.4, 0, 0.4, 0, 0, 0, 0.4])*1e-12
        Ke = np.array([0.289, 0.504, 0.835, 1.174, 1.46, 1.754, 2.359, 2.881, 2.879, 3.401, 3.446, 3.928, 4.369, 4.372, 4.814, 5.13])
        
        einf = einf[Concentration==c]
        e1 = e1[Concentration==c]
        tau1 = tau1[Concentration==c]
        e2 = e2[Concentration==c]
        tau2 = tau2[Concentration==c]
        e3 = e3[Concentration==c]
        tau3 = tau3[Concentration==c]
        e4 = e4[Concentration==c]
        tau4 = tau4[Concentration==c]
        e5 = e5[Concentration==c]
        tau5 = tau5[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + (e1-e2)/(1+1j*omega*tau1)
                           + (e2-e3)/(1+1j*omega*tau2)
                           + (e3-e4)/(1+1j*omega*tau3)
                           + (e4-e5)/(1+1j*omega*tau4)
                           + (e5-einf)/(1+1j*omega*tau5))
        return data

def CuSO4Aqueous_Akilan(frequency,temperature=25,c=0.0201):
    data = {}
    data['minFREQ'] = 0.2e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Akilan, C., Hefter, G., Rohman, N., & Buchner, R. (2006). Ion association and hydration in aqueous solutions of copper (II) sulfate from 5 to 65 C by dielectric spectroscopy. The Journal of Physical Chemistry B, 110(30), 14961-14970.'
    data['minTemperature (C)'] = 5
    data['maxTemperature (C)'] = 65    
    data['minConcentration (x)'] = 0.0201
    data['maxConcentration (x)'] = 1.4183
    data['Remarks'] = 'c in mol/kg of the substance in water (molality). No interpolation for temperature and concentration is allowed (See the original article).'
    
    Temperature = np.array([5,25,45,65])
    
    if temperature not in Temperature:
        print('No interpolation is allowed in this substance at this stage. Available temperatures are:'+str(Temperature)+'. See the cited article.')
        return None
    else:
        if temperature == 5:
            Concentration = np.array([0.0201, 0.0501, 0.0703, 0.1003, 0.1506, 0.2005, 0.2507, 0.3517, 0.4013, 0.5019, 0.6539, 0.8048, 1.008])
            if c not in Concentration:
                print('No interpolation is allowed in this substance at this stage. Available concentrations  at '+str(temperature)+' C are: '+str(Concentration)+'. See the cited article.')
                return None
            else:
                einf = np.array([6.64, 7, 6.47, 6.5, 6.87, 6.13, 6.56, 7.09, 5.04, 6.87, 5.98, 7.62, 3.26])
                e1 = np.array([87.31, 90.09, 91.3, 93.64, 94.04, 94.58, 94.87, 95.07, 94.4, 92.98, 91.33, 89.02, 85.21])
                tau1 = np.array([520, 489, 474, 604, 443, 416, 420, 444, 370, 0, 0, 0, 0])*1e-12
                e2 = np.array([85.36, 86.74, 88.53, 88.72, 89.65, 91.28, 92.07, 93.98, 93.02, 92.98, 91.33, 89.02, 85.21])
                tau2 = np.array([220, 220, 275, 220, 222, 219, 215, 218, 200, 189, 181, 172, 157])*1e-12
                e3 = np.array([84.55, 84.07, 83.7, 83.07, 81.74, 80.75, 80, 78.41, 77.25, 75.8, 73.21, 70.97, 67.97])
                tau3 = np.array([0, 0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30.5, 30])*1e-12
                e4 = np.array([84.55, 84.07, 81.33, 79.08, 77.9, 76.41, 74.77, 70.84, 70.76, 67.99, 65.79, 63.99, 59.53])
                tau4 = np.array([14.3, 14.7, 14.3, 14.2, 14.1, 14, 14, 14, 14, 14.1, 14.3, 14.5, 14.8])*1e-12
                e5 = np.array([6.64, 7, 6.47, 6.5, 6.87, 6.59, 6.56, 7.09, 6.96, 6.87, 7.29, 7.62, 7.93])
                tau5 = np.array([0, 0, 0, 0, 0, 0.5, 0, 0, 0.5, 0, 0.5, 0, 0.5])*1e-12
                Ke = np.array([0.182, 0.36, 0.464, 0.599, 0.88, 1.05, 1.24, 1.61, 1.79, 2.13, 2.53, 2.89, 3.32])
        elif temperature == 25:
            Concentration = np.array([0.0201, 0.0501, 0.0702, 0.1003, 0.1504, 0.2005, 0.2507, 0.3511, 0.4013, 0.5019, 0.6531, 0.8048, 1.0079, 1.2123, 1.418])
            if c not in Concentration:
                print('No interpolation is allowed in this substance at this stage. Available concentrations at '+str(temperature)+' C are: '+str(Concentration)+'. See the cited article.')
                return None
            else:
                einf = np.array([6.06, 5.8, 6.71, 6.2, 7.22, 2.31, 6.77, 6.83, 2.03, 6.79, 2.66, 7.32, 4.2, 7.69, 7.95])
                e1 = np.array([81.11, 83.16, 83.91, 87.59, 87.83, 89.07, 88.93, 90.08, 89.41, 87.94, 86.43, 85.06, 81.44, 78.77, 75.82])
                tau1 = np.array([455, 425, 356, 380, 241, 407, 423, 346, 405, 288, 380, 413, 0, 0, 0])*1e-12
                e2 = np.array([78.6, 80.56, 82.38, 81.35, 81.79, 85.79, 84.4, 86.34, 87.67, 86.27, 85.49, 83.24, 81.44, 78.77, 75.82])
                tau2 = np.array([96.1, 171, 188, 150, 165, 153, 127, 136, 144, 132, 127, 120, 122, 124, 122])*1e-12
                e3 = np.array([77.68, 77.1, 77.05, 76.23, 76.18, 75.05, 72.94, 72.81, 72.46, 71.06, 69.27, 67.82, 65.51, 64.15, 62.26])
                tau3 = np.array([0, 0, 25, 30.8, 25, 25, 24.7, 23.5, 25.4, 25.6, 24, 27, 24.7, 26.3, 26.5])*1e-12
                e4 = np.array([77.68, 77.1, 75.48, 75.39, 71.99, 71.75, 71.21, 67.82, 66.91, 65.05, 61.4, 59.5, 55.38, 51.79, 49.08])
                tau4 = np.array([8.3, 8.28, 8.27, 8.25, 8.23, 8.29, 8.2, 8.19, 8.13, 8.18, 8.14, 8.17, 8.15, 8.17, 8.16])*1e-12
                e5 = np.array([6.06, 5.8, 6.71, 6.2, 7.22, 6.74, 6.77, 6.83, 6.97, 6.79, 6.99, 7.32, 7.47, 7.69, 7.95])
                tau5 = np.array([0, 0, 0, 0, 0, 0.3, 0, 0, 0.3, 0, 0.3, 0, 0.3, 0, 0])*1e-12
                Ke = np.array([0.25, 0.49, 0.62, 0.827, 1.18, 1.47, 1.79, 2.26, 2.51, 2.93, 3.51, 4.07, 4.63, 5.13, 5.57])
        elif temperature == 45:
            Concentration = np.array([0.0201, 0.0508, 0.0703, 0.1003, 0.1506, 0.2005, 0.2515, 0.3517, 0.4013, 0.4989, 0.6539, 0.8048, 1.0078, 1.2123, 1.4183])
            if c not in Concentration:
                print('No interpolation is allowed in this substance at this stage. Available concentrations at '+str(temperature)+' C are: '+str(Concentration)+'. See the cited article.')
                return None
            else:
                einf = np.array([5.12, 5.5, 8.59, 4.96, 4.95, 5.35, 5.21, 6.71, 5.41, 5.49, 6.26, 6.88, 6.63, 6, 6.15])
                e1 = np.array([73.99, 76.35, 77.57, 78.97, 81.28, 80.14, 81.31, 82.04, 81.81, 81.51, 80.46, 79.27, 77.03, 74.9, 71.5])
                tau1 = np.array([339, 301, 304, 425, 456, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])*1e-12
                e2 = np.array([72.01, 72.86, 74.22, 77.07, 78.75, 80.14, 81.31, 82.04, 81.81, 81.51, 80.46, 79.27, 77.03, 74.9, 71.5])
                tau2 = np.array([82.1, 86.1, 104, 131, 122, 123, 127, 121, 116, 120, 104, 108, 110, 120, 115])*1e-12
                e3 = np.array([70.98, 69.99, 69.9, 70.12, 69.72, 69.33, 69.63, 68.89, 68.48, 68.5, 66.82, 66.8, 65.42, 65.49, 62.89])
                tau3 = np.array([0, 0, 0, 24.1, 22.9, 29.3, 30.1, 25.4, 26.4, 29.7, 27, 29.6, 30, 33.1, 31.9])*1e-12
                e4 = np.array([70.98, 69.99, 69.9, 68.46, 67.08, 66.44, 65.19, 62.29, 61.49, 60.15, 57.76, 55.89, 53.21, 51.43, 48.53])
                tau4 = np.array([5.27, 5.23, 5.51, 5.18, 5.14, 5.15, 5.1, 5.08, 5.02, 5.08, 5.13, 5.14, 5.18, 5.24, 5.28])*1e-12
                e5 = np.array([5.12, 5.5, 8.59, 4.96, 4.95, 5.35, 5.21, 6.71, 5.41, 5.49, 6.26, 6.88, 6.63, 6, 6.15])
                tau5 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])*1e-12
                Ke = np.array([0.317, 0.639, 0.816, 1.1, 1.49, 1.87, 2.26, 2.84, 3.17, 3.77, 4.48, 5.3, 5.99, 6.86, 7.44])
        elif temperature == 65:
            Concentration = np.array([0.0201, 0.0508, 0.0703, 0.1003, 0.1506, 0.2005, 0.2515, 0.3517, 0.4013, 0.4989, 0.6539, 0.8048, 1.0078, 1.2123, 1.4183])
            if c not in Concentration:
                print('No interpolation is allowed in this substance at this stage. Available concentrations at '+str(temperature)+' C are: '+str(Concentration)+'. See the cited article.')
                return None
            else:
                einf = np.array([6.91, 4.07, 5.21, 5.47, 6.34, 5.3, 5.85, 6.87, 5.16, 6.68, 7.65, 7.64, 7.14, 5.31, 7.08])
                e1 = np.array([68.03, 71.1, 71.73, 73.8, 74.47, 75.99, 76.68, 76.9, 77.4, 76.97, 76.68, 75.13, 73.08, 71.5, 69.41])
                tau1 = np.array([465, 438, 429, 613, 293, 407, 458, 0, 0, 0, 0, 0, 0, 0, 0])*1e-12
                e2 = np.array([65.95, 68, 69.75, 71.64, 73.19, 73.09, 75.15, 76.9, 77.4, 76.97, 76.68, 75.13, 73.08, 71.5, 69.41])
                tau2 = np.array([65.8, 78.5, 86, 109, 101, 91.9, 95.6, 137, 103, 103, 91, 95.1, 99.7, 121, 129])*1e-12
                e3 = np.array([64.25, 64.06, 64.08, 65.27, 65.11, 64.6, 65.66, 67.8, 66.64, 66.89, 65.66, 65.6, 64.96, 64.9, 64.49])
                tau3 = np.array([0, 0, 0, 30.7, 21.8, 30.6, 24.4, 35.4, 29.4, 28.3, 25.8, 27.9, 27.8, 30.3, 30.1])*1e-12
                e4 = np.array([64.25, 64.06, 64.08, 62.99, 60.99, 60.63, 59.52, 57.85, 58.07, 55.95, 53.56, 51.75, 48.61, 47.71, 45.94])
                tau4 = np.array([3.7, 3.7, 3.7, 3.7, 3.7, 3.57, 3.7, 3.7, 3.7, 3.7, 3.85, 3.7, 3.7, 3.7, 3.7])*1e-12
                e5 = np.array([6.91, 4.07, 5.21, 5.47, 6.34, 5.3, 5.85, 6.87, 5.16, 6.68, 7.65, 7.64, 7.14, 5.31, 7.08])
                tau5 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])*1e-12
                Ke = np.array([0.408, 0.772, 1.01, 1.31, 1.81, 2.3, 2.74, 3.55, 3.87, 4.66, 5.57, 6.34, 7.66, 8.31, 9.35])
                                
        einf = einf[Concentration==c]
        e1 = e1[Concentration==c]
        tau1 = tau1[Concentration==c]
        e2 = e2[Concentration==c]
        tau2 = tau2[Concentration==c]
        e3 = e3[Concentration==c]
        tau3 = tau3[Concentration==c]
        e4 = e4[Concentration==c]
        tau4 = tau4[Concentration==c]
        e5 = e5[Concentration==c]
        tau5 = tau5[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + (e1-e2)/(1+1j*omega*tau1)
                           + (e2-e3)/(1+1j*omega*tau2)
                           + (e3-e4)/(1+1j*omega*tau3)
                           + (e4-e5)/(1+1j*omega*tau4)
                           + (e5-einf)/(1+1j*omega*tau5))
        return data
    
def LaCl3Aqueous_Friesen(frequency,temperature=25,c=0.03972):
    data = {}
    data['minFREQ'] = 0.05e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Friesen, S., Krickl, S., Luger, M., Nazet, A., Hefter, G., & Buchner, R. (2018). Hydration and ion association of La 3+ and Eu 3+ salts in aqueous solution. Physical Chemistry Chemical Physics, 20(13), 8812-8821.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.03972
    data['maxConcentration (x)'] = 0.9664
    data['Remarks'] = 'c in mol/L of the substance in water. No interpolation is allowed (See the original article).'
    
    Concentration = np.array([0.03972, 0.06953, 0.1486, 0.2957, 0.4917, 0.732, 0.9664])
    if c not in Concentration:
        print('No interpolation is allowed in this substance at this stage. Available concentrations are:'+str(Concentration)+'. See the cited article.')
        return None
    else:
        einf = np.array([5.43, 5.37, 5.12, 5.22, 6.02, 6.25, 6.58])
        de1 = np.array([2.89, 3.91, 5.08, 4.82, 5.53, 5.14, 5.46])
        tau1 = np.array([268, 171, 221, 279, 436, 366, 387])*1e-12
        de2 = np.array([0, 0, 0.84, 2.15, 5.5, 7.42, 11.33])
        tau2 = np.array([0, 0, 30, 30, 29.8, 22.1, 19.6])*1e-12
        de3 = np.array([70.54, 69.57, 66.69, 61.19, 52.53, 43.97, 35.07])
        tau3 = np.array([8.26, 8.2, 8.14, 8.01, 7.66, 7.23, 6.49])*1e-12
        a3 = np.array([0.011, 0.012, 0.024, 0.036, 0.028, 0.033, 0.031])
        
        einf = einf[Concentration==c]
        de1 = de1[Concentration==c]
        tau1 = tau1[Concentration==c]
        de2 = de2[Concentration==c]
        tau2 = tau2[Concentration==c]
        de3 = de3[Concentration==c]
        tau3 = tau3[Concentration==c]
        a3 = a3[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + de1/(1+1j*omega*tau1)
                           + de2/(1+1j*omega*tau2)
                           + de3/(1+(1j*omega*tau3)**(1-a3)))
        return data
    
def LaNO33Aqueous_Friesen(frequency,temperature=25,c=0.004979):
    data = {}
    data['minFREQ'] = 0.05e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Friesen, S., Krickl, S., Luger, M., Nazet, A., Hefter, G., & Buchner, R. (2018). Hydration and ion association of La 3+ and Eu 3+ salts in aqueous solution. Physical Chemistry Chemical Physics, 20(13), 8812-8821.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.004979
    data['maxConcentration (x)'] = 0.838
    data['Remarks'] = 'c in mol/L of the substance in water. No interpolation is allowed (See the original article).'
    
    Concentration = np.array([0.004979, 0.01381, 0.0279, 0.05179, 0.1002, 0.2083, 0.2611, 0.4327, 0.605, 0.778, 0.838])
    if c not in Concentration:
        print('No interpolation is allowed in this substance at this stage. Available concentrations are:'+str(Concentration)+'. See the cited article.')
        return None
    else:
        einf = np.array([3.78, 3.52, 3.52, 2.82, 2.47, 3.52, 4.23, 2.36, 3.51, 3.84, 3.48])
        de1 = np.array([0.97, 2.08, 2.11, 0.71, 1.79, 1.28, 0.76, 0, 0, 0, 0])
        tau1 = np.array([1298, 1103, 999, 618, 641, 369, 415, 0, 0, 0, 0])*1e-12
        de2 = np.array([0.38, 0.99, 0.86, 3.31, 4.31, 5.35, 6.26, 7.51, 9.87, 12.55, 13.37])
        tau2 = np.array([479, 283, 258, 240, 145, 122, 127, 150, 158, 169, 173])*1e-12
        de3 = np.array([0, 0, 0, 0, 0, 0, 0.49, 2.96, 4.52, 5.08, 8.61])
        tau3 = np.array([0, 0, 0, 0, 0, 0, 37.2, 31.1, 24.4, 29.7, 21.4])*1e-12
        de4 = np.array([70.85, 69.94, 69.94, 68.54, 65.44, 61.56, 58.95, 53.53, 47.83, 43.82, 39.68])
        tau4 = np.array([8.38, 8.38, 8.38, 8.25, 8.27, 8.09, 8.07, 7.69, 7.35, 7.2, 6.77])*1e-12
        de5 = np.array([2.4, 3.05, 3.05, 3.71, 4.76, 4.05, 4.08, 5.58, 4.41, 4.57, 4.46])
        tau5 = np.array([0.86, 1.06, 1.06, 0.77, 0.88, 1.09, 1.47, 0.66, 0.64, 0.68, 0.58])*1e-12
        
        einf = einf[Concentration==c]
        de1 = de1[Concentration==c]
        tau1 = tau1[Concentration==c]
        de2 = de2[Concentration==c]
        tau2 = tau2[Concentration==c]
        de3 = de3[Concentration==c]
        tau3 = tau3[Concentration==c]
        de4 = de4[Concentration==c]
        tau4 = tau4[Concentration==c]
        de5 = de5[Concentration==c]
        tau5 = tau5[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + de1/(1+1j*omega*tau1)
                           + de2/(1+1j*omega*tau2)
                           + de3/(1+(1j*omega*tau3))
                           + de4/(1+(1j*omega*tau4))
                           + de5/(1+(1j*omega*tau5)) )
        
        return data

def La2SO43Aqueous_Friesen(frequency,temperature=25,c=0.00356):
    data = {}
    data['minFREQ'] = 0.05e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Friesen, S., Krickl, S., Luger, M., Nazet, A., Hefter, G., & Buchner, R. (2018). Hydration and ion association of La 3+ and Eu 3+ salts in aqueous solution. Physical Chemistry Chemical Physics, 20(13), 8812-8821.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.00356
    data['maxConcentration (x)'] = 0.0357
    data['Remarks'] = 'c in mol/L of the substance in water. No interpolation is allowed (See the original article).'
    
    Concentration = np.array([0.00356, 0.00705, 0.0106, 0.018, 0.0212, 0.0283, 0.0357])
    if c not in Concentration:
        print('No interpolation is allowed in this substance at this stage. Available concentrations are:'+str(Concentration)+'. See the cited article.')
        return None
    else:
        einf = np.array([5.42, 5.13, 5.51, 5.46, 5.55, 5.5, 5.6])
        de1 = np.array([0.94, 1.8, 2.33, 3.32, 3.61, 4.31, 5.13])
        tau1 = np.array([613, 539, 697, 666, 800, 649, 602])*1e-12
        de2 = np.array([1.01, 1.92, 3.23, 4.98, 6.23, 7.41, 8.58])
        tau2 = np.array([146, 131, 142, 141, 148, 138, 132])*1e-12
        de3 = np.array([72.06, 71.62, 71.15, 70.54, 70.06, 69.38, 68.84])
        tau3 = np.array([8.19, 8.14, 8.17, 8.15, 8.18, 8.16, 8.15])*1e-12
        
        einf = einf[Concentration==c]
        de1 = de1[Concentration==c]
        tau1 = tau1[Concentration==c]
        de2 = de2[Concentration==c]
        tau2 = tau2[Concentration==c]
        de3 = de3[Concentration==c]
        tau3 = tau3[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + de1/(1+1j*omega*tau1)
                           + de2/(1+1j*omega*tau2)
                           + de3/(1+(1j*omega*tau3)) )
        
        return data

def EuNO33Aqueous_Friesen(frequency,temperature=25,c=0.01551):
    data = {}
    data['minFREQ'] = 0.05e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Friesen, S., Krickl, S., Luger, M., Nazet, A., Hefter, G., & Buchner, R. (2018). Hydration and ion association of La 3+ and Eu 3+ salts in aqueous solution. Physical Chemistry Chemical Physics, 20(13), 8812-8821.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.01551
    data['maxConcentration (x)'] = 0.2276
    data['Remarks'] = 'c in mol/L of the substance in water. No interpolation is allowed (See the original article).'
    
    Concentration = np.array([0.01551, 0.03083, 0.05379, 0.1147, 0.2276])
    if c not in Concentration:
        print('No interpolation is allowed in this substance at this stage. Available concentrations are:'+str(Concentration)+'. See the cited article.')
        return None
    else:
        einf = np.array([3.52, 3.52, 3.52, 3.52, 3.52])
        de1 = np.array([2.45, 3.19, 4.19, 6.07, 8.39])
        tau1 = np.array([161, 143, 141, 126, 113])*1e-12
        de2 = np.array([66.45, 65.47, 64.66, 61.59, 57.6])
        tau2 = np.array([8.4, 8.31, 8.29, 8.23, 8.03])*1e-12
        de3 = np.array([2.88, 2.74, 2.79, 3.2, 3.19])
        tau3 = np.array([0.28, 0.28, 0.28, 0.28, 0.28])*1e-12
        
        einf = einf[Concentration==c]
        de1 = de1[Concentration==c]
        tau1 = tau1[Concentration==c]
        de2 = de2[Concentration==c]
        tau2 = tau2[Concentration==c]
        de3 = de3[Concentration==c]
        tau3 = tau3[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + de1/(1+1j*omega*tau1)
                           + de2/(1+1j*omega*tau2)
                           + de3/(1+(1j*omega*tau3)) )
        
        return data

def Eu2SO43Aqueous_Friesen(frequency,temperature=25,c=0.004001):
    data = {}
    data['minFREQ'] = 0.05e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Friesen, S., Krickl, S., Luger, M., Nazet, A., Hefter, G., & Buchner, R. (2018). Hydration and ion association of La 3+ and Eu 3+ salts in aqueous solution. Physical Chemistry Chemical Physics, 20(13), 8812-8821.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0.004001
    data['maxConcentration (x)'] = 0.0249
    data['Remarks'] = 'c in mol/L of the substance in water. No interpolation is allowed (See the original article).'
    
    Concentration = np.array([0.004001, 0.007966, 0.01194, 0.01492, 0.01787, 0.02086, 0.0249])
    if c not in Concentration:
        print('No interpolation is allowed in this substance at this stage. Available concentrations are:'+str(Concentration)+'. See the cited article.')
        return None
    else:
        einf = np.array([5.66, 5.69, 5.61, 5.67, 5.76, 5.62, 5.41])
        de1 = np.array([1.4, 2.37, 3.02, 3.95, 4.25, 4.91, 5.37])
        tau1 = np.array([540, 503, 510, 418, 433, 409, 382])*1e-12
        de2 = np.array([1.2, 2.15, 3.16, 3.44, 4.15, 4.6, 5.07])
        tau2 = np.array([90, 90.5, 99.3, 84.6, 94.7, 88.2, 89.9])*1e-12
        de3 = np.array([71.41, 70.95, 70.07, 70.12, 69.94, 69.66, 68.85])
        tau3 = np.array([8.22, 8.19, 8.13, 8.11, 8.13, 8.07, 8])*1e-12
        
        einf = einf[Concentration==c]
        de1 = de1[Concentration==c]
        tau1 = tau1[Concentration==c]
        de2 = de2[Concentration==c]
        tau2 = tau2[Concentration==c]
        de3 = de3[Concentration==c]
        tau3 = tau3[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + de1/(1+1j*omega*tau1)
                           + de2/(1+1j*omega*tau2)
                           + de3/(1+(1j*omega*tau3)) )
        
        return data

def NdCl3Aqueous_Yoon(frequency,temperature=23,c=0.027):
    data = {}
    data['minFREQ'] = 0.2e9
    data['maxFREQ'] = 40e9
    data['Citation'] = 'Yoon, T. J., Vigil, M. J., Raby, E. Y., Singh, R. P., Maerzke, K. A., Currier, R. P., & Findikoglu, A. T. (2020). Dielectric relaxation of neodymium chloride in water and in methanol. Journal of Molecular Liquids, 308, 112981.'
    data['minTemperature (C)'] = 23
    data['maxTemperature (C)'] = 23    
    data['minConcentration (x)'] = 0.027
    data['maxConcentration (x)'] = 0.789
    data['Remarks'] = 'c in mol/kg of the substance in water. No interpolation is allowed (See the original article).'
    
    Concentration = np.array([0.027, 0.052, 0.077, 0.104, 0.129, 0.232, 0.473, 0.789])
    if c not in Concentration:
        print('No interpolation is allowed in this substance at this stage. Available concentrations are:'+str(Concentration)+'. See the cited article.')
        return None
    else:
        einf = np.array([5.259, 5.255, 5.34, 5.46, 5.42, 6.88, 7.45, 6.36])
        de1 = np.array([73.23, 70.76, 68.74, 66.92, 65.34, 59.53, 50.6, 43.66])
        tau1 = np.array([8.924, 8.778, 8.675, 8.57, 8.455, 8.025, 8.11, 7.97])*1e-12
        de2 = np.array([0, 1.22, 1.8, 2.1, 2.48, 4.543, 5.9, 4.79])
        tau2 = np.array([0, 10.9, 19.5, 22.8, 24.5, 26.9, 30.3, 35.5])*1e-12
        de3 = np.array([2.17, 2.45, 3.32, 3.08, 2.86, 2.96, 1.45, 0.67])
        tau3 = np.array([203, 202, 180, 119, 84.2, 139, 146, 158])*1e-12
        Ke = np.array([0.7711, 1.4153, 2.033, 2.6713, 3.2296, 5.4152, 8.9818, 12.5059]) # S/m
        
        einf = einf[Concentration==c]
        de1 = de1[Concentration==c]
        tau1 = tau1[Concentration==c]
        de2 = de2[Concentration==c]
        tau2 = tau2[Concentration==c]
        de3 = de3[Concentration==c]
        tau3 = tau3[Concentration==c]
        
        omega = 2*np.pi*frequency
        data['epsilon'] = (einf + de1/(1+1j*omega*tau1)
                           + de2/(1+1j*omega*tau2)
                           + de3/(1+(1j*omega*tau3)) )
        
        return data
    
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Binary mixture (non-electrolytes)
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
def GlycerolEthanol_Kaatze(frequency,temperature=25,c=1.0):
    data = {}
    data['minFREQ'] = 400e3
    data['maxFREQ'] = 3e9
    data['Citation'] = 'Kaatze, U. (2012). Hydrogen network fluctuations: Dielectric spectra of glycerol–ethanol mixtures. Chemical Physics, 403, 74-80.'
    data['minTemperature (C)'] = 10
    data['maxTemperature (C)'] = 90
    data['frequency'] = frequency
    data['Temperature (C)'] = temperature
    data['Remarks'] = 'Concentration is in mole fraction of glycerol in ethanol. The high-frequency limit of the dielectric constant epsinf is fixed in this work.'
        
    einf = 4.3 # Fixed.
    omega = 2*np.pi*frequency
    Temperature = np.array([10,20,25,30,40,50])
    if c >=0 and c<=1.0:
        Concentration = np.array([0, 0.25, 0.5, 0.75, 0.9, 1])
        Concentration = np.tile(Concentration,len(Temperature))
        Temperature = np.repeat(Temperature,int(len(Concentration)/6))
        epss = np.array([26.7, 29.4, 32.6, 38.7, 43.2, 48.1, 25.2, 27.4, 31, 37, 41.9, 45.5, 24.5, 26.7, 30.2, 35.3, 38.2, 43, 23.9, 25.9, 29.4, 34.6, 37.9, 42.1, 22.4, 24.5, 27.7, 32.6, 36.2, 39.7, 21.2, 22.8, 26, 31, 33.9, 37.8])
        tau1 = np.array([233, 207, 372, 926, 3313, 4100, 184, 178, 258, 583, 1573, 1700, 162, 169, 241, 540, 1026, 1400, 143, 150, 203, 386, 835, 1240, 105, 112, 142, 262, 348, 740, 82, 91, 119, 206, 286, 495])*1e-12
        b1 = np.array([1, 0.98, 0.78, 0.72, 0.66, 0.65, 1, 0.89, 0.81, 0.72, 0.63, 0.69, 1, 0.88, 0.8, 0.65, 0.7, 0.69, 1, 0.9, 0.79, 0.74, 0.66, 0.7, 1, 0.9, 0.82, 0.73, 0.72, 0.69, 1, 0.82, 0.75, 0.77, 0.71, 0.64])
        
        f1 = interp2d(Temperature,Concentration,epss)
        f2 = interp2d(Temperature,Concentration,tau1)
        f3 = interp2d(Temperature,Concentration,b1)
        
        epss = f1(temperature,c)
        tau1 = f2(temperature,c)
        b1 = f3(temperature,c)
        
        de1 = epss - einf
        
        data['epsilon'] = einf + de1/(1+(1j*omega*tau1)**b1)
        return data
    else:
        print('out of range. See the remarks.')
        return None

def EthanolAqueous_Sato(frequency,temperature=25,c=1.0):
    data = {}
    data['minFREQ'] = 0.1e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Sato, T., & Buchner, R. (2004). Dielectric relaxation processes in ethanol/water mixtures. The Journal of Physical Chemistry A, 108(23), 5007-5015.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0
    data['maxConcentration (x)'] = 1
    data['Remarks'] = 'c in mole fraction of the substance in water. The mole fraction between 0.18 and 0.3 should be avoided.'
    
    if c >=0 and c<=0.18:
        Concentration = np.array([0, 0.04, 0.08, 0.11, 0.18])
        einf = np.array([3.96, 3.7, 3.6, 3.55, 3.4])
        de1 = np.array([72.15, 67.58, 62.84, 59.34, 51.76])
        tau1 = np.array([8.32, 11.97, 15.85, 18.97, 25.48])*1e-12
        b1 = np.array([1, 0.978, 0.958, 0.954, 0.954])
        de2 = np.array([0,0,0,0,0])
        tau2 = np.array([0,0,0,0,0])*1e-12
        de3 = np.array([2.14, 1.92, 1.86, 1.9, 2.16])
        tau3 = np.array([0.39, 0.56, 1.12, 0.65, 1.46])*1e-12
    elif c>=0.3 and c<=1:
        Concentration = np.array([0.3, 0.5, 0.7, 0.9, 1])
        einf = np.array([3.14, 2.87, 2.7, 2.68, 2.6])
        de1 = np.array([42.45, 31.07, 25.39, 21.73, 19.94])
        tau1 = np.array([35.98, 59.57, 94.41, 142.2, 164.9])*1e-12
        b1 = np.array([0.94, 0.934, 0.956, 0.988, 1])
        de2 = np.array([0.86, 1.79, 1.44, 0.77, 0.74])
        tau2 = np.array([7.44, 9.5, 8.91, 9.57, 10.4])*1e-12
        de3 = np.array([1.49, 1.21, 1.28, 1.32, 1.19])
        tau3 = np.array([1.25, 1.08, 1.39, 1.97, 1.69])*1e-12
    else:
        print('out of range. See the remarks.')
        return None
    
    f1 = interp1d(Concentration,einf)
    f2 = interp1d(Concentration,de1)
    f3 = interp1d(Concentration,tau1)        
    f4 = interp1d(Concentration,b1)
    f5 = interp1d(Concentration,de2)        
    f6 = interp1d(Concentration,tau2)        
    f7 = interp1d(Concentration,de3)
    f8 = interp1d(Concentration,tau3)
        
    einf,de1,tau1,b1,de2,tau2,de3,tau3 = f1(c), f2(c), f3(c), f4(c), f5(c), f6(c), f7(c), f8(c)
    
    omega = 2*np.pi*frequency
    data['epsilon'] = einf + de1/(1+(1j*omega*tau1)**b1) + de2/(1+1j*omega*tau2) + de3/(1+1j*omega*tau3)
    
    return data

def Propanol2Aqueous_Sato(frequency,temperature=25,c=1.0):
    data = {}
    data['minFREQ'] = 0.1e9
    data['maxFREQ'] = 89e9
    data['Citation'] = 'Sato, T., & Buchner, R. (2003). Dielectric relaxation spectroscopy of 2-propanol–water mixtures. The Journal of chemical physics, 118(10), 4606-4613.'
    data['minTemperature (C)'] = 25
    data['maxTemperature (C)'] = 25    
    data['minConcentration (x)'] = 0
    data['maxConcentration (x)'] = 1
    data['Remarks'] = 'c in mole fraction of the substance in water. The mole fraction between 0.14 and 0.3 should be avoided.'
    
    if c >=0 and c<=0.14:
        Concentration = np.array([0, 0.03, 0.065, 0.14])
        einf = np.array([3.96, 3.8, 3.7, 3.4])
        de1 = np.array([72.15, 67.05, 60.95, 48.04])
        tau1 = np.array([8.32, 12.41, 17.89, 27.26])*1e-12
        b1 = np.array([1, 0.983, 0.955, 0.939])
        de2 = np.array([0,0,0,0])
        tau2 = np.array([0,0,0,0])*1e-12
        de3 = np.array([2.14, 2.15, 1.91, 1.74])
        tau3 = np.array([0.39, 0.81, 0.73, 1.48])*1e-12
    elif c>=0.3 and c<=1:
        Concentration = np.array([0.3, 0.5, 0.7, 0.9, 1])
        einf = np.array([2.95, 2.8, 2.7, 2.54, 2.48])
        de1 = np.array([31.43, 20.63, 16.64, 15.53, 15.68])
        tau1 = np.array([49.7, 90.57, 149.4, 254.5, 354.6])*1e-12
        b1 = np.array([0.912, 0.962, 0.981, 0.991, 1])
        de2 = np.array([1.47, 2.02, 1.3, 0.78, 0.55])
        tau2 = np.array([9.63, 14.1, 17.7, 20.6, 23.4])*1e-12
        de3 = np.array([1.17, 0.97, 0.86, 0.72, 0.63])
        tau3 = np.array([0.89, 1.86, 2.57, 2.2, 2.12])*1e-12
    else:
        print('out of range. See the remarks.')
        return None
    
    f1 = interp1d(Concentration,einf)
    f2 = interp1d(Concentration,de1)
    f3 = interp1d(Concentration,tau1)        
    f4 = interp1d(Concentration,b1)
    f5 = interp1d(Concentration,de2)        
    f6 = interp1d(Concentration,tau2)        
    f7 = interp1d(Concentration,de3)
    f8 = interp1d(Concentration,tau3)
        
    einf,de1,tau1,b1,de2,tau2,de3,tau3 = f1(c), f2(c), f3(c), f4(c), f5(c), f6(c), f7(c), f8(c)
    
    omega = 2*np.pi*frequency
    data['epsilon'] = einf + de1/(1+(1j*omega*tau1)**b1) + de2/(1+1j*omega*tau2) + de3/(1+1j*omega*tau3)
    
    return data

def Propanol1Aqueous_Sato(frequency,temperature=25,c=1.0):
    data = {}
    data['minFREQ'] = 0.1e9
    data['maxFREQ'] = 25e9
    data['Citation'] = 'Sato, T., Chiba, A., & Nozaki, R. (2000). Composition-dependent dynamical structures of 1-propanol–water mixtures determined by dynamical dielectric properties. The Journal of Chemical Physics, 113(21), 9748-9758.'
    data['minTemperature (C)'] = 20
    data['maxTemperature (C)'] = 30    
    data['minConcentration (x)'] = 0
    data['maxConcentration (x)'] = 1
    data['Remarks'] = 'c in mole fraction of the substance in water. The mole fraction in the range of 0.30 - 0.35 should be avoided, since the fitting model is changed. For 30 C data, UnivariateSpline function was used to estimate the parameters that are not given in the original paper.'
    
    omega = 2*np.pi*frequency
    
    Temperature = np.array([20,25,30])
    if c >=0 and c<=0.30:
        Concentration = np.array([0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.22, 0.25, 0.3])
        Concentration = np.tile(Concentration,3)
        Temperature = np.repeat(Temperature,int(len(Concentration)/3))
        einf = np.array([5.3, 5.3, 5.3, 5.2, 5.2, 5.2, 5.1, 5.1, 5, 4.9, 4.9, 4.9, 4.8, 4.7, 4.7, 4.7, 4.7, 4.7, 4.6, 4.6, 4.6, 4.6, 4.6, 4.6, 4.6, 4.4, 4.5, 4.4, 4.4, 4.3, 4.2, 4.2, 4.2, 3.9, 5.2, 5.2, 5.2, 5.1, 5.1, 5, 5, 4.9, 4.8, 4.8, 4.8, 4.7, 4.7, 4.6, 4.6, 4.6, 4.6, 4.6, 4.6, 4.5, 4.5, 4.5, 4.4, 4.4, 4.4, 4.4, 4.4, 4.4, 4.3, 4.2, 4.2, 4.1, 4, 3.9, 5.1, 5.1, 5.1, 5.1, 5, 4.9, 4.9, 4.9, 4.8, 4.7, 4.7, 4.7, 4.7, 4.7, 4.6, 4.6, 4.6, 4.6, 4.6, 4.5, 4.5, 4.5, 4.5, 4.4, 4.4, 4.4, 4.3, 4.3, 4.3, 4.3, 4.2, 4.1, 4, 3.9])
        de1 = np.array([75.1, 74.2, 73.36, 72.43, 71.55, 70.83, 70.1, 69.34, 68.55, 67.3, 66.87, 65.94, 64.92, 63.73, 62.9, 61.93, 61.12, 60.27, 59.55, 58.82, 57.71, 56.18, 54.71, 52.94, 51.74, 49.51, 49.51, 47.86, 46.71, 45.77, 44.51, 42.16, 39.27, 35.29, 73.12, 72.19, 71.44, 70.41, 69.58, 68.89, 67.94, 67.26, 66.31, 65.49, 64.79, 63.93, 62.94, 61.92, 61.1, 60.05, 59.09, 58.22, 57.43, 56.67, 55.79, 54.02, 52.56, 50.97, 49.56, 48.45, 47.15, 45.97, 44.5, 43.55, 42.74, 40.64, 37.41, 34.1, 71.1, 70.5, 69.58, 68.54, 67.58, 66.84, 66.12, 65.24, 64.35, 63.6, 62.86, 62, 61.07, 60.13, 59.17, 58.17, 57.26, 56.54, 55.82, 54.92, 54.02, 52.64, 50.91, 49.2, 47.97, 47.03, 45.43, 44.08, 43.05, 42.07, 41.09, 39.13, 36.35, 32.85])
        log10tau1 = np.array([-11.032, -11.009, -10.983, -10.956, -10.93, -10.904, -10.883, -10.863, -10.844, -10.823, -10.803, -10.783, -10.765, -10.748, -10.733, -10.718, -10.704, -10.69, -10.678, -10.666, -10.654, -10.632, -10.612, -10.593, -10.575, -10.558, -10.54, -10.523, -10.506, -10.489, -10.473, -10.44, -10.392, -10.313, -11.083, -11.061, -11.036, -11.011, -10.987, -10.964, -10.945, -10.929, -10.911, -10.892, -10.873, -10.855, -10.839, -10.824, -10.81, -10.797, -10.783, -10.77, -10.757, -10.745, -10.733, -10.71, -10.688, -10.668, -10.649, -10.632, -10.615, -10.598, -10.581, -10.564, -10.548, -10.515, -10.468, -10.39, -11.134, -11.112, -11.089, -11.066, -11.044, -11.025, -11.008, -10.994, -10.979, -10.961, -10.942, -10.925, -10.911, -10.899, -10.887, -10.875, -10.863, -10.85, -10.837, -10.824, -10.812, -10.787, -10.764, -10.742, -10.723, -10.705, -10.688, -10.671, -10.655, -10.638, -10.622, -10.59, -10.543, -10.467])
        b1 = np.array([1, 0.995, 0.989, 0.985, 0.981, 0.983, 0.972, 0.967, 0.965, 0.958, 0.955, 0.949, 0.945, 0.948, 0.947, 0.939, 0.941, 0.939, 0.934, 0.932, 0.929, 0.921, 0.919, 0.918, 0.913, 0.8, 0.902, 0.893, 0.895, 0.894, 0.893, 0.886, 0.876, 0.856, 1, 0.997, 0.993, 0.988, 0.985, 0.983, 0.983, 0.973, 0.971, 0.97, 0.961, 0.955, 0.951, 0.947, 0.947, 0.944, 0.942, 0.935, 0.936, 0.935, 0.935, 0.932, 0.927, 0.923, 0.915, 0.909, 0.908, 0.903, 0.897, 0.895, 0.894, 0.89, 0.884, 0.867, 1, 1, 0.997, 0.993, 0.991, 0.99, 0.986, 0.974, 0.962, 0.96, 0.961, 0.959, 0.955, 0.952, 0.949, 0.945, 0.943, 0.943, 0.942, 0.935, 0.928, 0.927, 0.925, 0.92, 0.908, 0.905, 0.897, 0.893, 0.895, 0.898, 0.9, 0.894, 0.872, 0.87])
        
        f1 = interp2d(Temperature,Concentration,einf)
        f2 = interp2d(Temperature,Concentration,de1)
        f3 = interp2d(Temperature,Concentration,log10tau1)
        f4 = interp2d(Temperature,Concentration,b1)
        
        einf = f1(temperature,c)
        de1 = f2(temperature,c)
        log10tau1 = f3(temperature,c)
        b1 = f4(temperature,c)
        tau1 = 10**log10tau1
        data['epsilon'] = einf + de1/(1+(1j*omega*tau1)**b1)
        return data
    
    elif c>=0.35 and c<=1:
        Concentration = np.array([0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1])
        Concentration = np.tile(Concentration,3)
        Temperature = np.repeat(Temperature,int(len(Concentration)/3))
        einf = np.array([3.8, 3.7, 3.6, 3.6, 3.4, 3.4, 3.4, 3.3, 3.3, 3.3, 3.7, 3.6, 3.6, 3.5, 3.4, 3.3, 3.3, 3.2, 3.2, 3.2, 3.7, 3.6, 3.6, 3.5, 3.4, 3.3, 3.3, 3.2, 3.2, 3.2])
        de1 = np.array([31.17, 27.8, 25.5, 23.78, 20.93, 19.25, 18.37, 17.59, 17.32, 17.01, 30.39, 27.27, 24.89, 23.07, 20.42, 18.95, 17.82, 17.12, 16.84, 16.6, 28.76, 26.6, 24.13, 22.33, 19.61, 18.06, 17.11, 16.54, 16.21, 16.04])
        log10tau1 = np.array([-10.234, -10.155, -10.074, -9.997, -9.866, -9.748, -9.632, -9.515, -9.457, -9.398, -10.312, -10.235, -10.155, -10.078, -9.948, -9.832, -9.717, -9.601, -9.544, -9.486, -10.392, -10.313, -10.234, -10.158, -10.03, -9.915, -9.801, -9.687, -9.63, -9.573])
        b1 = np.array([0.885, 0.911, 0.934, 0.932, 0.943, 0.965, 0.975, 0.987, 0.994, 1, 0.89, 0.916, 0.935, 0.939, 0.946, 0.967, 0.976, 0.986, 0.993, 1, 0.902, 0.91, 0.932, 0.941, 0.948, 0.968, 0.973, 0.988, 0.995, 1])
        de2 = np.array([0.76, 1.15, 1.33, 1.66, 1.85, 1.69, 1.33, 0.99, 0.87, 0.71, 0.76, 1.12, 1.32, 1.58, 1.83, 1.6, 1.29, 1.13, 0.81, 0.7, 0.75, 1.1, 1.3, 1.58, 1.7, 1.32, 1.35, 1.12, 0.83, 0.72])
        log10tau2 = np.array([-11.04, -11.02, -11.03, -10.99, -10.96, -10.93, -10.89, -10.855, -10.82, -10.805, -11.055, -11.04, -11.03, -11.01, -10.98, -10.95, -10.905, -10.87, -10.83, -10.82, -11.07, -11.055, -11.04, -11.02, -10.995, -10.97, -10.92, -10.885, -10.85, -10.835])
        
        f1 = interp2d(Temperature,Concentration,einf)
        f2 = interp2d(Temperature,Concentration,de1)
        f3 = interp2d(Temperature,Concentration,log10tau1)
        f4 = interp2d(Temperature,Concentration,b1)
        f5 = interp2d(Temperature,Concentration,de2)
        f6 = interp2d(Temperature,Concentration,log10tau2)
        
        einf = f1(temperature,c)
        de1 = f2(temperature,c)
        log10tau1 = f3(temperature,c)
        b1 = f4(temperature,c)
        de2 = f5(temperature,c)
        log10tau2 = f6(temperature,c)
        tau1 = 10**log10tau1
        tau2 = 10**log10tau2        
        data['epsilon'] = einf + de1/(1+(1j*omega*tau1)**b1) + de2/(1+(1j*omega*tau2))
        return data        
    else:
        print('out of range. See the remarks.')
        return None