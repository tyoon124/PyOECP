# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 20:19:11 2021

@author: tae-jun_yoon
"""

from setuptools import setup

setup(
      name = 'PyOECP',
      version = '0.0.1',
      author = 'Tae Jun Yoon',
      author_email = 'oheim124@gmail.com',
      packages=['PyOECP'],
      py_modules=['PyOECP/Transform','PyOECP/References','PyOECP/Models'],
      url='To be declared',
      description='A flexible open-source software library for estimating and modeling the complex permittivity based on the open-ended coaxial probe technique',
)
