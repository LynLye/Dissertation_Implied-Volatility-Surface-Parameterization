#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 16:55:38 2017

@author: lyn
"""

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup



setup(
     name='py_implied_vol',
     version='0.1.1',
     description='Python implied vol parameterization based on SVI and SABR',
     author='Ye,Lin',
     author_email='yelin_lyn@163.com',
     url='https://github.com/LynLye/Dissertation_Implied-Volatility-Surface-Parameterization',
    package_dir={'Dissertation_Implied-Volatility-Surface-Parameterization':'py_implied_vol'},
    include_package_data=True,
     )
