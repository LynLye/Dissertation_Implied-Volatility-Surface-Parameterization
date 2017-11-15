#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 08:00:48 2017

@author: lyn
"""
import numpy as np
class SVI():
    def __init__(self, a, p, b, s, m):
        self.a = a
        self.p = p
        self.b = b
        self.s = s
        self.m = m
        
    def SVI_vol(self, x):
        return np.sqrt(self.a + self.b * (self.p * (x - self.m) + np.sqrt((x - self.m) * (x - self.m) + self.s * self.s)))