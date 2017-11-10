#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 19:20:54 2017

@author: lyn
"""

def svi(A, P, B, S, M, x):
    return np.sqrt(A + B * (P * (x - M) + np.sqrt((x - M) * (x - M) + S * S)))
