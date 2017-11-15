#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 23:02:46 2017

@author: lyn
"""

 # CALL = 1 PUT = -1
import math  
import sys 
sys.path.append("..")
import py_lets_be_rational as lets_be_rational
 # CALL = 1 PUT = -1
 
        
def cal_implied_vol(price, F, K, r, T, opt_type):
    _undisc_price = price * math.exp(r * T)
    _opt_type = 1 if opt_type == 'C' else -1  
    try :
        return lets_be_rational.implied_volatility_from_a_transformed_rational_guess(_undisc_price, F, K, T, _opt_type)
    except:
        return math.nan   
if __name__ == "__main__":   
    print(cal_implied_vol(0.0005,2.5, 2.2, 0.04,0.0556, "P"))
    print(cal_implied_vol(0.0005,2.5, 2.2, 0.04,0.0556, "C"))

