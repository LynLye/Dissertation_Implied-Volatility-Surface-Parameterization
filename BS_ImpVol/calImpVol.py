# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 01:30:31 2017

@author: Lyn Ye
"""
import py_lets_be_rational as lets_be_rational
 # CALL = 1 PUT = -1
import math
        
def calImpVol(price, F, K, r, T, opt_type):
    _undisc_price = price * math.exp(r * T)
    _opt_type = 1 if opt_type == 'C' else -1  
    try :
        return lets_be_rational.implied_volatility_from_a_transformed_rational_guess(_undisc_price, F, K, T, _opt_type)
    except:
        return math.nan