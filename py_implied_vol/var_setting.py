#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 23:04:22 2017

@author: lyn
"""



def enum(**enums):
    return type('Enum', (), enums)

option_type = enum(CALL='C', PUT='P')
ParamModel = enum(SABR_MODEL='sabr', SVI_PARAM='svi')
greeks_measure = enum(VALUE='value', DELTA='delta', THETA='theta', RHO='rho', VEGA='vega', GAMMA='gamma')
price_type = enum(ASK='ask', BID='bid', MID='bid')

DEFAULT_RISK_FREE_RATE = 0.025
DEFAULT_TRANS_DAYS = 252

