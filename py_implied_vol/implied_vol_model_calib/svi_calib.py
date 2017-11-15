#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 20:22:33 2017

@author: ocs3
"""

import numpy as np
import pandas as pd
from scipy import optimize
import sys
sys.path.append("..")
import cal_implied_vol


INTEREST_RATE = 0.025

import sys
# parameter boundary check
def acceptable(sigma, a, d, c, vT):
    eps = sys.float_info.epsilon
#    print(eps)
    return  -eps < c and \
            c < 4 * sigma +eps and \
            abs(d)-eps < min(c, 4 * sigma - c) +eps and \
            -eps < a and \
            a < min(vT.max() ,10) + eps
            
#object function
def sum_of_squares(x, a, d, c, vT,vega):
    diff = (a + d * x + c * np.sqrt(x * x + 1) - vT) 
    return (vega * diff * diff).mean()


# explicit-Qausi method by Zeiled
#step 1: slove the gradiant(the parameterizarion meet the free-arbitrage condition)
def solve_grad( S, M, x, vT,vega):
    # var @ S, M:paramters,
    #     @ x: Moneyness
    #     @ vT: maturity
    #     @ vega: option vega, performing as the LSM weights
    ys = (x - M) / S 
    w = vega.mean()
    y = (vega * ys ).mean()
    y2 = (vega *ys * ys).mean()
    y2one =(vega * (ys * ys + 1)).mean()
    ysqrt =( vega *np.sqrt(ys * ys + 1)).mean()
    y2sqrt = (vega *(ys * np.sqrt(ys * ys + 1))).mean()
    v = (vega *vT).mean()
    vy = (vega *vT * ys).mean()
    vsqrt = (vega *vT * np.sqrt(ys * ys + 1)).mean()

    matrix = [
        [w, y, ysqrt],
        [y, y2, y2sqrt],
        [ysqrt, y2sqrt, y2one]
    ]
    
    vector = [v, vy, vsqrt]
    # opmtimization get from  Largranie
    _a, _d, _c = np.linalg.solve(np.array(matrix), np.array(vector))
    # check if parameters in free-arbitrage constraints
    if acceptable(S, _a, _d, _c, vT):
        #print(1)
#        print(acceptable(S, _a, _d, _c, vT))
        _cost = sum_of_squares(ys, _a, _d, _c, vT, vega)
        #print(_cost)
        #return if it's ture
        return _a, _d, _c, _cost
    
    # optimization on the contraints boundary.
    _a, _d, _c, _cost = None, None, None, None
    for matrix, vector, clamp_params in [
        ([[1,0,0],[y,y2,y2sqrt],[ysqrt,y2sqrt,y2one]], [0, vy, vsqrt], False), # a = 0
        ([[1,0,0],[y,y2,y2sqrt],[ysqrt,y2sqrt,y2one]], [vT.max(),vy,vsqrt], False), # a = _vT.max()
        ([[1,y,ysqrt],[0,-1,1],[ysqrt,y2sqrt,y2one]], [v,0,vsqrt], False), # d = c
        ([[1,y,ysqrt],[0,1,1],[ysqrt,y2sqrt,y2one]], [v,0,vsqrt], False), # d = -c
        ([[1,y,ysqrt],[0,1,1],[ysqrt,y2sqrt,y2one]], [v,4*S,vsqrt], False), # d <= 4*s-c
        ([[1,y,ysqrt],[0,-1,1],[ysqrt,y2sqrt,y2one]], [v,4*S,vsqrt], False), # -d <= 4*s-c
        ([[1,y,ysqrt],[y,y2,y2sqrt],[0,0,1]], [v,vy,0], False), # c = 0
        ([[1,y,ysqrt],[y,y2,y2sqrt],[0,0,1]], [v,vy,4*S], False), # c = 4*S

        ([[1,y,ysqrt],[0,1,0],[0,0,1]], [v,0,0], True), # c = 0, implies d = 0, find optimal a
        ([[1,y,ysqrt],[0,1,0],[0,0,1]], [v,0,4*S ], True), # c = 4s, implied d = 0, find optimal a
        ([[1,0,0],[0,-1,1],[ysqrt,y2sqrt,y2one]], [0,0,vsqrt], True), # a = 0, d = c, find optimal c
        ([[1,0,0],[0,1,1],[ysqrt,y2sqrt,y2one]], [0,0,vsqrt], True), # a = 0, d = -c, find optimal c
        ([[1,0,0],[0,1,1],[ysqrt,y2sqrt,y2one]], [0,4*S,vsqrt], True), # a = 0, d = 4s-c, find optimal c
        ([[1,0,0],[0,-1,1],[ysqrt,y2sqrt,y2one]], [0,4*S,vsqrt], True) # a = 0, d = c-4s, find optimal c
    ]:
        a, d, c = np.linalg.solve(np.array(matrix), np.array(vector))
#        _a, _d, _c, cost  = a, d, c, sum_of_squares(ys, a, d, c, vT)
#        print(acceptable(S, _a, _d, _c, vT))
       
        if clamp_params:
            dmax = min(c, 4 * S - c)
            a = min(max(a, 0), vT.max())
            d = min(max(d, -dmax), dmax)
            c = min(max(c, 0), 4 * S)
        
        cost = sum_of_squares(ys, a, d, c, vT, vega)
        if acceptable(S, a, d, c, vT) and (_cost is None or cost < _cost):
            _a, _d, _c, _cost = a, d, c, cost
            #print(2)
            #print(cost)
    a_margin = [0, min(vT.max() ,10) ]
    c_margin = [0, 4 * S]
    for _a_ in  a_margin:
        for _c_ in c_margin:
            cost = sum_of_squares(ys, _a_, d, _c_, vT, vega)
            if _cost is None or cost < _cost:
                _a, _d, _c, _cost = a, d, c, cost
                #print(3)
        
    assert _cost is not None, "_cost is None, S=%s, M=%s" % (S, M)
    #print(_a, _d, _c)
    return _a, _d, _c, _cost



def solve_grad_get_score(p,x, vT,vega):
    S, M = p
    return solve_grad(S, M, x, vT,vega)[3]

#initial guess
def inital_guess(df):
    vT = df.MktImpliedVol * df.MktImpliedVol * df.Maturity
    abs_moneyness = abs(df.Moneyness)
    n = len(df.Maturity)
    if n <3:
        mInitial = 0.1
        sigmaInitial = 0.0
    else:
              
        aL = (df.Moneyness.iloc[0]*vT.iloc[1] - vT.iloc[0]*df.Moneyness.iloc[1])/(df.Moneyness.iloc[0]-df.Moneyness.iloc[1])
        aR = (df.Moneyness.iloc[n-1]*vT.iloc[n-2] - vT.iloc[n-1]*df.Moneyness.iloc[n-2])/(df.Moneyness.iloc[n-1]-df.Moneyness.iloc[n-2])
        
        bL = -(vT.iloc[0] - vT.iloc[1]) / (df.Moneyness.iloc[0] - df.Moneyness.iloc[1])
        bR = (vT.iloc[n-1] - vT.iloc[n-2]) / (df.Moneyness.iloc[n-1] - df.Moneyness.iloc[n-2])
        
        recurring = 0 if bL == bR else bL*(aR-aL)/(bL-bR)
        aInitial = aL + recurring 
        rhoInitial = 0 if bL == bR else (bL+bR)/(bR-bL)
        bInitial = 0.5*(bR-bL)
        vol_min = df.MktImpliedVol.min()
        mInitial = df[df.MktImpliedVol == vol_min]['Moneyness'].iloc[0]  
#        mInitial = recurring/(bInitial * (rhoInitial-1.0))
        sigmaInitial_ = abs((-(vT.min()) + + aL + recurring) / (bInitial * np.sqrt( abs (1.0 - rhoInitial * rhoInitial))) )
        sigmaInitial = min(max(sigmaInitial_, 0.0),10)
    return sigmaInitial, mInitial
#print(inital_guess(item))              
        
    # perform calibration
def calibrate(df):
    sigmaInitial, mInitial =  inital_guess(df)
    vT = df.MktImpliedVol * df.MktImpliedVol * df.Maturity
    res = optimize.minimize(solve_grad_get_score, [sigmaInitial, mInitial], args=(df.Moneyness, vT,df.vega), bounds=[(0.001, None), (None, None)])
    assert res.success
    S, M = res.x
    a, d, c, _ = solve_grad(S, M, df.Moneyness, vT,df.vega)
    #print(a, d, c)
    T = df.Maturity.max() # should be the same for all rows
    if c != 0 :
        A, P, B = a / T, d / c, c /(S*T)
    else:
        A, P, B = a / T, 0, 0
#    print(A, P, B)
    assert T >= 0 and S >= 0 and abs(P) <= 1

    return A, P, B, S, M


if __name__ == '__main__':
    import numpy as np
    from data_process import *
    import implied_vol_model.svi as svi
    cal_implied_vol
    data_read()
    TIME_TEST = pd.datetime(2017,6,14,10,30,0,0)
    data_process(TIME_TEST)
    mkt_Option_Equity,option_set_by_maturity  = data_process_kind('ask', TIME_TEST)
    item = option_set_by_maturity[0][['Maturity', 'F', 'Moneyness','MktImpliedVol']]
    item['vega'] = 1
    A, P, B, S, M = calibrate(item) 
    svi_impl_vol = item.apply( lambda x: svi.SVI(A, P, B, S, M).SVI_vol(x.Moneyness), axis = 1)

