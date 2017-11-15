#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 06:08:24 2017

@author: lyn
"""

import numpy as np
import scipy
import math
import sys
sys.path.append("..")
import cal_implied_vol
import implied_vol_model.sabr as sabr
_B_ = 0



def ATMVolStats(z_d1, sigma_d1, z_0, sigma_0, z_u1, sigma_u1):
    w_d1 = 1/(z_d1 - z_0)/(z_d1 - z_u1)
    w_0 = 1/(z_0 - z_d1)/(z_0 - z_u1)
    w_u1 = 1/(z_u1 - z_0)/(z_u1 - z_d1)
    sigma = z_0 * z_u1 * w_d1 * sigma_d1 + z_d1 * z_u1 * w_0 * sigma_0 \
            + z_d1 * z_0 * w_u1 * sigma_u1
    sigmaSkew = -(z_0 + z_u1) * w_d1 * sigma_d1 -  (z_d1 + z_u1) * w_0 * sigma_0 \
                 - (z_d1 + z_0) * w_u1 * sigma_u1
    sigmaCur =  2 * w_d1 * sigma_d1 + 2 * w_0 * sigma_0 + 2 * w_u1 * sigma_u1
    return sigma, sigmaSkew, sigmaCur
def sgn(x):
    return (-1 + 2/(1+math.exp(-(x+1)/2)))

def recalATMVolStats(alpha, beta, rho, vu, F, b):
    sigma = alpha * (F + b) ** (beta -1)
    sigmaSkew = (rho * vu - ( 1- beta) * sigma)/2
    sigmaCur = vu **2 / 3/ sigma +\
                ((1 - beta)**2 * sigma**2 - 3 * rho**2 * vu**2)/6/sigma
    return sigma, sigmaSkew, sigmaCur

           
def is_acceptable(a, b, c, F, b_, beta):
#    alpha = sigma *　(F + b)**( 1- beta)
    nuSquare = 3 * a * c - 1 / 2 * ( 1- beta) **2 *a**2 + 3 * ( 2 * b + (1 - beta) *a)**2 / 2
    if  nuSquare < 0:
        rho = sgn(( 2 * b + ( 1 - beta) * a))
        nu = ( 2 * b + ( 1 - beta) * a) / rho
    else:
        nu = math.sqrt(nuSquare)
        rho = ( 2 * b + ( 1 - beta) * a)/nu
    alpha = a * (F + b_) **(1 - beta)
    return  alpha, rho, nu


def get_terms(alpha, rho, nu,F, b_, beta):
    a = alpha* pow((F+b_),(beta-1))
    b = (rho * nu - (1-beta)*a)/2
    c = nu**2/(3 * a) + (pow((1-beta),2)*pow(a,2) - 3*pow(rho,2)*pow(nu,2)/(6*a))
    return a,b,c

# waiting for futher check
def get_vega(S,K,T,r,sigma):
    d1 = (math.log(S/K) +(r+(sigma**2)/2)*(T))/(sigma*math.sqrt(T))
    n1 = math.exp(-d1*d1/2)/(2 * math.pi)
    vega = S*math.sqrt(T)*n1
    return vega


def sum_of_sqaures(alpha, rho, nu, weight, moneyness,mkt_imp,F, b_, beta,T):
#    _a, _b, _c  = get_terms(alpha, rho, nu,F, b_, beta)
    diff = sabr.SABR(alpha, beta, rho,nu).SABR_vol_lognormal_by_moneyness(moneyness, T, F) - mkt_imp
    return (weight*diff * diff).mean()



# guess rho and nu 
# output: alpha, rho, nu
def guessParam(sigma, sigmaSkew, sigmaCur, F, b, beta):
#    alpha = sigma *　(F + b)**( 1- beta)
    nuSquare = 3 * sigma * sigmaCur - 1 / 2 * ( 1- beta) **2 * sigma**2 + 3 * ( 2 * sigmaSkew + (1 - beta) * sigma)**2 / 2
    if  nuSquare < 0:
        rho = sgn(( 2 * sigmaSkew + ( 1 - beta) * sigma))
        nu = ( 2 * sigmaSkew + ( 1 - beta) * sigma) / rho
    else:
        nu = math.sqrt(nuSquare)
        rho = ( 2 * sigmaSkew + ( 1 - beta) * sigma)/nu
    alpha = sigma * (F + b) **(1 - beta)
    return  alpha, rho, nu


def findAlpha(sigma, rho, nu, beta, T, F, b):
    def findAlphaFuncSub(sigma, rho, nu, beta, T, F, b):
        func1 = T * (1 - beta)**2 / 24 / ((F + b) ** (2 - 2 * beta))
        func2 = rho * beta * nu * T / ((F + b) ** (1 - beta)) / 4
        func3 = 1 + (2 - 3 * rho **2)* nu**2 * T /24
        func4 = sigma * (F + b)** (1 - beta)
        return func1, func2, func3, func4

    func1, func2, func3, func4 = findAlphaFuncSub(sigma, rho, nu, beta, T, F, b)
    
    def findAlpha(alpha):
        return func1* alpha ** 3 + func2 * alpha **2 + func3 * alpha - func4
    
    
    def difAlpha(alpha):
        return 3 * func1* alpha ** 2 + 2 *func2 * alpha + func3 
    
    def sloveAlpha():
        roots = scipy.optimize.fsolve(findAlpha, 0, fprime = difAlpha) 
        return roots[0]
#        if roots > 0:
#            return roots
#        else:
#            return 0.2

    return sloveAlpha()
 
    
def get_initial_guess(sigma, sigmaSkew, sigmaCur,T ,F, b_, beta):
    _alpha, _rho, _nu = guessParam(sigma, sigmaSkew, sigmaCur, F, b_, beta)
    _alpha = findAlpha(sigma, _rho, _nu, beta, T, F, b_)
    return _alpha, _rho, _nu

def get_cost_scores(parm,  weight, moneyness,mkt_imp,F, b_, beta,T):
    alpha, rho, nu = parm
    return sum_of_sqaures(alpha, rho, nu,  weight, moneyness,mkt_imp,F, b_, beta,T)
    
    

    
from scipy import optimize


def calibrate(df, alpha_init, rho_init,nu_init):
    res = optimize.minimize(get_cost_scores, [alpha_init, rho_init,nu_init], args=(df.vega, df.Moneyness, df.MktImpliedVol,df.F, _B_, 0.5, df.Maturity), bounds=[(0,None),(-1,1),(0, None)])
    assert res.success
    alpha, rho, nu = res.x
    #print(alpha, rho, nu) 
    return alpha, rho, nu







if __name__ == '__main__':
    import numpy as np
    from data_process import *
    data_read()
    TIME_TEST = pd.datetime(2017,6,14,10,30,0,0)
    data_process(TIME_TEST)
    mkt_Option_Equity,option_set_by_maturity  = data_process_kind('ask', TIME_TEST)
    item = option_set_by_maturity[0][['Maturity', 'F', 'Moneyness','MktImpliedVol']]
    item['SelectedOption'] = item.apply(lambda x: 'P' if x.Moneyness < 0 else 'C', axis=1)
    Num_Put = item.SelectedOption.value_counts()[0]
    item.loc[item.index,'No'] = range(len(item.index))- Num_Put
    a,b,c = ATMVolStats(item[item.No == 1]['Moneyness'].iloc[0], item[item.No ==1]['MktImpliedVol'].iloc[0], \
    item[item.No == 0]['Moneyness'].iloc[0], item[item.No == 0]['MktImpliedVol'].iloc[0], \
    item[item.No == -1]['Moneyness'].iloc[0], item[item.No ==-1]['MktImpliedVol'].iloc[0])
    item['vega'] = 1
    alpha_init, rho_init,vu_init = get_initial_guess(a,b,c,item[item.No == 0]['Maturity'].iloc[0],item[item.No == 0]['F'].iloc[0],0.5, _B_)
    alpha, rho, nu = calibrate(item, alpha_init, rho_init,vu_init)    
    sabr_implied_vol = sabr.SABR(alpha, 0.5, rho, nu).SABR_vol_lognormal_by_moneyness(item.Moneyness, item.Maturity, item.F)


