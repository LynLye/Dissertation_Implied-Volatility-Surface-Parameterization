#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 05:37:13 2017

@author: lyn
"""
import abc
import pandas as pd 
import implied_vol_model.sabr
import implied_vol_model_calib.sabr_calib
import implied_vol_model.svi
import implied_vol_model_calib.svi_calib


class ParamModel():
    SABR='sabr'
    SVI='svi'
ParamModels = ParamModel()







DEFAULT_RISK_FREE_RATE = 0.025
DEFAULT_TRANS_DAYS = 252
_b_ = 0


class ImplVolSurfCndt(object):
    @abc.abstractmethod
    def run_model(self):
        """perform implied volatility surface parameterization"""

class Parameterization(ImplVolSurfCndt):
    def __init__(self, mkt_impl_vol, F, moneyness,maturity, No=None, SelectedOption=None):
        self.mkt_impl_vol = mkt_impl_vol
        self.F = F
        self.monyness = moneyness
        self.maturity = maturity
        self.df = pd.DataFrame({'MktImpliedVol': mkt_impl_vol, 'Moneyness': moneyness})
        self.df['F'] = F
        self.df['Maturity'] = maturity
        self.model_cache = {}
        self.model_cache_param_hashes = {}
        self.No = No
        self.SelectedOption = SelectedOption
    
    def copy(self):
        return Parameterization(self.mkt_impl_vol, self.F, self.monyness, self.No, self.SelectedOption)
    
    def param_hash(self):
        return hash((self.F, self.maturity))

    def run_model(self, model=ParamModels.SABR, **kwargs):
        curr_param_hash = self.param_hash()
        if model in self.model_cache:
            prev_param_hash = self.model_cache_param_hashes[model]
            if self.param_hash() == prev_param_hash:
                return self.model_cache[model]
        
        if model == ParamModels.SABR:
            #print(ParamModels.SABR)
            self.model_cache[model] = self.SABR_param_calib(**kwargs)
        elif model == ParamModels.SVI:
            #print(ParamModels.SVI)
            
            self.model_cache[model] = self.SVI_param_calib(**kwargs)
          
        result = self.model_cache[model]
        self.model_cache_param_hashes[model] = curr_param_hash
        return result
    
    
    def SABR_param_calib(self, set_beta = 0.5):
        item = self.df
        check_list = list(item.columns)
        if self.No is not None :
            pass
        else:
            if self.SelectedOption is not None :
                pass
            else:
                item['SelectedOption'] = item.apply(lambda x: 'P' if x.Moneyness < 0 else 'C', axis=1)
            Num_Put = item.SelectedOption.value_counts()[0]
            item.loc[item.index,'No'] = range(len(item.index))- Num_Put  
        a,b,c = sabr_calib.ATMVolStats(item[item.No == 1]['Moneyness'].iloc[0], item[item.No ==1]['MktImpliedVol'].iloc[0], \
        item[item.No == 0]['Moneyness'].iloc[0], item[item.No == 0]['MktImpliedVol'].iloc[0], \
        item[item.No == -1]['Moneyness'].iloc[0], item[item.No ==-1]['MktImpliedVol'].iloc[0])
        item['vega'] = 1
        alpha_init, rho_init,nu_init = sabr_calib.get_initial_guess(a,b,c,item[item.No == 0]['Maturity'].iloc[0],item[item.No == 0]['F'].iloc[0],0.5, 0)
        alpha, rho, nu = sabr_calib.calibrate(item, alpha_init, rho_init,nu_init)
        beta = set_beta
        
        sabr_implied_vol = item.apply(lambda x: \
            SABR_model.SABR(alpha, beta, rho, nu).SABR_vol_lognormal_by_moneyness(x.Moneyness, x.Maturity, x.F),axis=1) 
        self.df['SABR_impl_vol'] = sabr_implied_vol
        return {'alpha':alpha,
                'beta': beta,
                'rho': rho,
                'nu': nu}
       
    def SVI_param_calib(self):
        item = self.df
        item['vega'] = 1
        a,p,b,s,m= svi_calib.calibrate(item)
        svi_implied_vol = svi_param.SVI(a,p,b,s,m).SVI_vol(item.Moneyness) 
        self.df['SVI_impl_vol'] = svi_implied_vol
        return {'a': a,
                'p': p,
                'b': b,
                's': s,
                'm': m}

        
if __name__ == '__main__':
    import numpy as np
    from data_process import *
    data_read()
    TIME_TEST = pd.datetime(2017,6,14,10,30,0,0)
    data_process(TIME_TEST)
    mkt_Option_Equity,option_set_by_maturity  = data_process_kind('ask', TIME_TEST)
    item1 = option_set_by_maturity[0][['Maturity', 'F', 'Moneyness','MktImpliedVol']]
    item1['vega'] = 1
    params = Parameterization(item1.MktImpliedVol,item1.F.iloc[0], item1.Moneyness, item1.Maturity.iloc[0])

    for i in ['svi', 'sabr', '2']:
        try:
            params.run_model(i)
            print(params.model_cache)
        except:
            print('input an unsupported parameterization model, and please input a new one : "sabr" or "svi" ')
