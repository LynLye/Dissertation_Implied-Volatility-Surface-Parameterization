#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 03:13:01 2017

@author: lyn
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 00:12:11 2017

@author: lyn
"""
import pandas as pd
import numpy as np
import math
import sys 
sys.path.append("..")
#from cal_implied_vol import cal_implied_vol, calImpVol
from var_setting import *

import math
        
import py_lets_be_rational as lets_be_rational
 # CALL = 1 PUT = -1
import math
        
def cal_implied_vol(price, F, K, r, T, opt_type):
    _undisc_price = price * math.exp(r * T)
    _opt_type = 1 if opt_type == 'C' else -1  
    try :
        return lets_be_rational.implied_volatility_from_a_transformed_rational_guess(_undisc_price, F, K, T, _opt_type)
    except:
        return math.nan
def get_vega(S,K,T,r,sigma):
    d1 = (math.log(S/K) +(r+(sigma**2)/2)*(T))/(sigma*math.sqrt(T))
    n1 = math.exp(-d1*d1/2)/(2 * math.pi)
    vega = S*math.sqrt(T)*n1
    return vega

    
    
def data_read():
    global raw_mkt_data, instruments
    instruments = pd.read_csv('/Users/lyn/Desktop/untitled folder 2/untitled folder 2/Disseration0831/SABR/md1706/instruments.csv', sep = ',',  \
                             usecols = ['Type', 'Symbol', 'UnderlyingSymbol','Expiry',\
                             'ContractMultplierNominator','ContractMultplierDenominator','Strike','OptionType'],\
                            dtype = {'UnderlyingSymbol':str})
    
    raw_mkt_data = pd.read_csv('/Users/lyn/Desktop/untitled folder 2/untitled folder 2/Disseration0831/SABR/md1706/marketdata-20170614-qsar1.csv', sep = ',',\
                          header = None, usecols = [ 0, 3, 6, 8], names = ['Time', 'Symbol', 'Ask1', 'Bid1'],\
                        dtype = {'Symbol':str})
    
    
    instruments['Expiry'] = pd.to_datetime(instruments['Expiry'], format = '%Y%m%d')
    instruments = instruments[instruments['ContractMultplierNominator'] != 10220]
    raw_mkt_data['Bid1'] = raw_mkt_data[raw_mkt_data.Bid1 != 0.0].Bid1.astype(np.float64)
    raw_mkt_data['Ask1'] = raw_mkt_data[raw_mkt_data.Ask1 != 0.0].Ask1.astype(np.float64)
    try:
        raw_mkt_data['Time'] = pd.to_datetime(raw_mkt_data['Time'], format = '%Y-%b-%d %H:%M:%S.%f')
    except:
        raw_mkt_data['Time'] = pd.to_date5time(raw_mkt_data['Time'], format = '%Y-%b-%d %H:%M:%S')

#data_read()

def data_process(test_time):
    global Call, Put, mkt_Equity,use_data
    Equity = instruments[(instruments.Type == 'Equity')]
    Option = instruments[(instruments.Type == 'Option')&(instruments.Expiry >= test_time)]
    Option = Option.sort_values(by = 'Expiry', ascending = False)
    Option = Option.reset_index(drop=True)
    Call = Option[Option.OptionType == 'C']
    Put = Option[Option.OptionType == 'P']
    
    use_data = raw_mkt_data[raw_mkt_data.Time <= test_time].groupby('Symbol').tail(1)  
    
    mkt_Equity = pd.merge(use_data, Equity, 'right', on = 'Symbol')
    mkt_Equity['EquityPrice'] = (mkt_Equity.Ask1 + mkt_Equity.Bid1) / 2
    mkt_Equity = mkt_Equity[['EquityPrice','Symbol']]
    mkt_Equity.rename(columns={'Symbol':'UnderlyingSymbol'}, inplace = True)

   
    

def data_process_kind(_price_kind,test_date):
    if _price_kind == price_type.ASK:
        use_data['cal_price'] = use_data.Ask1 
    elif _price_kind == price_type.BID:
        use_data['cal_price'] = use_data.Bid1 
    elif _price_kind == price_type.MID:
        use_data['cal_price'] = (use_data.Ask1 + use_data.Bid1) / 2
    else:
         print('please input the right price kind!')
        
    
    mkt_Call = pd.merge(Call, use_data, 'left', on = 'Symbol')
    mkt_Call = mkt_Call[['Expiry','Strike','cal_price','UnderlyingSymbol']]
    mkt_Call.rename(columns={'cal_price':'call_cal_price'}, inplace = True)
#    return mkt_Call


    mkt_Put = pd.merge(Put, use_data, 'left', on = 'Symbol')
    mkt_Put = mkt_Put[['Expiry','Strike','cal_price','UnderlyingSymbol']]
    mkt_Put.rename(columns={'cal_price':'put_cal_price'}, inplace = True)
    
    mkt_Option_Pairs = pd.merge(mkt_Call, mkt_Put, how = 'inner', \
                                on = ['Strike','Expiry','UnderlyingSymbol'])
    mkt_Option_Pairs = mkt_Option_Pairs.dropna()
    

    #select Option to calculate
    mkt_Option_Pairs['OptionCalPrice'] = mkt_Option_Pairs.apply(lambda x: \
                    min(x.call_cal_price,x.put_cal_price), axis = 1)
    mkt_Option_Pairs['SelectedOption'] = mkt_Option_Pairs.apply(lambda x: \
                    'C' if x.call_cal_price < x.put_cal_price else'P', axis = 1)
    mkt_Option_Pairs['Maturity'] = ((mkt_Option_Pairs['Expiry']-test_date).astype('timedelta64[D]')+1)/ DEFAULT_TRANS_DAYS
                
    # merge selected option and Equity
    mkt_Option_Equity = pd.merge(mkt_Option_Pairs,mkt_Equity, 'left', on = 'UnderlyingSymbol')
    mkt_Option_Equity = mkt_Option_Equity[['Maturity','Strike','UnderlyingSymbol',\
                                           'EquityPrice','OptionCalPrice','SelectedOption','put_cal_price','call_cal_price']]
    mkt_Option_Equity = mkt_Option_Equity.sort_values(by =['Maturity', 'Strike'])
   
       
    mkt_Option_Equity_Group= mkt_Option_Equity.groupby('Maturity')
    mkt_Option_Equity_Set = []

    for key in mkt_Option_Equity_Group.groups.keys():
       mkt_Option_Equity_Set.append(mkt_Option_Equity_Group.get_group(key))
     
    Num_Set = len(mkt_Option_Equity_Set)
   
    for i in range(Num_Set):
        item = mkt_Option_Equity_Set[i]
        
        item = item.sort_values(by = ['Maturity','Strike'], \
                                                          ascending = True)
        Num_Put = item.SelectedOption.value_counts()[0]
        item.loc[item.index,'No'] = range(len(item.index))- Num_Put
        def calImpliedForward(C1, P1, K1, C2, P2, K2):
            forward = K1 + (K2 - K1)*(C1 - P1)/((C1 - P1) - (C2 - P2))
            return forward
           
        item ['F'] =  calImpliedForward(
                                                     item[item.No == 0]['call_cal_price'].iloc[0],\
                                                     item[item.No == 0]['put_cal_price'].iloc[0],\
                                                     item[item.No == 0]['Strike'].iloc[0],\
                                                     item[item.No == -1]['call_cal_price'].iloc[0],\
                                                     item[item.No == -1]['put_cal_price'].iloc[0],\
                                                     item[item.No == -1]['Strike'].iloc[0])
                
        # dealing with moneyness
        item['Moneyness'] = np.log(item['Strike'] / item['F'])
        item = item.sort_values(by = ['Maturity','Strike','Moneyness'], \
                                                      ascending = True)                
    # calculate market implied volatility
     # cal_implied_vol(price, F, K, r, T, opt_type)   
        item['MktImpliedVol'] = item.apply(lambda x: \
                          cal_implied_vol(x.OptionCalPrice, x.F, x.Strike,  0.025, x.Maturity,x.SelectedOption),axis=1)
        #S,K,T,r,sigma
        item['vega'] = item.apply(lambda x: \
                          get_vega(x.OptionCalPrice,  x.Strike, x.Maturity, 0.025, x.MktImpliedVol ),axis=1)
        mkt_Option_Equity_Set[i] = item
    return  mkt_Option_Equity, mkt_Option_Equity_Set

 

        
if __name__ == "__main__":
    data_read()
    TIME_TEST = pd.datetime(2017,6,14,10,30,0,0)
    data_process(TIME_TEST)
    mkt_Option_Equity,option_set_by_maturity  = data_process_kind('ask', TIME_TEST)
    
        