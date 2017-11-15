#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 23:44:13 2017

@author: lyn
"""


import numpy as np
class SABR():
    def __init__(self, alpha, beta, rho, nu):
        self.alpha = alpha
        self.beta = beta
        self.rho = rho
        self.nu = nu
        
    def SABR_vol_lognormal(self, K, T, F, b=0):
        if self.beta ==1:
            g_k = 0
            zeta_k = self.nu * np.log((F + b) / (K + b)) / self.alpha
        elif self.beta == 0:
            g_k = self.alpha * self.alpha / (24 * (F + b) * (K + b))
            zeta_k = self.nu * (F - K) / self.alpha
        else: 
            g_k = 1/ 24 * pow(self.beta-1, 2)* pow((F + b), beta -1) * pow((K + b), self.beta -1)* self.alpha * self.alpha
            zeta_k = self.nu * (pow(F + b,1-self.beta) - pow(K + b, 1 - self.beta)) / self.alpha / (1 - self.beta)        
        if F == K:
            sigmaLognormal = self.alpha * pow((F + b), self.beta -1) * (1 + (g_k + 1/4 * 
                                        self.rho * self.nu * self.alpha * self.beta\
                                        *  pow((F + b), self.beta -1) \
                                        + 1/24 *(2 - 3 * self.rho * self.rho)*\
                                               self.nu * self.nu) *  T)
        else:
            if self.nu  < 1e-6:
                # problem            
                x_zeta = zeta_k / self.nu
            elif self.rho == 1:
                x_zeta = np.log(K) * F **self.beta -K ** self.beta * F * self.nu\
                + K ** self.beta * self.alpha * self.beta -\
                                 K** self.beta * self.alpha * F **self.beta + \
                                 np.log(K)*self.beta * self.beta - K ** self.beta * F ** self.beta-\
                                                                       np.log(self.alpha)
            elif self.rho == -1:
                x_zeta =  np.log(K)*self.beta * self.beta- K ** self.beta * F ** self.beta + np.log(self.alpha)-\
                                  np.log(K) * F **self.beta + K ** self.beta * F * self.nu \
                                  + K ** self.beta * self.alpha * self.beta-\
                                   K** self.beta * self.alpha * F **self.beta       
            else:
                x_zeta = np.log((np.sqrt( 1 - 2 * self.rho * zeta_k+ zeta_k **2) - self.rho + zeta_k))/self.nu
            sub1 = 1/4 * self.rho * self.nu * self.alpha * self.beta * pow((F + b), (self.beta/2 - 1/2)) \
            * pow((K + b), (self.beta - 1)/2)
            sub2 = 1/24 *(2 - 3 * self.rho * self.rho) * self.nu * self.nu
            sigma_lognormal = 1 / x_zeta * np.log((F + b)/ (K + b)) * (1 +(g_k + sub1 +sub2 ) * T )
            
        return sigma_lognormal
    def SABR_vol_lognormal_by_moneyness(self, z, T, F, b=0):
        if self.beta == 1:
            g_k = 0
        elif self.beta == 0:
            g_k = self.alpha * self.alpha / (24 * (F + b) * (np.exp(z) + b))
        else:
            g_k = 1 / 24 * pow(self.beta - 1, 2) * pow((F + b), self.beta - 1) * pow((np.exp(z) + b), self.beta - 1) * self.alpha * self.alpha
        if np.all(z) != 0:
            sigma_by_z = self.alpha * pow((F + b), self.beta - 1) + 1 / 2 * (self.rho * self.nu - \
                                       (1 - self.beta) * self.alpha * pow(F + b, self.beta - 1)) * z + \
                       ((1 - self.beta) ** 2 * ((self.alpha * pow((F + b), self.beta - 1))) ** 2 + \
                        self.nu ** 2 * (2 - 3 * self.rho * self.rho)) * \
                       z * z / (12 * self.alpha * pow((F + b), self.beta - 1))
        else:
            # ==============================================================================
            #     # problem
            # ==============================================================================
            sigma_by_z = self.alpha * pow((F + b), self.beta - 1) * (1 + (g_k + 1 / 4 *
                                                              self.rho * self.nu * self.alpha * self.beta *\
                                                              pow((F + b), self.beta - 1) + 1 / 24 * (
                                                              2 - 3 * self.rho * self.rho) * \
                                                              self.nu * self.nu) * T)
        return sigma_by_z
    
    
                
            
            