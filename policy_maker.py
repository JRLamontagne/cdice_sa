# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 23:07:04 2017

@author: jrl276
"""
def make_policy_single(func_form,a,year):
    import numpy as np
    a = float(a)
    BAU = [0.0299,0.0323,0.0349,0.0377,0.0408,0.0441,0.0476,0.0515,0.0556,0.0601,0.0650,0.0702,0.0759,0.0821,0.0887,0.0959,0.1036,0.1120]
    YEARS={'2020':1,'2030':3,'2040':5,'2050':7,'2060':9,'2070':11,'2075':12}
    beta0 = BAU[0]
    bau = BAU[YEARS[str(year)]]
    year = float(YEARS[str(year)])
    policy = np.zeros([99])
    
    for t in range(1,100):
        if func_form == 'linear':
            a1 = (1/year)*(a+bau*(1-a)-beta0)
            policy[t-1] = np.min([a1*t+beta0,1.0])
        elif func_form == 'quad':
            a1 = (1/(year**2.0))*(a+bau*(1-a)-beta0)
            policy[t-1] = np.min([a1*(t**2.0)+beta0,1.0])
        elif func_form == 'exp':
            a1 = (1/year)*np.log((1/beta0)*(a+bau*(1-a)))
            policy[t-1] = np.min([beta0*np.exp(a1*t),1.0])
        else:
            raise ValueError('Function Form not Recognized')
        
        if t >= 29 and t<49:
            policy[t-1] = 1.2
        elif t >= 49:
            policy[t-1] = 1.0
    
    return policy