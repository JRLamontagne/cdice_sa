# -*- coding: utf-8 -*-
"""
Script tests policy maker
"""
import numpy as np
from policy_maker import make_pol_exp
import matplotlib.pyplot as plt
import seaborn as sns

a = np.zeros([100,1])
b = np.zeros([100,1])
c = np.zeros([3,1])
c[0] = 25
c[1] = 29
c[2] = 33
params = np.zeros([3,1])
policy = np.zeros([100*100*3,99])
cnt = 0
for i in range(100):
    a[i] = 0.4*float(i)/100
    params[0] = a[i]
    for j in range(100):
        b[j] = (0.405-0.075)*float(j)/100+0.075
        #3b[j] = np.exp((1+2.59)*float(j)/100-2.59)
        params[1] = b[j]
        for k in range(3):
            params[2] = c[k]
            
            policy[cnt,:] = make_pol_exp(params)
            cnt = cnt + 1

plt.plot(np.transpose(policy))
'''
year = range(99)
abate = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.15,1.2,1.25]
X = np.zeros([99,25])
XX = np.zeros([99,25])
for i in range(99):
    for j in range(25):
        X[i,j] = sum(policy[:,i]<abate[j])

    XX[:,0] = X[:,0]
    for j in range(24,0,-1):
        XX[i,j] = X[i,j] - X[i,j-1]
    print i
'''
        