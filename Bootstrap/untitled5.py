# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 14:00:21 2017

@author: jrl276
"""

import numpy as np

S1 = np.zeros([29,2,5])
ST = np.zeros([29,2,5])
S2 = np.zeros([29,29,2,5])

S1_conf = np.zeros([29,2,5])
ST_conf = np.zeros([29,2,5])
S2_conf = np.zeros([29,29,2,5])

l_10=np.load('Dam_10.txt.npz')
l_100=np.load('Dam_100.txt.npz')
l_1000=np.load('Dam_1000.txt.npz')

S1[:,:,0] = l_10['S1']
S1[:,:,1] = l_100['S1']
S1[:,:,2] = l_1000['S1']

S1_conf[:,:,0] = l_10['S1_conf']
S1_conf[:,:,1] = l_100['S1_conf']
S1_conf[:,:,2] = l_1000['S1_conf']

ST[:,:,0] = l_10['ST']
ST[:,:,1] = l_100['ST']
ST[:,:,2] = l_1000['ST']

ST_conf[:,:,0] = l_10['ST_conf']
ST_conf[:,:,1] = l_100['ST_conf']
ST_conf[:,:,2] = l_1000['ST_conf']

S2[:,:,:,0] = l_10['S2']
S2[:,:,:,1] = l_100['S2']
S2[:,:,:,2] = l_1000['S2']

S2_conf[:,:,:,0] = l_10['S2_conf']
S2_conf[:,:,:,1] = l_100['S2_conf']
S2_conf[:,:,:,2] = l_1000['S2_conf']