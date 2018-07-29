# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 17:06:19 2015

@author: jrl276
"""

import numpy as np
#NPV_Dam
S1 = np.zeros((29,60))
S2 = np.zeros((29,29,60))
ST = np.zeros((29,60))

for i in range(0,60):
    filename = "NPV_Dam_%i.txt.npz" % i
    npzfile=np.load(filename)
    hold_S1 = npzfile['S1']
    hold_S2 = npzfile['S2']
    hold_ST = npzfile['ST']
    for j in range(0,29):
        S1[j,i] = hold_S1[j]
        ST[j,i] = hold_ST[j]
        for k in range(0,29):
            S2[j,k,i] = hold_S2[j,k]

np.savetxt("NPV_Dam_S1.txt", S1, fmt = '%10.5f', delimiter=' ')
np.savetxt("NPV_Dam_ST.txt", ST, fmt = '%10.5f', delimiter=' ')

for i in range(0,60):
    filename = "NPV_Dam_2_%i.txt" %i
    np.savetxt(filename, S2[:,:,i], fmt = '%10.5f', delimiter=' ')
#Dam
S1 = np.zeros((29,60))
S2 = np.zeros((29,29,60))
ST = np.zeros((29,60))

for i in range(0,60):
    filename = "Dam_%i.txt.npz" % i
    npzfile=np.load(filename)
    hold_S1 = npzfile['S1']
    hold_S2 = npzfile['S2']
    hold_ST = npzfile['ST']
    for j in range(0,29):
        S1[j,i] = hold_S1[j]
        ST[j,i] = hold_ST[j]
        for k in range(0,29):
            S2[j,k,i] = hold_S2[j,k]

np.savetxt("Dam_S1.txt", S1, fmt = '%10.5f', delimiter=' ')
np.savetxt("Dam_ST.txt", ST, fmt = '%10.5f', delimiter=' ')

for i in range(0,60):
    filename = "Dam_2_%i.txt" %i
    np.savetxt(filename, S2[:,:,i], fmt = '%10.5f', delimiter=' ')
#Abat
S1 = np.zeros((29,60))
S2 = np.zeros((29,29,60))
ST = np.zeros((29,60))

for i in range(0,60):
    filename = "Abat_%i.txt.npz" % i
    npzfile=np.load(filename)
    hold_S1 = npzfile['S1']
    hold_S2 = npzfile['S2']
    hold_ST = npzfile['ST']
    for j in range(0,29):
        S1[j,i] = hold_S1[j]
        ST[j,i] = hold_ST[j]
        for k in range(0,29):
            S2[j,k,i] = hold_S2[j,k]

np.savetxt("Abat_S1.txt", S1, fmt = '%10.5f', delimiter=' ')
np.savetxt("Abat_ST.txt", ST, fmt = '%10.5f', delimiter=' ')

for i in range(0,60):
    filename = "Abat_2_%i.txt" %i
    np.savetxt(filename, S2[:,:,i], fmt = '%10.5f', delimiter=' ')
#NPV_Abat
S1 = np.zeros((29,60))
S2 = np.zeros((29,29,60))
ST = np.zeros((29,60))

for i in range(0,60):
    filename = "NPV_Abat_%i.txt.npz" % i
    npzfile=np.load(filename)
    hold_S1 = npzfile['S1']
    hold_S2 = npzfile['S2']
    hold_ST = npzfile['ST']
    for j in range(0,29):
        S1[j,i] = hold_S1[j]
        ST[j,i] = hold_ST[j]
        for k in range(0,29):
            S2[j,k,i] = hold_S2[j,k]

np.savetxt("NPV_Abat_S1.txt", S1, fmt = '%10.5f', delimiter=' ')
np.savetxt("NPV_Abat_ST.txt", ST, fmt = '%10.5f', delimiter=' ')

for i in range(0,60):
    filename = "NPV_Abat_2_%i.txt" %i
    np.savetxt(filename, S2[:,:,i], fmt = '%10.5f', delimiter=' ')
# SCC
S1 = np.zeros((29,60))
S2 = np.zeros((29,29,60))
ST = np.zeros((29,60))

for i in range(0,60):
    filename = "SCC_%i.txt.npz" % i
    npzfile=np.load(filename)
    hold_S1 = npzfile['S1']
    hold_S2 = npzfile['S2']
    hold_ST = npzfile['ST']
    for j in range(0,29):
        S1[j,i] = hold_S1[j]
        ST[j,i] = hold_ST[j]
        for k in range(0,29):
            S2[j,k,i] = hold_S2[j,k]

np.savetxt("SCC_S1.txt", S1, fmt = '%10.5f', delimiter=' ')
np.savetxt("SCC_ST.txt", ST, fmt = '%10.5f', delimiter=' ')

for i in range(0,60):
    filename = "SCC_2_%i.txt" %i
    np.savetxt(filename, S2[:,:,i], fmt = '%10.5f', delimiter=' ')