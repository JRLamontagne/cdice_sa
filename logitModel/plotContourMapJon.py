# -*- coding: utf-8 -*-
"""
Created on Thu Aug 03 14:08:34 2017

@author: jrl276
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm

def plotContourMapJon(result, LHsamples, dta, variable, contour_cmap, dot_cmap, xgrid, y1grid, levels, \
    xvar, y1var, xlabel, ylabel1, prb, title, filename):
    
    # find probability of success for x=xgrid, y1var=ygrid and y2var=1
    X, Y = np.meshgrid(xgrid, y1grid)
    x = X.flatten()
    y = Y.flatten()
    grid = np.column_stack([np.ones(len(x)),np.ones(len(x))*prb,x,y])
    z = result.predict(grid)
    Z = np.reshape(z, np.shape(X))
            
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(1,1,1)
    contourset = ax.contourf(X, Y, Z, levels, cmap=contour_cmap)
    x = LHsamples[xvar].values
    y = LHsamples[y1var].values
    #ax.scatter(LHsamples[xvar].values,LHsamples[y1var].values,\
    #    facecolor=dta[variable].values,edgecolor='none',cmap=dot_cmap,s=10)
    ax.set_xlim(np.min(X),np.max(X))
    ax.set_ylim(np.min(Y),np.max(Y))
    ax.set_xlabel(xlabel,fontsize=24)
    ax.set_ylabel(ylabel1,fontsize=24)
    ax.tick_params(axis='both',labelsize=18)
    
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(contourset, cax=cbar_ax)
    cbar_ax.set_ylabel('Probability of Success',fontsize=20)
    yticklabels = cbar.ax.get_yticklabels()
    cbar.ax.set_yticklabels(yticklabels,fontsize=18)
    
    fig.suptitle(title,fontsize=24)
    fig.set_size_inches([14.4375, 6.7])
    fig.savefig(filename)
    fig.clf()
    
    return None