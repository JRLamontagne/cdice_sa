import numpy as np
from scipy.stats import norm
import pandas as pd
import statsmodels.api as sm
import matplotlib as mpl
from plotContourMap import plotContourMap
from plotContourMapJon import plotContourMapJon
# uncertain factors
LHsamples = pd.read_csv('CDICE_SAvar.txt', sep=' ', \
    names=["elasmu","prstp","pop0","popadj","popasym","q0","ga0","dela","gsigma1","dsig","eland0","deland","e0","mat0","b12","b23","t2xco2","fex0","fex1","fco22x","a1","a2","a3","expcost2","pback","gback","a","b","c"])
    
paramBounds = np.loadtxt('uncertain_params.txt',usecols=[1,2])
    
# normalize LHsamples
#for i, col in enumerate(LHsamples.columns):
#    LHsamples[col] = (LHsamples[col] - paramBounds[i,0])/(paramBounds[i,1] - paramBounds[i,0])
    
# color maps for plotting
dot_cmap = mpl.colors.ListedColormap(np.array([[227,26,28],[166,206,227]])/255.0)
contour_cmap = mpl.cm.get_cmap('RdBu')
class_cmap = mpl.colors.ListedColormap(np.array([[251,154,153],[31,120,180]])/255.0)

#############################################################################################
# fit logistic regression to flood successes of best flood solution
dta = pd.read_csv('NPV_Abat_100.txt',\
    sep='   ', names=['NPV_Abat'])
dta['NPV_Abat'] = 100*np.loadtxt('NPV_Abat_100.txt')/np.loadtxt('GWP_DISC_100.txt')
dta['NPV_Dam'] = 100*np.loadtxt('NPV_Dam_100.txt')/np.loadtxt('GWP_DISC_100.txt')
dta['Tatm_100'] = np.loadtxt('Tatm_100.txt')

dta['ExpPVAbat'] = (((dta.NPV_Abat + dta.NPV_Dam) < 5.00) & (dta.Tatm_100 < 2.00)).astype(int)
# concatenate values of uncertain factors and intercept column of 1s
dta = pd.concat([LHsamples,dta],axis=1)
dta['Intercept'] = np.ones(np.shape(LHsamples)[0])
R2 = np.zeros([29,1])
# get columns of predictors
for i in range(29):
    cols = dta.columns.tolist()
    cols = cols[-1:] + cols[i:i+1]

    #fit logistic regression
    logit = sm.Logit(dta['ExpPVAbat'], dta[cols])
    result = logit.fit()
    R2[i,0] = result.prsquared  
#############################################################################################
# fit logistic regression to flood successes of best flood solutio
dta = pd.read_csv('NPV_Abat_100.txt',\
    sep='   ', names=['NPV_Abat'])
dta['NPV_Abat'] = 100*np.loadtxt('NPV_Abat_100.txt')/np.loadtxt('GWP_DISC_100.txt')
dta['NPV_Dam'] = 100*np.loadtxt('NPV_Dam_100.txt')/np.loadtxt('GWP_DISC_100.txt')
dta['Tatm_100'] = np.loadtxt('Tatm_100.txt')

dta['ExpPVAbat'] = (((dta.NPV_Abat + dta.NPV_Dam) < 5.00) & (dta.Tatm_100 < 2.00)).astype(int)
# concatenate values of uncertain factors and intercept column of 1s
dta = pd.concat([LHsamples,dta],axis=1)
dta['Intercept'] = np.ones(np.shape(LHsamples)[0])
R2 = np.zeros([29,1])
# get columns of predictors
for i in range(29):
    if i != 16:
        cols = dta.columns.tolist()
        cols = cols[-1:] + cols[16:17] + cols[i:i+1]

        #fit logistic regression
        logit = sm.Logit(dta['ExpPVAbat'], dta[cols])
        result = logit.fit()
        R2[i,0] = result.prsquared
#############################################################################################
# fit logistic regression to flood successes of best flood solutio
dta = pd.read_csv('NPV_Abat_100.txt',\
    sep='   ', names=['NPV_Abat'])
dta['NPV_Abat'] = 100*np.loadtxt('NPV_Abat_100.txt')/np.loadtxt('GWP_DISC_100.txt')
dta['NPV_Dam'] = 100*np.loadtxt('NPV_Dam_100.txt')/np.loadtxt('GWP_DISC_100.txt')
dta['Tatm_100'] = np.loadtxt('Tatm_100.txt')

dta['ExpPVAbat'] = (((dta.NPV_Abat + dta.NPV_Dam) < 5.00) & (dta.Tatm_100 < 2.00)).astype(int)
# concatenate values of uncertain factors and intercept column of 1s
dta = pd.concat([LHsamples,dta],axis=1)
dta['Intercept'] = np.ones(np.shape(LHsamples)[0])
R2 = np.zeros([29,1])
# get columns of predictors
for i in range(29):
    if i != 16 and i != 28:
        cols = dta.columns.tolist()
        cols = cols[-1:] + cols[16:17] + cols[28:29] + cols[i:i+1]

        #fit logistic regression
        logit = sm.Logit(dta['ExpPVAbat'], dta[cols])
        result = logit.fit()
        R2[i,0] = result.prsquared  
############################################################################################
cols = dta.columns.tolist()
cols = cols[-1:] + cols[16:17] + cols[27:28] + cols[28:29]
logit = sm.Logit(dta['ExpPVAbat'], dta[cols])
result = logit.fit()
R2 = result.prsquared

# countour map
xgrid = np.arange(np.min(LHsamples['t2xco2']), np.max(LHsamples['t2xco2'])+0.001,0.001)
y1grid = np.arange(np.min(LHsamples['b']), np.max(LHsamples['b'])+0.001,0.001)
y2grid = np.arange(np.min(LHsamples['c']), np.max(LHsamples['c'])+0.001,0.001)
levels=np.arange(0,1.05,0.1)
#levels = [0,0.95,1.00]
plotContourMap(result, LHsamples, dta, 'ExpPVAbat', contour_cmap, dot_cmap, xgrid, y1grid, y2grid, levels, \
    't2xco2', 'b', 'c', 't2xco2', 'b', 'c', \
    'PV Abatement + PV Damages > 3% and Tatm < 2', 'PV_Abat_Dam_GT5_Tarm_LT2.png')
############################################################################################
cols = dta.columns.tolist()
cols = cols[-1:] + cols[16:17] + cols[27:28] + cols[28:29]
logit = sm.Logit(dta['ExpPVAbat'], dta[cols])
result = logit.fit()
R2 = result.prsquared

# countour map
xgrid = np.arange(np.min(LHsamples['b']), np.max(LHsamples['b'])+0.001,0.001)
ygrid = np.arange(np.min(LHsamples['c']), np.max(LHsamples['c'])+0.001,0.001)
levels=np.arange(0,1.05,0.10)
#levels = [0,0.95,1.00]
plotContourMapJon(result, LHsamples, dta, 'ExpPVAbat', contour_cmap, dot_cmap, xgrid, ygrid, levels, \
    'b', 'c', 'b', 'c', 0.01, \
    'PV Abatement + PV Damages > 3% and Tatm < 2', 'PV_Abat_Dam_GT5_Tarm_LT2.png')
############################################################################################
## fit logistic regression to flood successes of best flood solutio
#dta = pd.read_csv('NPV_Abat_100.txt',\
#    sep='   ', names=['NPV_Abat'])
#dta['NPV_Abat'] = 100*np.loadtxt('NPV_Abat_100.txt')/np.loadtxt('GWP_DISC_100.txt')
#dta['NPV_Dam'] = 100*np.loadtxt('NPV_Dam_100.txt')/np.loadtxt('GWP_DISC_100.txt')
#dta['Tatm_100'] = np.loadtxt('Tatm_100.txt')
#
#dta['ExpPVAbat'] = (((dta.NPV_Abat + dta.NPV_Dam) < 3.00) & (dta.Tatm_100 < 2.00)).astype(int)
## concatenate values of uncertain factors and intercept column of 1s
#dta = pd.concat([LHsamples,dta],axis=1)
#dta['Intercept'] = np.ones(np.shape(LHsamples)[0])
#R2 = np.zeros([29,1])
## get columns of predictors
#for i in range(29):
#    if i != 16:
#        cols = dta.columns.tolist()
#        cols = cols[-1:] + cols[16:17] + cols[i:i+1]
#
#        #fit logistic regression
#        logit = sm.Logit(dta['ExpPVAbat'], dta[cols])
#        result = logit.fit()
#        R2[i,0] = result.prsquared  
## fit logistic regression to flood successes of best flood solutio
#dta = pd.read_csv('NPV_Abat_100.txt',\
#    sep='   ', names=['NPV_Abat'])
#dta['NPV_Abat'] = 100*np.loadtxt('NPV_Abat_100.txt')/np.loadtxt('GWP_DISC_100.txt')
#dta['NPV_Dam'] = 100*np.loadtxt('NPV_Dam_100.txt')/np.loadtxt('GWP_DISC_100.txt')
#dta['Tatm_100'] = np.loadtxt('Tatm_100.txt')
#
#dta['ExpPVAbat'] = (((dta.NPV_Abat + dta.NPV_Dam) < 4.00) & (dta.Tatm_100 < 2.00)).astype(int)
## concatenate values of uncertain factors and intercept column of 1s
#dta = pd.concat([LHsamples,dta],axis=1)
#dta['Intercept'] = np.ones(np.shape(LHsamples)[0])
#R2 = np.zeros([29,1])
## get columns of predictors
#for i in range(29):
#    if i != 16 and i != 27 and i != 8:
#        cols = dta.columns.tolist()
#        cols = cols[-1:] + cols[16:17] + cols[27:28] + cols[8:9] + cols[i:i+1]
#
#        #fit logistic regression
#        logit = sm.Logit(dta['ExpPVAbat'], dta[cols])
#        result = logit.fit()
#        R2[i,0] = result.prsquared  
##############################################################################################
## fit logistic regression to flood successes of best flood solutio
#dta = pd.read_csv('NPV_Abat_100.txt',\
#    sep='   ', names=['NPV_Abat'])
#dta['NPV_Abat'] = 100*np.loadtxt('NPV_Abat_100.txt')/np.loadtxt('GWP_DISC_100.txt')
#dta['NPV_Dam'] = 100*np.loadtxt('NPV_Dam_100.txt')/np.loadtxt('GWP_DISC_100.txt')
#dta['Tatm_100'] = np.loadtxt('Tatm_100.txt')
#
#dta['ExpPVAbat'] = (((dta.NPV_Abat + dta.NPV_Dam) < 4.00) & (dta.Tatm_100 < 2.00)).astype(int)
## concatenate values of uncertain factors and intercept column of 1s
#dta = pd.concat([LHsamples,dta],axis=1)
#dta['Intercept'] = np.ones(np.shape(LHsamples)[0])
#R2 = np.zeros([29,1])
## get columns of predictors
#for i in range(29):
#    if i != 16:
#        cols = dta.columns.tolist()
#        cols = cols[-1:] + cols[16:17] + cols[i:i+1]
#
#        #fit logistic regression
#        logit = sm.Logit(dta['ExpPVAbat'], dta[cols])
#        result = logit.fit()
#        R2[i,0] = result.prsquared  
## fit logistic regression to flood successes of best flood solutio
#dta = pd.read_csv('NPV_Abat_100.txt',\
#    sep='   ', names=['NPV_Abat'])
#dta['NPV_Abat'] = 100*np.loadtxt('NPV_Abat_100.txt')/np.loadtxt('GWP_DISC_100.txt')
#dta['NPV_Dam'] = 100*np.loadtxt('NPV_Dam_100.txt')/np.loadtxt('GWP_DISC_100.txt')
#dta['Tatm_100'] = np.loadtxt('Tatm_100.txt')
#
#dta['ExpPVAbat'] = (((dta.NPV_Abat + dta.NPV_Dam) < 4.00) & (dta.Tatm_100 < 2.00)).astype(int)
## concatenate values of uncertain factors and intercept column of 1s
#dta = pd.concat([LHsamples,dta],axis=1)
#dta['Intercept'] = np.ones(np.shape(LHsamples)[0])
#R2 = np.zeros([29,1])
## get columns of predictors
#for i in range(29):
#    if i != 16 and i != 27 and i != 8 and i != 28:
#        cols = dta.columns.tolist()
#        cols = cols[-1:] + cols[16:17] + cols[27:28] + cols[8:9] + cols[28:29] + cols[i:i+1]
#
#        #fit logistic regression
#        logit = sm.Logit(dta['ExpPVAbat'], dta[cols])
#        result = logit.fit()
#        R2[i,0] = result.prsquared  
##################################
#        