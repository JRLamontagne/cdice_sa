import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib as mpl

# uncertain factors
LHsamples = pd.read_csv('../Final_Data/1000/CDICE_SAvar.txt', sep=' ', \
    names=["pop0","popadj","popasym","q0","ga0","dela","gsigma1","dsig","eland0","deland","e0","mat0","b12","b23","t2xco2","fex0","fex1","fco22x","a1","a2","a3","Expcost2","pback","gback","a"])
    
paramBounds = np.loadtxt('uncertain_params.txt',usecols=[1,2])
    
#normalize LHsamples
for i, col in enumerate(LHsamples.columns):
    LHsamples[col] = (LHsamples[col] - paramBounds[i,0])/(paramBounds[i,1] - paramBounds[i,0])
    
# color maps for plotting
dot_cmap = mpl.colors.ListedColormap(np.array([[227,26,28],[166,206,227]])/255.0)
contour_cmap = mpl.cm.get_cmap('RdBu')
class_cmap = mpl.colors.ListedColormap(np.array([[251,154,153],[31,120,180]])/255.0)

#############################################################################################
# fit logistic regression to tolerable window attainment
dta = pd.read_csv('../Final_Data/1000/NPV_Abat_150.txt',\
    sep='   ', names=['NPV_Abat'])
dta['NPV_Abat'] = 100*np.loadtxt('../Final_Data/1000/NPV_Abat_150.txt')/np.loadtxt('../Final_Data/1000/GWP_DISC_150.txt')
dta['NPV_Dam'] = 100*np.loadtxt('../Final_Data/1000/NPV_Dam_150.txt')/np.loadtxt('../Final_Data/1000/GWP_DISC_150.txt')
dta['Tatm_100'] = np.loadtxt('../Final_Data/1000/Tatm_100.txt')

dta['LinearPVAbat'] = ((dta.NPV_Abat < 2.00) & (dta.NPV_Dam < 2.00) & (dta.Tatm_100 < 2.00)).astype(int)
#dta['LinearPVAbat'] = (((dta.NPV_Abat + dta.NPV_Dam) < 5.00) & (dta.Tatm_100 < 2.00)).astype(int)
# concatenate values of uncertain factors and intercept column of 1s
dta = pd.concat([LHsamples,dta],axis=1)
dta['Intercept'] = np.ones(np.shape(LHsamples)[0])
R2 = np.zeros([25,1])
# get columns of predictors
for i in range(25):
    cols = dta.columns.tolist()
    cols = cols[-1:] + cols[i:i+1]

    #fit logistic regression
    logit = sm.Logit(dta['LinearPVAbat'], dta[cols])
    result = logit.fit()
    R2[i,0] = result.prsquared
#############################################################################################