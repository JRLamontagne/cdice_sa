print 'Loading in numpy...'
import csv
import numpy as np
from scipy.stats import norm
from policy_maker import make_opt2
from policy_maker import make_2C

#Load Data in
print 'Loading Data in...'
control_opt = np.loadtxt("control_rate.txt")
control_opt2 = np.loadtxt("control_rate_2C_full.txt")
savings_rate = np.loadtxt("savings_rate.txt")
SAvars = np.loadtxt("./SALib-master/CDICE_SAvar.txt")

#Transform uniform (0,1) variates to log-normal for the climate sensitivity
SAvars[:,16] = np.exp(norm.ppf(SAvars[:,16],1.098001424,0.265206276))

npzfile=np.load("configuration.npz")
nominal = npzfile['nominal']
var_index = npzfile['var_index']

#Now assemble default params for input to CDICE
N = len(SAvars)
d = 0#len(var_index)
default = np.hstack((nominal,control_opt,savings_rate))
M = len(default);

params = np.zeros((N,M))

hd = np.zeros((N,8))
hd[:,6] = 57
hd[:,7] = 0.956942624
hd[:,0] = SAvars[:,(d+0)]
hd[:,1] = SAvars[:,(d+1)]
hd[:,2] = SAvars[:,(d+2)]
hd[:,3] = 9
hd[:,4] = 0.2
hd[:,5] = 1.0

np.savetxt("test.txt",hd, fmt = '%10.5f', delimiter=' ')

print 'Running over N...'
for i in range(0,N):
    params[i,:] = default
    control_rate = make_2C(hd[i,:],control_opt2)
    params[i,len(nominal):len(nominal)+59] = control_rate

# Include variable parameters
print 'Running over d...'
#if d > 1:
#    for i in range(0,d):
#        params[:,var_index[i]] = SAvars[:,i]
#else:
#    params[:,var_index[0]] = SAvars

# Save to txt file
np.savetxt("CDICE_input.txt", params, fmt = '%10.5f', delimiter=' ')