print 'Loading in numpy...'
import csv
import numpy as np
from scipy.stats import norm
#Load Data in
print 'Loading Data in...'
control_rate = np.loadtxt("control_rate.txt")
savings_rate = np.loadtxt("savings_rate.txt")
SAvars = np.loadtxt("./SALib-master/CDICE_SAvar.txt")

#Transform uniform (0,1) variates to log-normal for the climate sensitivity
SAvars[:,16] = np.exp(norm.ppf(SAvars[:,16],1.098001424,0.265206276))

npzfile=np.load("configuration.npz")
nominal = npzfile['nominal']
var_index = npzfile['var_index']

#Now assemble default params for input to CDICE
N = len(SAvars)
d = len(var_index)
default = np.hstack((nominal,control_rate,savings_rate))
M = len(default);

params = np.zeros((N,M))

print 'Running over N...'
for i in range(0,N):
    params[i,:] = default

# Include variable parameters
print 'Running over d...'
if d > 1:
    for i in range(0,d):
        params[:,var_index[i]] = SAvars[:,i]
else:
    params[:,var_index[0]] = SAvars

# Save to txt file
np.savetxt("CDICE_input.txt", params, fmt = '%10.5f', delimiter=' ')