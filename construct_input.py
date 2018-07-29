#This code reads the configuration.npz file and generates the files needed for
#the SALib Saltelli sampling code
import numpy as np

npzfile=np.load("configuration.npz")

fix_index = npzfile['fix_index']
fix_nominal = npzfile['fix_nominal']
fix_names = npzfile['fix_names']
num_fix = npzfile['num_fix']

var_index = npzfile['var_index']
var_lo_bounds = npzfile['var_lo_bounds']
var_hi_bounds = npzfile['var_hi_bounds']
var_names = npzfile['var_names']
num_vars = npzfile['num_vars']

f=open("./SALib-master/SALib_Input.txt","w")
for i in range(0,num_vars):
    f.write(var_names[i])
    f.write(" ")
    f.write(str(var_lo_bounds[i]))
    f.write(" ")
    f.write(str(var_hi_bounds[i]))
    f.write("\n")
#        print(i)
f.close()
print(float(num_vars))