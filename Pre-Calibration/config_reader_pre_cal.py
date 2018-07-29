#This code reads the configuration.txt file, and pulls out which factors are
    #variables and which are fixed
import csv
import numpy as np

cnt = 0;

fix_index = []
fix_names = []
fix_nominal = []
num_fix = 0;

var_index =[]
var_names = []
var_lo_bounds = []
var_hi_bounds = []
num_vars = 0

nominal = []
    
filename = "configuration.txt"
fieldnames=["Name","Include","Nominal","Low","High"]
with open(filename) as csvfile:
    dialect = csv.Sniffer().sniff(csvfile.read(1024))
    csvfile.seek(0)
    reader = csv.DictReader(csvfile, fieldnames=fieldnames, dialect=dialect)
    for row in reader:
        if row['Name'].strip().startswith('#'):
            pass
        else:
            if float(row['Include']) == 0:
                fix_index.append(cnt-1)
                fix_names.append(row['Name'])
                fix_nominal.append(float(row['Nominal']))
                num_fix = num_fix+1
            else:
                var_index.append(cnt-1)
                var_names.append(row['Name'])
                var_lo_bounds.append(float(row['Low']))
                var_hi_bounds.append(float(row['High']))
                num_vars = num_vars+1
            nominal.append(float(row['Nominal']))
        cnt = cnt+1
np.savez("configuration.npz",fix_index=fix_index,fix_names=fix_names,fix_nominal=fix_nominal,num_fix=num_fix,var_index=var_index,var_names=var_names,var_lo_bounds=var_lo_bounds,var_hi_bounds=var_hi_bounds,num_vars=num_vars,nominal=nominal)