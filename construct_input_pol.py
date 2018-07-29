#This code reads the configuration.npz file and generates the files needed for
#the SALib Saltelli sampling code
def construct_input(type_analysis):
    #type_analysis is an indicator of what type of anlysis we are doing.
    #1 is both DICE and policy
    #2 is only DICE
    #3 is only policy
    print type_analysis
    import numpy as np
    
    f=open("./SALib-master/SALib_Input.txt","w")
    
    if type_analysis == 1 or type_analysis == 2:
        npzfile=np.load("configuration.npz")

        var_lo_bounds = npzfile['var_lo_bounds']
        var_hi_bounds = npzfile['var_hi_bounds']
        var_names = npzfile['var_names']
        num_vars = npzfile['num_vars']
        
        for i in range(0,num_vars):
            f.write(var_names[i])
            f.write(" ")
            f.write(str(var_lo_bounds[i]))
            f.write(" ")
            f.write(str(var_hi_bounds[i]))
            f.write("\n")
    else:
        npzfile=np.load("configuration_pol.npz")

        var_lo_bounds = npzfile['var_lo_bounds']
        var_hi_bounds = npzfile['var_hi_bounds']
        var_names = npzfile['var_names']
        num_vars_pol = npzfile['num_vars']

        for i in range(0,num_vars_pol):
            f.write(var_names[i])
            f.write(" ")
            f.write(str(var_lo_bounds[i]))
            f.write(" ")
            f.write(str(var_hi_bounds[i]))
            f.write("\n")
        print 'shit'
    if type == 1:
        npzfile=np.load("configuration_pol.npz")

        var_lo_bounds = npzfile['var_lo_bounds']
        var_hi_bounds = npzfile['var_hi_bounds']
        var_names = npzfile['var_names']
        num_vars_pol = npzfile['num_vars']

        for i in range(0,num_vars_pol):
            f.write(var_names[i])
            f.write(" ")
            f.write(str(var_lo_bounds[i]))
            f.write(" ")
            f.write(str(var_hi_bounds[i]))
            f.write("\n")
        print 'fuck'
    f.close()
    print(float(num_vars)+float(num_vars_pol))