def cdice_input_maker(type_analysis,func_form,year):
    #type_analysis is an indicator of what type of anlysis we are doing.
    #1 is both DICE and policy
    #2 is only DICE
    #3 is only policy
    #func_form indicates what policy approximator to use, options are linear, quad, or exp
    #year indicates the year in which we have uniform sampling between BAU and 1.

    import numpy as np
    from scipy.stats import norm
    from policy_maker import make_policy_single

    #Load Data in
    
    control_opt = np.loadtxt("control_rate.txt")
    savings_rate = np.loadtxt("savings_rate.txt")
    SAvars = np.loadtxt("./SALib-master/CDICE_SAvar.txt")
	
    #Transform uniform (0,1) variates to log-normal for the climate sensitivity
    if type_analysis == 1 or type_analysis == 2:
        SAvars[:,16] = np.exp(norm.ppf(SAvars[:,16],1.098001424,0.265206276)) #confirm 16 index still right.

    npzfile=np.load("configuration.npz")
    nominal = npzfile['nominal']#need to transform txco2 to log-normal
    var_index = npzfile['var_index']

    #Now assemble default params for input to CDICE
    N = len(SAvars)
    if type_analysis == 3:
        d = 0
    else:
        d = len(var_index)
        npzfile = np.load("configuration_pol.npz")
        var_index2 = npzfile['var_index']
        k = len(var_index2)
        #print k
    
    default = np.hstack((nominal,control_opt,savings_rate))
    M = len(default);

    params = np.zeros((N,M))
    print 'Running over N...'
    if type_analysis == 2:
        for i in range(0,N):
            params[i,:] = default    
    else:
        policy_params = np.zeros((N,k))
        print np.shape(SAvars)[0]
        print np.shape(SAvars)[1]
        print d
        print k
        print np.shape(policy_params)[0]
        print np.shape(policy_params)[1]
        
        for i in range(k):
            print i
            policy_params[:,i] = SAvars[:,(d+i)]
    
        for i in range(0,N):
            params[i,:] = default
            control_rate = make_policy_single(func_form,policy_params[i,0],year)
            params[i,len(nominal):len(nominal)+99] = control_rate
    
    print 'Running over d...'    
    if type_analysis == 1 or type_analysis == 2:
        # Include variable parameters
        if d > 1:
            print d
            for i in range(0,d):
                params[:,var_index[i]] = SAvars[:,i]
        else:
            params[:,var_index[0]] = SAvars

    # Save to txt file
    np.savetxt("CDICE_input.txt", params, fmt = '%10.5f', delimiter=' ')