def cdice_input_maker(type_analysis):
    #type_analysis is an indicator of what type of anlysis we are doing.
    #1 is both DICE and policy
    #2 is only DICE
    #3 is only policy

    import numpy as np
    from scipy.stats import norm
    from policy_maker import make_pol_exp
    from policy_maker import make_2C

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
    
    default = np.hstack((nominal,control_opt,savings_rate))
    M = len(default);

    params = np.zeros((N,M))
    print 'Running over N...'
    if type_analysis == 2:
        for i in range(0,N):
            params[i,:] = default
    else:
        policy_params = np.zeros((N,3))
        policy_params[:,0] = SAvars[:,(d+0)]
        policy_params[:,1] = SAvars[:,(d+1)]
        policy_params[:,2] = SAvars[:,(d+2)]
    
        for i in range(0,N):
            params[i,:] = default
            #control_rate = make_pol_exp(policy_params[i,:])
            control_rate = make_2C(policy_params[i,:])
            params[i,len(nominal):len(nominal)+99] = control_rate
    
    print 'Running over d...'    
    if type_analysis == 1 or type_analysis == 2:
        # Include variable parameters
        if d > 1:
            for i in range(0,d):
                params[:,var_index[i]] = SAvars[:,i]
        else:
            params[:,var_index[0]] = SAvars

    # Save to txt file
    np.savetxt("CDICE_input.txt", params, fmt = '%10.5f', delimiter=' ')