# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 16:21:53 2017

@author: jrl276
"""
# -*- coding: utf-8 -*-
def sobol_index(input_filename,N,B):
    import sys
    import numpy as np
    import time

    sys.path.append('./SALib-master')

    from SALib.analyze import sobol
    from SALib.util import read_param_file
    start = time.time()
    Y_names = ['NPV_Abat.txt','NPV_Dam.txt']

    #Load in data
    filename = "SALib/%s" % input_filename
    problem = read_param_file(filename,' ')

    S1 = np.zeros((29,2))
    S1_conf = np.zeros((29,2))
    S2 = np.zeros((29,29,2))
    S2_conf = np.zeros((29,29,2))
    ST = np.zeros((29,2))
    ST_conf = np.zeros((29,2))

    cnt = 0

    for i in range(2):
        filename = Y_names[i]
        Y = np.loadtxt(filename)
        Y = Y[:,17] #2200
        Si = sobol.analyze(problem, Y, calc_second_order=True, conf_level=0.95, num_resamples=B, print_to_console=False,parallel=True,
            n_processors=15)

        S1[:,cnt] = Si.get("S1")
        S1_conf[:,cnt] = Si.get("S1_conf")
        S2[:,:,cnt] = Si.get("S2")
        S2_conf[:,:,cnt] = Si.get("S2_conf")
        ST[:,cnt] = Si.get("ST")
        ST_conf[:,cnt] = Si.get("ST_conf")
        cnt = cnt+1
        
    filename = "Dam_%i.txt" % N
    np.savez(filename,S1=S1,S1_conf=S1_conf,S2=S2,S2_conf=S2_conf,ST=ST,ST_conf=ST_conf)
    end = time.time()
    print(end-start)
def sobol_sample(input_filename,output_filename,N):
    import sys
    import numpy as np
    
    sys.path.append('./SALib-master')
    
    from SALib.sample import saltelli
    from SALib.util import read_param_file
    
    filename = "SALib-master/%s" % input_filename
    problem = read_param_file(filename,' ')
    
    param_values = saltelli.sample(problem, N, calc_second_order=True)
    filename = "SALib-master/%s" % output_filename
    np.savetxt(filename,param_values,delimiter=' ')
    