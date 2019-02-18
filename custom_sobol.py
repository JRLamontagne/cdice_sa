# -*- coding: utf-8 -*-
"""
Created on Sun Jul 02 23:16:01 2017

@author: jrl276
"""

def custom_sobol(input_filename,n_resample,start_ind,end_ind,N,D,output_filename):
    import numpy as np
    H=np.loadtxt(input_filename)
    Y_hld = H[:,start_ind:end_ind+1]
    del H
    
    step=2*D+2
    tsteps = end_ind-start_ind+1
    
    ST = np.zeros([D,tsteps])
    S1 = np.zeros([D,tsteps])
    S2 = np.zeros([D,D,tsteps])
    AB = np.zeros((N, D))
    BA = np.zeros((N, D))
    
    if n_resample > 0:
        ST_conf = np.zeros([D,tsteps])
        S1_conf = np.zeros([D,tsteps])
        S2_conf = np.zeros([D,D,tsteps])
    
    for i in range(0,tsteps):
        Y = Y_hld[:,i]
        Y = (Y - Y.mean())/Y.std()
        A=Y[0:Y.size:step]
        B = Y[(step - 1):Y.size:step]

        for j in range(D):
            AB[:, j] = Y[(j + 1):Y.size:step]
            BA[:, j] = Y[(j + 1 + D):Y.size:step]

        for j in range(D):
            S1[j,i]=np.mean(B * (AB[:,j] - A), axis=0) / np.var(np.r_[A, B], axis=0)
            ST[j,i]=0.5 * np.mean((A - AB[:,j]) ** 2, axis=0) / np.var(np.r_[A, B], axis=0)

            for k in range(j+1,D):
                Vjk = np.mean(BA[:,j] * AB[:,k] - A * B, axis=0) / np.var(np.r_[A, B], axis=0)
                Sj = S1[j,i]
                Sk = np.mean(B * (AB[:,k] - A), axis=0) / np.var(np.r_[A, B], axis=0)
                S2[j,k,i]=Vjk - Sj - Sk            
        
        if n_resample > 0:
            r = np.random.randint(N, size=(N, n_resample))

            ST_re = np.zeros([D,n_resample])
            S1_re = np.zeros([D,n_resample])
            S2_re = np.zeros([D,D,n_resample])
            
            
            for j in range(n_resample):
                A_re = A[r[:,j]]
                B_re = B[r[:,j]]
                AB_re = AB[r[:,j],:]
                BA_re = BA[r[:,j],:]
                
                for k in range(D):
                    S1_re[k,j]=np.mean(B_re * (AB_re[:,k] - A_re), axis=0) / np.var(np.r_[A_re, B_re], axis=0)
                    ST_re[k,j]=0.5 * np.mean((A_re - AB_re[:,k]) ** 2, axis=0) / np.var(np.r_[A_re, B_re], axis=0)
                    
                    for l in range(k+1,D):
                        Vjk = np.mean(BA_re[:,k] * AB_re[:,l] - A_re * B_re, axis=0) / np.var(np.r_[A_re, B_re], axis=0)
                        Sj = S1_re[k,j]
                        Sk = np.mean(B_re * (AB_re[:,l] - A_re), axis=0) / np.var(np.r_[A_re, B_re], axis=0)
                        S2_re[k,l,j]=Vjk - Sj - Sk
            for j in range(D):
                S1_conf[j,i]=S1_re[j,:].std(ddof=1)
                ST_conf[j,i]=ST_re[j,:].std(ddof=1)
                for k in range(j,D):                    
                    S2_conf[j,k,i]=S2_re[j,k,:].std(ddof=1)
    filename = "%s_%i_%i" % (output_filename,N,start_ind)
    np.savez(filename,S1=S1,S1_conf=S1_conf,S2=S2,S2_conf=S2_conf,ST=ST,ST_conf=ST_conf)