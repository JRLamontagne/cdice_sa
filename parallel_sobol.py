# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 09:30:38 2017

@author: jrl276
"""
def parallel_sobol(input_filename,n_resample,N,D,output_filename):
    import sys
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    sys.path.append('./SALib-master')

    from custom_sobol import custom_sobol

    #Get the number of processors and the rank of processors
    rank = comm.rank
    nprocs = comm.size

    #Determine the chunk which each processor will neeed to do
    count = 64/nprocs

    #Use the processor rank to determine the chunk of work each processor will do
    start = rank*count
    stop = start+count
    
    for i in range(start,stop):
        print i
        custom_sobol(input_filename,n_resample,i,i,N,D,output_filename)