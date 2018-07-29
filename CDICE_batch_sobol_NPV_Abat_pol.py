# -*- coding: utf-8 -*-
import sys
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
sys.path.append('./SALib-master')

from SALib.analyze import sobol
from SALib.util import read_param_file

#Get the number of processors and the rank of processors
rank = comm.rank
nprocs = comm.size

#Determine the chunk which each processor will neeed to do
count = 60/nprocs

#Use the processor rank to determine the chunk of work each processor will do
start = rank*count
stop = start+count

#Initialize variables
S1 = np.zeros((29,stop-start))
S1_conf = np.zeros((29,stop-start))
S2 = np.zeros((29,29,stop-start))
S2_conf = np.zeros((29,29,stop-start))
ST = np.zeros((29,stop-start))
ST_conf = np.zeros((29,stop-start))

#Load in data
problem = read_param_file('SALib-master/SALib_Input.txt',' ')
Y_matrix = np.loadtxt('Abat.txt')
N = np.size(Y_matrix,0)
cnt = 0
for i in range(start,stop):
    #Y = np.zeros((N))
    Y = Y_matrix[:,i]
    Si = sobol.analyze(problem,Y,calc_second_order=True,num_resamples=2,conf_level=0.95,print_to_console=False)#Need to confirm
    S1[:,cnt] = Si.get("S1")
    S1_conf[:,cnt] = Si.get("S1_conf")
    S2[:,:,cnt] = Si.get("S2")
    S2_conf[:,:,cnt] = Si.get("S2_conf")
    ST[:,cnt] = Si.get("ST")
    ST_conf[:,cnt] = Si.get("ST_conf")
    cnt = cnt+1

filename = "Abat_%i.txt" % rank
np.savez(filename,S1=S1,S1_conf=S1_conf,S2=S2,S2_conf=S2_conf,ST=ST,ST_conf=ST_conf)

#Initialize variables
S1 = np.zeros((29,stop-start))
S1_conf = np.zeros((29,stop-start))
S2 = np.zeros((29,29,stop-start))
S2_conf = np.zeros((29,29,stop-start))
ST = np.zeros((29,stop-start))
ST_conf = np.zeros((29,stop-start))

#Load in data
problem = read_param_file('SALib-master/SALib_Input.txt',' ')
Y_matrix = np.loadtxt('NPV_Abat.txt')
N = np.size(Y_matrix,0)
cnt = 0
for i in range(start,stop):
    #Y = np.zeros((N))
    Y = Y_matrix[:,i]
    Si = sobol.analyze(problem,Y,calc_second_order=True,num_resamples=2,conf_level=0.95,print_to_console=False)#Need to confirm
    S1[:,cnt] = Si.get("S1")
    S1_conf[:,cnt] = Si.get("S1_conf")
    S2[:,:,cnt] = Si.get("S2")
    S2_conf[:,:,cnt] = Si.get("S2_conf")
    ST[:,cnt] = Si.get("ST")
    ST_conf[:,cnt] = Si.get("ST_conf")
    cnt = cnt+1

filename = "NPV_Abat_%i.txt" % rank
np.savez(filename,S1=S1,S1_conf=S1_conf,S2=S2,S2_conf=S2_conf,ST=ST,ST_conf=ST_conf)