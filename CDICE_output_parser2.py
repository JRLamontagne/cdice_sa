# This code reads the CDICE output, collapses it to a vector of output for a 
# single variable then saves it to a txt file for SALib.

def parsers(nProcs): 
    print 'nProcs %i' % nProcs    
    print 'loading np...'
    import numpy as np
    from mpi4py import MPI
    
    comm = MPI.COMM_WORLD
    
    var_name =['util','pop','util_disc','real_intr','NPV_Dam','NPV_Abat','Abat','Dam','GWP','SCC','Tatm','Mat','forc','caa','caa_up_tstep','ygross']
    file_name =['util.txt','pop.txt','util_disc.txt','real_intr.txt','NPV_Dam.txt','NPV_Abat.txt','Abat.txt','Dam.txt','GWP.txt','SCC.txt','Tatm.txt','Mat.txt','forc.txt','caa.txt','caa_up_tstep.txt','ygross.txt']

    #Get the number of processors and the rank of processors
    rank = comm.rank
    nProcs = comm.size
    
    a = nProcs-1
    filename = "CDICE_output_%03d.txt" % a
    hold = np.loadtxt(filename)
    N = int(hold[np.size(hold,0)-1,0])
    T = 100
    
    if rank <= 15:
	print 'proc %i' % rank
        print 'running %s...' % var_name[rank]
        hold2 = np.zeros((N+1,T))
    
        for i in range(0,nProcs):
            filename = "CDICE_output_%03d.txt" % i
            hold = np.loadtxt(filename)
        
            for j in range(0,np.size(hold,0)):
                hold2[hold[j,0],hold[j,1]] = hold[j,rank+2]
			
            
        np.savetxt(file_name[rank],hold2,delimiter=' ')
    return

def make_scalar():
    print 'loading np...'
    import numpy as np
    # abat
    print 'Running Abat...'
    filename = "NPV_Abat.txt"
    hold = np.loadtxt(filename)
    ri = np.loadtxt("Mon_Disc.txt")
    N = int(hold[np.size(hold,0)-1,0])
    T = 60
    abat = np.zeros([N,T])
    for i in range(0,N):
        abat[i,0] = hold[i,0]
        for j in range(1,T):
            ri_disc = 1/((1+ri[i,j])**(5*j))
            abat[i,j] = (1/ri_disc)*(hold[i,j]-hold[i,j-1])
    np.savetxt('NPV_Abat.txt',abat,delimiter=' ')
    abat_scalar = abat[:,58]
    np.savetxt('NPV_Abat_scalar.txt',abat_scalar,delimiter=' ')
    
    # dam
    print 'Running Dam...'
    filename = "NPV_Dam.txt"
    hold = np.loadtxt(filename)
    N = int(hold[np.size(hold,0)-1,0])
    T = 60
    dam = np.zeros([N,T])
    for i in range(0,N):
        dam[i,0] = hold[i,0]
        for j in range(1,T):
            ri_disc = 1/((1+ri[i,j])**(5*j))
            dam[i,j] = (1/ri_disc)*(hold[i,j]-hold[i,j-1])
    np.savetxt('NPV_Dam.txt',dam,delimiter=' ')
    dam_scalar = dam[:,58]
    np.savetxt('NPV_Dam_scalar.txt',dam_scalar,delimiter=' ')