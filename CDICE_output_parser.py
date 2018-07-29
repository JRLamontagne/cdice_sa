# This code reads the CDICE output, collapses it to a vector of output for a 
# single variable then saves it to a txt file for SALib.

def parsers(nProcs): 
    print 'nProcs %i' % nProcs    
    print 'loading np...'
    import numpy as np
    a = nProcs-1
    filename = "CDICE_output_%i.txt" % a
    hold = np.loadtxt(filename)
    N = int(hold[np.size(hold,0)-1,0])
    T = 60
# util
    print 'running util...'    
    hold2 = np.zeros((N+1,T))
    
    for i in range(0,nProcs):
        filename = "CDICE_output_%i.txt" % i
        hold = np.loadtxt(filename)
        
        for j in range(0,np.size(hold,0)):
            hold2[hold[j,0],hold[j,1]] = hold[j,2]
			
            
    np.savetxt('util.txt',hold2,delimiter=' ')
# pop    
    print 'running pop...'    
    hold2 = np.zeros((N+1,T))
    
    for i in range(0,nProcs):
        filename = "CDICE_output_%i.txt" % i
        hold = np.loadtxt(filename)
        
        for j in range(0,np.size(hold,0)):
            hold2[hold[j,0],hold[j,1]] = hold[j,3]

    np.savetxt('pop.txt',hold2,delimiter=' ')
# dis_fact    
    print 'running dis_fact...'        
    hold2 = np.zeros((N+1,T))
    
    for i in range(0,nProcs):
        filename = "CDICE_output_%i.txt" % i
        hold = np.loadtxt(filename)
        
        for j in range(0,np.size(hold,0)):
            hold2[hold[j,0],hold[j,1]] = hold[j,4]			
            
    np.savetxt('dis_fact.txt',hold2,delimiter=' ')
# NPV_Dam    
    print 'running NPV_Dam...'        
    hold2 = np.zeros((N+1,T))
    
    for i in range(0,nProcs):
        filename = "CDICE_output_%i.txt" % i
        hold = np.loadtxt(filename)
        
        for j in range(0,np.size(hold,0)):
            hold2[hold[j,0],hold[j,1]] = hold[j,5]
            
    np.savetxt('NPV_Dam.txt',hold2,delimiter=' ')
 # NPV_Abat   
    print 'running NPV_Abat...'        
    hold2 = np.zeros((N+1,T))
    
    for i in range(0,nProcs):
        filename = "CDICE_output_%i.txt" % i
        hold = np.loadtxt(filename)
        
        for j in range(0,np.size(hold,0)):
            hold2[hold[j,0],hold[j,1]] = hold[j,6]

    np.savetxt('NPV_Abat.txt',hold2,delimiter=' ')
# SCC
    print 'running SCC...'        
    hold2 = np.zeros((N+1,T))
    
    for i in range(0,nProcs):
        filename = "CDICE_output_%i.txt" % i
        hold = np.loadtxt(filename)
        
        for j in range(0,np.size(hold,0)):
            hold2[hold[j,0],hold[j,1]] = hold[j,7]			
            
    np.savetxt('SCC.txt',hold2,delimiter=' ')
# Tatm
    print 'running Tatm...'    
    hold2 = np.zeros((N+1,T))
    
    for i in range(0,nProcs):
        filename = "CDICE_output_%i.txt" % i
        hold = np.loadtxt(filename)
        
        for j in range(0,np.size(hold,0)):
            hold2[hold[j,0],hold[j,1]] = hold[j,8]			
            
    np.savetxt('Tatm.txt',hold2,delimiter=' ')
# mat
    print 'running Mat...'    
    hold2 = np.zeros((N+1,T))
    
    for i in range(0,nProcs):
        filename = "CDICE_output_%i.txt" % i
        hold = np.loadtxt(filename)
        
        for j in range(0,np.size(hold,0)):
            hold2[hold[j,0],hold[j,1]] = hold[j,9]			
            
    np.savetxt('Mat.txt',hold2,delimiter=' ')
# GWP
    print 'running GWP...'    
    hold2 = np.zeros((N+1,T))
    
    for i in range(0,nProcs):
        filename = "CDICE_output_%i.txt" % i
        hold = np.loadtxt(filename)
        
        for j in range(0,np.size(hold,0)):
            hold2[hold[j,0],hold[j,1]] = hold[j,10]
			
    np.savetxt('GWP.txt',hold2,delimiter=' ')
# Monetary Discount Rate
    print 'running ri...'    
    hold2 = np.zeros((N+1,T))
    
    for i in range(0,nProcs):
        filename = "CDICE_output_%i.txt" % i
        hold = np.loadtxt(filename)
        
        for j in range(0,np.size(hold,0)):
            hold2[hold[j,0],hold[j,1]] = hold[j,11]
			
    np.savetxt('RI.txt',hold2,delimiter=' ')
# caa
    print 'running caa...'    
    hold2 = np.zeros((N+1,T))
    
    for i in range(0,nProcs):
        filename = "CDICE_output_%i.txt" % i
        hold = np.loadtxt(filename)
        
        for j in range(0,np.size(hold,0)):
            hold2[hold[j,0],hold[j,1]] = hold[j,12]
			
    np.savetxt('caa.txt',hold2,delimiter=' ')
# caa_up_tstep
    print 'running caa_up_tstep...'    
    hold2 = np.zeros((N+1,T))
    
    for i in range(0,nProcs):
        filename = "CDICE_output_%i.txt" % i
        hold = np.loadtxt(filename)
        
        for j in range(0,np.size(hold,0)):
            hold2[hold[j,0],hold[j,1]] = hold[j,13]
			
    np.savetxt('caa_up_tstep.txt',hold2,delimiter=' ')
    return

def make_scalar():
    print 'loading np...'
    import numpy as np
    # abat
    print 'Running Abat...'
    filename = "NPV_Abat.txt"
    hold = np.loadtxt(filename)
    ri = np.loadtxt('Mon_Disc.txt')
    N = np.size(hold,0)
    T = 60
    abat = np.zeros([N,T])
    print 'N...%i' % N
    print 'T...%i' % T
    for i in range(0,N):
        abat[i,0] = hold[i,0]
        for j in range(1,T):
            ri_disc = 1/((1+ri[i,j])**(5*j))
            abat[i,j] = (1/ri_disc)*(hold[i,j]-hold[i,j-1])
    np.savetxt('abat.txt',abat,delimiter=' ')
    abat_scalar = abat[:,58]
    np.savetxt('NPV_Abat_scalar.txt',abat_scalar,delimiter=' ')
    
    # dam
    print 'Running Dam...'
    filename = "NPV_Dam.txt"
    hold = np.loadtxt(filename)
    N = np.size(hold,0)
    T = 60
    dam = np.zeros([N,T])
    for i in range(0,N):
        dam[i,0] = hold[i,0]
        for j in range(1,T):
            ri_disc = 1/((1+ri[i,j])**(5*j))
            dam[i,j] = (1/ri_disc)*(hold[i,j]-hold[i,j-1])
    np.savetxt('dam.txt',dam,delimiter=' ')
    dam_scalar = dam[:,58]
    np.savetxt('NPV_Dam_scalar.txt',dam_scalar,delimiter=' ')