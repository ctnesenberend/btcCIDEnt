import random
import numpy as np 
import matplotlib.pyplot as plt
import math
import gzip
import zlib 
import struct
import lz4.frame
import lzma


def YakovenkoSim (N = 500, NSweeps=[1,10,50,200,500,1000], moneyint = 10, DeltaMax = 1, binsize = None, discrete = False, seed = 2222):
    # The simple simulation used for the data
    
    random.seed(seed)
    if binsize == None:
        binsize = DeltaMax
    money = [moneyint]*N
    print(f"Doing {NSweeps[0]} sweeps for plot 0")
    YakovenkoIters(NSweeps[0]*N, money, DeltaMax, discrete)
    for i in range(0,len(NSweeps)):
        YakovenkoPlot(money,N,NSweeps[i],binsize,moneyint,DeltaMax,i)
        if i < len(NSweeps) - 1:
            print(f"Doing {NSweeps[i+1]-NSweeps[i]} sweeps for plot {i + 1}")
            YakovenkoIters((NSweeps[i+1] - NSweeps[i])*N, money,DeltaMax, discrete)
    plt.show()
    
def CIDOverTime(N = 500, NSweeps = 10000, moneyint = 10, DeltaMax = 1, binsize = None, stepsize = None,
zliblevel = -1, gziplevel = 9, discrete = False, seed = 2222):
    # A script that plots the CID over time of our simulation
    
    random.seed(seed)
    if binsize == None:
        binsize = moneyint/10
    if stepsize == None:
        stepsize = math.floor(NSweeps/1000)
    money = [moneyint]*N
    zlibCIDs = [0]*(math.ceil(NSweeps/stepsize)+1)
    gzipCIDs = [0]*(math.ceil(NSweeps/stepsize)+1)
    lz4CIDs = [0]*(math.ceil(NSweeps/stepsize)+1)
    lzmaCIDs = [0]*(math.ceil(NSweeps/stepsize)+1)
    nsweeparray = range(0,NSweeps+1,stepsize)
    sweepstotal = 0
    istep = 1
    print("Doing sweeps", end = "")
    while(sweepstotal < NSweeps):
        print("\r                                                                         ", end = "")
        print(f"\rDoing sweeps {sweepstotal} to {sweepstotal + stepsize} out of {NSweeps}\r", end = "")
        YakovenkoIters(stepsize*N,money,DeltaMax, discrete)
        zlibCIDs[istep], gzipCIDs[istep], lz4CIDs[istep], lzmaCIDs[istep] = CIDs(money,zliblevel,gziplevel)
        sweepstotal += stepsize
        istep += 1
    plt.figure(1)
    #plt.subplot(1,2,1)
    #YakovenkoPlot(money,N,sweepstotal,binsize,moneyint,DeltaMax)
    #plt.subplot(1,2,2)
    plt.plot(nsweeparray,zlibCIDs,label="zlib")
    plt.plot(nsweeparray,gzipCIDs,label="gzip")
    plt.plot(nsweeparray,lz4CIDs,label="lz4")
    plt.plot(nsweeparray,lzmaCIDs,label="lzma")
    plt.title("CID over time of the model")
    plt.xlabel(f"Number of Yakovenko iterations/{N}")
    plt.ylabel("CID")
    plt.legend()
    plt.show()
        
def CIDOverM0(N = 500, DeltaMax = 1, stepsize = 10, minsteps = 50, CIDConstraint = 0.02, zliblevel = -1, gziplevel = 9, discrete = False, seed = 2222, m0scale = 1,maxsteps= 100000):
    # A script that plots the CID over m0 of our simulation
    
    random.seed(seed)
    basic = list(range(1,5,1)) + list(range(5,25,5)) + list(range(25,100,25)) + list(range(100,450,50))
    m0s = [i* m0scale for i in basic]
    zlibCIDm0s = [0]*len(m0s)
    gzipCIDm0s = [0]*len(m0s)
    lz4CIDm0s = [0]*len(m0s)
    lzmaCIDm0s = [0]*len(m0s)
    shannonentr = [0]*len(m0s)
    for i in range(len(m0s)):
        zlibCIDm0s[i], gzipCIDm0s[i], lz4CIDm0s[i], lzmaCIDm0s[i], shannonentr[i] = CIDGivenM0(N, m0s[i], DeltaMax, stepsize, minsteps, CIDConstraint, zliblevel, gziplevel, discrete,maxsteps)
        lz4CIDm0s[i] = lz4CIDm0s[i]*.6
        lzmaCIDm0s[i] = lzmaCIDm0s[i]*1.15
    plt.figure(1)
    plt.plot(m0s,zlibCIDm0s, label = "zlib")
    plt.plot(m0s,gzipCIDm0s, label = "gzip")
    plt.plot(m0s,lz4CIDm0s, label = ".6*lz4")
    plt.plot(m0s,lzmaCIDm0s, label = "1.15*lzma")
    plt.plot(m0s,shannonentr, label = "Shannon Entropy/17")
    plt.title(f"CID vs $m_0$, $\Delta_{{max}} =$ {DeltaMax}")
    plt.xlabel("$m_0$")
    plt.ylabel("CID")
    plt.legend()
    plt.figure(2)
    zlibdiv = [x/y for x,y in zip(zlibCIDm0s,shannonentr)]
    gzipdiv = [x/y for x,y in zip(gzipCIDm0s,shannonentr)]
    lz4div = [x/y for x,y in zip(lz4CIDm0s,shannonentr)]
    lzmadiv = [x/y for x,y in zip(lzmaCIDm0s,shannonentr)]
    plt.plot(m0s,zlibdiv, label = "zlib/Shannon Entropy")
    plt.plot(m0s,gzipdiv, label = "gzip/Shannon Entropy")
    plt.plot(m0s,lz4div, label = ".6*lz4/Shannon Entropy")
    plt.plot(m0s,lzmadiv, label = "1.15*lzma/Shannon Entropy")
    plt.legend()
    plt.show()
    
def CIDGivenM0(N = 500, moneyint = 10, DeltaMax = 1, stepsize = 10, minsteps = 50, CIDConstraint = 0.02, 
zliblevel = -1, gziplevel = 9, discrete = False,maxsteps = 100000):
    # A script that calculates the CID for particular parameters (called CIDGivenM0 because I initially used it for the m0 graphs)
    
    money = [moneyint]*N
    # zlibCIDs = False
    # gzipCIDs = False
    # nsweeps = 0
    YakovenkoIters((maxsteps-stepsize*minsteps)*N,money,DeltaMax,discrete)
    zlibCIDs = [0]*(minsteps)
    gzipCIDs = [0]*(minsteps)
    lz4CIDs = [0]*(minsteps)
    lzmaCIDs = [0]*(minsteps)
    shannonentr = [0]*(minsteps)
    istep = 0
    while(istep < minsteps):
        YakovenkoIters(stepsize*N,money,DeltaMax, discrete)
        zlibCIDs[istep], gzipCIDs[istep], lz4CIDs[istep], lzmaCIDs[istep] = CIDs(money,zliblevel,gziplevel, discrete)
        shannonentr[istep] = ShannonEntropy(money)
        istep += 1  
    # while(nsweeps < 100000):#CIDsNotConverged(zlibCIDs, gzipCIDs, CIDConstraint)):
        # zlibCIDs = [0]*(minsteps)
        # gzipCIDs = [0]*(minsteps)
        # shannonentr = [0]*(minsteps)
        # istep = 0
        # while(istep < minsteps):
            # zlibCIDs[istep], gzipCIDs[istep] = CIDs(money,zliblevel,gziplevel, discrete)
            # shannonentr[istep] = ShannonEntropy(money)
            # YakovenkoIters(stepsize*N,money,DeltaMax, discrete)
            # istep += 1
            # nsweeps += stepsize
    print(f"Got CID for m0 = {moneyint} after {maxsteps} sweeps of size {N}")
    return np.mean(zlibCIDs), np.mean(gzipCIDs), np.mean(lz4CIDs), np.mean(lzmaCIDs), np.mean(shannonentr)/17

def CIDOverN(moneyint = 10, DeltaMax = 1, stepsize = 10, minsteps = 50, CIDConstraint = 0.02, zliblevel = -1, gziplevel = 9, discrete = False, seed = 2222, m0scale = 1):
    # a script that plots the CID over different N
    
    random.seed(seed)
    basic = list(range(10,50,10)) + list(range(50,250,50)) + list(range(250,1000,250)) + list(range(1000,5000,1000))
    Ns = [i* m0scale for i in basic]
    zlibCIDNs = [0]*len(Ns)
    gzipCIDNs = [0]*len(Ns)
    lz4CIDNs = [0]*len(Ns)
    lzmaCIDNs = [0]*len(Ns)
    shannonentr = [0]*len(Ns)
    for i in range(len(Ns)):
        zlibCIDNs[i], gzipCIDNs[i], lz4CIDNs[i], lzmaCIDNs[i], shannonentr[i] = CIDGivenM0(Ns[i], moneyint, DeltaMax, stepsize, minsteps, CIDConstraint, zliblevel, gziplevel, discrete, 5500)
    plt.figure(1)
    plt.plot(Ns,zlibCIDNs, label = "zlib")
    plt.plot(Ns,gzipCIDNs, label = "gzip")
    plt.plot(Ns,lz4CIDNs, label = "lz4")
    plt.plot(Ns,lzmaCIDNs, label = "lzma")
    plt.plot(Ns,shannonentr, label = "Shannon Entropy/17")
    plt.title(f"CID vs $N$, $\Delta_{{max}} =$ {DeltaMax}, $m_{{0}} =$ {moneyint}")
    plt.xlabel("$N$")
    plt.ylabel("CID")
    plt.legend()
    plt.figure(2)
    zlibdiv = [x/y for x,y in zip(zlibCIDNs,shannonentr)]
    gzipdiv = [x/y for x,y in zip(gzipCIDNs,shannonentr)]
    lz4div = [x/y for x,y in zip(lz4CIDNs,shannonentr)]
    lzmadiv = [x/y for x,y in zip(lzmaCIDNs,shannonentr)]
    plt.plot(Ns,zlibdiv, label = "zlib/Shannon Entropy")
    plt.plot(Ns,gzipdiv, label = "gzip/Shannon Entropy")
    plt.plot(Ns,lz4div, label = "lz4/Shannon Entropy")
    plt.plot(Ns,lzmadiv, label = "lzma/Shannon Entropy")
    plt.legend()
    plt.show()    
    
def CIDOverTimeMean(N = 500, NSweeps = 100000, moneyint = 10, DeltaMax = 1, binsize = None, stepsize = None, stepstepsize = None,
zliblevel = -1, gziplevel = 9, discrete = False, seed = 2222):
    # A script which plots the CID over different N but takes the average over a few steps to get less noise
    
    random.seed(seed)
    if binsize == None:
        binsize = moneyint/10
    if stepsize == None:
        stepsize = math.floor(NSweeps/40)
    if stepstepsize == None:
        stepstepsize = math.floor(stepsize/50)
    money = [moneyint]*N
    zlibCIDs = [0]*(math.ceil(NSweeps/stepsize)+1)
    gzipCIDs = [0]*(math.ceil(NSweeps/stepsize)+1)
    lz4CIDs = [0]*(math.ceil(NSweeps/stepsize)+1)
    lzmaCIDs = [0]*(math.ceil(NSweeps/stepsize)+1)
    nsweeparray = range(0,NSweeps+1,stepsize)
    sweepstotal = 0
    istep = 1
    print("Doing sweeps", end = "")
    while(sweepstotal < NSweeps):
        print("\r                                                                         ", end = "")
        print(f"\rDoing sweeps {sweepstotal} to {sweepstotal + stepsize} out of {NSweeps}\r", end = "")
        zlibCIDs[istep], gzipCIDs[istep], lz4CIDs[istep], lzmaCIDs[istep] = CIDSweepMean(money,DeltaMax,math.floor(stepsize/stepstepsize),stepstepsize)
        sweepstotal += stepsize
        istep += 1
    plt.figure(1)
    plt.subplot(1,2,1)
    YakovenkoPlot(money,N,sweepstotal,binsize,moneyint,DeltaMax)
    plt.subplot(1,2,2)
    plt.plot(nsweeparray,zlibCIDs,label="zlib")
    plt.plot(nsweeparray,gzipCIDs,label="gzip")
    plt.plot(nsweeparray,lz4CIDs,label="lz4")
    plt.plot(nsweeparray,lzmaCIDs,label="lzma")
    plt.title("CID over time of the model")
    plt.xlabel(f"Number of Yakovenko iterations/{N}")
    plt.ylabel(f"CID averaged over {math.floor(stepsize/stepstepsize)} steps of size {stepstepsize}")
    plt.legend()
    plt.show()
    
def CIDSweepMean(money, DeltaMax, steps, stepsize, zliblevel = -1, gziplevel = 9, discrete = False):
    # Applies steps * stepsize sweeps of size len(money) and calculates the average CID of all steps
    
    N = len(money)
    zlibCIDs = [0]*steps
    gzipCIDs = [0]*steps
    lz4CIDs = [0]*steps
    lzmaCIDs = [0]*steps
    istep  = 0
    while(istep < steps):
        YakovenkoIters(stepsize*N, money, DeltaMax, discrete)
        zlibCIDs[istep], gzipCIDs[istep], lz4CIDs[istep], lzmaCIDs[istep] = CIDs(money, zliblevel, gziplevel, discrete)
        istep += 1
    return sum(zlibCIDs)/steps, sum(gzipCIDs)/steps, sum(lz4CIDs)/steps, sum(lzmaCIDs)/steps
    

def CIDsNotConverged(zlibCIDs, gzipCIDs, lz4CIDs, lzmaCIDs, CIDConstraint = 0.02):
    # an attempt at trying a constraint on the CID being converged
    
    if zlibCIDs == False or gzipCIDs == False or lz4CIDs == False or lzmaCIDs == False:
        return True
    return ((max(zlibCIDs) - min(zlibCIDs)) > CIDConstraint) or ((max(gzipCIDs) - min(gzipCIDs)) > CIDConstraint) or ((max(lz4CIDs) - min(lz4CIDs)) > CIDConstraint) or ((max(lzmaCIDs) - min(lzmaCIDs)) > CIDConstraint)
    
def CIDs (money, zliblevel = -1, gziplevel = 9, discrete = False):
    # the actual CID calculation happens here
    
    if discrete:
        buf = struct.pack(f'{len(money)}I', *money)
    else:
        money2 = [round(x) for x in money]
        buf = struct.pack(f'{len(money)}I', *money2) # needs to be int for a nice compression
        # if I really want to be exact I can do this instead of rounding, I don't think it's worthwhile though:
        # moneynumerator = [0]*len(money)
        # moneydivisor = [0]*len(money)
        # for i in len(money)
            # moneynumerator[i], moneydivisor[i] = money[i].as_integer_ratio()
        # div = max(moneydivisor)
        # moneydiscr = [(div * moneynumerator[i])//moneydivisor[i] for i in len(money)]
    zlibCIDv = len(zlib.compress(buf,zliblevel))/len(buf)
    gzipCIDv = len(gzip.compress(buf,gziplevel))/len(buf)
    lz4CIDv = len(lz4.frame.compress(buf))/len(buf)
    lzmaCIDv = len(lzma.compress(buf))/len(buf)
    return zlibCIDv, gzipCIDv, lz4CIDv, lzmaCIDv

            
            
def zlibCID (bbytes, level = -1):
    return len(zlib.compress(bbytes,level))/len(bbytes)
    
def gzipCID (bbytes,compresslevel = 9):
    return len(gzip.compress(bbytes,compresslevel))/len(bbytes)
    
def lz4CID (bbytes):
    return len(lz4.frame.compress(bbytes))/len(bbytes)
                
def lzmaCID (bbytes):
    return len(lzma.compress(bbytes))/len(bbytes)


def BoltzmannHist(x,binsize,moneyint,N):
    return N*(1 - np.exp(-binsize/moneyint))*np.exp(-(x - binsize/2)/moneyint)

def BoltzmannHistPlot(moneyint,maxmoney,N,binsize,figure=None):
    if figure != None:  
        plt.figure(figure)
    x = np.arange(0.0,maxmoney+binsize,binsize/4)
    y = BoltzmannHist(x,binsize,moneyint,N)
    plt.plot(x,y)
    
    

def YakovenkoPlot(money,N,NSweepstotal,binsize,moneyint,DeltaMax,figure=None,boltzmanndistr=True):
    #A script that plots the basic Yakovenko distribution
    
    if figure != None:  
        plt.figure(figure)
    plt.hist(money, bins=np.arange(0,max(money)+ binsize,binsize))
    plt.title(f"Histogram of N = {N} agents after {NSweepstotal}$\cdot N$ Yakovenko iterations")
    plt.xlabel(f"Money ($m_0 = ${moneyint}, $\Delta_{{max}} = ${DeltaMax})")
    plt.ylabel("Amount of agents")    
    if boltzmanndistr == True:
        BoltzmannHistPlot(moneyint,max(money),N,binsize)
    
    
    
def YakovenkoIter(money, DeltaMax, discrete = False):
    #A script that does a single Yakovenko itration on money
    
    N = len(money)
    i = random.randrange(N)
    j = random.randrange(N)
    while i == j:
        j = random.randrange(N)
    if discrete and DeltaMax == 1:
        Delta = 1
    elif discrete:
        Delta = random.randint(1,DeltaMax)
    else:
        Delta = random.uniform(0,DeltaMax)
    if money[i] >= Delta:
        money[i] = money[i] - Delta
        money[j] = money[j] + Delta
        

def YakovenkoIters(NIters, money, DeltaMax, discrete = False):
    #A script that does multiple yakovenko iterations on money 
    
    for i in range(NIters):
        YakovenkoIter(money, DeltaMax, discrete)


def ShannonEntropy(money):
    #A simple shannon entropy calculation 
    
    N = math.ceil(max(money))
    S = 0
    for i in range(N):
        l = 0
        for x in money:
            if i <= x and x < i + 1:
                l += 1
        if l != 0:
            p = l / len(money)
            S += -p * math.log2(p)
    return S

# def ShannonEntropy (money, bins = 25):
    # hist, edges = np.histogram(money, bins)
    # hist = hist/hist.sum()
    # ent = 0
    # for bin in hist:
        # ent += -bin*math.log2(bin)


    



# def sweepsim (N = 1000, Nsweeps = 10, moneyint = 100, DeltaMax = 1):
    # money = [moneyint]*N
    # sweeps (money, Nsweeps, DeltaMax)
    # plt.hist(money, bins=25)
    # plt.title("Histogram of N = "+str(N)+" agents after "+str(Nsweeps)+" sweeps")
    # plt.xlabel("Money ($m_i = $"+str(moneyint)+", $\Delta_{max} = $"+str(DeltaMax)+")")
    # plt.ylabel("Amount of agents")
    # plt.show()
    
    
    
    

# def sweep (money, DeltaMax):
    # N = len(money)
    # givers = list(range(N))
    # random.shuffle(givers)
    ## there's a possibility of money being exchanged from and to the same agent
    ## the thesis is unclear on whether this is supposed to happen
    # for i in range(N):
        # Delta = random.uniform(0,DeltaMax)
        # if money[givers[i]] >= Delta:
            # money[givers[i]] = money[givers[i]] - Delta
            # money[i] = money[i] + Delta
    
    
    

# def sweeps (money, Nsweeps, DeltaMax):
    # for n in range(Nsweeps):
        # sweep (money, DeltaMax)
        
