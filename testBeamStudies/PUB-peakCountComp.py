import ROOT
from ROOT import TCanvas, TH1D, TGraph, TFile, TTree, TH2D, TLegend
import pickle
import numpy as np
from array import array
import math
import glob

""" Dylan Mead's Peak Counting algorithm rewritten in Python by Jason Barkeloo-Univ of Oregon"""
sim=True
test=False #W-first
revTest=True #Si-first
######################################################
def listoflist(length):
    a = []
    for i in range(length):
        a.append([])
    return a

def distance(x1,y1,x2,y2):
    return math.sqrt((x2-x1)**2+(y2-y1)**2)

#######################################################


########################################################
##  Creation of Pixel X- and Y-Positions  ##############
########################################################
if test or revTest:
    inFile = '/home/jason/jas-assembly-3.1.2/extensions/BarkelooCode/resources/geom+centers2.txt'
    f = open(inFile,'r')
    pos = np.zeros((1024,2))
    xTest, yTest, n  = 0., 0. ,0
    for i in range(1024):
        for line in f:
            if float(line.split()[0]) == i:
                xTest += float(line.split()[1])
                yTest += float(line.split()[2])
                n+=1
            else:
                break
        pos[i][0] = xTest/n
        pos[i][1] = yTest/n
        xTest,yTest,n =0,0,0
    del xTest, yTest, n
    print "X and Y positions to Pixels Calculated"

########################################################
###########   FOR SIMULATION DATA
########################################################
if sim:
    inFile = '/home/jason/data/amanda/smear9FinalOuts/DENS23-SiWNi/withPos/shifted/12.1GeV-smear9-SiWNi(0.8725)_shifted23.endataIMP.pkl'
    inFile1 = '/home/jason/data/amanda/smear9FinalOuts/DENS23-SiWNi/withPos/shifted/12.1GeV-smear9-SiWNi(0.8725)_shifted23_elinfoIMP.pkl'
    unpicklefile = open(inFile, 'r')
    endataSim = pickle.load(unpicklefile)
    unpicklefile = open(inFile1,'r')
    elinfo = pickle.load(unpicklefile)
    NmbEventsSim=len(endataSim)
    print NmbEventsSim
    #positionArray1=pickle.load(open('/home/jason/data/amanda/smear9FinalOuts/12.1GeV-smear9-SiWNi(0.8725)_1shifted23.positions.pkl','rb'))
    positionArray2=pickle.load(open('/home/jason/data/amanda/smear9FinalOuts/DENS23-SiWNi/withPos/shifted/12.1GeV-smear9-SiWNi(0.8725)_2shifted23.positions.pkl','rb'))
    #positionArray3=pickle.load(open('/home/jason/data/amanda/smear9FinalOuts/12.1GeV-smear9-SiWNi(0.8725)_3shifted23.positions.pkl','rb'))
    #positionArray4=pickle.load(open('/home/jason/data/amanda/smear9FinalOuts/12.1GeV-smear9-SiWNi(0.8725)_4shifted23.positions.pkl','rb'))
    #positionArray5=pickle.load(open('/home/jason/data/amanda/smear9FinalOuts/12.1GeV-smear9-SiWNi(0.8725)_5shifted23.positions.pkl','rb'))
    multEArray=[8000,3490,1015,221,38] #SiWNi
#    multEArray=[8000,3092,796,153,23]#WNiSi
    
    print inFile
    print inFile1
########################################################
###########   FOR TEST BEAM DATA
########################################################
if test or revTest:
    inFile = '/home/jason/data/TestBeam/TestBeamROOT/2013_07_26_14_13_43.bin.endata.pkl'
    #inFile1 = '/home/jason/data/TestBeam/TestBeamROOT/2013_07_26_14_13_43.bin.NmbEvents.pkl'    
    unpicklefile = open(inFile, 'r')
    endataTest = pickle.load(unpicklefile)
    NmbEventsTest=len(endataTest)
    sep2=[]

print inFile
#print inFile1
########################################################
###########   FOR DEALING WITH NEIGHBORS
########################################################
#lets go through a file and find which pixels neighbor others
neighborFile = '/home/jason/jas-assembly-3.1.2/extensions/BarkelooCode/resources/neighborsNewest.txt'
##neighbors = []
##for i in range(1024):
##    neighbors.append([])
neighbors = listoflist(1024)

testpnt = 0
for i in open(neighborFile,'r'):
    neighbors[testpnt]=i.split()
    testpnt+=1
##  changing neighbors to int as it should be.
for i in range(len(neighbors)):
	for j in range(len(neighbors[i])):
		neighbors[i][j]=int(neighbors[i][j])

##### neighbors is now a list with index of a kpixChannel and each entry is a list of its' neighboring cells
# Now some very basic memory management
del testpnt
del neighborFile

def flatten(l):
  out = []
  for item in l:
    if isinstance(item, (list, tuple)):
      out.extend(flatten(item))
    else:
      out.append(item)
  return out


##### clean test beam
if test or revTest:
    cleanedEvents=[]
    for i in range(NmbEventsTest):
        if i%100==0:
            print "Test Beam Event: "+str(i)+" of "+str(NmbEventsTest)
        totalEnergy=endataTest[i].sum()*100000000000000
        statNum=0
        statDenom=totalEnergy
        nullCount=0
        for j in range(9):
            for k in range(1024):
                if test:
                    statNum+=endataTest[i][j][k]*(j+1)*(j+1)*100000000000000
                if revTest:
                    statNum+=endataTest[i][j][k]*(9-j)*(9-j)*100000000000000
        for j in range(9):
            if endataTest[i][j].sum()==0:
                nullCount+=1
                if test:
                    statNum+=5.*(j+1)*(j+1)
                if revTest:
                    statNum+=5.*(9-j)*(9-j)
                statDenom+=5.
        if nullCount==8 or statNum/statDenom<44:
            cleanedEvents.append(i) #contaminant events


def peakCount(kpixArray, neighbor, NmbEvents):
#    NmbEvents = len(kpixArray)
    """Peak Counting Method.  Coded by Jason Barkeloo.  Takes Array of NmbEvents x layers(0-8 for current iteration) x kpixChannel filled with amplitude info (adc/energy/etc.)"""
    countedPeaks=np.zeros(NmbEvents)
    countedPeaks1=np.zeros(NmbEvents)
    countedPeaks2=np.zeros(NmbEvents)
    localMaxes =listoflist(NmbEvents)
    arrayMaxes=[]
    for i in range(NmbEvents):
        if i%100==0:
            print "Event: "+str(i)+" of "+str(NmbEvents)
        detector = kpixArray[i] #should be 9x1024 array accordingly
        ##### EVENS and ODDS Layer summed list
        evens = np.zeros(1024) #sum of even layers
        odds = np.zeros(1024) #sum of odd layers
        full = np.zeros(1024) #sum of all layers
        for j in range(9):
            if revTest:
                weight=(0.2+0.1*j)
            else:
                weight = (1-0.1*j)
            full += kpixArray[i][j]*weight #weighting full towards front
            if j%2==0:#if even add to evens list( Layers 0,2,4,6,8)
                evens += kpixArray[i][j]
            else:
                odds += kpixArray[i][j]
        ############ evens and odds are created
        #Lets find some local maxes in each layer!
        arrayMaxes =[[],[],[],[],[],[],[],[],[]] #get maxes for each layer individually 
        evenMaxes = [] #get maxes for the sum of all even layers
        oddMaxes = [] #get maxes for the sum of all odd layers
        fullMaxes =[] #get maxes for the sum of all layers
        threshold =0.1* detector.max()
        for j in range(len(detector)):#j is layer
            for k in range(len(detector[0])):#k is channels
                if detector[j][k]==0:
                    pass
                elif detector[j][k]>0:#threshold=0
                    neighborCheck = 0
                    for l in range(len(neighbor[k])):
                        if detector[j][k] > detector[j][neighbor[k][l]]:
                            continue
                        else:
                            neighborCheck +=1
                    if neighborCheck == 0: #if no neighbors are higher than it, it is candidate
                        arrayMaxes[j].append(k)#appending channels that are local maxes
            ###### Evens and Odd deals start here:
                if evens[k]==0:
                    pass
                elif evens[k] > threshold:
                    neighborCheck = 0
                    for l in range(len(neighbor[k])):
                        if evens[k] > evens[neighbor[k][l]]:
                            continue
                        else:
                            neighborCheck+=1
                    if neighborCheck ==0:
                        evenMaxes.append(k)
                #### Odds
                if odds[k]==0:
                    pass
                elif odds[k] > threshold:
                    neighborCheck = 0
                    for l in range(len(neighbor[k])):
                        if odds[k] >odds[neighbor[k][l]]:
                            continue
                        else:
                            neighborCheck+=1
                    if neighborCheck ==0:
                        oddMaxes.append(k)
                #### full
                if full[k]==0:
                    pass
                elif full[k] > 2*threshold:
                    neighborCheck = 0
                    for l in range(len(neighbor[k])):
                        if odds[k] >full[neighbor[k][l]]:
                            continue
                        else:
                            neighborCheck+=1
                    if neighborCheck ==0:
                        fullMaxes.append(k)

        #### Now using evenMaxes, oddMaxes and arrayMaxes to count the number of peaks
       ####### Dylan's method ######
        testlayer=7
        for j in range(len(arrayMaxes[testlayer])):
            a = 0
            for k in range(len(arrayMaxes)):
                if arrayMaxes[k].count(arrayMaxes[testlayer][j])>0:
                    a+=1 #find how many layers have a max at the same pixel
            if a >2: #require this many layers have a max for an electron
                countedPeaks[i] +=1
       ####### ammended Dylan method (Amanda's method) #####
        arrayMaxesTot=[] 
        eventPeaks=[] #counts the peaks recorded in each event
        prev=[-1] #stores past list elements after they've been analyzed
        count=0
        for m in range(len(arrayMaxes)):
            arrayMaxesTot.extend(arrayMaxes[m]) #combine all arrayMaxes into one large array
        for n in range(len(arrayMaxesTot)):
            if arrayMaxesTot[n] not in prev:
                count=arrayMaxesTot.count(arrayMaxesTot[n])#count how many layers a certain pixel is a local max
                if count>3:
                    countedPeaks2[i]+=1 #get electron hit if at least 4 layers share the same local max
                    eventPeaks.append(arrayMaxesTot[n]) # add that max to the array of event peaks
                count=0
                prev.append(arrayMaxesTot[n]) #put each entry in prev to not double-count
        lenEPeaks=len(eventPeaks)
        for p in range(lenEPeaks):
            for q in range(lenEPeaks-1):
                if p+q<lenEPeaks-1:
                    if eventPeaks[p] in neighbors[eventPeaks[p+q+1]]:#check to see if >1 peaks are neighbors, where an electron hits near the border and deposits evenly between two pixels. If so, remove an electron count
                        countedPeaks2[i]-=1
                        eventPeaks[p]=-1#remove the neighboring pixel - NOT ACCURATE
        if test or revTest:
            if countedPeaks2[i]==2:
                sep2.append(distance(pos[eventPeaks[0]][0],pos[eventPeaks[0]][1],pos[eventPeaks[1]][0],pos[eventPeaks[1]][1]))#distance from center of pixel maxes
       ####### Jason's Method ######
        evenNeighbors, oddNeighbors, ePeaks = [],[],[]
        for j in range(len(evenMaxes)):
            evenNeighbors.append(neighbors[evenMaxes[j]])
        for j in range(len(oddMaxes)):
            oddNeighbors.append(neighbors[oddMaxes[j]])
        evenNeighbors = set(flatten(evenNeighbors))
        oddNeighbors = set(flatten(oddNeighbors))
        for j in set(flatten(arrayMaxes)):
            if j in set(fullMaxes):
                if j in set(evenMaxes) or j in evenNeighbors:
                    if j in set(oddMaxes) or j in oddNeighbors:
                        countedPeaks1[i] +=1 #need max in all 4 array sets to be electron event
                        ePeaks.append(j)
        lenePeaks=len(ePeaks)
        for p in range(lenePeaks):
            for q in range(lenePeaks-1):
                if p+q<lenePeaks-1:
                    if ePeaks[p] in neighbors[ePeaks[p+q+1]]:
                        countedPeaks1[i]-=1

    print "Counted: "+ str(sum(countedPeaks)) +" electrons"
    print "Counted: "+ str(sum(countedPeaks1)) +" electrons"
    print "Counted: "+ str(sum(countedPeaks2)) +" electrons"
    return countedPeaks, countedPeaks1, countedPeaks2

if sim:
    countedPeaksSim,countedPeaks1Sim,countedPeaks2Sim=peakCount(endataSim,neighbors,NmbEventsSim)

if test or revTest:
    countedPeaksTest,countedPeaks1Test,countedPeaks2Test=peakCount(endataTest,neighbors,NmbEventsTest)

##print countedPeaks for sim
truecount=0
truecount1=0
truecount2=0
tag2=0
if sim:
    for i in range(len(elinfo)):#update??
        if countedPeaksSim[i]==elinfo[i]:
            truecount+=1.
        if countedPeaks1Sim[i]==elinfo[i]:
            truecount1+=1.
        if countedPeaks2Sim[i]==elinfo[i]:
            truecount2+=1.
        if countedPeaks2Sim[i]==2:
            tag2+=1
    print "Counted "+str(truecount)+" of "+str(len(elinfo))+" events correctly"
    print "Counted "+str(truecount1)+" of "+str(len(elinfo))+" events correctly"    
    print "Counted "+str(truecount2)+" of "+str(len(elinfo))+" events correctly"
    print "2e tagged events = "+str(tag2)

########################################################
###### For Nick - events with 1.5x mean single electron energy and electron number distribution info
########################################################
"""
totalE=0

for i in range(NmbEvents):
    for j in range(9):
        for k in range(1024):
            totalE+=endata[i][j][k]
    if 70<totalE<75:
        print i, elinfo[i], countedPeaks[i], countedPeaks2[i], countedPeaks1[i]
    totalE=0
    i+=1

oneE=[]
twoE=[]
threeE=[]
fourE=[]
fiveE=[]
oneECount=[]
twoECount=[]
threeECount=[]
fourECount=[]
fiveECount=[]

for i in range(NmbEvents):
    if elinfo[i]==1:
        oneE.append(i)
        oneECount.append(countedPeaks2[i])
    elif elinfo[i]==2:
        twoE.append(i)
        twoECount.append(countedPeaks2[i])
    elif elinfo[i]==3:
        threeE.append(i)
        threeECount.append(countedPeaks2[i])
    elif elinfo[i]==4:
        fourE.append(i)
        fourECount.append(countedPeaks2[i])
    else:
        fiveE.append(i)
        fiveECount.append(countedPeaks2[i])
np.savetxt('/home/jason/data/amanda/smear9FinalOuts/oneE_shifted2.txt',np.c_[oneE,oneECount])
np.savetxt('/home/jason/data/amanda/smear9FinalOuts/twoE_shifted2.txt',np.c_[twoE,twoECount])
np.savetxt('/home/jason/data/amanda/smear9FinalOuts/threeE_shifted2.txt',np.c_[threeE,threeECount])
np.savetxt('/home/jason/data/amanda/smear9FinalOuts/fourE_shifted2.txt',np.c_[fourE,fourECount])
np.savetxt('/home/jason/data/amanda/smear9FinalOuts/fiveE_shifted2.txt',np.c_[fiveE,fiveECount])

badEvents=[]
badTruth=[]
badCount=[]

for i in range(NmbEvents):
    if elinfo[i]!=countedPeaks2[i]:
        badEvents.append(i)
        badTruth.append(elinfo[i])
        badCount.append(countedPeaks2[i])
np.savetxt('/home/jason/data/amanda/smear9FinalOuts/miscounted_shifted2.txt',np.c_[badEvents,badTruth,badCount])
"""
########################################################
###########   2-electron event separation analysis - for simulation
########################################################
separation=np.zeros(NmbEventsSim)

if sim:
    for i in range(NmbEventsSim):
        if elinfo[i]==2:
            separation[i]=distance(positionArray2[i-multEArray[0]][0][0],positionArray2[i-multEArray[0]][0][1],positionArray2[i-multEArray[0]][1][0],positionArray2[i-multEArray[0]][1][1])

########################################################
###########   FOR Graphing
########################################################

#### Lets do some graphing now that we supposedly know the counts
#Set up graphs

PeakGraphsSim =[0,0,0,0,0,0,0,0,0]
PeakGraphs1Sim =[0,0,0,0,0,0,0,0,0]
PeakGraphs2Sim=[0,0,0,0,0,0,0,0,0]
PeakGraphsTest =[0,0,0,0,0,0,0,0,0]
PeakGraphs1Test =[0,0,0,0,0,0,0,0,0]
PeakGraphs2Test=[0,0,0,0,0,0,0,0,0]

bins=100
xmin=0
xmax=1000

if sim:
    for i in range(9):
        PeakGraphsSim[i] = TH1D(str(i)+" Electron Events - SumEnergy",str(i)+" Electron Events - SumEnergy",bins,xmin,xmax)
        PeakGraphs1Sim[i] = TH1D(str(i)+" Electron Events - SumEnergy",str(i)+" Electron Events - SumEnergy",bins,xmin,xmax)
        PeakGraphs2Sim[i] = TH1D(str(i)+" Electron Events - SumEnergy",str(i)+" Electron Events - SumEnergy",bins,xmin,xmax)
        PeakGraphsSim[i].GetXaxis().SetTitle("Total Deposited Charge (x10e-14 C)")
        PeakGraphs2Sim[i].GetXaxis().SetTitle("Total Deposited Charge (x10e-14 C)")
    peakLegendSim=TLegend(0.5,0.5,0.8,0.8)
    OverCount =  TH1D("Overcount","Overcount - All Events",bins/2,0,xmax/3)
    UnderCount =  TH1D("Undercount","Undercount - All Events",bins/2,0,xmax/3)
    CorrCount =  TH1D("Correct Count","All Events",bins/2,0,xmax/3)
    CorrCount.GetXaxis().SetTitle("Total Deposited Energy (MeV)")
    totLegend=TLegend(0.5,0.75,0.8,0.9)
    OverCount2 =  TH1D("Overcount","Overcount - 2-electron events",bins/2,0,xmax/3)
    UnderCount2 =  TH1D("Undercount","Undercount - 2-electron events",bins/2,0,xmax/3)
    CorrCount2 =  TH1D("Correct count"," 2-electron events",bins/2,0,xmax/3)
    CorrCount2.GetXaxis().SetTitle("Total Deposited Energy (MeV)")
    twoLegend=TLegend(0.4,0.75,0.9,0.9)
    SepHistU = TH1D("Separation-Under","Separation-Under",bins/12.5,xmin,xmax/20)
    SepHistO = TH1D("Separation-Over","Separation-Over",bins/12.5,xmin,xmax/20)
    SepHistC = TH1D("Separation Correct","2 Electron Separation",bins/12.5,xmin,xmax/20)
    SepHistC.GetXaxis().SetTitle("Separation [mm]")
    sepLegend=TLegend(0.5,0.75,0.7,0.9)
#    sep2TotHist=TH1D("Separation","Truth Separation of 2-electron events; Separation (mm)",bins/5.2,xmin,xmax/13)
    sep2TotEffHist=TH1D("Separation","Truth Separation of 2-electron events; Separation (mm)",bins/2,xmin,xmax/20)
    sep2CorrNormHist=TH1D("Separation","Separation of 2-electron correctly tagged events - Normalized to 100 Events",bins/5.2,xmin,xmax/13)
    sep2CorrNormHist.GetXaxis().SetTitle("Separation [mm]")
    sep2TagHist=TH1D("Separation","Separation of 2-electron tagged events;Separation (mm)",bins/5.2,xmin,xmax/13)
    sep2TruthNormHist=TH1D("Separation","Separation of 2-electron truth events; Separation [mm]; Entries (Normalized to 100)",bins/5.2,xmin,xmax/13)
    sep2TagNormHist=TH1D("Separation","Separation of 2-electron tagged events; Separation [mm]; Entries (Normalized to 100)",bins/5.2,xmin,xmax/13)
    sepEfficiency=TH1D("Correct count","2-electron Counting Efficiency",bins/2,xmin,xmax/20)
    sepEfficiency.GetXaxis().SetTitle("Separation [mm]")
    sepEfficiency1=TH1D("Correct count","2-electron Counting Efficiency",bins/2,xmin,xmax/20)
    sep2SimLegend=TLegend(0.5,0.3,0.8,0.4)
#    sep2Legend=TLegend(0.5,0.3,0.8,0.4)
    multEHistSim=[0,0,0,0,0,0,0,0,0]
    for i in range(9):
        multEHistSim[i]=TH1D(str(i)+" Electrons","Electron Events - Simulation Counted",2.5*bins,xmin,xmax)
        multEHistSim[i].GetXaxis().SetTitle("Deposited Total Charge (x10e-14 C)")
    multELegendSim=TLegend(0.65,0.5,0.85,0.7)
if test or revTest:
    for i in range(9):
        PeakGraphsTest[i] = TH1D(str(i)+" Electron Events - SumEnergy",str(i)+" Electron Events - SumEnergy",bins,xmin,xmax)
        PeakGraphs1Test[i] = TH1D(str(i)+" Electron Events - SumEnergy",str(i)+" Electron Events - SumEnergy",bins,xmin,xmax)
        PeakGraphs2Test[i] = TH1D(str(i)+" Electron Events - SumEnergy",str(i)+" Electron Events - SumEnergy",bins,xmin,xmax)
        PeakGraphsTest[i].GetXaxis().SetTitle("Total Deposited Charge (x10e-14 C)")
        PeakGraphs2Test[i].GetXaxis().SetTitle("Total Deposited Charge (x10e-14 C)")
    peakLegendTest=TLegend(0.5,0.5,0.8,0.8)
    sep2Hist=TH1D("Separation","Separation of Tagged 2-electron Test Beam Events",bins/5.2,xmin,xmax/13)
    sep2Hist.GetXaxis().SetTitle("Separation (mm)")
    sep2NormHist=TH1D("Separation","Separation of 2-electron Events; Separation [mm]; Entries (Normalized to 100)",bins/5.2,xmin,xmax/13)
#    sep2NormHist.GetXaxis().SetTitle("Separation (mm)")
    multEHistTest=[0,0,0,0,0,0,0,0,0]
    for i in range(9):
        multEHistTest[i]=TH1D(str(i)+" Electrons","Electron Events - Test Beam Counted",2.5*bins,xmin,xmax)
        multEHistTest[i].GetXaxis().SetTitle("Deposited Total Charge (x10e-14 C)")
    multELegendTest=TLegend(0.65,0.5,0.85,0.7)

sep2ComboLegend=TLegend(0.5,0.6,0.8,0.7)

#start variables for sim statistics
if sim:
    overc=0.
    underc=0.
    corrc=0.
    overc2=0.
    underc2=0.
    corrc2=0.

#sort events into appropriate histograms
if sim:
    for i in range(NmbEventsSim):
        PeakGraphsSim[int(countedPeaksSim[i])].Fill(endataSim[i].sum()*2.9575)
        PeakGraphs1Sim[int(countedPeaks1Sim[i])].Fill(endataSim[i].sum()*2.9575)
        PeakGraphs2Sim[int(countedPeaks2Sim[i])].Fill(endataSim[i].sum()*2.9575)
        if countedPeaks2Sim[i] != elinfo[i]:
            if countedPeaks2Sim[i] > elinfo[i]:
                OverCount.Fill(endataSim[i].sum())
                overc+=1
            else:
                UnderCount.Fill(endataSim[i].sum())
                underc+=1
        else:
            CorrCount.Fill(endataSim[i].sum())
            corrc+=1
        if elinfo[i]==2:
#            sep2TotHist.Fill(separation[i])
            sep2TotEffHist.Fill(separation[i])
            sep2TruthNormHist.Fill(separation[i],1/34.90)
            if countedPeaks2Sim[i]>2:
                OverCount2.Fill(endataSim[i].sum())
                SepHistO.Fill(separation[i])
                overc2+=1
            elif countedPeaks2Sim[i]<2:
                UnderCount2.Fill(endataSim[i].sum())
                SepHistU.Fill(separation[i])
                underc2+=1
            else:
                CorrCount2.Fill(endataSim[i].sum())
                SepHistC.Fill(separation[i])
                sepEfficiency.Fill(separation[i])
                #sep2CorrNormHist.Fill(separation[i],1/28.84)
                sep2CorrNormHist.Fill(separation[i],1/34.9)
                corrc2+=1
            if countedPeaks1Sim[i]==2:
                sepEfficiency1.Fill(separation[i])
        if countedPeaks2Sim[i]==0:
            multEHistSim[0].Fill(endataSim[i].sum()*2.9575,1/127.64)
        elif countedPeaks2Sim[i]==1:
            multEHistSim[1].Fill(endataSim[i].sum()*2.9575,1/127.64)
        elif countedPeaks2Sim[i]==2:
            multEHistSim[2].Fill(endataSim[i].sum()*2.9575,1/127.64)
            sep2TagHist.Fill(separation[i])
            sep2TagNormHist.Fill(separation[i],1/32.96)
        elif countedPeaks2Sim[i]==3:
            multEHistSim[3].Fill(endataSim[i].sum()*2.9575,1/127.64)
        elif countedPeaks2Sim[i]==4:
            multEHistSim[4].Fill(endataSim[i].sum()*2.9575,1/127.64)
        elif countedPeaks2Sim[i]==5:
            multEHistSim[5].Fill(endataSim[i].sum()*2.9575,1/127.64)
        else:
            multEHistSim[6].Fill(endataSim[i].sum()*2.9575,1/127.64)
    print "correct tagged 2e events = "+str(corrc2)

"""
cOther=0
c10c=0
c11c=0
c12c=0
c13c=0
c14c=0
c20c=0
c21c=0
c22c=0
c23c=0
c24c=0
c25c=0
c30c=0
c31c=0
c32c=0
c33c=0
c34c=0
c35c=0
c36c=0
c40c=0
c41c=0
c42c=0
c43c=0
c44c=0
c45c=0
c46c=0
c50c=0
c51c=0
c52c=0
c53c=0
c54c=0
c55c=0
c56c=0
if sim:
    for i in range(NmbEvents):
        if elinfo[i]==1:
            if countedPeaks2[i]==0:
                c10c+=1
            elif countedPeaks2[i]==1:
                c11c+=1
            elif countedPeaks2[i]==2:
                c12c+=1
            elif countedPeaks2[i]==3:
                c13c+=1
            elif countedPeaks2[i]==4:
                c14c+=1
        elif elinfo[i]==2:
            if countedPeaks2[i]==0:
                c20c+=1
            elif countedPeaks2[i]==1:
                c21c+=1
            elif countedPeaks2[i]==2:
                c22c+=1
            elif countedPeaks2[i]==3:
                c23c+=1
            elif countedPeaks2[i]==4:
                c24c+=1
            elif countedPeaks2[i]==5:
                c25c+=1
        elif elinfo[i]==3:
            if countedPeaks2[i]==0:
                c30c+=1
            elif countedPeaks2[i]==1:
                c31c+=1
            elif countedPeaks2[i]==2:
                c32c+=1
            elif countedPeaks2[i]==3:
                c33c+=1
            elif countedPeaks2[i]==4:
                c34c+=1
            elif countedPeaks2[i]==5:
                c35c+=1
            elif countedPeaks2[i]==6:
                c36c+=1
        elif elinfo[i]==4:
            if countedPeaks2[i]==0:
                c40c+=1
            elif countedPeaks2[i]==1:
                c41c+=1
            elif countedPeaks2[i]==2:
                c42c+=1
            elif countedPeaks2[i]==3:
                c43c+=1
            elif countedPeaks2[i]==4:
                c44c+=1
            elif countedPeaks2[i]==5:
                c45c+=1
            elif countedPeaks2[i]==6:
                c46c+=1
        elif elinfo[i]==5:
            if countedPeaks2[i]==0:
                c50c+=1
            elif countedPeaks2[i]==1:
                c51c+=1
            elif countedPeaks2[i]==2:
                c52c+=1
            elif countedPeaks2[i]==3:
                c53c+=1
            elif countedPeaks2[i]==4:
                c54c+=1
            elif countedPeaks2[i]==5:
                c55c+=1
            elif countedPeaks2[1]==6:
                c56c+=1
        else:
            cOther+=1
print "10 :"+str(c10c)
print "11 :"+str(c11c)
print "12 :"+str(c12c)
print "13 :"+str(c13c)
print "14 :"+str(c14c)
print "20 :"+str(c20c)
print "21 :"+str(c21c)
print "22 :"+str(c22c)
print "23 :"+str(c23c)
print "24 :"+str(c24c)
print "25 :"+str(c25c)
print "30 :"+str(c30c)
print "31 :"+str(c31c)
print "32 :"+str(c32c)
print "33 :"+str(c33c)
print "34 :"+str(c34c)
print "35 :"+str(c35c)
print "36 :"+str(c36c)
print "40 :"+str(c40c)
print "41 :"+str(c41c)
print "42 :"+str(c42c)
print "43 :"+str(c43c)
print "44 :"+str(c44c)
print "45 :"+str(c45c)
print "46 :"+str(c46c)
print "50 :"+str(c50c)
print "51 :"+str(c51c)
print "52 :"+str(c52c)
print "53 :"+str(c53c)
print "54 :"+str(c54c)
print "55 :"+str(c55c)
print "56 :"+str(c56c)
print "other :"+str(cOther)
"""

if test or revTest:
    for i in range(len(sep2)):
        sep2Hist.Fill(sep2[i])
        sep2NormHist.Fill(sep2[i],1/74.96)
    for i in range(NmbEventsTest):
        if i not in cleanedEvents:
            PeakGraphsTest[int(countedPeaksTest[i])].Fill(endataTest[i].sum()*100000000000000)
            PeakGraphs1Test[int(countedPeaks1Test[i])].Fill(endataTest[i].sum()*100000000000000)
            PeakGraphs2Test[int(countedPeaks2Test[i])].Fill(endataTest[i].sum()*100000000000000)
            if countedPeaks2Test[i]==0:
                multEHistTest[0].Fill(endataTest[i].sum()*100000000000000,1/176.16)
            elif countedPeaks2Test[i]==1:
                multEHistTest[1].Fill(endataTest[i].sum()*100000000000000,1/176.16)
            elif countedPeaks2Test[i]==2:
                multEHistTest[2].Fill(endataTest[i].sum()*100000000000000,1/176.16)
            elif countedPeaks2Test[i]==3:
                multEHistTest[3].Fill(endataTest[i].sum()*100000000000000,1/176.16)
            elif countedPeaks2Test[i]==4:
                multEHistTest[4].Fill(endataTest[i].sum()*100000000000000,1/176.16)
            elif countedPeaks2Test[i]==5:
                multEHistTest[5].Fill(endataTest[i].sum()*100000000000000,1/176.16)
            else:
                multEHistTest[6].Fill(endataTest[i].sum()*100000000000000,1/176.16)


##draw histograms
canvas6=TCanvas('canvas6','2e Separation Comparison',1200,800)
canvas6.Divide(2,2)
canvas6.cd(1)
#sep2TotHist.Draw()
#sep2TotHist.SetStats(0)
#sep2TagHist.Draw("same")
#sep2TagHist.SetStats(0)
SepHistC.Draw()
#SepHistC.SetStats(0)
#sep2TagHist.SetLineColor(ROOT.kGreen+3)
#SepHistC.SetLineColor(ROOT.kRed)
#sep2SimLegend.AddEntry(sep2TotHist,"True 2e events","l")
#sep2SimLegend.AddEntry(sep2TagHist,"Tagged 2e events","l")
#sep2SimLegend.AddEntry(SepHistC,"Correctly tagged 2e events","l")
#sep2SimLegend.Draw()
canvas6.cd(2)
sep2Hist.Draw()
sep2Hist.SetLineColor(ROOT.kRed)
canvas6.cd(3)
sep2NormHist.Draw()
sep2NormHist.SetStats(0)
sep2CorrNormHist.Draw("same")
sep2CorrNormHist.SetStats(0)
sep2TruthNormHist.Draw("same")
sep2TruthNormHist.SetStats(0)
#sep2TagNormHist.Draw("same")
#sep2TagNormHist.SetStats(0)
sep2CorrNormHist.SetLineColor(ROOT.kBlue)
sep2NormHist.SetLineColor(ROOT.kRed)
sep2TruthNormHist.SetLineColor(ROOT.kBlue)
sep2TruthNormHist.SetLineStyle(7)
#sep2TagNormHist.SetLineColor(ROOT.kBlue)
#sep2TagNormHist.SetLineStyle(3)
sep2ComboLegend.AddEntry(sep2TruthNormHist,"Simulation Truth","l")
sep2ComboLegend.AddEntry(sep2CorrNormHist,"Simulation Correctly Tagged","l")
#sep2ComboLegend.AddEntry(sep2TagNormHist,"Simulation Tagged","l")
sep2ComboLegend.AddEntry(sep2NormHist,"43 Test Beam Tagged","l")
sep2ComboLegend.Draw()

if sim:
    canvas = TCanvas( 'aCanvas', 'distribution', 1350, 700 )
    canvas.Divide(3,3)
    for i in range(9):
        canvas.cd(i+1)
        PeakGraphsSim[i].SetLineColor(ROOT.kRed)
        PeakGraphs1Sim[i].SetLineColor(ROOT.kOrange)
        PeakGraphs2Sim[i].SetLineColor(ROOT.kGreen)
        if i==0:
            PeakGraphsSim[0].Draw()
            PeakGraphsSim[0].SetStats(0)
            PeakGraphs2Sim[0].Draw("same")
            PeakGraphs2Sim[0].SetStats(0)
            PeakGraphs1Sim[0].Draw("same")
            PeakGraphs1Sim[0].SetStats(0)
        else:
            PeakGraphs2Sim[i].Draw()
            PeakGraphs2Sim[i].SetStats(0)
            PeakGraphsSim[i].Draw("same")
            PeakGraphsSim[i].SetStats(0)
            PeakGraphs1Sim[i].Draw("same")
            PeakGraphs1Sim[i].SetStats(0)
    peakLegendSim.AddEntry(PeakGraphsSim[6],"Dylan's algorithm","l")
    peakLegendSim.AddEntry(PeakGraphs1Sim[6],"Jason's algorithm","l")
    peakLegendSim.AddEntry(PeakGraphs2Sim[6],"Amanda's algorithm","l")
    peakLegendSim.Draw()
    
    canvas1 = TCanvas( 'bCanvas', 'all event counts', 700, 500 )
    CorrCount.Draw()
    CorrCount.SetStats(0)
    OverCount.SetLineColor(ROOT.kRed)
    UnderCount.SetLineColor(ROOT.kGreen-2)
    OverCount.Draw("same")
    OverCount.SetStats(0)
    UnderCount.Draw("same")
    UnderCount.SetStats(0)
    corrPerc="{0:.2f}".format(round(corrc/NmbEventsSim*100.,2))#to display percentages in legend
    overPerc="{0:.2f}".format(round(overc/NmbEventsSim*100.,2))
    underPerc="{0:.2f}".format(round(underc/NmbEventsSim*100.,2))
    totLegend.AddEntry(CorrCount,"Correctly counted - "+str(corrPerc)+"%","l")
    totLegend.AddEntry(OverCount,"Overcounted - "+str(overPerc)+"%","l")
    totLegend.AddEntry(UnderCount,"Undercounted - "+str(underPerc)+"%","l")
    totLegend.Draw()

    canvas2=TCanvas('cCanvas','2e analysis',900,600)
    canvas2.Divide(2,1)
    canvas2.cd(1)
    CorrCount2.Draw()
    CorrCount2.SetStats(0)
    OverCount2.SetLineColor(ROOT.kRed)
    UnderCount2.SetLineColor(ROOT.kGreen-2)
    OverCount2.Draw("same")
    OverCount2.SetStats(0)
    UnderCount2.Draw("same")
    UnderCount2.SetStats(0)
    corrPerc2="{0:.2f}".format(round(corrc2/multEArray[1]*100.,2))#to display percentages in legend
    overPerc2="{0:.2f}".format(round(overc2/multEArray[1]*100.,2))
    underPerc2="{0:.2f}".format(round(underc2/multEArray[1]*100.,2))
    twoLegend.AddEntry(CorrCount2,"Correctly counted - "+str(corrPerc2)+"%","l")
    twoLegend.AddEntry(OverCount2,"Overcounted - "+str(overPerc2)+"%","l")
    twoLegend.AddEntry(UnderCount2,"Undercounted - "+str(underPerc2)+"%","l")
    twoLegend.Draw()
    canvas2.cd(2)
    SepHistU.SetLineColor(ROOT.kGreen-2)
    SepHistO.SetLineColor(ROOT.kRed)
    SepHistC.Draw()
#    SepHistC.SetStats(0)
    SepHistU.Draw("same")
    SepHistU.SetStats(0)
    SepHistO.Draw("same")
    SepHistO.SetStats(0)
    sepLegend.AddEntry(SepHistC,"Correctly counted","l")
    sepLegend.AddEntry(SepHistO,"Overcounted","l")
    sepLegend.AddEntry(SepHistU,"Undercounted","l")
    sepLegend.Draw()

    canvas4=TCanvas('canvas4','Separation',900,600)
    canvas4.Divide(2,1)
    canvas4.cd(1)
    sep2TotEffHist.Draw()
    sep2TotEffHist.SetFillColor(ROOT.kBlue-2)
    canvas4.cd(2)
    #sepEfficiency1.Divide(sep2Hist)
    #sepEfficiency1.Draw()
    #sepEfficiency1.SetLineColor(ROOT.kRed)
    sepEfficiency.Divide(sep2TotEffHist)
    sepEfficiency.Draw()
    sepEfficiency.SetStats(0)
    sepEfficiency.SetLineColor(ROOT.kBlue)
    #sep2Legend.AddEntry(sepEfficiency,"Amanda's Algorithm","l")
    #sep2Legend.AddEntry(sepEfficiency1,"Jason's Algorithm","l")
    #sep2Legend.Draw()
    canvas3=TCanvas('canvas3','Total',600,350)
#canvas3.SetLogy()
    multELegendSim.AddEntry(multEHistSim[0],"0 electrons","f")
    multEHistSim[1].Draw()#want biggest hist in the back
    multEHistSim[1].SetStats(0)
    multEHistSim[1].SetFillColor(ROOT.kGreen+3)
    multELegendSim.AddEntry(multEHistSim[1],"1 electron","f")
    elec=2
    for i in range(2,7):
        multEHistSim[i].Draw("same")
        multEHistSim[i].SetStats(0)
        multEHistSim[i].SetFillColor(2*elec-2)
        multELegendSim.AddEntry(multEHistSim[i],str(elec)+" electrons","f")
        elec+=1
    multEHistSim[0].Draw("same")
    multEHistSim[0].SetStats(0)
    multEHistSim[0].SetFillColor(ROOT.kGreen)
    multELegendSim.Draw()

if test or revTest:
    canvas10 = TCanvas( 'aCanvas10', 'distribution', 1350, 700 )
    canvas10.Divide(3,3)
    for i in range(9):
        canvas10.cd(i+1)
        PeakGraphsTest[i].SetLineColor(ROOT.kRed)
        PeakGraphs1Test[i].SetLineColor(ROOT.kOrange)
        PeakGraphs2Test[i].SetLineColor(ROOT.kGreen)
        if i==0:
            PeakGraphsTest[0].Draw()
            PeakGraphsTest[0].SetStats(0)
            PeakGraphs2Test[0].Draw("same")
            PeakGraphs2Test[0].SetStats(0)
            PeakGraphs1Test[0].Draw("same")
            PeakGraphs1Test[0].SetStats(0)
        else:
            PeakGraphs2Test[i].Draw()
            PeakGraphs2Test[i].SetStats(0)
            PeakGraphsTest[i].Draw("same")
            PeakGraphsTest[i].SetStats(0)
            PeakGraphs1Test[i].Draw("same")
            PeakGraphs1Test[i].SetStats(0)
    peakLegendTest.AddEntry(PeakGraphsTest[6],"Dylan's algorithm","l")
    peakLegendTest.AddEntry(PeakGraphs1Test[6],"Jason's algorithm","l")
    peakLegendTest.AddEntry(PeakGraphs2Test[6],"Amanda's algorithm","l")
    peakLegendTest.Draw()

#    canvas5=TCanvas('canvas5',"Test Beam 2e Separation",600,400)
#    sep2Hist.Draw()
    canvas7=TCanvas('canvas7','Total',600,350) 
    multELegendTest.AddEntry(multEHistTest[0],"0 electrons","f")
    multEHistTest[1].Draw()#want biggest hist in the back 
    multEHistTest[1].SetStats(0)
    multEHistTest[1].SetFillColor(ROOT.kGreen+3)
    multELegendTest.AddEntry(multEHistTest[1],"1 electron","f")
    elec=2
    for i in range(2,7):
        multEHistTest[i].Draw("same")
        multEHistTest[i].SetStats(0)
        multEHistTest[i].SetFillColor(2*elec-2)
        multELegendTest.AddEntry(multEHistTest[i],str(elec)+" electrons","f")
        elec+=1
    multEHistTest[0].Draw("same")
    multEHistTest[0].SetStats(0)
    multEHistTest[0].SetFillColor(ROOT.kGreen)
    multELegendTest.Draw()

"""
canvas6=TCanvas('canvas6','2e Separation Comparison',1200,400)
canvas6.Divide(2,2)
canvas6.cd(1)
CorrCount2.Draw()
canvas6.cd(2)
sep2Hist.Draw()
sep2Hist.SetLineColor(ROOT.kRed)
canvas6.cd(3)
sep2NormHist.Draw()
sep2NormHist.SetStats(0)
sep2TotNormHist.Draw("same")
sep2TotNormHist.SetStats(0)
sep2TotNormHist.SetLineColor(ROOT.kBlue)
sep2NormHist.SetLineColor(ROOT.kRed)
sep2ComboLegend.AddEntry(sep2TotNormHist,"Simulation Tagged")
sep2ComboLegend.AddEntry(sep2NormHist,"43 Test Beam Tagged")
sep2ComboLegend.Draw()
"""
