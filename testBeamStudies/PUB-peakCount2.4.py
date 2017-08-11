import ROOT
from ROOT import TCanvas, TH1D, TGraph, TFile, TTree, TH2D, TLegend
import pickle
import numpy as np
from array import array
import math
import glob

""" Dylan Mead's Peak Counting algorithm rewritten in Python by Jason Barkeloo-Univ of Oregon"""
sim=False
if sim:
    siFirst=True #if False, W-First
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

##~~~~~~~~~~~~~~~~~~~~~~~
# Input Simulation Data
##~~~~~~~~~~~~~~~~~~~~~~~
if sim:
    inFile = '/home/jason/data/amanda/smear9FinalOuts/DENS23-SiWNi/withPos/shifted/12.1GeV-smear9-SiWNi(0.8725)_shifted23.endataIMP.pkl'
    inFile1 = '/home/jason/data/amanda/smear9FinalOuts/DENS23-SiWNi/withPos/shifted/12.1GeV-smear9-SiWNi(0.8725)_shifted23_elinfoIMP.pkl'
    unpicklefile = open(inFile, 'r')
    endata = pickle.load(unpicklefile)
    unpicklefile = open(inFile1,'r')
    elinfo = pickle.load(unpicklefile)
    NmbEvents=len(endata)
    print NmbEvents
    positionArray2=pickle.load(open('/home/jason/data/amanda/smear9FinalOuts/DENS23-SiWNi/withPos/shifted/12.1GeV-smear9-SiWNi(0.8725)_2shifted23.positions.pkl','rb'))
    multEArray=[8000,3490,1015,221,38] #SiWNi
#    multEArray=[8000,3092,796,153,23]#WNiSi
    
    print inFile
    print inFile1

##~~~~~~~~~~~~~~~~~~~~~~~~
# Input Test Beam Data
##~~~~~~~~~~~~~~~~~~~~~~~~
if test or revTest:
#    inFile = '/home/jason/data/TestBeam/TestBeamROOT/2013_07_26_14_13_43.bin.endata.pkl'
    inFile='/home/jason/data/TestBeam/TestBeamROOT/2013_07_26_10_41_28.bin.endata.pkl'
    unpicklefile = open(inFile, 'r')
    endata = pickle.load(unpicklefile)
    NmbEvents=len(endata)

print inFile

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Input information about pixel neighbors
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#this text file lists each pixel and all its neighboring pixels with GUI ID number
neighborFile = '/home/jason/jas-assembly-3.1.2/extensions/BarkelooCode/resources/neighborsNewest.txt'
neighbors = listoflist(1024)

testpnt = 0
for i in open(neighborFile,'r'):
    neighbors[testpnt]=i.split()
    testpnt+=1
for i in range(len(neighbors)):
	for j in range(len(neighbors[i])):
		neighbors[i][j]=int(neighbors[i][j])#change to type int
# neighbors is now a list with index of a kpixChannel and each entry is a list of its' neighboring cells

# Basic memory management
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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clean low-energy photon contamination from test beam data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if test or revTest:
    if test:
        cut=44
    else:
        cut=40
    cleanedEvents=[]
    for i in range(NmbEvents):
        if i%1000==0:
            print "Test Beam Event: "+str(i)+" of "+str(NmbEvents)
        totalEnergy=endata[i].sum()*1000000000000000
        statNum=0
        statDenom=totalEnergy
        nullCount=0
        for j in range(9):
            for k in range(1024):
                if test:
                    statNum+=endata[i][j][k]*(j+1)*(j+1)*1000000000000000
                if revTest:
                    statNum+=endata[i][j][k]*(9-j)*(9-j)*1000000000000000
        for j in range(9):
            if endata[i][j].sum()==0:
                nullCount+=1
                if test:
                    statNum+=4.*(j+1)*(j+1)
                if revTest:
                    statNum+=4.*(9-j)*(9-j)
                statDenom+=4.
        if nullCount==8 or statNum/statDenom<cut:
            cleanedEvents.append(i)

######################################################################################################################
#method to find local maxima of charge deposits in each layer and determine using those how many incident electrons were in the event
def peakCount(kpixArray, neighbor):
    countedPeaks2=np.zeros(NmbEvents)
    arrayMaxes=[]
    for i in range(NmbEvents):#i is event number
        if i%1000==0:
            print "Event: "+str(i)+" of "+str(NmbEvents)
        detector = kpixArray[i] #should be 9x1024 array accordingly
        #Find local maxes of deposited charge in each layer
        arrayMaxes =[[],[],[],[],[],[],[],[],[]] #get maxes for each layer individually 
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
                        arrayMaxes[j].append(k) #append channels that are local maxes
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

    print "Counted: "+ str(sum(countedPeaks2)) +" electrons"
    return countedPeaks2
################################################################################################################################

#submit the data to be counted
countedPeaks2= peakCount(endata, neighbors)

#check accuracy of counting method
truecount2=0
if sim:
    for i in range(len(countedPeaks2)):
        if countedPeaks2[i]==elinfo[i]:
            truecount2+=1.    
    print "Counted "+str(truecount2)+" of "+str(len(elinfo))+" events correctly"


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute true separation of simulated 2-electron events (for simulation only)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if sim:
    separation=np.zeros(NmbEvents)
    for i in range(NmbEvents):
        if elinfo[i]==2:
            separation[i]=distance(positionArray2[i-multEArray[0]][0][0],positionArray2[i-multEArray[0]][0][1],positionArray2[i-multEArray[0]][1][0],positionArray2[i-multEArray[0]][1][1])

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create histograms and analyze output to inform plotting
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PeakGraphs2=[0]*9

bins=100
xmin=0
xmax=1000
if sim and siFirst:
    convert=29.575
elif sim:
    convert=29.39

if sim:
    OverCount =  TH1D("Overcount","Overcount - All Events",bins,0,xmax*10)
    UnderCount =  TH1D("Undercount","Undercount - All Events",bins,0,xmax*10)
    CorrCount =  TH1D("Correct Count","Counting of All Simulated Electron Events;Total Measured Charge [fC];Entries",bins,0,xmax*10)
    totLegend=TLegend(0.5,0.75,0.8,0.9)
    OverCount2 =  TH1D("Overcount","Overcount - 2-electron events",bins*4/5,0,xmax*8)
    UnderCount2 =  TH1D("Undercount","Undercount - 2-electron events",bins*4/5,0,xmax*8)
    CorrCount2 =  TH1D("Correct count","Counting of 2-Electron Simulated Events;Total Measured Charge [fC];Entries",bins*4/5,0,xmax*8)
    twoLegend=TLegend(0.4,0.75,0.9,0.9)
    SepHistU = TH1D("Separation-Under","Separation-Under",bins/2,xmin,xmax/20)
    SepHistO = TH1D("Separation-Over","Separation-Over",bins/2,xmin,xmax/20)
    SepHistC = TH1D("Separation Correct","Counting of Simulated 2-Electron Events by Separation;Separation [mm];Entries",bins/2,xmin,xmax/20)
    sepLegend=TLegend(0.5,0.75,0.9,0.9)
    sep2Hist=TH1D("Separation","Separation of 2-electron events",bins/2,xmin,xmax/20)
    sepEfficiency=TH1D("Correct count","Simulated 2-Electron Event Counting Efficiency;Separation [mm];Entries",bins/2,xmin,xmax/20)
    sepEfficiency1=TH1D("Correct count","Simulated 2-Electron Event Counting Efficiency",bins/2,xmin,xmax/20)
    sep2Legend=TLegend(0.5,0.3,0.8,0.4)
multEHist=[0,0,0,0,0,0,0,0,0]
for i in range(9):
    if sim:
        multEHist[i]=TH1D(str(i)+" Electrons","Electron Events - Simulation Counted;Total Measured Charge [fC];Entries (normalized to 100)",100,0,10000)
    if test or revTest:
        multEHist[i]=TH1D(str(i)+" Electrons","Electron Events - Test Beam Counted; Total Measured Charge; Entries (normalized to 100)",100,0,10000)
multELegend=TLegend(0.65,0.5,0.85,0.7)

#calculate statistics regarding counting performance
if sim:
    overc=0.
    underc=0.
    corrc=0.
    overc2=0.
    underc2=0.
    corrc2=0.

#sort events into appropriate histograms
for i in range(NmbEvents):
    if sim:
        if countedPeaks2[i] != elinfo[i]:
            if countedPeaks2[i] > elinfo[i]:
                OverCount.Fill(endata[i].sum()*convert)
                overc+=1
            else:
                UnderCount.Fill(endata[i].sum()*convert)
                underc+=1
        else:
            CorrCount.Fill(endata[i].sum()*convert)
            corrc+=1
if sim:
    for i in range(NmbEvents):
        if elinfo[i]==2:
            sep2Hist.Fill(separation[i])
            if countedPeaks2[i]>2:
                OverCount2.Fill(endata[i].sum()*convert)
                SepHistO.Fill(separation[i])
                overc2+=1
            elif countedPeaks2[i]<2:
                UnderCount2.Fill(endata[i].sum()*convert)
                SepHistU.Fill(separation[i])
                underc2+=1
            else:
                CorrCount2.Fill(endata[i].sum()*convert)
                SepHistC.Fill(separation[i])
                sepEfficiency.Fill(separation[i])
                corrc2+=1

for i in range(NmbEvents):
    if sim:
        if countedPeaks2[i]==0:
            multEHist[0].Fill(endata[i].sum()*convert,1/127.64)
        elif countedPeaks2[i]==1:
            multEHist[1].Fill(endata[i].sum()*convert,1/127.64)
        elif countedPeaks2[i]==2:
            multEHist[2].Fill(endata[i].sum()*convert,1/127.64)
        elif countedPeaks2[i]==3:
            multEHist[3].Fill(endata[i].sum()*convert,1/127.64)
        elif countedPeaks2[i]==4:
            multEHist[4].Fill(endata[i].sum()*convert,1/127.64)
        elif countedPeaks2[i]==5:
            multEHist[5].Fill(endata[i].sum()*convert,1/127.64)
        else:
            multEHist[6].Fill(endata[i].sum()*convert,1/127.64)
    if test or revTest:
        if i not in cleanedEvents:
            if countedPeaks2[i]==0:
                multEHist[0].Fill(endata[i].sum()*1000000000000000,1/176.16)
            elif countedPeaks2[i]==1:
                multEHist[1].Fill(endata[i].sum()*1000000000000000,1/176.16)
            elif countedPeaks2[i]==2:
                multEHist[2].Fill(endata[i].sum()*1000000000000000,1/176.16)
            elif countedPeaks2[i]==3:
                multEHist[3].Fill(endata[i].sum()*1000000000000000,1/176.16)
            elif countedPeaks2[i]==4:
                multEHist[4].Fill(endata[i].sum()*1000000000000000,1/176.16)
            elif countedPeaks2[i]==5:
                multEHist[5].Fill(endata[i].sum()*1000000000000000,1/176.16)
            else:
                multEHist[6].Fill(endata[i].sum()*1000000000000000,1/176.16)

##~~~~~~~~~~~~~~~~~~~~
# Draw Histograms
##~~~~~~~~~~~~~~~~~~~~
if sim:
    allCountedCanvas = TCanvas( 'allCountedCanvas', 'Counting Accuracy of All Electron Events', 800, 700)
    CorrCount.Draw()
    CorrCount.SetFillColor(ROOT.kBlue)
    CorrCount.SetFillStyle(3005)
    CorrCount.SetStats(0)
    UnderCount.SetLineColor(ROOT.kGreen-2)
    UnderCount.SetFillColor(ROOT.kGreen-2)
    UnderCount.SetFillStyle(3004)
    OverCount.SetFillColor(ROOT.kRed)
    OverCount.Draw("same")
    OverCount.SetStats(0)
    UnderCount.Draw("same")
    UnderCount.SetStats(0)
    corrPerc="{0:.2f}".format(round(corrc/NmbEvents*100.,2))#to display percentages in legend
    overPerc="{0:.2f}".format(round(overc/NmbEvents*100.,2))
    underPerc="{0:.2f}".format(round(underc/NmbEvents*100.,2))
    totLegend.AddEntry(CorrCount,"Correctly counted - "+str(corrPerc)+"%","f")
    totLegend.AddEntry(UnderCount,"Undercounted - "+str(underPerc)+"%","f")
    totLegend.AddEntry(OverCount,"Overcounted - "+str(overPerc)+"%","f")
    totLegend.Draw()
    allCountedCanvas.SetLogy()
    CorrCount.GetYaxis().SetTitleOffset(1.3)

    eeCanvas=TCanvas('2eCanvas','Counting of 2-electron events',900,600)
    eeCanvas.Divide(2,1)
    eeCanvas.cd(1)
    CorrCount2.Draw()
    CorrCount2.SetFillColor(ROOT.kBlue)
    CorrCount2.SetStats(0)
    OverCount2.SetFillColor(ROOT.kRed)
    UnderCount2.SetFillColor(ROOT.kGreen-2)
    UnderCount2.Draw("same")
    UnderCount2.SetStats(0)
    OverCount2.Draw("same")
    OverCount2.SetStats(0)
    corrPerc2="{0:.2f}".format(round(corrc2/multEArray[1]*100.,2))#to display percentages in legend
    overPerc2="{0:.2f}".format(round(overc2/multEArray[1]*100.,2))
    underPerc2="{0:.2f}".format(round(underc2/multEArray[1]*100.,2))
    twoLegend.AddEntry(CorrCount2,"Correctly counted - "+str(corrPerc2)+"%","f")
    twoLegend.AddEntry(UnderCount2,"Undercounted - "+str(underPerc2)+"%","f")
    twoLegend.AddEntry(OverCount2,"Overcounted - "+str(overPerc2)+"%","f")
    twoLegend.Draw()
    CorrCount2.GetYaxis().SetTitleOffset(1.3)
    eeCanvas.cd(2)
    SepHistC.SetFillColor(ROOT.kBlue)
    SepHistC.SetFillStyle(3005)
    SepHistU.SetLineColor(ROOT.kGreen-2)
    SepHistU.SetFillColor(ROOT.kGreen-2)
    SepHistU.SetFillStyle(3004)
    SepHistO.SetFillColor(ROOT.kRed)
    SepHistC.Draw()
    SepHistC.SetStats(0)
    SepHistU.Draw("same")
    SepHistU.SetStats(0)
    SepHistO.Draw("same")
    SepHistO.SetStats(0)
    SepHistC.GetYaxis().SetTitleOffset(1.3)
    sepLegend.AddEntry(SepHistC,"Correctly counted","f")
    sepLegend.AddEntry(SepHistU,"Undercounted","f")
    sepLegend.AddEntry(SepHistO,"Overcounted","f")
    sepLegend.Draw()

    efficiencyCanvas=TCanvas('efficiencyCanvas','Efficiency of counting 2-electron events',900,600)
    sepEfficiency.Divide(sep2Hist)
    sepEfficiency.Draw()
    sepEfficiency.SetStats(0)
    sepEfficiency.SetFillColor(ROOT.kBlue)    

countedCanvas=TCanvas('countedCanvas','Distribution of Counted Electron Events',800,700)
multELegend.AddEntry(multEHist[0],"0 electrons","f")
multEHist[1].Draw()#want biggest hist in the back
multEHist[1].SetStats(0)
multEHist[1].SetFillColor(ROOT.kGreen+3)
multELegend.AddEntry(multEHist[1],"1 electron","f")
elec=2
for i in range(2,7):
    multEHist[i].Draw("same")
    multEHist[i].SetStats(0)
    multEHist[i].SetFillColor(2*elec-2)
    multELegend.AddEntry(multEHist[i],str(elec)+" electrons","f")
    elec+=1
multEHist[0].Draw("same")
multEHist[0].SetStats(0)
multEHist[0].SetFillColor(ROOT.kOrange)
multELegend.Draw()
countedCanvas.SetLogy()
multEHist[1].GetYaxis().SetTitleOffset(1.3)
