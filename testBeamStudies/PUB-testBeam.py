import ROOT
from ROOT import *
import numpy as np
from array import array
import math

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Input orientation of test beam
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
siFirst=True
#if set to false, then W-first is run


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import root data file
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if siFirst:
    fileIn=ROOT.TFile('/home/jason/data/amanda/TestBeam/TestBeamROOT/2013_07_26_14_13_43.binv2.root')
#    fileIn=ROOT.TFile('/home/jason/data/amanda/TestBeam/TestBeamROOT/2013_07_26_10_41_28.binv2.root')
else:
    fileIn=ROOT.TFile('/home/jason/data/amanda/TestBeam/TestBeamROOT/2013_07_29_16_12_24.binv2.root')

fileIn.ls()
testBeam=fileIn.Get("TestBeam")

##~~~~~~~~~~~~~~~~~~~~~~~~
# Create histograms
##~~~~~~~~~~~~~~~~~~~~~~~~
sensorsHist=[0]*9
nullHist=[0]*9
null8Hist=[0]*9
statHist=[0]*9

energySum=TH1D("Total","Total Measured Charge per Electron Event;Total Measured Charge [fC]; Entries",100,0,10000)
positionSum=TH2D("Position","Transverse Distribution - Sum of all Hits;x [mm]; y[mm]; Measured Charge [fC]",35,-100,100,35,-100,100)
for i in range(9):
    nullHist[i]=TH1D("Null plots "+str(i),"Events with "+str(i)+" empty layers; Total Measured Charge [fC]; Entries",50,0,10000)
    nullHist[i].SetLabelSize(0.05,"xy")
    nullHist[i].SetTitleSize(0.06,"xy")
    nullHist[i].SetTitleSize(0.2,"t")
    null8Hist[i]=TH1D("8 null layers "+str(i),"Energy in Only Layer "+str(i)+"; Total Measured Charge [fC]; Entries",50,0,500)
    null8Hist[i].SetLabelSize(0.05,"xy")
    null8Hist[i].SetTitleSize(0.06,"xy")
    null8Hist[i].SetTitleSize(0.1,"")
nullHist[6].GetYaxis().SetTitle("Entries (log scale)")
nullHist[7].GetYaxis().SetTitle("Entries (log scale)")
nullHist[8].GetYaxis().SetTitle("Entries (log scale)")
for i in range(9):
    statHist[i]=TH1D("Weighted"+str(i),str(i)+" Empty Layers; Weighted Energy [fC]; Entries",50,0,100)
    statHist[i].SetLabelSize(0.05,"xy")
    statHist[i].SetTitleSize(0.06,"xy")
energySumCleaned=TH1D("Cleaned Total","Total measured Charge After Cleaning; Total Measured Charge [fC];Entries",100,0,10000)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set universal variables
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
entry0=0
eventCount=0
eventCount0=0
for entry in testBeam:
    if entry0==0:
        eventCount=entry.eventID-1 #subtract 1 so that the first event is included in the following loop
        eventCount0=entry.eventID-1
        break
totalEnergy=0.0
sensorEnergy=[0.0]*9
nullCount=0
statNum=0.0
statDenom=0.0

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cycle through entries and sum energy deposits for each electron event. When a new electron even begins, fill the histograms and clear all variables
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for entry in testBeam:
    #move to next event after all hits from the previous event have been counted up
    if entry.eventID!=eventCount:
        energySum.Fill(totalEnergy)
        statDenom=totalEnergy
        for i in range(9):
            if sensorEnergy[i]==0: #insert MIP in empty layers
                nullCount+=1
                if siFirst:
                    statNum+=4.*(9-i)*(9-i) 
                else:
                    statNum+=4.*(i+1)*(i+1)
                statDenom+=4.
        for i in range(9):
            if nullCount==i:
                nullHist[i].Fill(totalEnergy)
                statHist[i].Fill(statNum/statDenom)
        if nullCount==8:
            for i in range(9):
                if sensorEnergy[i]!=0:#record hits from only nonzero layer
                    null8Hist[i].Fill(totalEnergy)
        if nullCount!=8: #Events with deposits in at least 2 layers and with a statistics value of >44 are not contamination and retained in cleaned data version
            if siFirst and (statNum/statDenom)>44:
                energySumCleaned.Fill(totalEnergy)
            elif not siFirst and (statNum/statDenom)>40:
                energySumCleaned.Fill(totalEnergy)
        nullCount=0
        statNum=0.0
        statDenom=0.0
        totalEnergy=0.0
        sensorEnergy=[0.0]*9
        eventCount+=1
        #count through loop to display progress
        if (eventCount-eventCount0)%1000==0:
            print eventCount-eventCount0
    #sum energies in hits to get totals for the whole event
    sensorEnergy[int(entry.SensorID)]+=entry.Charge*1000000000000000
    totalEnergy+=entry.Charge*1000000000000000
    if siFirst:
        statNum+=(entry.Charge*1000000000000000)*(9-entry.SensorID)*(9-entry.SensorID)
    else:
        statNum+=(entry.Charge*1000000000000000)*(entry.SensorID+1)*(entry.SensorID+1)
    #fill transverse distribution histogram with all spatial positions of all hits
    positionSum.Fill(entry.XPosition,entry.YPosition,entry.Charge)

##~~~~~~~~~~~~~~~~~~
# Print plots
##~~~~~~~~~~~~~~~~~~
totCanvas=TCanvas('totCanvas','Total',800,700)
totCanvas.SetLogy()
energySum.Draw()
energySum.SetFillColor(ROOT.kBlue)
energySum.SetStats(0)

distCanvas=TCanvas('distCanvas','Total Transverse Distribution',700,700)
ROOT.gStyle.SetPalette(56)
distCanvas.GetPad(0).SetRightMargin(2.)
positionSum.Draw("COLZ")
positionSum.SetStats(0)
positionSum.GetYaxis().SetTitleOffset(1.2)
positionSum.GetZaxis().SetTitleOffset(1.2)

nullCanvas=TCanvas('nullCanvas','Events with Empty Layers',1500,1100)
nullCanvas.Divide(3,3)
nullCanvas.GetPad(7).SetLogy()
nullCanvas.GetPad(8).SetLogy()
nullCanvas.GetPad(9).SetLogy()
for i in range(9):
    nullCanvas.cd(i+1)
    nullHist[i].Draw()
    nullHist[i].SetStats(0)
    nullHist[i].SetFillColor(ROOT.kRed)
    nullHist[i].GetYaxis().SetTitleOffset(0.8)
    nullHist[i].GetXaxis().SetTitleOffset(0.8)

null8Canvas=TCanvas('null8Canvas','Events with 8 Empty Layers',1500,1100)
null8Canvas.Divide(3,3)
for i in range(9):
    null8Canvas.cd(i+1)
    if siFirst:
        index=8-i
    else:
        index=i
    null8Hist[index].Draw()
    null8Hist[index].SetStats(0)
    null8Hist[index].SetFillColor(ROOT.kMagenta)
    null8Hist[index].GetYaxis().SetTitleOffset(0.9)
    null8Hist[index].GetXaxis().SetTitleOffset(0.8)


statCanvas=TCanvas('stat1Canvas','Statistics Ratio R',1500,1100)
statCanvas.Divide(3,3)
for i in range(9):
    statCanvas.cd(i+1)
    statHist[i].Draw()
    statHist[i].SetStats(0)
    statHist[i].SetFillColor(ROOT.kGreen-3)
    statHist[i].GetYaxis().SetTitleOffset(0.8)
    statHist[i].GetXaxis().SetTitleOffset(0.8)

cleanedCanvas=TCanvas('cleanedCanvas','Cleaned Total',800,700)
cleanedCanvas.SetLogy()
energySumCleaned.Draw()
energySumCleaned.SetStats(0)
energySumCleaned.SetFillColor(ROOT.kBlue)
