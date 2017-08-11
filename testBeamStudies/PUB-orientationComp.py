import ROOT
from ROOT import TCanvas, TH1D, TGraph, TFile, TTree, TH2D, TBranch, TLine, TMath, TPaveText, TLegend, TStyle
import numpy as np
from array import array
import pickle
import math

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import test beam and simulation data files
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fileIn24=ROOT.TFile('/home/jason/data/amanda/TestBeam/TestBeamROOT/2013_07_29_16_12_24.binv2.root')
fileIn24.ls()
testBeam24=fileIn24.Get("TestBeam")

fileIn43=ROOT.TFile('/home/jason/data/amanda/TestBeam/TestBeamROOT/2013_07_26_14_13_43.binv2.root')
fileIn43.ls()
testBeam43=fileIn43.Get("TestBeam")

fileIn28=ROOT.TFile('/home/jason/data/amanda/TestBeam/TestBeamROOT/2013_07_26_10_41_28.binv2.root')
fileIn28.ls()
testBeam28=fileIn28.Get("TestBeam")

inFileSi='/home/jason/data/amanda/smear9FinalOuts/DENS23-SiWNi/withPos/shifted/12.1GeV-smear9-SiWNi(0.8725)_shifted23.endataIMP.pkl'
unpicklefile=open(inFileSi,'r')
matrixSi=pickle.load(unpicklefile)
unpicklefile.close()

inFileW='/home/jason/data/amanda/smear9FinalOuts/DENS23-WNiSi/12.1GeV-smear9-WNiSi(0.773).endataIMP.pkl'
unpicklefile=open(inFileW,'r')
matrixW=pickle.load(unpicklefile)
unpicklefile.close()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create histograms and legends
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sensorsHist24=[0]*9
sensorsHist43=[0]*9
sensorsHist28=[0]*9
sensorsHistSi=[0]*9
sensorsHistW=[0]*9
totalLegend=TLegend(0.65,0.7,0.95,0.9)
layerLegend=TLegend(0.65,0.7,0.95,0.9)
totalLegendSim=TLegend(0.65,0.7,0.95,0.9)
layerLegendSim=TLegend(0.65,0.7,0.95,0.9)
energySum24=TH1D("Total24","Total Measured Charge (All layers of Si-First, first 8 layers of W-First); Total Measured Charge [fC]; Entries (normalized to 100)",100,0.,10000.)
energySum43=TH1D("Total43","Total Measured Charge [fC]",100,0.,10000.)
energySum28=TH1D("Total28","Total Measured Charge [fC]",100,0.,10000.)
energySumW=TH1D("TotalSi","Total Measured Charge (All layers of Si-First, first 8 layers of W-First Simulation); Total Measured Charge [fC]; Entries (normalized to 100)",100,0.,10000.)
energySumSi=TH1D("TotalW","Total Measured Charge [fC]",100,0.,10000.)
for i in range(9):
    sensorsHist43[i]=TH1D("Layer "+str(i),"Total Measured Charge Deposits Per Layer; Total Measured Charge [fC]; Entries (normalized to 100)",30,0.,3000.)
    sensorsHist24[i]=TH1D("Layer "+str(i)+"24","Sensor Energy Deposits - Layer "+str(i),30,0.,3000.)
    sensorsHist28[i]=TH1D("Layer "+str(i)+"28","Sensor Energy Deposits - Layer "+str(i),30,0.,3000.)
    sensorsHistW[i]=TH1D("Layer "+str(i)+"W","Total Measured Charge Deposits Per Layer; Total Measured Charge [fC]; Entries (normalized to 100)",30,0.,3000.)
    sensorsHistSi[i]=TH1D("Layer "+str(i)+"Si","Sensor Energy Deposits - Layer "+str(i),30,0.,3000.)
sensorsHist43[0].SetTitle("Total Measured Charge Deposits Per Layer - Upstream Layer")
sensorsHist43[8].SetTitle("Total Measured Charge Deposits Per Layer - Downstream Layer")

##~~~~~~~~~~~~~~~~~~~~~~~~
# Set universal variables
##~~~~~~~~~~~~~~~~~~~~~~~~
i=0
m=0
for entry in testBeam24:
    if i==0:
        eventCount24=entry.eventID-1
        break
for entry in testBeam43:
    if i==0:
        eventCount43=entry.eventID-1
        break
for entry in testBeam28:
    if i==0:
        eventCount28=entry.eventID-1
        break
eventNoSi=len(matrixSi)
eventNoW=len(matrixW)
pixelSi=len(matrixSi[0][0])
pixelW=len(matrixW[0][0])
totalEnergy=0.0
sensorEnergy=[0.0]*9
statNum=0.0
statDenom=0.0
statTotal=0.0
nullCount=0


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cycle through entries and sum energy deposits for each electron event. When a new electron even begins, fill the histograms and clear all variables
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# W-First Test Beam 24 (Upstream layer = 0)
print "Counting Test Beam 24"

tempTot=0.0 #is identical to "totalEnergy" but only sums over first 8 layers
for entry in testBeam24:
    if entry.eventID!=eventCount24:
        statDenom=totalEnergy
        for i in range(9):
            if sensorEnergy[i]==0:
                nullCount+=1
                statNum+=4.*(1+i)*(1+i)
                statDenom+=4.
        statTotal=statNum/statDenom
        if nullCount!=8 and statTotal>40: 
            for i in range(8):
                sensorsHist24[i].Fill(sensorEnergy[i],1/68.32)
                tempTot+=sensorEnergy[i]
            energySum24.Fill(tempTot,1/68.32)
        totalEnergy=0.0
        tempTot=0.0
        sensorEnergy=[0.0]*9
        nullCount=0
        statNum=0.0
        statDenom=0.0
        statTotal=0.0
        eventCount24+=1
    sensorEnergy[int(entry.SensorID)]+=entry.Charge*1000000000000000
    totalEnergy+=entry.Charge*1000000000000000
    statNum+=(entry.Charge*1000000000000000)*(1+entry.SensorID)*(1+entry.SensorID)
    m+=1
    if m%100000==0:
        print m

# Si-First Test Beam 43 (Upstream layer = 8)
m=0
i=0
totalEnergy=0.0
sensorEnergy=[0.0]*9
nullCount=0
statNum=0.0
statDenom=0.0
statTotal=0.0

print "Counting Test Beam 43"

for entry in testBeam43:
    if entry.eventID!=eventCount43:
        statDenom=totalEnergy
        for i in range(9):
            if sensorEnergy[i]==0:
                nullCount+=1
                statNum+=4.*(9-i)*(9-i)
                statDenom+=4.
        statTotal=statNum/statDenom
        if nullCount!=8 and statTotal>44:
            for i in range(9):
                sensorsHist43[i].Fill(sensorEnergy[8-i],1/151.06)
            energySum43.Fill(totalEnergy,1/151.06)
        totalEnergy=0.0
        sensorEnergy=[0.0]*9
        nullCount=0
        statNum=0
        statDenom=0
        statTotal=0
        eventCount43+=1
    sensorEnergy[int(entry.SensorID)]+=entry.Charge*1000000000000000
    totalEnergy+=entry.Charge*1000000000000000
    statNum+=(entry.Charge*1000000000000000)*(9-entry.SensorID)*(9-entry.SensorID)
    m+=1
    if m%100000==0:
        print m

# Si-First Test Beam 28 (Upstream layer = 8)
m=0
i=0
totalEnergy=0.0
sensorEnergy=[0.0]*9
nullCount=0
statNum=0.0
statDenom=0.0
statTotal=0.0

print "Counting Test Beam 28"

for entry in testBeam28:
    if entry.eventID!=eventCount28:
        statDenom=totalEnergy
        for i in range(9):
            if sensorEnergy[i]==0:
                nullCount+=1
                statNum+=4.*(9-i)*(9-i)
                statDenom+=4.
        statTotal=statNum/statDenom
        if nullCount!=8 and statTotal>44:
            for i in range(9):
                sensorsHist28[i].Fill(sensorEnergy[8-i],1/121.84)
            energySum28.Fill(totalEnergy,1/121.84)
        totalEnergy=0.0
        sensorEnergy=[0.0]*9
        nullCount=0
        statNum=0.0
        statDenom=0.0
        statTotal=0.0
        eventCount28+=1
    sensorEnergy[int(entry.SensorID)]+=entry.Charge*1000000000000000
    totalEnergy+=entry.Charge*1000000000000000
    statNum+=(entry.Charge*1000000000000000)*(9-entry.SensorID)*(9-entry.SensorID)
    m+=1
    if m%100000==0:
        print m

# W-First Simulation (Upstream layer = 0)
count=0
totalEnergy=0.0
sensorEnergy=[0.0]*9
                                                                                          
print "Counting W-first simulation"

for i in range(eventNoW):
    for j in range(8):#only use first 8 layers to compare                                      
        for k in range(pixelW):
            sensorEnergy[j]+=matrixW[i][j][k]*29.39 #conversion factor from MeV to fC
            totalEnergy+=matrixW[i][j][k]*29.39
    energySumW.Fill(totalEnergy)
    for q in range(9):
        sensorsHistW[q].Fill(sensorEnergy[q])
    sensorEnergy=[0.0]*9
    totalEnergy=0.0
    if count%1000==0:
        print count
    count+=1

# Si-First Simulation (Upstream layer = 0)                                                                                                                                        
count=0
totalEnergy=0.0
sensorEnergy=[0.0]*9

print "Counting Si-first simulation"

for i in range(eventNoSi):
    for j in range(9):
        for k in range(pixelSi):
            totalEnergy+=matrixSi[i][j][k]*29.575
            sensorEnergy[j]+=matrixSi[i][j][k]*29.575
    energySumSi.Fill(totalEnergy)
    for q in range(9):
        sensorsHistSi[q].Fill(sensorEnergy[q])
    sensorEnergy=[0.0]*9
    totalEnergy=0.0
    if count%1000==0:
        print count
    count+=1


##~~~~~~~~~~~~~~~~  
# Print plots
##~~~~~~~~~~~~~~~
layerCanvas= TCanvas('layerCanvas', 'Layer Comparison', 1500,1100)
layerCanvas.Divide(3,3)
for i in range(9):#plot all layers for Si-First runs
    layerCanvas.cd(i+1)
    layerCanvas.GetPad(i+1).SetLogy()
    sensorsHist43[i].Draw()
    sensorsHist43[i].SetLineColor(ROOT.kGreen)
    sensorsHist43[i].SetStats(0)
    sensorsHist28[i].Draw("same")
    sensorsHist28[i].SetLineColor(ROOT.kBlue)
    sensorsHist28[i].SetStats(0)
for j in range(8):#plot only first 8 layers of W-First runs
    layerCanvas.cd(j+2)
    sensorsHist24[j].Draw("same")
    sensorsHist24[j].SetLineColor(ROOT.kRed)
    sensorsHist24[j].SetStats(0)
layerCanvas.Update()
layerCanvas.cd(2)
layerLegend.AddEntry(sensorsHist24[1],"W-First Test Beam 24","l")
layerLegend.AddEntry(sensorsHist43[1],"Si-First Test Beam 43","l")
layerLegend.AddEntry(sensorsHist28[1],"Si-First Test Beam 28","l")
layerLegend.Draw()   

totalCanvas=TCanvas('totalCanvas','Total Charge Comparison',800,700)
totalCanvas.SetLogy()
energySum24.Draw()
energySum24.SetLineColor(ROOT.kRed)
energySum24.SetStats(0)
energySum43.Draw("same")
energySum43.SetLineColor(ROOT.kGreen)
energySum43.SetStats(0)
energySum28.Draw("same")
energySum28.SetLineColor(ROOT.kBlue)
energySum28.SetStats(0)
totalCanvas.Update()
totalLegend.AddEntry(energySum24,"W-First Test Beam 24","l")
totalLegend.AddEntry(energySum43,"Si-First Test Beam 43","l")
totalLegend.AddEntry(energySum28,"Si-First Test Beam 28","l")
totalLegend.Draw()

layerCanvasSim= TCanvas('layerCanvasSim', 'Simulation Layer Comparison', 1500,1100)
layerCanvasSim.Divide(3,3)
for i in range(9):#plot all layers for Si-First runs
    layerCanvasSim.cd(i+1)
    layerCanvasSim.GetPad(i+1).SetLogy()
    sensorsHistSi[i].Draw()
    sensorsHistSi[i].SetLineColor(ROOT.kGreen)
    sensorsHistSi[i].SetStats(0)
for j in range(8):#plot only first 8 layers of W-First runs
    layerCanvasSim.cd(j+2)
    sensorsHistW[j].Draw("same")
    sensorsHistW[j].SetLineColor(ROOT.kRed)
    sensorsHistW[j].SetStats(0)
layerCanvasSim.Update()
layerCanvasSim.cd(2)
layerLegendSim.AddEntry(sensorsHistW[1],"W-First Simulation","l")
layerLegendSim.AddEntry(sensorsHistSi[1],"Si-First Simulation","l")
layerLegendSim.Draw()

totalCanvasSim=TCanvas('totalCanvasSim','Total Charge Comparison of Simulation',800,700)
totalCanvasSim.SetLogy()
energySumW.Draw()
energySumW.SetLineColor(ROOT.kRed)
energySumW.SetStats(0)
energySumSi.Draw("same")
energySumSi.SetLineColor(ROOT.kGreen)
energySumSi.SetStats(0)
totalCanvasSim.Update()
totalLegendSim.AddEntry(energySumW,"W-First Simulation","l")
totalLegendSim.AddEntry(energySumSi,"Si-First Simulation","l")
totalLegendSim.Draw()
