import ROOT
from ROOT import TCanvas, TH1D, TGraph, TFile, TTree, TH2D, TBranch, TLine, TMath, TPaveText, TLegend, TStyle
import pickle
import numpy as np
from array import array
import math

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Choose what orientation to compare data to simulation
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
siFirst=False #If False, then Tungsten layer first
if not siFirst:
    fluctuation=True #If True, use simulation file that partially accounts for fluctuations. If False, use simulation file that has no fluctuation correction

##~~~~~~~~~~~~~~~~~~~
# Import data files
##~~~~~~~~~~~~~~~~~~~
if siFirst:
    fileIn=ROOT.TFile('/home/jason/data/amanda/TestBeam/TestBeamROOT/2013_07_26_14_13_43.binv2.root')
    fileIn2=ROOT.TFile('/home/jason/data/amanda/TestBeam/TestBeamROOT/2013_07_26_10_41_28.binv2.root')
    inFile='/home/jason/data/amanda/smear9FinalOuts/DENS23-SiWNi/withPos/shifted/12.1GeV-smear9-SiWNi(0.8725)_shifted23.endataIMP.pkl'
    fileIn2.ls()
    testBeam2=fileIn2.Get("TestBeam")
else:
    fileIn=ROOT.TFile('/home/jason/data/amanda/TestBeam/TestBeamROOT/2013_07_29_16_12_24.binv2.root')
    if fluctuation:
        inFile='/home/jason/data/amanda/smear9FinalOuts/DENS23-WNiSi/12.1GeV-smear9-WNiSi(0.773)_fluctDeposit6.endataIMP.pkl'
    else:
        inFile='/home/jason/data/amanda/smear9FinalOuts/DENS23-WNiSi/12.1GeV-smear9-WNiSi(0.773).endataIMP.pkl'

fileIn.ls()
testBeam=fileIn.Get("TestBeam")
unpicklefile=open(inFile,'r')
matrix=pickle.load(unpicklefile)
unpicklefile.close()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create histograms and legends
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sensorsHist=[0]*9
sensorsHistS=[0]*9
layerLegend=[0]*9
energySum=TH1D("Total","Total Measured Charge per Cleaned or Simulated Electron Events;Total  Measured Charge [fC];Entries",100,0.,10000.)
energySumS=TH1D("Total","Total Measured Charge per Cleaned or Simulated Electron Events;Total Measured Charge [fC];Entries",100,0.,10000.)
if siFirst:
    sensorsHist2=[0]*9
    energySum2=TH1D("Total","Total Measured Charge per Cleaned or Simulated Electron Events;Total Measured Charge [fC];Entries",100,0.,10000.)
    for i in range(9):
        sensorsHist[i]=TH1D("Layer "+str(8-i)+"test1","Total Measured Charge - Layer "+str(9-i)+";Total Measured Charge [fC];Entries",30,0.,3000.)
        sensorsHist2[i]=TH1D("Layer "+str(8-i)+"test2","Total Measured Charge - Layer "+str(9-i)+";Total Measured Charge [fC];Entries",30,0.,3000.)
        sensorsHistS[i]=TH1D("Layer "+str(8-i)+"sim","Total Measured Charge - Layer "+str(9-i)+";Total Measured Charge [fC];Entries",30,0.,3000)
else:
    for i in range(9):
        sensorsHist[i]=TH1D("Layer "+str(i)+"test","Total Measured Charge - Layer "+str(i)+";Total Measured Charge [fC];Entries",30,0.,3000.)
        sensorsHistS[i]=TH1D("Layer "+str(i)+"sim","Total Measured Charge - Layer "+str(i)+";Total Measured Charge [fC];Entries",30,0.,3000.)
for i in range(9):
    sensorsHist[i].SetLabelSize(0.05,"xy")
    sensorsHist[i].SetTitleSize(0.06,"xy")
    layerLegend[i]=TLegend(0.5,0.7,0.75,0.9)
totLegend=TLegend(0.5,0.7,0.75,0.9)

##~~~~~~~~~~~~~~~~~~~~~~
# Set global variables
##~~~~~~~~~~~~~~~~~~~~~~
i=0
m=0
for entry in testBeam:
    if i==0:
        eventCount=entry.eventID-1
        break
if siFirst:
    for entry in testBeam2:
        if i==0:
            eventCount2=entry.eventID-1
            break
totalEnergy=0.0
sensorEnergy=[0.0]*9
eventNo=len(matrix)
pixel=len(matrix[0][0])
nullCount=0
statNum=0.0
statDenom=0.0
statTotal=0.0

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cycle through entries and sum energy deposits for each electron event. When a new electron even begins, fill the histograms and clear all variables
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Test beam data file
print "Counting test beam data file"

for entry in testBeam:
    if entry.eventID!=eventCount:
        statDenom=totalEnergy
        for i in range(9):
            if sensorEnergy[i]==0:
                nullCount+=1
                statDenom+=4.
                if siFirst:
                    statNum+=4.*(9-i)*(9-i)
                else:
                    statNum+=4.*(1+i)*(1+i)
        statTotal=statNum/statDenom
        if siFirst:
            if nullCount!=8 and statTotal>44:
                for i in range(9):
                    sensorsHist[i].Fill(sensorEnergy[i],1/151.06)#normalize all to 100 events
                energySum.Fill(totalEnergy,1/151.06)
        else:
            if nullCount!=8 and statTotal>40:
                for i in range(9):
                    sensorsHist[i].Fill(sensorEnergy[i],1/68.32)
                energySum.Fill(totalEnergy,1/68.32)

        nullCount=0
        statNum=0.0
        statDenom=0.0
        statTotal=0.0
        totalEnergy=0.0
        sensorEnergy=[0.0]*9
        eventCount+=1
    sensorEnergy[int(entry.SensorID)]+=entry.Charge*1000000000000000
    totalEnergy+=entry.Charge*1000000000000000
    if siFirst:
        statNum+=(entry.Charge*1000000000000000)*(9-entry.SensorID)*(9-entry.SensorID)
    else:
        statNum+=(entry.Charge*1000000000000000)*(1+entry.SensorID)*(1+entry.SensorID)
    m+=1
    if m%100000==0:
        print m

if siFirst:
    #Second test beam data file - for Si-First runs only (no two good W-First runs)
    m=0
    totalEnergy=0.0
    nullCount=0
    statNum=0.0
    statDenom=0.0
    statTotal=0.0
    sensorEnergy=[0.0]*9

    print "Counting second test beam data file"
 
    for entry in testBeam2:
        if entry.eventID!=eventCount2:
            statDenom=totalEnergy
            for i in range(9):
                if sensorEnergy[i]==0:
                    nullCount+=1
                    statNum+=4.*(9-i)*(9-i)
                    statDenom+=4.
                    statTotal=statNum/statDenom
            if nullCount!=8 and statTotal>44:
                for i in range(9):
                    sensorsHist2[i].Fill(sensorEnergy[i],1/121.84)
                    energySum2.Fill(totalEnergy,1/121.84)
            nullCount=0
            statNum=0.0
            statDenom=0.0
            statTotal=0.0
            totalEnergy=0.0
            sensorEnergy=[0.0]*9
            eventCount2+=1
        sensorEnergy[int(entry.SensorID)]+=entry.Charge*1000000000000000
        totalEnergy+=entry.Charge*1000000000000000
        statNum+=(entry.Charge*1000000000000000)*(9-entry.SensorID)*(9-entry.SensorID)
        m+=1
        if m%100000==0:
            print m

#Simulation file
m=0
totalEnergy=0.0
sensorEnergy=[0.0]*9

print "Counting simulation file"

#Conversion factor from fC to MeV from fit to test beam data and differs for Si-First and W-First runs
if siFirst:
    norm=1/127.64
    convert=29.575
else:
    norm=1/120.64
    convert=29.39

for i in range(eventNo):
    for j in range(9):
        for k in range(pixel):
            totalEnergy+=convert*matrix[i][j][k]
            sensorEnergy[j]+=convert*matrix[i][j][k]
    for q in range(9):
        sensorsHistS[q].Fill(sensorEnergy[q],norm)
    energySumS.Fill(totalEnergy,norm)
    sensorEnergy=[0.0]*9
    totalEnergy=0.0
    if m%1000==0:
        print m
    m+=1

##~~~~~~~~~~~~~
# Print plots
##~~~~~~~~~~~~~
canvas= TCanvas('canvas', 'Layers High', 1400, 1000)
canvas.Divide(3,3)
for i in range(9):
    canvas.cd(i+1)
    canvas.GetPad(i+1).SetLogy()
    if siFirst:
        sensorsHist[8-i].Draw()
        sensorsHist[8-i].SetLineColor(ROOT.kRed)
        sensorsHist[8-i].SetStats(0)
        sensorsHist2[8-i].Draw("same")
        sensorsHist2[8-i].SetLineColor(ROOT.kOrange+5)
        sensorsHist2[8-i].SetStats(0)
    else:
        sensorsHist[i].Draw()
        sensorsHist[i].SetLineColor(ROOT.kRed)
        sensorsHist[i].SetStats(0)
    sensorsHistS[i].Draw("same")
    sensorsHistS[i].SetLineColor(ROOT.kBlack)
    sensorsHistS[i].SetStats(0)
    sensorsHist[8-i].GetXaxis().SetTitleOffset(0.8)
    sensorsHist[8-i].GetYaxis().SetTitleOffset(0.7)
    canvas.Update()
for a in range(9):
    canvas.cd(a+1)
    if siFirst:
        layerLegend[a].AddEntry(sensorsHist[a],"Beam Test Run 43","l")
        layerLegend[a].AddEntry(sensorsHist2[a],"Beam Test Run 28","l")
    else:
        layerLegend[a].AddEntry(sensorsHist[a],"Beam Test Run 24","l")
    layerLegend[a].AddEntry(sensorsHistS[a],"Simulation","l")
    layerLegend[a].Draw()
canvas.Update()

canvas1=TCanvas('canvas1','Total High',800,700)
canvas1.SetLogy()
energySum.Draw()
energySum.SetLineColor(ROOT.kGreen+2)
energySum.SetStats(0)
if siFirst:
    energySum2.Draw("sames")
    energySum2.SetLineColor(ROOT.kBlue)
    energySum2.SetStats(0)
energySumS.Draw("sames")
energySumS.SetLineColor(ROOT.kMagenta)
energySumS.SetStats(0)
canvas1.Update()
energySumS.GetYaxis().SetRangeUser(0.001,100.)
if siFirst:
    totLegend.AddEntry(energySum,"Beam Test Run 43","l")
    totLegend.AddEntry(energySum2,"Beam Test Run 28","l")
else:
    totLegend.AddEntry(energySum,"Beam Test Run 24","l")
totLegend.AddEntry(energySumS,"Simulation","l")
totLegend.Draw()
