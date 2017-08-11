import ROOT
from ROOT import *
import pickle
import numpy as np
from array import array
import math

##~~~~~~~~~~~~~~~~~~~~~~~~
# Indicate orientation
##~~~~~~~~~~~~~~~~~~~~~~~~
siFirst=False #if False, then tungsten layer is first

##~~~~~~~~~~~~~~~~~~~~~~~~
# Import .pkl  data file
##~~~~~~~~~~~~~~~~~~~~~~~~
if siFirst:
    inFile='/home/jason/data/amanda/smear9FinalOuts/DENS23-SiWNi/withPos/shifted/12.1GeV-smear9-SiWNi(0.8725)_shifted23.endataIMP.pkl'
else:
    inFile='/home/jason/data/amanda/smear9FinalOuts/DENS23-WNiSi/12.1GeV-smear9-WNiSi(0.773).endataIMP.pkl'

unpicklefile=open(inFile,'r')
matrix=pickle.load(unpicklefile)
unpicklefile.close()

eventNo= len(matrix) 
pixel= len(matrix[0][0])

##~~~~~~~~~~~~~~~~~~~
# Create histograms
##~~~~~~~~~~~~~~~~~~~
multEHist=[0,0,0,0,0]    
energySumHist=TH1D("Total","Total Measured Charge per Simulated Electron Event;Total Measured Charge [fC];Entries",100,0,10000)
singleEHist=TH1D("Total 1 electron","Total Measured Charge Per Single-Electron Simulated Event; Total Measured Energy [fC];Entries",50,0,3000)
for i in range(5):
    multEHist[i]=TH1D(str(i)+" Electrons","Electron Events - Simulation Truth; Total Measured Charge [fC]; Entries (normalized)",100,0,10000)
multELegend=TLegend(0.65,0.5,0.85,0.7)

##~~~~~~~~~~~~~~~~~~~
# Set variables
##~~~~~~~~~~~~~~~~~~~
count=0
totalEnergy=0.0
sensorEnergy=[0.0]*9
if siFirst:
    multEArray=[8000,3490,1015,221,38] #SiWNi
else:
    multEArray=[8000,3092,796,153,23]#WNiSi

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cycle through entries in matrix and sum energy deposits in each layer for each electron event. When a new electron event begins, fill the histograms and clear all variables
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for i in range(eventNo):
    for j in range(9):
        for k in range(pixel):
            sensorEnergy[j]+=matrix[i][j][k]
            totalEnergy+=matrix[i][j][k]*29.57 #Factor of acts as a conversion factor from the simulation's MeV to fC, found by fitting the simulation data
    if i<multEArray[0]:
        multEHist[0].Fill(totalEnergy,1/127.64) #Normalize with respect to 12764 initial events
        singleEHist.Fill(totalEnergy)
    elif i< multEArray[0]+multEArray[1]:
        multEHist[1].Fill(totalEnergy,1/127.64)
    elif i<multEArray[0]+multEArray[1]+multEArray[2]:
        multEHist[2].Fill(totalEnergy,1/127.64)
    elif i<multEArray[0]+multEArray[1]+multEArray[2]+multEArray[3]:
        multEHist[3].Fill(totalEnergy,1/127.64)
    else:
        multEHist[4].Fill(totalEnergy,1/127.64)
    #fill histograms with each event's total energy after summing over all deposits
    energySumHist.Fill(totalEnergy)
    sensorEnergy=[0.0]*9
    totalEnergy=0.0
    #count throguh loop to display progress
    if count%4000==0:
        print count
    count+=1

##~~~~~~~~~~~~~~~~
# Print plots
##~~~~~~~~~~~~~~~~

canvas=TCanvas('canvas1','Total',800,700)
canvas.SetLogy()
energySumHist.Draw()
energySumHist.SetStats(0)
energySumHist.SetFillColor(ROOT.kBlue)

elec=1
multECanvas=TCanvas('multECanvas','Totals Separated by Number of Electrons per Event',800,700)
multECanvas.SetLogy()
for i in range(5):
    if i==0:
        multEHist[0].Draw()
        multEHist[0].SetStats(0)
        multEHist[0].SetFillColor(ROOT.kGreen+3)
        multELegend.AddEntry(multEHist[0],str(elec)+" electrons","f")
        elec+=1
    else:
        multEHist[i].Draw("same")
        multEHist[i].SetStats(0)
        multEHist[i].SetFillColor(2*elec-2)
        multELegend.AddEntry(multEHist[i],str(elec)+" electrons","f")
        elec+=1
multELegend.Draw()

singleECanvas=TCanvas('singleECanvas','Total - Single Electron Events',800,700)
singleEHist.Fit("gaus")
singleEHist.Draw()
singleEHist.SetFillColor(ROOT.kBlue)
singleEHist.SetStats(0)
singleEHist.GetYaxis().SetTitleOffset(1.2)

