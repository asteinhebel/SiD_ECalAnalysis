from pyLCIO import IOIMPL, EVENT,UTIL
from pyLCIO.io.LcioReader import LcioReader
from array import array
from ROOT import *
import math

##~~~~~~~~~~~~~~~~~~~~~~~~~
#Define helpful methods
##~~~~~~~~~~~~~~~~~~~~~~~~~
def solidAngle(hit):
    angle=(math.acos((math.cos(phi*math.pi/180.)*hit.getPosition()[0]+math.sin(phi*math.pi/180.)*hit.getPosition()[1])/math.sqrt(hit.getPosition()[0]*hit.getPosition()[0]+hit.getPosition()[1]*hit.getPosition()[1]+hit.getPosition()[2]*hit.getPosition()[2])))
    if hit.getPosition()[1]<0:
        angle=-angle 
    return angle

##~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define global variables
##~~~~~~~~~~~~~~~~~~~~~~~~~
phi=0
phiRad=phi*math.pi/180.
evtsPerEn=5000
inFile='reco_'+str(evtsPerEn)+'a.GeVscan.'+str(phi)+'phi.slcio'
nmbEvents=7*evtsPerEn

#create a reader
readerL = LcioReader(inFile )

#useful variables
energyArray=[1,2,5,10,20,50,100]
energyArray_offset=[x*1.05 for x in energyArray] #for nicer looking plot
ecalLayers=31
hcalLayers=40
#input values from Julia toy calibration code
calibrationLinear=58.464354462350954
calibrationNonlinear=0.5941123907726221

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create histograms, legends, and graphs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hitEnergyHist=[0]*7 #full distribution of all 5000 events before correction
hitEnergyHist1=[0]*7 #full distribution of all 5000 events after shower shape correction
hitEnergyHistCorr=[0]*7 #full distribution of all 5000 events after correction
evtGraph=[0]*7 #graph the shower profile in ECal and HCal for each event
enRadLenGraph=[0]*7
enResLegend1=TLegend(0.6,0.75,0.8,0.9)
enResLegend2=TLegend(0.6,0.75,0.8,0.9)

for i in range(7):
    hitEnergyHist[i] = TH1D( 'TotalDepositedEnergy'+str(i), 'Event Energy Deposit ('+str(energyArray[i])+' GeV photons, phi='+str(phi)+', theta=90);Total Measured Deposited Energy [# MIPs];Entries',500,0.,20000.)
    hitEnergyHist1[i] = TH1D( 'TotalDepositedEnergyCalib'+str(i), 'Event Energy Deposit ('+str(energyArray[i])+' GeV photons, phi='+str(phi)+', theta=90);Total Measured Deposited Energy [# MIPs];Entries',100000,0.,400000000.)
    hitEnergyHistCorr[i] = TH1D( 'TotalDepositedEnergyCorr'+str(i), 'Corrected Event Energy Deposit ('+str(energyArray[i])+' GeV photons, phi='+str(phi)+', theta=90) ;Total Measured Deposited Energy [# MIPs];Entries', 500,0.,20000.)
    evtGraph[i]=[0]*evtsPerEn

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Decode Layer Information
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~
layerNos=[0]*nmbEvents
for i in range(nmbEvents):
    layerNos[i]=[]
print "Layer arrays created"

eventNo=-1
# loop over the events
for event in readerL:
    eventNo+=1
    if eventNo%100==0:
        print "recording layers of event "+str(eventNo)
    # get a hit collection
    ecalHits = event.getCollection( 'ECalBarrelHits' )
    # get the cell ID encoding string from the collection parameters
    cellIdEncoding = ecalHits.getParameters().getStringVal( EVENT.LCIO.CellIDEncoding )
    # define a cell ID decoder for the collection
    idDecoder = UTIL.BitField64( cellIdEncoding )
    # loop over all hits in the collection
    for caloHit in ecalHits:
        # combine the two 32 bit cell IDs of the hit into one 64 bit integer
        cellID = long( caloHit.getCellID0() & 0xffffffff ) | ( long( caloHit.getCellID1() ) << 32 )
        # set up the ID decoder for this cell ID
        idDecoder.setValue( cellID )
        # access the field information using a valid field from the cell ID encoding string
        layerNos[eventNo].append(idDecoder['layer'].value())
        
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a reader, open an LCIO file, and initialize arrays for data/analysis
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open( inFile )
count=-1
evNo=0 #tracks what energy the event has for histogram filling loop
eventLayerArray=[0]*7 #record total deposits in each layer of each event of ECal (weighted for sampling fraction correction)
#eventLayerArray2=[0]*7 #record total deposits in each layer of each event of ECal (weighted for sampling fraction correction) after shower shape calibration
corrTotalArray=[0]*7 #record mean of each event after leakage correction
evtFitFnArray=[0]*7 #create Gaussian function to fit to each event's profile
for i in range(7):
    corrTotalArray[i]=[0.]*evtsPerEn
    evtFitFnArray[i]=[0.]*evtsPerEn
    eventLayerArray[i]=[0.]*evtsPerEn

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loop over all events in file
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for dep in reader:
    count+=1
    if count%100==0:
	print "Summing energy of event "+str(count)
    hitTotal=0.
    hitTotal2=0.	
    layerTotal=[0.]*31
    # get the collection from the event
    hitCollection = dep.getCollection( 'ECalBarrelHits' )     
    # loop over all hits in the ecal collection and record total deposits for each layer and whole (weighted) event after backscatter rejection
    hitNo=0
    for hit in hitCollection:
        if solidAngle(hit)<phiRad+0.2 and solidAngle(hit)>phiRad-0.2: #hits within cone of 0.2 radians of original particle gun trajectory (rejects backscatter)
            if layerNos[count][hitNo]>0 and layerNos[count][hitNo]<21: #exclude tracking layer0
                hitTotal+=hit.getEnergy()
                layerTotal[layerNos[count][hitNo]]+=hit.getEnergy()
            elif layerNos[count][hitNo]>0:#exclude tracking layer0
                hitTotal+=hit.getEnergy()*2
                layerTotal[layerNos[count][hitNo]]+=hit.getEnergy()
            hitNo+=1
 
    # fill arrays with deposit information in units of MIPs (in silicon 1 MIP=0.124 MeV)
    if count<evtsPerEn*(evNo+1):
        hitEnergyHist[evNo].Fill(hitTotal/0.000124)
        hitEnergyHist1[evNo].Fill(hitTotal/0.000124*(calibrationLinear+calibrationNonlinear*math.sqrt(hitTotal)))
        eventLayerArray[evNo][count-evNo*evtsPerEn]=[x/0.000124 for x in layerTotal] 
	corrTotalArray[evNo][count-evNo*evtsPerEn]=hitTotal/0.000124
    if count==(evtsPerEn-1)*(evNo+1)+evNo:
        evNo+=1

reader.close()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fill Arrays for Graphs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
layerMeanArray=[0]*7 #find mean deposit for each ECal layer
radLengthArray=[0.]*ecalLayers #position of each ECal layer in radiation lengths

for i in range(7):
    layerMeanArray[i]=[0.]*ecalLayers
    for m in range(ecalLayers):
        for n in range(evtsPerEn):
            layerMeanArray[i][m]+=eventLayerArray[i][n][m]/float(evtsPerEn) #mean of all deposits for each layer

for z in range(ecalLayers+hcalLayers):
    if z==0:
        radLengthArray[z]=0.
    elif z<21:
        radLengthArray[z]=radLengthArray[z-1]+(2.5/math.cos(math.pi/180.*phi)*0.26) #X0=3.85 mm of tungsten alloy
    elif z<ecalLayers:
        radLengthArray[z]=radLengthArray[z-1]+(5./math.cos(math.pi/180.*phi)*0.26)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fit Distributions and Estimate Leakage for each Event
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Fit full distributions (depth vs measured energy) with Gamma Distribution to get starting point for parameters
fullDistFitFnArray=[0]*7
params=[0.]*7

for i in range(7):
    enRadLenGraph[i]=TGraph(31,array('d',radLengthArray),array('d',layerMeanArray[i]))
    fullDistFitFnArray[i]=TF1("fullDistFitFnArray"+str(i),"[0]*TMath::Power(x,[1])*TMath::Exp(-[2]*x)",radLengthArray[0],70)
    params[i]=[0.]*3 #for parameters [0], [1], and [2] at each energy
    #a few first guesses to convince fits to converge
    if i==6:
        fullDistFitFnArray[i].SetParameter(1,3.)
    elif i==3 or i==4:
        fullDistFitFnArray[i].SetParameter(0,2.)
        fullDistFitFnArray[i].SetParameter(1,3.)
    else:
        fullDistFitFnArray[i].SetParameter(1,4.)
    enRadLenGraph[i].Fit(fullDistFitFnArray[i],"R")
    for k in range(3):
        params[i][k]=fullDistFitFnArray[i].GetParameter(k)

valid=[0]*7
#Fit each event independently with a Gamma Distribution, using full distribution parameters as first guess of fit parameters
xInt=(-1/0.3)*math.log(1e-11) #where y=10e-11 assuming exp(-0.3x)
for i in range(7):
    for z in range(evtsPerEn):
        evtGraph[i][z]=TGraph(31,array('d',radLengthArray),array('d',eventLayerArray[i][z]))
        evtFitFnArray[i][z]=TF1("evtFit"+str(z)+str(i),"[0]*TMath::Power(x,[1])*TMath::Exp(-[2]*x)",radLengthArray[0],radLengthArray[30])
        evtFitFnArray[i][z].SetParameters(params[i][0],params[i][1],params[i][2])
        evtGraph[i][z].Fit(evtFitFnArray[i][z],"R")
        fitPtr=evtGraph[i][z].Fit(evtFitFnArray[i][z],"S")
        if fitPtr.IsValid():#confirm that fit converges
            corrTotalArray[i][z]+=evtFitFnArray[i][z].Integral(radLengthArray[30],xInt,1.e-12)*50./57. #integrate up to x~84, convert HCal MIPS to GeV
            hitEnergyHistCorr[i].Fill(corrTotalArray[i][z])
            #if fit doesn't converge, no leakage correction is added and original measured energy sum remains
            valid[i]+=1./evtsPerEn*100.
print "Percentage of converging fits for each energy: "+str(valid)
 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# Calculate energy resolution with mean and RMS of the full distribution before and after leakage correction
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#In these 2d arrays of 7x2, 0 holds fit info from raw event distribution and 1 holds fit info from leakage corrected event distribution
fitFunctionGauss=[0]*7
scaledEnRes=[0]*7
error=[0.]*7
zeros=[0.]*7

for i in range(7):
    fitFunctionGauss[i]=[0.]*3
    scaledEnRes[i]=[0.]*3
    error[i]=[0.]*3
    for j in range(3):
        #fit a Gaussian to the total distribution (measured energy vs entries) to find mean and RMS
        fitFunctionGauss[i][j]=TF1('fitFnGauss'+str(i),'gaus',0.,400000000.)
        fitFunctionGauss[i][j].SetLineStyle( kDashed )
        if j==0:
            fitFunctionGauss[i][0].SetLineColor( kBlue )
            hitEnergyHist[i].Fit( fitFunctionGauss[i][0],'R')
        elif j==1:
            fitFunctionGauss[i][1].SetLineColor( kRed )
            hitEnergyHist1[i].Fit( fitFunctionGauss[i][1],'R')
        else:
            fitFunctionGauss[i][2].SetLineColor( kRed )
            hitEnergyHistCorr[i].Fit( fitFunctionGauss[i][2],'R')
        #calculate scaled energy resolution and errors
        scaledEnRes[i][j]=math.sqrt(float(energyArray[i]))*fitFunctionGauss[i][j].GetParameter(2)/fitFunctionGauss[i][j].GetParameter(1)
        error[i][j]=scaledEnRes[i][j]/math.sqrt(evtsPerEn)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Display graphs/histograms
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Display a few event profiles and their fit
evtFitCheckCanvas=TCanvas('evtFitCheckCanvas','100 GeV Shower Profiles and their Fit Gaussian',800,700)
evtFitCheckCanvas.Divide(2,2)
for i in range(4):
    evtFitCheckCanvas.cd(i+1)
    evtGraph[6][i].Draw()
    evtGraph[6][i].SetTitle("100 GeV - Event "+str(i)+";Depth [X_{0}];Total Measured Charge [# MIPs]")
    evtFitFnArray[6][i].Draw("same")
    evtGraph[6][i].GetYaxis().SetTitleOffset(1.5)

scaledEnRes=zip(*scaledEnRes) #transpose
error=zip(*error) #transpose
scaledEnResGraph=TGraphErrors(7,array('d',energyArray),array('d',scaledEnRes[0]),array('d',zeros),array('d',error[0]))
scaledEnResGraphCalib=TGraphErrors(7,array('d',energyArray_offset),array('d',scaledEnRes[1]),array('d',zeros),array('d',error[1]))
scaledEnResGraphCorr=TGraphErrors(7,array('d',energyArray_offset),array('d',scaledEnRes[2]),array('d',zeros),array('d',error[2]))

enResGraph1=TMultiGraph() #uncorrected and calibrated
enResGraph1.SetTitle("Scaled Energy Resolution; Initial Energy [GeV]; Energy~Resolution * #\sqrt{E}")
scaledEnResGraph.SetLineColor(kBlue)
scaledEnResGraph.SetMarkerColor(kBlue)
scaledEnResGraph.SetMarkerStyle(20)
scaledEnResGraph.SetMarkerSize(1.)
scaledEnResGraph.SetLineWidth(3)
scaledEnResGraphCalib.SetLineColor(kGreen+2)
scaledEnResGraphCalib.SetMarkerColor(kGreen+2)
scaledEnResGraphCalib.SetMarkerStyle(22)
scaledEnResGraphCalib.SetMarkerSize(1.)
scaledEnResGraphCalib.SetLineWidth(3)
enResGraph1.Add(scaledEnResGraph)
enResGraph1.Add(scaledEnResGraphCalib)

enResGraph2=TMultiGraph() #uncorrected and event-by-event corrected
enResGraph2.SetTitle("Scaled Energy Resolution; Initial Energy [GeV]; Energy~Resolution * #\sqrt{E}")
scaledEnResGraphCorr.SetLineColor(kRed)
scaledEnResGraphCorr.SetMarkerColor(kRed)
scaledEnResGraphCorr.SetMarkerStyle(21)
scaledEnResGraphCorr.SetMarkerSize(0.75)
scaledEnResGraphCorr.SetLineWidth(2)
enResGraph2.Add(scaledEnResGraph)
enResGraph2.Add(scaledEnResGraphCorr)

canvasUncorr=TCanvas('canvasUncorr','Uncorreced Scaled Energy Resolution',800,700)
scaledEnResGraph.Draw("AP")
scaledEnResGraph.GetYaxis().SetRangeUser(0.0,0.3)
canvasUncorr.SetLogx()
scaledEnResGraph.SetTitle("Scaled Energy Resolution; Initial Energy [GeV]; Energy~Resolution * #\sqrt{E}")
scaledEnResGraph.GetYaxis().SetTitleOffset(1.5)

canvasCalib=TCanvas('canvasCalib','Calibrated Scaled Energy Resolution',800,700)
canvasCalib.SetLogx()
enResGraph1.Draw("AP")
enResLegend1.AddEntry(scaledEnResGraph,"Uncalibrated","p")
enResLegend1.AddEntry(scaledEnResGraphCalib,"Calibrated by Shower Shape","p")
enResLegend1.Draw()
enResGraph1.GetYaxis().SetTitleOffset(1.5)

canvasCorr=TCanvas('canvasCorr','Corrected Scaled Energy Resolution',800,700)
canvasCorr.SetLogx()
enResGraph2.Draw("AP")
enResLegend2.AddEntry(scaledEnResGraph,"Uncalibrated","p")
enResLegend2.AddEntry(scaledEnResGraphCorr,"Corrected Event-by-event","p")
enResLegend2.Draw()
enResGraph2.GetYaxis().SetTitleOffset(1.5)

