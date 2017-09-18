from pyLCIO import IOIMPL, EVENT,UTIL
from pyLCIO.io.LcioReader import LcioReader
from array import array
from ROOT import *
import math

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set global variables and import data files
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
phi=[0,3.75,7.5,9.3,11.25,15,18.75,22.5,30] #Each input data file contains showers at one phi location. This array defines which of those files are pulled into the code for analysis. Here, 14 files are pulled in - 7 phi locations for both 100 GeV and 10 GeV EM showers
phi1=[x+0.5 for x in phi] #to offset the data points for easy visual comparison

nmbEvents=500 #number of showers in each file
nmbLayers=32 #number of ECal layers
inFile=[0]*len(phi)
readerL=[0]*len(phi)
inFile1=[0]*len(phi)
readerL1=[0]*len(phi)
for i in range(len(phi)):
    inFile[i]='reco_500a.10GeV.'+str(phi[i])+'phi.slcio'
    readerL[i] = LcioReader(inFile[i])
    inFile1[i]='reco_500a.100GeV.'+str(phi[i])+'phi.slcio'
    readerL1[i] = LcioReader(inFile1[i])
#All objects associated with 100 GeV runs has an additional "1" on their name, to differentiate from the 10 GeV-associated objects

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create histograms and legends to compare ECal and HCal hits as a function of phi at 100 GeV and 10 GeV
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sumHistE=[0]*len(phi)
sumHistH=[0]*len(phi)
sumHistE1=[0]*len(phi)
sumHistH1=[0]*len(phi)
for i in range(len(phi)):
    sumHistE[i] = TH1D( 'totalEcal'+str(i), 'Event Energy Deposit (10 GeV photons, phi='+str(phi[i])+', theta=90);ECAL Barrel Hit Energy [GeV];Entries', 200, 0., 2. )
    sumHistH[i] = TH1D( 'totalHcal'+str(i), 'Event Energy Deposit (10 GeV photons, phi='+str(phi[i])+', theta=90);HCAL Barrel Hit Energy [GeV];Entries', 20, 0., 0.1 )
    sumHistE1[i] = TH1D( 'totalEcal'+str(i)+str(i), 'Event Energy Deposit (100 GeV photons, phi='+str(phi[i])+', theta=90);ECAL Barrel Hit Energy [GeV];Entries', 200, 0., 2. )
    sumHistH1[i] = TH1D( 'totalHcal'+str(i)+str(i), 'Event Energy Deposit (100 GeV photons, phi='+str(phi[i])+', theta=90);HCAL Barrel Hit Energy [GeV];Entries', 20, 0., 0.1 )

ecalLegend=TLegend(0.8,0.8,0.95,0.95)
hcalLegend=TLegend(0.5,0.5,0.7,0.7)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Decode ECal layer information for proper sampling fraction weighting of energy deposits
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
layerNos=[0]*len(phi)
layerNos1=[0]*len(phi)
for i in range(len(phi)):
    layerNos[i]=[0]*nmbEvents
    layerNos1[i]=[0]*nmbEvents
    for j in range(nmbEvents):
        layerNos[i][j]=[]
        layerNos1[i][j]=[]
print "Layer array created"

for i in range(len(phi)):
    ###for 10 GeV files###
    eventNo=-1
    # loop over the events
    for event in readerL[i]:
        eventNo+=1
        if eventNo%100==0:
            print "10 GeV, phi "+str(phi[i])+", event "+str(eventNo)        
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
            layerNos[i][eventNo].append(idDecoder['layer'].value())

    ###for 100 GeV files###
    eventNo=-1
    for event in readerL1[i]:
        eventNo+=1
        if eventNo%100==0:
            print "100 GeV, phi "+str(phi[i])+", event "+str(eventNo)
        ecalHits = event.getCollection( 'ECalBarrelHits' )
        cellIdEncoding = ecalHits.getParameters().getStringVal( EVENT.LCIO.CellIDEncoding )
        idDecoder = UTIL.BitField64( cellIdEncoding )
        for caloHit in ecalHits:
            cellID = long( caloHit.getCellID0() & 0xffffffff ) | ( long( caloHit.getCellID1() ) << 32 )
            idDecoder.setValue( cellID )
            layerNos1[i][eventNo].append(idDecoder['layer'].value())


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Collect and sum up energy deposits in ECal and HCal
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()

for i in range(len(phi)):
    ###for 10 GeV files###
    reader.open( inFile[i] )
    j=-1
    # loop over all events in the file
    for dep in reader:
        j+=1
        if j%100==0:
            print "Summing energy of event "+str(j)+", 10 GeV, phi = "+str(phi[i])
        hitTotalE=0.	
        hitTotalH=0.
        # get uncalibrated collections from the event
        hitCollectionH = dep.getCollection( 'HCalBarrelHits' )
        hitCollectionE = dep.getCollection( 'ECalBarrelHits' )

        # loop over all hits in the collections and weight those in ECal to take geometry and sampling fraction into account
        hitNo=0
        for hit in hitCollectionE:
            if layerNos[i][j][hitNo]>0 and layerNos[i][j][hitNo]<21:
                hitTotalE+=(hit.getEnergy())	
            elif layerNos[i][j][hitNo]>0 and layerNos[i][j][hitNo]<nmbLayers:
                hitTotalE+=(hit.getEnergy()*2)
            hitNo+=1
        for hit in hitCollectionH:
            hitTotalH+=hit.getEnergy()

        #fill histograms
        sumHistH[i].Fill(hitTotalH)   
        sumHistE[i].Fill(hitTotalE)
    reader.close()

    ###for 100 GeV files###
    reader.open( inFile1[i] )
    j=-1
    # loop over all events in the file
    for dep in reader:
        j+=1
        if j%100==0:
            print "Summing energy of event "+str(j)+", 100 GeV, phi = "+str(phi[i])
        hitTotalE=0.	
        hitTotalH=0.
        # get uncalibrated collections from the event
        hitCollectionH = dep.getCollection( 'HCalBarrelHits' )
        hitCollectionE = dep.getCollection( 'ECalBarrelHits' )

        # loop over all hits in the collections and weight those in ECal to take geometry and sampling fraction into account
        hitNo=0
        for hit in hitCollectionE:
            if layerNos1[i][j][hitNo]>0 and layerNos1[i][j][hitNo]<21:
                hitTotalE+=(hit.getEnergy())	
            elif layerNos1[i][j][hitNo]>0 and layerNos1[i][j][hitNo]<nmbLayers:
                hitTotalE+=(hit.getEnergy()*2)
            hitNo+=1
        for hit in hitCollectionH:
            hitTotalH+=hit.getEnergy()

        #fill histograms
        sumHistH1[i].Fill(hitTotalH)   
        sumHistE1[i].Fill(hitTotalE)
    reader.close()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fit ECal distributions with Gaussian curve. Pull out mean and standard deviation of these Gaussian functions, as well as histogram mean and RMS for HCal histograms.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fitE=[0]*len(phi)
fitE1=[0]*len(phi)
meanE=[0]*len(phi)
meanE1=[0]*len(phi)
sigmaE=[0]*len(phi)
sigmaE1=[0]*len(phi)
meanH=[0]*len(phi)
meanH1=[0]*len(phi)
sigmaH=[1/math.sqrt(nmbEvents)]*len(phi)
zeros=[0]*len(phi)
for i in range(len(phi)):
    fitE[i]=TF1('fitFnGauss'+str(i),'gaus',0.,2.)
    fitE1[i]=TF1('fitFnGauss1'+str(i),'gaus',0.,2.)
    sumHistE[i].Fit( fitE[i], 'R' )
    sumHistE1[i].Fit( fitE1[i], 'R' )
    meanE[i]=fitE[i].GetParameter(1)*10.
    meanE1[i]=fitE1[i].GetParameter(1)
    sigmaE[i]=fitE[i].GetParameter(2)*math.sqrt(10.)
    sigmaE1[i]=fitE1[i].GetParameter(2)
    meanH[i]=sumHistH[i].GetMean()*10./meanE[i]*100.
    meanH1[i]=sumHistH1[i].GetMean()/meanE1[i]*100.

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define graphs for ECal and HCal comparisons, and draw graphs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
graphE=TGraphErrors(len(phi),array('d',phi),array('d',meanE),array('d',zeros),array('d',sigmaE))
graphE.SetMarkerStyle(21)
graphE.SetMarkerColor(kRed)
graphE.SetMarkerSize(2)
graphE.SetLineColor(kRed)
graphE.SetLineWidth(3)
graphE1=TGraphErrors(len(phi),array('d',phi1),array('d',meanE1),array('d',zeros),array('d',sigmaE1))
graphE1.SetMarkerStyle(8)
graphE1.SetMarkerColor(kBlue)
graphE1.SetMarkerSize(2)
graphE1.SetLineColor(kBlue)
graphE1.SetLineWidth(3)
graphH=TGraphErrors(len(phi),array('d',phi),array('d',meanH),array('d',zeros),array('d',sigmaH))
graphH.SetMarkerStyle(21)
graphH.SetMarkerColor(kRed)
graphH.SetMarkerSize(1.5)
graphH.SetLineColor(kRed)
graphH1=TGraphErrors(len(phi),array('d',phi),array('d',meanH1),array('d',zeros),array('d',sigmaH))
graphH1.SetMarkerStyle(8)
graphH1.SetMarkerColor(kBlue)
graphH1.SetMarkerSize(1.5)
graphH1.SetLineColor(kBlue)

totGraphE=TMultiGraph()
totGraphE.SetTitle("ECal Hits;#varphi [deg]; Total Measured Charge [GeV]")
totGraphE.Add(graphE)
totGraphE.Add(graphE1)
ecalLegend.AddEntry(graphE,"10 GeV (x10)","lep")
ecalLegend.AddEntry(graphE1,"100 GeV","lep")
totGraphH=TMultiGraph()
totGraphH.SetTitle("HCal Hits;#varphi [deg]; % Shower Energy in HCal ")
totGraphH.Add(graphH)
totGraphH.Add(graphH1)
hcalLegend.AddEntry(graphH,"10 GeV","lep")
hcalLegend.AddEntry(graphH1,"100 GeV","lep")

#Draw graphs
ecalCanvas = TCanvas( 'ECal Canvas', 'ECalCanvas', 800, 700 )
totGraphE.Draw("AP")
totGraphE.GetYaxis().SetTitleOffset(1.3)
ecalLegend.Draw()

hcalCanvas=TCanvas('HCal Canvas','HCalCanvas',800,700)
totGraphH.Draw("AP")
totGraphH.GetYaxis().SetTitleOffset(1.3)
hcalLegend.Draw()
