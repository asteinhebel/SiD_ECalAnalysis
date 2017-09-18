from pyLCIO import IOIMPL, EVENT,UTIL
from pyLCIO.io.LcioReader import LcioReader
from array import array
from ROOT import *
import math

##~~~~~~~~~~~~~~~~~~~~
# Define methods
##~~~~~~~~~~~~~~~~~~~~
############################################################################
def solidAngle(phi,hit):
    angle=(math.acos((math.cos(phi*math.pi/180.)*hit.getPosition()[0]+math.sin(phi*math.pi/180.)*hit.getPosition()[1])/math.sqrt(hit.getPosition()[0]*hit.getPosition()[0]+hit.getPosition()[1]*hit.getPosition()[1]+hit.getPosition()[2]*hit.getPosition()[2])))
    if hit.getPosition()[1]<0:
        angle=-angle 
    return angle
############################################################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Input global variables and data file
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
phi=0 #angle of incidence - from particle gun
phir=phi*math.pi/180. #phi in radians
pRange=0.1 #Specify backscatter cut. Hits with solid angle outside phir +/- pRange are considered backscatter
evtsPerEn=500
inFile='reco_'+str(evtsPerEn)+'a.GeVScan.'+str(phi)+'phi.slcio'
nmbEvents=7*evtsPerEn
readerL = LcioReader(inFile ) #create a reader

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create histograms and legends
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
energyArray=[1,2,5,10,20,50,100] #simulated energies in data file
angleHistogram=[0]*7
xHistogram=[0]*7
yHistogram=[0]*7
zHistogram=[0]*7
backscatterEHist=[0]*7
centerscatterEHist=[0]*7
scatterLegend=[0]*7
for i in range(7):
    angleHistogram[i]=TH1D('showerSpread'+str(i),'Angle Between Beam and Shower Hits (500 photon showers, phi='+str(phi)+', theta=90);Angle [rad];Entries', 32,-math.pi,math.pi)
    xHistogram[i]=TH1D('showerSpreadx'+str(i),'x position of hits ( phi='+str(phi)+', theta=90);Position [mm];Entries',240,1265.,1405.)    
    yHistogram[i]=TH1D('showerSpready'+str(i),'y position of hits ( phi='+str(phi)+', theta=90);Position [mm];Entries', 200,-100.,100.)    
    zHistogram[i]=TH1D('showerSpreadz'+str(i),'z position of hits ( phi='+str(phi)+', theta=90);Position [mm];Entries', 200,-100.,100.)
    backscatterEHist[i]=TH1D('Backscatter'+str(i),'Backscattered Hits ('+str(evtsPerEn)+' '+str(energyArray[i])+' GeV photons, phi='+str(phi)+', theta=90);Entries;Total Measured Deposits in ECal [MeV]',56, 0.,6.944)#bins of 1 MIP
    centerscatterEHist[i]=TH1D('centerscatter'+str(i),'Measured Energy of Hits ('+str(evtsPerEn)+' '+str(energyArray[i])+' GeV photons, phi='+str(phi)+', theta=90, bins = 1 MIP);Measured Hit Deposits in ECal Barrel [MeV];Entries',56, 0.,6.944)#bins of 1 MIP
    scatterLegend[i]=TLegend(0.4,0.5,0.9,0.7)
xLegend=TLegend(0.7,0.7,0.8,1)
yLegend=TLegend(0.7,0.7,0.8,1)
zLegend=TLegend(0.7,0.7,0.8,1)
angleLegend=TLegend(0.7,0.7,0.9,0.9)

##~~ get layer ~~###############################################
layerNos=[0]*nmbEvents
for i in range(nmbEvents):
    layerNos[i]=[]
print "Layer array created"

eventNo=-1
# loop over the events
for event in readerL:
    eventNo+=1
    if eventNo%100==0:
        print "recording layers of event "+str(eventNo)
#    print eventNo
    # get a hit collection
    ecalHits = event.getCollection( 'ECalBarrelHits' )
    # get the cell ID encoding string from the collection parameters
    cellIdEncoding = ecalHits.getParameters().getStringVal( EVENT.LCIO.CellIDEncoding )
#    print cellIdEncoding
    # define a cell ID decoder for the collection
    idDecoder = UTIL.BitField64( cellIdEncoding )
    # loop over all hits in the collection
    for caloHit in ecalHits:
        # combine the two 32 bit cell IDs of the hit into one 64 bit integer
        cellID = long( caloHit.getCellID0() & 0xffffffff ) | ( long( caloHit.getCellID1() ) << 32 )
        # set up the ID decoder for this cell ID
        idDecoder.setValue( cellID )
        # access the field information using a valid field from the cell ID encoding string
#        print 'layer:', idDecoder['layer'].value()
        layerNos[eventNo].append(idDecoder['layer'].value())

###############################################################################

# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open( inFile )
count=-1
evNo=0 #tracks what energy the event has for histogram filling loop
backEn=[0.0]*7
centerEn=[0.0]*7
totLayerSum=[0.0]*32
totHitNo=[0]*32
# loop over all events in the file
for dep in reader:
    count+=1
    if count%100==0:
	print "Summing energy of event "+str(count)
    # get the collection from the event
    hitCollection = dep.getCollection( 'ECalBarrelHits' ) 
    hitNo=0
    for hit in hitCollection:
        angle=solidAngle(phi,hit)
        if count<evtsPerEn*(evNo+1):
            angleHistogram[evNo].Fill(angle)
            if angle<phir+pRange and angle>phir-pRange:
                if layerNos[count][hitNo]<21:
                    centerEn[evNo]+=hit.getEnergy()*100
                else:
                    centerEn[evNo]+=hit.getEnergy()*200
                centerscatterEHist[evNo].Fill(hit.getEnergy()*100)
            else:
                if layerNos[count][hitNo]<21:
                    backEn[evNo]+=hit.getEnergy()*100
                else:
                    backEn[evNo]+=hit.getEnergy()*200
                backscatterEHist[evNo].Fill(hit.getEnergy()*100)
            xHistogram[evNo].Fill(hit.getPosition()[0])
            yHistogram[evNo].Fill(hit.getPosition()[1])
            zHistogram[evNo].Fill(hit.getPosition()[2])
        hitNo+=1
    if count==(evtsPerEn-1)*(evNo+1)+evNo:
        evNo+=1
reader.close()


#plot histograms
angleCanvas=[0]*7
angleCanvas=TCanvas('angleCanvas','Solid Angle Distribution of Hits',800,700)
angleCanvas.SetLogy()
for i in range(7):
    angleHistogram[6-i].Draw("same")
    angleHistogram[6-i].SetStats(0)
    angleHistogram[6-i].SetLineColor(7-i)
    angleHistogram[6-i].SetStats(0)
    angleLegend.AddEntry(angleHistogram[6-i],str(energyArray[6-i])+' GeV','l')
angleLegend.Draw()

spatialCanvas=TCanvas('spatialCanvas','Spatial Distribution of Hits',1000,1100)
spatialCanvas.Divide(1,3)
spatialCanvas.cd(1)
for i in range(7):
    xHistogram[6-i].Draw("same")
    xHistogram[6-i].SetFillColor(7-i)
    xHistogram[6-i].SetStats(0)
    xLegend.AddEntry(xHistogram[6-i],str(energyArray[6-i])+' GeV','f')
xLegend.Draw()
spatialCanvas.cd(2)
for i in range(7):
    yHistogram[6-i].Draw("same")
    yHistogram[6-i].SetFillColor(7-i)
    yHistogram[6-i].SetStats(0)
    yLegend.AddEntry(yHistogram[6-i],str(energyArray[6-i])+' GeV','f')
yLegend.Draw()
spatialCanvas.cd(3)
for i in range(7):
    zHistogram[6-i].Draw("same")
    zHistogram[6-i].SetFillColor(7-i)
    zHistogram[6-i].SetStats(0)
    zLegend.AddEntry(zHistogram[6-i],str(energyArray[6-i])+' GeV','f')
zLegend.Draw()

backscatterCanvas=[0]*7
for i in range(7):
    backscatterCanvas[i]=TCanvas('backscatterCanvas'+str(i),'Backscattered vs Centered Hits'+str(i),800,700)
    backscatterCanvas[i].SetLogy()
    centerscatterEHist[i].Draw()
    centerscatterEHist[i].SetStats(0)
    centerscatterEHist[i].SetFillColor(kBlue)
    scatterLegend[i].AddEntry(centerscatterEHist[i],"Centered Hits ("+str(phir-pRange)+" rad < #Omega < "+str(phir+pRange)+" rad)","f")
    backscatterEHist[i].Draw("same")
    backscatterEHist[i].SetStats(0)
    backscatterEHist[i].SetFillColor(kRed)
    scatterLegend[i].AddEntry(backscatterEHist[i],"Backscattered Hits ("+str(phir-pRange)+" rad >#Omega or #Omega > "+str(phir+pRange)+" rad)","f")
    scatterLegend[i].AddEntry("",str(backEn[i]/(centerEn[i]+backEn[i])*100.)+"% of energy in backscatter","")
    scatterLegend[i].Draw()
