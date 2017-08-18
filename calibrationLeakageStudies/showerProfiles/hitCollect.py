from pyLCIO import IOIMPL, EVENT,UTIL
from pyLCIO.io.LcioReader import LcioReader
from array import array
from ROOT import *
import math
#my modules
import methods

#######################################################################################################
def hits(evtsPerEn,ecalLayers,hcalLayers,inFile,layerNos,layerNosH,hitEnergyHistogram,eventLayerArray,eventLayerArrayH):
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open( inFile )
    count=-1
    evNo=0
    # loop over all events in the file
    for dep in reader:
        count+=1
        if count%100==0:
            print "Summing energy of event "+str(count)
        hitTotal=0.	
        layerTotal=[0.]*ecalLayers
        layerTotalH=[0.]*hcalLayers
        # get the collection from the event
        hitCollection = dep.getCollection( 'ECalBarrelHits' )   
        hitCollectionH = dep.getCollection( 'HCalBarrelHits' ) 

        hitNo=0
        for hit in hitCollection: #exclude layer0
            if layerNos[count][hitNo]>0 and layerNos[count][hitNo]<21:
                hitTotal+=(hit.getEnergy())	
                layerTotal[layerNos[count][hitNo]]+=hit.getEnergy()	
            elif layerNos[count][hitNo]>0:
                hitTotal+=(hit.getEnergy()*2)
                layerTotal[layerNos[count][hitNo]]+=2*hit.getEnergy()
            hitNo+=1

        hitNo4=0
        for hit2 in hitCollectionH:
            layerTotalH[layerNosH[count][hitNo4]]+=hit2.getEnergy()
            hitNo4+=1    
        
        #record sums of deposits at each energy for later use
        if count<evtsPerEn*(evNo+1):
            hitEnergyHistogram[evNo].Fill(hitTotal)
            eventLayerArray[evNo][count-evNo*evtsPerEn]=[x for x in layerTotal]
            eventLayerArrayH[evNo][count-evNo*evtsPerEn]=[y for y in layerTotalH]
        if count==(evtsPerEn-1)*(evNo+1)+evNo:
            evNo+=1

    reader.close()

    return hitEnergyHistogram,eventLayerArray,eventLayerArrayH 

########################################################################################################
#cuts on backscattered hits  >= 0.2 rad from beam 
def hits2(phi,evtsPerEn,ecalLayers,hcalLayers,inFile,layerNos,layerNosH,hitEnergyHistogram,eventLayerArray,eventLayerArrayH):
    phiRad=phi*math.pi/180.
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open( inFile )
    count=-1
    evNo=0
    # loop over all events in the file
    for dep in reader:
        count+=1
        if count%100==0:
            print "Summing energy of event "+str(count)
        hitTotal=0.	
        layerTotal=[0.]*ecalLayers
        layerTotalH=[0.]*hcalLayers
        # get the collection from the event
        hitCollection = dep.getCollection( 'ECalBarrelHits' )   
        hitCollectionH = dep.getCollection( 'HCalBarrelHits' ) 

        hitNo=0
        for hit in hitCollection: #exclude layer0
            if methods.solidAngle(phi,hit)<phiRad+0.2 and methods.solidAngle(phi,hit)>phiRad-0.2: #hits within cone of 0.2 radians of original particle gun trajectory
                if layerNos[count][hitNo]>0 and layerNos[count][hitNo]<21:
                    hitTotal+=(hit.getEnergy())	
                    layerTotal[layerNos[count][hitNo]]+=hit.getEnergy()
                elif layerNos[count][hitNo]>0:
                    hitTotal+=(hit.getEnergy()*2)
                    layerTotal[layerNos[count][hitNo]]+=2*hit.getEnergy()
                hitNo+=1

        hitNo2=0
        for hit2 in hitCollectionH:
            if methods.solidAngle(phi,hit)<phiRad+0.2 and methods.solidAngle(phi,hit)>phiRad-0.2: #hits within cone of 0.2 radians of original particle gun trajectory
                layerTotalH[layerNosH[count][hitNo2]]+=hit2.getEnergy()
                hitNo2+=1    

        #record sums of deposits at each energy for later use
        if count<evtsPerEn*(evNo+1):
            hitEnergyHistogram[evNo].Fill(hitTotal)
            eventLayerArray[evNo][count-evNo*evtsPerEn]=[x for x in layerTotal]
            eventLayerArrayH[evNo][count-evNo*evtsPerEn]=[y for y in layerTotalH]
        if count==(evtsPerEn-1)*(evNo+1)+evNo:
            evNo+=1

    reader.close()

    return hitEnergyHistogram,eventLayerArray,eventLayerArrayH
#######################################################################################################
