from pyLCIO import IOIMPL, EVENT,UTIL
from pyLCIO.io.LcioReader import LcioReader
from array import array
from ROOT import *
import math

##############################################################################
def layerNos(nmbEvents,inFile):
    layerNosH=[0]*nmbEvents
    layerNos=[0]*nmbEvents
    for i in range(nmbEvents):
        layerNos[i]=[]
        layerNosH[i]=[]

    #create a reader
    readerL = LcioReader(inFile )
    eventNo=-1
    # loop over the events
    for event in readerL:
        eventNo+=1
        if eventNo%100==0:
            print "recording layers of event "+str(eventNo)
        # get a hit collection
        hcalHits = event.getCollection( 'HCalBarrelHits' )
        ecalHits = event.getCollection( 'ECalBarrelHits' )
        # get the cell ID encoding string from the collection parameters
        cellIdEncoding = ecalHits.getParameters().getStringVal( EVENT.LCIO.CellIDEncoding )
        cellIdEncodingH = hcalHits.getParameters().getStringVal( EVENT.LCIO.CellIDEncoding )
        # define a cell ID decoder for the collection
        idDecoder = UTIL.BitField64( cellIdEncoding )
        idDecoderH = UTIL.BitField64( cellIdEncodingH )
        # loop over all hits in the collection
        for caloHit in ecalHits:
            # combine the two 32 bit cell IDs of the hit into one 64 bit integer
            cellID = long( caloHit.getCellID0() & 0xffffffff ) | ( long( caloHit.getCellID1() ) << 32 )
            # set up the ID decoder for this cell ID
            idDecoder.setValue( cellID )
            # access the field information using a valid field from the cell ID encoding string
            layerNos[eventNo].append(idDecoder['layer'].value())
        for caloHitH in hcalHits:
            cellIDH=long(caloHitH.getCellID0() & 0xffffffff)|(long(caloHitH.getCellID1())<<32)
            idDecoderH.setValue(cellIDH)
            layerNosH[eventNo].append(idDecoderH['layer'].value())

    return layerNosH, layerNos
##############################################################################
