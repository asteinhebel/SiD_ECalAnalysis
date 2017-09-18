from pyLCIO import IOIMPL, EVENT,UTIL
from pyLCIO.io.LcioReader import LcioReader
from array import array
from ROOT import *
import math

#my modules
import methods
import decodeLayerNos
import hitCollect
import fitting

##~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set initial variables
##~~~~~~~~~~~~~~~~~~~~~~~~~
#This code is robust for any phi that results in a shower encountering all 30 ECal layers fully (not in the overlap region) -> -15[deg]<phi<4 [deg] and all multiples
phi=0
evtsPerEn=5000
ecalLayers=31
hcalLayers=40

#This code can be run on data files with or without the 5T solenoid magnetic field, however the seeded fit parameters may need to be altered for fit convergence
inFile='reco_'+str(evtsPerEn)+'a.GeVscan.'+str(phi)+'phi.slcio'
nmbEvents=7*evtsPerEn
energyArray=[1,2,5,10,20,50,100]

###choose axis values###
xPos=False #if xPos, then GeV vs x position [mm]
           #if not xPos, then energy vs length (energy and length variables set below)

if not xPos:
    critEn=False #if critEn, then energy plotted as critical energy
                 #if not critEn, then energy plotted as average # of MIPs (based upon GeV sums)
    radiationLength=True #if radiationLength, then length plotted as radiation length
                          #if not radiationLength, then length plotted as "scaled radiation length" with ECal radiation lengths plotted as usual and HCal radiation lengths scaled by (dE/dx)_W/(dE/dx)_Fe

###choose cuts###
backscatterCut=False #if backscatterCut, then all hits outside a code of radius 0.2 radians from the original particle gun are removed in final plots
                     #if not backscatterCut, then all hits considered

###choose fitting method###
#linFit is the default. Gamma distributions can ONLY FIT MIP VS X0 PLOTS
linFit=False #if linFit, then deposits in the last 5 ECal layers are fit to an exponential and extrapolated into the HCal to estimate leakage
            #if not linFit, then the entire ECal shower is fit to a gamma distribution and extrapolated into the HCal to estimate leakage

if not linFit and not critEn and not radiationLength:
    print "ERROR - This fit can only be used for MIP vs radiation length plots. Changing fit option to 'linear fit'"
    linFit=True

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create histograms and graphs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
showerProfileGraph=[0]*7
showerProfileGraphLong=[0]*7

fullEnergyDistHist=[0]*7
fitLegend=[0]*7
for i in range(7):
    fullEnergyDistHist[i] = TH1D( 'TotalDepositedEnergy'+str(i), 'Event Energy Deposit ('+str(energyArray[i])+' GeV photons, phi='+str(phi)+', theta=90);Measured Deposited Energy [GeV];Entries', 5*(i+1)*(i+1), 0.,0.075*(i+1)*(i+1))
    fitLegend[i]=TLegend(0.15,0.2,0.35,0.4)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Record layer/module arrays
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
layerNosH, layerNos = decodeLayerNos.layerNos(nmbEvents,inFile)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cycle through data, store raw hit information, and fill and fit ECal energy distribution histogram
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
eventLayerArray=[0]*7
eventLayerArrayH=[0]*7
for i in range(7):
    eventLayerArray[i]=[0.]*evtsPerEn #stores sum of measured depopsits in each ECal layer for each event
    eventLayerArrayH[i]=[0.]*evtsPerEn #stores sum of measured depopsits in each HCal layer for each event

if not backscatterCut:
    hitCollect.hits(evtsPerEn,ecalLayers,hcalLayers,inFile,layerNos,layerNosH,fullEnergyDistHist,eventLayerArray,eventLayerArrayH)
else:
    hitCollect.hits2(phi,evtsPerEn,ecalLayers,hcalLayers,inFile,layerNos,layerNosH,fullEnergyDistHist,eventLayerArray,eventLayerArrayH)

#fit ECal distribution to a Gaussian
fitFunctionGauss=[0]*7
totEnArray=[0.]*7 #stores mean of each ECal distribution
for i in range(7):
    fitFunctionGauss[i]=TF1('fitFnGauss'+str(i),'gaus',0.,5.)
    fitFunctionGauss[i].SetLineColor( kRed )
    fitFunctionGauss[i].SetLineStyle( kDashed )
    fullEnergyDistHist[i].Fit( fitFunctionGauss[i], 'R' )
    totEnArray[i]=fitFunctionGauss[i].GetParameter(1)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create arrays to graph (x position/radiation length depth, measured energy deposits in each layer)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
layerMeanArray=[0]*7 #stores mean measured energy for each ECal layer
layerMeanArrayLong=[0]*7 #stores mean measured energy for each ECal and HCal layer

if xPos:
    layerPos=[0]*ecalLayers
    layerPosLong=[0]*(ecalLayers+hcalLayers)
else:
    radLen=[0]*ecalLayers
    radLenLong=[0]*(ecalLayers+hcalLayers)

#fill energy array
for i in range(7):
    layerMeanArray[i]=[0]*ecalLayers
    layerMeanArrayLong[i]=[0]*(ecalLayers+hcalLayers)
    for z in range(ecalLayers+hcalLayers):
        if z<21:
            for y in range(evtsPerEn):
                layerMeanArray[i][z]+=eventLayerArray[i][y][z]/float(evtsPerEn)
            if not xPos:
                if critEn:
                    layerMeanArray[i][z]=layerMeanArray[i][z]/0.0078
                else:
                    layerMeanArray[i][z]=layerMeanArray[i][z]/0.000124 #In Silicon, 1 MIP = 0.124 MeV
            layerMeanArrayLong[i][z]=layerMeanArray[i][z]
        elif z<ecalLayers:
            for y in range(evtsPerEn):
                layerMeanArray[i][z]+=eventLayerArray[i][y][z]/2./float(evtsPerEn) #eliminate sampling fraction correction from earlier
            if not xPos:
                if critEn:
                    layerMeanArray[i][z]=layerMeanArray[i][z]/0.0078
                else:
                    layerMeanArray[i][z]=layerMeanArray[i][z]/0.000124 #In Silicon, 1 MIP = 0.124 MeV
            layerMeanArrayLong[i][z]=layerMeanArray[i][z]
        else: #add on HCal info
            for y in range(evtsPerEn):
                layerMeanArrayLong[i][z]+=eventLayerArrayH[i][y][z-ecalLayers]/float(evtsPerEn)
            if not xPos:
                if critEn:
                    layerMeanArrayLong[i][z]=layerMeanArrayLong[i][z]/0.0213
                else:
                    layerMeanArrayLong[i][z]=layerMeanArrayLong[i][z]/0.000604 #In polystyrene (HCal scintillator), 1 MIP = 0.604 MeV

#fill position array
for z in range(ecalLayers+hcalLayers):
    if not xPos:
        if z==0:
            radLen[z]=0.
            radLenLong[z]=0.
        elif z<21:
            radLen[z]=radLen[z-1]+(2.5/math.cos(math.pi/180.*phi)*0.26)
            radLenLong[z]=radLenLong[z-1]+(2.5/math.cos(math.pi/180.*phi)*0.26)
        elif z<ecalLayers:
            radLen[z]=radLen[z-1]+(5./math.cos(math.pi/180.*phi)*0.26)
            radLenLong[z]=radLenLong[z-1]+(5./math.cos(math.pi/180.*phi)*0.26) # x0=3.85 mm of alloy
        else:
            if radiationLength:
                radLenLong[z]=radLenLong[z-1]+(20./math.cos(math.pi/180.*phi)*0.057) # X0=17.6 mm for stainless steel. Using Steel235 (99.8% Fe, 0.2% C) where for pure iron X0=17.57 mm
            elif not radiationLength:
                radLenLong[z]=radLenLong[z-1]+(20./math.cos(math.pi/180.*phi)*0.057)*20.2665/12.1338 #typical radiation length scaled by (dE/dx)_W/(dE/dx)_Fe
    else:
        if z==0:
            layerPos[z]=1264./math.cos(math.pi/180.*phi)
            layerPosLong[z]=layerPos[z]
        elif z<21:
            layerPos[z]=layerPos[z-1]+3.75/math.cos(math.pi/180.*phi)
            layerPosLong[z]=layerPos[z]
        elif z<ecalLayers:
            layerPos[z]=layerPos[z-1]+6.25/math.cos(math.pi/180.*phi)
            layerPosLong[z]=layerPos[z]
        else:
            layerPosLong[z]=layerPosLong[z-1]+27./math.cos(math.pi/180.*phi)/6. #approximate for radiation length difference

##~~~~~~~~~~~~~~~~~~~~~~~
#Display Final Graphs
##~~~~~~~~~~~~~~~~~~~~~~~
if xPos:
    for i in range(7):
        showerProfileGraph[i]=TGraph(31,array('d',layerPos),array('d',layerMeanArray[i]))

    #fit the last 5 layers to an exponential
    linFitFn,integral,params=fitting.linearFit1(layerPos[26],layerPos[30],showerProfileGraph,totEnArray)

    linFitCanvas=[0]*7
    for i in range(7):
        linFitCanvas[i]=TCanvas('linFitCanvas'+str(i),str(energyArray[i])+' GeV Shower Profile in [mm] vs [GeV] with Linear Fit',800,700)
        showerProfileGraphLong[i]=TGraph(71,array('d',layerPosLong),array('d',layerMeanArrayLong[i]))
        showerProfileGraphLong[i].Draw("AC*")
        showerProfileGraphLong[i].SetTitle("Mean Energy per Layer vs x position in ECal and HCal ("+str(evtsPerEn)+" "+str(energyArray[i])+" GeV photons, theta=90, phi="+str(phi)+");x position [mm];Mean Measured Deposits per Layer [GeV]")
        linFitFn[i].Draw("same")
        fitLegend[i].AddEntry(showerProfileGraphLong[i],"Fit function = b+mx","")
        fitLegend[i].AddEntry(showerProfileGraphLong[i],"b = "+str(params[i][0]),"")
        fitLegend[i].AddEntry(showerProfileGraphLong[i],"m = "+str(params[i][1]),"")
        fitLegend[i].Draw()
        showerProfileGraphLong[i].GetYaxis().SetTitleOffset(1.5)
        linFitCanvas[i].SetLogy()
        showerProfileGraphLong[i].GetYaxis().SetRangeUser(0.00000001,1.)
 
    print "Integral values: "+str(integral)


if not xPos:
    for i in range(7):
        #convert GeV to appropriate units
        if critEn:
            totEnArray[i]=totEnArray[i]/0.0078
        else:
            totEnArray[i]=totEnArray[i]/0.000124
        showerProfileGraph[i]=TGraph(31,array('d',radLen),array('d',layerMeanArray[i]))

    if linFit:
        #fit the last 5 layers to an exponential
        fitFn,integral,params=fitting.linearFit2(radLen[26],radLen[30],showerProfileGraph,totEnArray) 
    else: 
        #gamma fit - plotting MIPs vs radiation lengths
        fitFn,integral,params=fitting.gammaFit(radLen[30],showerProfileGraph,totEnArray)

    showerCanvas=[0]*7
    calLine=TLine(26.5,0,26.5,1000)
    for i in range(7):           
        showerProfileGraphLong[i]=TGraph(71,array('d',radLenLong),array('d',layerMeanArrayLong[i]))
        showerCanvas[i]=TCanvas('showerCanvas'+str(i),str(energyArray[i])+' GeV Shower Profile with Leakage Extrapolation Fit',800,700)
        showerProfileGraphLong[i].Draw("AC*")
        if not radiationLength and not critEn:
            showerProfileGraphLong[i].SetTitle("Mean Energy per Layer vs Scaled Radiation Length X_{0} in ECal and HCal ("+str(evtsPerEn)+" "+str(energyArray[i])+" GeV photons, theta=90, phi="+str(phi)+");(Scaled HCal) Radiation Lengths; Average Number of MIPs")
            showerProfileGraphLong[i].GetYaxis().SetRangeUser(0.9,1000.)
        elif not radiationLength and critEn:
            showerProfileGraphLong[i].SetTitle("Mean Energy per Layer vs Scaled Radiation Length X_{0} in ECal and HCal ("+str(evtsPerEn)+" "+str(energyArray[i])+" GeV photons, theta=90, phi="+str(phi)+");(Scaled HCal) Radiation Lengths;Critical Energies [E_{C}]")
            showerProfileGraphLong[i].GetYaxis().SetRangeUser(0.0001,20.)
        elif radiationLength and not critEn:
            showerProfileGraphLong[i].SetTitle("Mean Energy per Layer vs Radiation Length X_{0} in ECal and HCal ("+str(evtsPerEn)+" "+str(energyArray[i])+" GeV photons, theta=90, phi="+str(phi)+");Radiation Lengths [X_{0}]; Average Number of MIPs")
            showerProfileGraphLong[i].GetYaxis().SetRangeUser(0.0001,1000.)
        elif radiationLength and critEn:
            showerProfileGraphLong[i].SetTitle("Mean Energy per Layer vs Radiation Length X_{0} in ECal and HCal ("+str(evtsPerEn)+" "+str(energyArray[i])+" GeV photons, theta=90, phi="+str(phi)+");Radiation Lengths [X_{0}]; Critical Energies [E_{C}]")
            showerProfileGraphLong[i].GetYaxis().SetRangeUser(0.0001,20.)
        if linFit:
            fitLegend[i].AddEntry(showerProfileGraphLong[i],"Fit function = b+mx","")
            fitLegend[i].AddEntry(showerProfileGraphLong[i],"b = "+str(params[i][0]),"")
            fitLegend[i].AddEntry(showerProfileGraphLong[i],"m = "+str(params[i][1]),"")
            fitLegend[i].Draw()
        else:
            fitLegend[i].AddEntry(showerProfileGraphLong[i],"Fit function = A x^{#alpha} e^{-bx}","")
            fitLegend[i].AddEntry(showerProfileGraphLong[i],"A = "+str(params[i][0]),"")
            fitLegend[i].AddEntry(showerProfileGraphLong[i],"#alpha = "+str(params[i][1]),"")
            fitLegend[i].AddEntry(showerProfileGraphLong[i],"b = "+str(params[i][2]),"")
            fitLegend[i].Draw()
        showerProfileGraphLong[i].GetYaxis().SetTitleOffset(1.3)
        showerCanvas[i].SetLogy()
        fitFn[i].Draw("same")
        calLine.Draw("same")

    print "Integral values: "+str(integral)
