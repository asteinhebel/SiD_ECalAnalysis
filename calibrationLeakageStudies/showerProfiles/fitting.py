from pyLCIO import IOIMPL, EVENT,UTIL
from pyLCIO.io.LcioReader import LcioReader
from array import array
from ROOT import *
import math

#################################################################################
# for x in [mm]
def linearFit1(ecalMid,ecalEnd,graph,totEnArray):
    fitFn=[0]*7
    integral=[0.]*7
    for en in range(7):
        fitFn[en]=TF1('linFitFn'+str(en),'expo(0)',ecalMid,1700) #ecalMid defines the 5th to last ECal layer where fitting begins
        graph[en].Fit(fitFn[en], "R" )
        linB=fitFn[en].GetParameter(0)
        slope=fitFn[en].GetParameter(1)
        xInt=(1/slope)*(math.log(0.000000001)-linB) #until y=10e-10
        integral[en]=fitFn[en].Integral(ecalEnd,xInt,1.e-12)/totEnArray[en]/(6.25/2.) #integrate up to y=10e-10 and divide by 6.25/2 mm from extrapolating from thick layers but doubling their deposits in the energy sum

    return fitFn,integral

##################################################################################
# for x in [X_0]
def linearFit2(ecalMid,ecalEnd,graph,totEnArray):
    fitFn=[0]*7
    integral=[0.]*7
    for en in range(7):
        fitFn[en]=TF1('linFitFn'+str(en),'expo(0)',ecalMid,1700) #ecalMid defines the 5th to last ECal layer where fitting begins
        graph[en].Fit(fitFn[en], "R" )
        linB=fitFn[en].GetParameter(0)
        slope=fitFn[en].GetParameter(1)
        xInt=(1/slope)*(math.log(0.000000001)-linB) #until y=10e-10
        integral[en]=fitFn[en].Integral(ecalEnd,xInt,1.e-12)/totEnArray[en]/(6.25/2.) #integrate up to y=10e-10 and divide by 6.25/2 mm from extrapolating from thick layers but doubling their deposits in the energy sum

    return fitFn,integral

##################################################################################
# for x in [X_0] and y in MIP *only* - gamma distribution fit of entire ECal shower
def gammaFit(ecalEnd,graph):
    fitFn=[0]*7
    integral=[0.]*7
    xInt=(-1/0.3)*math.log(1e-11) #where y=10e-11 assuming exp(-0.3x) (x~84)
    for i in range(7):
        fitFn[i]=TF1("gammaFitFn"+str(i),"[0]*TMath::Power(x,[1])*TMath::Exp(-[2]*x)",0,70)
        #a few first guesses to convince fits to converge (optimized for runs WITH 5T magnetic field)
        if i==6:
            fitFn[i].SetParameter(1,3.)
        elif i==3 or i==4:
            fitFn[i].SetParameter(0,2.)
            fitFn[i].SetParameter(1,3.)
        else:
            fitFn[i].SetParameter(1,4.)
        graph[i].Fit(fitFn[i])
        integral[i]=fitFn[i].Integral(ecalEnd,xInt,1.e-12)*50./57. #integrate up to y=10e-11 and multiply by HCal X_0/layer to get idea of HCal energy that's unaccounted for

    return fitFn, integral
#################################################################################
