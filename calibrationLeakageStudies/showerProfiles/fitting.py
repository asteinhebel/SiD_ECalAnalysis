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
    params=[0]*7
    roundParams=[0]*7
    for i in range(7):
        fitFn[i]=TF1('linFitFn'+str(i),'expo(0)',ecalMid,1700) #ecalMid defines the 5th to last ECal layer where fitting begins
        graph[i].Fit(fitFn[i], "R" )
        params[i]=fitFn[i].GetParameters()
        xInt=(1/params[i][1])*(math.log(0.000000001)-params[i][0]) #until y=10e-10
        integral[i]=fitFn[i].Integral(ecalEnd,xInt,1.e-12)/totEnArray[i]/(6.25/2.) #integrate up to y=10e-10 and divide by 6.25/2 mm from extrapolating from thick layers but doubling their deposits in the energy sum
        roundParams[i]=[0.]*2
        for j in range(2):
            roundParams[i][j]="{0:.4f}".format(round(params[i][j],4))

    return fitFn,integral,roundParams

##################################################################################
# for x not in [mm]
def linearFit2(ecalMid,ecalEnd,graph,totEnArray):
    fitFn=[0]*7
    integral=[0.]*7
    params=[0]*7
    roundParams=[0]*7
    for i in range(7):
        fitFn[i]=TF1('linFitFn'+str(i),'expo(0)',ecalMid,175) #ecalMid defines the 5th to last ECal layer where fitting begins
        graph[i].Fit(fitFn[i], "R" )
        params[i]=fitFn[i].GetParameters()
        xInt=(1/params[i][1])*(math.log(0.000000001)-params[i][0]) #until y=10e-10
        integral[i]=fitFn[i].Integral(ecalEnd,xInt,1.e-12)/totEnArray[i]/(5.*0.26) #integrate up to y=10e-10 and divide by 5 mm*0.26 from extrapolating from radiation lengths of thick layers
        roundParams[i]=[0.]*2
        for j in range(2):
            roundParams[i][j]="{0:.4f}".format(round(params[i][j],4))

    return fitFn,integral,roundParams

##################################################################################
# for x in [X_0] and y in MIP *only* - gamma distribution fit of entire ECal shower
def gammaFit(ecalEnd,graphn,totEnArray):
    fitFn=[0]*7
    integral=[0.]*7
    params=[0]*7
    roundParams=[0]*7
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
        params[i]=fitFn[i].GetParameters()
        integral[i]=fitFn[i].Integral(ecalEnd,xInt,1.e-12)/totEnArray[i]/(5.*0.26) #integrate up to y=10e-10 and divide by 5 mm*0.26 from extrapolating from radiation lengths of thick layers
        roundParams[i]=[0.]*3
        for j in range(3):
            roundParams[i][j]="{0:.4f}".format(round(params[i][j],4))

    return fitFn, integral,roundParams
#################################################################################
