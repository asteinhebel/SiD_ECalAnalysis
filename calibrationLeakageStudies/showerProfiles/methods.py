from pyLCIO import IOIMPL, EVENT,UTIL
from pyLCIO.io.LcioReader import LcioReader
from array import array
from ROOT import *
import math


################################################
#calculate solid angle from hit to initial particle gun trajectory
def solidAngle(phi,hit):
    angle=(math.acos((math.cos(phi*math.pi/180.)*hit.getPosition()[0]+math.sin(phi*math.pi/180.)*hit.getPosition()[1])/math.sqrt(hit.getPosition()[0]*hit.getPosition()[0]+hit.getPosition()[1]*hit.getPosition()[1]+hit.getPosition()[2]*hit.getPosition()[2])))
    #correct for arccos range
    if hit.getPosition()[1]<0:
        angle=-angle 

    return angle

################################################
#calculate phi angle
def phiAngle(phi,hit):
    angle=math.atan(hit.getPosition()[1]/hit.getPosition()[0])*(180./math.pi)
    #correct for arctan range
    if hit.getPosition()[0]<0 and hit.getPosition()[1]>0:
        angle=angle+180.
    elif hit.getPosition()[0]<0 and hit.getPosition()[1]<0:
        angle=angle-180.

    return angle

################################################
