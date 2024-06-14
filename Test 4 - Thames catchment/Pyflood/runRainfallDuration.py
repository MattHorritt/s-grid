#!/usr/bin/python

import numpy
import time
import fileIO
import timeStep
import random
import os

import ctypes
from numpy.ctypeslib import ndpointer
import pyfloodUtils
import scsCnUtils
from subprocess import call
import tempfile
import sys

#############################################################################
# Control Panel

n=0.045
g=9.81

# Arguments: xmin xmax ymin ymax duration rainfall result_root
xmin=float(sys.argv[1])
xmax=float(sys.argv[2])
ymin=float(sys.argv[3])
ymax=float(sys.argv[4])

rainfallDuration=float(sys.argv[5])
rainfallDepth=float(sys.argv[6])

resultsPrefix=sys.argv[7]
diagnosticFileName=sys.argv[8]


dtmFileName=r"../DTM.tiff"
manningsN=0.06

noDataValue=-9999.
noDataReplacement=1e6

bcBurnValue=2.0
bcBurnAll=True

# Boundary conditions
flowLineShp=None # "/owl3/S-Grid vs Pyflood/GIS/lincsFlowLineInland.shp"
flowMultiplier=None # 0.225

# Rainfall
#rainfallDepth=100. # In mm
# pcRunoff=35

# Parameters
totalTime=rainfallDuration*3*3600.
#totalTime=16*60.

minTimeStep=1.0
defaultTimeStep=10.0
maxTimeStep=100.0

dt=minTimeStep

depthThresh=0.01

convDepthThresh=0.1
convTime=600.0
convTestOff=True

tsFactor=0.7

displayInterval=10.
#displayInterval=-1.
lastDisplayTime=0.

wlInit=-100.

#############################################################################

cnwOld=-1
prevConvTestTime=0

useExternalTimeStep=True

# Topography
z,dx,xsz,ysz,xll,yll=fileIO.readScalarGrid(dtmFileName,window=[xmin,xmax,ymin,ymax])

# No non null cells, so exit nicely
if z.max()==noDataValue:

    print "NULL TILE"

    f=open(diagnosticFileName,'a')
    f.write("%s\tNULL TILE"%(resultsPrefix))

    for arg in sys.argv:
        f.write("\t%s"%arg)

    f.write("\n")
    f.close()


    sys.exit()


# Roughness - assume already sampled to same extents/resolution as DTM
manningsGrid=numpy.zeros((xsz,ysz))+manningsN


zNotNullList=numpy.where(z!=noDataValue)
z[numpy.where(z==noDataValue)]=noDataReplacement

###########################################################################
# Calculate total rainfall and runoff for info and mass balance
print "Calculating total runoff..."

runoffGrid=numpy.zeros((xsz,ysz))+rainfallDepth
additionalRunoff=numpy.zeros((xsz,ysz))
totalRainfallSoFarGrid=numpy.zeros((xsz,ysz))

print "Average rainfall=",runoffGrid.sum()/(xsz*ysz), "mm"
print "Average runoff=", runoffGrid.sum()/(xsz*ysz), "mm"
print "Total runoff input= %e m3"%(runoffGrid.sum()*dx*dx/1000.)
###########################################################################


# Set up grids
xCentre=numpy.array(range(xsz))*dx+dx/2
yCentre=numpy.array(range(ysz))*dx+dx/2
x=numpy.array(range(xsz))*dx+dx/2
y=numpy.array(range(ysz))*dx+dx/2

h=numpy.zeros((xsz,ysz))
hMax=numpy.zeros((xsz,ysz))
absVel=numpy.zeros((xsz,ysz))
u=numpy.zeros((xsz+1,ysz))
v=numpy.zeros((xsz,ysz+1))
hnew=numpy.zeros((xsz,ysz))
unew=numpy.zeros((xsz+1,ysz))
vnew=numpy.zeros((xsz,ysz+1))
courantNumber=numpy.zeros((xsz,ysz))

cellActive=numpy.zeros((xsz,ysz),dtype=numpy.bool)

# Initial values
h[:,:]=wlInit-z[:,:]
h[numpy.where(h<=0)]=0.

cellActive[numpy.where(h>depthThresh)]=True
initialVolume=h.sum()*dx*dx

print "Initial Volume=", initialVolume

prevTotalVolume=initialVolume


integratedQ=0.
integratedQout=0.
rainfallQ=0.

oldWetCount=-1

u[:,:]=0.0
v[:,:]=0.0
hnew[:]=h

if useExternalTimeStep:
    libPath="./C++/"
    libFileName="timeStep.so"
    tmpLibFileName="%i"%int(random.uniform(1000000,9999999))+".so"
    os.system("cp \""+libPath+libFileName+"\" \""+libPath+tmpLibFileName+"\"")

    lib=ctypes.cdll.LoadLibrary(libPath+tmpLibFileName)
    timeStepFunc=lib.timeStep

    timeStepFunc.argtypes=[ndpointer(ctypes.c_double),\
        ndpointer(ctypes.c_double),\
        ndpointer(ctypes.c_double),\
        ndpointer(ctypes.c_double),\
        ndpointer(ctypes.c_double),\
        ndpointer(ctypes.c_double),\
        ndpointer(ctypes.c_double),\
        ndpointer(ctypes.c_double),\
        ndpointer(ctypes.c_int),\
        ndpointer(ctypes.c_int),\
        ctypes.c_int,\
        ndpointer(ctypes.c_double),\
        ndpointer(ctypes.c_int),\
        ndpointer(ctypes.c_int),\
        ctypes.c_int,\
        ndpointer(ctypes.c_double),\
        ndpointer(ctypes.c_int),\
        ndpointer(ctypes.c_int),\
        ctypes.c_int,\
        ndpointer(ctypes.c_double),\
        ndpointer(ctypes.c_int),\
        ndpointer(ctypes.c_int),\
        ctypes.c_int,\
        ctypes.c_int,\
        ctypes.c_int,\
        ctypes.c_double,\
        ctypes.c_double,\
        ctypes.CFUNCTYPE(ctypes.c_double),\
        ctypes.CFUNCTYPE(None,ctypes.c_double,ctypes.c_double),\
        ctypes.CFUNCTYPE(None,ctypes.c_double,ctypes.c_double),\
        ctypes.CFUNCTYPE(None,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double),\
        ctypes.CFUNCTYPE(ctypes.c_bool,ctypes.c_double),\
        ctypes.c_double,\
        ndpointer(ctypes.c_double),\
        ctypes.c_int]

    lib.totaliser.restype=ctypes.c_double
    lib.totaliser.argtypes=[ndpointer(ctypes.c_double),ctypes.c_int,ctypes.c_int]

    lib.courantNumber.argtypes=[ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),\
        ctypes.c_int,ctypes.c_int,\
        ctypes.c_double,ctypes.c_double,
        ndpointer(ctypes.c_double)]

#    print "Lib speed test:"
#    randomArray=numpy.random.random((xsz,ysz))
#    t1=time.time()
#    vol=lib.totaliser(randomArray,xsz,ysz)*dx*dx
#    print "External library: ", vol, " in ",time.time()-t1,"s"
#    t1=time.time()
#    vol=randomArray.sum()*dx*dx
#    print "Numpy library: ", vol, " in ",time.time()-t1,"s"

else:
    timeStepFunc=timeStep.timeStep

crZeroTime=0.

def calcTimeStep():
    global dt, crZeroTime

    lib.courantNumber(hnew,unew,vnew,xsz,ysz,dt,dx,courantNumber)
    cmm=courantNumber.max()
    if cmm!=0:
        dt=tsFactor*dt/cmm
    else:
        dt=defaultTimeStep

    if numpy.isinf(dt):
        dt=defaultTimeStep

    dt=min(dt,maxTimeStep)
    dt=max(dt,minTimeStep)


    return float(dt)

def updateBCs(t,dt):
    global sbcVals, integratedQ

    # Calculate overtopping inflows
    f=numpy.sin(2*numpy.pi*t/(12.5*3600.))
    if t/3600.>6.25:
        f=0.

    Q=0.

    for i, pot in enumerate(potrList):
        sbcVals[i]=f*pot
        Q+=f*pot

    integratedQ+=sum(sbcVals)*dt


    return

#tInScs=0.0

def rainfall(t,dt):
    global h, integratedQ, rainfallQ, zNotNullList, tInScs

    rainfallQ=0.

    if t<=(3600.*rainfallDuration):

        totalRainfallSoFarGrid=runoffGrid*(t/3600)/rainfallDuration
        additionalRainfallGrid=(runoffGrid/rainfallDuration)*(dt/3600.)

#        runoffGrid=scsCnUtils.scsAdditionalRunoffGrid(\
#                    totalRainfallSoFarGrid,additionalRainfallGrid,cnGrid)

#        for i in range(xsz):
#            for j in range(ysz):
#                runoffGrid[i,j]=scsCnUtils.scsAdditionalRunoff(\
#                    totalRainfallSoFarGrid[i,j],additionalRainfallGrid[i,j],cnGrid[i,j])

#        tInScs+=(time.time()-tTmp)

#        print runoffGrid[0,0]

#        sys.exit()

        dDepth=dt*(runoffGrid/(3600.*rainfallDuration))/1000.

#        print "In rainfall():", dDepth.min(), dDepth.max(), dDepth.sum()/(xsz*ysz)

        h[:,:][zNotNullList]+=dDepth[zNotNullList]
        rainfallQ=dDepth.sum()*dx*dx/dt
        integratedQ+=rainfallQ*dt

    return

def convergenceTest(currentTime):
    global cnwOld, hMax, prevConvTestTime, convTestOff, totalTime

    if convTestOff:
        return False

    if (currentTime-prevConvTestTime)<convTime:
        return False

    prevConvTestTime=currentTime

    hMax=numpy.maximum(hMax,h)

    cnw=(hMax>convDepthThresh).sum()

    stopNow=(cnw==cnwOld and currentTime>(3600.*rainfallDuration))

#    print cnw, cnwOld, currentTime, rainfallDuration, stopNow

    cnwOld=cnw

    if stopNow:
        totalTime=currentTime

    return stopNow

#    global hMax, oldWetCount
#
#
#    updateCount=(h>hMax).sum()
#
#    updateDepth=(h-hMax).sum()
#
#    hMax=numpy.maximum(hMax,h)
#
#    wetCount=(hMax>depthThresh).sum()
#
#    updateDepth/=wetCount
#
#    stopNow=float(updateCount)/max(0.1,wetCount)<(pcStopCriterion/100.) and currentTime>rainfallDuration*3600.
#
#    if stopNow or True:
#        print "Model stopping:", updateCount,wetCount, currentTime, updateDepth
#
#    oldWetCount=wetCount

    return False

def display(ts,currentTime,Qout,intQout): # Use this to keep track of maximums too
    global totalVolume, prevTotalVolume, lastDisplayTime, hMax, integratedQout
    totalVolume=h.sum()*dx*dx

    Qin=sum(sbcVals)+rainfallQ

    if currentTime>lastDisplayTime:
        dVdt=(totalVolume-prevTotalVolume)/(currentTime-lastDisplayTime)
    else:
        dVdt=0.

    wetCount=(h>depthThresh).sum()
#    cnw1=(h>convDepthThresh).sum()
#    cnw2=(hMax>convDepthThresh).sum()

#    print ts,pyfloodUtils.formatTime(currentTime),dt,h.max(),totalVolume,Qin,Qout,dVdt,wetCount
#    assert False

    print "Time step: %5i Time: %s dt: %5.2f Max H: %5.2f Vol: %5.3e Qin: %6.1f Qout: %6.1f dV/dt: %6.1f nWet=%i"\
        %(ts,pyfloodUtils.formatTime(currentTime),dt,h.max(),totalVolume,Qin,Qout,dVdt,wetCount)

#    print cnw1, cnw2


#    mbError=((totalVolume-prevTotalVolume)/(currentTime-lastDisplayTime)+Qout-rainfallQ)
#    print totalVolume
#    print prevTotalVolume
#    print rainfallQ
#    print dt
#    print currentTime-lastDisplayTime
#    print intQout
#    print "mbError=", mbError

    prevTotalVolume=totalVolume
    lastDisplayTime=currentTime



    hMax=numpy.maximum(hMax,h)

    integratedQout=intQout

    return


# Iterate
nIt=int(totalTime/dt)
maxIt=1000000
lastDisplayTime=0.

t1=time.time()

lastDisplayTime=0.
currentTime=0.

# Set up BC arrays and fill before first call to C++ lib
ubcVals=numpy.array([])
ubcLocX=numpy.array([],numpy.int32)
ubcLocY=numpy.array([],numpy.int32)
ubcCount=0

vbcVals=numpy.array([])
vbcLocX=numpy.array([],numpy.int32)
vbcLocY=numpy.array([],numpy.int32)
vbcCount=0

hbcVals=numpy.array([])
hbcLocX=numpy.array([],numpy.int32)
hbcLocY=numpy.array([],numpy.int32)
hbcCount=0

# Source boundary conditions - peak overtopping rate
flowPointsN,flowPointsXi,flowPointsYi,flowPointsMaxQ=\
    pyfloodUtils.processFlowLineBCs(flowLineShp,xll,yll,xsz,ysz,dx)

potrList=[]
ptList=[]
sbcCount=0

# Remove duplicates and apply flow multiplier
for i in range(flowPointsN):
    if (flowPointsXi[i],flowPointsYi[i]) not in ptList:
        ptList.append((flowPointsXi[i],flowPointsYi[i]))
        potrList.append(flowPointsMaxQ[i]*flowMultiplier)
        sbcCount+=1
    else:
        idx=ptList.index((flowPointsXi[i],flowPointsYi[i]))
        potrList[idx]+=flowPointsMaxQ[i]*flowMultiplier

if flowPointsN>0:
    xiList,yiList=zip(*ptList)
else:
    xiList=[]
    yiList=[]

sbcVals=numpy.array([0.]*sbcCount)
sbcLocX=numpy.array(xiList,numpy.int32)
sbcLocY=numpy.array(yiList,numpy.int32)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Check z levels for BCs - burn in to make sure we're not adding flow to
# NaN or high area
for i in range(sbcCount):
    if bcBurnAll:
        z[xiList[i], yiList[i]]=bcBurnValue
    elif z[xiList[i], yiList[i]]==noDataValue or z[xiList[i], yiList[i]]==noDataReplacement:
        z[xiList[i], yiList[i]]=bcBurnValue

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

updateBCs(0.,dt) # Make sure our BCs are filled in before we send to C++ lib

print "Entering time step code..."
timeStepFunc(u,v,h,unew,vnew,hnew,z,\
    ubcVals,ubcLocX,ubcLocY,ubcCount,\
    vbcVals,vbcLocX,vbcLocY,vbcCount,\
    hbcVals,hbcLocX,hbcLocY,hbcCount,\
    sbcVals,sbcLocX,sbcLocY,sbcCount,\
    xsz,ysz,totalTime,displayInterval,\
    ctypes.CFUNCTYPE(ctypes.c_double)(calcTimeStep),\
    ctypes.CFUNCTYPE(None,ctypes.c_double,ctypes.c_double)(updateBCs),\
    ctypes.CFUNCTYPE(None,ctypes.c_double,ctypes.c_double)(rainfall),\
    ctypes.CFUNCTYPE(None,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double)(display),\
    ctypes.CFUNCTYPE(ctypes.c_bool,ctypes.c_double)(convergenceTest),\
    dx,manningsGrid,0)

t2=time.time()

wl=h+z
wl[numpy.where(h<depthThresh)]=-9999.

for i in range(xsz):
    for j in range(ysz):
        absVel[i,j]=numpy.sqrt((0.5*(unew[i,j]+unew[i+1,j]))**2+(0.5*(vnew[i,j]+vnew[i,j+1]))**2)

totalHead=wl+absVel**2/(2*g)
totalHead[numpy.where(h<depthThresh)]=-9999.

#fileIO.saveScalarGrid(h,xll,yll,dx,resultsPrefix+".depth.tiff")
fileIO.saveScalarGrid(hMax,xll,yll,dx,resultsPrefix+".maxdepth.tiff")
#fileIO.saveScalarGrid(z,xll,yll,dx,resultsPrefix+".topo.tiff")
#fileIO.saveScalarGrid(wl,xll,yll,dx,resultsPrefix+".wl.tiff")
#fileIO.saveScalarGrid(totalHead,xll,yll,dx,resultsPrefix+".th.tiff")
#fileIO.saveScalarGrid(courantNumber,xll,yll,dx,resultsPrefix+".Cr.tiff")
#
#fileIO.saveScalarCSV(h,xll,yll,dx,resultsPrefix+".depth.csv")
#fileIO.saveVectorCSV(u,v,xll,yll,dx,resultsPrefix+".vel.csv")
#fileIO.saveVectorCSV(advCrossTermX,advCrossTermY,0,0,dx,resultsPrefix+".adv.csv")

finalVolume=h.sum()*dx*dx
print "Initial Volume:    ", initialVolume
print "Final Volume:      ", finalVolume
print "Change in volume:  ", finalVolume-initialVolume
print "Total Flow in:     ", integratedQ
print "Total Flow out:    ", integratedQout
print "Mass balance error:", (finalVolume-initialVolume-integratedQ+integratedQout)
print "Mass balance error:", 100*(finalVolume-initialVolume-integratedQ+integratedQout)/max(finalVolume,integratedQ),"%"
print "Run for:           ", totalTime/3600., "h"

#print tInScs

f=open(diagnosticFileName,'a')
f.write("%s\t%f\t%f\t%f"%(resultsPrefix,t2-t1, \
    100*(finalVolume-initialVolume-integratedQ+integratedQout)/max(finalVolume,integratedQ),\
    totalTime/3600))

for arg in sys.argv:
    f.write("\t%s"%arg)

f.write("\n")
f.close()

#for i in range(2,40,2):
#    print "%f,%f,%f,%f"%(x[i],z[i,20],wl[i,20],absVel[i,20])


print "Completed in %0.2fs"%((t2-t1))

#del timeStepFunc
if useExternalTimeStep:
    os.system("rm \""+libPath+tmpLibFileName+"\"")
