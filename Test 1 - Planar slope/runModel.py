#!/usr/bin/python

import cPickle as pickle
import numpy
import time
import os
from subprocess import call
import tempfile
import sys
import fileIO


nExpansionIts=None; nSmoothingIts=None; expansionSlope=None;
goring.setPrecision32()

#############################################################################
# Control Panel

# Provide paths to Python modules and C++ library
sys.path.append('../') # sgrid.py
extLibName=r"../hydraulics_f32.so"# C++ library



parametersFile="params1km.pck"
dtmFileName="GB_SE.tiff"

outputDirectory="results/" # Must have final /
outputPrefix="100mm_12h"

rainfallDuration=24. # In hours
rainfallStart=120. # In hours
duration=rainfallStart+rainfallDuration*3

rainfallDepth=100. # In mm


initialWlFile=None  # This can be used to specify initial water depths from a
                    # csv file output by a previos run - use None to turn off

pcRunoff=35.       # Percentage runoff

noDataValue=-9999
noDataReplacement=-9999

saveMax=True # Set to true to save max water levels, flows etc
saveEnd=False # Set to true to save final water levels, flows etc

defaultDepth=1.0
flowThreshold=10. # Use this to switch off cells with flow below this value

# Flow points
flowPointsShp='Flow BC.shp'
flowAttr='flow'
flowMultiplier=1.0 # Easy way to adjust all values by this factor

# Baseflow - useful for groundwater contributions etc
baseFlow=0.1 # in m3/s/km2, introduced into all cells

# Water level boundary conditions
wlShp="Water Level BC.shp" # Water level BCs in shapefile - steady state only - use None to turn off
wlAttr='wlbc'              # Attribute holding water level

# Initial conditions
initialWL=2.0   # Initial water level applied everywhere - use None to turn off

initialTimeStep=30.  # Time step at start of run
minTimeStep=30.
maxTimeStep=3600.

#############################################################################

import sgrid

# Channel parameters
channel=False


flowPathOutput="new"
extendWlGrid=False

# Create output folder if it doesn't exist already
if not os.path.isdir(outputDirectory):
    assert not os.path.exists(outputDirectory)
    os.makedirs(outputDirectory)


arrayType=goring.getPrecision()

file=open(parametersFile)
(xll,yll,cellSize,xsz,ysz,convParX,convParY,storagePar)=pickle.load(file)
file.close()

(cppCalcFlow,cppCalcFlowGrid,cppDryCheck,cppTimeStep,cppWlFromVolGrid,\
    cppConveyanceParameters,cppMaxVolGrid,cppResample2,cppResample3, \
    cppFlowPaths,cppSum,cppCalcStorageParameters,cppLazyFlowPaths, \
    cppWlFill,cppBurnFlowPaths,cppMakeWlGrid,cppClipZero,cppDryCheckDiagnostic,
    cppScsAdditionalRunoff,cppCalcFlowEdges)=\
    goring.loadCppLib("./hydraulics_f32.so")

nNonNullCells=(storagePar[:,:,0]!=-9999.).sum()

runoffGrid=numpy.zeros((xsz,ysz),dtype=arrayType)
runoffGrid[:,:]=rainfallDepth*pcRunoff/100. # Can introduce spatial variation


wlGrid=numpy.zeros((xsz,ysz),dtype=arrayType)
wlGrid[:,:]=storagePar[:,:,0]
volGrid=numpy.zeros((xsz,ysz),dtype=arrayType)

# Initial water levels
wlGiven=False
if initialWL is not None:
    wlGrid[numpy.where(wlGrid<initialWL)]=initialWL
    wlGiven=True

if initialWlFile is not None: # Initial conditions in file
    wlGrid[:,:]=fileIO.readCSV(initialWlFile,"WL",xsz,ysz)
    wlGiven=True


if wlGiven:
    for i in range(xsz):
        for j in range(ysz):
            if storagePar[i,j,0]!=-9999.:
                volGrid[i,j]=goring.volFromWl(wlGrid[i,j],i,j,storagePar)

maxVolGrid=numpy.zeros((xsz,ysz),dtype=arrayType)

flowX=numpy.zeros((xsz+1,ysz),dtype=arrayType)
flowY=numpy.zeros((xsz,ysz+1),dtype=arrayType)
maxFlowX=numpy.zeros((xsz+1,ysz),dtype=arrayType)
maxFlowY=numpy.zeros((xsz,ysz+1),dtype=arrayType)

dryMask=numpy.zeros((xsz,ysz),dtype=arrayType)+1
dryThresh=0.1

print "Initial Volume=%e"%volGrid.sum()

# Read sources and process
if flowPointsShp is not None:
    flowPointsRaw=fileIO.readPointShapefile(flowPointsShp)

    flowPoints=[]
    flowPointsN=0
    for p in flowPointsRaw:
        xi=int((p[0]-xll)/cellSize)
        yi=int((p[1]-yll)/cellSize)
        qi=p[2][flowAttr]*flowMultiplier

        if xi>=0 and xi<xsz and yi>=0 and yi<ysz:
            flowPointsN+=1
            flowPoints.append((xi,yi,qi))

    if flowPointsN>0:
        flowPointsXi=numpy.zeros(flowPointsN,dtype=numpy.int32)
        flowPointsYi=numpy.zeros(flowPointsN,dtype=numpy.int32)
        flowPointsQ=numpy.zeros(flowPointsN,dtype=arrayType)

        for i in range(flowPointsN):
            flowPointsXi[i]=flowPoints[i][0]
            flowPointsYi[i]=flowPoints[i][1]
            flowPointsQ[i]=flowPoints[i][2]

else: # Null
    flowPoints=None
    flowPointsN=0
    flowPointsXi=numpy.array([0],dtype=numpy.int32)
    flowPointsYi=numpy.array([0],dtype=numpy.int32)
    flowPointsQ=numpy.array([0],dtype=arrayType)

# Read downstream BCs and process
if wlShp is not None:
    bcPointsRaw=fileIO.readPolylineShapefile(wlShp)
    nl=len(bcPointsRaw)
    wlPoints=[]
    wlPointsN=0
    for i in range(nl):
        xl=bcPointsRaw[i][0]
        yl=bcPointsRaw[i][1]
        wl=bcPointsRaw[i][2][wlAttr]

        for j in range(len(xl)-1):
            x1=xl[j]
            y1=yl[j]
            x2=xl[j+1]
            y2=yl[j+1]

            segmentLength=numpy.sqrt((x2-x1)**2+(y2-y1)**2)
            nSteps=int(2*segmentLength/cellSize)

            for k in range(nSteps):
                x=x1+(x2-x1)*k/nSteps
                y=y1+(y2-y1)*k/nSteps

                xi=int((x-xll)/cellSize)
                yi=int((y-yll)/cellSize)

#                print (xi,yi,wl)

                if xi>=0 and xi<xsz and yi>=0 and yi<ysz:
                    if (xi,yi,wl) not in wlPoints:
                        wlPointsN+=1
                        wlPoints.append((xi,yi,wl))
else:
    wlPointsN=0
    wlPoints=[]


# Perform timesteps
currentTime=0.
timeStep=initialTimeStep

nextDisplayTime=0.
displayInterval=3600.*duration/100.

totalRunoffInputVolume=0.
t1=time.time()

while currentTime<(duration*3600.):

    timeStep=goring.calcTimeStep(wlGrid,storagePar[:,:,0],cellSize)

    if numpy.isnan(timeStep) or numpy.isinf(timeStep):
        timeStep=initialTimeStep
    timeStep=max(timeStep,minTimeStep)
    timeStep=min(timeStep,maxTimeStep)

    # Add baseflow
    bf=baseFlow*timeStep
    volGrid+=bf

    # Add rainfall
    if currentTime>=rainfallStart*3600. and currentTime<(rainfallStart+rainfallDuration)*3600.:
        totalRainfallSoFarGrid=rainfallDepth*(currentTime/3600.-rainfallStart)/rainfallDuration

        totalRunoffInputVolume+=runoffGrid.sum()*\
            timeStep/(rainfallDuration*3600.)*cellSize*cellSize/1000.

        volGrid+=runoffGrid*timeStep/(rainfallDuration*3600.)*\
            cellSize*cellSize/1000.

    # Mask dry cells
    depthGrid=wlGrid-storagePar[:,:,0]
    dryMask[:,:]=1
    dryMask[numpy.where(depthGrid<dryThresh)]=0
    dryMask[numpy.where(storagePar[:,:,0]==-9999.)]=0

    totalActiveVol=volGrid[numpy.where(dryMask>0)].sum()

    cppCalcFlowGrid(wlGrid,convParX,convParY,\
        cellSize,timeStep,xsz,ysz,flowX,flowY,dryMask)

    # Calculate flow out of edge cells
    Qout=0.
    Qout=cppCalcFlowEdges(wlGrid,convParX,convParY,storagePar,cellSize,xsz,ysz,\
        flowX,flowY,dryMask)


    # Check for drying cells
    nDCI=cppDryCheck(volGrid,flowX,flowY,timeStep,xsz,ysz,\
        flowPointsXi,flowPointsYi,flowPointsQ,flowPointsN)

    # Update cell volumes
    cppTimeStep(volGrid,flowX,flowY,timeStep,xsz,ysz,\
        flowPointsXi,flowPointsYi,flowPointsQ,flowPointsN)

    volGrid[numpy.where(storagePar[:,:,0]==-9999)]=0.
    cppWlFromVolGrid(volGrid,wlGrid,storagePar,xsz,ysz,channel,cellSize)

    # Modify cell values with boundary water levels
    for wlPt in wlPoints:
        if storagePar[wlPt[0],wlPt[1],0]==-9999:
            continue

        wlGrid[wlPt[0],wlPt[1]]=wlPt[2]
        newV=goring.volFromWl(wlPt[2],wlPt[0],wlPt[1],storagePar)
        Qout+=(volGrid[wlPt[0],wlPt[1]]-newV)/timeStep

        volGrid[wlPt[0],wlPt[1]]=newV

    Qin=sum(flowPointsQ)+baseFlow*nNonNullCells


# Track maximum volumes
    if (currentTime/3600.)>=rainfallDuration:
        cppMaxVolGrid(volGrid,maxVolGrid,flowX,maxFlowX,flowY,maxFlowY,xsz,ysz)



    nActiveCells=(dryMask==1).sum()
    if (currentTime>=nextDisplayTime):

        if currentTime>0:
            pcComplete=float(currentTime)/(3600.*duration)
            pcToGo=1.-pcComplete
            t2=time.time()-t1
            projectedFinish=time.time()+pcToGo*(t2/pcComplete)
            projectedFinishString=time.strftime("%H:%M:%S",time.localtime(projectedFinish))
        else:
            projectedFinishString=' - '

        print "t=%s dt=%0.2f V=%e aV=%e rV=%e nDC=%i nActive=%i Qin=%0.1f Qout=%0.1f Finish=%s"\
            %(goring.formatTime(currentTime),timeStep,volGrid.sum(),
              totalActiveVol,totalRunoffInputVolume,nDCI,nActiveCells,Qin,Qout,
              projectedFinishString)

        nextDisplayTime+=displayInterval

    currentTime+=timeStep

print "Final Volume= %e"%volGrid.sum()

#cppWlFromVolGrid(maxVolGrid,wlGrid,storagePar,xsz,ysz,channel)
cppWlFromVolGrid(volGrid,wlGrid,storagePar,xsz,ysz,channel,cellSize)


t2=time.time()
print "Completed simulation in %0.2fs"%((t2-t1))

###########################################################################
# Save results

tmpFileName=tempfile._get_candidate_names().next()+'.tiff'

ogrCommand=['gdal_translate']
ogrCommand+=['-co','BIGTIFF=YES']
ogrCommand+=['-co','TILED=YES']
ogrCommand+=[dtmFileName]
ogrCommand+=[tmpFileName]


call(ogrCommand)

if saveEnd:
    print "Resampling and saving end depths/flows to file..."

    maskList=numpy.where((wlGrid-storagePar[:,:,0])<dryThresh)
    wlGrid[maskList]=storagePar[:,:,0][maskList]

    goring.saveResults(volGrid,wlGrid,flowX,flowY,storagePar,xsz,ysz,cellSize,xll,yll,
                       defaultDepth,flowThreshold,channel,flowPathOutput,
                       tmpFileName,noDataValue,noDataReplacement,
                       outputDirectory,outputPrefix,cppResample3,cppLazyFlowPaths,
                       extendWlGrid,cppBurnFlowPaths,cppMakeWlGrid,cppWlFill,cppClipZero)

if saveMax:
    print "Resampling and saving max depths/flows to file..."
    cppWlFromVolGrid(maxVolGrid,wlGrid,storagePar,xsz,ysz,channel,cellSize)

    maskList=numpy.where((wlGrid-storagePar[:,:,0])<dryThresh)
    wlGrid[maskList]=storagePar[:,:,0][maskList]

#    f=open("/owl4/Ambiental Sgrid/Australia/PYTHON/AWS bundle/t1.pck",'w')
#    pickle.dump((maxVolGrid,wlGrid,maxFlowX,maxFlowY,storagePar,xsz,ysz,cellSize,xll,yll,
#                       defaultDepth,flowThreshold,channel,flowPathOutput,
#                       tmpFileName,noDataValue,noDataReplacement,
#                       outputDirectory,outputPrefix+'_max',
#                       extendWlGrid),f)
#    f.close()

    goring.saveResults(maxVolGrid,wlGrid,maxFlowX,maxFlowY,storagePar,xsz,ysz,cellSize,xll,yll,
                       defaultDepth,flowThreshold,channel,flowPathOutput,
                       tmpFileName,noDataValue,noDataReplacement,
                       outputDirectory,outputPrefix+'_max',cppResample3,cppLazyFlowPaths,
                       extendWlGrid,cppBurnFlowPaths,cppMakeWlGrid,cppWlFill,cppClipZero)

os.remove(tmpFileName)

t2=time.time()

print "Completed all in %0.2fs"%((t2-t1))
