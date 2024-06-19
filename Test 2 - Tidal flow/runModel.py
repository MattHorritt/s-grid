#!/usr/bin/python

#############################################################################
# Control Panel

# Provide paths to Python modules and dll etc
sgridPath = "../"


parametersFile="params.pck"
dtmFileName="dtm.tiff"

outputDirectory="results"
outputPrefix="output"

# Rainfall and runoff
rainfallDuration=0. # In hours
rainfallStart=0. # In hours
duration=72
rainfallDepth=0. # In mm
pcRunoff=35.       # Percentage runoff

# Initial conditions
initialWlFile=None  # This can be used to specify initial water depths from a
                    # csv file output by a previos run - use None to turn off
initialWL=0.0

# Water level boundary conditions
wlShp="wlBCs.shp" # Water level BCs in shapefile
wlAttr='wl'              # Attribute holding water level

def tidalWl(time):
    return math.sin(2*math.pi*time/43200.)

wlTimeSeriesFunction=tidalWl

# Flow points
flowPointsShp=None
flowAttr='flow'
flowMultiplier=1.0 # Easy way to adjust all values by this factor

noDataValue=None
noDataReplacement=None

saveMax=True # Set to true to save max water levels, flows etc
saveEnd=True # Set to true to save final water levels, flows etc

defaultDepth=1.0
flowThreshold=10. # Use this to switch off cells with flow below this value


# Baseflow - useful for groundwater contributions etc
baseFlow=0.0 # in m3/s/km2, introduced into all cells



initialTimeStep=30.  # Time step at start of run
minTimeStep=30.
maxTimeStep=3600.

#############################################################################
import sys
sys.path.append(sgridPath) # sgrid.py etc

import pickle
import numpy
import time
import os

import sgrid
import fileIO
import platform

system = platform.uname().system
if system == 'Linux':
    extLibName=r"../sgridHydraulics.so"# C++ library
elif system == 'Windows':
    extLibName=r"../sgridHydraulics.dll"# C++ library
else:
    sys.exit(f"S-Grid not implemented for operating system {system}")



sgrid.setPrecision32()

# Channel parameters
channel=False


flowPathOutput="new"
extendWlGrid=False

# Create output folder if it doesn't exist already
if not os.path.isdir(outputDirectory):
    assert not os.path.exists(outputDirectory)
    os.makedirs(outputDirectory)


arrayType=sgrid.getPrecision()

file=open(parametersFile,'rb')
(xll,yll,cellSize,xsz,ysz,convParX,convParY,storagePar)=pickle.load(file)
file.close()

(cppCalcFlow,cppCalcFlowGrid,cppDryCheck,cppTimeStep,cppWlFromVolGrid,\
    cppConveyanceParameters,cppMaxVolGrid,cppResample2,cppResample3, \
    cppFlowPaths,cppSum,cppCalcStorageParameters,cppLazyFlowPaths, \
    cppWlFill,cppBurnFlowPaths,cppMakeWlGrid,cppClipZero,cppDryCheckDiagnostic,
    cppScsAdditionalRunoff,cppCalcFlowEdges,cppCheckLicence)=\
    sgrid.loadCppLib(extLibName)

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
                volGrid[i,j]=sgrid.volFromWl(wlGrid[i,j],i,j,storagePar,cellSize)

maxVolGrid=numpy.zeros((xsz,ysz),dtype=arrayType)

flowX=numpy.zeros((xsz+1,ysz),dtype=arrayType)
flowY=numpy.zeros((xsz,ysz+1),dtype=arrayType)
maxFlowX=numpy.zeros((xsz+1,ysz),dtype=arrayType)
maxFlowY=numpy.zeros((xsz,ysz+1),dtype=arrayType)

dryMask=numpy.zeros((xsz,ysz),dtype=arrayType)+1
dryThresh=0.1

print("Initial Volume=%e"%volGrid.sum())

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

    timeStep=sgrid.calcTimeStep(wlGrid,storagePar[:,:,0],cellSize)

    if numpy.isnan(timeStep) or numpy.isinf(timeStep):
        timeStep=initialTimeStep
    timeStep=max(timeStep,minTimeStep)
    timeStep=min(timeStep,maxTimeStep)

    # Add baseflow
    bf=baseFlow*timeStep*cellSize*cellSize/1e6
    volGrid+=bf

    ############################################################################################
    # Add rainfall
    # Edit this to add different rainfall profiles, spatial variation, runoff etc
    if currentTime>=rainfallStart*3600. and currentTime<(rainfallStart+rainfallDuration)*3600.:
        # totalRunoffInputVolume (cumulative) used for mass balance tracking
        totalRunoffInputVolume+=runoffGrid.sum()*\
            timeStep/(rainfallDuration*3600.)*cellSize*cellSize/1000.

        # Rainfall is added as a volume to each cell
        volGrid+=runoffGrid*timeStep/(rainfallDuration*3600.)*\
            cellSize*cellSize/1000.

    ############################################################################################

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

    ############################################################################################
    # Modify flow boundary
    # Insert code here to modify values in the flowPointsQ array to represent a flow hydrograph
    ############################################################################################


    # Check for drying cells
    nDCI=cppDryCheck(volGrid,flowX,flowY,timeStep,xsz,ysz,\
        flowPointsXi,flowPointsYi,flowPointsQ,flowPointsN)

    # Update cell volumes
    cppTimeStep(volGrid,flowX,flowY,timeStep,xsz,ysz,\
        flowPointsXi,flowPointsYi,flowPointsQ,flowPointsN)

    volGrid[numpy.where(storagePar[:,:,0]==-9999)]=0.
    cppWlFromVolGrid(volGrid,wlGrid,storagePar,xsz,ysz,channel,cellSize)

    ############################################################################################
    # Apply water level boundary
    # wlPoints is list of [i,j,wl] values processed from shapefile for entering into model grid
    for wlPt in wlPoints:
        if storagePar[wlPt[0],wlPt[1],0]==-9999: # This indicates NULL cell so do nothing
            continue

        wlGrid[wlPt[0],wlPt[1]]=wlPt[2] # Modify this for time varying boundary conditions

        newV=sgrid.volFromWl(wlPt[2],wlPt[0],wlPt[1],storagePar,cellSize)
        Qout+=(volGrid[wlPt[0],wlPt[1]]-newV)/timeStep

        volGrid[wlPt[0],wlPt[1]]=newV
    ############################################################################################



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

        print("t=%s dt=%0.2f V=%e aV=%e rV=%e nDC=%i nActive=%i Qin=%0.1f Qout=%0.1f Finish=%s"\
            %(sgrid.formatTime(currentTime),timeStep,volGrid.sum(),
              totalActiveVol,totalRunoffInputVolume,nDCI,nActiveCells,Qin,Qout,
              projectedFinishString))

        nextDisplayTime+=displayInterval

    currentTime+=timeStep

print("Final Volume= %e"%volGrid.sum())

#cppWlFromVolGrid(maxVolGrid,wlGrid,storagePar,xsz,ysz,channel)
cppWlFromVolGrid(volGrid,wlGrid,storagePar,xsz,ysz,channel,cellSize)


t2=time.time()
print("Completed simulation in %0.2fs"%((t2-t1)))

###########################################################################
# Save results

if saveEnd:
    fileIO.saveVectorCSV(flowX,flowY,xll,yll,cellSize,\
        os.path.join(outputDirectory,outputPrefix+"_flow.csv"),thresholdVal=1e-3)

    fileIO.saveScalarCSV(wlGrid,xll,yll,cellSize,\
        os.path.join(outputDirectory,outputPrefix+"_wl.csv"), headerList=['WL'])


if saveMax:
    cppWlFromVolGrid(maxVolGrid,wlGrid,storagePar,xsz,ysz,channel,cellSize)

    maskList=numpy.where((wlGrid-storagePar[:,:,0])<dryThresh)
    wlGrid[maskList]=storagePar[:,:,0][maskList]

    fileIO.saveVectorCSV(maxFlowX,maxFlowY,xll,yll,cellSize,\
        os.path.join(outputDirectory,outputPrefix+"_max_flow.csv"),thresholdVal=1e-3)

    fileIO.saveScalarCSV(wlGrid,xll,yll,cellSize,\
        os.path.join(outputDirectory,outputPrefix+"_max_wl.csv"), headerList=['WL'])

t2=time.time()

print("Completed all in %0.2fs"%((t2-t1)))
