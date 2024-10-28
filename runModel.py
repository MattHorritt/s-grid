#!/usr/bin/python

import sys

#############################################################################
# Control Panel

# Get return period and SLR from arguments
seaLevelRise = float(sys.argv[1])
returnPeriod = float(sys.argv[2])

parametersFile="50m_params.pck" # Output from buildModel.py

outputDirectory="50m_results" # Folder for results

slrStr = f"{seaLevelRise:.1f}"
slrStr = slrStr.replace('.', 'p')
outputPrefix=f"output_slr{slrStr}_rp{returnPeriod:.0f}"  # Filename for results is based on return period and sea level

# Rainfall information and length of simulation - in hours
rainfallDuration=0.
rainfallStart=0.
duration=31.05 # Takes us to ebb tide after first peak
rainfallDepth=0. # In mm

initialWlFile=None  # This can be used to specify initial water depths from a
                    # csv file output by a previos run - use None to turn off

initialWL = 0.0

pcRunoff=100.       # Percentage runoff

saveMax=True # Set to true to save max water levels, flows etc
saveEnd=True # Set to true to save final water levels, flows etc

# Baseflow - useful for groundwater contributions etc
baseFlow=0.0 # in m3/s/km2, introduced into all cells

initialTimeStep=30.  # Time step at start of run
minTimeStep=30.
maxTimeStep=3600.

seaLevelRise = 0.0
defenceLineStr = r"PG: host=localhost dbname=ltis2025 active_schema=slr user=postgres password=postgres"
layerName = "defence_toe_levels"
# See comment "Add rainfall" for where to edit rainfall/runoff code
# See comment "Apply water level boundary" for where to edit water level boundary conditions
# See comment "Modify flow boundary" for where to edit flow boundary conditions


#############################################################################

import pickle
import numpy
import time
import os

import sgrid
import ltis_slr
import fileIO

# Path to C++ library
extLibName=r"./sgridHydraulics.so"# C++ library


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

file=open(parametersFile, 'rb')
(xll,yll,cellSize,xsz,ysz,convParX,convParY,storagePar)=pickle.load(file)
file.close()

(cppCalcFlow,cppCalcFlowGrid,cppDryCheck,cppTimeStep,cppWlFromVolGrid,\
    cppConveyanceParameters,cppMaxVolGrid,cppResample2,cppResample3, \
    cppFlowPaths,cppSum,cppCalcStorageParameters,cppLazyFlowPaths, \
    cppWlFill,cppBurnFlowPaths,cppMakeWlGrid,cppClipZero,cppDryCheckDiagnostic,
    cppScsAdditionalRunoff,cppCalcFlowEdges,cppLicenceCheck)=\
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

# Tide level points for LTIS SLR
wlPoints = ltis_slr.slr_tide_points(seaLevelRise, returnPeriod, defenceLineStr, xll, yll, xsz, ysz, cellSize,
                                    outputDirectory, outputPrefix, layerName = layerName)

# Some dummy arrays because we have no flow points
flowPointsN=0
flowPointsXi=numpy.array([0],dtype=numpy.int32)
flowPointsYi=numpy.array([0],dtype=numpy.int32)
flowPointsQ=numpy.array([0],dtype=arrayType)

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
    Qin, Qout = ltis_slr.applyTideLevels(wlPoints, currentTime, timeStep, seaLevelRise,
                             storagePar, wlGrid, volGrid, cellSize)

    ############################################################################################

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
