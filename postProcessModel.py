#!/usr/bin/python

#############################################################################
# Control Panel
parametersFile="params.pck"

# Topography
dtmFileName=r"DTM.tiff"
useTempTopoFile=False # Use this to create uncompressed, tiled topo file to speed up access

# These values can be used to replace NULLs (e.g. at sea) with sensible values
noDataValue=None
noDataReplacement=None

# Folder and first part of filename where CSV outputs from runModel are stored
resultsDirectory="results"
resultsPrefix="outputs"

processMax=True # Set to true to save max water levels, flows etc
processEnd=True # Set to true to save final water levels, flows etc

defaultDepth=1.0    # Depth burnt into flow paths
flowThreshold=10.   # Use this to switch off interpolation between cells with
                    # flows below this value

# Reservoirs present?
reservoirsPresent=False
###########################################################################

import os
import cPickle as pickle
import numpy

import reservoirs
import fileIO
import sgrid

# Path to C++ library
extLibName=r"sgridHydraulics.dll"# C++ library

sgrid.setPrecision32()

file=open(parametersFile)
(xll,yll,cellSize,xsz,ysz,convParX,convParY,storagePar)=pickle.load(file)
file.close()

(cppCalcFlow,cppCalcFlowGrid,cppDryCheck,cppTimeStep,cppWlFromVolGrid,\
    cppConveyanceParameters,cppMaxVolGrid,cppResample2,cppResample3, \
    cppFlowPaths,cppSum,cppCalcStorageParameters,cppLazyFlowPaths, \
    cppWlFill,cppBurnFlowPaths,cppMakeWlGrid,cppClipZero,cppDryCheckDiagnostic,
    cppScsAdditionalRunoff,cppCalcFlowEdges,cppLicenceCheck)=\
    sgrid.loadCppLib(extLibName)

if useTempTopoFile:
    tmpDtmFileName=sgrid.uncompressGeoTiff(dtmFileName,tiled=True)
else:
    tmpDtmFileName=dtmFileName

dryThresh=0.1
channel=False

flowPathOutput="new"
extendWlGrid=False

if reservoirsPresent:
    fileName=os.path.join(resultsDirectory,resultsPrefix+"_reservoirs.pck")
    with open(fileName) as f:
        (fillCellDict,zeroPolyList)=pickle.load(f)
else:
    fillCellDict=None
    zeroPolyList=None

if processEnd:
    wlFileName=os.path.join(resultsDirectory,resultsPrefix+'_wl.csv')
    flowFileName=os.path.join(resultsDirectory,resultsPrefix+'_flow.csv')

    wlGrid=fileIO.readCSV(wlFileName,"WL",xsz,ysz,dataType=sgrid.getPrecision())
    flowX,flowY=fileIO.readFlowCsv(flowFileName,xsz,ysz,dataType=sgrid.getPrecision())

    print "Resampling and saving end depths/flows to file..."

    maskList=numpy.where((wlGrid-storagePar[:,:,0])<dryThresh)
    wlGrid[maskList]=storagePar[:,:,0][maskList]

    if fillCellDict is not None:
        reservoirs.fillEdgeCells(fillCellDict,wlGrid)

    sgrid.saveResults(wlGrid,wlGrid,flowX,flowY,storagePar,xsz,ysz,cellSize,xll,yll,
                       defaultDepth,flowThreshold,channel,flowPathOutput,
                       tmpDtmFileName,noDataValue,noDataReplacement,
                       resultsDirectory,resultsPrefix,cppResample3,cppLazyFlowPaths,
                       extendWlGrid,cppBurnFlowPaths,cppMakeWlGrid,cppWlFill,cppClipZero,
                       saveCsv=False,zeroPolyList=zeroPolyList)

if processMax:
    wlFileName=os.path.join(resultsDirectory,resultsPrefix+'_max_wl.csv')
    flowFileName=os.path.join(resultsDirectory,resultsPrefix+'_max_flow.csv')

    wlGrid=fileIO.readCSV(wlFileName,"WL",xsz,ysz,dataType=sgrid.getPrecision())
    flowX,flowY=fileIO.readFlowCsv(flowFileName,xsz,ysz,dataType=sgrid.getPrecision())

    print "Resampling and saving max depths/flows to file..."

    maskList=numpy.where((wlGrid-storagePar[:,:,0])<dryThresh)
    wlGrid[maskList]=storagePar[:,:,0][maskList]

    if fillCellDict is not None:
        reservoirs.fillEdgeCells(fillCellDict,wlGrid)

    sgrid.saveResults(wlGrid,wlGrid,flowX,flowY,storagePar,xsz,ysz,cellSize,xll,yll,
                       defaultDepth,flowThreshold,channel,flowPathOutput,
                       tmpDtmFileName,noDataValue,noDataReplacement,
                       resultsDirectory,resultsPrefix+'_max',cppResample3,cppLazyFlowPaths,
                       extendWlGrid,cppBurnFlowPaths,cppMakeWlGrid,cppWlFill,cppClipZero,
                       saveCsv=False,zeroPolyList=zeroPolyList)

if useTempTopoFile:
    os.remove(tmpDtmFileName)


