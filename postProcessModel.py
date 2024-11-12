#!/usr/bin/python

#############################################################################
# Control Panel
parametersFile="2m_EA_params.pck"

# Topography
dtmFileName=r"/merlin1/Projects/LTIS SLR/GIS/DTM/All_clip_range.tif"
useTempTopoFile=False # Use this to create uncompressed, tiled topo file to speed up access

# These values can be used to replace NULLs (e.g. at sea) with sensible values
noDataValue=None
noDataReplacement=None

# Folder and first part of filename where CSV outputs from buildModel are stored
resultsDirectory="2m_EA_results"
resultsPrefix="output_slr1p0_rp200"

processMax=True # Set to true to save max water levels, flows etc
processEnd=False # Set to true to save final water levels, flows etc

defaultDepth=1.0    # Depth burnt into flow paths
flowThreshold=1.0   # Use this to switch off interpolation between cells with
                    # flows below this value

###########################################################################

import os
import pickle
import numpy

import fileIO
import sgrid

# Path to C++ library
extLibName=r"./sgridHydraulics.so"# C++ library

sgrid.setPrecision32()

file=open(parametersFile, 'rb')
(xll,yll,cellSize,xsz,ysz,convParX,convParY,storagePar)=pickle.load(file)
file.close()

(cppCalcFlow,cppCalcFlowGrid,cppDryCheck,cppTimeStep,cppWlFromVolGrid,\
    cppConveyanceParameters,cppMaxVolGrid,cppResample2,cppResample3, \
    cppFlowPaths,cppSum,cppCalcStorageParameters,cppLazyFlowPaths, \
    cppWlFill,cppBurnFlowPaths,cppMakeWlGrid,cppClipZero,cppDryCheckDiagnostic,
    cppScsAdditionalRunoff,cppCalcFlowEdges,cppCheckLicence)=\
    sgrid.loadCppLib(extLibName)

if useTempTopoFile:
    tmpDtmFileName=sgrid.uncompressGeoTiff(dtmFileName,tiled=True)
else:
    tmpDtmFileName=dtmFileName

dryThresh=0.1
channel=False

flowPathOutput=None
extendWlGrid=False

if processEnd:
    wlFileName=os.path.join(resultsDirectory,resultsPrefix+'_wl.csv')
    flowFileName=os.path.join(resultsDirectory,resultsPrefix+'_flow.csv')

    wlGrid=fileIO.readCSV(wlFileName,"WL",xsz,ysz,dataType=sgrid.getPrecision())
    flowX,flowY=fileIO.readFlowCsv(flowFileName,xsz,ysz,dataType=sgrid.getPrecision())

    print("Resampling and saving end depths/flows to file...")

    maskList=numpy.where((wlGrid-storagePar[:,:,0])<dryThresh)
    wlGrid[maskList]=storagePar[:,:,0][maskList]

    sgrid.saveResults(wlGrid,wlGrid,flowX,flowY,storagePar,xsz,ysz,cellSize,xll,yll,
                       flowThreshold,
                       noDataValue,noDataReplacement,
                       resultsDirectory,resultsPrefix+'_final',cppResample3,
                       saveCsv=True)


if processMax:
    wlFileName=os.path.join(resultsDirectory,resultsPrefix+'_max_wl.csv')
    flowFileName=os.path.join(resultsDirectory,resultsPrefix+'_max_flow.csv')

    wlGrid=fileIO.readCSV(wlFileName,"WL",xsz,ysz,dataType=sgrid.getPrecision())
    flowX,flowY=fileIO.readFlowCsv(flowFileName,xsz,ysz,dataType=sgrid.getPrecision())

    print("Resampling and saving max depths/flows to file...")

    maskList=numpy.where((wlGrid-storagePar[:,:,0])<dryThresh)
    wlGrid[maskList]=storagePar[:,:,0][maskList]

    sgrid.saveResults(wlGrid,wlGrid,flowX,flowY,storagePar,xsz,ysz,cellSize,xll,yll,
                       flowThreshold,
                       noDataValue,noDataReplacement,
                       resultsDirectory,resultsPrefix+'_max',cppResample3,
                       saveCsv=True, threads=None)


if useTempTopoFile:
    os.remove(tmpDtmFileName)


