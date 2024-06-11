#!/usr/bin/python

import cPickle as pickle
import os
import sys

###############################################################################
# CONTROL PANEL
# Provide paths to Python modules and C++ library
sys.path.append('../') # sgrid.py
extLibName=r"../hydraulics_f32.so"# C++ library

# Topography
dtmFileName=r"dtm.tiff"
clipPolyName=None # Provide polygon to clip catchment
useTempTopoFile=False # Use this to create uncompressed, tiled topo file to speed up access

# This can be used to replace NULLs at sea with sensible values
noDataValue=None
noDataReplacement=None

# Extent and resolution of model
xll=-5000.
yll=-5000.
cellSize=1000.
xsz=110
ysz=20

# Manning's n
nFloodplain=0.045 # Uniform value
nFloodplainFile=None # Raster file of roughness, with same size/res as DTM

# Output options
outputFile="./params.pck"
gridFileName="grid.csv"

###############################################################################

import sgrid

sgrid.setPrecision32()


# Channel parameters
channel=False

arrayType=sgrid.getPrecision()

cellSize=arrayType(cellSize)
nFloodplain=arrayType(nFloodplain)
nChannel=nFloodplain

(cppCalcFlow,cppCalcFlowGrid,cppDryCheck,cppTimeStep,cppWlFromVolGrid,\
    cppConveyanceParameters,cppMaxVolGrid,cppResample2,cppResample3, \
    cppFlowPaths,cppSum,cppCalcStorageParameters,cppLazyFlowPaths, \
    cppWlFill,cppBurnFlowPaths,cppMakeWlGrid,cppClipZero,cppDryCheckDiagnostic,
    cppScsAdditionalRunoff,cppCalcFlowEdges)=\
    sgrid.loadCppLib(extLibName)


print "Parameterising topography..."
if useTempTopoFile:
    tmpDtmFileName=sgrid.uncompressGeoTiff(dtmFileName,tiled=True)
else:
    tmpDtmFileName=dtmFileName

if clipPolyName is not None:
	sgrid.clipRasterPoly(tmpDtmFileName,clipPolyName)

convParX, convParY, storagePar=sgrid.gridFlowSetupTiled(tmpDtmFileName,\
    xll, yll, cellSize, xsz, ysz, nChannel, nFloodplain, \
    nFileName=nFloodplainFile,
    ndv=noDataValue,ndr=noDataReplacement,conveyanceFunc=cppConveyanceParameters,\
    storageFunc=cppCalcStorageParameters,outputPrefix='')

if useTempTopoFile:
    os.remove(tmpDtmFileName)

storagePar[:,0,0]=-9999.
storagePar[:,-1,0]=-9999.
storagePar[0,:,0]=-9999.
storagePar[-1,:,0]=-9999.

if gridFileName is not None:
     sgrid.writeGridCSV(xll,yll,cellSize,xsz,ysz,storagePar,gridFileName)

file=open(outputFile,"w")
pickle.dump((xll,yll,cellSize,xsz,ysz,convParX, convParY, storagePar),file)
file.close()
