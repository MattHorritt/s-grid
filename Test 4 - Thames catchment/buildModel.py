#!/usr/bin/python

###############################################################################
# CONTROL PANEL

# Provide paths to Python modules and dll etc
sgridPath = "../"

# Topography
dtmFileName=r"DTM.tiff"
clipPolyName=None # Provide polygon to clip catchment
useTempTopoFile=False # Use this to create uncompressed, tiled topo file to speed up access

# This can be used to replace NULLs at sea with sensible values
noDataValue=None
noDataReplacement=None
addNullEdges=True # Use this to add NULL cells around edge - allows water to fall out of model

# Extent and resolution of model
xll=395000.
yll=143000.
cellSize=1000.
xsz=125
ysz=71

# Manning's n
nFloodplain=0.06
nFloodplainFile=None # Raster file of roughness, with same size/res as DTM

# Output options
outputFile="./params.pck"
gridFileName="grid.csv"

###############################################################################
import sys
sys.path.append(sgridPath) # sgrid.py etc

import pickle
import os
import sgrid
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

arrayType=sgrid.getPrecision()

cellSize=arrayType(cellSize)

nFloodplain=arrayType(nFloodplain)
nChannel=nFloodplain

(cppCalcFlow,cppCalcFlowGrid,cppDryCheck,cppTimeStep,cppWlFromVolGrid,\
    cppConveyanceParameters,cppMaxVolGrid,cppResample2,cppResample3, \
    cppFlowPaths,cppSum,cppCalcStorageParameters,cppLazyFlowPaths, \
    cppWlFill,cppBurnFlowPaths,cppMakeWlGrid,cppClipZero,cppDryCheckDiagnostic,
    cppScsAdditionalRunoff,cppCalcFlowEdges,cppCheckLicence)=\
    sgrid.loadCppLib(extLibName)


print("Parameterising topography...")
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

if addNullEdges:
    storagePar[:,0,0]=-9999.
    storagePar[:,-1,0]=-9999.
    storagePar[0,:,0]=-9999.
    storagePar[-1,:,0]=-9999.

if gridFileName is not None:
     sgrid.writeGridCSV(xll,yll,cellSize,xsz,ysz,storagePar,gridFileName)

file=open(outputFile,"wb")
pickle.dump((xll,yll,cellSize,xsz,ysz,convParX, convParY, storagePar),file)
file.close()
