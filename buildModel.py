#!/usr/bin/python

###############################################################################
# CONTROL PANEL

# Topography
dtmFileName=r"/merlin1/Projects/LTIS SLR/GIS/DTM/OS Panorama 50m.tif"
clipPolyName="/merlin1/Projects/LTIS SLR/GIS/England_2km_buffer.gpkg" # Provide polygon to clip catchment etc
useTempTopoFile=False # Use this to create uncompressed, tiled topo file to speed up access for large grids

# These values can be used to replace NULLs (e.g. at sea) with sensible values
noDataValue=None
noDataReplacement=None

# Use this to add NULL cells around edge - allows water to fall out of model
addNullEdges=False

# Extent and resolution of model
xll=80000.    # Lower left corner
yll=3000.
cellSize=1000.
xsz=600
ysz=700

# Manning's n
nFloodplain=0.06    # Can omit this if grid data supplied
nFloodplainFile=None # Raster file of roughness, with same size/res as DTM

# Output options
outputFile="./params.pck"
gridFileName="grid.csv"

###############################################################################

import pickle
import os
import sgrid

# Path to C++ library
extLibName=r"./sgridHydraulics.so"# C++ library

sgrid.setPrecision32()

# Channel parameters
channel=False

arrayType=sgrid.getPrecision()

cellSize=arrayType(cellSize)

# If nFlodplain not defined - use default value
if 'nFloodplain' not in locals():
    nFloodplain=0.03

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
pickle.dump((xll, yll, cellSize, xsz, ysz, convParX, convParY, storagePar), file)
file.close()
