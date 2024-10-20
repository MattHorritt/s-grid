#!/usr/bin/python

###############################################################################
# CONTROL PANEL

# Topography
dtmFileName=r"/merlin1/Projects/LTIS SLR/GIS/DTM/OS Panorama 50m clip range.tif"
clipPolyName="/merlin1/Projects/LTIS SLR/GIS/England_buffer.gpkg:1km" # Provide polygon to clip catchment etc
useTempTopoFile=False # Use this to create uncompressed, tiled topo file to speed up access for large grids

# These values can be used to replace NULLs (e.g. at sea) with sensible values
replacement_values = {-10:20, 20:20}

# Use this to add NULL cells around edge - allows water to fall out of model
addNullEdges=False

# Extent and resolution of model
xll=80000.    # Lower left corner
yll=0.
cellSize=1000.
xsz=600
ysz=700

threads = None

# Manning's n
nFloodplain=0.06    # Can omit this if grid data supplied
nFloodplainFile=None # Raster file of roughness, with same size/res as DTM

# Output options
outputPrefix = '50m_'
outputFile= outputPrefix+"params.pck"
gridFileName=outputPrefix+"grid.csv"

###############################################################################

import pickle
import os
import sgrid
import time


t1=time.time()


# Path to C++ library
extLibName=r"./sgridHydraulics.so"# C++ library

sgrid.setPrecision32()

# Channel parameters
channel=False

arrayType=sgrid.getPrecision()

cellSize=arrayType(cellSize)

# If nFloodplain not defined - use default value
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

convParX, convParY, storagePar=sgrid.gridFlowSetupTiled(tmpDtmFileName,\
    xll, yll, cellSize, xsz, ysz, nChannel, nFloodplain, \
    nFileName=nFloodplainFile,
    rvs = replacement_values,conveyanceFunc=cppConveyanceParameters,\
    storageFunc=cppCalcStorageParameters,outputPrefix=outputPrefix, clipRasterPoly=clipPolyName, threads = threads)

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

t2=time.time()
print("Completed in %0.2fs"%((t2-t1)))
