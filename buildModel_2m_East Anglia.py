#!/usr/bin/python

###############################################################################
# CONTROL PANEL

# Topography
dtmFileName=r"/merlin1/Projects/LTIS SLR/GIS/DTM/All_clip_range.tif"
clipPolyName=r"/merlin1/Projects/LTIS SLR/GIS/England_buffer.gpkg:1km" # Provide polygon to clip catchment etc
maskGridName=r"/merlin1/Projects/LTIS SLR/GIS/mask_1km.tif"
highWaterMaskFileName=r"/merlin1/Projects/LTIS SLR/GIS/high_water_grid_2m_100m.tif"
useTempTopoFile=False # Use this to create uncompressed, tiled topo file to speed up access for large grids

# These values can be used to replace NULLs (e.g. at sea) with sensible values
replacement_values = {-10:-10, 20:20}

# Use this to add NULL cells around edge - allows water to fall out of model
addNullEdges=False

# Extent and resolution of model
xll=500000.    # Lower left corner
yll=200000.
cellSize=1000.
xsz=160
ysz=160

threads = 8

# Manning's n
nFloodplain=0.06    # Can omit this if grid data supplied
nFloodplainFile=None # Raster file of roughness, with same size/res as DTM

# Output options
outputPrefix = '2m_EA_'
outputFile= outputPrefix+"params.pck"
gridFileName=outputPrefix+"grid.csv"
saveDtmTiles = True

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
    storageFunc=cppCalcStorageParameters,outputPrefix=outputPrefix, clipRasterPolyName=clipPolyName,
    threads = threads, saveDtmTiles=saveDtmTiles, maskGridName=maskGridName,
    highWaterMaskFileName = highWaterMaskFileName)

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
