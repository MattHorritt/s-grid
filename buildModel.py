#!/usr/bin/python

import cPickle as pickle
import os
import goring_obf as goring

goring.setPrecision32()
extLibName=r"./hydraulics_f32.so"     # Don't change this for demo version

###############################################################################
# CONTROL PANEL
# Output options

dtmFileName=r"GB_SE.tiff"

clipPolyName=None # Provide polygon to clip catchment

# This can be used to replace NULLs at sea with sensible values
noDataValue=None
noDataReplacement=None

outputFile="./params.pck"
gridFileName="grid.csv"


# Extent and resolution of model
xll=395000.
yll=131000.
cellSize=1000.
xsz=215
ysz=102

nFloodplain=0.045 # Uniform value
nFloodplainFile=None # Raster file of roughness, with same size/res as DTM

###############################################################################

# Channel parameters
channel=False

arrayType=goring.getPrecision()

cellSize=arrayType(cellSize)
nFloodplain=arrayType(nFloodplain)
nChannel=nFloodplain

(cppCalcFlow,cppCalcFlowGrid,cppDryCheck,cppTimeStep,cppWlFromVolGrid,\
    cppConveyanceParameters,cppMaxVolGrid,cppResample2,cppResample3, \
    cppFlowPaths,cppSum,cppCalcStorageParameters,cppLazyFlowPaths, \
    cppWlFill,cppBurnFlowPaths,cppMakeWlGrid,cppClipZero,cppDryCheckDiagnostic,
    cppScsAdditionalRunoff,cppCalcFlowEdges)=\
    goring.loadCppLib(extLibName)


print "Parameterising topography..."
tmpDtmFileName=goring.uncompressGeoTiff(dtmFileName,tiled=True)

if clipPolyName is not None:
	goring.clipRasterPoly(tmpDtmFileName,clipPolyName)

convParX, convParY, storagePar=goring.gridFlowSetupTiled(tmpDtmFileName,\
    xll, yll, cellSize, xsz, ysz, nChannel, nFloodplain, \
    nFileName=nFloodplainFile,
    ndv=noDataValue,ndr=noDataReplacement,conveyanceFunc=cppConveyanceParameters,\
    storageFunc=cppCalcStorageParameters,outputPrefix='')

os.remove(tmpDtmFileName)

storagePar[:,0,0]=-9999.
storagePar[:,-1,0]=-9999.
storagePar[0,:,0]=-9999.
storagePar[-1,:,0]=-9999.

if gridFileName is not None:
     goring.writeGridCSV(xll,yll,cellSize,xsz,ysz,storagePar,gridFileName)

file=open(outputFile,"w")
pickle.dump((xll,yll,cellSize,xsz,ysz,convParX, convParY, storagePar),file)
file.close()
