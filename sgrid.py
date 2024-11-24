import fileIO
import numpy
from bresenhamalgorithm import bresenham
import scipy.stats
import scipy.optimize
import matplotlib.pyplot as plt
import os
import sys
import ctypes
from numpy.ctypeslib import ndpointer
import tempfile
from subprocess import call
import reservoirs
import shapely, shapely.wkb
from osgeo import ogr
from multiprocessing.dummy import Pool as ThreadPool
from itertools import product
import time
from pathlib import Path
import time
from scipy.interpolate import RegularGridInterpolator

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
arrayType=numpy.float64

def setPrecision32():
    global arrayType
    arrayType=numpy.float32

def setPrecision64():
    global arrayType
    arrayType=numpy.float64

def getPrecision():
    return arrayType

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Utility which replaces one value in an array with another
# Done line-by-line to avoid memory hogging

def replaceArrayVals(arr,val1,val2):
    xsz,ysz=arr.shape

    for row in range(ysz):
        rList=numpy.where(arr[:,row]==val1)
        arr[:,row][rList]=val2

    return

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Use ReFH PDM model to calculate runoff (mm) in a time step
Ct=0
previousRainfall=0.
rainfallProfileSummer=[0.024,0.036,0.054,0.087,0.154,0.291,0.154,0.087,0.054,
    0.036,0.024]
rainfallProfileWinter=[0.022,0.037,0.061,0.101,0.165,0.229,0.165,0.101,0.061,
    0.037,0.022]

def refhRunoff(t,stormDuration,rainfallDepth,bfihost,propwet,summerProfile=False):
    global Ct, previousRainfall, rainfallProfile, Cmax

    # First call, so reset initial and max soil moisture etc
    if t<=0.:
        Cmax=596.7*(bfihost**0.95)*(propwet**-0.24)
        Ct=(Cmax/2)*(0.90-0.82*bfihost-0.43*propwet)
        if summerProfile is not None and summerProfile:
            rainfallProfile=rainfallProfileSummer
        else:
            rainfallProfile=rainfallProfileWinter
        previousRainfall=0.

    if t<0 or t>stormDuration:
        return 0

    rainfallStep=int(11.*t/stormDuration)

    if rainfallStep>=1 and t<stormDuration:
        totalRainfall=sum(rainfallProfile[:rainfallStep])
    elif t>=stormDuration:
        totalRainfall=1.0
    else:
        totalRainfall=0.

    if t<stormDuration:
        totalRainfall+=(t-rainfallStep*stormDuration/11.)/(stormDuration/11.)\
            *rainfallProfile[rainfallStep]

    totalRainfall*=rainfallDepth

    totalRainfall=min(totalRainfall,rainfallDepth)

    additionalRainfall=totalRainfall-previousRainfall

    runoff=additionalRainfall*(Ct/Cmax+0.5*additionalRainfall/Cmax)

    Ct+=additionalRainfall
    previousRainfall=totalRainfall


    return runoff


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Use ReFH PDM model to calculate runoff (mm) in a time step
rainfallProfileSummer=[0.024,0.036,0.054,0.087,0.154,0.291,0.154,0.087,0.054,\
    0.036,0.024]
rainfallProfileWinter=[0.022,0.037,0.061,0.101,0.165,0.229,0.165,0.101,0.061,\
    0.037,0.022]

def pRunoff(t,stormDuration,rainfallDepth,pr,summerProfile=False):

    global previousRainfall

    if t<=0:
        previousRainfall=0.

    if summerProfile is not None and summerProfile:
        rainfallProfile=rainfallProfileSummer
    else:
        rainfallProfile=rainfallProfileWinter

    if t<0 or t>stormDuration:
        return 0

    rainfallStep=int(11.*t/stormDuration)

    if rainfallStep>=1 and t<stormDuration:
        totalRainfall=sum(rainfallProfile[:rainfallStep])
    elif t>=stormDuration:
        totalRainfall=1.0
    else:
        totalRainfall=0.

    if t<stormDuration:
        totalRainfall+=(t-rainfallStep*stormDuration/11.)/(stormDuration/11.)\
            *rainfallProfile[rainfallStep]

    totalRainfall*=rainfallDepth

    totalRainfall=min(totalRainfall,rainfallDepth)

    additionalRainfall=totalRainfall-previousRainfall

    runoff=additionalRainfall*pr/100.

    previousRainfall=totalRainfall

    return runoff

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Writes a CSV with WKT containing the calculation grid - useful for display
def writeGridCSV(xll,yll,cellSize,xsz,ysz,sp,fileName): #,rain,runoff,cn,fileName):


    outCsvt=open(fileName+"t","w")
    outCsvt.write("\"Integer\",\"Integer\",\"Integer\",\"String\",\"Real\"\n")
    outCsvt.close()

    csvFile=open(fileName,"w")

    id=0

#    csvFile.write("id;wkt;elevation;rainfall;runoff;cn\n")
    csvFile.write("id;i;j;wkt;elevation\n")

    for i in range(xsz):
        for j in range(ysz):
            x0=xll+i*cellSize
            x1=x0+cellSize
            y0=yll+j*cellSize
            y1=y0+cellSize

            wktString="POLYGON (("
            wktString+="%0.3f %0.3f,"%(x0,y0)
            wktString+="%0.3f %0.3f,"%(x1,y0)
            wktString+="%0.3f %0.3f,"%(x1,y1)
            wktString+="%0.3f %0.3f,"%(x0,y1)
            wktString+="%0.3f %0.3f"%(x0,y0)
            wktString+="))"

            if sp[i,j,0] != -9999:
                csvFile.write("%i;%i;%i;%s;%f\n"%(id,i,j,wktString,sp[i,j,0]))
#               csvFile.write("%i;%s;%f;%f;%f;%f\n"%(id,wktString,sp[i,j,0],rain[i,j],runoff[i,j],cn[i,j]))
                id+=1
    csvFile.close()

    return

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def __extractProfilePolyline(xyList,raster,cellSize,strict=False):
    profile=[]

#    xsz=raster.shape[0]
#    ysz=raster.shape[1]

    numPoints=len(xyList)
    crossSectionLength=0.
    for i in range(numPoints-1):
        x0=xyList[i][0]
        y0=xyList[i][1]
        x1=xyList[i+1][0]
        y1=xyList[i+1][1]

        if x0==x1:
            crossSectionLength+=(y1-y0)
            profile=list(raster[x0,y0:y1])
            profile.append(raster[x1,y1])

        elif y0==y1:
            crossSectionLength+=(x1-x0)
            profile=list(raster[x0:x1,y0])
            profile.append(raster[x1,y1])

        else:
            crossSectionLength+=numpy.sqrt((x1-x0)**2+(y1-y0)**2)

            if x0==x1 and y0==y1:
                profile.append(raster[x0,y0])
            else:
                if i==(numPoints-2):
                    for x,y in bresenham([x0,y0],[x1,y1]).path[:]:
                        profile.append(raster[x,y])
                else:
                    for x,y in bresenham([x0,y0],[x1,y1]).path[:-1]:
                        profile.append(raster[x,y])

    # What's the average distance between profile points?
    dl=cellSize*crossSectionLength/len(profile)

    # If we're being strict, don't include last cell - avoids double counting
    if strict:
        profile=profile[0:-1]

    return dl, profile

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def __qf(x,a,b):
    return a*x*x+b*x

def quadFit(xList,yList):
    res=scipy.optimize.curve_fit(__qf, numpy.array(xList), numpy.array(yList))

    return res[0]


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def calcStorageParametersChannel(rect,z,xll,yll,dx,cellSize,\
    cag,cagXll,cagYll,cagDx,chanExp,chanMult,chanAR,chanMaxD):

    xi0=int((rect[0][0]-xll)/dx)
    xi1=int((rect[1][0]-xll)/dx)

    yi0=int((rect[0][1]-yll)/dx)
    yi1=int((rect[1][1]-yll)/dx)

    cell=z[xi0:(xi1+1),yi0:(yi1+1)]

    cell=cell.reshape((-1))


    cell.sort()


# Get the area grid tile
    xi0=int((rect[0][0]-cagXll)/cagDx)
    xi1=int((rect[1][0]-cagXll)/cagDx)

    yi0=int((rect[0][1]-cagYll)/cagDx)
    yi1=int((rect[1][1]-cagYll)/cagDx)

#    cagCell=cag[xi0:(xi1+1),yi0:(yi1+1)]
    cagXsz,cagYsz=cag.shape

    if xi0>0 and xi0<cagXsz and \
        xi1>0 and xi1<cagXsz and \
        yi0>0 and yi0<cagYsz and \
        yi1>0 and yi1<cagYsz:

        cagCell=cag[xi0:xi1,yi0:yi1]

        drainageArea=cagCell.max()
    else:
        drainageArea=0.


    width=chanMult*(drainageArea**chanExp)
    depth=min(width/chanAR,chanMaxD)

    zBank=cell.min()
    chanVol=width*depth*cellSize # Fix this with cellSize
    chanArea=width*cellSize
    numChanCells=int(chanArea/(dx*dx))

    #Remove lowest cells and replace with bed
    if numChanCells>0:
        cell=cell[numChanCells:]
        cell=list(cell)+[zBank-depth]*numChanCells

    zMin=min(cell)
    zMax=max(cell)

    trueVol=[] # Actually volume per unit area


    if zMax>(zBank+5):
        wlIP=[zMin,zBank,zBank+1,zBank+5,zMax]
    else:
        wlIP=[zMin,zBank,zBank+1,zBank+5,zBank+5]
        zMax=zBank+5

    volIP=[]

    for wl in wlIP:
        hi=sum([max(0,wl-zi) for zi in cell])
        volIP.append(hi/len(cell))

    return [zMin, zBank, zMax]+volIP[1:]



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def calcStorageParameters(rect,z,xll,yll,dx,plotName=None,csvOutput=None):

    xi0=int((rect[0][0]-xll)/dx)
    xi1=int((rect[1][0]-xll)/dx)

    yi0=int((rect[0][1]-yll)/dx)
    yi1=int((rect[1][1]-yll)/dx)

    cell=z[xi0:(xi1+1),yi0:(yi1+1)]

    cell=cell.reshape((-1))

    zMin=cell.min()

    zMax=cell.max()

    if numpy.isnan(zMin) or numpy.isnan(zMax):
        print(cell)

    trueVol=[] # Actually volume per unit area


    if zMax>(zMin+5) and zMax>(zMin+1):
        wlIP=[zMin,zMin+1,zMin+5,zMax]
    else:
        wlIP=[zMin,zMin+1,zMin+5,zMin+5]
        zMax=zMin+5

    volIP=[]

    for wl in wlIP:
        hi=sum([max(0,wl-zi) for zi in cell])
        volIP.append(hi/len(cell))



    return zMin, zMax, volIP[1], volIP[2], volIP[3]




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def conveyanceParameters(profile,dx,n,plotName=None,csvOutput=None, nList=None):
    maxDepth=10.

    minZ=min(profile)
    maxZ=max(profile)

    nConveyanceLevels=5
    wlList=[]
    for i in range(nConveyanceLevels):
        wlf=(10**(-2.+i*2./(nConveyanceLevels)))
        wlList.append(minZ+wlf*maxDepth)

    conveyanceList=[]
    conveyanceUsingMinList=[]

    if nList is None:
        nList=[n]*len(wlList)

    for wl in wlList:
        Ki=0.

        conveyanceUsingMinList.append(dx*len(profile)*((wl-minZ)**1.666667)/n)

#        for z in profile:
#            Ki+=dx*(max(0,wl-z)**1.666667)/n

#        Ki=sum([dx*(max(0,wl-z)**1.6667)/n for z in profile])

        for i,z in enumerate(profile):
            Ki+=dx*(max(0,wl-z)**1.6667)/nList[i]

        conveyanceList.append(Ki)

    # Transform to logs
    conveyanceUsingMinList=[numpy.log(k) for k in conveyanceUsingMinList]
    conveyanceList=[numpy.log(k) for k in conveyanceList]

    res=scipy.stats.linregress(conveyanceUsingMinList,conveyanceList)
    slope=res[0]
    intercept=res[1]


    if plotName is not None: # Produce some graphical output to check fit
        plt.figure(1,facecolor='w',edgecolor='w')

        conveyanceList=[numpy.exp(k) for k in conveyanceList]
        try:
            plt.semilogx(conveyanceList,wlList,'+',color='k')
        except:
            print(conveyanceList,wlList)
            assert False

        conveyanceApprox=[numpy.exp(slope*km+intercept) for km in conveyanceUsingMinList]

        plt.semilogx(conveyanceApprox,wlList,'-',color='k')

        if csvOutput:
              f=open(plotName+'.csv',"w")
              f.write("WL,K,Kapprox\n")
              for i, wl in enumerate(wlList):
                  f.write("%f,%f,%f\n"%(wl,conveyanceList[i],conveyanceApprox[i]))
              f.close()

        rmsHeightError=0.
        c=0
        for i, Ki in enumerate(conveyanceList):
            wl=wlList[i]

#            if wl-minZ>5:
#                continue

            kmin=numpy.exp((numpy.log(Ki)-intercept)/slope)

            wlApprox=minZ+(kmin*n/(dx*len(profile)))**0.6

            rmsHeightError+=(wlApprox-wl)**2

            c+=1

        rmsHeightError=numpy.sqrt(rmsHeightError/c)

        plt.title(plotName+" RMS Error=%fm"%rmsHeightError)
        plt.xlabel('Conveyance (m3s-1)')
        plt.ylabel('WL (m)')

        plt.savefig(plotName+'.png', bbox_inches=0)
        plt.close('all')


    return minZ, slope, intercept

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Test if any of one geometry in a vector layer intersects with a shapely polygon
def vectors_intersect(l1, p2):
    l1.ResetReading()
    f1 = l1.GetNextFeature()

    while f1 is not None:
        g1 = f1.GetGeometryRef()
        p1 = shapely.wkb.loads(bytes(g1.ExportToIsoWkb()))

        e1 = shapely.envelope(p1).bounds
        x11 = e1[0]
        x12 = e1[2]
        y11 = e1[1]
        y12 = e1[3]

        e2 = shapely.envelope(p2).bounds
        x21 = e2[0]
        x22 = e2[2]
        y21 = e2[1]
        y22 = e2[3]

        if x12 >= x21 and x11 <= x22 and y12 >= y21 and y11 <= y22:
            if p1.intersects(p2):
                return True

        f1 = l1.GetNextFeature()

    return False

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Wrapper for processCell so we can pass single argument tuple for threading.
# Bit of a hack - but works.
def processCellWrapper(argTuple):
    # print(f"In processCellWrapper processing cell [{argTuple[0]},{argTuple[1]}]...")

    dtm = fileIO.geoGrid(argTuple[5], objOnly = True)

    retTuple = processCell(argTuple[0], argTuple[1], argTuple[2],
                argTuple[3], argTuple[4], dtm,
                argTuple[6], argTuple[7], argTuple[8],
                nFileObj = argTuple[9], clipLyr = argTuple[10], rvs = argTuple[11],
                saveDtmTiles = argTuple[12], maskGrid = argTuple[13], highWaterMaskFileName = argTuple[14])

    return argTuple, retTuple # Return arguments as well - useful for tracking the results

def processCell(i, j, xll, yll, cellSize, dtm, nFp, conveyanceFunc,storageFunc, nFileObj = None, clipLyr = None,
                rvs = None, saveDtmTiles = None, maskGrid = None, highWaterMaskFileName = None):
    # These are the square extents in map coordinates
    x0 = xll + i * cellSize
    x1 = x0 + cellSize
    y0 = yll + j * cellSize
    y1 = y0 + cellSize

    # Return arrays for results
    conveyanceValuesX = numpy.zeros(7,dtype=arrayType) - 9999
    conveyanceValuesY = numpy.zeros(7,dtype=arrayType) - 9999
    storageValues = numpy.zeros(5,dtype=arrayType) - 9999

    # These are the square extents in array coordinates. All these are INSIDE the box.
    xi0 = int((x0 - dtm.xll) / dtm.dx)
    xi1 = int((x1 - dtm.xll) / dtm.dx) - 1
    yi0 = int((y0 - dtm.yll) / dtm.dx)
    yi1 = int((y1 - dtm.yll) / dtm.dx) - 1

    windowXsz = xi1 - xi0 + 1
    windowYsz = yi1 - yi0 + 1

    if xi0 >= 0 and xi1 <= dtm.xsz - 1 and yi0 >= 0 and yi1 <= dtm.ysz - 1:  # Within DTM extent - grab window

        if maskGrid is not None:
            xc = 0.5 * (x0 + x1)
            yc = 0.5 * (y0 + y1)
            if maskGrid.getPointValueXy(xc, yc) == 0:
                return conveyanceValuesX, conveyanceValuesY, storageValues

        dtmWindow = dtm.obj.ReadAsArray(xoff=xi0, yoff=dtm.ysz - 1 - yi1, xsize=windowXsz,
                                           ysize=windowYsz).transpose().copy()
        dtmWindow = numpy.array(dtmWindow[:, ::-1], arrayType)

        # Use -9999 as no data value as this is what's used in C++ code
        dtmWindow[numpy.where(dtmWindow == dtm.noDataValue)] = -9999

        # Replace values according to replacement values dict
        if rvs is not None:
            for v1, v2 in rvs.items():
                if v1 != v2:
                    dtmWindow[numpy.where(dtmWindow == v1)] = v2

        # TODO Might need to remove this
        # if dtmWindow.max() < 1.0:
        #     return conveyanceValuesX, conveyanceValuesY, storageValues

        # Set values below high water to -10
        if highWaterMaskFileName is not None:
            highWaterMask = fileIO.geoGrid(highWaterMaskFileName, objOnly=True)
            highWaterMaskWindow = highWaterMask.obj.ReadAsArray(xoff=xi0, yoff=dtm.ysz - 1 - yi1, xsize=windowXsz,
                                               ysize=windowYsz).transpose().copy()
            highWaterMaskWindow = highWaterMaskWindow[:, ::-1]
            dtmWindow[numpy.where(highWaterMaskWindow == 1)] = 20



        # Decide whether to include this cell - if all NaNs, or all replacement values, skip
        if numpy.all(dtmWindow == -9999):
            return conveyanceValuesX, conveyanceValuesY, storageValues
        if rvs is not None:
            for v1 in rvs.keys():
                if numpy.all(dtmWindow == v1):
                    return conveyanceValuesX, conveyanceValuesY, storageValues

        # Check if this square intersect clip polygon if given - this test is potentially slow - so do last
        if clipLyr is not None:
            if isinstance(clipLyr, str):
                if ':' in clipLyr:
                    fileName = clipLyr.split(':')[0]
                    lyrName = clipLyr.split(':')[1]
                    clipRasterPolyDataSource = ogr.Open(fileName)
                    clipRasterPolyLayer = clipRasterPolyDataSource.GetLayer(lyrName)
                else:
                    clipRasterPolyDataSource = ogr.Open(clipLyr)
                    clipRasterPolyLayer = clipRasterPolyDataSource.GetLayerByIndex(0)

                sq = shapely.Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])

                if not vectors_intersect(clipRasterPolyLayer, sq):
                    return conveyanceValuesX, conveyanceValuesY, storageValues

            else:
                sq = shapely.Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])

                if not vectors_intersect(clipLyr, sq):
                    return conveyanceValuesX, conveyanceValuesY, storageValues

        # Save tile
        if saveDtmTiles:
            if '.' in dtm.fileName:
                tileRoot = dtm.fileName.split('.')[0]
            else:
                tileRoot = dtm.fileName

            tileRoot += '_tiles'

            # If folder doesn't exist, create it
            Path(tileRoot).mkdir(parents=True, exist_ok=True)

            fileName = Path(tileRoot) / f"{i:03d}_{j:03d}.tif"
            tile = fileIO.geoGrid(dtmWindow, windowXsz, windowYsz, x0, y0, dtm.dx)
            tile.save(str(fileName))

        # And get landcover window
        if nFileObj is not None:
            nWindow = nFileObj.ReadAsArray(xoff=xi0, yoff=dtm.ysz - yi1, xsize=windowXsz,
                                           ysize=windowYsz).transpose().copy()
            nWindow = nWindow[:, ::-1]
        else:
            nWindow = dtmWindow.copy()
            nWindow[:, :] = nFp

        # X-direction
        conveyanceFunc(0, 0, 0, yi1 - yi0, \
                       dtmWindow, dtm.dx, windowXsz, windowYsz, nFp, \
                       conveyanceValuesX, False, \
                       0., dtmWindow, 0, 0, 0, 0, \
                       0, 0, 0, 0, 0, 0, nWindow)

        # Y-direction
        conveyanceFunc(0, 0, xi1 - xi0, 0, \
                       dtmWindow, dtm.dx, xi1 - xi0 + 1, yi1 - yi0 + 1, nFp, \
                       conveyanceValuesY, False, \
                       0., dtmWindow, 0, 0, 0, 0, \
                       0, 0, 0, 0, 0, 0, nWindow)

        # Cell storage
        if storageFunc is not None:
            storageFunc(arrayType(x0), arrayType(y0), arrayType(x1), arrayType(y1), dtmWindow, \
                        windowXsz, windowYsz, arrayType(x0), arrayType(y0), \
                        arrayType(dtm.dx), storageValues, \
                        False, dtmWindow, 0, 0, \
                        0, 0, 0, \
                        0, 0, 0, 0, 0)

    return conveyanceValuesX, conveyanceValuesY, storageValues

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def gridFlowSetupTiled(dtmFileName,xll,yll,cellSize,xsz,ysz,nChan,nFP,
    nFileName=None,
    plotNamePrefix=None, outputPrefix=None,
    rvs=None,ndr=None,conveyanceFunc=None,storageFunc=None,
    clipRasterPolyName = None, threads = None, saveDtmTiles = None, maskGridName = None,
    highWaterMaskFileName = None):

    if plotNamePrefix is None:
        plotNamePrefix=""

    if outputPrefix is None:
        outputPrefix=""



    # dtmFileObj,dtmCellSize,dtmXsz,dtmYsz,dtmXll,dtmYll=fileIO.readScalarGridObj(dtmFileName)
    dtm = fileIO.geoGrid(dtmFileName, objOnly = True)
    # fileNdv = dtmFileObj.GetRasterBand(1).GetNoDataValue()

    if maskGridName is not None:
        maskGrid = fileIO.geoGrid(maskGridName)
    else:
        maskGrid = None

    if nFileName is not None:
        nFileObj,lcCellSize,lcXsz,lcYsz,lcXll,lcYll=fileIO.readScalarGridObj(nFileName)

    # # Might need to fettle no data value
    # ndv=arrayType(ndv)

    convParX=numpy.zeros((xsz+1,ysz,7),dtype=arrayType)-9999.
    convParY=numpy.zeros((xsz,ysz+1,7),dtype=arrayType)-9999.

    storagePar=numpy.zeros((xsz,ysz,5),dtype=arrayType)-9999.

    plotName=None

    ticker=0

    tickerStep=int(xsz * ysz / 100)

    t1 = time.time()

    if threads is None:

        if clipRasterPolyName is not None:
            if ':' in clipRasterPolyName:
                fileName = clipRasterPolyName.split(':')[0]
                lyrName = clipRasterPolyName.split(':')[1]
                clipRasterPolyDataSource = ogr.Open(fileName)
                clipRasterPolyLayer = clipRasterPolyDataSource.GetLayer(lyrName)
            else:
                clipRasterPolyDataSource = ogr.Open(clipRasterPolyName)
                clipRasterPolyLayer = clipRasterPolyDataSource.GetLayerByIndex(0)

            if clipRasterPolyDataSource is None:
                raise ValueError("Can't open clip polygon")

            clipRasterPolyLayer.ResetReading()
        else:
            clipRasterPolyLayer = None

        for i in range(xsz):
            for j in range(ysz):

                if (ticker%tickerStep)==0:
                    pc = int(100.*ticker/(xsz * ysz))
                    print(f"{pc}%% {i}/{xsz},{j}/{ysz} "%(), end='')

                    if ticker > 0:
                        pcComplete = float(ticker) / (xsz * ysz)
                        pcToGo = 1. - pcComplete
                        t2 = time.time() - t1
                        projectedFinish = time.time() + pcToGo * (t2 / pcComplete)

                        if pcToGo * (t2 / pcComplete) < 86400: # <1 day, report time only
                            projectedFinishString = time.strftime("%H:%M:%S", time.localtime(projectedFinish))
                        else: # Report date too
                            projectedFinishString = time.strftime("%d/%m/%y %H:%M:%S", time.localtime(projectedFinish))

                        print(projectedFinishString, end='')

                    print("...", end='\r', flush=True)
                    sys.stdout.flush()

                ticker += 1


                cX, cY, st = processCell(i, j, xll, yll, cellSize, dtm, nFP, conveyanceFunc, storageFunc,
                                         clipLyr = clipRasterPolyLayer, rvs = rvs, saveDtmTiles = saveDtmTiles,
                                         maskGrid = maskGrid, highWaterMaskFileName = highWaterMaskFileName)

                convParX[i, j, :] = cX
                convParY[i, j, :] = cY
                storagePar[i, j, :] = st
    else:
        pool = ThreadPool(threads)
        funcArgList = []
        returnValues = []
        for i, j in product(range(xsz), range(ysz)):
            funcArgList.append((i, j, xll, yll, cellSize, dtmFileName, nFP, conveyanceFunc, storageFunc,
                                         None, clipRasterPolyName, rvs, saveDtmTiles, maskGrid, highWaterMaskFileName))

        counter = 0

        tickerStep=int(xsz * ysz / 100)
        ticker = 0

        for arg, ret in pool.imap(processCellWrapper, funcArgList, chunksize=10):

            if (ticker%tickerStep)==0:
                print("%i%% "%(100.*ticker/(xsz * ysz)), end='')

                if ticker > 0:
                    pcComplete = float(ticker) / (xsz * ysz)
                    pcToGo = 1. - pcComplete
                    t2 = time.time() - t1
                    projectedFinish = time.time() + pcToGo * (t2 / pcComplete)
                    projectedFinishString = time.strftime("%H:%M:%S", time.localtime(projectedFinish))
                    print(projectedFinishString, end='')

                print("...", end='')
                sys.stdout.flush()

            ticker += 1

            i = arg[0]
            j = arg[1]
            convParX[i, j, :] = ret[0]
            convParY[i, j, :] = ret[1]
            storagePar[i, j, :] = ret[2]

        pass


    print("Done.")

    fileIO.saveConveyanceParametersCSV(convParX,convParY,xll,yll,cellSize,\
        outputPrefix+"conveyanceParams.csv")
    fileIO.saveStorageParametersCSV(storagePar,xll,yll,cellSize,\
        outputPrefix+"storageParams.csv")

    return convParX, convParY, storagePar #, bankLevel


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def gridFlowSetup(dtmFileName,xll,yll,cellSize,xsz,ysz,nChan,nFP,
    plotNamePrefix=None, outputPrefix=None,
    ndv=None,ndr=None,conveyanceFunc=None,storageFunc=None, \
    catchmentAreaGridFile=None,chanExp=None,chanMult=None,chanAR=None,chanMaxD=None):

    if plotNamePrefix is None:
        plotNamePrefix=""

    if outputPrefix is None:
        outputPrefix=""

    dtm,dtmCellSize,dtmXsz,dtmYsz,dtmXll,dtmYll=fileIO.readScalarGrid(dtmFileName,\
        dataType=arrayType)

    if catchmentAreaGridFile is not None:
        catchmentAreaGrid,cagCellSize,cagXsz,cagYsz,cagXll,cagYll=fileIO.readScalarGrid(catchmentAreaGridFile,\
            dataType=arrayType)

        catchmentAreaGrid[numpy.where(numpy.isnan(catchmentAreaGrid))]=0.


    # Might need to fettle no data value
    ndv=arrayType(ndv)

    if ndv is not None and ndr is not None:
        replaceArrayVals(dtm,ndv,ndr)
#        dtm[numpy.where(dtm==ndv)]=ndr

    dtm[numpy.where(numpy.isnan(dtm))]=ndr

    convParX=numpy.zeros((xsz+1,ysz,6),dtype=arrayType)
    convParY=numpy.zeros((xsz,ysz+1,6),dtype=arrayType)

    if catchmentAreaGridFile is None:
        storagePar=numpy.zeros((xsz,ysz,5),dtype=arrayType)
    else:
        storagePar=numpy.zeros((xsz,ysz,7),dtype=arrayType)

    bankLevel=numpy.zeros((xsz,ysz),dtype=arrayType)

    count=0
    plotName=None

    returnArray=numpy.zeros(6,dtype=arrayType)

    # Parameters for flow in x-direction
    ticker=0
    print("Calculating X-conveyance",end='')
    for i in range(xsz+1):
        if (ticker%int(xsz/10))==0: print(".", end='')
        ticker+=1
        for j in range(ysz):
            if i==0 or i==xsz:
                convParX[i,j,:]=-9999. # Signalled as blockage
            else:
                x=xll+i*cellSize
                y0=yll+j*cellSize
                y1=y0+cellSize

                xi=int((x-dtmXll)/dtmCellSize)
                yi0=int((y0-dtmYll)/dtmCellSize)
                yi1=int((y1-dtmYll)/dtmCellSize)

                if xi<0 or xi>=dtmXsz or yi0<0 or yi1>dtmYsz:
                    convParX[i,j,:]=-9999.
                    continue

                if conveyanceFunc is None:
                    dl,profile=__extractProfilePolyline([(xi,yi0),(xi,yi1)],dtm,dtmXll,dtmYll,dtmCellSize)

                    if catchmentAreaGridFile is not None:
                        cag_xi=int((x-cagXll)/cagCellSize)
                        cag_yi0=int((y0-cagYll)/cagCellSize)
                        cag_yi1=int((y1-cagYll)/cagCellSize)

                        if cag_xi>0 and cag_xi<cagXsz and \
                            cag_yi0>0 and cag_yi1<cagYsz and \
                            cag_yi1>0 and cag_yi1<cagYsz:

                            cagdl,cagProfile=__extractProfilePolyline(\
                                [(cag_xi,cag_yi0),(cag_xi,cag_yi1)],catchmentAreaGrid,cagXll,cagYll,cagCellSize,strict=True)

                            width=chanMult*(max(cagProfile)**chanExp)
                        else:
                            width=0.

                        depth=min(width/chanAR,chanMaxD)

                        oldProfile=profile[:]
                        nList=__lowerChannelCells(profile,dl,width,depth,nChan,nFP)


                    else:
                        nList=None

                    convParX[i,j,:]=conveyanceParameters(profile,dl,nFP,\
                        plotName=plotName,csvOutput=False,nList=nList)

                else:
                    if catchmentAreaGridFile is not None:

                        cag_xi=int((x-cagXll)/cagCellSize)
                        cag_yi0=int((y0-cagYll)/cagCellSize)
                        cag_yi1=int((y1-cagYll)/cagCellSize)

                        conveyanceFunc(xi,yi0,xi,yi1,\
                            dtm,dtmCellSize,dtmXsz,dtmYsz,nFP,\
                            returnArray,
                            True,nChan,catchmentAreaGrid,cagXsz,cagYsz,
                            cag_xi,cag_yi0,cag_xi,cag_yi1,
                            chanMult,chanExp,chanAR,chanMaxD)

                    else:

                        conveyanceFunc(xi,yi0,xi,yi1,\
                            dtm,dtmCellSize,dtmXsz,dtmYsz,nFP,\
                            returnArray,False,\
                            0.,dtm,0,0,0,0, \
                            0,0,0,0,0,0)

                    convParX[i,j,:]=returnArray[:]

    print(" Done")

#    assert False

    # Parameters for flow in y-direction
    count=0
    ticker=0
    print("Calculating Y-conveyance",end='')
    for i in range(xsz):
        if (ticker%int(xsz/10))==0: print(".",end='')
        ticker+=1
        for j in range(ysz+1):
            if j==0 or j==ysz:
                convParY[i,j,:]=-9999. # Signalled as blockage
            else:
                x0=xll+i*cellSize
                x1=x0+cellSize
                y=yll+j*cellSize

                xi0=int((x0-dtmXll)/dtmCellSize)
                xi1=int((x1-dtmXll)/dtmCellSize)
                yi=int((y-dtmYll)/dtmCellSize)

                if yi<0 or yi>=dtmYsz or xi0<0 or xi1>dtmXsz:
                    convParY[i,j,:]=-9999.
                    continue

                if conveyanceFunc is None:
                    dl,profile=__extractProfilePolyline([(xi0,yi),(xi1,yi)],dtm,dtmXll,dtmYll,dtmCellSize)

                    if catchmentAreaGridFile is not None:
                        cag_xi0=int((x0-cagXll)/cagCellSize)
                        cag_xi1=int((x1-cagXll)/cagCellSize)
                        cag_yi=int((y-cagYll)/cagCellSize)

                        if cag_xi0>0 and cag_xi0<cagXsz and \
                            cag_xi1>0 and cag_xi1<cagXsz and \
                            cag_yi>0 and cag_yi<cagYsz:

                            cagdl,cagProfile=__extractProfilePolyline(\
                                [(cag_xi0,cag_yi),(cag_xi1,cag_yi)],catchmentAreaGrid,cagXll,cagYll,cagCellSize,strict=True)

                            width=chanMult*(max(cagProfile)**chanExp)
                        else:
                            width=0.

                        depth=min(width/chanAR,chanMaxD)
                        nList=__lowerChannelCells(profile,dl,width,depth,nChan,nFP)
                    else:
                        nList=None

#                    if count%100==0:
#                        plotName=plotNamePrefix+'conveyancePlotY%06i'%count
#                    else:
#                        plotName=None
#                    count+=1

                    convParY[i,j,:]=conveyanceParameters(profile,dl,nFP,\
                        plotName=plotName,csvOutput=False,nList=nList)

                else:

                    if catchmentAreaGridFile is not None:

                        cag_xi0=int((x0-cagXll)/cagCellSize)
                        cag_xi1=int((x1-cagXll)/cagCellSize)
                        cag_yi=int((y-cagYll)/cagCellSize)

                        conveyanceFunc(xi0,yi,xi1,yi,\
                            dtm,arrayType(dtmCellSize),dtmXsz,dtmYsz,\
                            nFP,\
                            returnArray,\
                            True,nChan,catchmentAreaGrid,cagXsz,cagYsz,\
                            cag_xi0,cag_yi,cag_xi1,cag_yi,\
                            chanMult,chanExp,chanAR,chanMaxD)


                    else:

                        conveyanceFunc(xi0,yi,xi1,yi,\
                            dtm,arrayType(dtmCellSize),dtmXsz,dtmYsz,\
                            nChan,\
                            returnArray,False,\
                            0,dtm,0,0,0,0,\
                            0,0,0,0,0,0)

                    convParY[i,j,:]=returnArray[:]

    print(" Done")

    # Cell storage curve
    count=0
    ticker=0
    print("Calculating storage",end='')

    dtmXtr=dtmXll+dtmXsz*dtmCellSize;
    dtmYtr=dtmYll+dtmYsz*dtmCellSize;
    tmp5=numpy.zeros(5,dtype=arrayType)
    tmp7=numpy.zeros(7,dtype=arrayType)

    for i in range(xsz):
        if (ticker%int(xsz/10))==0: print(".",end='')
        ticker+=1
        for j in range(ysz):
            x0=xll+i*cellSize
            x1=x0+cellSize

            y0=yll+j*cellSize
            y1=y0+cellSize

            if x0<dtmXll or x1>dtmXtr or y0<dtmYll or y1>dtmYtr:
                storagePar[i,j,:]=-9999.
                continue

            if storageFunc is not None:


                if catchmentAreaGridFile is not None:
                    storageFunc(arrayType(x0),arrayType(y0),arrayType(x1),arrayType(y1),dtm,\
                        dtmXsz,dtmYsz,arrayType(dtmXll),arrayType(dtmYll), \
                        arrayType(dtmCellSize),tmp7,\
                        True,catchmentAreaGrid,cagXsz,cagYsz,\
                        cagXll,cagYll,cagCellSize,\
                        chanMult,chanExp,chanAR,chanMaxD,cellSize)
                    storagePar[i,j,:]=tmp7
                else:
                    storageFunc(arrayType(x0),arrayType(y0),arrayType(x1),arrayType(y1),dtm,\
                        dtmXsz,dtmYsz,arrayType(dtmXll),arrayType(dtmYll), \
                        arrayType(dtmCellSize),tmp5,\
                        False,dtm,0,0,\
                        0,0,0,\
                        0,0,0,0,0)
                    storagePar[i,j,:]=tmp5
            else:
                if catchmentAreaGridFile is not None:
                    storagePar[i,j,:]=calcStorageParametersChannel([(x0,y0),(x1,y1)],\
                        dtm,dtmXll,dtmYll,dtmCellSize,\
                        catchmentAreaGrid,cagXll,cagYll,cagCellSize,\
                        chanExp,chanMult,chanAR,chanMaxD)
                else:
                    storagePar[i,j,:]=calcStorageParameters([(x0,y0),(x1,y1)],\
                        dtm,dtmXll,dtmYll,dtmCellSize,plotName=plotName,csvOutput=False)

    print("Done.")

    fileIO.saveConveyanceParametersCSV(convParX,convParY,xll,yll,cellSize,\
        outputPrefix+"conveyanceParams.csv")
    fileIO.saveStorageParametersCSV(storagePar,xll,yll,cellSize,\
        outputPrefix+"storageParams.csv")

    return convParX, convParY, storagePar #, bankLevel

def saveMonitoringPoints(f,pts,t,wlGrid,sp,qX=None,qY=None):
    f.write("%f"%t)
    for p in pts:
        i=p[0]
        j=p[1]
        f.write(",%f,%f"%(wlGrid[i,j],wlGrid[i,j]-sp[i,j,0]))

        if qX is not None:
            f.write(",%f,%f,%f,%f"%(qX[i,j],qX[i+1,j],qY[i,j],qY[i,j+1]))


    f.write("\n")
    return


def dryCheck(v,Qx,Qy,sp,dt,verbose=False,sources=None):

    maxIt=10

    for i in range(maxIt):
        nc=__dc(v,Qx,Qy,sp,dt,verbose,sources)
        if verbose:
            print("In dryCheck, number of flows changed=", nc)
        if nc==0: break

    return i+1


def __dc(v,Qx,Qy,sp,dt,verbose=False,sources=None):
    xsz,ysz=v.shape

    nChanged=0

    if sources is not None: # Don't forget to remove these later
        for s in sources:
            v[s[0],s[1]]+=dt*s[2]

    for i in range(xsz):
        for j in range(ysz):
            Q1t=-Qx[i,j]*dt # Flows out
            Q2t=Qx[i+1,j]*dt
            Q3t=-Qy[i,j]*dt
            Q4t=Qy[i,j+1]*dt

            if (Q1t+Q2t+Q3t+Q4t)>v[i,j]:
                alpha=v[i,j]/(Q1t+Q2t+Q3t+Q4t)
#                qPosSum=Q1t*(Q1t>0)+Q2t*(Q2t>0)+Q3t*(Q3t>0)+Q4t*(Q4t>0)
#                qNegSum=Q1t*(Q1t<0)+Q2t*(Q2t<0)+Q3t*(Q3t<0)+Q4t*(Q4t<0)
#
#                alpha=v[i,j]*qPosSum/qNegSum

#                if Q1t<0: Qx[i,j]*=alpha
#                if Q2t<0: Qx[i+1,j]*=alpha
#                if Q3t<0: Qy[i,j]*=alpha
#                if Q4t<0: Qy[i,j+1]*=alpha

                if Q1t>0:
                    Qx[i,j]*=alpha
                    nChanged+=1
                if Q2t>0:
                    Qx[i+1,j]*=alpha
                    nChanged+=1
                if Q3t>0:
                    Qy[i,j]*=alpha
                    nChanged+=1
                if Q4t>0:
                    Qy[i,j+1]*=alpha
                    nChanged+=1

                if verbose:

                    Q1=Qx[i,j]
                    Q2=-Qx[i+1,j]
                    Q3=Qy[i,j]
                    Q4=-Qy[i,j+1]

                    vNew=v[i,j]+dt*(Q1+Q2+Q3+Q4)

                    if vNew<-1:
                        print("vNew=",vNew)
                        print(i,j, alpha, v[i,j])
                        print(Q1t, Q2t, Q3t, Q4t)
                        print(Q1*dt,Q2*dt,Q3*dt,Q4*dt)


    if sources is not None: # Removing these now
        for s in sources:
            v[s[0],s[1]]-=dt*s[2]

    return nChanged

def timeStep(v,Qx,Qy,dt,sources=None):
    xsz,ysz=v.shape

    if sources is not None:
        for s in sources:
            v[s[0],s[1]]+=dt*s[2]

    for i in range(xsz):
        for j in range(ysz):
            Q1=Qx[i,j]
            Q2=-Qx[i+1,j]
            Q3=Qy[i,j]
            Q4=-Qy[i,j+1]

            v[i,j]+=dt*(Q1+Q2+Q3+Q4)

def trappedVol(sp,cpx,cpy):
    xsz,ysz,dum=sp.shape

    tVol=0.

    for i in range(xsz):
        for j in range(ysz):
            zMin=sp[i,j,0]
            zx1=cpx[i,j,0]
            zx2=cpx[i+1,j,0]
            zy1=cpy[i,j,0]
            zy2=cpy[i,j+1,0]

            if zMin<min(zx1,zx2,zy1,zy2):
                tVol+=volFromWl(min(zx1,zx2,zy1,zy2),i,j,sp)

    return tVol


def volFromWl(wl,i,j,sp,dx):

    if sp.shape[2]==5: # No channel
        zMin=sp[i,j,0]
        zMax=sp[i,j,1]

        vips=sp[i,j,2:]

        if wl>=zMax:
            v=dx*dx*(vips[-1]+(wl-zMax))
        else:
            v=dx*dx*numpy.interp(wl,[zMin,zMin+1,zMin+5,zMax],numpy.concatenate(([0],vips)))

    else: # With channel
        zMin=sp[i,j,0]
        zBank=sp[i,j,1]
        zMax=sp[i,j,2]

        vips=sp[i,j,3:]

        if wl>=zMax:
            v=vips[-1]+dx*dx*(wl-zMax)
        else:
            v=dx*dx*numpy.interp(wl,[zMin,zBank,zBank+1,zBank+5,zMax],numpy.concatenate(([0],vips)))

    return v

def wlFromVol(v,i,j,sp,dx):

    if len(sp)==5: # No channel

        zMin=sp[i,j,0]
        zMax=sp[i,j,1]
        vip1=sp[i,j,2]
        vip2=sp[i,j,3]
        vip3=sp[i,j,4]

        if v/(dx*dx)>=vip3: # Cell size
            wl=(v/(dx*dx)-vip3)+zMax
        else:
            wl=numpy.interp(v/(dx*dx),[0,vip1,vip2,vip3],[zMin,zMin+1,zMin+5,zMax])
    else:
        zMin=sp[i,j,0]
        zBank=sp[i,j,1]
        zMax=sp[i,j,2]

        vips=sp[i,j,3:]

        if v/(dx*dx)>=vips[-1]:
            wl=(v/(dx*dx)-vips[-1])+zMax
        else:
            wl=numpy.interp(v/(dx*dx),numpy.concatenate(([0],vips)),[zMin,zBank,zBank+1,zBank+5,zMax])

    return wl





# Add a volume to a given cell
def addVol(wlg,i,j,sp,dv):
    v=volFromWl(wlg[i,j],i,j,sp)
    v+=dv
    wlg[i,j]=wlFromVol(v,i,j,sp)

def resample(wl,xll,yll,dx,dtmFileName):
    xsz,ysz=wl.shape

    dtm,dtmCellSize,dtmXsz,dtmYsz,dtmXll,dtmYll=fileIO.readScalarGrid(dtmFileName)

    depth=numpy.zeros((dtmXsz,dtmYsz),dtype=arrayType)
    wlGrid=numpy.zeros((dtmXsz,dtmYsz),dtype=arrayType)-9999.

    for i in range(dtmXsz):
        for j in range(dtmYsz):
            xc=dtmXll+i*dtmCellSize+0.5*dtmCellSize
            yc=dtmYll+j*dtmCellSize+0.5*dtmCellSize

            iwl=int((xc-xll)/dx)
            jwl=int((yc-yll)/dx)

            if xc<xll or yc<yll:
                continue

            if iwl>=0 and iwl<xsz and jwl>=0 and jwl<ysz:
                depth[i,j]=max(0,wl[iwl,jwl]-dtm[i,j])
                if wl[iwl,jwl]>dtm[i,j]:
                    wlGrid[i,j]=wl[iwl,jwl]
            else:
                depth[i,j]=0.

    depth[numpy.where(depth<=0)]=-9999.

    return depth,wlGrid,dtmXll,dtmYll,dtmCellSize

def nint(f):
    return int(f+0.5)

def resample2(wl,v,xll,yll,dx,dtmFileName,ndv=None,ndr=None):
    xsz,ysz=wl.shape

    dtm,dtmCellSize,dtmXsz,dtmYsz,dtmXll,dtmYll=fileIO.readScalarGrid(dtmFileName)

    if ndv is not None and ndr is not None:
        dtm[numpy.where(dtm==ndv)]=ndr

    dtm[numpy.where(numpy.isnan(dtm))]=ndr

    depth=numpy.zeros((dtmXsz,dtmYsz),dtype=arrayType)
    wlGrid=numpy.zeros((dtmXsz,dtmYsz),dtype=arrayType)-9999.

    dtmWindow=numpy.zeros((nint(dx/dtmCellSize),nint(dx/dtmCellSize)),dtype=arrayType)
    depthWindow=numpy.zeros((nint(dx/dtmCellSize),nint(dx/dtmCellSize)),dtype=arrayType)
    wlWindow=numpy.zeros((nint(dx/dtmCellSize),nint(dx/dtmCellSize)),dtype=arrayType)

    for i in range(xsz):
        for j in range(ysz):
            if v[i,j]<1e-3:
                continue

            xCell1=xll+i*dx
            yCell1=yll+j*dx
            xCell2=xCell1+dx
            yCell2=yCell1+dx

            xi0=int((xCell1-dtmXll)/dtmCellSize)
            xi1=xi0+nint(dx/dtmCellSize)

            yi0=int((yCell1-dtmYll)/dtmCellSize)
            yi1=yi0+nint(dx/dtmCellSize)

            try:
                dtmWindow[:,:]=dtm[xi0:xi1,yi0:yi1]
            except:
                print(xCell1, xCell2)
                print(yCell1, yCell2)
                print(xi0, yi0, xi1, yi1)
                print()
                print((xCell1-dtmXll)/dtmCellSize)
                print((xCell2-dtmXll)/dtmCellSize)

                print((yCell1-dtmYll)/dtmCellSize)
                print((yCell2-dtmYll)/dtmCellSize)

            depthWindow=wl[i,j]-dtmWindow
            depthWindow[numpy.where(depthWindow<0)]=-9999.
            wlWindow[:,:]=wl[i,j]
            wlWindow[numpy.where(depthWindow==-9999.)]=-9999.

            depth[xi0:xi1,yi0:yi1]=depthWindow
            wlGrid[xi0:xi1,yi0:yi1]=wlWindow


    return depth,wlGrid,dtmXll,dtmYll,dtmCellSize

def formatTime(t):
    nd=int(t/86400.)
    nh=int((t-nd*86400.)/3600.)
    nm=int((t-nd*86400-nh*3600)/60.)
    ns=int(t-nd*86400-nh*3600-nm*60.)

    return "%02id:%02ih:%02im:%02is"%(nd,nh,nm,ns)


def loadCppLib(libPath):

    if arrayType==numpy.float64:
        arrayArgType=ndpointer(ctypes.c_double)
        scalarType=ctypes.c_double
        staticLibName=libPath
    else: # Must be float 32
        arrayArgType=ndpointer(ctypes.c_float)
        scalarType=ctypes.c_float
        staticLibName=libPath

    lib=ctypes.cdll.LoadLibrary(staticLibName)

    cppSum=lib.sum
    cppSum.restype=scalarType
    cppSum.argtypes=[arrayArgType, \
        ctypes.c_int, ctypes.c_int]

    cppCalcFlow=lib.calcFlow
    cppCalcFlow.restype=scalarType
    cppCalcFlow.argtypes=[scalarType,scalarType,scalarType,\
        scalarType,scalarType,scalarType,scalarType]

    cppCalcFlowGrid=lib.calcFlowGrid
    cppCalcFlowGrid.argtypes=[arrayArgType, \
        arrayArgType, \
        arrayArgType, \
        scalarType, scalarType, ctypes.c_int, ctypes.c_int,
        arrayArgType, \
        arrayArgType,\
        arrayArgType]

    cppDryCheck=lib.dryCheck
    cppDryCheck.restype=ctypes.c_int
    cppDryCheck.argtypes=[arrayArgType,\
        arrayArgType,arrayArgType,\
        scalarType,ctypes.c_int,ctypes.c_int,\
        ndpointer(ctypes.c_int),ndpointer(ctypes.c_int),arrayArgType,\
        ctypes.c_int]

    cppDryCheckDiagnostic=lib.dryCheckDiagnostic
    cppDryCheckDiagnostic.restype=ctypes.c_int
    cppDryCheckDiagnostic.argtypes=[arrayArgType,\
        arrayArgType,arrayArgType,\
        scalarType,ctypes.c_int,ctypes.c_int,\
        ndpointer(ctypes.c_int),ndpointer(ctypes.c_int),arrayArgType,\
        ctypes.c_int,arrayArgType]


    cppTimeStep=lib.timeStep
    cppTimeStep.argtypes=[arrayArgType,\
        arrayArgType,arrayArgType,
        scalarType,ctypes.c_int,ctypes.c_int,
        ndpointer(ctypes.c_int),ndpointer(ctypes.c_int),arrayArgType,\
        ctypes.c_int]

    cppWlFromVolGrid=lib.wlFromVolGrid
    cppWlFromVolGrid.argtypes=[arrayArgType,\
        arrayArgType,arrayArgType,\
        ctypes.c_int,ctypes.c_int,ctypes.c_bool,scalarType]

    cppConveyanceParameters=lib.conveyanceParameters
    cppConveyanceParameters.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,\
        arrayArgType,scalarType,ctypes.c_int,ctypes.c_int,\
        scalarType,\
        arrayArgType,\
        ctypes.c_bool,scalarType,arrayArgType,ctypes.c_int,ctypes.c_int,\
        ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,\
        scalarType,scalarType,scalarType,scalarType,arrayArgType]


    cppMaxVolGrid=lib.maxVolGrid
    cppMaxVolGrid.argtypes=[arrayArgType, \
                            arrayArgType,\
                            arrayArgType,\
                            arrayArgType,\
                            arrayArgType,\
                            arrayArgType,\
                            ctypes.c_int,ctypes.c_int]


    cppResample2=lib.resample2
    cppResample2.argtypes=[arrayArgType, \
        arrayArgType, \
        scalarType,scalarType,scalarType,\
        ctypes.c_int,ctypes.c_int,\
        arrayArgType,scalarType,\
        ctypes.c_int,ctypes.c_int,\
        scalarType,scalarType,\
        arrayArgType,arrayArgType]

    cppCalcStorageParameters=lib.calcStorageParameters
    cppCalcStorageParameters.argtypes=[scalarType,scalarType,\
        scalarType,scalarType, \
        arrayArgType,ctypes.c_int,ctypes.c_int,\
        scalarType,scalarType,scalarType,\
        arrayArgType, \
        ctypes.c_bool,arrayArgType,ctypes.c_int,ctypes.c_int,\
        scalarType,scalarType,scalarType,scalarType,scalarType,\
        scalarType,scalarType,scalarType]

    cppFlowPaths=lib.flowPaths
    cppFlowPaths.argtypes=[scalarType,scalarType,scalarType,\
        ctypes.c_int,ctypes.c_int,\
        arrayArgType,scalarType,\
        ctypes.c_int,ctypes.c_int,\
        scalarType,scalarType,\
        ctypes.c_char_p,\
        arrayArgType,arrayArgType,\
        ctypes.c_int]

    cppResample3=lib.resample3
    cppResample3.argtypes=[arrayArgType, \
        arrayArgType, \
        arrayArgType, \
        arrayArgType, \
        arrayArgType, arrayArgType, scalarType, \
        scalarType,scalarType,scalarType,\
        ctypes.c_int,ctypes.c_int,\
        arrayArgType,scalarType,\
        ctypes.c_int,ctypes.c_int,\
        scalarType,scalarType,\
        arrayArgType,arrayArgType]


    cppLazyFlowPaths=lib.lazyFlowPaths
    cppLazyFlowPaths.argtypes=[\
        scalarType,scalarType,scalarType,
        ctypes.c_int,ctypes.c_int,
        arrayArgType,scalarType,ctypes.c_int,ctypes.c_int,
        scalarType,scalarType,
        arrayArgType,arrayArgType,arrayArgType,scalarType,scalarType,scalarType]

    cppWlFill=lib.fillWlGrid
    cppWlFill.argtypes=[arrayArgType,arrayArgType,ctypes.c_int,ctypes.c_int,scalarType, \
        ctypes.c_int,ctypes.c_int,scalarType]


    cppBurnFlowPaths=lib.burnFlowPaths
    cppBurnFlowPaths.argtypes=[arrayArgType,arrayArgType,ctypes.c_int,ctypes.c_int]

    cppMakeWlGrid=lib.makeWlGrid
    cppMakeWlGrid.argtypes=[arrayArgType,arrayArgType,arrayArgType,\
        ctypes.c_int,ctypes.c_int,scalarType]

    cppClipZero=lib.clipZero
    cppClipZero.argtypes=[arrayArgType,ctypes.c_int,ctypes.c_int]


    cppScsAdditionalRunoff=lib.scsAdditionalRunoff
    cppScsAdditionalRunoff.argtypes=[arrayArgType,arrayArgType,arrayArgType,arrayArgType,\
        ctypes.c_int,ctypes.c_int]


    cppCalcFlowEdges=lib.calcFlowEdges
    cppCalcFlowEdges.restype=scalarType
    cppCalcFlowEdges.argtypes=[arrayArgType,arrayArgType,arrayArgType,arrayArgType,\
        scalarType,ctypes.c_int,ctypes.c_int,arrayArgType,arrayArgType,arrayArgType]

    cppCheckLicence=lib.checkLicence
    return cppCalcFlow, cppCalcFlowGrid, cppDryCheck, cppTimeStep, \
        cppWlFromVolGrid, cppConveyanceParameters, cppMaxVolGrid, \
        cppResample2,cppResample3,cppFlowPaths,cppSum,cppCalcStorageParameters,\
        cppLazyFlowPaths, cppWlFill, cppBurnFlowPaths, cppMakeWlGrid, cppClipZero, cppDryCheckDiagnostic,\
        cppScsAdditionalRunoff,cppCalcFlowEdges,cppCheckLicence




def __lowerChannelCells(topoProfile,dxt,width,depth,nChan,nFP):
    numCellsLower=int(width/dxt)
    minZ=min(topoProfile)



    topoProfile.sort()

    nList=[nFP]*len(topoProfile)

    # Find lowest cells to change
    for i in range(numCellsLower):
        topoProfile[i]=minZ-depth
        nList[i]=nChan

    return nList

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Utility function to tidy up a load of if/else code
def getNeighbouringWl(i1, j1, i2, j2, qx, qy, flowThreshold, g):

    xsz, ysz = g.shape

    # If neighbouring point falls outside grid, use cell value
    if i2 <0 or i2 >= xsz or j2 < 0 or j2 >= ysz:
        return g[i1,j1]

    # If neighbouring cell is inactive, use cell value
    if g[i2,j2] == -9999:
        return g[i1,j1]

    # x-neighbour
    if i2 != i1:
        if i2 < i1:
            q = qx[i1,j1]
        else:
            q = qx[i1 + 1,j1]

        if abs(q) > flowThreshold:
            return g[i2,j2]
        else:
            return g[i1,j1]

    # Must be j-neighbour
    if j2 < j1:
        q = qy[i1, j1]
    else:
        q = qy[i1, j1 + 1]

    if abs(q) > flowThreshold:
        return g[i2, j2]
    else:
        return g[i1, j1]

    return # Should never get here

# Wrapper for interpolateTile for threading; unpacks arguments, saves to TIFF
def interpolateTileWrapper(argTuple):

    i = argTuple[0]
    j = argTuple[1]
    wlGrid = argTuple[2]
    dtmFileName = argTuple[3]
    flowX = argTuple[4]
    flowY = argTuple[5]
    flowThreshold = argTuple[6]
    method = argTuple[7]
    outputFilePathRoot = argTuple[8]
    saveWl = argTuple[9]
    xsz = argTuple[10]
    ysz = argTuple[11]
    cellSize = argTuple[12]
    xll = argTuple[13]
    yll = argTuple[14]

    x0 = xll + i * cellSize
    x1 = x0 + cellSize
    y0 = yll + j * cellSize
    y1 = y0 + cellSize

    dtmFileGeoGrid = fileIO.geoGrid(str(dtmFileName), objOnly=True)
    dx = dtmFileGeoGrid.dx

    xi0 = int((x0 - dtmFileGeoGrid.xll) / dtmFileGeoGrid.dx)
    xi1 = int((x1 - dtmFileGeoGrid.xll) / dtmFileGeoGrid.dx) - 1
    yi0 = int((y0 - dtmFileGeoGrid.yll) / dtmFileGeoGrid.dx)
    yi1 = int((y1 - dtmFileGeoGrid.yll) / dtmFileGeoGrid.dx) - 1

    windowXsz = xi1 - xi0 + 1
    windowYsz = yi1 - yi0 + 1

    dtmTile = dtmFileGeoGrid.obj.ReadAsArray(xoff=xi0, yoff=dtmFileGeoGrid.ysz - 1 - yi1, xsize=windowXsz,
                                             ysize=windowYsz).transpose().copy()
    dtmTile = numpy.array(dtmTile[:, ::-1], arrayType)


    depth, wl =  interpolateTile(i, j, wlGrid, dtmTile, xi0, yi0, dx, flowX, flowY, flowThreshold, method = method)

    tileString = f"{i:03d}_{j:03d}"
    fileIO.saveScalarGrid(depth, x0, y0, dtmFileGeoGrid.dx,
                          outputFilePathRoot + "_depth_" + tileString + ".tif")

    if saveWl:
        fileIO.saveScalarGrid(wl, x0, y0, dtmFileGeoGrid.dx,
                              outputFilePathRoot + "_wl_" + tileString + ".tif")

    return True


def  interpolateTile(i, j, wlg, dtmTile, xll, yll, dx, qx, qy, flowThreshold, method = 1):

    dxsz, dysz = dtmTile.shape

    if method == 1:
        depth = wlg[i, j] - dtmTile
        depth[depth < 0] = -9999

        wl = numpy.zeros((dxsz, dysz)) - 9999
        wl[wlg[i, j] > dtmTile] = wlg[i, j]

        return depth, wl

    xsz, ysz = wlg.shape

    dxsz2 = int(dxsz / 2)
    dysz2 = int(dysz / 2)

    wl = numpy.zeros((dxsz, dysz))
    depth = numpy.zeros((dxsz, dysz))

    y = numpy.vstack([numpy.linspace(0,0.5,dxsz2+1)[1:]] * dxsz2)
    x = y[:,:].transpose()

    z0 = wlg[i,j]

    # Need to process in 4 quadrants
    # North east
    z1 = getNeighbouringWl(i, j, i + 1, j, qx, qy, flowThreshold, wlg)
    z2 = getNeighbouringWl(i, j, i, j + 1, qx, qy, flowThreshold, wlg)
    wl[dxsz2:, dysz2:] = z0 + (z1 - z0) * x + (z2 - z0) * y

    # South east
    z1 = getNeighbouringWl(i, j, i + 1, j, qx, qy, flowThreshold, wlg)
    z2 = getNeighbouringWl(i, j, i, j - 1, qx, qy, flowThreshold, wlg)
    wl[dxsz2:, 0 : dysz2] = z0 + (z1 - z0) * x + (z2 - z0) * y[:,::-1]

    # North west
    z1 = getNeighbouringWl(i, j, i - 1, j, qx, qy, flowThreshold, wlg)
    z2 = getNeighbouringWl(i, j, i, j + 1, qx, qy, flowThreshold, wlg)
    wl[0:dxsz2, dysz2:] = z0 + (z1 - z0) * x[::-1,:] + (z2 - z0) * y

    # South west
    z1 = getNeighbouringWl(i, j, i - 1, j, qx, qy, flowThreshold, wlg)
    z2 = getNeighbouringWl(i, j, i, j - 1, qx, qy, flowThreshold, wlg)
    wl[0:dxsz2, 0:dysz2] = z0 + (z1 - z0) * x[::-1, :] + (z2 - z0) * y[:,::-1]

    # Calculate depth from water level and assign no data cells
    depth = wl - dtmTile
    depth[depth<=0] = -9999
    wl[depth == -9999] = -9999

    return depth, wl

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def saveResults(wlGrid,flowX,flowY,xsz,ysz, xll, yll, cellSize,
                flowThreshold, dtmFileName,
                       outputDirectory,outputPrefix,
                       threads = None, method = 1, saveWl = False):

    outputFilePathRoot=os.path.join(outputDirectory,outputPrefix)

    ticker=0

    if xsz>100:
        tickerStep=int(xsz * ysz / 100)
    else:
        tickerStep=1

    t1 = time.time()

    tileList = [] # Use this to track which TIFFs we need to merge later

    if threads is None:
        for i in range(xsz):
            for j in range(ysz):

                if (ticker % tickerStep) == 0:
                    pc = int(100. * ticker / (xsz * ysz))
                    print(f"{pc}%% {i}/{xsz},{j}/{ysz} " % (), end='')

                    if ticker > 0:
                        pcComplete = float(ticker) / (xsz * ysz)
                        pcToGo = 1. - pcComplete
                        t2 = time.time() - t1
                        projectedFinish = time.time() + pcToGo * (t2 / pcComplete)

                        if pcToGo * (t2 / pcComplete) < 86400:  # <1 day, report time only
                            projectedFinishString = time.strftime("%H:%M:%S", time.localtime(projectedFinish))
                        else:  # Report date too
                            projectedFinishString = time.strftime("%d/%m/%y %H:%M:%S", time.localtime(projectedFinish))

                        print(projectedFinishString, end='')

                    print("...", end='', flush=True)
                    sys.stdout.flush()

                ticker += 1

                if wlGrid[i,j] == -9999: # This should skip inactive and dry cells
                    continue

                x0 = xll + i * cellSize
                x1 = x0 + cellSize
                y0 = yll + j * cellSize
                y1 = y0 + cellSize

                dtmFileGeoGrid = fileIO.geoGrid(str(dtmFileName), objOnly=True)

                xi0 = int((x0 - dtmFileGeoGrid.xll) / dtmFileGeoGrid.dx)
                xi1 = int((x1 - dtmFileGeoGrid.xll) / dtmFileGeoGrid.dx) - 1
                yi0 = int((y0 - dtmFileGeoGrid.yll) / dtmFileGeoGrid.dx)
                yi1 = int((y1 - dtmFileGeoGrid.yll) / dtmFileGeoGrid.dx) - 1

                windowXsz = xi1 - xi0 + 1
                windowYsz = yi1 - yi0 + 1

                dtmTile = dtmFileGeoGrid.obj.ReadAsArray(xoff=xi0, yoff=dtmFileGeoGrid.ysz - 1 - yi1, xsize=windowXsz,
                                                ysize=windowYsz).transpose().copy()
                dtmTile = numpy.array(dtmTile[:, ::-1], arrayType)

                depthGrid, wl = interpolateTile(i, j, wlGrid,
                                                dtmTile, xi0, yi0, dtmFileGeoGrid.dx,
                                                flowX, flowY, flowThreshold, method)

                tileString = f"{i:03d}_{j:03d}"
                fileIO.saveScalarGrid(depthGrid, x0, y0, dtmFileGeoGrid.dx,
                                      outputFilePathRoot + "_depth_" + tileString + ".tif")

                tileList.append((i,j))

                if saveWl:
                    fileIO.saveScalarGrid(wl, x0, y0, dtmFileGeoGrid.dx,
                                          outputFilePathRoot + "_wl_" + tileString + ".tif")
    else: # Threaded
        pool = ThreadPool(threads)
        funcArgList = []
        returnValues = []
        for i, j in product(range(xsz), range(ysz)):
            if wlGrid[i, j] == -9999:  # This should skip inactive and dry cells
                continue

            tileList.append((i, j))

            funcArgList.append((i, j, wlGrid, dtmFileName, flowX, flowY, flowThreshold, method, outputFilePathRoot,
                                saveWl, xsz, ysz, cellSize, xll, yll))

        counter = 0

        tickerStep=int(len(tileList) / 100)
        ticker = 0


        for ret in pool.imap(interpolateTileWrapper, funcArgList, chunksize=10):

            if (ticker%tickerStep)==0:
                print("%i%% "%(100.*ticker/(len(tileList))), end='')

                if ticker > 0:
                    pcComplete = float(ticker) / (len(tileList))
                    pcToGo = 1. - pcComplete
                    t2 = time.time() - t1
                    projectedFinish = time.time() + pcToGo * (t2 / pcComplete)
                    projectedFinishString = time.strftime("%H:%M:%S", time.localtime(projectedFinish))
                    print(projectedFinishString, end='')

                print("...", end='')
                sys.stdout.flush()

            ticker += 1


    # Build VRTs for grid outputs for converting to TIFF later - this is much quicker than
    # adding TIFFs individually using gdal_merge
    print()
    print("Generating VRT(s)...",)

    # Depths
    vrtCommand=['gdalbuildvrt']
    vrtCommand += ['-vrtnodata', '-9999']
    vrtCommand += [outputFilePathRoot+'_depth.vrt']

    for i, j in tileList:
        tileString=f"{i:03d}_{j:03d}"
        vrtCommand.append(outputFilePathRoot+"_depth_"+tileString+".tif")

    call(vrtCommand)

    if saveWl:
        vrtCommand = ['gdalbuildvrt']
        vrtCommand += ['-vrtnodata', '-9999']
        vrtCommand += [outputFilePathRoot + '_wl.vrt']

        for i, j in tileList:
            tileString = f"{i:03d}_{j:03d}"
            vrtCommand.append(outputFilePathRoot + "_wl_" + tileString + ".tif")

        call(vrtCommand)

    print("Converting to TIFF(s)...",)
    mergeCommand = ['gdal_translate']
    mergeCommand+= ['-of', 'GTiff']
    mergeCommand+= ['-co', 'COMPRESS=LZW']
    mergeCommand+= ['-co', 'BIGTIFF=YES']
    mergeCommand+= [outputFilePathRoot+"_depth.vrt"]
    mergeCommand+= [outputFilePathRoot+"_depth.tif"]

    call(mergeCommand)

    if saveWl:
        mergeCommand = ['gdal_translate']
        mergeCommand += ['-of', 'GTiff']
        mergeCommand += ['-co', 'COMPRESS=LZW']
        mergeCommand += ['-co', 'BIGTIFF=YES']
        mergeCommand += [outputFilePathRoot + "_wl.vrt"]
        mergeCommand += [outputFilePathRoot + "_wl.tif"]

        call(mergeCommand)

    # Remove temporary geotiffs
    print("Deleting temporary files...",)
    for i, j in tileList:
        tileString=f"{i:03d}_{j:03d}"
        os.remove(outputFilePathRoot+"_depth_"+tileString+".tif")
        if saveWl:
            os.remove(outputFilePathRoot + "_wl_" + tileString + ".tif")

    os.remove(outputFilePathRoot+'_depth.vrt')
    if saveWl:
        os.remove(outputFilePathRoot + '_wl.vrt')

    return
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def deleteFileComplete(progName,str):
    fileName=os.path.splitext(os.path.basename(progName))[0]

    fileName+='_'+str+'.txt'

    try:
        os.remove(fileName)
    except:
        pass

    return

def notifyFileComplete(progName,str):
    fileName=os.path.splitext(os.path.basename(progName))[0]

    fileName+='_'+str+'.txt'

    f=open(fileName,'w')
    f.close()

    return


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def uncompressGeoTiff(fName,tiled=False,extents=None):
    tmpFileName=tempfile._get_candidate_names().next()+'.tiff'

    ogrCommand=['gdal_translate']
    ogrCommand+=['-co','BIGTIFF=YES']
#    ogrCommand+=['-co','COMPRESS=NO']


    if tiled:
        ogrCommand+=['-co','TILED=YES']

    if extents is not None:
        ogrCommand+=['-projwin']
        ogrCommand+=["%f"%extents[0]]
        ogrCommand+=["%f"%extents[3]]
        ogrCommand+=["%f"%extents[2]]
        ogrCommand+=["%f"%extents[1]]

    ogrCommand+=[fName]
    ogrCommand+=[tmpFileName]

    call(ogrCommand)

    return tmpFileName


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def clipRasterPoly(rasterFileName,polyFileName):
    ogrCommand=['gdal_rasterize']
    ogrCommand+=['-burn','-9999']
    ogrCommand+=['-i']
    ogrCommand+=[polyFileName]
    ogrCommand+=[rasterFileName]

    call(ogrCommand)

    return

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def calcTimeStep(wlGrid,zGrid,dx):
    depth=wlGrid-zGrid
    depth[numpy.where(zGrid==-9999)]=0.

    if depth.max() > 0:
        return 0.7*dx/numpy.sqrt(9.81*depth.max())
    else: # Assume depth = 1 if no depth yet
        return 0.7*dx/numpy.sqrt(9.81)





# Test linear interpolation function
if __name__ == "__main__":

        dtm = numpy.random.rand(500,500)
        dtmGeoGrid = fileIO.geoGrid(dtm, 500, 500, 0, 0, 2)

        wlGrid = numpy.array([
            [1, 2, 3],
            [2, 4, 4],
            [3, 5, 5],
        ]).transpose()[:,::-1]

        depth, wl =  interpolateTile(0, 1, wlGrid, dtmGeoGrid)

        fileIO.saveScalarGrid(depth,0,0,2,"depth_test.tif")
        fileIO.saveScalarGrid(wl,0,0,2,"wl_test.tif")

        pass
