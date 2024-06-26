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
    outCsvt.write("\"Integer\",\"String\",\"Real\"\n")
    outCsvt.close()

    csvFile=open(fileName,"w")

    id=0

#    csvFile.write("id;wkt;elevation;rainfall;runoff;cn\n")
    csvFile.write("id;wkt;elevation\n")

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

            csvFile.write("%i;%s;%f\n"%(id,wktString,sp[i,j,0]))
#            csvFile.write("%i;%s;%f;%f;%f;%f\n"%(id,wktString,sp[i,j,0],rain[i,j],runoff[i,j],cn[i,j]))
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
def gridFlowSetupTiled(dtmFileName,xll,yll,cellSize,xsz,ysz,nChan,nFP,
    nFileName=None,
    plotNamePrefix=None, outputPrefix=None,
    ndv=None,ndr=None,conveyanceFunc=None,storageFunc=None):

    if plotNamePrefix is None:
        plotNamePrefix=""

    if outputPrefix is None:
        outputPrefix=""

    dtmFileObj,dtmCellSize,dtmXsz,dtmYsz,dtmXll,dtmYll=fileIO.readScalarGridObj(dtmFileName)

    if nFileName is not None:
        nFileObj,lcCellSize,lcXsz,lcYsz,lcXll,lcYll=fileIO.readScalarGridObj(nFileName)

    # Might need to fettle no data value
    ndv=arrayType(ndv)

    convParX=numpy.zeros((xsz+1,ysz,7),dtype=arrayType)-9999.
    convParY=numpy.zeros((xsz,ysz+1,7),dtype=arrayType)-9999.

    storagePar=numpy.zeros((xsz,ysz,5),dtype=arrayType)-9999.

    plotName=None

    returnArray=numpy.zeros(7,dtype=arrayType)
    tmp5=numpy.zeros(5,dtype=arrayType)

    ticker=0

    if xsz>100:
        tickerStep=int(xsz/100)
    else:
        tickerStep=1

    for i in range(xsz):
        if (ticker%tickerStep)==0: print("%i%% ..."%(100.*ticker/xsz), end='')
        sys.stdout.flush()
        ticker+=1
        for j in range(ysz):

            x0=xll+i*cellSize
            x1=x0+cellSize
            y0=yll+j*cellSize
            y1=y0+cellSize

            xi0=int((x0-dtmXll)/dtmCellSize)
            xi1=int((x1-dtmXll)/dtmCellSize)
            yi0=int((y0-dtmYll)/dtmCellSize)
            yi1=int((y1-dtmYll)/dtmCellSize)

            windowXsz=xi1-xi0+1
            windowYsz=yi1-yi0+1

            if xi0>=0 and xi1<dtmXsz-1 and yi0>=0 and yi1<dtmYsz-1: # Within DTM extent - grab window
                dtmWindow=dtmFileObj.ReadAsArray(xoff=xi0,yoff=dtmYsz-yi1,xsize=windowXsz,ysize=windowYsz).transpose().copy()
                dtmWindow=numpy.array(dtmWindow[:,::-1],arrayType)

                dtmWindow[numpy.where(dtmWindow==ndv)]=ndr
                dtmWindow[numpy.where(dtmWindow<-1e6)]=ndr
                dtmWindow[numpy.where(numpy.isnan(dtmWindow))]=ndr

                if dtmWindow.max()==ndr:
                    continue

                # And get landcover window
                if nFileName is not None:
                    nWindow=nFileObj.ReadAsArray(xoff=xi0,yoff=dtmYsz-yi1,xsize=windowXsz,ysize=windowYsz).transpose().copy()
                    nWindow=nWindow[:,::-1]
                else:
                    nWindow=dtmWindow.copy()
                    nWindow[:,:]=nFP

                # X-direction

                conveyanceFunc(0,0,0,yi1-yi0,\
                    dtmWindow,dtmCellSize,windowXsz,windowYsz,nFP,\
                    returnArray,False,\
                    0.,dtmWindow,0,0,0,0, \
                    0,0,0,0,0,0,nWindow)

                convParX[i,j,:]=returnArray[:]

                # Y-direction
                conveyanceFunc(0,0,xi1-xi0,0,\
                    dtmWindow,dtmCellSize,xi1-xi0+1,yi1-yi0+1,nFP,\
                    returnArray,False,\
                    0.,dtmWindow,0,0,0,0, \
                    0,0,0,0,0,0,nWindow)

                convParY[i,j,:]=returnArray[:]

                # Cell storage
                if storageFunc is not None:
                    storageFunc(arrayType(x0),arrayType(y0),arrayType(x1),arrayType(y1),dtmWindow,\
                        windowXsz,windowYsz,arrayType(x0),arrayType(y0), \
                        arrayType(dtmCellSize),tmp5,\
                        False,dtmWindow,0,0,\
                        0,0,0,\
                        0,0,0,0,0)
                    storagePar[i,j,:]=tmp5
                else:
                    storagePar[i,j,:]=calcStorageParameters([(x0,y0),(x1,y1)],\
                        dtmWindow,x0,y0,dtmCellSize,plotName=plotName,csvOutput=False)

    print("Done.")

    fileIO.saveConveyanceParametersCSV(convParX,convParY,xll,yll,cellSize,\
        outputPrefix+"conveyanceParams.csv")
    fileIO.saveStorageParametersCSV(storagePar,xll,yll,cellSize,\
        outputPrefix+"storageParams.csv")

    return convParX, convParY, storagePar #, bankLevel


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def gridFlowSetupTiled2(dtmFileName,xll,yll,cellSize,xsz,ysz,nChan,nFP,
    plotNamePrefix=None, outputPrefix=None,
    ndv=None,ndr=None,conveyanceFunc=None,storageFunc=None):

    if plotNamePrefix is None:
        plotNamePrefix=""

    if outputPrefix is None:
        outputPrefix=""

    dtmFileObj,dtmCellSize,dtmXsz,dtmYsz,dtmXll,dtmYll=fileIO.readScalarGridObj(dtmFileName)

    # Might need to fettle no data value
    ndv=arrayType(ndv)


    convParX=numpy.zeros((xsz+1,ysz,6),dtype=arrayType)-9999.
    convParY=numpy.zeros((xsz,ysz+1,6),dtype=arrayType)-9999.

    storagePar=numpy.zeros((xsz,ysz,5),dtype=arrayType)-9999.

    count=0
    plotName=None

    returnArray=numpy.zeros(6,dtype=arrayType)
    tmp5=numpy.zeros(5,dtype=arrayType)
    tmp7=numpy.zeros(7,dtype=arrayType)

    ticker=0

    stepSize=10 # Load bigger tiles - should speed things up

    if xsz>100:
        tickerStep=int(xsz/100)
    else:
        tickerStep=1

    for i in range(0,xsz,stepSize):
        print(i, xsz)
        if (ticker%tickerStep)==0: print("%i%% ..."%(100.*ticker/xsz),end='')
        sys.stdout.flush()
        ticker+=1
        for j in range(0,ysz,stepSize):

            # Dimensions of tile to load from file
            x0=xll+i*cellSize
            x1=x0+cellSize*stepSize
            y0=yll+j*cellSize
            y1=y0+cellSize*stepSize

            xi0=int((x0-dtmXll)/dtmCellSize)
            xi1=int((x1-dtmXll)/dtmCellSize)
            yi0=int((y0-dtmYll)/dtmCellSize)
            yi1=int((y1-dtmYll)/dtmCellSize)

            # Make sure doesn't extend beyond file
            xi1=min(xi1,dtmXsz-1)
            yi1=min(yi1,dtmYsz-1)

            windowXsz=xi1-xi0+1
            windowYsz=yi1-yi0+1

            # Now read in array
            dtmWindow=dtmFileObj.ReadAsArray(xoff=xi0,yoff=dtmYsz-yi1,xsize=windowXsz,ysize=windowYsz).transpose().copy()

            dtmWindow=numpy.array(dtmWindow[:,::-1],arrayType)
            dtmWindow[numpy.where(dtmWindow==ndv)]=ndr
            dtmWindow[numpy.where(numpy.isnan(dtmWindow))]=ndr

            # And calculate parameters for model cells
            for ii in range(i,i+stepSize):
                for jj in range(j,j+stepSize):

                    if ii>=xsz or jj>=ysz:
                        continue

                    xx0=x0+(ii-i)*cellSize
                    xx1=xx0+cellSize
                    yy0=y0+(jj-j)*cellSize
                    yy1=yy0+cellSize

                    xxi0=int((xx0-x0)/dtmCellSize)
                    xxi1=int((xx1-x0)/dtmCellSize)
                    yyi0=int((yy0-y0)/dtmCellSize)
                    yyi1=int((yy1-y0)/dtmCellSize)


                    cellWindow=dtmWindow[xxi0:xxi1+1,yyi0:yyi1+1]

                    cwXsz,cwYsz=cellWindow.shape

                    # X-direction
                    conveyanceFunc(0,0,0,cwYsz-1,\
                        cellWindow,dtmCellSize,cwXsz,cwYsz,nFP,\
                        returnArray,False,\
                        0.,cellWindow,0,0,0,0, \
                        0,0,0,0,0,0)

                    convParX[ii,jj,:]=returnArray[:]

                    # Y-direction
                    conveyanceFunc(0,0,cwXsz-1,0,\
                        cellWindow,dtmCellSize,cwXsz,cwYsz,nFP,\
                        returnArray,False,\
                        0.,cellWindow,0,0,0,0, \
                        0,0,0,0,0,0)

                    convParY[ii,jj,:]=returnArray[:]

                    # Cell storage

                    storageFunc(arrayType(xx0),arrayType(yy0),arrayType(xx1),arrayType(yy1),cellWindow,\
                        cwXsz,cwYsz,arrayType(xx0),arrayType(yy0), \
                        arrayType(dtmCellSize),tmp5,\
                        False,cellWindow,0,0,\
                        0,0,0,\
                        0,0,0,0,0)
                    storagePar[ii,jj,:]=tmp5



#                    storagePar[i,j,:]=calcStorageParameters([(xx0,yy0),(xx1,yy1)],\
#                        cellWindow,xx0,yy0,dtmCellSize,plotName=plotName,csvOutput=False)

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
def saveResults(volGrid,wlGrid,flowX,flowY,storagePar,xsz,ysz,cellSize,xll,yll,
                       defaultDepth,flowThreshold,channel,flowPathOutput,
                       dtmFileName,noDataValue,noDataReplacement,
                       outputDirectory,outputPrefix,
                       resampleFunction,lfpFunction,
                       extendWlGrid,cppBurnFlowPaths,cppMakeWlGrid,cppWlFill,cppClipZero,
                       saveCsv=True,zeroPolyList=None):

    if flowPathOutput is not None:
        defaultDepth=arrayType(defaultDepth)
        flowThreshold=arrayType(flowThreshold)
        dWeight=arrayType(0.25)

    # Save zeroPolyList as csv file, and get poly mask extents - don't
    # need to process these for all tiles
    if zeroPolyList is not None:
        polyXmin=1e20
        polyXmax=-1e20
        polyYmin=1e20
        polyYmax=-1e20

        for poly in zeroPolyList:
            b=poly.bounds
            polyXmin=min(b[0],polyXmin)
            polyYmin=min(b[1],polyYmin)

            polyXmax=max(b[2],polyXmax)
            polyYmax=max(b[3],polyYmax)

            tmpWktFileName=tempfile._get_candidate_names().next()+'.csv'

            wktFile=open(tmpWktFileName,"w")
            wktFile.write("id;wkt\n")
            count=0

            for poly in zeroPolyList:
                wktFile.write("%i;%s\n"%(count,poly.wkt))
                count+=1

            wktFile.close()

    outputFilePathRoot=os.path.join(outputDirectory,outputPrefix)

#    nExpansionIts=None;
#    nSmoothingIts=None;
#    expansionSlope=None;

    nExpansionIts=5;
    nSmoothingIts=5;
    expansionSlope=0.01;

    zMin=numpy.zeros((xsz,ysz),dtype=numpy.float32,order='C') # Should use array type
    zMax=numpy.zeros((xsz,ysz),dtype=numpy.float32,order='C')

    if channel:
        zMin[:,:]=storagePar[:,:,1]
        zMax[:,:]=storagePar[:,:,2]
    else:
        zMin[:,:]=storagePar[:,:,0]
        zMax[:,:]=storagePar[:,:,1]

    dtmObj,ddx,dtmXsz,dtmYsz,dxll,dyll=fileIO.readScalarGridObj(dtmFileName)

    tileSize=10000 # Tile to maximum of this size square
    numTilesX=int(dtmXsz/tileSize)+1
    numTilesY=int(dtmYsz/tileSize)+1

    for iTile in range(numTilesX):
        for jTile in range(numTilesY):

            print("Processing tile %i/%i,%i/%i"%(iTile+1,numTilesX,jTile+1,numTilesY))

            x0=dxll+iTile*tileSize*ddx
            x1=min(x0+tileSize*ddx,x0+ddx*dtmXsz)

            y1=dyll+dtmYsz*ddx-jTile*tileSize*ddx
            y0=max(y1-tileSize*ddx,dyll)

            xoff=iTile*tileSize
            yoff=jTile*tileSize

            tileXsize=min(tileSize,dtmXsz-iTile*tileSize)
            tileYsize=min(tileSize,dtmYsz-jTile*tileSize)

            dtmTile=dtmObj.ReadAsArray(xoff=xoff,yoff=yoff,xsize=tileXsize,ysize=tileYsize).transpose().copy()
            dtmTile=numpy.array(dtmTile[:,::-1],arrayType)

            if noDataValue is not None and noDataReplacement is not None:
                replaceArrayVals(dtmTile,arrayType(noDataValue),noDataReplacement)

            wlGrid2=numpy.zeros((tileXsize,tileYsize),dtype=arrayType)-9999.
            depthGrid=numpy.zeros((tileXsize,tileYsize),dtype=arrayType)-9999.

            resampleFunction(wlGrid,volGrid,flowX,flowY,zMin,zMax, \
                flowThreshold,\
                xll,yll,cellSize,xsz,ysz,dtmTile,ddx,tileXsize,tileYsize,\
                x0,y0,wlGrid2,depthGrid)

            tileString="tile%02i%02i"%(iTile,jTile)

            fileIO.saveScalarGrid(depthGrid,x0,y0,ddx,\
                outputFilePathRoot+"_depth_"+tileString+".tif")
            ###################################################################
            if zeroPolyList is not None and \
                polyXmin<x1 and polyXmax>x0 and polyYmin<y1 and polyYmax>y0:
                print("Zeroing polys...")
                reservoirs.maskZeroPoly(tmpWktFileName,\
                    os.path.join(outputDirectory,\
                    outputPrefix+"_depth_"+tileString+".tif"))
            ###################################################################
            fileIO.saveScalarGrid(wlGrid2,x0,y0,ddx,\
                outputFilePathRoot+"_wl_"+tileString+".tif")

            lfpGrid=numpy.zeros((tileXsize,tileYsize),dtype=arrayType)

            lfpFunction(xll,yll,cellSize,xsz,ysz,dtmTile,ddx,tileXsize,tileYsize,x0,y0,\
                lfpGrid,flowX,flowY,flowThreshold,defaultDepth,dWeight)

            fileIO.saveScalarGrid(lfpGrid,x0,y0,ddx,\
                outputFilePathRoot+"_lfp_"+tileString+".tif")

            # Merge depth and flowpath grids
            repList=numpy.where((depthGrid<defaultDepth) & (lfpGrid>0))
            depthGrid[repList]=defaultDepth

            fileIO.saveScalarGrid(depthGrid,x0,y0,ddx,\
                outputFilePathRoot+"_merge_"+tileString+".tif")

           # Extend WL grid
            if extendWlGrid:
                # This is memory critical - use 3 arrays only
#                arr1=dtm
#                arr2=depthGrid
#                arr3=lfpGrid

                # Insert flow paths, if present
                if lfpGrid is not None:
                    cppBurnFlowPaths(depthGrid,lfpGrid,tileXsize,tileYsize)

                wlGrid2[:,:]=0

                cppMakeWlGrid(dtmTile,depthGrid,wlGrid2,dtmXsz,dtmYsz,0.1)

                depthGrid[:,:]=0

                # Default parameters for expanding/filling water level
                if nExpansionIts is None:
                    nExpansionIts=int(0.5*cellSize/ddx)

                if nSmoothingIts is None:
                    nSmoothingIts=5

                if expansionSlope is None:
                    expansionSlope=2.*defaultDepth/cellSize

                wlGrid3=numpy.zeros((tileXsize,tileYsize),dtype=arrayType)-9999.
                depthGrid3=numpy.zeros((tileXsize,tileYsize),dtype=arrayType)-9999.

                cppWlFill(wlGrid2,wlGrid3,dtmXsz,dtmYsz,ddx,nExpansionIts,nSmoothingIts,expansionSlope)

                fileIO.saveScalarGrid(wlGrid3,dxll,dyll,ddx,\
                    outputFilePathRoot+"_wl_fill_"+tileString+".tif")


                depthGrid3[:,:]=wlGrid3-dtmTile

                cppClipZero(depthGrid3,dtmXsz,dtmYsz)

                fileIO.saveScalarGrid(depthGrid3,dxll,dyll,ddx,\
                    outputFilePathRoot+"_depth_fill_"+tileString+".tif")


    # Build VRTs for grid outputs
    print("Generating VRTs...",)

    # Depths
    vrtCommand=['gdalbuildvrt']
    vrtCommand.append(outputFilePathRoot+'_depth.vrt')

    for iTile in range(numTilesX):
        for jTile in range(numTilesY):
            tileString="tile%02i%02i"%(iTile,jTile)
            vrtCommand.append(outputFilePathRoot+"_depth_"+tileString+".tif")

    call(vrtCommand)

    # Water levels
    vrtCommand=['gdalbuildvrt']
    vrtCommand.append(outputFilePathRoot+'_wl.vrt')

    for iTile in range(numTilesX):
        for jTile in range(numTilesY):
            tileString="tile%02i%02i"%(iTile,jTile)
            vrtCommand.append(outputFilePathRoot+"_wl_"+tileString+".tif")

    call(vrtCommand)

    # Flow paths
    vrtCommand=['gdalbuildvrt']
    vrtCommand.append(outputFilePathRoot+'_lfp.vrt')

    for iTile in range(numTilesX):
        for jTile in range(numTilesY):
            tileString="tile%02i%02i"%(iTile,jTile)
            vrtCommand.append(outputFilePathRoot+"_lfp_"+tileString+".tif")

    call(vrtCommand)

    # Merged depth/flow paths
    vrtCommand=['gdalbuildvrt']
    vrtCommand.append(outputFilePathRoot+'_merge.vrt')

    for iTile in range(numTilesX):
        for jTile in range(numTilesY):
            tileString="tile%02i%02i"%(iTile,jTile)
            vrtCommand.append(outputFilePathRoot+"_merge_"+tileString+".tif")

    call(vrtCommand)

    # Filled depth/wl
    if extendWlGrid:
        vrtCommand=['gdalbuildvrt']
        vrtCommand.append(outputFilePathRoot+'_depth_fill.vrt')

        for iTile in range(numTilesX):
            for jTile in range(numTilesY):
                tileString="tile%02i%02i"%(iTile,jTile)
                vrtCommand.append(outputFilePathRoot+"_depth_fill_"+tileString+".tif")

        call(vrtCommand)

        vrtCommand=['gdalbuildvrt']
        vrtCommand.append(outputFilePathRoot+'_wl_fill.vrt')

        for iTile in range(numTilesX):
            for jTile in range(numTilesY):
                tileString="tile%02i%02i"%(iTile,jTile)
                vrtCommand.append(outputFilePathRoot+"_wl_fill_"+tileString+".tif")

        call(vrtCommand)

    print("done.")

    if saveCsv:
        fileIO.saveVectorCSV(flowX,flowY,xll,yll,cellSize,\
            outputFilePathRoot+"_flow.csv",thresholdVal=flowThreshold)

        fileIO.saveScalarCSV(wlGrid,xll,yll,cellSize,\
            outputFilePathRoot+"_wl.csv", headerList=['WL'])

    # Delete temporary CSV/WKT file for reservoirs
    if zeroPolyList is not None:
        os.remove(tmpWktFileName)

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
