# -*- coding: utf-8 -*-
import fileIO
import numpy

def sampleRainfallFromGrid(rainfallGrid,xll,yll,cellSize,xsz,ysz):
    rGrid,rGridCellSize,rGridXsz,rGridYsz,rGridXll,rGridYll=fileIO.readScalarGrid(rainfallGrid)

    rGridOut=numpy.zeros((xsz,ysz),dtype=numpy.float32)

    for i in range(xsz):
        for j in range(ysz):
            x=xll+i*cellSize+cellSize/2
            y=yll+j*cellSize+cellSize/2

            xi=int((x-rGridXll)/rGridCellSize)
            yi=int((y-rGridYll)/rGridCellSize)

            rGridOut[i,j]=rGrid[xi,yi]

    return rGridOut


def sampleCnFromGrid(cnGridFileName,xll,yll,cellSize,xsz,ysz):

    cnGrid,cnGridCellSize,cnGridXsz,cnGridYsz,cnGridXll,cnGridYll=fileIO.readScalarGrid(cnGridFileName)

    cnGridOut=numpy.zeros((xsz,ysz),dtype=numpy.float32)

    ticker=0
    print "Sampling curve number to cells",

    cnGridXtr=cnGridXll+cnGridXsz*cnGridCellSize;
    cnGridYtr=cnGridYll+cnGridYsz*cnGridCellSize;

    for i in range(xsz):
        if (ticker%int(xsz/10))==0: print ".",
        ticker+=1
        for j in range(ysz):
            x0=xll+i*cellSize
            x1=x0+cellSize

            y0=yll+j*cellSize
            y1=y0+cellSize

            if x0<cnGridXll or x1>cnGridXtr or y0<cnGridYll or y1>cnGridYtr:
                cnGridOut[i,j]=0.
                continue

            xi0=int((x0-cnGridXll)/cnGridCellSize)
            xi1=int((x1-cnGridXll)/cnGridCellSize)

            yi0=int((y0-cnGridYll)/cnGridCellSize)
            yi1=int((y1-cnGridYll)/cnGridCellSize)

            subWindow=cnGrid[xi0:(xi1+1),yi0:(yi1+1)]

            if subWindow.min()==0:
                cnGridOut[i,j]=1e-6
            else:
                cnGridOut[i,j]=numpy.average(subWindow[numpy.where(subWindow>0)])

#            print i,j,subWindow.min(),subWindow.max(),cnGridOut[i,j]


    return cnGridOut


def scsRunoff(rainfall,cn):
    rainfallInches=rainfall/25.4

    S=1000./cn-10.

#    print rainfallInches, S


    if rainfallInches<=0.2*S:
        return 0

    return 25.4*(rainfallInches-0.2*S)**2/(rainfallInches+0.8*S)

def scsRunoffGrid(rainfall,cn):
    returnGrid=numpy.zeros(rainfall.shape)

    rainfallInches=rainfall/25.4

    S=1000./cn-10.

    returnGrid[numpy.where(rainfallInches<=0.2*S)]=0.0

    nzCells=numpy.where(rainfallInches>0.2*S)

    returnGrid[nzCells]=25.4*(rainfallInches[nzCells]-0.2*S[nzCells])**2/(rainfallInches[nzCells]+0.8*S[nzCells])

    return returnGrid

def scsAdditionalRunoff(r,dr,cn):

    S=25.4*(1000./cn-10.)

    if r<=0.2*S:
        return 0

    dQdI=2.*(r-0.2*S)/(r+0.8*S)-((r-0.2*S)/(r+0.8*S))**2

    return dQdI*dr

def scsAdditionalRunoffGrid(r,dr,cn):
    returnGrid=numpy.zeros(r.shape)

    S=25.4*(1000./cn-10.)

#    print (r+0.8*S).min()
#    print r.min()
#    print S.min()

    returnGrid=dr*(2.*(r-0.2*S)/(r+0.8*S)-((r-0.2*S)/(r+0.8*S))**2)
    returnGrid[numpy.where(r<=0.2*S)]=0.0

    return returnGrid

#def scsRunoffI(rainfall,cn):
#    rainfallInches=rainfall
#    S=1000./cn-10.
#
#    if rainfallInches<=0.2*S:
#        return 0
#
#    return (rainfallInches-0.2*S)**2/(rainfall+0.8*S)
#
#def scsAdditionalRunoffI(r,dr,cn):
#
#    S=(1000./cn-10.)
#
#    if r<=0.2*S:
#        return 0
#
#    dQdI=2.*(r-0.2*S)/(r+0.8*S)-((r-0.2*S)/(r+0.8*S))**2
#
#    return dQdI*dr