# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 09:52:48 2016

@author: mattyh
"""
import shapely
import shapely.wkt

from osgeo import gdal
from osgeo import ogr
import numpy
import csv

import fileIO

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def formatTime(t):
    nd=int(t/86400.)
    nh=int((t-nd*86400.)/3600.)
    nm=int((t-nd*86400-nh*3600)/60.)
    ns=int(t-nd*86400-nh*3600-nm*60.)

    return "%02id:%02ih:%02im:%02is"%(nd,nh,nm,ns)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def processFlowLineBCs(flowLineShp,xll,yll,xsz,ysz,cellSize,flowAttr=None):

    if flowLineShp is None:
        return 0, None, None, None


    flowLine=fileIO.readPolylineShapefile(flowLineShp)

    # Generate list of shapely shapes corresponding to grid
    gridSquareList=[]
    for i in range(xsz):
        for j in range(ysz):
            x0=xll+i*cellSize
            x1=x0+cellSize
            y0=yll+j*cellSize
            y1=y0+cellSize

            wkt="POLYGON((%f %f,%f %f,%f %f,%f %f,%f %f))"%(x0,y0,x1,y0,x1,y1,x0,y1,x0,y0)

            gridSquareList.append((i,j,shapely.wkt.loads(wkt)))


    flowPointsN=0
    flowPoints=[]

    for xl,yl,attrDict in flowLine:
        wkt="LINESTRING("
        for xpt,ypt in zip(xl,yl):
            wkt+="%f %f,"%(xpt,ypt)
        wkt=wkt[:-1]+")" # Remove last comma and close brackets

#        print wkt
        line=shapely.wkt.loads(wkt)

        for i,j,gridSq in gridSquareList:
            if gridSq.intersects(line):
                qi=gridSq.intersection(line).length
                if flowAttr is not None:
                    qi*=attrDict[flowAttr]

                flowPointsN+=1
                flowPoints.append((i,j,qi))


    if flowPointsN>0:
        flowPointsXi=numpy.zeros(flowPointsN)
        flowPointsYi=numpy.zeros(flowPointsN)
        flowPointsMaxQ=numpy.zeros(flowPointsN)

        for i in range(flowPointsN):
            flowPointsXi[i]=flowPoints[i][0]
            flowPointsYi[i]=flowPoints[i][1]
            flowPointsMaxQ[i]=flowPoints[i][2]

    return flowPointsN,flowPointsXi,flowPointsYi,flowPointsMaxQ


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
