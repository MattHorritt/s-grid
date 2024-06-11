# -*- coding: utf-8 -*-

import shapely
import shapely.wkt
from osgeo import ogr
import numpy
from math import sqrt
import fileIO
from subprocess import call

###########################################################################
def fillEdgeCells(cellDict,wl):
    for i,j in cellDict.keys():
        wl[i,j]=wl[cellDict[(i,j)][0],cellDict[(i,j)][1]]

    return

###########################################################################
def maskZeroPoly(wktFileName,rasterFileName):

    print(wktFileName,rasterFileName)

    gdalCmd=['gdal_rasterize']
    gdalCmd+=['-burn','0']
#    gdalCmd+=['-l',tmpWktFileName.split('.')[0]]
    gdalCmd+=[wktFileName]
    gdalCmd+=[rasterFileName]

    call(gdalCmd)

    return

###########################################################################
def weirFlowCalc(wl1,wl2,cl,len,wc,ml):

    h1=wl1-cl
    h2=wl2-cl

    if h1<0 and h2<0: # Zero flow
        return 0

    if wl1>wl2 and (h2/h1)<=ml: # Forward free flow
        return wc*len*(h1**1.5)

    if wl1>wl2 and (h2/h1)>ml: # Forward drowned flow
        return wc*len*h1*(((h1-h2)/(1-ml))**0.5)

    if wl1<wl2 and (h1/h2)<=ml: # Reverse free flow
        return -wc*len*(h2**1.5)

    if wl1<wl2 and (h1/h2)>ml: # Reverse drowned flow
        return -wc*len*h2*(((h2-h1)/(1-ml))**0.5)

    return 0.0

def calcSpills(wlGrid,spills,qX,qY):
    modularLimit=0.9

    for (i,j,dir) in spills.keys():

        crestLevel=spills[(i,j,dir)][0]
        crestLength=spills[(i,j,dir)][1]
        weirCoeff=spills[(i,j,dir)][2]

        if dir=='x':
            wl1=wlGrid[i-1,j]
            wl2=wlGrid[i,j]

            qX[i,j]+=weirFlowCalc(wl1,wl2,crestLevel,crestLength,weirCoeff,modularLimit)

            # Add spillway flow if present
            if len(spills[(i,j,dir)])>3:
                crestLevel=spills[(i,j,dir)][3]
                crestLength=spills[(i,j,dir)][4]
                weirCoeff=spills[(i,j,dir)][5]

                qX[i,j]+=weirFlowCalc(wl1,wl2,crestLevel,crestLength,weirCoeff,modularLimit)

        if dir=='y':
            wl1=wlGrid[i,j]
            wl2=wlGrid[i,j-1]

            qY[i,j]+=weirFlowCalc(wl1,wl2,crestLevel,crestLength,weirCoeff,modularLimit)

            # Add spillway flow if present
            if len(spills[(i,j,dir)])>3:
                crestLevel=spills[(i,j,dir)][3]
                crestLength=spills[(i,j,dir)][4]
                weirCoeff=spills[(i,j,dir)][5]

                qY[i,j]+=weirFlowCalc(wl1,wl2,crestLevel,crestLength,weirCoeff,modularLimit)

    return

###########################################################################

###########################################################################
def spillways(xll,yll,cellSize,xsz,ysz,spills,shpFileName):
    pts=fileIO.readPointShapefile(shpFileName)

    for x, y, attr in pts:
        i=int((x-xll)/cellSize)
        j=int((y-yll)/cellSize)

        x0=i*cellSize+xll
        y0=j*cellSize+yll

        # How many edges do we have?
        edges=[]
        if (i,j,'x') in spills.keys():
            edges.append('w')
        if (i+1,j,'x') in spills.keys():
            edges.append('e')
        if (i,j,'y') in spills.keys():
            edges.append('s')
        if (i,j+1,'y') in spills.keys():
            edges.append('n')

        # If none - something wrong
        if len(edges)==0:
            print("***** Warning - spillway at %f %f cannot be assigned to edge - ignoring ***"%(x,y))
            continue

        # If multiple edges, use nearest
        if len(edges)>1:
            minDist=1e20
            if 'n' in edges:
                dist=y0+cellSize-y
                if dist<minDist:
                    minDist=dist
                    edge='n'
            if 's' in edges:
                dist=y-y0
                if dist<minDist:
                    minDist=dist
                    edge='s'
            if 'e' in edges:
                dist=x0+cellSize-x
                if dist<minDist:
                    minDist=dist
                    edge='e'
            if 'w' in edges:
                dist=x-x0
                if dist<minDist:
                    minDist=dist
                    edge='w'
        else:
            edge=edges[0]

        if edge=='n':
            spills[(i,j+1,'y')]+=[attr['cl'],attr['width'],attr['wc']]
            spills[(i,j+1,'y')][1]-=attr['width']
        if edge=='s':
            spills[(i,j,'y')]+=[attr['cl'],attr['width'],attr['wc']]
            spills[(i,j,'y')][1]-=attr['width']
        if edge=='e':
            spills[(i+1,j,'x')]+=[attr['cl'],attr['width'],attr['wc']]
            spills[(i+1,j,'x')][1]-=attr['width']
        if edge=='w':
            spills[(i,j,'x')]+=[attr['cl'],attr['width'],attr['wc']]
            spills[(i,j,'x')][1]-=attr['width']

    return


###########################################################################
def overlaps(xmin1,xmax1,ymin1,ymax1,xmin2,xmax2,ymin2,ymax2):
    return (xmin1<xmax2 and ymin1<ymax2 and xmax1>xmin2 and ymax1>ymin2)

def almostSame(x1,x2):
    return abs((x1-x2))/(0.5*(x1+x2))<1e-6

def edgeSame(e1,e2):
    return e1==e2 or e1==[e2[2],e2[3],e2[0],e2[1]]

def addEdgeArray(arrX,arrY,clX,clY,edge,xll,yll,dx):
    x1=edge[0]
    y1=edge[1]
    x2=edge[2]
    y2=edge[3]

    if x1==x2: # X Edge
        xi=int((x1-xll)/dx)
        yi=int((0.5*(y1+y2)-yll)/dx)

        arrX[xi,yi]=True

        clX[xi,yi]=dx

    if y1==y2: # Y Edge
        xi=int((0.5*(x1+x2)-xll)/dx)
        yi=int((y1-yll)/dx)

        arrY[xi,yi]=True
        clY[xi,yi]=dx

    return

def rmEdgeArray(arrX,arrY,edge,xll,yll,dx):
    x1=edge[0]
    y1=edge[1]
    x2=edge[2]
    y2=edge[3]

    if x1==x2: # X Edge
        xi=int((x1-xll)/dx)
        yi=int((0.5*(y1+y2)-yll)/dx)

        arrX[xi,yi]=False

    if y1==y2: # Y Edge
        xi=int((0.5*(x1+x2)-xll)/dx)
        yi=int((y1-yll)/dx)

        arrY[xi,yi]=False

    return

def vectorProduct(x1,y1,x2,y2):
    return x1*y2-y1*x2

def modifyConvPar(xll,yll,cellSize,xsz,ysz,convParX,convParY,shpFileName,csvEdgeFileName=None):

    shpSource = ogr.Open(shpFileName,0)

    # Generate list of crest geometries, crest levels and weir coefficients
    lyr=shpSource.GetLayer()
    polylineGeoms=[shapely.wkt.loads(f.GetGeometryRef().ExportToWkt()) for f in lyr]

    lyr.ResetReading()
    crestLevels=[f.GetField("cl") for f in lyr]

    lyr.ResetReading()
    weirCoeffs=[f.GetField("wc") for f in lyr]

#    combEdgeList=[]
    edgeArrX=numpy.zeros((xsz+1,ysz),dtype=bool)
    edgeArrY=numpy.zeros((xsz,ysz+1),dtype=bool)

    crestLengthX={}
    crestLengthY={}
    crestSpillList={}

    totalPolylineLength=0

    # Lists of polygons to be filled or zeroed later to tidy up reservoir edge
    fillCellDict={}
    zeroPolyList=[]

    for pi, pline in enumerate(polylineGeoms):

        totalPolylineLength+=pline.length

        xList=[]
        yList=[]

        for point in pline.coords:
            xList.append(point[0])
            yList.append(point[1])

        i1=int((xList[0]-xll)/cellSize)
        j1=int((yList[0]-yll)/cellSize)

        first=True
        last=False
        pointList=[(xList[0],yList[0])]
        edgeList=[]
#        ijEdgeList=[]

        count=0

        # The logic in this loop is mangled - need to refactor
        while True:
            count+=1
            if count>100:
                break

            x0=xll+i1*cellSize
            x1=x0+cellSize
            y0=yll+j1*cellSize
            y1=y0+cellSize

            xSqCentre=0.5*(x0+x1)
            ySqCentre=0.5*(y0+y1)

            wkt="LINESTRING(%f %f,%f %f,%f %f,%f %f,%f %f)"%(x0,y0,x1,y0,x1,y1,x0,y1,x0,y0)
            gridSquare=shapely.wkt.loads(wkt)


            intersectionPoints=gridSquare.intersection(pline)
            if intersectionPoints.type=='MultiPoint':
                nextPoints=[]
                for i in range(len(intersectionPoints)):
                    nextPoints.append(intersectionPoints[i].coords[0])

                for nx,ny in nextPoints:
                    if (nx,ny) not in pointList:
                        pointList.append((nx,ny))

            else:
#                nextPoints=[c for c in intersectionPoints.coords]
                nextPoints=intersectionPoints.coords[0]
                if first:
                    pointList.append(nextPoints)
                else:
                    pointList.append((xList[-1],yList[-1]))
                    last=True


            nextPoints=nextPoints[:2] # Only keep first two points (e.g. if shape is really wiggly)


            # Calculate crest length
            wkt="POLYGON((%f %f,%f %f,%f %f,%f %f,%f %f))"%(x0,y0,x1,y0,x1,y1,x0,y1,x0,y0)
            gridSquare=shapely.wkt.loads(wkt)
            clippedLine=gridSquare.intersection(pline)
            crestLength=clippedLine.length

            # Calculate fill/zero polygons
            clippedLineBuff=clippedLine.buffer(1.0) # Buffer line by 1m
            mp=gridSquare.difference(clippedLineBuff) # Should be multipolygon with 2 parts

            if mp.type=='MultiPolygon' and len(mp)==2:
                c1x,c1y=mp[0].centroid.coords[0]
                c2x,c2y=mp[1].centroid.coords[0]

                p1x,p1y=clippedLine.coords[0]
                p2x,p2y=clippedLine.coords[-1]

                vp=vectorProduct(p2x-p1x,p2y-p1y,c2x-c1x,c2y-c1y)

                if vp>0:
#                    fillPolyList.append(mp[0])
                    zeroPolyList.append(mp[1])
                else:
#                    fillPolyList.append(mp[1])
                    zeroPolyList.append(mp[0])

#                fillCellList.append((i1,j1))


            # Find midpoint of line - use this to work out where the edge goes
            xm=0.5*(pointList[-1][0]+pointList[-2][0])
            ym=0.5*(pointList[-1][1]+pointList[-2][1])


            if almostSame(pointList[-1][1],y0): # South

                # Eastern edge
                if not first and xm>xSqCentre:
                    edgeList.append([x1,y1,x1,y0])
                    edgeArrX[i1+1,j1]=True
                    crestLengthX[(i1+1,j1)]=crestLength

                # Western edge
                if not first and xm<=xSqCentre:
                    edgeList.append([x0,y1,x0,y0])
                    edgeArrX[i1,j1]=True
                    crestLengthX[(i1,j1)]=crestLength

                    fillCellDict[(i1,j1)]=(i1-1,j1)


                j1-=1


            if almostSame(pointList[-1][1],y1): # North

                # Eastern edge
                if not first and xm>xSqCentre:
                    edgeList.append([x1,y0,x1,y1])
                    edgeArrX[i1+1,j1]=True
                    crestLengthX[(i1+1,j1)]=crestLength

                    fillCellDict[(i1,j1)]=(i1+1,j1)

                # Western edge
                if not first and xm<=xSqCentre:
                    edgeList.append([x0,y0,x0,y1])
                    edgeArrX[i1,j1]=True
                    crestLengthX[(i1,j1)]=crestLength

                j1+=1


            if almostSame(pointList[-1][0],x0): # West

                # Northern edge
                if not first and ym>ySqCentre:
                    edgeList.append([x1,y1,x0,y1])
                    edgeArrY[i1,j1+1]=True
                    crestLengthY[(i1,j1+1)]=crestLength

                    fillCellDict[(i1,j1)]=(i1,j1+1)

                # Southern edge
                if not first and ym<=ySqCentre:
                    edgeList.append([x1,y0,x0,y0])
                    edgeArrY[i1,j1]=True
                    crestLengthY[(i1,j1)]=crestLength

                i1-=1

            if almostSame(pointList[-1][0],x1): # East

                # Northern edge
                if not first and ym>ySqCentre:
                    edgeList.append([x0,y1,x1,y1])
                    edgeArrY[i1,j1+1]=True
                    crestLengthY[(i1,j1+1)]=crestLength
                # Southern edge
                if not first and ym<=ySqCentre:
                    edgeList.append([x0,y0,x1,y0])
                    edgeArrY[i1,j1]=True
                    crestLengthY[(i1,j1)]=crestLength

                    fillCellDict[(i1,j1)]=(i1,j1-1)

                i1+=1


            first=False
            if last:
                break

        # Link gaps in edge list
        for i, edge in enumerate(edgeList[:-1]):
            if edge[2]!=edgeList[i+1][0] or edge[3]!=edgeList[i+1][1]:
                edgeList.append([edge[2],edge[3],edgeList[i+1][0],edgeList[i+1][1]])
                addEdgeArray(edgeArrX,edgeArrY,crestLengthX,crestLengthY,edgeList[-1],xll,yll,cellSize)

        # And remove duplicated edges - these are always wrong, so remove both
        for i, edge1 in enumerate(edgeList[:]):
            for edge2 in edgeList[i+1:]:
                if edgeSame(edge1,edge2):
                    rmEdgeArray(edgeArrX,edgeArrY,edge1,xll,yll,cellSize)
                    edgeList.remove(edge1)
                    edgeList.remove(edge2)

    # Block convPar arrays and add to crest spill list
    totalCrestLength=0
    for i in range(xsz):
        for j in range(ysz):
            if edgeArrX[i,j]:
                convParX[i,j,:]=-9999
                crestSpillList[(i,j,'x')]=[crestLevels[pi],crestLengthX[(i,j)],weirCoeffs[pi]]
                totalCrestLength+=crestLengthX[(i,j)]
            if edgeArrY[i,j]:
                convParY[i,j,:]=-9999
                crestSpillList[(i,j,'y')]=[crestLevels[pi],crestLengthY[(i,j)],weirCoeffs[pi]]
                totalCrestLength+=crestLengthY[(i,j)]

    # Save edges as diagnostic
    if csvEdgeFileName is not None:
        wktFile=open(csvEdgeFileName,"w")
        wktFile.write("id;wkt;level;length;wc\n")
        count=0
        for i in range(xsz):
            for j in range(ysz):
                if edgeArrX[i,j]:
                    x1=i*cellSize+xll
                    y1=j*cellSize+yll
                    x2=x1
                    y2=y1+cellSize

                    wktFile.write("%i;LINESTRING(%f %f, %f %f);%f;%f;%f\n"%\
                        (count,x1,y1,x2,y2,\
                        crestSpillList[(i,j,'x')][0],crestSpillList[(i,j,'x')][1],crestSpillList[(i,j,'x')][2]))
                    count+=1

                if edgeArrY[i,j]:
                    x1=i*cellSize+xll
                    y1=j*cellSize+yll
                    x2=x1+cellSize
                    y2=y1

                    wktFile.write("%i;LINESTRING(%f %f, %f %f);%f;%f;%f\n"%\
                        (count,x1,y1,x2,y2,\
                        crestSpillList[(i,j,'y')][0],crestSpillList[(i,j,'y')][1],crestSpillList[(i,j,'y')][2]))
                    count+=1
        wktFile.close()

###############################################################################

#    csvFile=open("fillCellList.csv",'w')
#    csvFile.write("wkt\n")
#    for i,j in fillCellDict.keys():
#
#        wkt='linestring(%f %f, %f %f)'%(xll+fillCellDict[(i,j)][0]*cellSize+cellSize/2,
#            yll+fillCellDict[(i,j)][1]*cellSize+cellSize/2,
#            xll+i*cellSize+cellSize/2,yll+j*cellSize+cellSize/2)
#
#        csvFile.write(wkt+"\n")
#    csvFile.close()
#    assert False


###############################################################################
#    wktFile=open("fillZeroPoly.csv","w")
#    wktFile.write("id;wkt;type\n")
#    count=0
#
#    for poly in fillPolyList:
#        wktFile.write("%i;%s;fill\n"%(count,poly.wkt))
#        count+=1
#
#    for poly in zeroPolyList:
#        wktFile.write("%i;%s;zero\n"%(count,poly.wkt))
#        count+=1
#
#    wktFile.close()
#    assert False
###############################################################################

    print("Total polyline length=",totalPolylineLength)
    print("Total crest length added=",totalCrestLength)

    return crestSpillList, fillCellDict, zeroPolyList