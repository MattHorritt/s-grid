from osgeo import gdal
from osgeo import ogr
import numpy
import csv
import os
import psycopg2
import math
import tempfile

# Hack to avoid proj issue throwing exception
gdal.DontUseExceptions()
os.environ['PROJ_LIB'] = '/home/mattyh/miniconda3/envs/sgrid/share/proj'

class geoGrid:
    def __init__(self,*args,**kwargs):

#        print(args, len(args))

        if len(args)==6:
            self.grid=args[0]
            self.xsz=args[1]
            self.ysz=args[2]
            self.xll=args[3]
            self.yll=args[4]
            self.dx=args[5]

        elif len(args)==1:
            f=gdal.Open(args[0])

            self.fileName = args[0]

            if "objOnly" in kwargs and kwargs["objOnly"]:
                self.obj = f
            else:
                self.grid=f.ReadAsArray().transpose().copy()

                if "dataType" in kwargs:
                    self.grid=numpy.array(self.grid[:,::-1],kwargs["dataType"])
                else:
                    self.grid=numpy.array(self.grid[:,::-1])

            geo=f.GetGeoTransform()

            self.xsz = f.RasterXSize
            self.ysz = f.RasterYSize

            self.dx = geo[1]
            self.xll = geo[0]
            self.yll = geo[3] - self.dx * self.ysz

            self.noDataValue = f.GetRasterBand(1).GetNoDataValue()

    def subArea(self,xmin,xmax,ymin,ymax):
        new_grid=self.copy()

        ximin=int((xmin-self.xll)/self.dx)
        ximax=int((xmax-self.xll)/self.dx)
        yimin=int((ymin-self.yll)/self.dx)
        yimax=int((ymax-self.yll)/self.dx)

        new_grid.xsz=ximax-ximin
        new_grid.ysz=yimax-yimin
        new_grid.xll=xmin
        new_grid.yll=ymin

        new_grid.grid=new_grid.grid[ximin:ximax,yimin:yimax]

        return new_grid

    def getPointValueXy(self,x,y):
        i=int((x-self.xll)/self.dx)
        j=int((y-self.yll)/self.dx)

        if i>=0 and i<self.xsz and j>=0 and j<self.ysz:
            return self.grid[i,j]
        else:
            return None

    def save(self,fileName):
        geotiffDriver=gdal.GetDriverByName("GTiff")

        outputTIFF=geotiffDriver.Create(fileName,self.xsz,self.ysz,\
            1,gdal.GDT_Float32,options=['COMPRESS=LZW'])
        outputTIFF.SetGeoTransform([self.xll,self.dx,0,self.yll+self.ysz*self.dx,0,-self.dx])
        outputTIFF.GetRasterBand(1).WriteArray(self.grid.transpose()[::-1,:])
        outputTIFF.FlushCache()
        geotiffDriver=None

    def copy(self):
        cp=geoGrid(None,self.xsz,self.ysz,self.xll,self.yll,self.dx)
        cp.grid=self.grid.copy()

        return cp




def readPolylineShapefile(fileName,layerName=None):
    dataSource=ogr.Open(fileName)
    if dataSource is None:
        return None

    if layerName is None:
        layer=dataSource.GetLayerByIndex(0)
    else:
        layer=dataSource.GetLayer(layerName)


    layer.ResetReading()

    layerDefn=layer.GetLayerDefn()
    fieldNames=[layerDefn.GetFieldDefn(i).GetName() for i in range(layerDefn.GetFieldCount())]

    rds=[]

    layer.ResetReading()
    feature=layer.GetNextFeature()

    while feature is not None:
        xl=[]
        yl=[]
        geom=feature.GetGeometryRef()

        for gr in range(geom.GetGeometryCount()):
            g=geom.GetGeometryRef(gr)

            for i in range(g.GetPointCount()):
                pt=g.GetPoint(i)
                xl.append(pt[0])
                yl.append(pt[1])

        tmpDict={}
        for fn in fieldNames:
            val=feature.GetField(fn)
            tmpDict[fn]=val

        rds.append((xl,yl,tmpDict))

        feature=layer.GetNextFeature()

    return rds


def readPointShapefile(fileName):
    dataSource=ogr.Open(fileName)
    if dataSource is None:
        return None

    pointsLayer=dataSource.GetLayerByIndex(0)
    pointsLayer.ResetReading()

    layerDefn=pointsLayer.GetLayerDefn()
    fieldNames=[layerDefn.GetFieldDefn(i).GetName() for i in range(layerDefn.GetFieldCount())]

    rds=[]

    # featureCount = pointsLayer.GetFeatureCount()
    # pointsLayer.ResetReading()

    for pointFeature in pointsLayer:
        # pointFeature = pointsLayer.GetFeature(fi)

        X=pointFeature.geometry().GetPoint()[0]
        Y=pointFeature.geometry().GetPoint()[1]

        tmpDict={}
        for fn in fieldNames:
            val=pointFeature.GetField(fn)
            tmpDict[fn]=val

        rds.append((X,Y,tmpDict))

        # pointFeature=pointsLayer.GetNextFeature()

    return rds



def readScalarGrid(fileName,dataType=numpy.float64):
    f=gdal.Open(fileName)
    arr=f.ReadAsArray().transpose().copy()

    arr=numpy.array(arr[:,::-1],dataType)

    geo=f.GetGeoTransform()

    xsz,ysz=arr.shape
    dx=geo[1]
    xll=geo[0]
    yll=geo[3]-dx*ysz

    return arr,dx,xsz,ysz,xll,yll

def readScalarGridObj(fileName,dataType=numpy.float64):
    f=gdal.Open(fileName)

    geo=f.GetGeoTransform()

    xsz=f.RasterXSize
    ysz=f.RasterYSize

    dx=geo[1]
    xll=geo[0]
    yll=geo[3]-dx*ysz

    return f,dx,xsz,ysz,xll,yll

def saveDiagCSV(u45,u135,xll,yll,dx,fileName):
    xsz=u45.shape[0]
    ysz=u45.shape[1]

    f=open(fileName,"w")
    f.write("X,Y,Val,Dir\n")

    for i in range(xsz):
        for j in range(ysz):
            if u45[i,j]>0:
                dir=45
            else:
                dir=225
            f.write("%i,%i,%f,%f,%f,%f\n"%(i,j,xll+dx*i,yll+dx*j,abs(u45[i,j]),dir))

    for i in range(xsz):
        for j in range(ysz):
            if u135[i,j]>0:
                dir=-45
            else:
                dir=135
            f.write("%i,%i,%f,%f,%f,%f\n"%(i,j,xll+dx*i,yll+dx*j,abs(u135[i,j]),dir))

    f.close()


def saveConveyanceParametersCSV(parX,parY,xll,yll,dx,fileName):
    xsz=parY.shape[0]
    ysz=parX.shape[1]

    f=open(fileName,"w")
    f.write("I,J,X,Y,Zmin,Zbank,slope1,int1,slope2,int2\n")

    # X directions
    for i in range(xsz+1):
        for j in range(ysz):
            x=xll+i*dx
            y=yll+j*dx+dx/2

#            x=xll+i*dx+dx/2
#            y=yll+j*dx

            if parX[i, j, 0] != -9999:
                f.write("%i,%i,%f,%f,%f,%f,%f,%f,%f,%f\n"%\
                    (i,j,x,y,parX[i,j,0],parX[i,j,1],parX[i,j,2],parX[i,j,3],parX[i,j,4],parX[i,j,5]))

    # Y directions
    for i in range(xsz):
        for j in range(ysz+1):
            x=xll+i*dx+dx/2
            y=yll+j*dx

#            x=xll+i*dx
#            y=yll+j*dx+dx/2
            if parY[i, j, 0] != -9999:
                f.write("%i,%i,%f,%f,%f,%f,%f,%f,%f,%f\n"%\
                    (i,j,x,y,parY[i,j,0],parY[i,j,1],parY[i,j,2],parY[i,j,3],parY[i,j,4],parY[i,j,5]))

    f.close()
    return

def saveStorageParametersCSV(par,xll,yll,dx,fileName):
    xsz=par.shape[0]
    ysz=par.shape[1]

    f=open(fileName,"w")

    if par.shape[-1]==7:
        f.write("I,J,X,Y,zMin,zBank,zMax,vip1,vip2,vip3,vip4\n")

        for i in range(xsz):
            for j in range(ysz):
                x=xll+i*dx+dx/2
                y=yll+j*dx+dx/2

                if par[i,j,0] != -9999:
                    f.write("%i,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f\n"\
                        %(i,j,x,y,par[i,j,0],par[i,j,1],par[i,j,2],par[i,j,3],par[i,j,4],\
                            par[i,j,5],par[i,j,6]))
    else:
        f.write("I,J,X,Y,zMin,zMax,vip1,vip2,vip3\n")

        for i in range(xsz):
            for j in range(ysz):
                x=xll+i*dx+dx/2
                y=yll+j*dx+dx/2

                if par[i,j,0] != -9999:
                    f.write("%i,%i,%f,%f,%f,%f,%f,%f,%f\n"\
                        %(i,j,x,y,par[i,j,0],par[i,j,1],par[i,j,2],par[i,j,3],par[i,j,4]))


    f.close()
    return

def saveVectorCSV(vx,vy,xll,yll,dx,fileName, thresholdVal=None):
    xsz=vx.shape[0]
    ysz=vy.shape[1]

    f=open(fileName,"w")

    f.write("i,j,X,Y,Val,Dir\n")


    if thresholdVal is None:
        thresholdVal=-100.

    for i in range(xsz):
        for j in range(ysz-1):
            if abs(vx[i,j])>thresholdVal:
                if vx[i,j]>0:
                    dir=90
                else:
                    dir=270
                f.write("%i,%i,%f,%f,%f,%f\n"%(i,j,xll+dx*i,yll+dx*j+dx/2,abs(vx[i,j]),dir))

    for i in range(xsz-1):
        for j in range(ysz):
            if abs(vy[i,j])>thresholdVal:
                if vy[i,j]>0:
                    dir=0
                else:
                    dir=180
                f.write("%i,%i,%f,%f,%f,%f\n"%(i,j,xll+dx*i+dx/2,yll+dx*j,abs(vy[i,j]),dir))

    f.close()


def readFlowCsv(fileName,xsz,ysz,dataType=numpy.float64):

    retArrX=numpy.zeros((xsz+1,ysz),dtype=dataType)
    retArrY=numpy.zeros((xsz,ysz+1),dtype=dataType)

    csvFile=open(fileName)

    csvReader=csv.reader(csvFile)

    headerRow=next(csvReader)
    magIdx=headerRow.index("Val")
    dirIdx=headerRow.index("Dir")

    for row in csvReader:
        i=int(row[0])
        j=int(row[1])

        mag=float(row[magIdx])
        dir=float(row[dirIdx])

        if dir==0.0: # North
            retArrY[i,j]=mag
        elif dir==180.0: # South
            retArrY[i,j]=-mag
        if dir==90.0: # East
            retArrX[i,j]=mag
        if dir==270.0: # West
            retArrX[i,j]=-mag

    return retArrX, retArrY


def readCSV(fileName,columnHeader,xsz,ysz,dataType=numpy.float64):
    retArr=numpy.zeros((xsz,ysz),dtype=dataType)

    csvFile=open(fileName)

    csvReader=csv.reader(csvFile)

    headerRow=next(csvReader)
    dataIdx=headerRow.index(columnHeader)

    for row in csvReader:
        i=int(row[2])
        j=int(row[3])

        val=float(row[dataIdx])

        retArr[i,j]=val

    return retArr



def saveScalarCSV(sArg,xll,yll,dx,fileName,headerList=None):

    if type(sArg)==numpy.ndarray:
        s=[sArg]
        numGrids=1
    else:
        numGrids=len(sArg)
        s=sArg

    xsz=s[0].shape[0]
    ysz=s[0].shape[1]

    f=open(fileName,"w")
    f.write("X,Y,i,j")
    if headerList is None:
        for i in range(numGrids):
            f.write(",Val%i"%i)
        f.write("\n")
    else:
        for header in headerList:
            f.write(",%s"%header)
        f.write("\n")

    for i in range(xsz):
        for j in range(ysz):
            f.write("%f,%f,%i,%i"%(xll+dx*i+dx/2,yll+dx*j+dx/2,i,j))
            for k in range(numGrids):
                f.write(",%f"%s[k][i,j])
            f.write("\n")

    f.close()


def saveScalarGrid(s,xll,yll,dx,fileName):
    xsz=s.shape[0]
    ysz=s.shape[1]


    geotiffDriver=gdal.GetDriverByName("GTiff")

    outputTIFF=geotiffDriver.Create(fileName,xsz,ysz,1,gdal.GDT_Float32,options=['COMPRESS=LZW'])
    outputTIFF.SetGeoTransform([xll,dx,0,yll+ysz*dx,0,-dx])
#    outputTIFF.SetGeoTransform([xll+dx/2,dx,0,yll+ysz*dx+dx/2,0,-dx])
    outputTIFF.GetRasterBand(1).WriteArray(s.transpose()[::-1,:])
    outputTIFF.FlushCache()
    geotiffDriver=None

################################################################################
def uploadGridToDb(gridFileName, outputTableName, depthThreshold = None, dropTable = False):
    pgStr='PG:"dbname=ltis2025 host=localhost password=''postgres'' port=5432 user=postgres ACTIVE_SCHEMA=slr"'

    try:
        conn = psycopg2.connect("dbname='ltis2025' port=5432 user='postgres' host='localhost' password='postgres'")
    except:
        print
        "Unable to connect to the database"

    conn.autocommit = True

    cur = conn.cursor()

    gGrid = geoGrid(gridFileName, objOnly = True)

    tileSize = 100000
    pondAreaThreshold = 1e4  # 10 ha
    islandAreaThreshold = 5e4  # 50 ha

    nTilesX=int(math.ceil(gGrid.xsz*gGrid.dx/tileSize))
    nTilesY=int(math.ceil(gGrid.ysz*gGrid.dx/tileSize))

    tableNameList=[]

    if depthThreshold is None:
        depthThreshold = 0

    # Run through tiles
    for i in range(nTilesX):
        for j in range(nTilesY):

            print(f"Processing tile {i+1}/{nTilesX} {j+1}/{nTilesY}... ", end='')

            tileXmin=gGrid.xll + i * tileSize
            tileXmax=tileXmin + tileSize
            tileYmin=gGrid.yll + j * tileSize
            tileYmax=tileYmin + tileSize

            # Generate depth map tile
            gdalCmdStr = 'gdal_translate -ot Float32 -of GTiff '
            gdalCmdStr += f'-projwin {tileXmin} {tileYmax} {tileXmax} {tileYmin} '
            gdalCmdStr += gridFileName + ' '

            tmpDepthFileName=next(tempfile._get_candidate_names())+'.tiff'

            gdalCmdStr+=tmpDepthFileName
            os.system(gdalCmdStr + ' > /dev/null 2>&1')

            # Threshold into extent
            tmpExtentFileName=next(tempfile._get_candidate_names())+'.tiff'

            gdalCmdStr='gdal_calc.py -A %s --calc="A>%f" --outfile=%s --type=Byte'\
                %(tmpDepthFileName,depthThreshold,tmpExtentFileName)
            os.system(gdalCmdStr + ' > /dev/null 2>&1')

            # Convert to vector and load to PostGIS table - use temporary table for now
            tableName='tmp_%02i%02i'%(i,j)
            tableNameList.append(tableName)
            cur.execute("drop table if exists slr.tmp")

            gdalCmdStr='gdal_polygonize.py %s -f PostgreSQL  %s tmp'%(tmpExtentFileName,pgStr)
            os.system(gdalCmdStr + ' > /dev/null 2>&1')

            os.remove(tmpExtentFileName)
            os.remove(tmpDepthFileName)

            # Remove dn=0 polygons
            cur.execute("delete from slr.tmp where dn=0 or dn is NULL")

            # Remove ponds
            cur.execute("delete from slr.tmp where st_area(wkb_geometry)<%s"%pondAreaThreshold)

            # Remove islands by writing outer rings to new table
            cur.execute('drop table if exists slr.%s'%tableName)

            sqlStr='create table slr.%s as '%tableName
            sqlStr+='WITH rings AS ( '
            sqlStr+='SELECT ogc_fid, (ST_DumpRings((st_dump(st_makevalid(wkb_geometry))).geom)).geom '
            sqlStr+='FROM slr.tmp) '
            sqlStr+='SELECT -ogc_fid as ogc_fid, ST_BuildArea(ST_Collect(geom)) as geom '
            sqlStr+='from rings WHERE ST_Area(geom) > %f GROUP BY ogc_fid'%islandAreaThreshold
            cur.execute(sqlStr)

    # Merge outputs to single table
    if dropTable:
        sqlStr='drop table if exists slr.%s;'%outputTableName
        sqlStr+='create table slr.%s as '%outputTableName

    sqlStr+='select * from slr.%s '%tableNameList[0]
    for tableName in tableNameList[1:]:
        sqlStr+='union select * from slr.%s '%tableName
    cur.execute(sqlStr)

    # Add depth values
    sqlStr = 'alter table slr.%s add column depth real' % tableName
    cur.execute(sqlStr)

    sqlStr = 'update slr.%s set depth = %s' % (tableName, depthThreshold)
    cur.execute(sqlStr)

    # And delete tmp files
    for tableName in tableNameList:
        sqlStr='drop table if exists slr.%s '%tableName
        cur.execute(sqlStr)

    return

################################################################################
def uploadGridToDb2(depthGridFileName, outputTableName, depthThreshold = None, dropTable = False):
    pgStr='PG:"dbname=ltis2025 host=localhost password=''postgres'' port=5432 user=postgres ACTIVE_SCHEMA=slr"'

    try:
        conn = psycopg2.connect("dbname='ltis2025' port=5432 user='postgres' host='localhost' password='postgres'")
    except:
        print
        "Unable to connect to the database"

    conn.autocommit = True

    cur = conn.cursor()

    # gGrid = geoGrid(gridFileName, objOnly = True)

    if depthThreshold is None:
        depthThreshold = 0


    print(f"Processing depth grid {depthGridFileName}... ", end='')
    tmpName = next(tempfile._get_candidate_names())


    # Threshold into extent
    tmpExtentFileName = tmpName + '.tiff'

    gdalCmdStr='gdal_calc.py -A %s --calc="A>%f" --outfile=%s --type=Byte'\
        %(depthGridFileName,depthThreshold,tmpExtentFileName)
    os.system(gdalCmdStr) # + ' > /dev/null 2>&1')

    # Convert to vector and load to PostGIS table - use temporary table for now
    cur.execute("drop table if exists slr.%s"%(tmpName))

    gdalCmdStr='gdal_polygonize.py %s -f PostgreSQL  %s %s'%(tmpExtentFileName,pgStr,tmpName)
    os.system(gdalCmdStr) # + ' > /dev/null 2>&1')

    os.remove(tmpExtentFileName)

    # Remove dn=0 polygons
    cur.execute("delete from slr.%s where dn=0 or dn is NULL"%(tmpName))

    # Merge outputs to output table
    if dropTable:
        sqlStr='drop table if exists slr.%s;'%outputTableName
        sqlStr+='create table slr.%s (depth real, geom geometry) '%outputTableName

        cur.execute(sqlStr)

    # Copy rows from slr.tmp into output table
    sqlStr = 'insert into slr.%s select %s, wkb_geometry from slr.%s' % (outputTableName, depthThreshold, tmpName)
    cur.execute(sqlStr)

    # And delete tmp files
    sqlStr='drop table if exists slr.%s '%(tmpName)
    cur.execute(sqlStr)

    return

def saveScalarGridByte(s,xll,yll,dx,fileName):
    xsz=s.shape[0]
    ysz=s.shape[1]


    geotiffDriver=gdal.GetDriverByName("GTiff")

    outputTIFF=geotiffDriver.Create(fileName,xsz,ysz,1,gdal.GDT_Byte,options=['COMPRESS=LZW'])
    outputTIFF.SetGeoTransform([xll,dx,0,yll+ysz*dx,0,-dx])
#    outputTIFF.SetGeoTransform([xll+dx/2,dx,0,yll+ysz*dx+dx/2,0,-dx])
    outputTIFF.GetRasterBand(1).WriteArray(s.transpose()[::-1,:])
    outputTIFF.FlushCache()
    geotiffDriver=None

if __name__ == "__main__":

    # Test upload to DB
    gridFileName = r"2m_EA_results/output_slr1p0_rp200_max_depth.tif"

    for d in [0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
        print("Processing threshold ", d)
        depthLabel = f"{d:.1f}"
        depthLabel = depthLabel.replace('.', 'p')
        uploadGridToDb2(gridFileName, f"slr1p0_rp200_d{depthLabel}", depthThreshold = d, dropTable=True)
