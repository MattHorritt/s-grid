# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 17:43:30 2018

@author: mattyh
"""

import fileIO
import math
import numpy
import tempfile
import os
import psycopg2
import time


# Connect to database
try:
    conn = psycopg2.connect(\
        "dbname='awash' port=5433 user='postgres' host='rabbit' password='pg_1langarr0n_Hr9'")
except:
    print "Unable to connect to the database"

conn.autocommit=True

cur = conn.cursor()


tileSize=100000 # 100km tiles

pondAreaThreshold=1e4      # 10 ha or approximately 100 cells
islandAreaThreshold=5e4    # 50 ha or approximately 100 cells

depthThreshold=0.1

rpList=[10,100,1000] #,2,5,10,20,50,100,200,500,1000,2000,5000,10000]

t1=time.time()
count=0

for rp in rpList:
    depthFileName="maxDepth6h_%iyr_cc30.tif"%rp
    outputExtentFileName="extent%iyr_cc30.tiff"%rp
    outputTableName="s%iyr_cc30"%rp

    # Load raster information
    f,dx,xsz,ysz,xmin,ymin=fileIO.readScalarGridObj(depthFileName)

    xmax=xmin+dx*xsz
    ymax=ymin+dx*ysz

    nTilesX=int(math.ceil(xsz*dx/tileSize))
    nTilesY=int(math.ceil(ysz*dx/tileSize))

    nTiles=nTilesX*nTilesY*len(rpList)

    # Create blank raster - this will have extents burned in later
    fileIO.saveScalarGridByte(numpy.zeros((xsz,ysz),dtype=numpy.byte),xmin,ymin,dx,outputExtentFileName)

    # Parameters for postgis
    pgStr='PG:"dbname=awash host=rabbit port=5433 user=postgres ACTIVE_SCHEMA=hazard_sw"'

    tableNameList=[]

    # Run through tiles
    for i in range(nTilesX):
        for j in range(nTilesY):

            print "Processing tile %i/%i %i/%i"%((i+1),nTilesX,(j+1),nTilesY),

            tileXmin=xmin+i*tileSize
            tileXmax=tileXmin+tileSize
            tileYmin=ymin+j*tileSize
            tileYmax=tileYmin+tileSize

            # Generate depth map tile
            gdalCmdStr='gdal_translate -ot Float32 -of GTiff '
            gdalCmdStr+='-projwin %f %f %f %f '%(tileXmin,tileYmax,tileXmax,tileYmin)
            gdalCmdStr+=depthFileName+' '

            tmpDepthFileName=tempfile._get_candidate_names().next()+'.tiff'

            gdalCmdStr+=tmpDepthFileName

#            print gdalCmdStr

#            assert False

#            print "CMD1:", gdalCmdStr
            os.system(gdalCmdStr+' > /dev/null 2>&1')

            # Threshold into extent
            tmpExtentFileName=tempfile._get_candidate_names().next()+'.tiff'

            gdalCmdStr='gdal_calc.py -A %s --calc="A>%f" --outfile=%s --type=Byte'\
                %(tmpDepthFileName,depthThreshold,tmpExtentFileName)
#            print "CMD2:", gdalCmdStr
            os.system(gdalCmdStr+' > /dev/null 2>&1')

            # Convert to vector and load to PostGIS table - use temporary table for now
            tableName='tmp_%02i%02i'%(i,j)
            tableNameList.append(tableName)
            cur.execute("drop table if exists hazard_sw.tmp")

            gdalCmdStr='gdal_polygonize.py %s -f PostgreSQL  %s tmp'%(tmpExtentFileName,pgStr)
#            print "CMD3:", gdalCmdStr
            os.system(gdalCmdStr+' > /dev/null 2>&1')

            os.remove(tmpExtentFileName)
            os.remove(tmpDepthFileName)

#            assert False

            # Remove dn=0 polygons
            cur.execute("delete from hazard_sw.tmp where dn=0 or dn is NULL")

            # Remove ponds
            cur.execute("delete from hazard_sw.tmp where st_area(wkb_geometry)<%s"%pondAreaThreshold)

            # Remove islands by writing outer rings to new table
            cur.execute('drop table if exists hazard_sw.%s'%tableName)

            sqlStr='create table hazard_sw.%s as '%tableName
            sqlStr+='WITH rings AS ( '
            sqlStr+='SELECT ogc_fid, (ST_DumpRings((st_dump(st_makevalid(wkb_geometry))).geom)).geom '
            sqlStr+='FROM hazard_sw.tmp) '
            sqlStr+='SELECT -ogc_fid as ogc_fid, ST_BuildArea(ST_Collect(geom)) as geom '
            sqlStr+='from rings WHERE ST_Area(geom) > %f GROUP BY ogc_fid'%islandAreaThreshold
            cur.execute(sqlStr)

            # burn into binary extent raster
            gdalCmdStr='gdal_rasterize -burn 1 -l %s %s %s'%(tableName,pgStr,outputExtentFileName)
            os.system(gdalCmdStr+' > /dev/null 2>&1')

            # Predict finish time
            count+=1
            pcComplete=float(count)/nTiles
            pcToGo=1.-pcComplete
            t2=time.time()-t1
            projectedFinish=time.time()+pcToGo*(t2/pcComplete)
#            projectedFinishString=time.strftime("%H:%M:%S",time.localtime(projectedFinish))
            projectedFinishString=time.strftime("%H:%M:%S %A %-d %B",time.localtime(projectedFinish))


            print 'Finish=',projectedFinishString

    # Merge outputs to single table
#    print "Merging to single table..."
    sqlStr='drop table if exists hazard_sw.%s;'%outputTableName
    sqlStr+='create table hazard_sw.%s as '%outputTableName
    sqlStr+='select * from hazard_sw.%s '%tableNameList[0]
    for tableName in tableNameList[1:]:
        sqlStr+='union select * from hazard_sw.%s '%tableName
    cur.execute(sqlStr)

    # And delete tmp files
    for tableName in tableNameList:
        sqlStr='drop table if exists hazard_sw.%s '%tableName
        cur.execute(sqlStr)
