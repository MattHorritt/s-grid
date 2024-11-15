#!/usr/bin/env python

import sys, os
import fileIO
from multiprocessing.dummy import Pool as ThreadPool

def wrapper(argTuple):
    fileIO.uploadGridToDb2(argTuple[0], argTuple[1], depthThreshold=argTuple[2], dropTable=True)

#############################################################################
# Control Panel

resultsDirectory="2m_EA_results"
resultsPrefix=sys.argv[1]

outputDirectory="2m_EA_results" # Folder for results

outputFilePathRoot = os.path.join(outputDirectory, resultsPrefix)
gridFileName = outputFilePathRoot+"_max_depth.tif"

threads = 8

if threads is None:
    for d in [0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
        print("Processing threshold ", d)
        depthLabel = f"{d:.1f}"
        depthLabel = depthLabel.replace('.', 'p')
        tableName = resultsPrefix+f"_d{depthLabel}"
        fileIO.uploadGridToDb2(gridFileName, tableName, d, True)
else:
    pool = ThreadPool(threads)
    funcArgList = []

    for d in [0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
        depthLabel = f"{d:.1f}"
        depthLabel = depthLabel.replace('.', 'p')
        tableName = resultsPrefix+f"_d{depthLabel}"

        funcArgList.append((gridFileName, tableName, d))

    pool.map(wrapper, funcArgList)