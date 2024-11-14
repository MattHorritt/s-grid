
#!/usr/bin/env python

import sys, os
import fileIO

#############################################################################
# Control Panel

resultsDirectory="2m_EA_results"
resultsPrefix=sys.argv[1]

outputDirectory="2m_EA_results" # Folder for results

outputFilePathRoot = os.path.join(outputDirectory, resultsPrefix)
gridFileName = outputFilePathRoot+"_max_depth.tif"

for d in [0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
    print("Processing threshold ", d)
    depthLabel = f"{d:.1f}"
    depthLabel = depthLabel.replace('.', 'p')
    tableName = resultsPrefix+f"_d{depthLabel}"
    fileIO.uploadGridToDb2(gridFileName, tableName, depthThreshold=d, dropTable=True)