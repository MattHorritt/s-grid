###############################################################################
# Various utilities for LTIS sea level rise and coastal risk project

import sgrid
import fileIO
import numpy as np
from pathlib import Path

# Some useful constants
minWaveHeight = 0.2  # Height required above toe level required to damage defence
tidePeriod = 12.42  # Period between high tides


###############################################################################
# Read defence polyline, identify defences at risk of erosion, and generate
# points as list of (i, j, t1, t200) tuples. Also save to CSV as diagnostic.
def slr_tide_points(slr, defenceLineStr, xll, yll, xsz, ysz, cellSize, outputDirectory, outputPrefix, layerName=None):
    tidePointsCsvFile = open(str(Path(outputDirectory)/ (outputPrefix + '_tidePoints.csv')), "w")
    tidePointsCsvFile.write("X,Y,XC,YC,T1,T200,ToeLevel\n")

    xyiList = []

    tidePointsRaw = fileIO.readPolylineShapefile(defenceLineStr, layerName=layerName)
    nl = len(tidePointsRaw)
    tidePoints = []
    tidePointsN = 0
    for i in range(nl):
        xl = tidePointsRaw[i][0]
        yl = tidePointsRaw[i][1]

        t1 = tidePointsRaw[i][2]["t1"]
        t200 = tidePointsRaw[i][2]["t200"]

        toeLevel = tidePointsRaw[i][2]["toe_level_avg"]

        if toeLevel > (slr - minWaveHeight):
            continue

        for j in range(len(xl) - 1):
            x1 = xl[j]
            y1 = yl[j]
            x2 = xl[j + 1]
            y2 = yl[j + 1]

            segmentLength = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            nSteps = max(2, int(2 * segmentLength / cellSize))  # Always 2 points

            for k in range(nSteps):
                x = x1 + (x2 - x1) * k / nSteps
                y = y1 + (y2 - y1) * k / nSteps

                xi = int((x - xll) / cellSize)
                yi = int((y - yll) / cellSize)

                if xi >= 0 and xi < xsz and yi >= 0 and yi < ysz:
                    if (xi, yi) not in xyiList:
                        tidePointsN += 1
                        tidePoints.append((xi, yi, t1, t200))
                        xyiList.append((xi, yi))

                        # Save centres of cells too
                        xc = xll + xi * cellSize + cellSize / 2
                        yc = yll + yi * cellSize + cellSize / 2

                        tidePointsCsvFile.write("%f,%f,%f,%f,%f,%f,%f\n" % (x, y, xc, yc, t1, t200, toeLevel))

    tidePointsCsvFile.close()

    return tidePoints


###############################################################################
# Calculate tide level across multiple high tides; peak is applied to 2nd high
# tide.
def applyTideLevels(tidePoints, t, dt, slr, storagePar, wlGrid, volGrid):
    # For tracking flows in and out
    Qin = 0.
    Qout = 0.

    for tidePt in tidePoints:
        if storagePar[tidePt[0], tidePt[1], 0] == -9999:
            continue

        theta = 2 * np.pi * t / (3600. * tidePeriod)

        if theta <= 2 * np.pi or theta > 3 * np.pi:
            wl = np.sin(theta) * tidePt[2]
        else:
            wl = np.sin(theta) * tidePt[3]

        wl += slr

        wlGrid[tidePt[0], tidePt[1]] = wl
        newV = sgrid.volFromWl(wl, tidePt[0], tidePt[1], storagePar)

        Qtide = (newV - volGrid[tidePt[0], tidePt[1]]) / dt

        if Qtide < 0:
            Qout -= Qtide
        else:
            Qin += Qtide

        volGrid[tidePt[0], tidePt[1]] = newV

    return Qin, Qout

