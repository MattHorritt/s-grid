import fileIO
import numpy as np

arr = np.zeros((2000, 1000))
arr[1750:1950, 50:950] = 100
arr[0:1700, 650:] = 100
arr[0:1700, 0:350] = 100
for i in range(0, 1700):
    arr[i, :] = arr[i, :] + (1700 - i) * 10 * 0.001

fileIO.saveScalarGrid(arr, 0, 0, 10, "test_grid_100m.tif")

# NaN version
arr[np.where(arr > 50)] = -9999
fileIO.saveScalarGrid(arr, 0, 0, 10, "test_grid_NaN.tif")
