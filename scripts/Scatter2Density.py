import numpy as np
from scipy.stats import gaussian_kde

def Scatter2Density(arrX, arrY):

    XY = np.vstack([arrX, arrY])
    z_XY = gaussian_kde(XY)(XY)
    idx_Density_XY = z_XY.argsort()
    arrX, arrY, z_XY = arrX[idx_Density_XY], arrY[idx_Density_XY], z_XY[idx_Density_XY]

    return arrX, arrY, z_XY