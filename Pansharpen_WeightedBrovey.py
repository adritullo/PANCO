import numpy
from osgeo import gdal, ogr, osr, gdal_array
import pandas as pd
import numpy as np
import geopandas
import matplotlib.pyplot as plt


def pansharpenWeightedBrovey(hiri_input, cas_input, outputName = "outputWeightedBrovey.tif", band1 = 1, band2 = 2, band3 = 3):

    #1/3
    red_Weight = 0.166
    green_Weight = 0.167
    blue_Weight = 0.167
    NIR_Weight = 0.5

    hiriInfo = gdal.Open(hiri_input)

    trgPxColumns = hiriInfo.RasterXSize
    trgPxRows = hiriInfo.RasterYSize
    print(f"Columns: {trgPxColumns}, Rows: {trgPxRows}")

    temp_cas = "temp_out.tif"
    gdal.Translate(temp_cas, cas_input, format="GTiff", outputType=gdal.gdalconst.GDT_Float32, width=trgPxColumns,
                   height=trgPxRows)

    rasterArray1 = gdal_array.LoadFile(temp_cas, band_list=[band1])
    rasterArray1[rasterArray1 > 2] = np.nan
    rasterArray2 = gdal_array.LoadFile(temp_cas, band_list=[band2])
    rasterArray2[rasterArray2 > 2] = np.nan
    rasterArray3 = gdal_array.LoadFile(temp_cas, band_list=[band3])
    rasterArray3[rasterArray3 > 2] = np.nan

    rasterArrayHiri = gdal_array.LoadFile(hiri_input, band_list=[1])

    weighted_R = np.multiply(rasterArray1, red_Weight)
    weighted_G = np.multiply(rasterArray2, green_Weight)
    weighted_B = np.multiply(rasterArray3, blue_Weight)

    psuedo_pancro = np.nansum((weighted_R, weighted_G, weighted_B), axis=0, dtype=np.float32)

    ratio = np.divide(rasterArrayHiri, psuedo_pancro)

    #synthetic intensity component I

    corrected_r = np.multiply(rasterArray1, ratio)
    corrected_g = np.multiply(rasterArray2, ratio)
    corrected_b = np.multiply(rasterArray3, ratio)

    print(corrected_r[500, 500])
    print(corrected_g[500, 500])
    print(corrected_b[500, 500])

    # Detail image P - I

    stacked_array = np.stack([corrected_r, corrected_g, corrected_b], axis=0)

    driver = gdal.GetDriverByName('GTiff')
    n, rows, cols = stacked_array.shape
    dataset = driver.Create(outputName, cols, rows, n,
                            gdal.GDT_Float32)

    for b in range(1, n + 1):
        band = dataset.GetRasterBand(b)  # GetRasterBand is not zero indexed
        band.WriteArray(stacked_array[b - 1])  # Numpy is zero indexed

    dataset = None
    stacked_array = None

    temp_cas_Open = gdal.Open(temp_cas)
    out_cas_Open = gdal.Open(outputName)

    gdal_array.CopyDatasetInfo(temp_cas_Open, out_cas_Open)
    temp_cas_Open = out_cas_Open = None


if __name__ == '__main__':
    PAN_input = "your_panchromatic_path"
    MS_input = "your_multispectral_path"
    outputname = "your_output_path/Pansharpened.tif"
    band1, band2, band3 = 1, 2, 3  # The three bands you want to pansharpen. The output will have the same order you choose

    pansharpenWeightedBrovey(PAN_input, MS_input)
