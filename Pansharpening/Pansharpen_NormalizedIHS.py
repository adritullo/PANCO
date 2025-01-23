import numpy
from osgeo import gdal, ogr, osr, gdal_array
import pandas as pd
import numpy as np
import geopandas
import matplotlib.pyplot as plt
import sklearn.preprocessing as sk


def pansharpenNormIHS(hiri_input, cas_input, outputName = "outputNormIHS4.tif", band1 = 1, band2 = 2, band3 = 3):

    #1/3
    red_Weight = 0.3333
    green_Weight = 0.3333
    blue_Weight = 0.3333
    NIR_Weight = 0


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


    #synthetic intensity component I
    i_r = np.multiply(rasterArray1, red_Weight)
    i_g = np.multiply(rasterArray2, green_Weight)
    i_b = np.multiply(rasterArray3, blue_Weight)

    intensity_image = np.sum((i_r,i_g,i_b), axis=0, dtype=np.float32)

    i_min, i_max = numpy.nanmin(intensity_image), numpy.nanmax(intensity_image)

    print(i_min, i_max)


    #SOLUZIONE MANUALE
    #X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    #X_scaled = X_std * (max - min) + min

    hiri_min, hiri_max = numpy.nanmin(rasterArrayHiri), numpy.nanmax(rasterArrayHiri)

    hiri_std = (rasterArrayHiri - hiri_min) / (hiri_max - hiri_min)
    normalized_PAN = (hiri_std * (i_max-i_min)) + i_min

    print(hiri_min, hiri_max)

    print(rasterArrayHiri[500,500])
    print(normalized_PAN[500,500])

    # Detail image P - I
    detail_Image = np.subtract(normalized_PAN, intensity_image, dtype=np.float32)

    red_out = np.nansum((rasterArray1, detail_Image), axis=0, dtype=np.float32)
    gre_out = np.nansum((rasterArray2, detail_Image), axis=0, dtype=np.float32)
    blu_out = np.nansum((rasterArray3, detail_Image), axis=0, dtype=np.float32)

    stacked_array = np.stack([red_out, gre_out, blu_out], axis=0)

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

    pansharpenNormIHS(PAN_input, MS_input, outputname, band1, band2, band3)

