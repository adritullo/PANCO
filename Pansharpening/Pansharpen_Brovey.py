import numpy
from osgeo import gdal, ogr, osr, gdal_array
import pandas as pd
import numpy as np
import geopandas
import matplotlib.pyplot as plt


def pansharpenBrovey(hiri_input, cas_input, outputName = "outputBrovey.tif", band1 = 1, band2 = 2, band3 = 3):
    #1/3

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

    reduct_Image = np.nanmean((rasterArray1, rasterArray2, rasterArray3), axis=0, dtype=np.float32)

    detail_Image = np.subtract(rasterArrayHiri, reduct_Image, dtype=np.float32)

    # Synthetic intensity component I

    correction_r = np.multiply(detail_Image, numpy.divide(rasterArray1, reduct_Image))
    correction_g = np.multiply(detail_Image, numpy.divide(rasterArray2, reduct_Image))
    correction_b = np.multiply(detail_Image, numpy.divide(rasterArray3, reduct_Image))

    # Detail image P - I

    red_out = np.nansum((rasterArray1, correction_r), axis=0, dtype=np.float32)
    gre_out = np.nansum((rasterArray2, correction_g), axis=0, dtype=np.float32)
    blu_out = np.nansum((rasterArray3, correction_b), axis=0, dtype=np.float32)

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

    band1, band2, band3 = 1,2,3  # The three bands you want to pansharpen. The output will have the same order you choose
    pansharpenBrovey(PAN_input, MS_input, outputname, band1, band2, band3)
