from osgeo import gdal, ogr, osr, gdal_array
import pandas as pd
import numpy as np
import geopandas
import matplotlib.pyplot as plt
import scipy.ndimage as sn


def pansharpenHPF(hiri_input, cas_input, outputName = "outputHPF.tif", band1 = 1, band2 = 2, band3 = 3):


    # Low Reso HiRISE

    casInfo = gdal.Open(cas_input)
    low_trgPxColumns = casInfo.RasterXSize
    low_trgPxRows = casInfo.RasterYSize
    print(f"LowRes Columns: {low_trgPxColumns}, Rows: {low_trgPxRows}")
    '''
    tempLowHiri = "tempLowHiri.tif"
    gdal.Translate(tempLowHiri, hiri_input, format="GTiff", outputType=gdal.gdalconst.GDT_Float32, 
                   width=low_trgPxColumns, height=low_trgPxRows)
    '''

    hiriInfo = gdal.Open(hiri_input)
    hi_trgPxColumns = hiriInfo.RasterXSize
    hi_trgPxRows = hiriInfo.RasterYSize
    print(f"HighRes Columns: {hi_trgPxColumns}, Rows: {hi_trgPxRows}")
    increment_ratio = hi_trgPxRows / low_trgPxRows

    tempHiCas = "temp_out.tif"
    gdal.Translate(tempHiCas, cas_input, format="GTiff", outputType=gdal.gdalconst.GDT_Float32, width=hi_trgPxColumns,
                   height=hi_trgPxRows)


    rasterArray1 = gdal_array.LoadFile(tempHiCas, band_list=[band1])
    rasterArray1[rasterArray1 > 3] = 0
    rasterArray1[rasterArray1 < -3] = 0
    rasterArray2 = gdal_array.LoadFile(tempHiCas, band_list=[band2])
    rasterArray2[rasterArray2 > 3] = 0
    rasterArray1[rasterArray2 < -3] = 0
    rasterArray3 = gdal_array.LoadFile(tempHiCas, band_list=[band3])
    rasterArray3[rasterArray3 > 3] = 0
    rasterArray1[rasterArray3 < -3] = 0

    rasterArrayHiri = gdal_array.LoadFile(hiri_input, band_list=[1])


    #STIMARE UN MIGLIORE RAPPORTO TRA DIFFERENZA FREQUENCY RESPONSE CUTOFF
    boxcar_size = round(increment_ratio) + 1

    print(f"Rapporto {increment_ratio}, Boxcar Filter Size {boxcar_size}x{boxcar_size}")

    lowPassArrayHiRI = sn.uniform_filter(rasterArrayHiri, size=boxcar_size)
    '''
    driver = gdal.GetDriverByName('GTiff')
    rows, cols = lowPassArrayHiRI.shape
    dataset = driver.Create("lowPass9.tif", cols, rows, 1, gdal.GDT_Float32)
    band = dataset.GetRasterBand(1)  # GetRasterBand is not zero indexed
    band.WriteArray(lowPassArrayHiRI)  # Numpy is zero indexed
    dataset = None
    lowPassArrayHiRI = None
    out_cas_Open = gdal.Open("lowPass9.tif")
    gdal_array.CopyDatasetInfo(hiriInfo, out_cas_Open)
    '''

    # Detail image P - lowpass filtered P
    detail_Image = np.subtract(rasterArrayHiri, lowPassArrayHiRI, dtype=np.float32)

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

    temp_cas_Open = gdal.Open(tempHiCas)
    out_cas_Open = gdal.Open(outputName)

    gdal_array.CopyDatasetInfo(temp_cas_Open, out_cas_Open)

    temp_cas_Open = out_cas_Open = None



if __name__ == '__main__':
    PAN_input = "your_panchromatic_path"
    MS_input = "your_multispectral_path"
    outputname = "your_output_path/Pansharpened.tif"
    band1, band2, band3 = 1,2,3  # The three bands you want to pansharpen. The output will have the same order you choose

    pansharpenHPF(PAN_input, MS_input, outputname, band1, band2, band3)