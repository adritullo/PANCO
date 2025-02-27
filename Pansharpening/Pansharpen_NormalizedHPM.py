from osgeo import gdal, ogr, osr, gdal_array
import pandas as pd
import numpy as np
import geopandas
import matplotlib.pyplot as plt
import scipy.ndimage as sn


def PansharpenNormHPM(hiri_input, cas_input, outputName = "outputNormHPM.tif", band1 = 1, band2 = 2, band3 = 3):

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
    rasterArray1[rasterArray1 > 2] = np.nan
    rasterArray2 = gdal_array.LoadFile(tempHiCas, band_list=[band2])
    rasterArray2[rasterArray2 > 2] = np.nan
    rasterArray3 = gdal_array.LoadFile(tempHiCas, band_list=[band3])
    rasterArray3[rasterArray3 > 2] = np.nan

    rasterArrayHiri = gdal_array.LoadFile(hiri_input, band_list=[1])

    # synthetic intensity component I
    intensity_image = np.sum((rasterArray1, rasterArray2, rasterArray3), axis=0, dtype=np.float32)

    i_min, i_max = np.nanmin(intensity_image), np.nanmax(intensity_image)

    print(i_min, i_max)

    # SOLUZIONE MANUALE
    # X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    # X_scaled = X_std * (max - min) + min

    hiri_min, hiri_max = np.nanmin(rasterArrayHiri), np.nanmax(rasterArrayHiri)

    hiri_std = (rasterArrayHiri - hiri_min) / (hiri_max - hiri_min)
    normalized_PAN = (hiri_std * (i_max - i_min)) + i_min

    # STIMARE UN MIGLIORE RAPPORTO TRA DIFFERENZA FREQUENCY RESPONSE CUTOFF
    boxcar_size = round(increment_ratio) + 1

    print(f"Rapporto {increment_ratio}, Boxcar Filter Size {boxcar_size}x{boxcar_size}")

    lowPassArrayHiRI = sn.uniform_filter(normalized_PAN, size=boxcar_size)
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
    highpass_image = np.subtract(normalized_PAN, lowPassArrayHiRI, dtype=np.float32)
    ratio_r = np.divide(rasterArray1, lowPassArrayHiRI)
    ratio_g = np.divide(rasterArray2, lowPassArrayHiRI)
    ratio_b = np.divide(rasterArray3, lowPassArrayHiRI)

    detail_image_r = np.multiply(ratio_r, highpass_image)
    detail_image_g = np.multiply(ratio_g, highpass_image)
    detail_image_b = np.multiply(ratio_b, highpass_image)


    red_out = np.nansum((rasterArray1, detail_image_r), axis=0, dtype=np.float32)
    gre_out = np.nansum((rasterArray2, detail_image_g), axis=0, dtype=np.float32)
    blu_out = np.nansum((rasterArray3, detail_image_b), axis=0, dtype=np.float32)

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
    band1, band2, band3 = 1, 2, 3  # The three bands you want to pansharpen. The output will have the same order you choose

    PansharpenNormHPM(PAN_input, MS_input, outputname, band1, band2, band3)
