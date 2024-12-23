from osgeo import gdal, ogr, osr, gdal_array
import pandas as pd
import numpy as np
import geopandas
import matplotlib.pyplot as plt
import scipy.ndimage as sn
from scipy import integrate
import skimage
from sklearn.linear_model import LinearRegression


def pansharpenMMSEBrovey(hiri_input, cas_input, outputName = "outputMMSEBrovey.tif", band1=1, band2=2, band3=3):

    #MSE MEAN minimum mean-square error (mse)
    # Low Reso HiRISE

    casInfo = gdal.Open(cas_input)
    low_trgPxColumns = casInfo.RasterXSize
    low_trgPxRows = casInfo.RasterYSize
    print(f"LowRes Columns: {low_trgPxColumns}, Rows: {low_trgPxRows}")
    tempLowHiri = "tempLowHiri.tif"


    hiriInfo = gdal.Open(hiri_input)
    hi_trgPxColumns = hiriInfo.RasterXSize
    hi_trgPxRows = hiriInfo.RasterYSize
    print(f"HighRes Columns: {hi_trgPxColumns}, Rows: {hi_trgPxRows}")
    increment_ratio = hi_trgPxRows / low_trgPxRows

    tempHiCas = "temp_out.tif"
    gdal.Translate(tempHiCas, cas_input, format="GTiff", outputType=gdal.gdalconst.GDT_Float32, width=hi_trgPxColumns,
                   height=hi_trgPxRows)


    rasterArray1 = gdal_array.LoadFile(tempHiCas, band_list=[band1])
    rasterArray1[rasterArray1 < -3] = 0.0
    rasterArray2 = gdal_array.LoadFile(tempHiCas, band_list=[band2])
    rasterArray2[rasterArray2 < -3] = 0.0
    rasterArray3 = gdal_array.LoadFile(tempHiCas, band_list=[band3])
    rasterArray3[rasterArray3 < -3] = 0.0

    rasterArray1_low = gdal_array.LoadFile(cas_input, band_list=[band1])
    rasterArray1_low[rasterArray1_low < -3] = 0.0
    rasterArray2_low = gdal_array.LoadFile(cas_input, band_list=[band2])
    rasterArray2_low[rasterArray2_low < -3] = 0.0
    rasterArray3_low = gdal_array.LoadFile(cas_input, band_list=[band3])
    rasterArray3_low[rasterArray3_low < -3] = 0.0

    rasterArray_Hi_Hiri = gdal_array.LoadFile(hiri_input, band_list=[1])

    boxcar_size = round(increment_ratio) + 1
    lowPassArrayHiRI = sn.uniform_filter(rasterArray_Hi_Hiri, size=boxcar_size)

    rasterArray_Low_Hiri = skimage.measure.block_reduce(lowPassArrayHiRI, block_size=(round(increment_ratio), round(increment_ratio)), func=np.mean)

    '''
    driver = gdal.GetDriverByName('GTiff')
    rows, cols = rasterArray_Low_Hiri.shape
    dataset = driver.Create("Pcurvo.tif", cols, rows, 1, gdal.GDT_Float32)
    band = dataset.GetRasterBand(1)  # GetRasterBand is not zero indexed
    band.WriteArray(rasterArray_Low_Hiri)  # Numpy is zero indexed
    dataset = None
    lowPassArrayHiRI = None
    out_cas_Open = gdal.Open("Pcurvo.tif")
    gdal_array.CopyDatasetInfo(casInfo, out_cas_Open)
    '''
    stacked_array = np.stack([rasterArray1, rasterArray2, rasterArray3], axis=0)

    rasterArray1_low = np.nan_to_num(rasterArray1_low)
    rasterArray2_low = np.nan_to_num(rasterArray2_low)
    rasterArray3_low = np.nan_to_num(rasterArray3_low)
    rasterArray_Low_Hiri = np.nan_to_num(rasterArray3_low)

    rasterArray1_low[rasterArray1_low == np.inf] = 0.1
    rasterArray2_low[rasterArray2_low == np.inf] = 0.1
    rasterArray3_low[rasterArray3_low == np.inf] = 0.1
    rasterArray_Low_Hiri[rasterArray_Low_Hiri == np.inf] = 0.1

    model = LinearRegression()
    model.fit(rasterArray1_low, rasterArray_Low_Hiri)
    r_sq1 = model.score(rasterArray1_low, rasterArray_Low_Hiri)
    print(f"coefficient of determination: {r_sq1}")
    red_Weight = np.mean(model.intercept_)
    print(f"intercept: {model.intercept_}")
    print(f"slope: {model.coef_}")

    model2 = LinearRegression()
    model2.fit(rasterArray2_low, rasterArray_Low_Hiri)
    r_sq2 = model2.score(rasterArray2_low, rasterArray_Low_Hiri)
    print(f"coefficient of determination: {r_sq2}")
    green_Weight = np.mean(model2.intercept_)
    print(f"intercept: {model2.intercept_}")
    print(f"slope: {model2.coef_}")

    model3 = LinearRegression()
    model3.fit(rasterArray3_low, rasterArray_Low_Hiri)
    r_sq3 = model3.score(rasterArray3_low, rasterArray_Low_Hiri)
    print(f"coefficient of determination: {r_sq3}")
    blue_Weight = np.mean(model3.intercept_)
    print(f"intercept: {model3.intercept_}")
    print(f"slope: {model3.coef_}")

    print(red_Weight, green_Weight, blue_Weight)

    i_r = np.multiply(rasterArray1, red_Weight)
    i_g = np.multiply(rasterArray2, green_Weight)
    i_b = np.multiply(rasterArray3, blue_Weight)

    intensity_image = np.nansum((i_r,i_g,i_b), axis=0, dtype=np.float32)

    i_min, i_max = np.nanmin(intensity_image), np.nanmax(intensity_image)

    # SOLUZIONE MANUALE
    # X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    # X_scaled = X_std * (max - min) + min

    hiri_min, hiri_max = np.nanmin(rasterArray_Hi_Hiri), np.nanmax(rasterArray_Hi_Hiri)

    hiri_std = (rasterArray_Hi_Hiri - hiri_min) / (hiri_max - hiri_min)
    normalized_PAN = (hiri_std * (i_max - i_min)) + i_min


    ratio = np.divide(normalized_PAN, intensity_image)

    corrected_r = np.multiply(rasterArray1, ratio)
    corrected_g = np.multiply(rasterArray2, ratio)
    corrected_b = np.multiply(rasterArray3, ratio)

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

    temp_cas_Open = gdal.Open(tempHiCas)
    out_cas_Open = gdal.Open(outputName)

    gdal_array.CopyDatasetInfo(temp_cas_Open, out_cas_Open)
    temp_cas_Open = out_cas_Open = None



if __name__ == '__main__':
    PAN_input = "your_panchromatic_path"
    MS_input = "your_multispectral_path"
    outputname = "your_output_path/Pansharpened.tif"
    band1, band2, band3 = 1,2,3  # The three bands you want to pansharpen. The output will have the same order you choose

    pansharpenMMSEBrovey(PAN_input, MS_input, outputname, band1, band2, band3)
