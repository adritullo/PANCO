import numpy as np
import scipy.ndimage as sn
from osgeo import gdal, gdal_array


def pansharpenGramSchmidt(hiri_input, cas_input, outputName = "outputGramSchmidt.tif", band1 = 1, band2 = 2, band3 = 3):

    red_Weight = 0.5
    green_Weight = 0.5
    blue_Weight = 0.5
    NIR_Weight = 0.5


    #Computing weights

    rasterArray1_lowreso = gdal_array.LoadFile(cas_input, band_list=[band1])
    rasterArray1_lowreso[rasterArray1_lowreso > 5] = np.nan
    rasterArray1_lowreso[np.isnan(rasterArray1_lowreso)] = 0
    rasterArray1_lowreso[rasterArray1_lowreso < -1] = 0
    rasterArray2_lowreso = gdal_array.LoadFile(cas_input, band_list=[band2])
    rasterArray2_lowreso[rasterArray2_lowreso > 3] = np.nan
    rasterArray2_lowreso[np.isnan(rasterArray2_lowreso)] = 0
    rasterArray2_lowreso[rasterArray2_lowreso < -1] = 0
    rasterArray3_lowreso = gdal_array.LoadFile(cas_input, band_list=[band3])
    rasterArray3_lowreso[rasterArray3_lowreso > 3] = np.nan
    rasterArray3_lowreso[np.isnan(rasterArray3_lowreso)] = 0
    rasterArray3_lowreso[rasterArray3_lowreso < -1] = 0

    weighted_R = np.multiply(rasterArray1_lowreso, red_Weight)
    weighted_G = np.multiply(rasterArray2_lowreso, green_Weight)
    weighted_B = np.multiply(rasterArray3_lowreso, blue_Weight)

    sim_pancro_low = np.sum((weighted_R, weighted_G, weighted_B), axis=0, dtype=np.float32)

    sim_pancro_low[np.isnan(sim_pancro_low)] = 0

    sim_pancro_low = sim_pancro_low - np.nanmean(sim_pancro_low)

    ra1_mean = np.nanmean(rasterArray1_lowreso)
    ra2_mean = np.nanmean(rasterArray2_lowreso)
    ra3_mean = np.nanmean(rasterArray3_lowreso)

    ra1_lr_dcfree = np.subtract(rasterArray1_lowreso, ra1_mean)
    ra2_lr_dcfree = np.subtract(rasterArray2_lowreso, ra2_mean)
    ra3_lr_dcfree = np.subtract(rasterArray3_lowreso, ra3_mean)


    x1v  = sim_pancro_low.flatten()
    x1 = v1 = sim_pancro_low

    x2 = ra1_lr_dcfree
    x2v = ra1_lr_dcfree.flatten()

    v2 = x2 - (x2v @ x1v) / (x1v @ x1v) * x1
    v2v = v2.flatten()

    x3 = ra2_lr_dcfree
    x3v = ra2_lr_dcfree.flatten()

    v3 = x3 - (x3v @ x1v) / (x1v @ x1v) * x1 + (x3v @ v2v) / (v2v @ v2v) * v2
    v3v = v3.flatten()

    x4 = ra3_lr_dcfree
    x4v = ra3_lr_dcfree.flatten()

    v4 = x4 - (x4v @ x1v) / (x1v @ x1v) * x1 + (x4v @ v2v) / (v2v @ v2v) * v2 + (x4v @ v3v) / (v3v @ v3v) * v3

    hiriInfo = gdal.Open(hiri_input)
    trgPxColumns = hiriInfo.RasterXSize
    trgPxRows = hiriInfo.RasterYSize
    print(f"Columns: {trgPxColumns}, Rows: {trgPxRows}")

    temp_cas = "temp_out.tif"
    gdal.Translate(temp_cas, cas_input, format="GTiff", outputType=gdal.gdalconst.GDT_Float32, width=trgPxColumns,
                   height=trgPxRows)


    rasterArrayHiri = gdal_array.LoadFile(hiri_input, band_list=[1])
    '''
    rasterArrayHiri[rasterArrayHiri > 3] = np.nan
    rasterArrayHiri[rasterArrayHiri < -0.1] = np.nan
    rasterArrayHiri[np.isnan(rasterArrayHiri)] = 0
    rasterArrayHiri[np.isinf(rasterArrayHiri)] = 0
    '''


    '''
    hiri_min, hiri_max = np.nanmin(rasterArrayHiri), np.nanmax(rasterArrayHiri)

    i_max = np.max(v1)
    i_min = np.min(v1)

    hiri_std = (rasterArrayHiri - hiri_min) / (hiri_max - hiri_min)
    normalized_PAN = (hiri_std * (i_max-i_min)) + i_min
    '''

    v1_mean = np.mean(x1)
    pan_mean = np.mean(rasterArrayHiri)

    v1_std = np.std(x1)
    pan_std = np.std(rasterArrayHiri)

    print(v1_mean, pan_mean)
    print(v1_std, pan_std)

    gain = v1_std/pan_std

    bias = v1_mean - (gain * pan_mean)

    pan_stretched = (rasterArrayHiri * gain) + bias

    print("after: ", rasterArrayHiri.shape, pan_stretched.shape)

    pan_mean = np.mean(pan_stretched)
    pan_std = np.std(pan_stretched)

    print(v1_mean, pan_mean)
    print(v1_std, pan_std)

    x1u = pan_stretched
    x1vu = pan_stretched.flatten()

    '''
    fig = plt.figure()
    ax1 = fig.add_subplot(121)  # left side
    ax2 = fig.add_subplot(122)  # right side
    ax1.imshow(v1)
    ax2.imshow(x1u)
    plt.show()
    '''

    upscalefactor = ((pan_stretched.shape[0] / v1.shape[0]) + (pan_stretched.shape[1] / v1.shape[1])) /2
    print(upscalefactor)


    v2_up = sn.zoom(v2, zoom=((pan_stretched.shape[0] / v1.shape[0]), (pan_stretched.shape[1] / v1.shape[1])), order=1)
    v3_up = sn.zoom(v3, zoom=((pan_stretched.shape[0] / v1.shape[0]), (pan_stretched.shape[1] / v1.shape[1])), order=1)
    v4_up = sn.zoom(v4, zoom=((pan_stretched.shape[0] / v1.shape[0]), (pan_stretched.shape[1] / v1.shape[1])), order=1)


    x2_up = (v2_up + ra1_mean) + ((x1v @ x1v) / (x1v @ x1v) * x1u)

    x3_up = (v3_up + ra2_mean) + (x3v @ x1v) / (x1v @ x1v) * x1u + (x3v @ v2v) / (v2v @ v2v) * v2_up

    x4_up = (v4_up + ra3_mean) + (x4v @ x1v) / (x1v @ x1v) * x1u + (x4v @ v2v) / (v2v @ v2v) * v2_up + (x4v @ v3v) / (v3v @ v3v) * v3_up


    '''
    x2_up = (v2_up + ra1_mean) + ((v2_up.flatten() @ x1vu) / (x1vu @ x1vu) * x1u)

    x3_up = (v3_up + ra2_mean) + (v3_up.flatten() @ x1vu) / (x1vu @ x1vu) * x1u + (x3v @ v2v) / (v2v @ v2v) * v2_up

    x4_up = (v4_up + ra3_mean) + (v4_up.flatten() @ x1vu) / (x1vu @ x1vu) * x1u + (x4v @ v2v) / (v2v @ v2v) * v2_up + (x4v @ v3v) / (v3v @ v3v) * v3_up
    '''

    stacked_array = np.stack([x2_up, x3_up, x4_up], axis=0)

    driver = gdal.GetDriverByName('GTiff')
    n, rows, cols = stacked_array.shape
    dataset = driver.Create(outputName, cols, rows, n,
                            gdal.GDT_Float32)

    for b in range(1, n + 1):
        band = dataset.GetRasterBand(b)  # GetRasterBand is not zero indexed
        band.WriteArray(stacked_array[b - 1])  # Numpy is zero indexed

    temp_cas_Open = gdal.Open(temp_cas)

    temp_cas_Proj = temp_cas_Open.GetProjection()
    temp_cas_Geotr = temp_cas_Open.GetGeoTransform()

    dataset = None
    stacked_array = None

    out_cas_Open = gdal.Open(outputName)

    gdal_array.CopyDatasetInfo(temp_cas_Open, out_cas_Open)

    out_cas_Proj = out_cas_Open.GetProjection()
    out_cas_Geotr = out_cas_Open.GetGeoTransform()


    print(temp_cas_Proj, temp_cas_Geotr)
    print(out_cas_Proj, out_cas_Geotr)


    temp_cas_Open = out_cas_Open = None



if __name__ == '__main__':
    PAN_input = "your_panchromatic_path"
    MS_input = "your_multispectral_path"
    outputname = "your_output_path/Pansharpened.tif"

    band1, band2, band3 = 1,2,3  # The three bands you want to pansharpen. The output will have the same order you choose
    pansharpenGramSchmidt(PAN_input, MS_input, outputname, band1, band2, band3)
