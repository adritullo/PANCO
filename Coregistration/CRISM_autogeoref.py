import numpy as np
import cv2 as cv
from osgeo import gdal, gdal_array
import matplotlib.pyplot as plt
from skimage import transform, measure, exposure
import shutil


def Sottocampiona2(file_input, output_name, decrement_ratio_x, decrement_ratio_y, format):
    inputinfo = gdal_array.LoadFile(file_input, band_list=[1])
    inputinfo[inputinfo < -2.0] = 0
    minval = np.min(inputinfo)
    maxval = np.max(inputinfo)
    inputinfo = None
    gdal.Translate(output_name, file_input, xRes=decrement_ratio_x, yRes=decrement_ratio_y, outputType=gdal.GDT_Byte, scaleParams=[[minval,maxval,0,255]])


def crism_align(hiri_image, crism_img, crism_if_img, reference_comb, base_format, outname, contrastCrism, contrastCTX, gooddistance):
    standard = {"TRU" : ["R600", "R530", "R440"],
                "123" : ["R770", "RBR", "BD530_2"],
                "VNA" : ["R770", "R770", "R770"],
                "FEM" : ["BD530_2", "SH600_2", "BDI1000"],
                "FM2" : ["BD530_2", "BD920_2", "BDI1000"],
                "TAN" : ["R2529", "IRA", "R770"],
                "IRA" : ["R1300", "R1300", "R1300"],
                "FAL" : ["R2529", "R1506", "R1080"],
                "MAF" : ["OLINDEX3", "LCPINDEX2", "HCPINDEX2"],
                "HYD" : ["SINDEX2", "BD2100_2", "BD1900_2"],
                "PHY" : ["D2300", "D2200", "BD1900R2"],
                "PFM" : ["BD2355", "D2300", "BD2290"],
                "PAL" : ["BD2210_2", "BD2190", "BD2165"],
                "HYS" : ["MIN2250", "BD2250", "BD1900R2"],
                "ICE" : ["BD1900_2", "BD1500_2", "BD1435"],
                "IC2" : ["R3920", "BD1500_2", "BD1435"],
                "CHL" : ["ISLOPE1", "BD3000", "IRR2"],
                "CAR" : ["D2300", "BD2500_2", "BD1900_2"],
                "CR2" : ["MIN2295_2480", "MIN2345_2537", "CINDEX2"]}

    crism_img_lbl = crism_img.replace(".img", ".lbl").replace(".tif", ".lbl")

    if outname.endswith(".img"): outname = outname.replace(".img", ".tif")

    if isinstance(reference_comb, list):
        crism_bands = reference_comb
    else:
        with open(crism_img_lbl, "r") as lbl:
            lbl_read = lbl.readlines()

        for key, val in enumerate(lbl_read):
            if "BAND_NAME" in val:
                goodstart = key
            if " )" in val:
                goodend = key

        #print(goodstart, goodend)

        labellist = []
        for i in range(goodstart, goodend+1):
            labellist.append(lbl_read[i].split('"')[1])
            #print(lbl_read[i].split('"')[1])

        crism_bands = [labellist.index(standard[reference_comb][0])+1, labellist.index(standard[reference_comb][1])+1, labellist.index(standard[reference_comb][2])+1]
        print(reference_comb, ":", labellist.index(standard[reference_comb][0]) + 1, labellist.index(standard[reference_comb][1]) + 1,
              labellist.index(standard[reference_comb][2]) + 1)

    MIN_MATCH_COUNT = 30
    #img1 = cv.imread('outHir.tif', cv.IMREAD_GRAYSCALE)  # queryImage
    #img2 = cv.imread('outCas.tif', cv.IMREAD_GRAYSCALE)  # trainImage


    hiri_subcamp_temp = f"{hiri_image.split('/')[-1].split('.')[0]}_UINT16.tif"
    tempcrism = "tempcrism.tif"


    temp_Hir_Open = gdal.Open(hiri_image)
    hirProj = temp_Hir_Open.GetProjection()
    hirGeotr = temp_Hir_Open.GetGeoTransform()
    hirBands = temp_Hir_Open.RasterCount

    trgPxColumns = temp_Hir_Open.RasterXSize
    trgPxRows = temp_Hir_Open.RasterYSize
    hir_sizepx_x, hir_sizepx_y = hirGeotr[1], hirGeotr[5]


    temp_Hir_Open = None

    temp_CRISM_Open = gdal.Open(crism_img)

    numbands = temp_CRISM_Open.RasterCount

    inputProj = temp_CRISM_Open.GetProjection()
    inputGeotr = temp_CRISM_Open.GetGeoTransform()

    inputGeotr_corr = list(inputGeotr)
    cri_sizepx_x = inputGeotr[1]
    cri_sizepx_y = inputGeotr[5]

    temp_CRISM_Open = None

    #CRIHiri_ratio_x = cri_sizepx_x / hir_sizepx_x
    #CRIHiri_ratio_y = cri_sizepx_y / hir_sizepx_y

    Sottocampiona2(hiri_image, hiri_subcamp_temp, cri_sizepx_x, cri_sizepx_y, base_format)
    # SottocampionaNoRuota(hiri_image, hiri_subcamp_temp, CRIHiri_ratio_x, CRIHiri_ratio_y, base_format)
    # SottocampionaNoRuota(hiri_image, hiri_subcamp_temp, CRIHiri_ratio, gdal.GDT_UInt16)

    rasterArray1_lowreso = gdal_array.LoadFile(crism_img, band_list=[crism_bands[0]])
    rasterArray1_lowreso[rasterArray1_lowreso > 6] = np.nan
    rasterArray1_lowreso[np.isnan(rasterArray1_lowreso)] = 0
    rasterArray1_lowreso[rasterArray1_lowreso < -2] = 0
    rasterArray2_lowreso = gdal_array.LoadFile(crism_img, band_list=[crism_bands[1]])
    rasterArray2_lowreso[rasterArray2_lowreso > 6] = np.nan
    rasterArray2_lowreso[np.isnan(rasterArray2_lowreso)] = 0
    rasterArray2_lowreso[rasterArray2_lowreso < -2] = 0
    rasterArray3_lowreso = gdal_array.LoadFile(crism_img, band_list=[crism_bands[2]])
    rasterArray3_lowreso[rasterArray3_lowreso > 6] = np.nan
    rasterArray3_lowreso[np.isnan(rasterArray3_lowreso)] = 0
    rasterArray3_lowreso[rasterArray3_lowreso < -2] = 0

    rasterArray_pseudopan = rasterArray1_lowreso + rasterArray2_lowreso + rasterArray3_lowreso

    rasterArray_pseudopan_uint16 = (65535 * (rasterArray_pseudopan - np.min(rasterArray_pseudopan)) / np.ptp(rasterArray_pseudopan)).astype(int)

    driver = gdal.GetDriverByName('GTiff')
    rows, cols = rasterArray_pseudopan_uint16.shape
    dataset = driver.Create(tempcrism, cols, rows, 1, gdal.GDT_UInt16)
    dataset.SetProjection(inputProj)
    dataset.SetGeoTransform(inputGeotr_corr)
    dataset.GetRasterBand(1).WriteArray(rasterArray_pseudopan_uint16)
    dataset = None


    #inizia sift


    img1 = cv.imread(tempcrism, cv.IMREAD_GRAYSCALE)  # queryImage

    img1_realmin = np.min(img1[img1 != 0])
    img1_realmax = np.max(img1[img1 != 0])

    img1 = ((img1 - img1_realmin) / (img1_realmax - img1_realmin)) * 255
    img1 = img1.astype(np.uint8)
    img1[img1 > 255] = 0

    plt.imshow(img1), plt.show()

    verSubcamp = gdal.Open(hiri_subcamp_temp)
    verSubcamp_Proj = verSubcamp.GetProjection()
    verSubcamp_Geotr = verSubcamp.GetGeoTransform()
    verSubcamp_sizepx_x, verSubcamp_sizepx_y = verSubcamp_Geotr[1], verSubcamp_Geotr[5]


    if hirBands == 1:
        img2 = cv.imread(hiri_subcamp_temp, cv.IMREAD_GRAYSCALE)  # trainImage
    else:
        img2 = gdal_array.LoadFile(hiri_subcamp_temp, band_list=[1])
        #img2 = gdal_array.LoadFile(hiri_subcamp_temp)
        #img2 = np.sum(img2, axis=0)


        #contrast stretch
        img2_realmin = np.min(img2[img2 != 0])
        img2_realmax = np.max(img2[img2 != 0])

        img2 = ((img2 - img2_realmin) / (img2_realmax - img2_realmin)) * 255
        img2[img2>255] = 0
        img2 = img2.astype(dtype=np.uint8)

        print(img2_realmin, img2_realmax)


    plt.imshow(img2), plt.show()

    sift1 = cv.SIFT_create(contrastThreshold=contrastCrism)
    sift2 = cv.SIFT_create(contrastThreshold=contrastCTX)

    # find the keypoints and descriptors with SIFT

    img1_q01 = np.quantile(img1, 0.001)
    img1_q99 = np.quantile(img1, 0.999)

    img1_mask = np.full_like(img1,1)
    if reference_comb == "HYD":
        img1_mask[img1 == 32] = 0
        img1_mask[img1 <= 5] = 0
        img1_mask[img1 >= 240] = 0
    elif reference_comb == "IC2":
        img1 = 255 - img1
        img1_mask[img1 == img1[0, 0]] = 0
    else:
        img1_mask[img1 <= img1_q01] = 0
        img1_mask[img1 >= img1_q99] = 0



    plt.imshow(img1_mask), plt.show()

    #contrast stretch and remove extremes
    img2_mask = np.full_like(img2, 1)

    img2_q10 = np.quantile(img2, 0.002)
    img2_q90 = np.quantile(img2, 0.998)

    img2_mask[img2 <= img2_q10] = 0
    img2_mask[img2 >= img2_q90] = 0

    if panmask is not None:
        img2_mask[panmask[0][0]:panmask[0][1], panmask[1][0]:panmask[1][1]] = 0



    print(img2_q10, img2_q90)
    plt.imshow(img2_mask), plt.show()

    #plt.imshow(img1), plt.show()
    kp1, des1 = sift1.detectAndCompute(img1, img1_mask)

    #plt.imshow(img2), plt.show()
    kp2, des2 = sift2.detectAndCompute(img2, img2_mask)


    imgKeyp1 = cv.drawKeypoints(img1, kp1, img1, flags=cv.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
    imgKeyp2 = cv.drawKeypoints(img2, kp2, img2, flags=cv.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
    print(len(kp1))
    print(len(kp2))
    cv.imwrite('sift_keypoints1.jpg', imgKeyp1)
    cv.imwrite('sift_keypoints2.jpg', imgKeyp2)

    index_params = dict(algorithm=0, trees=5)
    search_params = dict(checks=500)
    flann = cv.FlannBasedMatcher(index_params, search_params)
    matches = flann.knnMatch(des1, des2, k=2)
    # store all the good matches as per Lowe's ratio test.
    good = []
    for m, n in matches:
        if m.distance < gooddistance * n.distance:
            good.append(m)


    if len(good) > MIN_MATCH_COUNT:
        src_pts = np.float32([kp1[m.queryIdx].pt for m in good]).reshape(-1, 1, 2)
        dst_pts = np.float32([kp2[m.trainIdx].pt for m in good]).reshape(-1, 1, 2)
        M, mask = cv.findHomography(src_pts, dst_pts, cv.RANSAC, 5.0)
        matchesMask = mask.ravel().tolist()
        h, w = img1.shape
        pts = np.float32([[0, 0], [0, h - 1], [w - 1, h - 1], [w - 1, 0]]).reshape(-1, 1, 2)
        dst = cv.perspectiveTransform(pts, M)


        img2_p = cv.polylines(img2, [np.int32(dst)], True, 255, 3, cv.LINE_AA)
        print(f"{len(good)} matches found!")


        draw_params = dict(matchColor=(0, 255, 0),  # draw matches in green color
                           singlePointColor=None,
                           matchesMask=matchesMask,  # draw only inliers
                           flags=2)
        img3 = cv.drawMatches(img1, kp1, img2_p, kp2, good, None, **draw_params)
        plt.imshow(img3, 'gray'), plt.show()
        img3 = None



        #PROIETTO SR

        tform = transform.ProjectiveTransform(matrix=M)
        driver = gdal.GetDriverByName('GTiff')

        rows, cols = img2.shape
        if cols < img1.shape[1]:
            cols = img1.shape[1]
        if rows < img1.shape[0]:
            rows = img1.shape[0]

        dataset = driver.Create(outname, cols, rows, numbands, gdal.GDT_Float32)
        dataset.SetProjection(hirProj)

        newGeotr = (hirGeotr[0], inputGeotr[1], hirGeotr[2], hirGeotr[3], hirGeotr[4], inputGeotr[5])
        dataset.SetGeoTransform(newGeotr)


        for i in range(1,numbands+1):
            rasterArray = gdal_array.LoadFile(crism_img, band_list=[i])
            rasterArray[np.isnan(rasterArray)] = 0
            rasterArray[np.isinf(rasterArray)] = 0
            newRastArray = np.zeros((rows,cols), np.float32)
            origx,origy = rasterArray.shape
            newRastArray[0:origx,0:origy] = rasterArray
            if (i % 10) == 0:
                print(f"{i}/{numbands} done!")
            '''
            if reshaped:
                temp_arr = np.full((rows, cols), 0, np.float32)
                act_x, act_y = rasterArray.shape
                temp_arr[0 : act_x, 0:act_y] = rasterArray
                rasterArray = temp_arr
            '''
            #print(np.min(rasterArray), np.max(rasterArray))
            tf_img = transform.warp(newRastArray, tform.inverse, preserve_range=True)
            #plt.imshow(tf_img), plt.show()
            tf_img[tf_img >= 100] = 65535
            tf_img[tf_img == 0] = 65535
            dataset.GetRasterBand(i).WriteArray(tf_img)
            dataset.GetRasterBand(i).SetNoDataValue(65535)

        shutil.copy(crism_img.replace(".img",".hdr"),outname.replace(".tif",".hdr"))
        shutil.copy(crism_img.replace(".img", ".lbl"), outname.replace(".tif", ".lbl"))
        dataset = None


        # PROIETTO IF

        temp_CRISM_Open = gdal.Open(crism_if_img)
        numbands = temp_CRISM_Open.RasterCount
        temp_Hir_Open = None

        dataset = driver.Create(outname.replace("_sr", "_if"), cols, rows, numbands, gdal.GDT_Float32)
        dataset.SetProjection(hirProj)

        newGeotr = (hirGeotr[0], inputGeotr[1], hirGeotr[2], hirGeotr[3], hirGeotr[4], inputGeotr[5])
        dataset.SetGeoTransform(newGeotr)

        for i in range(1, numbands + 1):
            rasterArray = gdal_array.LoadFile(crism_if_img, band_list=[i])
            rasterArray[np.isnan(rasterArray)] = 0
            rasterArray[np.isinf(rasterArray)] = 0
            newRastArray = np.zeros((rows,cols), np.float32)
            origx, origy = rasterArray.shape
            newRastArray[0:origx, 0:origy] = rasterArray

            if (i % 10) == 0:
                print(f"{i}/{numbands} done!")

            tf_img = transform.warp(newRastArray, tform.inverse, preserve_range=True)
            # plt.imshow(tf_img), plt.show()
            tf_img[tf_img >= 100] = 65535
            tf_img[tf_img == 0] = 65535

            dataset.GetRasterBand(i).WriteArray(tf_img)
            dataset.GetRasterBand(i).SetNoDataValue(65535)

        shutil.copy(crism_if_img.replace(".img", ".hdr"), outname.replace("_sr", "_if").replace(".tif", ".hdr"))
        shutil.copy(crism_if_img.replace(".img", ".lbl"), outname.replace("_sr", "_if").replace(".tif", ".lbl"))
        dataset = None



    else:
        print("Not enough matches are found - {}/{}".format(len(good), MIN_MATCH_COUNT))
        matchesMask = None




# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    pan_image = "your_panchromatic_path"
    crism_sr_img = r"E:\Pansharpening2\HRL_CASSIS\H0\hrl0000baba_07_sr182j_mtr3.img"
    crism_if_img = r"E:\Pansharpening2\HRL_CASSIS\H0\hrl0000baba_07_if182j_mtr3.img"

    # REFERENCE BANDS FOR ALIGNMENT (USUALLY True Color "TRU" IS GOOD FOR CRISM)
    # it can be set as specific band number list or by standard name as Viviano‐Beck et al., 2014. See below for the list
    #reference = [304, 122, 41]
    reference = "TRU"

    # FORMAT OF PAN
    # Usually gdal.GDT_Float32 CaSSIS cubes - gdal.GDT_Float32 ISIS stitched HiRISE cubes - gdal.GDT_Byte PDS CTX - gdal.GDT_UInt16 HRSC - gdal.GDT_Byte HiRISE ASU Mosaics
    base_format = gdal.GDT_Float32
    contrastCrism = 0.01
    contrastCTX = 0.01
    panmask = [[0,400],[0,1500]]
    #panmask = None

    outname = r"your_output file name"
    gooddistance = 0.85

    crism_align(pan_image, crism_sr_img, crism_if_img, reference, base_format, outname, contrastCrism, contrastCTX, gooddistance)


    '''
    standard = {"TRU" : ["R600", "R530", "R440"],
                "VNA" : ["R770", "R770", "R770"],
                "FEM" : ["BD530_2", "SH600_2", "BDI1000"],
                "FM2" : ["BD530_2", "BD920_2", "BDI1000"],
                "TAN" : ["R2529", "IRA", "R770"],
                "IRA" : ["R1300", "R1300", "R1300"],
                "FAL" : ["R2529", "R1506", "R1080"],
                "MAF" : ["OLINDEX3", "LCPINDEX2", "HCPINDEX2"],
                "HYD" : ["SINDEX2", "BD2100_2", "BD1900_2"],
                "PHY" : ["D2300", "D2200", "BD1900r2"],
                "PFM" : ["BD2355", "D2300", "BD2290"],
                "PAL" : ["BD2210_2", "BD2190", "BD2165"],
                "HYS" : ["MIN2250", "BD2250", "BD1900r2"],
                "ICE" : ["BD1900_2", "BD1500_2", "BD1435"],
                "IC2" : ["R3920", "BD1500_2", "BD1435"],
                "CHL" : ["ISLOPE", "BD3000", "IRR2"],
                "CAR" : ["D2300", "BD2500_2", "BD1900_2"],
                "CR2" : ["MIN2295_2480", "MIN2345_2537", "CINDEX2"]}
    '''