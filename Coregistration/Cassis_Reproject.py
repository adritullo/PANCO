import os
import numpy as np
import cv2 as cv
from skimage import transform
from osgeo import gdal, gdal_array
from numba import njit
from time import perf_counter, sleep
import matplotlib.pyplot as plt
import pandas as pd


def buildPyramids(image_input):
    image = gdal.Open(image_input, 0)  # 0 = read-only, 1 = read-write.
    gdal.SetConfigOption('COMPRESS_OVERVIEW', 'DEFLATE')
    image.BuildOverviews('NEAREST', [4, 8, 16, 32, 64, 128], gdal.TermProgress_nocb)
    del image  # close the dataset (Python object and pointers)


def createScaleMatrix(sx, sy):
    return np.array([
        [sx, 0, 0],
        [0, sy, 0],
        [0, 0, 1]
    ])


def createTranslationMatrix(tx, ty):
    return np.array([
        [1, 0, tx],
        [0, 1, ty],
        [0, 0, 1]
    ])


def applyGeoTransform(geotransform, x, y, h):
    return (geotransform[0] + geotransform[1] * x + geotransform[2] * y,
            geotransform[3] + geotransform[4] * x + geotransform[5] * y, h)


@njit()
def RuotaArray(inarray):
    return np.rot90(inarray, 2)


def SottocampionaMinMax(file_input, output_name, decrement_ratio_x, decrement_ratio_y, band):
    inputinfo = gdal_array.LoadFile(file_input, band_list=band).astype(float)
    inputinfo[inputinfo <= 0] = np.nan
    minval = np.nanmin(inputinfo)
    maxval = np.nanmax(inputinfo)
    print(f"Image range: {minval} - {maxval}")
    inputinfo = None
    gdal.Translate(output_name, file_input, xRes=decrement_ratio_x, yRes=decrement_ratio_y, outputType=gdal.GDT_Byte,
                   scaleParams=[[minval, maxval, 0, 255]])


def SottocampionaQuantile(file_input, output_name, decrement_ratio_x, decrement_ratio_y, band):
    inputinfo = gdal_array.LoadFile(file_input, band_list=band).astype(float)
    inputinfo[inputinfo <= 0] = np.nan
    minval = np.nanquantile(inputinfo, 0.05)
    maxval = np.nanquantile(inputinfo, 0.95)
    print(f"Image range: {minval} - {maxval}")
    inputinfo = None
    gdal.Translate(output_name, file_input, xRes=decrement_ratio_x, yRes=decrement_ratio_y, outputType=gdal.GDT_Byte,
                   scaleParams=[[minval, maxval, 0, 255]])


@njit()
def proiezioneOmografica(inparray, matrix):
    targetX, targetY = inparray.shape
    outarray = np.zeros((targetX, targetY), dtype=np.float32)

    for (x, y), el in np.ndenumerate(inparray):
        destx = ((matrix[0, 0] * x) + (matrix[0, 1] * y) + (matrix[0, 2])) / (
                (matrix[2, 0] * x) + (matrix[2, 1] * y) + matrix[2, 2])
        desty = ((matrix[1, 0] * x) + (matrix[1, 1] * y) + (matrix[1, 2])) / (
                (matrix[2, 0] * x) + (matrix[2, 1] * y) + matrix[2, 2])

        if (0 <= destx <= targetX) and (0 <= desty <= targetY):
            outarray[int(np.round(destx)), int(np.round(desty))] = el

    return outarray


def good_filter(matches, startingdistanceratio=0.7, maxiteration=4, mingood=100):
    ratio = startingdistanceratio
    iteration = maxiteration

    while True:
        templist = []
        for m, n in matches:
            if m.distance < ratio * n.distance:
                templist.append(m)
        if len(templist) < mingood:
            if iteration > 0:
                iteration -= 1
                ratio += 0.05
            else:
                return templist
        else:
            return templist


def georeference(hiri_input, cas_input, hiri_output, hiri_contrast_thresold, cas_contrast_thresold, minmatch=100,
                 slavepanband=3, panband=3, contrastSlv=False, contrastRef=False):
    timestart = perf_counter()

    temp_Hir_Open = gdal.Open(hiri_input)
    hirProj = temp_Hir_Open.GetProjection()
    hirGeotr = temp_Hir_Open.GetGeoTransform()
    hirBands = temp_Hir_Open.RasterCount
    hirPxColumns = temp_Hir_Open.RasterXSize
    hirPxRows = temp_Hir_Open.RasterYSize
    hir_sizepx_x, hir_sizepx_y = hirGeotr[1], hirGeotr[5]

    bandsNames = []
    for i in range(1, hirBands + 1):
        banddescription = temp_Hir_Open.GetRasterBand(i).GetDescription()
        bandsNames.append(banddescription)

    temp_Hir_Open = None

    temp_Cas_Open = gdal.Open(cas_input)
    casProj = temp_Cas_Open.GetProjection()
    casGeotr = temp_Cas_Open.GetGeoTransform()
    casBands = temp_Cas_Open.RasterCount
    casPxColumns = temp_Cas_Open.RasterXSize
    casPxRows = temp_Cas_Open.RasterYSize
    cas_sizepx_x, cas_sizepx_y = casGeotr[1], casGeotr[5]
    temp_Cas_Open = None

    if os.path.exists("temp") is False:
        os.mkdir("temp")

    hiri_subcamp_temp = "temp/temp_hiri.tif"
    ref_cas = "temp/temp_cas.tif"

    if contrastSlv:
        SottocampionaQuantile(hiri_input, hiri_subcamp_temp, cas_sizepx_x, cas_sizepx_y, band=[slavepanband])
    else:
        SottocampionaMinMax(hiri_input, hiri_subcamp_temp, cas_sizepx_x, cas_sizepx_y, band=[slavepanband])

    if contrastRef:
        SottocampionaQuantile(cas_input, ref_cas, cas_sizepx_x, cas_sizepx_y, band=[panband])
    else:
        SottocampionaMinMax(cas_input, ref_cas, cas_sizepx_x, cas_sizepx_y, band=[panband])

    img1 = gdal_array.LoadFile(hiri_subcamp_temp, band_list=[1])  # queryImage
    plt.imshow(img1), plt.show()
    img2 = gdal_array.LoadFile(ref_cas, band_list=[1])
    plt.imshow(img2), plt.show()

    img1_q01 = np.quantile(img1, 0.001)
    img1_q99 = np.quantile(img1, 0.999)

    img1_mask = np.full_like(img1, 1)
    img1_mask[img1 <= img1_q01] = 0
    img1_mask[img1 >= img1_q99] = 0

    img2_mask = np.full_like(img2, 1)

    img2_q10 = np.quantile(img2, 0.001)
    img2_q90 = np.quantile(img2, 0.999)

    img2_mask[img2 <= img2_q10] = 0
    img2_mask[img2 >= img2_q90] = 0

    plt.imshow(img1_mask), plt.show()
    plt.imshow(img2_mask), plt.show()

    sift1 = cv.SIFT_create(contrastThreshold=hiri_contrast_thresold)
    sift2 = cv.SIFT_create(contrastThreshold=cas_contrast_thresold)

    # find the keypoints and descriptors with SIFT
    kp1, des1 = sift1.detectAndCompute(img1, img1_mask)
    kp2, des2 = sift2.detectAndCompute(img2, img2_mask)

    imgKeyp1 = cv.drawKeypoints(img1, kp1, img1, flags=cv.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
    imgKeyp2 = cv.drawKeypoints(img2, kp2, img2, flags=cv.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
    print(f"SIFT Keypoints HiRISE: {len(kp1)}, CaSSIS {len(kp2)}")
    cv.imwrite('sift_keypoints1.jpg', imgKeyp1)
    cv.imwrite('sift_keypoints2.jpg', imgKeyp2)

    ''' FLANN ALGORITHMS CODES
    FLANN_INDEX_LINEAR = 0
    FLANN_INDEX_KDTREE = 1
    FLANN_INDEX_KMEANS = 2
    FLANN_INDEX_COMPOSITE = 3
    FLANN_INDEX_KDTREE_SINGLE = 4
    FLANN_INDEX_HIERARCHICAL = 5
    FLANN_INDEX_LSH = 6
    FLANN_INDEX_SAVED = 254
    FLANN_INDEX_AUTOTUNED = 255
    '''
    index_params = dict(algorithm=0, trees=5)
    search_params = dict(checks=500)

    # bruteforceMatcher = cv.BFMatcher()
    # matches = bruteforceMatcher.knnMatch(des1, des2, k=2)

    flann = cv.FlannBasedMatcher(index_params, search_params)
    matches = flann.knnMatch(des1, des2, k=2)

    # store all the good matches as per Lowe's ratio test.
    good = good_filter(matches, 0.7, 4, minmatch)

    if len(good) > minmatch:
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
        img_match = cv.drawMatches(img1, kp1, img2_p, kp2, good, None, **draw_params)
        plt.imshow(img_match, 'gray'), plt.show()
        img_match = None

        #### aggiunte da qui
        print("Do you want to try other points in a specific area of the CaSSIS? ")

        if input("Y or [N]") in ["Y", "y", "si", "fuck yeah", "yep"]:
            print("Select limits")
            xmin = int(input("X Min"))
            xmax = int(input("X Max:"))
            ymin = int(input("Y Min"))
            ymax = int(input("Y Max:"))

            mask1 = np.zeros_like(img2)
            mask1[ymin:ymax, xmin:xmax] = 1
            mask1[img2 == 0] = 0
            plt.imshow(mask1), plt.show()

            sift3 = cv.SIFT_create()
            kp3_2, des3_2 = sift3.detectAndCompute(img2, mask1)

            matches_2 = flann.knnMatch(des1, des3_2, k=2)

            good3 = good_filter(matches_2, 0.7, 3, 50)

            src_pts_2 = np.float32([kp1[m.queryIdx].pt for m in good3]).reshape(-1, 1, 2)
            # dst_pts = np.float32([kp3[m.trainIdx].pt for m in good2]).reshape(-1, 1, 2)
            dst_pts_2 = np.float32([kp3_2[m.trainIdx].pt for m in good3]).reshape(-1, 1, 2)

            src_pts_temp = np.concatenate((src_pts, src_pts_2), axis=0)
            dst_pts_temp = np.concatenate((dst_pts, dst_pts_2), axis=0)

            M_temp, mask = cv.findHomography(src_pts_temp, dst_pts_temp, cv.RANSAC, 5.0)
            matchesMask = mask.ravel().tolist()
            h, w = img1.shape
            pts = np.float32([[0, 0], [0, h - 1], [w - 1, h - 1], [w - 1, 0]]).reshape(-1, 1, 2)
            dst = cv.perspectiveTransform(pts, M_temp)

            img2 = cv.polylines(img2, [np.int32(dst)], True, 255, 3, cv.LINE_AA)
            print(f"{len(good) + len(good3)} matches found!")
            print(M)

            draw_params = dict(matchColor=(0, 255, 0),  # draw matches in green color
                               singlePointColor=None,
                               matchesMask=matchesMask,  # draw only inliers
                               flags=2)
            img4 = cv.drawMatches(img1, kp1, img2, kp2, good + good3, None, **draw_params)
            plt.imshow(img4, 'gray'), plt.show()

            print("Do you want to keep the new points?")
            if input("Y or [N]") in ["Y", "y", "si", "fuck yeah", "yep"]:
                M = M_temp

        print("Start projection...")

        hircas_ratio_x = cas_sizepx_x / hir_sizepx_x
        hircas_ratio_y = cas_sizepx_y / hir_sizepx_y

        listlong = np.int32([dst[0, 0, 0], dst[1, 0, 0], dst[2, 0, 0], dst[3, 0, 0]])
        listlat = np.int32([dst[0, 0, 1], dst[1, 0, 1], dst[2, 0, 1], dst[3, 0, 1]])

        left_x, top_y, width, height = np.round(np.min(listlong) * hircas_ratio_x), np.round(
            np.min(listlat) * hircas_ratio_y), \
            np.round((np.max(listlong) - np.min(listlong)) * hircas_ratio_x), \
            np.round((np.max(listlat) - np.min(listlat)) * hircas_ratio_y)

        #scalefactor = createScaleMatrix(hircas_ratio_x, hircas_ratio_y)

        #tform = transform.ProjectiveTransform(matrix= np.dot(scalefactor,M))

        scalingfactormatrix = createScaleMatrix(1 / hircas_ratio_x, 1 / hircas_ratio_y)
        invscaling = np.linalg.inv(scalingfactormatrix)

        M2 = np.dot(invscaling, np.dot(M, scalingfactormatrix))
        #translx = (M2[0, 2] * hircas_ratio_x) - M2[0, 2]
        #transly = (M2[1, 2] * hircas_ratio_y) - M2[1, 2]

        #M2 = np.dot(createTranslationMatrix(translx, transly), M)
        tform = transform.ProjectiveTransform(matrix=M2)
        # tform2 = transform.AffineTransform(matrix=m_translate)
        # tform2 = transform.AffineTransform(translation=(282,1570))

        targetX, targetY = img2.shape
        targetX, targetY = targetX * (hircas_ratio_x), targetY * (hircas_ratio_y)
        targetX, targetY = np.floor(targetX).astype(int), np.floor(targetY).astype(int)
        img1_post = np.zeros((targetX, targetY), dtype=np.uint8)
        img1OrigX, img1OrigY = img1.shape
        if img1OrigX > targetX or img1OrigY > targetY:
            img1_post[0:img1OrigX, 0:img1OrigY] = img1[0:targetX, 0:targetY]
        else:
            img1_post[0:img1OrigX, 0:img1OrigY] = img1
        img1 = None

        targetX, targetY = img1_post.shape
        img1_post = None

        temp_bloccone = "hirblocco.tif"
        inputGeotr_corr = list(casGeotr)
        inputGeotr_corr[1], inputGeotr_corr[5] = hir_sizepx_x, hir_sizepx_y

        driver = gdal.GetDriverByName('GTiff')
        rows, cols = targetX, targetY
        dataset = driver.Create(temp_bloccone, cols, rows, hirBands, gdal.GDT_Float32)
        dataset.SetProjection(casProj)
        dataset.SetGeoTransform(inputGeotr_corr)

        for nband in range(1, hirBands + 1):
            hiri_arr = np.zeros((targetX, targetY), dtype=np.float32)
            rasterArrayHiri = gdal_array.LoadFile(hiri_input, band_list=[nband])
            rasterArrayHiri[rasterArrayHiri < 0] = 65535
            rasAH_OrigX, rasAH_OrigY = rasterArrayHiri.shape
            if (targetX >= rasAH_OrigX) and (targetY >= rasAH_OrigY):
                hiri_arr[0:rasAH_OrigX, 0:rasAH_OrigY] = rasterArrayHiri
            else:
                hiri_arr = rasterArrayHiri
            rasterArrayHiri = None

            tf_img = transform.warp(hiri_arr, tform.inverse, preserve_range=True)
            tf_img[tf_img > 2] = 65535
            dataset.GetRasterBand(nband).WriteArray(tf_img)
            dataset.GetRasterBand(nband).SetDescription(bandsNames[nband - 1])

            dataset.GetRasterBand(nband).SetNoDataValue(65535)
        dataset = None

        gdal.Translate(hiri_output, temp_bloccone, format="GTiff", srcWin=[left_x, top_y, width, height])

        try:
            os.remove(temp_bloccone)
        except:
            sleep(10)
            os.remove(temp_bloccone)
        os.remove(hiri_subcamp_temp)
        if os.path.exists(hiri_subcamp_temp.replace(".tif", ".tif.msk")):
            os.remove(hiri_subcamp_temp.replace(".tif", ".tif.msk"))
        if os.path.exists(hiri_subcamp_temp.replace(".tif", ".tif.aux.xml")):
            os.remove(hiri_subcamp_temp.replace(".tif", ".tif.aux.xml"))
        os.remove(ref_cas)

        print("Creating pyramids")
        buildPyramids(hiri_output)

        # base image
        z = 0.
        index = np.where(np.array(matchesMask) == 1)
        #print(index)
        pixels_base = np.squeeze(src_pts[index, :, :])
        (coord_base_x, coord_base_y, coord_base_z) = applyGeoTransform(casGeotr,
                                                                       pixels_base[:, 0], pixels_base[:, 1], z)
        # wrap image
        z = 0.
        pixels_wrap = np.squeeze(dst_pts[index, :, :])
        (coord_wrap_x, coord_wrap_y, coord_wrap_z) = applyGeoTransform(casGeotr,
                                                                       pixels_wrap[:, 0], pixels_wrap[:, 1], z)

        # Aggregation of transformed coordinates
        my_array = np.vstack((coord_base_x, coord_base_y, coord_wrap_x, coord_wrap_y))
        df_out = pd.DataFrame(np.transpose(my_array), columns=
        ['mapX', 'mapY', 'pixelX', 'pixelY'])
        df_out.head()
        # writing qgis controle point file
        exten = hiri_output[-4:]
        pointpath = hiri_output.replace(exten, "_qgis.points")
        df_out.to_csv(pointpath, sep=',', index=False)
        pointpath = hiri_output.replace(exten, "_arcgis.txt")
        df_out.to_csv(pointpath, sep='	', index=False)


    else:
        print("Not enough matches are found - {}/{}".format(len(good), minmatch))
        matchesMask = None

    print(f"{hiri_output} ready! Time: {perf_counter() - timestart}s")


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    min_matches = 500
    casSlave_contrast_thresold = 0.034
    casMaster_contrast_thresold = 0.06
    cassisSlavepanband = 1
    cassisMasterpanband = 1
    constratStretchSlave = False
    constratStretchMaster = True

    casSlave_input = r"path\MY36_020142_020_3.cubeit.cub"
    casMaster_input = r"path\Mosaic.tif"
    casSlave_output = r"path\Aligned_MY36_020142_020_3.cubeit.cub"

    georeference(casSlave_input, casMaster_input, casSlave_output, casSlave_contrast_thresold,
                 casMaster_contrast_thresold, minmatch=min_matches, slavepanband=cassisSlavepanband,
                 panband=cassisMasterpanband, contrastSlv=constratStretchSlave, contrastRef=constratStretchMaster)
