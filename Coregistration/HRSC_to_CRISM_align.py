import numpy as np
import cv2 as cv
from skimage import transform
from osgeo import gdal, gdal_array
import matplotlib.pyplot as plt

def siftOper():
    MIN_MATCH_COUNT = 30
    #img1 = cv.imread('outHir.tif', cv.IMREAD_GRAYSCALE)  # queryImage
    #img2 = cv.imread('outCas.tif', cv.IMREAD_GRAYSCALE)  # trainImage

    mes_NAC_orig = "D:\Sincro\MercuryPansharpening\David Crater\EW0221931523D.cal.map.cub"
    mes_WAC_filter_orig = "D:\Sincro\MercuryPansharpening\David Crater\EW0221931519C.cal.map.cub"

    mes_NAC = "a.tif"
    mes_WAC_filter = "b.tif"

    rasterArray = gdal_array.LoadFile(mes_NAC_orig)
    minArr = np.nanmin(rasterArray)
    maxArr = np.nanmax(rasterArray)
    gdal.Translate(mes_NAC, mes_NAC_orig, scaleParams = [[0, maxArr, 0, 65535]], format="GTiff", outputType=gdal.gdalconst.GDT_UInt16)

    rasterArray = gdal_array.LoadFile(mes_WAC_filter_orig)
    minArr = np.nanmin(rasterArray)
    maxArr = np.nanmax(rasterArray)
    gdal.Translate(mes_WAC_filter, mes_WAC_filter_orig, scaleParams = [[0, maxArr, 0, 65535]], format="GTiff", outputType=gdal.gdalconst.GDT_UInt16)


    NAC_px_size = 825.8435538871799508
    WAC_px_size = 825.8435538871799508

    filter_Output = f"{mes_WAC_filter.split('.')[0]}.tif"
    # Original image scale range


    casHiri_ratio = WAC_px_size / NAC_px_size

    img1 = cv.imread(mes_WAC_filter, cv.IMREAD_GRAYSCALE)  # queryImage
    img2 = cv.imread(mes_NAC, cv.IMREAD_GRAYSCALE)  # trainImage


    #img1 = cv.imread('outCas.tif')  # queryImage
    #img2 = cv.imread('outHir_sot.tif')  # trainImage

    # Initiate SIFT detector
    #sift = cv.SIFT_create(0, 5, 0.15, 10, 3)
    sift1 = cv.SIFT_create(contrastThreshold=0.01)
    sift2 = cv.SIFT_create(contrastThreshold=0.01)
    #kaze = cv.KAZE.create(True)

    #orb = cv.ORB_create()

    # find the keypoints and descriptors with SIFT
    kp1, des1 = sift1.detectAndCompute(img1, None)
    kp2, des2 = sift2.detectAndCompute(img2, None)

    #kp1, des1 = kaze.detectAndCompute(img1, None)
    #kp2, des2 = kaze.detectAndCompute(img2, None)

    #kp1, des1 = orb.detectAndCompute(img1, None)
    #kp2, des2 = orb.detectAndCompute(img2, None)

    imgKeyp1 = cv.drawKeypoints(img1, kp1, img1, flags=cv.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
    imgKeyp2 = cv.drawKeypoints(img2, kp2, img2, flags=cv.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
    print(len(kp1))
    print(len(kp2))
    cv.imwrite('sift_keypoints1.jpg', imgKeyp1)
    cv.imwrite('sift_keypoints2.jpg', imgKeyp2)

    FLANN_INDEX_LINEAR = 0
    FLANN_INDEX_KDTREE = 1
    FLANN_INDEX_KMEANS = 2
    FLANN_INDEX_COMPOSITE = 3
    FLANN_INDEX_KDTREE_SINGLE = 4
    FLANN_INDEX_HIERARCHICAL = 5
    FLANN_INDEX_LSH = 6
    FLANN_INDEX_SAVED = 254
    FLANN_INDEX_AUTOTUNED = 255

    index_params = dict(algorithm=0, trees=5)
    search_params = dict(checks=500)
    flann = cv.FlannBasedMatcher(index_params, search_params)
    matches = flann.knnMatch(des1, des2, k=2)
    # store all the good matches as per Lowe's ratio test.
    good = []
    for m, n in matches:
        if m.distance < 0.7 * n.distance:
            good.append(m)


    if len(good) > MIN_MATCH_COUNT:
        src_pts = np.float32([kp1[m.queryIdx].pt for m in good]).reshape(-1, 1, 2)
        dst_pts = np.float32([kp2[m.trainIdx].pt for m in good]).reshape(-1, 1, 2)
        M, mask = cv.findHomography(src_pts, dst_pts, cv.RANSAC, 5.0)
        matchesMask = mask.ravel().tolist()
        h, w = img1.shape
        pts = np.float32([[0, 0], [0, h - 1], [w - 1, h - 1], [w - 1, 0]]).reshape(-1, 1, 2)
        dst = cv.perspectiveTransform(pts, M)
        img2 = cv.polylines(img2, [np.int32(dst)], True, 255, 3, cv.LINE_AA)
        print(f"{len(good)} matches found!")
        print(M)

        # print(M)

        draw_params = dict(matchColor=(0, 255, 0),  # draw matches in green color
                           singlePointColor=None,
                           matchesMask=matchesMask,  # draw only inliers
                           flags=2)
        img3 = cv.drawMatches(img1, kp1, img2, kp2, good, None, **draw_params)
        plt.imshow(img3, 'gray'), plt.show()
        img3 = None

        M2 = M * casHiri_ratio
        M2[2, 2] = 1.0
        tform = transform.ProjectiveTransform(matrix=M2)

        img1_orig_array = gdal_array.LoadFile(mes_WAC_filter_orig, band_list=[1])

        targetX, targetY = img2.shape
        targetX, targetY = targetX * casHiri_ratio, targetY * casHiri_ratio
        img1_post = np.zeros((np.floor(targetX).astype(int), np.floor(targetY).astype(int)), dtype=np.float32)
        img1OrigX, img1OrigY = img1_orig_array.shape
        img1_post[0:img1OrigX, 0:img1OrigY] = img1_orig_array
        img1 = None
        img1_orig_array = None


        # Riproiezione
        tf_img2 = transform.warp(img1_post, tform.inverse, preserve_range=True)

        driver = gdal.GetDriverByName('GTiff')
        rows, cols = tf_img2.shape
        dataset = driver.Create(filter_Output, cols, rows, 1,
                                gdal.GDT_Float32)
        dataset.GetRasterBand(1).WriteArray(tf_img2)

        temp_cas_Open = gdal.Open(mes_WAC_filter)
        out_Open = gdal.Open(filter_Output)

        gdal_array.CopyDatasetInfo(temp_cas_Open, out_Open)

        dataset = None
        img2 = None
        kp1 = None
        kp2 = None
        tf_img2 = None



    else:
        print("Not enough matches are found - {}/{}".format(len(good), MIN_MATCH_COUNT))
        matchesMask = None




# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    siftOper()
