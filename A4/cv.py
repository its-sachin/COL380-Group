import cv2
import numpy as np

def readImg(path):
    img = None
    data = open(path, "r")
    for l in data.readlines():
        line = l.strip('\n').split(' ')
        if(img is None):
            m = int(line[0])
            n = int(line[1])
            img = np.zeros((m, n, 3), dtype=np.uint8)
        else:
            for i in range(m):
                for j in range(n):
                    for k in range(3):
                        img[i][j][k] = int(line[3*(i*n+j)+k])
    return img

def writeImg(path):
    img = cv2.imread(path)
    m, n, _ = img.shape
    data = open(path[:-3] + "txt", "w")
    data.write(str(m)+" "+str(n)+"\n")
    for i in range(m):
        for j in range(n):
            for k in range(3):
                data.write(str(img[i][j][k])+" ")
    data.close()

#writeImg("test.jpg")
dataimg = readImg("test_case_1_small_image/data_image_a.txt")
# dataimg = cv2.imread("test.jpg", cv2.IMREAD_COLOR)
# writeImg("query_rotate.jpg")
queryimg = readImg("test_case_1_small_image/query_image_a.txt")

dataimg = cv2.circle(dataimg, (200,199), 2, 255, 2)
cv2.imshow("data", dataimg)
cv2.imshow("query", queryimg)

cv2.waitKey()

# w0 = 100
# h0 = 0
# w1 = 180
# h1 = 250
# crop_img = dataimg[w0:w1, h0:h1]
# # cv2.imshow("test", dataimg)
# # cv2.imshow('qimg', crop_img)
# # data = open('test2.txt', "w")
# # data.write(str(w1-w0)+" "+str(h1-h0)+"\n")
# # for i in range(w1-w0):
# #     for j in range(h1-h0):
# #         for k in range(3):
# #             data.write(str(crop_img[i][j][k])+" ")
# # data.close()
# # # cv2.imshow("queryimg", queryimg)
# # cv2.waitKey()
# cv2.imwrite("test2.jpg",crop_img)
# writeImg('test2.jpg')
# # q2 = readImg('test2.txt')
# # cv2.imshow('q2', q2)
# # cv2.waitKey()











