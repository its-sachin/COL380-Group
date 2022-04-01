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

dataimg = readImg("data_image.txt")
# queryimg = readImg("query_image.txt")

crop_img = dataimg[0:200, 0:200]
cv2.imshow("dataimage", dataimg)
cv2.imshow('qimg', crop_img)
data = open('data_image1.txt', "w")
data.write(str(200)+" "+str(200)+"\n")
for i in range(200):
    for j in range(200):
        for k in range(3):
            data.write(str(crop_img[i][j][k])+" ")
data.close()
# cv2.imshow("queryimg", queryimg)
cv2.waitKey()

# writeImg('query_image2.jpg')
# q2 = readImg('query_image2.txt')
# cv2.imshow('q2', q2)
# cv2.waitKey()