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

dataimg = readImg("data_image.txt")
queryimg = readImg("query_image.txt")

cv2.imshow("dataimage", dataimg)
cv2.imshow("queryimg", queryimg)
cv2.waitKey()