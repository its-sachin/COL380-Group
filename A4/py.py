from math import sin, cos, pi, ceil, floor,sqrt
import numpy as np
import cv2
from sys import argv

def getInterpolated( a,  b,  i,  j, theta,  M,  N, dataImg,  ind, top=False):
    xx = a + i*cos(theta) - j*sin(theta)
    yy = b + i*sin(theta) + j*cos(theta)
    x = xx - floor(xx)
    y = yy - floor(yy) 
    if(xx<0 or ceil(xx)>=M or yy<0 or ceil(yy)>=N):
        return xx,yy,0
    z00 = dataImg[floor(xx)][floor(yy)][ind]
    z01 = dataImg[floor(xx)][ceil(yy)][ind]
    z10 = dataImg[ceil(xx)][floor(yy)][ind]
    z11 = dataImg[ceil(xx)][ceil(yy)][ind]
    cx = 1-x
    return xx,yy,( (z00*cx + z10*x)*(1-y) + (z01*cx + z11*x)*y )

def setVal(arr, i, j, m, n):
    if(i<0 or j<0 or i>=m or j>= n):return 0
    return arr[i][j]

def checkGeneral( dataImg, queryImg, prefixSum,  M,  N,  m,  n,  queryAvg,  th1,  th2, datacv, querycv):

    queryAvg /= m*n
    angles = [45,0,-45]
    for an in range(3):
        theta = angles[an]*pi/180
        for a in range(49,50):
            for b in range(49,50):
                
                aa1 = int(a - n*sin(theta))
                bb2 = int(b + n*cos(theta) + m*sin(theta))
                aa2 = int(a + m*cos(theta))
                bb1 = b

                a1 = max(min(aa1,M-1),0)
                b1 = b
                b2 = max(min(bb2,N-1),0)
                a2 = max(min(aa2,M-1),0)

                sum = (prefix[a1][b1] + prefix[a2][b2] - prefix[a2][b1] - prefix[a1][b2])/(abs(aa2-aa1)*abs(bb2-bb1))

                # datacv = cv2.circle(datacv, (b1,a1), 2, 255, 2)
                # datacv = cv2.circle(datacv, (b1,a2), 2, 255, 2)
                # datacv = cv2.circle(datacv, (b2,a1), 2, 255, 2)
                # datacv = cv2.circle(datacv, (b2,a2), 2, 255, 2)
                # datacv = cv2.circle(datacv, (b,a), 2, 0, 2)
                # # querycv = cv2.circle(querycv, (int(j),int(i)), 2, 255, 2)
                # cv2.imshow('datacv',datacv)
                # # cv2.imshow('querycv',querycv)
                # cv2.waitKey()


                print("(a: ", a," b: ",b,"angle:",angles[an],") (queryavg: ",queryAvg,")(sum:",sum,") (absdiff : ",abs(queryAvg-sum),")",)
                print("a1",a1,"b1",b1,"a2",a2,"b2",b2)
                print("prefix[a1][b1]: ",prefix[a1][b1],"prefix[a2][b2]: ",prefix[a2][b2],"prefix[a2][b1]: ",prefix[a2][b1],"prefix[a1][b2]: ",prefix[a1][b2])
                if(abs(queryAvg-sum)<=th2):
                    sum = 0
                    for i in range(m):
                        for j in range(n):
                            for r in range(3):
                                x,y,v = getInterpolated(a,b,i,j,theta,M,N,dataImg,r,True)
                                # print('x:',x,'y:',y,'interp : ' ,v, 'query : ', queryImg[i][j][r], 'data:',dataImg[int(x)][int(y)][r])
                                sum+=pow(v-queryImg[i][j][r],2)/(m*n*3)
                            # datacv = cv2.circle(datacv, (int(y),int(x)), 2, 255, 2)
                            # querycv = cv2.circle(querycv, (int(j),int(i)), 2, 255, 2)
                            # cv2.imshow('datacv',datacv)
                            # cv2.imshow('querycv',querycv)
                            # cv2.waitKey()
                    print("   -> -> -> -> -> -> -> ",a,b ,sqrt(sum))
                    if(sqrt(sum)<=th1):
                        print("Res: ",M-round(a+m*cos(theta)),round(b+m*sin(theta)))
                        return
                

def readImg(path):
    img = None
    nump = None
    data = open(path, "r")
    for l in data.readlines():
        line = l.strip('\n').split(' ')
        if(img is None):
            m = int(line[0])
            n = int(line[1])
            img = [[[0,0,0,0] for i in range(n)] for j in range(m)]
            nump = np.zeros((m, n, 3), dtype=np.uint8)
        else:
            for i in range(m):
                for j in range(n):
                    s = 0
                    for k in range(3):
                        img[i][j][k] = int(line[3*(i*n+j)+k])
                        nump[i][j][k] = int(line[3*(i*n+j)+k])
                        s += img[i][j][k]
                    img[i][j][3] = s/3
    return nump,img,m,n

def prefixSum2D(a) :
    R = len(a)
    C = len(a[0])

    psa = [[0 for x in range(C)]
              for y in range(R)]
    psa[0][0] = a[0][0][3]
 
    for i in range(1, C) :
        psa[0][i] = (psa[0][i - 1] + a[0][i][3])
    for i in range(0, R) :
        psa[i][0] = (psa[i - 1][0] + a[i][0][3])
 
    for i in range(1, R) :
        for j in range(1, C) :
            psa[i][j] = (psa[i - 1][j] +
                         psa[i][j - 1] -
                         psa[i - 1][j - 1] +
                           a[i][j][3])
    return psa

datacv,dataImg,M,N = readImg("data_image.txt")
querycv,queryImg,m,n = readImg("query_image.txt")

prefix = prefixSum2D(dataImg)
th1 = 10
th2 = 0.5
queryAvg = 0

for i in range(m):
    for j in range(n):
        queryAvg+=queryImg[i][j][3]

checkGeneral(dataImg,queryImg,prefix,M,N,m,n,queryAvg,th1,th2,datacv,querycv)