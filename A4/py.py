from math import sin, cos, pi, ceil, floor,sqrt
import numpy as np
import cv2

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

def checkGeneral( dataImg, queryImg,  M,  N,  m,  n,  queryAvg,  th1,  th2, theta, datacv, querycv):

    for a in range(47,55):
        for b in range(47,55):
            
            sum = 0
            for i in range(m):
                for j in range(n):
                    s = getInterpolated(a,b,i,j,theta,M,N,dataImg,3)[2]
                    sum += s
                    
            print("a: ", a," b: ",b,abs(queryAvg-sum)/(m*n))
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
                    print("Res: ",a,b)
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

datacv,dataImg,M,N = readImg("data_image.txt")
querycv,queryImg,m,n = readImg("query_image.txt")

th1 = 10
th2 = 0.5*m*n
queryAvg = 0

for i in range(m):
    for j in range(n):
        queryAvg+=queryImg[i][j][3]

checkGeneral(dataImg,queryImg,M,N,m,n,queryAvg,th1,th2,45*pi/180,datacv,querycv)