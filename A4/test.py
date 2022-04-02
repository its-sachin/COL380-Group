from sys import prefix


a=[
    [1,2,3],
    [4,5,6],
    [7,8,9]
]

prefixSum = [0]*9
n=3
for i in range(3):
    for j in range(3):
        prefixSum[i*n+j] = a[i][j]
        if(i == 0):
            if(j!=0):
                prefixSum[j] += prefixSum[j-1]
        else:
            if(j==0):
                prefixSum[i*n+j] += prefixSum[(i-1)*n]
            else:
                prefixSum[i*n+j] += prefixSum[(i-1)*n+j] + prefixSum[i*n + j-1] - prefixSum[(i-1)*n+j-1]

for i in range(3):
    print(prefixSum[i*n:i*n+n])