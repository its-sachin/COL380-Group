#include <bits/stdc++.h>
// using namespace std;

__global__
void saxpy(int n, float a, float *x, float *y)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) y[i] = a*x[i] + y[i];
}

void readImage(int* &img, int &m, int &n, std::string fileName){

    std::ifstream file(fileName);
    file >> m >> n;

    img = new int[m*n*4];

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            int sum = 0;
            for (int k = 0; k < 3; k++){
                file>>img[(i*n+j)*4+k];
                sum+=img[(i*n+j)*4+k];
            }
            img[(i*n+j)*4+3] = sum/3;   
        }
    }
    file.close();
}

__device__
int d_floor(float x){
    return floor(x);
}

__device__
int d_ceil(float x){
    return ceil(x);
}

__device__
int d_round(float x){
    return llrint(x);
}

__device__
float getInterpolated(int a, int b, int i, int j, float theta, int M, int N, int* &dataImg, int ind){
    float xx = a + i*cos(theta) - j*sin(theta);
    float yy = b + i*sin(theta) + j*cos(theta);
    float x = xx - d_floor(xx);
    float y = yy - d_floor(yy);
    if(xx<0 || ceil(xx)>=M || yy<0 || ceil(yy)>=N){
        return 0;
    }
    int z00 = dataImg[((d_floor(xx))*N + (d_floor(yy)))*4 + ind];
    int z01 = dataImg[((d_floor(xx))*N + (d_ceil(yy)))*4 + ind];
    int z10 = dataImg[((d_ceil(xx))*N + (d_floor(yy)))*4 + ind];
    int z11 = dataImg[((d_ceil(xx))*N + (d_ceil(yy)))*4 + ind];
    float cx = 1-x;
    // cout << "x: " << x << " y: " << y << " cx: " << cx << " z01: "<< z00 << " " <<z01 << " " <<z10 << " " <<z11<< endl;
    return ( (z00*cx + z10*x)*(1-y) + (z01*cx + z11*x)*y );
}

__global__
void checkGeneral(int * dataImg, int * queryImg, int M, int N, int m, int n, int queryAvg, double th1, double th2, float theta, int* result){

    int a,b;
    int absi = blockIdx.x*256 + threadIdx.x;
    a = absi/N;
    b = absi%N;

    // printf("At start bid: %d tid: %d\n", blockIdx.x, threadIdx.x);

    result[absi] = 1;
    // if(absi > 20000)
    // printf("abs: %d a: %d b: %d\n",absi, a, b);
    float sum = 0;    

    // printf("Before interpol bid: %d tid: %d\n", blockIdx.x, threadIdx.x);

    for(int i =0; i<m; i++){
        for(int j=0; j<n; j++){
            sum += getInterpolated(a,b,i,j,theta,M,N,dataImg,3);
        }
    }


    // cout << "a: " << a << " b: " << b << " " << abs(queryAvg-sum)/(m*n) << endl;
    // printf("After interpol bid: %d tid: %d\n", blockIdx.x, threadIdx.x);


    if(abs(queryAvg-sum)<=th2){
        double sum = 0;
        for (int i = 0; i<m; i++){
            for (int j = 0; j<n; j++){
                for (int r = 0; r < 3; r++){
                    sum+=pow(getInterpolated(a,b,i,j,theta,M,N,dataImg,r)-queryImg[(i*n+j)*4+r],2)/(m*n*3);
                }
            }
        }
        // cout << "   -> " <<sqrt(sum) << endl;
        // printf("    -> %f\n",sqrt(sum));
        if(sqrt(sum)<=th1){
            int ansx = M-d_round(a + m*cos(theta) );
            int ansy = d_round(b + n*sin(theta) );
            int anst = (int)(theta*180/M_PI);
            printf("res = %d %d %d\n",ansx,ansy,anst);
            return;
        }
    }
}

int getAvg(int* &queryImg, int m, int n){
    int queryAvg = 0;

    for(int i=0; i<m; i++){        
        for(int j=0; j<n; j++){
            queryAvg+=queryImg[(i*n+j)*4+3];
        }
    }
    return queryAvg;
}


int main(int argc, char** argv)
{
    std::string dataImgPath = argv[1];
    std::string queryImgPath = argv[2];
    double th1 = std::stod(argv[3]);
    double th2 = std::stod(argv[4]);
    int maxN = std::stoi(argv[5]);

    int M,N,m,n;

    int *dataImg;
    int *queryImg;
    
    readImage(dataImg,M,N,dataImgPath);
    readImage(queryImg,m,n,queryImgPath);

    th2*=m*n;
    
    int *d_dataImg;
    int *d_queryImg;
    
   
    cudaMalloc(&d_dataImg, (M*N*4)*sizeof(int));
    cudaMalloc(&d_queryImg, (m*n*4)*sizeof(int));

    cudaMemcpy(d_dataImg, dataImg, (M*N*4)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_queryImg, queryImg, (m*n*4)*sizeof(int), cudaMemcpyHostToDevice);

    int *result = new int[N*M];
    memset(result, 0, M*N*sizeof(int));
    int *d_result;
    cudaMalloc(&d_result, M*N*sizeof(float));
    cudaMemcpy(d_result, result, M*N*sizeof(float), cudaMemcpyHostToDevice);

    int queryAvg = getAvg(queryImg, m,n);
    // checkGeneral(dataImg, queryImg, M,N,m,n,queryAvg,th1,th2,45*M_PI/180);

    checkGeneral<<<(N*M+255)/256, 256>>>(d_dataImg, d_queryImg, M,N,m,n,queryAvg,th1,th2,45*M_PI/180,d_result);

    cudaMemcpy(result, d_result, M*N*sizeof(int), cudaMemcpyDeviceToHost);

    for(int i=0; i<M; i++){
        for(int j=0; j<N; j++){
            if(result[i*N+j]==1){
                printf("%d %d\n",i,j);
            }
        }
    }

    cudaFree(d_dataImg);
    cudaFree(d_queryImg);
    delete(dataImg);
    delete(queryImg);
}
