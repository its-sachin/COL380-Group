#include <bits/stdc++.h>
using namespace std;

__global__
void saxpy(int n, float a, float *x, float *y)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) y[i] = a*x[i] + y[i];
}

void readImage(int*** &img, int &m, int &n, string fileName){

    ifstream file(fileName);
    file >> m >> n;

    img = new int**[m];

    for (int i = 0; i < m; i++)
    {
        img[i] = new int*[n];
        for (int j = 0; j < n; j++)
        {
            img[i][j] = new int[4];
            int sum = 0;
            for (int k = 0; k < 3; k++)
            {
                file>>img[i][j][k];
                sum+=img[i][j][k];
            }
            img[i][j][3] = sum/3;   
        }
    }
    file.close();
}



bool RMSD(int*** &dataImg,int*** &queryImg,int i,int j,int M,int N,int m,int n,double th1){
    double sum = 0;
    for (int p = i; p < i+m; p++)
    {
        for (int q = j; q < j+n; q++)
        {
            for (int r = 0; r < 3; r++)
            {
                sum+=pow(dataImg[p][q][r]-queryImg[p-i][q-j][r],2)/(m*n*3);
            }
            
        }
        
    }
    cout<<i<<" " <<j<<" "<<sqrt(sum)<<" "<<th1<<endl;
    if(sqrt(sum)<=th1){
        return true;
    }
    return false;
    
}
pair<int,int> templateSearchBasic(int*** &dataImg,int*** &queryImg,int** &dataImgAvg,int M,int N,int m,int n,double th1,double th2,double queryAvg){
    for (int i = 0; i <= M-m; i++)
    {
        for (int j = 0; j <= N-n; j++)
        {
            if(abs(queryAvg-dataImgAvg[i][j])<=th2){      
                //cout<<"hi "<<endl; 
                if(RMSD(dataImg,queryImg,i,j,M,N,m,n,th1)){
                    return make_pair(i,j);
                }
            }
        }
        
    }
    return make_pair(-1,-1);
}


template <typename T>
void initialize2Darray(T** &arr,int m,int n) {
    arr = new T*[m];
    for (int i = 0; i < m; i++)
    {
        arr[i] = new T[n];
    }
}

void checkZero(int *** &dataImg, int *** &queryImg, int M, int N, int m, int n, int queryAvg, double th1, double th2){
    int **dataImgTotalSum,**dataImgAvg;
    initialize2Darray<int>(dataImgTotalSum,M,N);
    initialize2Darray<int>(dataImgAvg,M,N);

    for (int i = 0; i <= M - m; i++)
    {
        for (int j = 0; j <= N - n; j++)
        {
            dataImgTotalSum[i][j] = 0;
            if(i == 0){
                if(j == 0){
                    for(int a = 0; a < m; a++){
                        for(int b = 0; b < n; b++){
                            dataImgTotalSum[i][j]+=dataImg[i+a][j+b][3];
                        }
                    }
                }
                else{
                    dataImgTotalSum[i][j] = dataImgTotalSum[i][j-1];
                    for(int a=0; a<m; a++){
                        dataImgTotalSum[i][j]+= dataImg[i+a][j+n-1][3] - dataImg[i+a][j-1][3];
                    }
                }
            }
            else{
                dataImgTotalSum[i][j] = dataImgTotalSum[i-1][j];
                for(int a=0; a<n; a++){
                    dataImgTotalSum[i][j]+= dataImg[i+m-1][j+a][3] - dataImg[i-1][j+a][3];
                }
            }
        }
    }
    pair<int,int> pos = templateSearchBasic(dataImg,queryImg,dataImgTotalSum,M,N,m,n,th1,th2,queryAvg);
    cout<<"Res: "<<pos.first<<" "<<pos.second<<endl;
}

float getInterpolated(int a, int b, int i, int j, float theta, int M, int N, int*** &dataImg, int ind){
    float xx = a + i*cos(theta) - j*sin(theta);
    float yy = b + i*sin(theta) + j*cos(theta);
    float x = xx - floor(xx);
    float y = yy - floor(yy);
    if(xx<0 || ceil(xx)>=M || yy<0 || ceil(yy)>=N){
        return 0;
    }
    int z00 = dataImg[(int)floor(xx)][(int)floor(yy)][ind];
    int z01 = dataImg[(int)floor(xx)][(int)ceil(yy)][ind];
    int z10 = dataImg[(int)ceil(xx)][(int)floor(yy)][ind];
    int z11 = dataImg[(int)ceil(xx)][(int)ceil(yy)][ind];
    float cx = 1-x;
    // cout << "x: " << x << " y: " << y << " cx: " << cx << " z01: "<< z00 << " " <<z01 << " " <<z10 << " " <<z11<< endl;
    return ( (z00*cx + z10*x)*(1-y) + (z01*cx + z11*x)*y );
}

void checkGeneral(int *** &dataImg, int *** &queryImg, int M, int N, int m, int n, int queryAvg, double th1, double th2, float theta){

    for(int a=48; a<52; a++){
        for(int b=48; b<52;b ++){
            float sum = 0;
            for(int i =0; i<m; i++){
                for(int j=0; j<n; j++){
                    // cout << sum << endl;
                    sum += getInterpolated(a,b,i,j,theta,M,N,dataImg,3);
                }
            }
            cout << "a: " << a << " b: " << b << " " << abs(queryAvg-sum)/(m*n) << endl;
            if(abs(queryAvg-sum)<=th2){
                double sum = 0;
                for (int i = 0; i<m; i++){
                    for (int j = 0; j<n; j++){
                        for (int r = 0; r < 3; r++){
                            sum+=pow(getInterpolated(a,b,i,j,theta,M,N,dataImg,r)-queryImg[i][j][r],2)/(m*n*3);
                        }
                    }
                }
                cout << "   -> " <<sqrt(sum) << endl;
                if(sqrt(sum)<=th1){
                    cout<<"Res: "<<a<<" "<<b<< " " << (int)(theta*180/M_PI) << endl;
                    return;
                }
            }
        }
    }
}



int main(int argc, char** argv)
{
    string dataImgPath = argv[1];
    string queryImgPath = argv[2];
    double th1 = stod(argv[3]);
    double th2 = stod(argv[4]);
    int maxN = stoi(argv[5]);

    int M,N,m,n;

    int ***dataImg;
    int ***queryImg;
    
    readImage(dataImg,M,N,dataImgPath);
    readImage(queryImg,m,n,queryImgPath);

    th2*=m*n;
           

    int queryAvg = 0;

    for(int i=0; i<m; i++){        
        for(int j=0; j<n; j++){
            queryAvg+=queryImg[i][j][3];
        }
    }
    
    int ***d_dataImg,**d_queryImg;
    //float *x, *y, *d_x, *d_y;
    //x = (float*)malloc(N*sizeof(float));
    //y = (float*)malloc(N*sizeof(float));

    //cudaMalloc(&d_x, N*sizeof(float)); 
    //cudaMalloc(&d_y, N*sizeof(float));
    cudaMalloc((void **) &d_dataImg, (M*N*4)*sizeof(int));
    cudaMalloc((void **) &d_queryImg, (m*n*4)*sizeof(int));


    // for (int i = 0; i < N; i++) {
    //     x[i] = 1.0f;
    //     y[i] = 2.0f;
    // }

    cudaMemcpy(d_dataImg, dataImg, (M*N*4)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_queryImg, queryImg, (m*n*4)*sizeof(int), cudaMemcpyHostToDevice);

    // Perform SAXPY on 1M elements
    // saxpy<<<(N+255)/256, 256>>>(N, 2.0f, d_x, d_y);

    // cudaMemcpy(y, d_y, N*sizeof(float), cudaMemcpyDeviceToHost);

    // float maxError = 0.0f;
    // for (int i = 0; i < N; i++)
    //     maxError = max(maxError, abs(y[i]-4.0f));
    // printf("Max error: %f\n", maxError);

    cudaFree(d_dataImg);
    cudaFree(d_queryImg);
    // free(x);
    // free(y);
}