#include <bits/stdc++.h>
using namespace std;

__global__
void saxpy(int n, float a, float *x, float *y)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) y[i] = a*x[i] + y[i];
}

void readImage(int* &img, int &m, int &n, string fileName){

    ifstream file(fileName);
    file >> m >> n;

    img = new int[m*n*4];

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int sum = 0;
            for (int k = 0; k < 3; k++)
            {
                file>>img[(i*n+j)*4+k];
                sum+=img[(i*n+j)*4+k];
            }
            img[(i*n+j)*4+3] = sum/3;   
        }
    }
    file.close();
}


float getInterpolated(int a, int b, int i, int j, float theta, int M, int N, int* &dataImg, int ind){
    float xx = a + i*cos(theta) - j*sin(theta);
    float yy = b + i*sin(theta) + j*cos(theta);
    float x = xx - floor(xx);
    float y = yy - floor(yy);
    if(xx<0 || ceil(xx)>=M || yy<0 || ceil(yy)>=N){
        return 0;
    }
    int z00 = dataImg[(((int)floor(xx))*N + ((int)floor(yy)))*4 + ind];
    int z01 = dataImg[(((int)floor(xx))*N + ((int)ceil(yy)))*4 + ind];
    int z10 = dataImg[(((int)ceil(xx))*N + ((int)floor(yy)))*4 + ind];
    int z11 = dataImg[(((int)ceil(xx))*N + ((int)ceil(yy)))*4 + ind];
    float cx = 1-x;
    // cout << "x: " << x << " y: " << y << " cx: " << cx << " z01: "<< z00 << " " <<z01 << " " <<z10 << " " <<z11<< endl;
    return ( (z00*cx + z10*x)*(1-y) + (z01*cx + z11*x)*y );
}

void checkGeneral(int * &dataImg, int * &queryImg, int M, int N, int m, int n, int queryAvg, double th1, double th2, float theta){

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
                            sum+=pow(getInterpolated(a,b,i,j,theta,M,N,dataImg,r)-queryImg[(i*n+j)*4+r],2)/(m*n*3);
                        }
                    }
                }
                cout << "   -> " <<sqrt(sum) << endl;
                if(sqrt(sum)<=th1){
                    cout<<"Res: "<<M-round(a + m*cos(theta) )<<" "<<round(b + n*sin(theta) )<< " " << (int)(theta*180/M_PI) << endl;
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

    int *dataImg;
    int *queryImg;
    
    readImage(dataImg,M,N,dataImgPath);
    readImage(queryImg,m,n,queryImgPath);

    th2*=m*n;
           

    int queryAvg = 0;

    for(int i=0; i<m; i++){        
        for(int j=0; j<n; j++){
            queryAvg+=queryImg[(i*n+j)*4+3];
        }
    }
    
    int *d_dataImg;
    int* d_queryImg;
   
    cudaMalloc(&d_dataImg, (M*N*4)*sizeof(int));
    cudaMalloc(&d_queryImg, (m*n*4)*sizeof(int));

    cudaMemcpy(d_dataImg, dataImg, (M*N*4)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_queryImg, queryImg, (m*n*4)*sizeof(int), cudaMemcpyHostToDevice);

    cudaFree(d_dataImg);
    cudaFree(d_queryImg);
    free(dataImg);
    free(queryImg);
}