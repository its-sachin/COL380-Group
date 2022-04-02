#include <bits/stdc++.h>
#include <chrono>
// using namespace std;

void readImage(int* &img, int &m, int &n, std::string fileName, int *&prefixSum, bool isData = false){

    std::ifstream file(fileName);
    file >> m >> n;

    img = new int[m*n*4];

    if(isData){
        prefixSum = new int[m*n];
    }

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            int sum = 0;
            for (int k = 0; k < 3; k++){
                file>>img[(i*n+j)*4+k];
                sum+=img[(i*n+j)*4+k];
            }
            img[(i*n+j)*4+3] = sum/3;  
            if(isData){
                prefixSum[i*n+j] = sum/3;
                if(i == 0){
                    if(j!=0){
                        prefixSum[j] += prefixSum[j-1];
                    }
                }
                else{
                    if(j==0)
                        prefixSum[i*n+j] += prefixSum[(i-1)*n];
                    else
                        prefixSum[i*n+j] += prefixSum[(i-1)*n+j] + prefixSum[i*n + j-1] - prefixSum[(i-1)*n+j-1];
                }
            } 
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
void checkGeneral(int * dataImg, int * queryImg, int * prefixSum, int M, int N, int m, int n, int queryAvg, double th1, double th2, float pi, float* result){

    int a,b;
    int absi = blockIdx.x*256 + threadIdx.x;
    a = absi/N;
    b = absi%N;

    int angles[3] = {45,0,-45};
    // printf("At start bid: %d tid: %d\n", blockIdx.x, threadIdx.x);
    // if(absi > 20000)
    // printf("abs: %d a: %d b: %d\n",absi, a, b);

    for(int t=0; t<3; t++){
        float theta = angles[t]*pi/180;

        int a1 = a - n*sin(theta);
        int b1 = b;
        int b2 = b + n*cos(theta) + m*sin(theta);
        int a2 = a + m*cos(theta);
        int sum = prefixSum[a1*N + b1] + prefixSum[a2*N + b2] - prefixSum[a1*N + b2] - prefixSum[a2*N + b1];

        printf("a: %d b: %d val: %d\n", a, b, abs(queryAvg-sum));
        

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
            float sq = sqrt(sum);
            if(sq<=th1){
                int ansx = M-d_round(a + m*cos(theta) );
                int ansy = d_round(b + m*sin(theta) );
                // printf("IRes: %d %d %d %f\n",ansx,ansy,t,sq);
                result[ansx*N*3 + ansy*3 + t] = sq;
                // result[(ansx*N + ansy)*3 = ansx;
                return;
            }
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

class container{

    public:
    int x,y,angle;

    container(int a,int b,int c){
        x = a;
        y = b;
        angle = c;
    }
};


int main(int argc, char** argv)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::string dataImgPath = argv[1];
    std::string queryImgPath = argv[2];
    double th1 = std::stod(argv[3]);
    double th2 = std::stod(argv[4]);
    int maxN = std::stoi(argv[5]);

    int M,N,m,n;

    int *dataImg;
    int *queryImg;
    int *dataPrefix;

    readImage(dataImg,M,N,dataImgPath,dataPrefix,true);
    readImage(queryImg,m,n,queryImgPath,dataPrefix);

    th2*=m*n;
    
    int *d_dataImg;
    int *d_queryImg;
    int *d_dataPrefix;
    
   
    cudaMalloc(&d_dataImg, (M*N*4)*sizeof(int));
    cudaMalloc(&d_queryImg, (m*n*4)*sizeof(int));
    cudaMalloc(&d_dataPrefix, (M*N)*sizeof(int));

    cudaMemcpy(d_dataImg, dataImg, (M*N*4)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_queryImg, queryImg, (m*n*4)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dataPrefix, dataPrefix, (M*N)*sizeof(int), cudaMemcpyHostToDevice);

    float *result = new float[N*M*3];
    for(int i=0; i<M; i++){
        for(int j=0; j<N; j++){
            for(int k=0; k<3; k++){
                result[i*N*3+j*3+k]=-1;
            }
        }
    }
    float *d_result;
    cudaMalloc(&d_result, M*N*3*sizeof(float));
    cudaMemcpy(d_result, result, M*N*3*sizeof(float), cudaMemcpyHostToDevice);

    int queryAvg = getAvg(queryImg, m,n);
    // checkGeneral(dataImg, queryImg, M,N,m,n,queryAvg,th1,th2,45*M_PI/180);

    auto mid = std::chrono::high_resolution_clock::now();

    std::cout << "Pre processing Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(mid - start).count() << " ms" << std::endl;

    checkGeneral<<<(N*M+255)/256, 256>>>(d_dataImg, d_queryImg, d_dataPrefix, M,N,m,n,queryAvg,th1,th2,M_PI,d_result);

    cudaMemcpy(result, d_result, M*N*3*sizeof(float), cudaMemcpyDeviceToHost);

    std::priority_queue <std::pair<float, container*> > pq;
    int angles[3] = {45,0,-45};
    
    for(int i=0; i<M; i++){
        for(int j=0; j<N; j++){
            for(int k=0; k<3; k++){
                if(result[i*N*3+j*3+k]!=-1){
                    printf("%d %d %d %f\n",i,j,k,result[i*N*3+j*3+k]);
                    container* c = new container(i,j,angles[k]);
                    pq.push({result[i*N*3+j*3+k], c});
                }
            }
        }
    }

    std::cout << std::endl;
    for(int i=0; i<maxN && pq.size() > 0; i++){
        std::pair<float, container*> p = pq.top();
        pq.pop();
        printf("Res[%d]: %d %d %d %f\n",i,p.second->x,p.second->y,p.second->angle,p.first);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Computation Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - mid).count() << " ms" << std::endl;
    std::cout << "Total Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

    std::cout << "CUDA error: " << cudaGetErrorString(cudaGetLastError()) << std::endl; // add

    cudaFree(d_dataImg);
    cudaFree(d_queryImg);
    cudaFree(d_result);
    delete(dataImg);
    delete(result);
    delete(queryImg);
}
