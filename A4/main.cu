#include <bits/stdc++.h>
#include <chrono>
// using namespace std;

void readImage(int* &img, int &m, int &n, std::string fileName){

    std::ifstream file(fileName);
    file >> m >> n;

    img = new int[m*n*3];

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            int sum = 0;
            for (int k = 0; k < 3; k++){
                file>>img[(i*n+j)*3+k];
                sum+=img[(i*n+j)*3+k];
            }
        }
    }
    file.close();
}

void readImage(int* &img, int &m, int &n, std::string fileName, float* &prefix){

    std::ifstream file(fileName);
    file >> m >> n;

    img = new int[m*n*3];
    prefix = new float[m*n];

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            float sum = 0;
            for (int k = 0; k < 3; k++){
                file>>img[(i*n+j)*3+k];
                sum+=img[(i*n+j)*3+k];
            }
            prefix[i*n +j]=sum/3;
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

    float xx = a - i*cos(theta) - j*sin(theta);
    float yy = b - i*sin(theta) + j*cos(theta);
    float x = xx - d_floor(xx);
    float y = yy - d_floor(yy);
    if(xx<0 || ceil(xx)>=M || yy<0 || ceil(yy)>=N){
        return 0;
    }
    int z00 = dataImg[((d_floor(xx))*N + (d_floor(yy)))*3 + ind];
    int z01 = dataImg[((d_floor(xx))*N + (d_ceil(yy)))*3 + ind];
    int z10 = dataImg[((d_ceil(xx))*N + (d_floor(yy)))*3 + ind];
    int z11 = dataImg[((d_ceil(xx))*N + (d_ceil(yy)))*3 + ind];
    float cx = 1-x;
    // cout << "x: " << x << " y: " << y << " cx: " << cx << " z01: "<< z00 << " " <<z01 << " " <<z10 << " " <<z11<< endl;
    return ( (z00*cx + z10*x)*(1-y) + (z01*cx + z11*x)*y );
}

__global__
void checkGeneral(int * dataImg, int * queryImg, float * prefix, int M, int N, int m, int n, float queryAvg, double th1, double th2, float pi, float* result){

    int a,b;
    int absi = blockIdx.x*256 + threadIdx.x;
    a = absi/N;
    b = absi%N;

    int angles[3] = {45,0,-45};
    int t = threadIdx.y;

    if(a >= M or a < 0 )return ;
    // printf("At start bid: %d tid: %d\n", blockIdx.x, threadIdx.x);
    // if(absi > 20000)
    // printf("abs: %d a: %d b: %d\n",absi, a, b);

    // for(int t=0; t<3; t++){
        // float sum = 0;    
        // float theta = angles[t]*pi/180;

        // // printf("Before interpol bid: %d tid: %d\n", blockIdx.x, threadIdx.x);

        // for(int i =0; i<m; i++){
        //     for(int j=0; j<n; j++){
        //         sum += getInterpolated(a,b,i,j,theta,M,N,dataImg,3);
        //     }
        // }

        
        float theta = angles[t]*pi/180;
        

        // if(M-a-1 == 840 and b == 900)
        // printf("before a: %d b: %d theta: %f th2: %f, queryAvg: %f\n", a, b, theta,th2,queryAvg);

        int a1,b1,a2,b2;

        if(theta < 0){
            a1 = a + m*sin(theta);
            b2 = b + m*cos(theta) - n*sin(theta);
            a2 = a + n*cos(theta);
            b1 = b;
        }

        else{
            a2 = a;
            b1 = b - m*sin(theta) ;
            b2 = b + n*cos(theta);
            a1 = a - n*cos(theta) - m*sin(theta) ;
        }

        int denom = abs((a2-a1)*(b2-b1));

        a1 = max(min(a1,M-1),0);
        b2 = max(min(b2,N-1),0);
        b1 = max(min(b1,N-1),0);
        a2 = max(min(a2,M-1),0);
    
        float sum = (prefix[a1*N + b1] + prefix[a2*N + b2] - prefix[a2*N + b1] - prefix[a1*N + b2])/denom;

        // //printf("p: %d q: %d r: %d , s: %d \n", p, q, r, s);
        // //printf("a1: %d b1: %d a2: %d b2: %d  P: %d Q: %d R: %d S: %d theta: %f \n ", a1, b1,a2,b2 ,a1*N+b1, a2*N+b2, a1*N+b2, a2*N+b1 ,theta);
        // printf("a1: %d b1: %d a2: %d b2: %d  P: %d Q: %d R: %d S: %d theta: %f \n ", a1, b1,a2,b2 ,p,q,r,s,theta);
        // //printf("a: %d b: %d val: %d\n", a, b, abs(queryAvg-sum));
        // //int sum = prefixSum[a1*N + b1] + prefixSum[a2*N + b2] - prefixSum[a1*N + b2] - prefixSum[a2*N + b1];

        // if(M-a-1 == 840 and b == 900){
        //     printf("(a: %d b: %d angle : %d), (sum: %f) , (queryAvg: %f) , (absDiff: %f), (a1 %d, b1 %d, a2 %d, b2 %d), (p11 %f, p22 %f, p21 %f, p12 %f) \n", a, b,angles[t],sum,queryAvg,abs(queryAvg-sum),a1,b,a2,b2,prefix[a1*N + b],prefix[a2*N + b2],prefix[a2*N + b],prefix[a1*N + b2]);
        // }

        if(abs(queryAvg-sum)<=th2){
            double sum = 0;
            for (int i = 0; i<m; i++){
                for (int j = 0; j<n; j++){
                    for (int r = 0; r < 3; r++){
                        sum+=pow(getInterpolated(a,b,i,j,theta,M,N,dataImg,r)-queryImg[((m-i-1)*n+j)*3+r],2)/(m*n*3);
                    }
                }
            }
            
            
            // cout << "   -> " <<sqrt(sum) << endl;
            // printf("    -> %f\n",sqrt(sum));
            float sq = sqrt(sum);
            if(M-a-1 == 840 and b == 900){
                printf("sqrt %f",sq);
            }
            if(sq<=th1){
                int ansx =  M- a-1 ;
                int ansy =  b ;
                // printf("IRes: %d %d %d %f\n",ansx,ansy,t,sq);
                result[ansx*N*3 + ansy*3 + t] = sq;
                // result[(ansx*N + ansy)*3 = ansx;
                return;
            // }
        }
    }
}

float getAvg(int* &queryImg, int m, int n){
    float queryAvg = 0;

    for(int i=0; i<m; i++){        
        for(int j=0; j<n; j++){
            float sum = 0;
            for(int k=0; k<3; k++)
                sum += queryImg[(i*n+j)*3+k];
            queryAvg += sum/3;
        }
    }
    return queryAvg/(m*n);
}

__global__
void rowsum(float* arr, int m, int n){
    int rownum = blockIdx.x*256 + threadIdx.x;
    if(rownum>=m)return;
    for(int i=1; i<n;i++){
        arr[rownum*n+i]+=arr[rownum*n+i-1];
    }
}

__global__
void colsum(float* arr, int m, int n){
    int colnum = blockIdx.x*256 + threadIdx.x;
    if(colnum>=n)return;
    for(int i=1; i<m;i++){
        arr[i*n+colnum]+=arr[(i-1)*n+colnum];
    }
}


void setSum(int *a, float *&sum, int m, int n){

    sum = new float[m*n];
    sum[0] = (a[0]+a[1]+a[2])/3;
 
    for (int i=1; i<n; i++){
        float s = 0;
        for(int k=0; k<3; k++)s+=a[i*3+k];
        sum[i] = sum[i-1] + s/3; 
    }
    for (int i=1; i<m; i++){
        float s = 0;
        for(int k=0; k<3; k++)s+=a[i*n*3+k];
        sum[i*n] = sum[(i-1)*n] + s/3;
    }
 
    for (int i=1; i<m; i++){
        for (int j=1; j<n; j++){
            float s = 0;
            for(int k=0; k<3; k++)s+=a[i*n*3+j*3+k];
            sum[i*n+j] = sum[(i - 1)*n+j] + sum[i*n + j - 1] - sum[(i - 1)*n + j - 1] + s/3;
        }
    }
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
    float *dataPrefix;

    // std::cout << "1" << std::endl;

    readImage(dataImg,M,N,dataImgPath,dataPrefix);
    // std::cout << "1.5" << std::endl;
    readImage(queryImg,m,n,queryImgPath);
    // std::cout << "1.75" << std::endl;

    
    int *d_dataImg;
    int *d_queryImg;
     
    cudaMalloc(&d_dataImg, (M*N*3)*sizeof(int));
    cudaMalloc(&d_queryImg, (m*n*3)*sizeof(int));

    cudaMemcpy(d_dataImg, dataImg, (M*N*3)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_queryImg, queryImg, (m*n*3)*sizeof(int), cudaMemcpyHostToDevice);

    float *d_dataPrefix;
    cudaMalloc(&d_dataPrefix, (M*N)*sizeof(float));
    cudaMemcpy(d_dataPrefix, dataPrefix, (M*N)*sizeof(float), cudaMemcpyHostToDevice);

    rowsum<<<(M+255)/256,256>>>(d_dataPrefix,M,N);
    cudaMemcpy(dataPrefix, d_dataPrefix, (M*N)*sizeof(float), cudaMemcpyDeviceToHost);
    colsum<<<(N+255)/256,256>>>(d_dataPrefix,M,N);

    cudaMemcpy(dataPrefix, d_dataPrefix, (M*N)*sizeof(float), cudaMemcpyDeviceToHost);

    float* temp;
    setSum(dataImg,temp,M,N);

    // cudaMemcpy(d_dataPrefix, temp,(M*N)*sizeof(float), cudaMemcpyHostToDevice);

    // for(int i=0; i<M;i++){
    //     for(int j=0;j<N;j++)
    //         std::cout<<"("<<temp[i*N+j]<<" " <<dataPrefix[i*N+j] << ")";
    //     std::cout<<std::endl;
    // }

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

    float queryAvg = getAvg(queryImg, m,n);
    // checkGeneral(dataImg, queryImg, M,N,m,n,queryAvg,th1,th2,45*M_PI/180);

    auto mid = std::chrono::high_resolution_clock::now();

    std::cout << "Pre processing Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(mid - start).count() << " ms" << std::endl;

    checkGeneral<<<(N*M+255)/256, {256,3,1}>>>(d_dataImg, d_queryImg, d_dataPrefix, M,N,m,n,queryAvg,th1,th2,M_PI,d_result);

    cudaMemcpy(result, d_result, M*N*3*sizeof(float), cudaMemcpyDeviceToHost);

    std::priority_queue <std::pair<float, container*> > pq;
    int angles[3] = {45,0,-45};
    
    for(int i=0; i<M; i++){
        for(int j=0; j<N; j++){
            for(int k=0; k<3; k++){
                if(result[i*N*3+j*3+k]!=-1){
                    // printf("%d %d %d %f\n",i,j,k,result[i*N*3+j*3+k]);
                    container* c = new container(i,j,angles[k]);
                    pq.push({result[i*N*3+j*3+k], c});
                }
                if(pq.size()>maxN)pq.pop();
            }
        }
    }
    //std::cout<<"pq size: "<<pq.size()<<std::endl;
    int pqSize = pq.size();
    std::vector<std::pair<float, container*>> vecRes;
    for(int i = 0;i<pqSize;i++){
        vecRes.push_back(pq.top());
        pq.pop();
    }
    std::ofstream outfile("output.txt");
    std::cout << std::endl;
    reverse(vecRes.begin(), vecRes.end());
    for(int i = 0;i<vecRes.size();i++){
        outfile << vecRes[i].second->x << " " << vecRes[i].second->y << " " << vecRes[i].second->angle << std::endl;
        printf("Res[%d]: %d %d %d %f\n",i,vecRes[i].second->x,vecRes[i].second->y,vecRes[i].second->angle,vecRes[i].first);
    }

    // for(int i=0; i<maxN && pq.size() > 0; i++){
    //     std::pair<float, container*> p = pq.top();
    //     pq.pop();
    //     printf("Res[%d]: %d %d %d %f\n",i,p.second->x,p.second->y,p.second->angle,p.first);
    // }

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Computation Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - mid).count() << " ms" << std::endl;
    std::cout << "Total Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

    std::cout << "CUDA error: " << cudaGetErrorString(cudaGetLastError()) << std::endl; // add

    // cudaFree(d_dataImg);
    // cudaFree(d_queryImg);
    // cudaFree(d_result);
    // delete(dataImg);
    // delete(result);
    // delete(queryImg);
}

  
// #include <bits/stdc++.h>
// #include <chrono>
// // using namespace std;

// void readImage(int* &img, int &m, int &n, std::string fileName){

//     std::ifstream file(fileName);
//     file >> m >> n;

//     img = new int[m*n*3];

//     for (int i = 0; i < m; i++){
//         for (int j = 0; j < n; j++){
//             int sum = 0;
//             for (int k = 0; k < 3; k++){
//                 file>>img[(i*n+j)*3+k];
//                 sum+=img[(i*n+j)*3+k];
//             }
//         }
//     }
//     file.close();
// }

// __device__
// int d_floor(float x){
//     return floor(x);
// }

// __device__
// int d_ceil(float x){
//     return ceil(x);
// }

// __device__
// int d_round(float x){
//     return llrint(x);
// }
// __device__
// float getInterpolated(int a, int b, int i, int j, float theta, int M, int N, int* &dataImg, int ind){

//     float xx = a - i*cos(theta) - j*sin(theta);
//     float yy = b - i*sin(theta) + j*cos(theta);
//     float x = xx - d_floor(xx);
//     float y = yy - d_floor(yy);
//     if(xx<0 || ceil(xx)>=M || yy<0 || ceil(yy)>=N){
//         return 0;
//     }
//     int z00 = dataImg[((d_floor(xx))*N + (d_floor(yy)))*3 + ind];
//     int z01 = dataImg[((d_floor(xx))*N + (d_ceil(yy)))*3 + ind];
//     int z10 = dataImg[((d_ceil(xx))*N + (d_floor(yy)))*3 + ind];
//     int z11 = dataImg[((d_ceil(xx))*N + (d_ceil(yy)))*3 + ind];
//     float cx = 1-x;
//     // cout << "x: " << x << " y: " << y << " cx: " << cx << " z01: "<< z00 << " " <<z01 << " " <<z10 << " " <<z11<< endl;
//     return ( (z00*cx + z10*x)*(1-y) + (z01*cx + z11*x)*y );
// }

// __global__
// void checkGeneral(int * dataImg, int * queryImg, float * prefix, int M, int N, int m, int n, float queryAvg, double th1, double th2, float pi, float* result){

//     int a,b;
//     int absi = blockIdx.x*256 + threadIdx.x;
//     a = absi/N;
//     b = absi%N;

//     int angles[3] = {45,0,-45};

//     if(a >= M or a < 0 )return ;
//     // printf("At start bid: %d tid: %d\n", blockIdx.x, threadIdx.x);
//     // if(absi > 20000)
//     // printf("abs: %d a: %d b: %d\n",absi, a, b);

//     for(int t=0; t<3; t++){
//         // float sum = 0;    
//         // float theta = angles[t]*pi/180;

//         // // printf("Before interpol bid: %d tid: %d\n", blockIdx.x, threadIdx.x);

//         // for(int i =0; i<m; i++){
//         //     for(int j=0; j<n; j++){
//         //         sum += getInterpolated(a,b,i,j,theta,M,N,dataImg,3);
//         //     }
//         // }

        
//         float theta = angles[t]*pi/180;

//         // if(a == 49 and b == 49)
//         // printf("before a: %d b: %d theta: %f th2: %f, queryAvg: %f\n", a, b, theta,th2,queryAvg);

//         int a1,b1,a2,b2;

//         if(theta < 0){
//             a1 = a + m*sin(theta);
//             b2 = b + m*cos(theta) - n*sin(theta);
//             a2 = a + n*cos(theta);
//             b1 = b;
//         }

//         else{
//             a2 = a;
//             b1 = b - m*sin(theta);
//             b2 = b + n*cos(theta);
//             a1 = a - n*cos(theta) - m*sin(theta);
//         }

//         int denom = abs((a2-a1)*(b2-b1));

//         a1 = max(min(a1,M-1),0);
//         b2 = max(min(b2,N-1),0);
//         b1 = max(min(b1,N-1),0);
//         a2 = max(min(a2,M-1),0);
    
//         float sum = (prefix[a1*N + b1] + prefix[a2*N + b2] - prefix[a2*N + b1] - prefix[a1*N + b2])/denom;

//         // //printf("p: %d q: %d r: %d , s: %d \n", p, q, r, s);
//         // //printf("a1: %d b1: %d a2: %d b2: %d  P: %d Q: %d R: %d S: %d theta: %f \n ", a1, b1,a2,b2 ,a1*N+b1, a2*N+b2, a1*N+b2, a2*N+b1 ,theta);
//         // printf("a1: %d b1: %d a2: %d b2: %d  P: %d Q: %d R: %d S: %d theta: %f \n ", a1, b1,a2,b2 ,p,q,r,s,theta);
//         // //printf("a: %d b: %d val: %d\n", a, b, abs(queryAvg-sum));
//         // //int sum = prefixSum[a1*N + b1] + prefixSum[a2*N + b2] - prefixSum[a1*N + b2] - prefixSum[a2*N + b1];


//         if(abs(queryAvg-sum)<=th2){
//             double sum = 0;
//             for (int i = 0; i<m; i++){
//                 for (int j = 0; j<n; j++){
//                     for (int r = 0; r < 3; r++){
//                         sum+=pow(getInterpolated(a,b,i,j,theta,M,N,dataImg,r)-queryImg[((m-i-1)*n+j)*3+r],2)/(m*n*3);
//                     }
//                 }
//             }
//             // cout << "   -> " <<sqrt(sum) << endl;
//             // printf("    -> %f\n",sqrt(sum));
//             float sq = sqrt(sum);
//             if(sq<=th1){
//                 int ansx =  M- a-1 ;
//                 int ansy =  b ;
//                 // printf("IRes: %d %d %d %f\n",ansx,ansy,t,sq);
//                 result[ansx*N*3 + ansy*3 + t] = sq;
//                 // result[(ansx*N + ansy)*3 = ansx;
//                 return;
//             }
//         }
//     }
// }

// float getAvg(int* &queryImg, int m, int n){
//     float queryAvg = 0;

//     for(int i=0; i<m; i++){        
//         for(int j=0; j<n; j++){
//             float sum = 0;
//             for(int k=0; k<3; k++)
//                 sum += queryImg[(i*n+j)*3+k];
//             queryAvg += sum/3;
//         }
//     }
//     return queryAvg/(m*n);
// }

// class container{

//     public:
//     int x,y,angle;

//     container(int a,int b,int c){
//         x = a;
//         y = b;
//         angle = c;
//     }
// };

// void setSum(int *a, float *&sum, int m, int n){

//     sum = new float[m*n];
//     sum[0] = (a[0]+a[1]+a[2])/3;
 
//     for (int i=1; i<n; i++){
//         float s = 0;
//         for(int k=0; k<3; k++)s+=a[i*3+k];
//         sum[i] = sum[i-1] + s/3; 
//     }
//     for (int i=1; i<m; i++){
//         float s = 0;
//         for(int k=0; k<3; k++)s+=a[i*n*3+k];
//         sum[i*n] = sum[(i-1)*n] + s/3;
//     }
 
//     for (int i=1; i<m; i++){
//         for (int j=1; j<n; j++){
//             float s = 0;
//             for(int k=0; k<3; k++)s+=a[i*n*3+j*3+k];
//             sum[i*n+j] = sum[(i - 1)*n+j] + sum[i*n + j - 1] - sum[(i - 1)*n + j - 1] + s/3;
//         }
//     }
// }


// int main(int argc, char** argv)
// {
//     auto start = std::chrono::high_resolution_clock::now();

//     std::string dataImgPath = argv[1];
//     std::string queryImgPath = argv[2];
//     double th1 = std::stod(argv[3]);
//     double th2 = std::stod(argv[4]);
//     int maxN = std::stoi(argv[5]);

//     int M,N,m,n;

//     int *dataImg;
//     int *queryImg;

//     readImage(dataImg,M,N,dataImgPath);
//     readImage(queryImg,m,n,queryImgPath);
    
//     int *d_dataImg;
//     int *d_queryImg;
     
//     cudaMalloc(&d_dataImg, (M*N*3)*sizeof(int));
//     cudaMalloc(&d_queryImg, (m*n*3)*sizeof(int));

//     cudaMemcpy(d_dataImg, dataImg, (M*N*3)*sizeof(int), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_queryImg, queryImg, (m*n*3)*sizeof(int), cudaMemcpyHostToDevice);

//     float *dataPrefix;
//     float *d_dataPrefix;
//     setSum(dataImg,dataPrefix,M,N);
//     cudaMalloc(&d_dataPrefix, (M*N)*sizeof(float));
//     cudaMemcpy(d_dataPrefix, dataPrefix, (M*N)*sizeof(float), cudaMemcpyHostToDevice);

//     float *result = new float[N*M*3];
//     for(int i=0; i<M; i++){
//         for(int j=0; j<N; j++){
//             for(int k=0; k<3; k++){
//                 result[i*N*3+j*3+k]=-1;
//             }
//         }
//     }
//     float *d_result;
//     cudaMalloc(&d_result, M*N*3*sizeof(float));
//     cudaMemcpy(d_result, result, M*N*3*sizeof(float), cudaMemcpyHostToDevice);

//     float queryAvg = getAvg(queryImg, m,n);
//     // checkGeneral(dataImg, queryImg, M,N,m,n,queryAvg,th1,th2,45*M_PI/180);

//     auto mid = std::chrono::high_resolution_clock::now();

//     std::cout << "Pre processing Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(mid - start).count() << " ms" << std::endl;

//     checkGeneral<<<(N*M+255)/256, 256>>>(d_dataImg, d_queryImg, d_dataPrefix, M,N,m,n,queryAvg,th1,th2,M_PI,d_result);

//     cudaMemcpy(result, d_result, M*N*3*sizeof(float), cudaMemcpyDeviceToHost);

//     std::priority_queue <std::pair<float, container*> > pq;
//     int angles[3] = {45,0,-45};
    
//     for(int i=0; i<M; i++){
//         for(int j=0; j<N; j++){
//             for(int k=0; k<3; k++){
//                 if(result[i*N*3+j*3+k]!=-1){
//                     // printf("%d %d %d %f\n",i,j,k,result[i*N*3+j*3+k]);
//                     container* c = new container(i,j,angles[k]);
//                     pq.push({result[i*N*3+j*3+k], c});
//                 }
//                 if(pq.size()>maxN)pq.pop();
//             }
//         }
//     }
//     //std::cout<<"pq size: "<<pq.size()<<std::endl;
//     int pqSize = pq.size();
//     std::vector<std::pair<float, container*>> vecRes;
//     for(int i = 0;i<pqSize;i++){
//         vecRes.push_back(pq.top());
//         pq.pop();
//     }
//     std::ofstream outfile("output.txt");
//     std::cout << std::endl;
//     reverse(vecRes.begin(), vecRes.end());
//     for(int i = 0;i<vecRes.size();i++){
//         outfile << vecRes[i].second->x << " " << vecRes[i].second->y << " " << vecRes[i].second->angle << std::endl;
//         // printf("Res[%d]: %d %d %d %f\n",i,vecRes[i].second->x,vecRes[i].second->y,vecRes[i].second->angle,vecRes[i].first);
//     }
//     // for(int i=0; i<maxN && pq.size() > 0; i++){
//     //     std::pair<float, container*> p = pq.top();
//     //     pq.pop();
//     //     printf("Res[%d]: %d %d %d %f\n",i,p.second->x,p.second->y,p.second->angle,p.first);
//     // }

//     auto end = std::chrono::high_resolution_clock::now();
//     std::cout << "Computation Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - mid).count() << " ms" << std::endl;
//     std::cout << "Total Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

//     std::cout << "CUDA error: " << cudaGetErrorString(cudaGetLastError()) << std::endl; // add

//     // cudaFree(d_dataImg);
//     // cudaFree(d_queryImg);
//     // cudaFree(d_result);
//     // delete(dataImg);
//     // delete(result);
//     // delete(queryImg);
// }