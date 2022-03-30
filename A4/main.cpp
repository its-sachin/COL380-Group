#include <bits/stdc++.h>



using namespace std;

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


template <typename T>
T** initialize2Darray(T** arr,int m,int n) {
    arr = new T*[m];
    for (int i = 0; i < m; i++)
    {
        arr[i] = new T[n];
    }
    return arr;
}


int main(int argc, char** argv){
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


           

    int queryAvg = 0;
    for(int i=0; i<m; i++){        
        for(int j=0; j<n; j++){
            //cout<<"Hekkoi "<<end;
            //cout<<"hello "<<queryImg[i][j][3]<<endl;
            queryAvg+=queryImg[i][j][3];
        }
    }


    queryAvg/=m*n;
    
    double **dataImgTotalSum;
    dataImgTotalSum = initialize2Darray<double>(dataImgTotalSum,M,N);
    //vector<vector<double>> dataImgTotalSum(m,vector<double>(n));
    for (int i = 0; i <= M - m; i++)
    {
        for (int j = 0; j < N - n; j++)
        {
            if(i==0&&j==0){
                dataImgTotalSum[i][j]  = 0;
                for (int p = i; p < i+m; p++)
                {
                    for (int q = j; q < j+n; q++)
                    {
                        dataImgTotalSum[i][j]+=dataImg[p][q][3];
                    }
                    
                }
            }else if(i==0){
                dataImgTotalSum[i][j]  = dataImgTotalSum[i][j-1];
                for (int p = i; p < i+n; p++)
                {
                    dataImgTotalSum[i][j]-=dataImg[p][j-1][3];
                }
                for (int p = i; p < i+n; p++)
                {
                    dataImgTotalSum[i][j]+=dataImg[p][(j+n-1)-1][3];
                }               
            }else if(j==0){
                dataImgTotalSum[i][j]  = dataImgTotalSum[i-1][j];
                for (int p = j; p < j+m; p++)
                {
                    dataImgTotalSum[i][j]-=dataImg[i-1][p][3];
                }
                for (int p = j; p < j+m; p++)
                {
                    dataImgTotalSum[i][j]+=dataImg[(i+(m))-1][p][3];
                }
            }else{
                dataImgTotalSum[i][j]  = dataImgTotalSum[i-1][j];
                for (int p = j; p < j+m; p++)
                {
                    dataImgTotalSum[i][j]-=dataImg[i-1][p][3];
                }
                for (int p = j; p < j+m; p++)
                {
                    dataImgTotalSum[i][j]+=dataImg[(i+(m))-1][p][3];
                }
            }
        }
        
    }
    
    


}