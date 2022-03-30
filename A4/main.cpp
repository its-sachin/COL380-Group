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
    cout<<sqrt(sum)<<" "<<th1<<endl;
    if(sqrt(sum)<=th1){
        return true;
    }
    return false;
    
}
pair<int,int> templateSearchBasic(int*** &dataImg,int*** &queryImg,double** &dataImgAvg,int M,int N,int m,int n,double th1,double th2,double queryAvg){
    for (int i = 0; i < M-m; i++)
    {
        for (int j = 0; j < N-n; j++)
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


           

    double queryAvg = 0;
    for(int i=0; i<m; i++){        
        for(int j=0; j<n; j++){
            //cout<<"Hekkoi "<<end;
            //cout<<"hello "<<queryImg[i][j][3]<<endl;
            queryAvg+=queryImg[i][j][3];
        }
    }


    queryAvg/=m*n;
    
    double **dataImgTotalSum,**dataImgAvg;
    initialize2Darray<double>(dataImgTotalSum,M,N);
    initialize2Darray<double>(dataImgAvg,M,N);
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
                //cout<<i<<" "<<j<<endl;
                //cout<<dataImgTotalSum[i][j-1]<<endl;
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

    for (int i = 0; i <= M - m; i++)
    {
        for (int j = 0; j < N - n; j++)
        {
            dataImgAvg[i][j] = dataImgTotalSum[i][j]/(m*n);
            //cout<<dataImgAvg[i][j]<<endl;
        }
    }    

    pair<int,int> pos = templateSearchBasic(dataImg,queryImg,dataImgAvg,M,N,m,n,th1,th2,queryAvg);
    cout<<"Res: "<<pos.first<<" "<<pos.second<<endl;
    


}