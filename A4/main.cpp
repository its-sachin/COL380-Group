#include <bits/stdc++.h>



using namespace std;

void readImage(int*** img, int &m, int &n, string fileName){

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
            queryAvg+=queryImg[i][j][3];
        }
    }

    queryAvg/=m*n;

    for(int i=0; i<M-m; i++){
        for(int j=0; j<N-n; j++){
            int sum = 0;
            for(int k=0; k<m; k++){
                for(int l=0; l<n; l++){
                    sum+=dataImg[i+k][j+l][3];
                }
            }
            sum = sum/(m*n);
            if(abs(sum-queryAvg) <= th2){
                // do detailed matching
            }
        }
    }

}