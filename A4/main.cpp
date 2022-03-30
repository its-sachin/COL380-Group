#include <bits/stdc++.h>



using namespace std;

// int*** initialize3Darray(int*** arr,int m,int n,int s){
//     arr = new int**[m];
//     for (int i = 0; i < m; i++)
//     {
//         arr[i] = new int*[n];
//         for (int j = 0; j < n; j++)
//         {
//             arr[i][j] = new int[s];
//         }
        
//     }
//     return arr;
// }
// int** initialize2Darray(int** arr,int m,int n){
//     arr = new int*[m];
//     for (int i = 0; i < m; i++)
//     {
//         arr[i] = new int[n];
//     }
//     return arr;
// }



template <typename T>
T** initialize2Darray(T** arr,int m,int n) {
    arr = new T*[m];
    for (int i = 0; i < m; i++)
    {
        arr[i] = new T[n];
    }
    return arr;
}

template <typename T>
T*** initialize3Darray(T*** arr,int m,int n,int s){
    arr = new T**[m];
    for (int i = 0; i < m; i++)
    {
        arr[i] = new T*[n];
        for (int j = 0; j < n; j++)
        {
            arr[i][j] = new T[s];
        }
        
    }
    return arr;
}

int main(int argc, char** argv){
    string dataImgPath = argv[1];
    string queryImgPath = argv[2];
    double th1 = stod(argv[3]);
    double th2 = stod(argv[4]);
    int maxN = stoi(argv[5]);

    ifstream datafile; 
    datafile.open(dataImgPath); 

    int mDataImg,nDataImg;   
    datafile>>mDataImg;
    datafile>>nDataImg;


    int ***dataImg;
    int ***queryImg;
    int **dataImgGrey;
    int **queryImgGrey;
    dataImg = initialize3Darray<int>(dataImg,mDataImg,nDataImg,3);
    dataImgGrey = initialize2Darray<int>(dataImgGrey,mDataImg,nDataImg);

    //vector<vector<vector<int>>> dataImg(m, vector<vector<int>>(n, vector<int>(3)));
    //vector<vector<int>> dataImgGrey(m,vector<int>(n));


    for (int i = 0; i < mDataImg; i++)
    {
        for (int j = 0; j < nDataImg; j++)
        {
            int sum = 0;
            for (int k = 0; k < 3; k++)
            {
                datafile>>dataImg[i][j][k];
                sum+=dataImg[i][j][k];
            }
            dataImgGrey[i][j] = sum/3; 
        }
    }

    datafile.close();

    ifstream queryfile; 
    queryfile.open(queryImgPath); 

    int mQueryImg,nQueryImg;   
    queryfile>>mQueryImg;
    queryfile>>nQueryImg;

    queryImg = initialize3Darray<int>(queryImg,mQueryImg,nQueryImg,3);
    queryImgGrey = initialize2Darray<int>(queryImgGrey,mQueryImg,nQueryImg);
    // vector<vector<vector<int>>> queryImg(m, vector<vector<int>>(n, vector<int>(3)));
    // vector<vector<int>> queryImgGrey(m,vector<int>(n));
    for (int i = 0; i < mQueryImg; i++)
    {
        for (int j = 0; j < nQueryImg; j++)
        {
            int sum = 0;
            for (int k = 0; k < 3; k++)
            {
                queryfile>>queryImg[i][j][k];
                sum+=queryImg[i][j][k];
            }
            queryImgGrey[i][j] = sum/3; 
        }
    }

    
    double **dataImgTotalSum;
    dataImgTotalSum = initialize2Darray<double>(dataImgTotalSum,mDataImg,nDataImg);
    //vector<vector<double>> dataImgTotalSum(m,vector<double>(n));
    // for (int i = 0; i <= dataImg.size() - queryImg.size(); i++)
    // {
    //     for (int j = 0; j < dataImg[0].size() - queryImg[0].size(); j++)
    //     {
    //         if(i==0&&j==0){

    //         }
    //     }
        
    // }
    
    


}