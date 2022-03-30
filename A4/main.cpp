#include <bits/stdc++.h>



using namespace std;

int main(int argc, char** argv){
    string dataImgPath = argv[1];
    string queryImgPath = argv[2];
    double th1 = stod(argv[3]);
    double th2 = stod(argv[4]);
    int maxN = stoi(argv[5]);

    ifstream datafile; 
    datafile.open(dataImgPath); 

    int m,n;   
    datafile>>m;
    datafile>>n;

    vector<vector<vector<int>>> dataImg(m, vector<vector<int>>(n, vector<int>(3)));
    
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                datafile>>dataImg[i][j][k];
            }
        }
    }
    datafile.close();

    ifstream queryfile; 
    queryfile.open(queryImgPath); 

    m,n;   
    queryfile>>m;
    queryfile>>n;

    vector<vector<vector<int>>> queryImg(m, vector<vector<int>>(n, vector<int>(3)));
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                queryfile>>queryImg[i][j][k];
            }
        }
    }
    


}