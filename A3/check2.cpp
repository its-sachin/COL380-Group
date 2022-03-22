#include<bits/stdc++.h>

using namespace std;


int main(int argc, char* argv[]){
    std::fstream fs("vect.bin", std::ios::in | std::ios::binary);
    for(int i=0; i<6; i++){
        for(int j=0; j<5; j++){
            double n;
            fs.read((char *)&n, sizeof(double));
            cout << n << " ";
        }
        cout << endl;
    }

    // for(int i=0; i<3; i++){
    //         double n;
    //         fs.read((char *)&n, sizeof(double));
    //         cout << n << " \n";
    // }
}