#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include<bits/stdc++.h>

using namespace std;


int main(int argc, char* argv[]){
    //cout << "hello"<<endl;

    ifstream inFile;
    inFile.open(argv[1], ios::in);
    inFile.seekg(2380970010, std::ios::beg);
    for(int i=0; i<10; i++){
        char a;
        inFile.read((char *)&a, 1);
        cout << a<<endl;
    }
}