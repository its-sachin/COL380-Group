#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include<bits/stdc++.h>

using namespace std;
    
int main(int argc, char **argv) {
    ifstream file;
    file.open("outFolder/level_offset.bin",ios::binary);
    for(int i=0; i<4; i++){
        int k;
        file.read((char*)&k, 4);
        cout << k << endl;
    }
    return 0;
}