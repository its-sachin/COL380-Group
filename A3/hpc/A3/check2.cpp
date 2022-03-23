#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include<bits/stdc++.h>

using namespace std;


int main(int argc, char* argv[]){
    //cout << "hello"<<endl;

    char* buff1 = new char[60];
    string s1 = "";
    for (int i = 0; i < 60; i++)
    {
        buff1[i] = '0'+i;
        // cout<<buff1[i]<<" ";
    }

    #pragma omp parallel for num_threads(4) shared(s1, buff1) reduction(+: s1)
    for (int i = 0; i < 60; i++)
    {
        s1+=buff1[i];
    }
    cout<<endl;
    cout<<s<<endl;

}