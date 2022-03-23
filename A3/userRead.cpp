#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include<bits/stdc++.h>

using namespace std;

void writeAll(int n, std::ofstream &out, int size){
    for(int i =0; i<n ;i++){
        out.write((char*)&i, size);
    }
}

    
int main(int argc, char **argv) {
    // cout<<argv[1]<<endl;
    

    
    //criticalData[6] = 0;
    int rank, size;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_File input;
    int maxTextSize = 5;
    
    string fullFilePath = string(argv[1]);
    
    MPI_File_open(MPI_COMM_WORLD, fullFilePath.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &input);
    // MPI_File_open(MPI_COMM_WORLD, "vect.bin", MPI_MODE_CREATE |MPI_MODE_WRONLY, MPI_INFO_NULL, &outFile);
    

    
    MPI_Offset toStart,toEnd,filesize,toEndNoOffset;
    int curRankSize;
    char *buff;
    bool status = false;
    int start = 0;
    int realEnd;
    while (!status)
    {
        //free(buff);
        status = false;
        //status = true;
        maxTextSize = maxTextSize*2;
        //cout<<maxTextSize<<endl;
        MPI_File_get_size(input, &filesize);
        //cout<<" filesize "<<filesize<<endl;
        curRankSize = filesize/size;

        toStart = rank*curRankSize;
        toEnd = (maxTextSize+((rank+1)*curRankSize)) - 1;
        toEndNoOffset =((rank+1)*curRankSize) - 1;;
        if(toEnd>=filesize-1){
            status = true;
            toEnd = filesize - 1;
        }
        

        buff = (char*)malloc( ((toEnd  -toStart + 1) + 1)*sizeof(char));

        
        MPI_File_read_at_all(input, toStart, buff, (toEnd  -toStart + 1) , MPI_CHAR, MPI_STATUS_IGNORE);

        
        buff[(toEnd  -toStart + 1)] = '\0';

        realEnd = (toEnd  -toStart + 1);
        for(int i = toEndNoOffset-toStart + 1;i<((toEnd  -toStart + 1));i++){
            if(buff[i]==' '||buff[i]=='\n'){
                //cout<<"is SPace "<<rank<<endl;
                realEnd = i;
                status = true;
                break;
            }
        }
        //cout<<realEnd<<endl;
        if(rank==size-1){
            status = true;
        }

        start = 0;
        if(rank!=0){
            for(int i = 0;i<((toEnd  -toStart + 1));i++){
                //cout<<i<<" "<<buff[i]<<endl;
                if(buff[i]==' '||buff[i]=='\n'){
                    start = i+1;
                    break;
                }
            }
        }
    }

    // cout<<"At rank "<<rank<<endl;
    string str = "";
    for(int i = start;i<realEnd;i++){
        str+=buff[i];
        //cout<<buff[i];
    }

    MPI_File_close(&input);

    std::istringstream is( str);
    string output =  string(argv[2]) + "/users.bin";
    std::ofstream outfile(output, std::ios::out | std::ios::binary);
    int writeOffset = 0;
    int *allSizes = new int[size];
    int sizeoftype = sizeof(double);
    vector<double> nums;
    double n;
    while( is >> n ) {
        nums.push_back(n);
    }

    allSizes[rank] = nums.size();

    MPI_Allgather(&allSizes[rank], 1, MPI_INT, allSizes,1, MPI_INT, MPI_COMM_WORLD);

    for (int i = 0; i < rank; i++)
    {
        writeOffset+=allSizes[i];
    }
    
    if(rank == size -1){
        writeAll(writeOffset, outfile, sizeoftype);
    }


    MPI_Barrier(MPI_COMM_WORLD);
    outfile.seekp(sizeoftype*writeOffset, std::ios::beg);

    for(int i=0; i<nums.size(); i++){
        outfile.write((char*)&nums[i], sizeoftype);
    }
    outfile.close();

    MPI_Finalize();

    return 0;
}
