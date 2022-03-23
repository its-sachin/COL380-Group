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
    
    
    int rank, size;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    string fullFilePath = string(argv[1]);
    MPI_File input;
    MPI_File_open(MPI_COMM_WORLD, fullFilePath.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &input);
    MPI_Offset filesize;
    MPI_File_get_size(input, &filesize);
    MPI_File_close(&input);


    ifstream infile;
    infile.open(fullFilePath, std::ios::in);

    long long fsize = (long long)filesize;
    long long start = (rank*fsize)/size;
    long long end = ((rank+1)*fsize)/size;

    // if(rank > 0)start += 1;
    infile.seekg(start, std::ios::beg);

    if(rank > 0){
        char curr = '\0';
        while(!(curr == '\n' or curr == ' ')){
            infile.read((char*)&curr, 1);
            start += 1;
        }
    }


    string output = string(argv[2])+"/users.bin";
    std::ofstream outfile(output, std::ios::out | std::ios::binary);
    int writeOffset = 0;
    int *allSizes = new int[size];

    int sizeoftype = sizeof(double);
    vector<double> nums;
    double a;
    while(start <= end && infile >> a){  
        start += to_string(a).size() + 1;
        nums.push_back(a);
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
